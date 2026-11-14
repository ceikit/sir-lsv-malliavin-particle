"""Shared synthetic setup, positive conservative PDE and feedback.

Prices are spot-option prices. No CME data or futures-to-spot approximation is
used here. Grid arrays store cell masses, so moment ratios have no cell-volume
factor. Tail integrals need an extra dx when divided by marginal cell masses.
"""
import ctypes
import dataclasses
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sirlsv.model import RateParams
from sirlsv.sensitivities import sensitivity_step

PARAMS = RateParams(.07, .02, .13, .025, .25, -.20, .15)
RD, RF, TMAX = .025, .015, 2.0
PAR = np.array([*dataclasses.astuple(PARAMS), RD, RF], dtype=float)
VMIN, VMAX = .06**2, .24**2


def B(k, t):
    return t if k == 0 else -np.expm1(-k*t)/k


def shifts(t):
    p = PARAMS
    return .5*p.sd**2*B(p.kd,t)**2, .5*p.sf**2*B(p.kf,t)**2


def true_vol(t, x):
    # A time-dependent, genuinely nonflat bounded coefficient, not supplied to
    # the calibrators. Level and skew both change with time.
    level=.16*(1+.08*np.sin(np.pi*t/TMAX))
    a=.35+.04*np.sin(np.pi*t/TMAX)
    th=np.tanh(2*x)
    return level*(1-a*th), -2*level*a*(1-th**2)


def native():
    src=Path(__file__).with_name('forward.cpp')
    lib=src.with_suffix('.so')
    if not lib.exists() or lib.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(src),'-o',str(lib)],check=True)
    fn=ctypes.CDLL(str(lib)).rhs
    arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    fn.argtypes=[arr,arr,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,
                 ctypes.c_double,ctypes.c_double,ctypes.c_double,arr,arr,arr,ctypes.c_double]
    fn.restype=ctypes.c_double
    return fn


class Forward:
    def __init__(self, nx=81, nr=41, nf=1, width=1.2):
        self.x=np.linspace(-width,width,nx)
        self.u=np.linspace(-.1*width,.1*width,nr)
        self.v=np.linspace(-.125*width,.125*width,nr)
        self.hx=self.x[1]-self.x[0]
        self.hu=self.u[1]-self.u[0]; self.hv=self.v[1]-self.v[0]
        self.y=np.zeros((nf,nx,nr,nr))
        self.y[0,nx//2,nr//2,nr//2]=1
        self.tmp=np.empty_like(self.y); self.out=np.empty_like(self.y)
        self.fn=native(); self.max_cfl=0.

    def rhs(self, y, t, sig, g):
        nf,nx,nu,nv=y.shape
        rate=self.fn(y,self.out,nf,nx,nu,nv,self.hx,self.hu,self.hv,
                     np.ascontiguousarray(sig),np.ascontiguousarray(g),PAR,t)
        if rate<0: raise ValueError('Grid violates positive diagonal-diffusion decomposition')
        return rate

    def step(self,t,dt,sig,g):
        # SSPRK2, coefficient frozen over the step, sources/curves at each stage.
        rate=self.rhs(self.y,t,sig,g)
        if dt*rate>.95: raise ValueError(f'CFL exceeded: {dt*rate}')
        self.max_cfl=max(self.max_cfl,dt*rate)
        np.add(self.y,dt*self.out,out=self.tmp)
        self.rhs(self.tmp,t+dt,sig,g)
        self.y*=.5
        self.y+=.5*(self.tmp+dt*self.out)
        if not np.isfinite(self.y).all(): raise FloatingPointError('Nonfinite PDE field')

    def marginal(self): return self.y[0].sum(axis=(1,2))

    def correction(self,t,method):
        p=self.marginal()
        if method=='moments':
            return 2*(self.y[1]-self.y[3]).sum(axis=(1,2))/np.maximum(p,1e-300)
        cd,cf=shifts(t)
        jd=(self.y[0]*(self.u[None,:,None]+cd)).sum(axis=(1,2))
        jf=(self.y[0]*(self.v[None,None,:]+cf)).sum(axis=(1,2))
        def tail(z): return np.cumsum(z[::-1])[::-1]-.5*z
        # Use the representation with the shorter tail on either side of zero.
        # Exact identities E[D c_d]=E[D S c_f]=0 justify lower-tail integration.
        def balanced_tail(z):
            upper=tail(z); lower=-np.cumsum(z)+.5*z
            return np.where(self.x<(RD-RF)*t,lower,upper)
        n1=balanced_tail(jd)
        n3=np.exp(-self.x)*balanced_tail(np.exp(self.x)*jf)
        return 2*self.hx*(n1-n3)/np.maximum(p,1e-300)

    def diagnostics(self,t):
        p=self.y[0]; px=self.marginal(); cd,cf=shifts(t)
        edge=np.zeros_like(p,dtype=bool)
        edge[:2]=True;edge[-2:]=True;edge[:,:2]=True;edge[:,-2:]=True;edge[:,:,:2]=True;edge[:,:,-2:]=True
        return dict(mass=float(px.sum()),mass_error=float(px.sum()-np.exp(-RD*t)),
                    stock_error=float(px@np.exp(self.x)-np.exp(-RF*t)),min_mass=float(p.min()),
                    boundary_mass=float(p[edge].sum()),max_cfl=self.max_cfl,
                    domestic_center_error=float((p*(self.u[None,:,None]+cd)).sum()),
                    foreign_center_error=float((p*np.exp(self.x[:,None,None])*(self.v[None,None,:]+cf)).sum()))


class Target:
    def __init__(self,path):
        z=np.load(path)
        self.t=z['t'];self.x=z['x'];self.q=z['q'];self.mass=z['mass']; self.calls=z['calls'];self.strikes=z['strikes']; self.price_t=z['price_t']

    def variance(self,t,x):
        j=min(np.searchsorted(self.t,t,side='right')-1,len(self.t)-2);j=max(j,0)
        a=(t-self.t[j])/(self.t[j+1]-self.t[j])
        return np.interp(x,self.x,(1-a)*self.q[j]+a*self.q[j+1])

    def density(self,t,x):
        j=np.argmin(abs(self.t-t))
        return np.interp(x,self.x,self.mass[j]/(self.x[1]-self.x[0]))


def safe_surface(q, correction, x, weights, t, rel_cut=1e-4):
    # Shared denominator-support rule for all methods. Outside the observed
    # support extrapolate the *correction*, rather than use the hidden truth.
    ok=weights>rel_cut*weights.max()
    if not ok.any(): raise ValueError('Empty feedback support')
    corr=np.interp(x,x[ok],correction[ok])
    raw=q-corr
    clipped=np.clip(raw,VMIN,VMAX)
    sig=np.sqrt(clipped)
    spline=PchipInterpolator(x,sig,extrapolate=False)
    g=spline.derivative()(x)
    return sig,g,dict(floor_nodes=int(np.count_nonzero(ok&(raw<VMIN))),
                      cap_nodes=int(np.count_nonzero(ok&(raw>VMAX))),support_nodes=int(ok.sum()))


def run_pde(method,target,nx,nr,dt,out,reference=False,width=1.2):
    clock=time.perf_counter()
    f=Forward(nx,nr,5 if method=='moments' else 1,width)
    n=int(round(TMAX/dt));dt=TMAX/n; times=np.linspace(0,TMAX,n+1)
    tab=[];qs=[];masses=[];floor=cap=0;prices=[];diags=[]
    price_t=np.array([.5,1.,1.5,2.]); strikes=np.exp(np.linspace(-.3,.3,13))
    for j,t in enumerate(times):
        px=f.marginal()
        if reference:
            sig,g=true_vol(t,f.x)
            corr=f.correction(t,'density') if j else np.zeros(nx)
            ok=px>1e-6*px.max()
            # Definition of continuous Dupire variance from the fine physical
            # density's rate-tail correction. No numerical time derivative of
            # noisy prices. This exposes the target-generator discretisation.
            q=sig**2+corr
            q=np.interp(f.x,f.x[ok],q[ok])
            if np.any(q<=0): raise ValueError('Synthetic Dupire variance nonpositive')
            qs.append(q);masses.append(px.copy())
        else:
            q=target.variance(t,f.x)
            corr=f.correction(t,method) if j else np.zeros(nx)
            sig,g,diag=safe_surface(q,corr,f.x,px,t)
            floor+=diag['floor_nodes'];cap+=diag['cap_nodes']
        tab.append(sig)
        if np.any(abs(price_t-t)<1e-10):
            prices.append(px@np.maximum(np.exp(f.x[:,None])-strikes,0))
            diags.append(dict(t=float(t),**f.diagnostics(t)))
        if j<n:f.step(t,dt,sig,g)
        if j and j%max(n//4,1)==0:print(f'  {method} {nx}x{nr}x{nr} t={t:.2f} elapsed={time.perf_counter()-clock:.1f}s',flush=True)
    elapsed=time.perf_counter()-clock
    np.savez_compressed(out,t=times,x=f.x,sigma=np.array(tab),q=np.array(qs),mass=np.array(masses),
                        calls=prices,price_t=price_t,strikes=strikes,runtime=elapsed)
    return dict(method='reference' if reference else method,nx=nx,nr=nr,dt=dt,
                runtime=elapsed,floor_nodes=floor,cap_nodes=cap,diagnostics=diags)


def deposit(x, values, grid):
    """Linear deposition, retaining total weight even outside the tabulated box."""
    z=np.clip((x-grid[0])/(grid[1]-grid[0]),0,len(grid)-1-1e-12)
    j=z.astype(int); a=z-j
    return (np.bincount(j,weights=values*(1-a),minlength=len(grid))+
            np.bincount(j+1,weights=values*a,minlength=len(grid)))


def run_particle(target,particles,dt,seed,out,bandwidth_scale=1.2):
    clock=time.perf_counter();p=PARAMS;n=int(round(TMAX/dt));dt=TMAX/n
    times=np.linspace(0,TMAX,n+1);grid=np.linspace(-1.2,1.2,401)
    rng=np.random.default_rng(seed);chol=np.linalg.cholesky(p.correlation())
    x=np.zeros(particles);u=x.copy();v=x.copy();logd=x.copy()
    ld=x.copy();zf=x.copy();lf=x.copy();zd=x.copy()
    tabs=[];floor=cap=0;ess_min=[]
    for j,t in enumerate(times):
        D=np.exp(logd)
        h=max(bandwidth_scale*.16*np.sqrt(max(t,1e-4))*particles**(-.2),grid[1]-grid[0])
        smooth=lambda z:gaussian_filter1d(z,h/(grid[1]-grid[0]),mode='constant',truncate=5)
        den=smooth(deposit(x,D,grid))
        num=smooth(deposit(x,D*(ld-lf),grid))
        corr=2*num/np.maximum(den,1e-300)
        sig,g,diag=safe_surface(target.variance(t,grid),corr,grid,den,t)
        floor+=diag['floor_nodes'];cap+=diag['cap_nodes'];tabs.append(sig)
        if any(abs(t-T)<1e-10 for T in [.5,1.,1.5,2.]):
            ess=[]
            for zz in [-2.,-1.,0.,1.,2.]:
                xx=(RD-RF)*t+zz*.16*np.sqrt(t)
                w=D*np.exp(-.5*((x-xx)/h)**2)
                ess.append(float(w.sum()**2/max(w@w,1e-300)))
            ess_min.append(dict(t=float(t),z=[-2.,-1.,0.,1.,2.],effective_particles=ess))
        if j==n:break
        sp=np.interp(x,grid,sig);gp=np.interp(x,grid,g)
        z=rng.standard_normal((particles//2,3))@chol.T
        dw=np.concatenate((z,-z))*np.sqrt(dt)
        ld,zf,lf,zd=sensitivity_step(ld,zf,lf,zd,gp,sp,dw[:,0],dt,
                                    p.sd**2*B(2*p.kd,t),p.sf**2*B(2*p.kf,t),p)
        cd,cf=shifts(t)
        logd-=(RD+cd+u)*dt
        x+=(RD-RF+cd-cf+u-v-.5*sp**2)*dt+sp*dw[:,0]
        u+=-p.kd*u*dt+p.sd*dw[:,1]
        v+=(-p.kf*v-p.rSf*p.sf*sp)*dt+p.sf*dw[:,2]
        if j and j%max(n//4,1)==0:print(f'  particle N={particles} seed={seed} t={t:.2f}',flush=True)
    elapsed=time.perf_counter()-clock
    np.savez_compressed(out,t=times,x=grid,sigma=np.array(tabs),runtime=elapsed)
    return dict(method='particle',particles=particles,dt=dt,seed=seed,bandwidth_scale=bandwidth_scale,
                runtime=elapsed,floor_nodes=floor,cap_nodes=cap,effective_sample_size=ess_min)


class Surface:
    def __init__(self,path):
        z=np.load(path);self.t=z['t'];self.x=z['x'];self.s=z['sigma']

    def __call__(self,t,x):
        j=min(np.searchsorted(self.t,t,side='right')-1,len(self.t)-2);j=max(0,j)
        a=(t-self.t[j])/(self.t[j+1]-self.t[j])
        return np.interp(x,self.x,(1-a)*self.s[j]+a*self.s[j+1])


def reprice(surfaces,particles=100000,dt=.002,seed=98181):
    """Independent paired-path repricing against the *generating* model.

    All frozen calibrated surfaces and the truth use identical random numbers,
    independent of calibration. Pairing estimates the error with much lower
    variance than comparing two separately noisy prices. No in-sample fitting.
    """
    from scipy.special import ndtr
    p=PARAMS;n=round(TMAX/dt);dt=TMAX/n;N=particles
    labels=['truth']+list(surfaces);count=len(labels)
    x=np.zeros((count,N));u=np.zeros(N);v=x.copy();logd=np.zeros(N)
    rng=np.random.default_rng(seed);chol=np.linalg.cholesky(p.correlation())
    mats=[.5,1.,1.5,2.]; rows=[]
    for j in range(n):
        t=j*dt
        sig=np.vstack([true_vol(t,x[0])[0]]+[s(t,x[k+1]) for k,s in enumerate(surfaces.values())])
        z=rng.standard_normal((N//2,3))@chol.T;dw=np.concatenate([z,-z])*np.sqrt(dt)
        cd,cf=shifts(t)
        logd-=(RD+cd+u)*dt
        x+=(RD-RF+cd-cf+u[None,:]-v-.5*sig**2)*dt+sig*dw[:,0]
        u+=-p.kd*u*dt+p.sd*dw[:,1]
        v+=(-p.kf*v-p.rSf*p.sf*sig)*dt+p.sf*dw[:,2]
        tt=(j+1)*dt
        if any(abs(tt-T)<1e-10 for T in mats):
            S=np.exp(x);D=np.exp(logd)
            F=np.exp((RD-RF)*tt)
            for zstrike in [-1.5,-.75,0.,.75,1.5]:
                K=F*np.exp(zstrike*.16*np.sqrt(tt))
                cp=1 if K>=F else -1
                payoff=D[None,:]*np.maximum(cp*(S-K),0)
                pairs=.5*(payoff[:,:N//2]+payoff[:,N//2:])
                # Vega scale only, not an assumed target implied volatility.
                d1=(np.log(F/K)+.5*.16**2*tt)/(.16*np.sqrt(tt))
                vega=np.exp(-RD*tt)*F*np.exp(-.5*d1*d1)/np.sqrt(2*np.pi)*np.sqrt(tt)
                for k,label in enumerate(labels):
                    diff=pairs[k]-pairs[0]
                    rows.append(dict(label=label,t=tt,K=K,zstrike=zstrike,price=float(pairs[k].mean()),
                        price_se=float(pairs[k].std(ddof=1)/np.sqrt(N//2)),truth_price=float(pairs[0].mean()),
                        paired_error=float(diff.mean()),paired_se=float(diff.std(ddof=1)/np.sqrt(N//2)),
                        error_vol_points=float(100*diff.mean()/vega),se_vol_points=float(100*diff.std(ddof=1)/np.sqrt(N//2)/vega)))
            print(f'  independent repricing t={tt:.1f}, {count} models, N={N}',flush=True)
    return rows
