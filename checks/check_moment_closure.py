"""Checks for the exact discounted first-moment closure.

The symbolic check starts with the full seven-dimensional generator. The
simulation checks integrated weak balances for a prescribed nonflat volatility;
it is not a PDE calibration solver. Dependencies: numpy and sympy.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp


def symbolic_checks():
    x, u, v = z = sp.symbols("x u v", real=True)
    y = sp.symbols("y1:5", real=True)
    kd, kf, nd, nf, bd, bf, c, vd, vf, ed, ef = sp.symbols(
        "kd kf nd nf beta gamma c vd vf eta_d eta_f", real=True
    )
    s = sp.Function("s")(x)
    g = sp.diff(s, x)
    b = sp.Matrix([ed + u - ef - v - s*s/2, -kd*u, -kf*v - bf*s])
    a = sp.Matrix([[-kd-s*g, -1, 0, 0], [-bf*g, -kd-kf, 0, 0],
                   [0, 0, -kf, 1], [0, 0, bd*g, -kd-kf]])
    source = sp.Matrix([vd+bd*s, c, -vf+bf*s, c])
    # Build joint quadratic variation directly from the three correlated drivers.
    drivers = sp.Matrix([[s, 0, 0], [0, nd, 0], [0, 0, nf],
                         [g*y[0], 0, 0], [0, 0, 0],
                         [g*y[2], 0, 0], [0, 0, 0]])
    corr = sp.Matrix([[1, bd/nd, bf/nf], [bd/nd, 1, c/(nd*nf)],
                      [bf/nf, c/(nd*nf), 1]])
    cov = drivers*corr*drivers.T
    joint_b = b.col_join(a*sp.Matrix(y)+source)
    coords = z+y
    phis = [sp.Function(f"phi{i}")(*z) for i in range(5)]
    test = phis[0]+sum(y[i]*phis[i+1] for i in range(4))
    full = sum(joint_b[i]*sp.diff(test, coords[i]) for i in range(7))
    full += sum(cov[i,j]*sp.diff(test, coords[i],coords[j])/2
                for i in range(7) for j in range(7))-(ed+u)*test

    def killed(phi):
        return (sum(b[i]*sp.diff(phi,z[i]) for i in range(3))
                +sum(cov[i,j]*sp.diff(phi,z[i],z[j])/2
                     for i in range(3) for j in range(3))-(ed+u)*phi)

    expected = killed(phis[0])+sum(source[i]*phis[i+1] for i in range(4))
    for j in range(4):
        coefficient = killed(phis[j+1])+sum(a[i,j]*phis[i+1] for i in range(4))
        if j in (0,2):
            coefficient += g*(s*sp.diff(phis[j+1],x)
                              +bd*sp.diff(phis[j+1],u)+bf*sp.diff(phis[j+1],v))
        expected += y[j]*coefficient
    assert sp.expand(full-expected) == 0
    assert all(sp.diff(full,yi,yj)==0 for yi in y for yj in y)

    # Independently integrate the physical density equation and eliminate the
    # rate marginals using the two own-measure integration-by-parts identities.
    p, n1, n3, star, q = [sp.Function(name)(x) for name in ("p","n1","n3","star","q")]
    rd, rf = sp.symbols("rd0 rf0")
    jd, jf = -sp.diff(n1,x), -sp.diff(n3,x)-n3
    integrated = (-sp.diff((rd-rf-s*s/2)*p+jd-jf,x)
                  +sp.diff(s*s*p,x,2)/2-rd*p-jd)

    def base(w):
        return -(rd-rf)*sp.diff(w,x)-rd*w

    def variance_term(w):
        return (sp.diff(w,x,2)+sp.diff(w,x))/2

    compact = base(p)+variance_term(s*s*p+2*(n1-n3))
    assert sp.expand(integrated-compact) == 0
    calibrated = compact.subs(n1,n3+(q-s*s)*p/2).doit()
    assert sp.simplify(calibrated-base(p)-variance_term(q*p)) == 0
    known = compact.subs(n1,n3+(q-s*s)*star/2).doit()
    difference = known-base(star)-variance_term(q*star)
    assert sp.simplify(difference-base(p-star)-variance_term(s*s*(p-star))) == 0
    return {"joint_generator_identity": True, "no_second_moments": True,
            "marginal_equation": True, "model_denominator_reduction": True,
            "market_denominator_difference_equation": True}


def monte_carlo(paths, steps, horizon, seed):
    rng = np.random.default_rng(seed)
    kd, kf, nd, nf = .07, .13, .04, .05
    rsd, rsf, rdf = .36, -.28, .25
    beta, gamma, c = rsd*nd, rsf*nf, rdf*nd*nf
    corr = np.array([[1,rsd,rsf],[rsd,1,rdf],[rsf,rdf,1]])
    root = np.linalg.cholesky(corr)
    dt = horizon/steps
    z = np.zeros((3,paths))
    y = np.zeros((4,paths))
    discount = np.ones(paths)
    integrated = np.zeros((5,paths))
    cross_integrated = np.zeros((2,paths))
    # Different affine state test functions activate every state/driver channel.
    grad = np.array([[1.,.3,-.2],[1.,2.,-3.],[1.,-1.,.5],
                     [1.,-2.,2.],[1.,1.,1.]])

    def loading(k,t):
        return -np.expm1(-k*t)/k

    for j in range(steps):
        t = j*dt
        tanh = np.tanh(2*z[0])
        s = .18*(1-.6*tanh)
        g = -.216*(1-tanh*tanh)
        rd = .025+.5*nd*nd*loading(kd,t)**2+z[1]
        rf = .015+.5*nf*nf*loading(kf,t)**2+z[2]
        b = np.array([rd-rf-.5*s*s, -kd*z[1], -kf*z[2]-gamma*s])
        yd = np.array([(-kd-s*g)*y[0]-y[1]+nd*nd*loading(2*kd,t)+beta*s,
                       -(kd+kf)*y[1]-gamma*g*y[0]+c,
                       -kf*y[2]+y[3]-nf*nf*loading(2*kf,t)+gamma*s,
                       -(kd+kf)*y[3]+beta*g*y[2]+c])
        phi = grad@z
        lphi = grad@b-rd*phi
        rhs = np.empty((5,paths))
        rhs[0] = lphi[0]
        rhs[1:] = yd*phi[1:]+y*lphi[1:]
        for jcross, i in enumerate((0,2)):
            cross = g*y[i]*(grad[i+1,0]*s+grad[i+1,1]*beta+grad[i+1,2]*gamma)
            rhs[i+1] += cross
            cross_integrated[jcross] += dt*discount*cross
        integrated += dt*discount*rhs
        dw = np.sqrt(dt)*(root@rng.standard_normal((3,paths)))
        y[[0,2]] += (g*y[[0,2]])*dw[0]
        y += dt*yd
        z += dt*b+np.array([s*dw[0],nd*dw[1],nf*dw[2]])
        discount *= np.exp(-rd*dt)

    phi = grad@z
    terminal = discount*np.vstack((phi[0],y*phi[1:]))
    residual = terminal-integrated
    rows = []
    for i, name in enumerate(("p","m1","m2","m3","m4")):
        mean = float(residual[i].mean())
        se = float(residual[i].std(ddof=1)/np.sqrt(paths))
        rows.append({"field":name, "weak_residual":mean, "standard_error":se,
                     "z_score":mean/se})
    omitted = []
    for jcross,i in enumerate((1,3)):
        wrong = residual[i]+cross_integrated[jcross]
        omitted.append({"field":f"m{i}","omitted_flux_z_score":float(
            wrong.mean()/(wrong.std(ddof=1)/np.sqrt(paths))),
            "mean_shared_noise_integral":float(cross_integrated[jcross].mean())})
    # Check the full-state Malliavin compatibility identities independently of
    # the time-integrated generator balances above.
    cd = z[1]+.5*nd*nd*loading(kd,horizon)**2
    cf = z[2]+.5*nf*nf*loading(kf,horizon)**2
    vd, vf = nd*nd*loading(2*kd,horizon), nf*nf*loading(2*kf,horizon)
    domestic = discount*(cd*phi[1]-(y[0]*grad[1,0]+vd*grad[1,1]+y[1]*grad[1,2]))
    foreign = discount*np.exp(z[0])*(cf*phi[3]-(y[2]*grad[3,0]+y[3]*grad[3,1]+vf*grad[3,2]))
    compatibility=[]
    for name, error in (("domestic",domestic),("foreign",foreign)):
        se=float(error.std(ddof=1)/np.sqrt(paths))
        compatibility.append({"identity":name,"mean_error":float(error.mean()),
                              "standard_error":se,"z_score":float(error.mean()/se)})
    # Monte Carlo bounds are a reproducibility check, not an exact discretization bound.
    assert max(abs(row["z_score"]) for row in rows)<5, rows
    assert max(abs(row["z_score"]) for row in compatibility)<5, compatibility
    return {"paths":paths,"steps":steps,"horizon":horizon,"seed":seed,
            "volatility":"0.18*(1-0.6*tanh(2*log(S)))",
            "weak_balances":rows,"omitted_flux_comparison":omitted,
            "full_state_compatibility":compatibility,
            "scope":"Prescribed nonflat Euler simulation; includes time discretization error."}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paths",type=int,default=120000)
    parser.add_argument("--steps",type=int,default=600)
    parser.add_argument("--horizon",type=float,default=3.)
    parser.add_argument("--seed",type=int,default=271828)
    parser.add_argument("--symbolic-only",action="store_true")
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result={"symbolic":symbolic_checks()}
    print("PASS: full generator, affine closure, marginal and both feedback reductions",flush=True)
    if not args.symbolic_only:
        result["monte_carlo"]=monte_carlo(args.paths,args.steps,args.horizon,args.seed)
    text=json.dumps(result,indent=2)
    print(text,flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text+"\n")


if __name__=="__main__":
    main()
