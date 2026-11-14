"""Check the compiled discretisation against the continuous joint generator."""
import numpy as np
from synthetic import Forward, PARAMS, RD, RF, shifts, B


def check():
    f=Forward(41,21,5)
    t=.7; s=.16;g=-.08
    vals=np.array([1.,.03,-.02,.01,.04])
    f.y[:]=0
    mid=(20,10,10)
    f.y[(slice(None),*mid)]=vals
    f.rhs(f.y,t,np.full(41,s),np.full(41,g))
    p=PARAMS; cd,cf=shifts(t)
    drift=np.array([RD+cd-RF-cf-.5*s*s,0.,-p.rSf*p.sf*s])
    q=np.diag([s,p.sd,p.sf])@p.correlation()@np.diag([s,p.sd,p.sf])
    xyz=np.broadcast_arrays(f.x[:,None,None],f.u[None,:,None],f.v[None,None,:])
    kill=RD+cd
    for i in range(3):
        np.testing.assert_allclose((f.out[0]*xyz[i]).sum(),drift[i],atol=1e-13)
        for j in range(3):
            np.testing.assert_allclose((f.out[0]*xyz[i]*xyz[j]).sum(),q[i,j],atol=1e-13)
    np.testing.assert_allclose(f.out[0].sum(),-kill,atol=1e-13)
    h=np.array([s,p.rSd*p.sd,p.rSf*p.sf])
    for field in (1,3):
        for i in range(3):
            np.testing.assert_allclose((f.out[field]*xyz[i]).sum(),vals[field]*(drift[i]+g*h[i]),atol=1e-13)
    vd=p.sd**2*B(2*p.kd,t);vf=p.sf**2*B(2*p.kf,t);c=p.rdf*p.sd*p.sf
    expected=np.array([-(p.kd+s*g)*vals[1]-vals[2]+vd+p.rSd*p.sd*s,
                       -(p.kd+p.kf)*vals[2]-p.rSf*p.sf*g*vals[1]+c,
                       -p.kf*vals[3]+vals[4]-vf+p.rSf*p.sf*s,
                       -(p.kd+p.kf)*vals[4]+p.rSd*p.sd*g*vals[3]+c])-kill*vals[1:]
    np.testing.assert_allclose(f.out[1:].sum(axis=(1,2,3)),expected,atol=1e-13)
    print('PASS: compiled drift, full covariance, killing, affine reactions and shared-noise flux')


if __name__=='__main__':check()
