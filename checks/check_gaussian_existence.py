"""Numerical cross-check of the constructive Gaussian existence proof.

Requires scipy in addition to numpy. The test checks the independently integrated
forward variance, not merely the defining ODE identity.

usage: python checks/check_gaussian_existence.py
"""
import numpy as np
from scipy.integrate import solve_ivp, quad

import _repo  # noqa: F401  (import path)
from sirlsv.gaussian_covariance import bond_loading


def run():
    kd, kf, sd, sf, rdf = .05, .1, .012, .015, .4
    horizon = 10.
    q = lambda t: .04 + .001*t
    target = lambda t: .04*t + .0005*t*t
    for rd, rf in [(.35,-.30),(0.,0.)]:
        def rhs(t,u):
            vd=sd*sd*bond_loading(2*kd,t)
            vf=sf*sf*bond_loading(2*kf,t)
            z=rdf*sd*sf*bond_loading(kd+kf,t)
            gap=q(t)-2*(u[0]-u[1])
            if gap<=0:
                raise ValueError('positive variance domain exited')
            a=np.sqrt(gap)
            return [-kd*u[0]+vd-z+rd*sd*a,-kf*u[1]+z-vf+rf*sf*a]
        solution=solve_ivp(rhs,[0,horizon],[0,0],rtol=1e-11,atol=1e-13,dense_output=True)
        assert solution.success
        def a(t):
            u=solution.sol(t)
            return np.sqrt(q(t)-2*(u[0]-u[1]))
        worst=0.
        for maturity in [1.,3.,5.,10.]:
            def integrand(s):
                ad=sd*bond_loading(kd,maturity-s)
                af=sf*bond_loading(kf,maturity-s)
                vol=a(s)
                return vol*vol+ad*ad+af*af-2*rdf*ad*af+2*vol*(rd*ad-rf*af)
            actual=quad(integrand,0,maturity,epsabs=1e-11,epsrel=1e-11)[0]
            worst=max(worst,abs(actual-target(maturity)))
            np.testing.assert_allclose(actual,target(maturity),atol=1e-10,rtol=1e-9)
            if rd==rf==0:
                ad=sd*bond_loading(kd,maturity);af=sf*bond_loading(kf,maturity)
                np.testing.assert_allclose(a(maturity)**2,q(maturity)-(ad*ad+af*af-2*rdf*ad*af),atol=1e-11)
        print(f'PASS Gaussian construction rhoSd={rd}, rhoSf={rf}: max variance error {worst:.3e}')
    # A flat 5% implied-vol surface is incompatible with these rates at 30y,
    # even allowing arbitrary spot-dependent volatility, under the log-moment assumptions.
    maturity=30.
    def rate_variance(s):
        ad=sd*bond_loading(kd,maturity-s);af=sf*bond_loading(kf,maturity-s)
        return ad*ad+af*af-2*rdf*ad*af
    lower=quad(rate_variance,0,maturity)[0]
    market=.05**2*maturity
    assert lower>market
    print(f'Log-contract obstruction (zero spot/rate correlations): target={market:.6f}, necessary minimum={lower:.6f}')

if __name__=='__main__':
    run()
