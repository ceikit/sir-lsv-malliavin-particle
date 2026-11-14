"""Deterministic covariance regression checks; run directly (NumPy only)."""
import numpy as np
from gaussian_covariance import step_coefficients, bond_loading, step_covariance_factor
from check_forward_measure import Sigma


def run():
    count = 0
    correlation = np.array([[1,.35,-.30],[.35,1,.40],[-.30,.40,1]])
    for kd, kf in [(0.,0.),(0.,.1),(1e-8,2e-8),(1e-4,2e-4),(.05,.1),(1.,2.)]:
        for dt in [.5,.01,.0001]:
            for t0 in [0.,4.,6.-dt]:
                c=step_coefficients(t0,t0+dt,6.,kd,kf)
                cov=Sigma(c)
                np.linalg.cholesky(cov)
                factor = step_covariance_factor(t0,t0+dt,6.,kd,kf,correlation)
                np.testing.assert_allclose(factor @ factor.T, cov, rtol=2e-14, atol=2e-14)
                assert np.isfinite(cov).all()
                # Independent dense trapezoidal integration of all kernel products.
                u=np.linspace(t0,t0+dt,10001)
                ed=np.exp(-kd*(t0+dt-u));ef=np.exp(-kf*(t0+dt-u))
                ad=bond_loading(kd,6.-u);af=bond_loading(kf,6.-u)
                kernels=[ed,ed*ed,ad,ad*ad,ed*ad,ef,ef*ef,af,af*af,ef*af,ed*ef,ad*ef,ed*af,ad*af]
                keys=['J1d','J2d','I1d','I2d','Ld','J1f','J2f','I1f','I2f','Lf','Jdf','L1df','L2df','Kdf']
                for key,y in zip(keys,kernels):
                    ref=np.trapezoid(y,u)/dt
                    np.testing.assert_allclose(c[key],ref,rtol=2e-8,atol=1e-13)
                count+=1
    c=step_coefficients(0.,.01,6.,0.,0.)
    np.testing.assert_allclose(c['I1d'],5.995,rtol=1e-14)
    np.testing.assert_allclose(c['I2d'],36.-.06+.0001/3,rtol=1e-14)
    # Tiny steps: factor the rectangular Gram matrix, not its rounded covariance.
    for dt in [1e-6, 1e-8]:
        factor = step_covariance_factor(4.,4.+dt,6.,0.,0.,correlation)
        assert np.isfinite(factor).all() and (np.diag(factor) > 0).all()
        # With zero kappa, Var(Hd | spot increment, Gd) = dt^2/12.
        actual_dt = (4.+dt)-4.
        np.testing.assert_allclose(factor[2,2],actual_dt/np.sqrt(12),rtol=2e-6)
    for args in [(0,0,6,.05,.1),(0,.1,6,-.1,.1),(0,.1,6,float('nan'),.1)]:
        try:
            step_coefficients(*args)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid input accepted')
    print(f'PASS: {count} covariance cases; all coefficients, Cholesky, zero-kappa limits')

if __name__=='__main__':
    run()
