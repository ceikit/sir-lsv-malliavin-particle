"""Three-factor model with a prescribed local volatility, simulated under Q^d.

Used by the checks: the local volatility is given (a tanh skew), nothing is calibrated.

    dS/S  = (r^d - r^f) dt + sigma(S) dW
    r^i   = f^i + s_i^2/(2 k_i^2) (1 - e^{-k_i t})^2 + x^i,   i = d, f   (flat forward curves)
    dx^d  = -k_d x^d dt + s_d dB^d
    dx^f  = (-k_f x^f - rho_Sf s_f sigma(S)) dt + s_f dB^f
"""
import numpy as np

from .hull_white import convexity, variance
from .model import RateParams, tanh_local_vol
from .sensitivities import sensitivity_step

PARAMS = RateParams(kd=0.05, sd=0.012, kf=0.10, sf=0.015, rSd=0.35, rSf=-0.30, rdf=0.40)
FD, FF = 0.03, 0.01          # flat domestic and foreign forward rates
S0 = 1.0
VOL_LEVEL, VOL_B = 0.12, 2.0   # sigma(S) = VOL_LEVEL (1 - skew tanh(VOL_B ln(S/S0)))


def simulate(N, spy, T, skew, volscale, seed=7):
    """Euler scheme with N antithetic pairs and spy steps per year, to maturity T.

    volscale multiplies both rate volatilities. Returns terminal values: spot S, discount
    factor D, centred rates cd and cf, the sensitivities ld and lf of the paper, and the
    'naive' variants ld_n and lf_n (coupling terms dropped, and g (dW - sigma dt) used for
    the foreign one), which the checks show to be wrong.
    """
    rng = np.random.default_rng(seed)
    p = RateParams(PARAMS.kd, PARAMS.sd * volscale, PARAMS.kf, PARAMS.sf * volscale,
                   PARAMS.rSd, PARAMS.rSf, PARAMS.rdf)
    kd, sd, kf, sf, rSd, rSf, rdf = p.kd, p.sd, p.kf, p.sf, p.rSd, p.rSf, p.rdf
    chol = np.linalg.cholesky(p.correlation())
    n = int(round(T * spy))
    dt = T / n
    sq = np.sqrt(dt)

    M = 2 * N                              # antithetic pairs
    z = lambda: np.zeros(M)
    lnS, xd, xf, intr = z(), z(), z(), z()
    ld, zf, lf, zd = z(), z(), z(), z()            # the paper's processes
    ld_n, zf_n, lf_n, zd_n = z(), z(), z(), z()    # naive variants

    for j in range(n):
        t = j * dt
        Z = rng.standard_normal((N, 3)) @ chol.T
        Z = np.concatenate([Z, -Z]) * sq
        dW, dBd, dBf = Z[:, 0], Z[:, 1], Z[:, 2]
        sig, g = tanh_local_vol(lnS, VOL_LEVEL, skew, VOL_B)
        rd = FD + convexity(t, kd, sd) + xd
        rf = FF + convexity(t, kf, sf) + xf
        vd, vf = variance(t, kd, sd), variance(t, kf, sf)

        new = sensitivity_step(ld, zf, lf, zd, g, sig, dW, dt, vd, vf, p)

        ldn_new = ld_n * (1 + g * (dW - sig * dt) - kd * dt) + (vd - zf_n + rSd * sd * sig) * dt
        zfn_new = zf_n + (-(kd + kf) * zf_n + rdf * sd * sf) * dt
        lfn_new = lf_n * (1 + g * (dW - sig * dt) - kf * dt) + (zd_n - vf + rSf * sf * sig) * dt
        zdn_new = zd_n + (-(kd + kf) * zd_n + rdf * sd * sf) * dt

        ld, zf, lf, zd = new
        ld_n, zf_n, lf_n, zd_n = ldn_new, zfn_new, lfn_new, zdn_new

        intr += rd * dt
        lnS += (rd - rf - 0.5 * sig * sig) * dt + sig * dW
        xd_new = xd - kd * xd * dt + sd * dBd
        xf = xf + (-kf * xf - rSf * sf * sig) * dt + sf * dBf
        xd = xd_new

    return dict(S=np.exp(lnS), D=np.exp(-intr), cd=xd + convexity(T, kd, sd), cf=xf + convexity(T, kf, sf),
                ld=ld, lf=lf, ld_n=ld_n, lf_n=lf_n)
