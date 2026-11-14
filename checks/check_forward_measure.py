"""
Check of the simulation scheme under the domestic T-forward measure Q^{d,T}
("Time Discretization under Q^{d^T}" in the paper) and of the particle weights.

 1. stable step coefficients (I, J, L, K, ...) against independent quadrature;
 2. exact Gaussian-vector scheme for (f, x^d, x^f) under Q^{d,T}: martingale /
    bond-price identities at an intermediate date t < T;
 3. the sensitivity identities with Q^{d,T} particles and weights 1/P^d_{tT}.

usage: python checks/check_forward_measure.py [steps_per_year] [pairs] [seed]
"""
import sys

import numpy as np

import _repo  # noqa: F401  (import path)
from sirlsv.gaussian_covariance import bond_loading, step_coefficients, step_covariance_factor
from sirlsv.model import tanh_local_vol
from sirlsv.montecarlo import antithetic_mean as pm
from sirlsv.prescribed_vol import FD as fd, FF as ff, PARAMS, S0, VOL_B as b, VOL_LEVEL as s0
from sirlsv.sensitivities import sensitivity_step

kd, sd, kf, sf = PARAMS.kd, PARAMS.sd, PARAMS.kf, PARAMS.sf
rSd, rSf, rdf = PARAMS.rSd, PARAMS.rSf, PARAMS.rdf
a = 0.6                      # skew of the local volatility
T, t_eval = 6.0, 4.0


def A(k, t, TT):
    return bond_loading(k, TT - t)


def coeffs(t0, t1):
    return step_coefficients(t0, t1, T, kd, kf)


def quad_check():
    t0, t1 = 2.3, 2.8
    D = t1 - t0
    u = np.linspace(t0, t1, 200001)
    integ = lambda y: float(np.sum((y[1:] + y[:-1]) * 0.5) * (u[1] - u[0])) / D
    ed, ef = np.exp(-kd * (t1 - u)), np.exp(-kf * (t1 - u))
    Ad, Af = A(kd, u, T), A(kf, u, T)
    ref = {"J1d": integ(ed), "J2d": integ(ed * ed), "I1d": integ(Ad), "I2d": integ(Ad * Ad), "Ld": integ(ed * Ad),
           "J1f": integ(ef), "J2f": integ(ef * ef), "I1f": integ(Af), "I2f": integ(Af * Af), "Lf": integ(ef * Af),
           "Jdf": integ(ed * ef), "L1df": integ(Ad * ef), "L2df": integ(ed * Af), "Kdf": integ(Ad * Af)}
    c = coeffs(t0, t1)
    worst = max(abs(c[k] - ref[k]) / abs(ref[k]) for k in ref)
    print(f"1. step coefficients vs quadrature: max relative error = {worst:.1e}")


def bond(k, s, f0, t, x):
    c = x + .5 * s * s * bond_loading(k, t) ** 2
    variance = s * s * bond_loading(2 * k, t)
    loading = A(k, t, T)
    return np.exp(-f0 * (T-t) - loading*c - .5*loading*loading*variance), c


def run(N=400_000, seed=3, spy=50):
    rng = np.random.default_rng(seed)
    n = int(round(t_eval * spy))
    D = 1.0 / spy
    M = 2 * N
    f0 = S0 * np.exp(-ff * T) / np.exp(-fd * T)
    lnf = np.full(M, np.log(f0)); xd = np.zeros(M); xf = np.zeros(M)
    ld = np.zeros(M); zf = np.zeros(M); lf = np.zeros(M); zd = np.zeros(M)
    var = lambda t, k, s: s * s * bond_loading(2 * k, t)
    for j in range(n):
        t0, t1 = j * D, (j + 1) * D
        c = coeffs(t0, t1)
        ch = step_covariance_factor(t0, t1, T, kd, kf, PARAMS.correlation())
        Z = rng.standard_normal((N, 5)) @ ch.T
        Z = np.concatenate([Z, -Z]) * np.sqrt(D)
        dW, Gd, Hd, Gf, Hf = Z.T
        Pd, _ = bond(kd, sd, fd, t0, xd)
        Pf, _ = bond(kf, sf, ff, t0, xf)
        S = np.exp(lnf) * Pd / Pf
        sig, g = tanh_local_vol(np.log(S / S0), s0, a, b)
        dWq = dW - rSd * sd * D * c["I1d"]          # increment of the Q^d Brownian motion
        vd, vf = var(t0, kd, sd), var(t0, kf, sf)
        ld, zf, lf, zd = sensitivity_step(ld, zf, lf, zd, g, sig, dWq, D, vd, vf, PARAMS)
        lnf += (sig * dW + sd * Hd - sf * Hf
                - 0.5 * D * (sig * sig + 2 * rSd * sig * sd * c["I1d"] - 2 * rSf * sig * sf * c["I1f"]
                             + sd * sd * c["I2d"] + sf * sf * c["I2f"] - 2 * rdf * sd * sf * c["Kdf"]))
        xd_n = np.exp(-kd * D) * xd + sd * Gd - sd * sd * D * c["Ld"]
        xf = np.exp(-kf * D) * xf + sf * Gf - rSf * sf * sig * D * c["J1f"] - rdf * sd * sf * D * c["L1df"]
        xd = xd_n

    Pd, cd = bond(kd, sd, fd, t_eval, xd)
    Pf, cf = bond(kf, sf, ff, t_eval, xf)
    S = np.exp(lnf) * Pd / Pf
    w = np.exp(-fd * T) / Pd                         # P^d_{0T} / P^d_{tT}
    print(f"2. Q^(d,T) scheme, T={T}, t={t_eval}, paths={M}, steps/yr={spy}")
    for nm, x, ref in (("E[f_t]/f_0", np.exp(lnf) / f0, 1.0),
                       ("P_0T E[1/P_tT] / P^d_0t", w / np.exp(-fd * t_eval), 1.0),
                       ("P_0T E[S_t/P_tT] / (S_0 P^f_0t)", w * S / (S0 * np.exp(-ff * t_eval)), 1.0),
                       ("P_0T E[c^d_t/P_tT]", w * cd, 0.0),
                       ("P_0T E[S_t c^f_t/P_tT]", w * S * cf, 0.0)):
        m, se = pm(x)
        print(f"   {nm:34s} = {m:+.6f} +- {se:.1e}   (exact {ref})")
    print("3. sensitivity identities with Q^(d,T) particles and weights 1/P^d_tT  (z = difference / s.e.)")
    for K in (0.7, 1.0, 1.3):
        ind = S > K
        pay = np.maximum(S - K, 0)
        a1, b1 = w * cd * pay, w * ind * S * ld
        a2, b2 = w * S * cf * pay, w * ind * S * S * lf
        d1, s1 = pm(a1 - b1)
        d2, s2 = pm(a2 - b2)
        print(f"   K={K:.2f}  domestic {pm(a1)[0]:+.4e} vs {pm(b1)[0]:+.4e} (z={d1 / s1:+.2f})"
              f"   foreign {pm(a2)[0]:+.4e} vs {pm(b2)[0]:+.4e} (z={d2 / s2:+.2f})")


if __name__ == "__main__":
    quad_check()
    run(spy=int(sys.argv[1]) if len(sys.argv) > 1 else 50,
        N=int(sys.argv[2]) if len(sys.argv) > 2 else 400_000,
        seed=int(sys.argv[3]) if len(sys.argv) > 3 else 3)
