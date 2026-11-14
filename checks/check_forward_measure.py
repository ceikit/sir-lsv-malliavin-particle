"""
Check of the simulation scheme under the domestic T-forward measure Q^{d,T}
(Section 3.4 of the merged paper) and of the particle weights (Section 3.6).

 1. closed-form step coefficients (I, J, L, K, ...) against numerical quadrature;
 2. exact Gaussian-vector scheme for (f, x^d, x^f) under Q^{d,T}: martingale /
    bond-price identities at an intermediate date t < T;
 3. the sensitivity identities with Q^{d,T} particles and weights 1/P^d_{tT}.

usage: python check_forward_measure.py [steps_per_year] [pairs] [seed]
"""
import sys
import numpy as np

kd, sd = 0.05, 0.012
kf, sf = 0.10, 0.015
rSd, rSf, rdf = 0.35, -0.30, 0.40
fd, ff = 0.03, 0.01
s0, a, b = 0.12, 0.6, 2.0
S0 = 1.0
T, t_eval = 6.0, 4.0
spy = int(sys.argv[1]) if len(sys.argv) > 1 else 50


def A(k, t, TT):
    return (1 - np.exp(-k * (TT - t))) / k


def coeffs(t0, t1):
    D = t1 - t0
    c = {}
    for nm, k in (("d", kd), ("f", kf)):
        c["J1" + nm] = (1 - np.exp(-k * D)) / (k * D)
        c["J2" + nm] = (1 - np.exp(-2 * k * D)) / (2 * k * D)
        c["I1" + nm] = (k * D - np.exp(-k * T) * (np.exp(k * t1) - np.exp(k * t0))) / (k * k * D)
        c["I2" + nm] = (k * D - 2 * np.exp(-k * T) * (np.exp(k * t1) - np.exp(k * t0))
                        + 0.5 * np.exp(-2 * k * T) * (np.exp(2 * k * t1) - np.exp(2 * k * t0))) / (k ** 3 * D)
        c["L" + nm] = c["J1" + nm] / k - np.exp(-k * (T + t1)) * (np.exp(2 * k * t1) - np.exp(2 * k * t0)) / (2 * k * k * D)
    ks = kd + kf
    c["Jdf"] = (1 - np.exp(-ks * D)) / (ks * D)
    c["L1df"] = c["J1f"] / kd - np.exp(-(kd * T + kf * t1)) * (np.exp(ks * t1) - np.exp(ks * t0)) / (kd * ks * D)
    c["L2df"] = c["J1d"] / kf - np.exp(-(kf * T + kd * t1)) * (np.exp(ks * t1) - np.exp(ks * t0)) / (kf * ks * D)
    c["Kdf"] = (D - np.exp(-kd * T) * (np.exp(kd * t1) - np.exp(kd * t0)) / kd
                - np.exp(-kf * T) * (np.exp(kf * t1) - np.exp(kf * t0)) / kf
                + (np.exp(-ks * (T - t1)) - np.exp(-ks * (T - t0))) / ks) / (kd * kf * D)
    return c


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


def Sigma(c):
    return np.array([
        [1.0,            rSd * c["J1d"], rSd * c["I1d"],  rSf * c["J1f"],  rSf * c["I1f"]],
        [rSd * c["J1d"], c["J2d"],       c["Ld"],         rdf * c["Jdf"],  rdf * c["L2df"]],
        [rSd * c["I1d"], c["Ld"],        c["I2d"],        rdf * c["L1df"], rdf * c["Kdf"]],
        [rSf * c["J1f"], rdf * c["Jdf"], rdf * c["L1df"], c["J2f"],        c["Lf"]],
        [rSf * c["I1f"], rdf * c["L2df"], rdf * c["Kdf"], c["Lf"],         c["I2f"]]])


def bond(k, s, f0, t, x):
    c = x + s * s / (2 * k * k) * (1 - np.exp(-k * t)) ** 2
    return np.exp(-f0 * (T - t)) * np.exp(-A(k, t, T) * c
                                          - s * s / (4 * k ** 3) * (np.exp(-k * T) - np.exp(-k * t)) ** 2 * (np.exp(2 * k * t) - 1)), c


def pm(x):
    h = len(x) // 2
    p = 0.5 * (x[:h] + x[h:])
    return p.mean(), p.std(ddof=1) / np.sqrt(h)


def run(N=400_000, seed=3):
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else seed
    rng = np.random.default_rng(seed)
    n = int(round(t_eval * spy))
    D = 1.0 / spy
    M = 2 * N
    f0 = S0 * np.exp(-ff * T) / np.exp(-fd * T)
    lnf = np.full(M, np.log(f0)); xd = np.zeros(M); xf = np.zeros(M)
    ld = np.zeros(M); zf = np.zeros(M); lf = np.zeros(M); zd = np.zeros(M)
    var = lambda t, k, s: s * s * (1 - np.exp(-2 * k * t)) / (2 * k)
    for j in range(n):
        t0, t1 = j * D, (j + 1) * D
        c = coeffs(t0, t1)
        ch = np.linalg.cholesky(Sigma(c))
        Z = rng.standard_normal((N, 5)) @ ch.T
        Z = np.concatenate([Z, -Z]) * np.sqrt(D)
        dW, Gd, Hd, Gf, Hf = Z.T
        Pd, _ = bond(kd, sd, fd, t0, xd)
        Pf, _ = bond(kf, sf, ff, t0, xf)
        S = np.exp(lnf) * Pd / Pf
        th = np.tanh(b * np.log(S / S0))
        sig = s0 * (1 - a * th)
        g = -s0 * a * b * (1 - th * th)
        dWq = dW - rSd * sd * D * c["I1d"]          # increment of the Q^d Brownian motion
        vd, vf = var(t0, kd, sd), var(t0, kf, sf)
        ld_n = ld * (1 + g * (dWq - sig * D) - kd * D) + (vd - zf + rSd * sd * sig) * D
        zf_n = zf + (-(kd + kf) * zf - rSf * sf * g * ld + rdf * sd * sf) * D
        lf_n = lf * (1 + g * dWq - kf * D) + (zd - vf + rSf * sf * sig) * D
        zd_n = zd + (-(kd + kf) * zd + rSd * sd * g * lf + rdf * sd * sf) * D
        ld, zf, lf, zd = ld_n, zf_n, lf_n, zd_n
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
    run(N=int(sys.argv[2]) if len(sys.argv) > 2 else 400_000)
