"""
Monte Carlo check of the integrated Malliavin sensitivities used in the
three-factor FX model (local volatility + Hull-White domestic/foreign rates).

Domestic risk-neutral measure Q^d:
    dS/S  = (r^d - r^f) dt + sigma(t,S) dW
    r^i   = f^i(0,t) + s_i^2/(2 k_i^2) (1 - e^{-k_i t})^2 + x^i,   i = d, f
    dx^d  = -k_d x^d dt + s_d dB^d
    dx^f  = (-k_f x^f - rho_Sf s_f sigma(t,S)) dt + s_f dB^f

Auxiliary processes (all start at 0), g = S dsigma/dS, v_i(t) = s_i^2 (1-e^{-2 k_i t})/(2 k_i):
    dl_d = l_d [ g (dW - sigma dt) - k_d dt ] + [ v_d - z_f + rho_Sd s_d sigma ] dt
    dz_f = [ -(k_d+k_f) z_f - rho_Sf s_f g l_d + rho_df s_d s_f ] dt
    dl_f = l_f [ g dW - k_f dt ]               + [ z_d - v_f + rho_Sf s_f sigma ] dt
    dz_d = [ -(k_d+k_f) z_d + rho_Sd s_d g l_f + rho_df s_d s_f ] dt

Identities tested (integrated in strike; c^i = r^i_T - f^i(0,T), D = exp(-int r^d)):
    (I)   E[ D c^d (S_T-K)^+ ]     = E[ D 1_{S_T>K} S_T   l_d ]
    (II)  E[ D S_T c^f (S_T-K)^+ ] = E[ D 1_{S_T>K} S_T^2 l_f ]

'naive' = same recursions with the coupling terms (g l_d in z_f, g l_f in z_d)
dropped and with g (dW - sigma dt) instead of g dW in l_f.

usage: python check_sensitivities.py [pairs] [steps_per_year] [T] [skew] [rate_vol_scale] [seed]
"""
import sys
import time
import numpy as np


def simulate(N, spy, T, skew, volscale, seed=7):
    rng = np.random.default_rng(seed)
    kd, sd = 0.05, 0.012 * volscale
    kf, sf = 0.10, 0.015 * volscale
    rSd, rSf, rdf = 0.35, -0.30, 0.40
    fd, ff = 0.03, 0.01                    # flat forward curves
    s0, a, b = 0.12, skew, 2.0             # sigma(S) = s0 (1 - a tanh(b ln S)), S0 = 1
    chol = np.linalg.cholesky(np.array([[1, rSd, rSf], [rSd, 1, rdf], [rSf, rdf, 1.0]]))
    n = int(round(T * spy))
    dt = T / n
    sq = np.sqrt(dt)
    conv = lambda t, k, s: s * s / (2 * k * k) * (1 - np.exp(-k * t)) ** 2
    var = lambda t, k, s: s * s * (1 - np.exp(-2 * k * t)) / (2 * k)

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
        th = np.tanh(b * lnS)
        sig = s0 * (1 - a * th)
        g = -s0 * a * b * (1 - th * th)
        rd = fd + conv(t, kd, sd) + xd
        rf = ff + conv(t, kf, sf) + xf
        vd, vf = var(t, kd, sd), var(t, kf, sf)

        ld_new = ld * (1 + g * (dW - sig * dt) - kd * dt) + (vd - zf + rSd * sd * sig) * dt
        zf_new = zf + (-(kd + kf) * zf - rSf * sf * g * ld + rdf * sd * sf) * dt
        lf_new = lf * (1 + g * dW - kf * dt) + (zd - vf + rSf * sf * sig) * dt
        zd_new = zd + (-(kd + kf) * zd + rSd * sd * g * lf + rdf * sd * sf) * dt

        ldn_new = ld_n * (1 + g * (dW - sig * dt) - kd * dt) + (vd - zf_n + rSd * sd * sig) * dt
        zfn_new = zf_n + (-(kd + kf) * zf_n + rdf * sd * sf) * dt
        lfn_new = lf_n * (1 + g * (dW - sig * dt) - kf * dt) + (zd_n - vf + rSf * sf * sig) * dt
        zdn_new = zd_n + (-(kd + kf) * zd_n + rdf * sd * sf) * dt

        ld, zf, lf, zd = ld_new, zf_new, lf_new, zd_new
        ld_n, zf_n, lf_n, zd_n = ldn_new, zfn_new, lfn_new, zdn_new

        intr += rd * dt
        lnS += (rd - rf - 0.5 * sig * sig) * dt + sig * dW
        xd_new = xd - kd * xd * dt + sd * dBd
        xf = xf + (-kf * xf - rSf * sf * sig) * dt + sf * dBf
        xd = xd_new

    return dict(S=np.exp(lnS), D=np.exp(-intr), cd=xd + conv(T, kd, sd), cf=xf + conv(T, kf, sf),
                ld=ld, lf=lf, ld_n=ld_n, lf_n=lf_n)


def pm(x):
    h = len(x) // 2
    p = 0.5 * (x[:h] + x[h:])
    return p.mean(), p.std(ddof=1) / np.sqrt(h)


def report(r, strikes):
    S, D = r["S"], r["D"]
    m, se = pm(D * r["cd"]);      print(f"  E[D c^d]   = {m:+.2e} +- {se:.1e}  (should be 0)")
    m, se = pm(D * S * r["cf"]);  print(f"  E[D S c^f] = {m:+.2e} +- {se:.1e}  (should be 0)")
    for name, w, sens, naive in (("(I)  domestic", D * r["cd"], D * S * r["ld"], D * S * r["ld_n"]),
                                 ("(II) foreign ", D * S * r["cf"], D * S * S * r["lf"], D * S * S * r["lf_n"])):
        print(f"  {name}:     K |      direct sensitivity      z |       naive       z")
        for K in strikes:
            ind = S > K
            a_ = w * np.maximum(S - K, 0.0)
            b_, c_ = sens * ind, naive * ind
            d1, s1 = pm(a_ - b_)
            d2, s2 = pm(a_ - c_)
            print(f"              {K:5.2f} | {pm(a_)[0]:+.4e} {pm(b_)[0]:+.4e} {d1 / s1:+6.2f} | {pm(c_)[0]:+.4e} {d2 / s2:+7.2f}")


if __name__ == "__main__":
    arg = lambda i, d, f: f(sys.argv[i]) if len(sys.argv) > i else d
    N, spy, T, skew, vs = arg(1, 400_000, int), arg(2, 100, int), arg(3, 5.0, float), arg(4, 0.6, float), arg(5, 1.0, float)
    t0 = time.time()
    res = simulate(N, spy, T, skew, vs, seed=arg(6, 7, int))
    print(f"T={T}, paths={2 * N}, steps/yr={spy}, skew={skew}, rate-vol scale={vs}, {time.time() - t0:.0f}s")
    report(res, [0.70, 0.85, 1.00, 1.15, 1.30, 1.50])
