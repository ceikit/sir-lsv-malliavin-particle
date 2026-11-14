"""
Monte Carlo check of the Malliavin (tangent-process) weights for the
3-factor FX model: local volatility + Hull-White domestic and foreign rates.

Under the domestic risk-neutral measure Q^d:
    dS/S  = (r^d - r^f) dt + sigma(t,S) dW
    r^d   = phi_d(t) + x^d,   dx^d = -k_d x^d dt + s_d dB^d
    r^f   = phi_f(t) + x^f,   dx^f = (-k_f x^f - rho_Sf s_f sigma(t,S)) dt + s_f dB^f

Identities tested (integrated-in-strike versions of the ones used in the
local-vol formula; phi(S) = (S-K)^+, phi'(S) = 1_{S>K}):

 (I)  E^d[ D (r^d_T - r^0_d) (S_T-K)^+ ]      = e^{-k_d T} E^d[ D 1_{S_T>K} S_T L_T  ]
 (II) E^d[ D S_T (r^f_T - r^0_f) (S_T-K)^+ ]  = e^{-k_f T} E^d[ D S_T 1_{S_T>K} S_T Lt_T ]

with D = exp(-int r^d), r^0_i = market forward rate f_i(0,T), and (L, Lt) the
log-tangent processes of the write-up.  "Naive" variants are run alongside to
show that the coupling terms / measure-change term are really needed.
"""
import sys
import time
import numpy as np


def simulate(N, steps_per_year, T, seed, skew=0.6, volscale=1.0):
    rng = np.random.default_rng(seed)

    # ---- parameters ---------------------------------------------------
    S0 = 1.0
    kd, sd = 0.05, 0.012 * volscale
    kf, sf = 0.10, 0.015 * volscale
    rSd, rSf, rdf = 0.35, -0.30, 0.40
    fd, ff = 0.03, 0.01                      # flat instantaneous forward curves
    s0, a, b = 0.12, skew, 2.0               # local vol: s0*(1 - a*tanh(b*ln(S/S0)))

    corr = np.array([[1.0, rSd, rSf], [rSd, 1.0, rdf], [rSf, rdf, 1.0]])
    chol = np.linalg.cholesky(corr)

    n = int(round(T * steps_per_year))
    dt = T / n
    sq = np.sqrt(dt)

    def phi(t, f0, k, s):                    # HW shift: r = phi + x under own measure
        return f0 + s * s / (2 * k * k) * (1 - np.exp(-k * t)) ** 2

    M = 2 * N                                # antithetic pairs
    lnS = np.zeros(M)
    xd = np.zeros(M)
    xf = np.zeros(M)
    intr = np.zeros(M)

    # correct tangents
    L = np.zeros(M);  Df = np.zeros(M)       # domestic bump: log-tangent of S, tangent of r^f
    Lt = np.zeros(M); Ddt = np.zeros(M)      # foreign bump:  log-tangent of S, tangent of r^d
    # naive variants
    L_nc = np.zeros(M);  Df_nc = np.zeros(M)     # domestic bump, coupling term dropped
    Lt_nv = np.zeros(M); Ddt_nv = np.zeros(M)    # foreign bump, GHL's V-dynamics + no coupling

    for j in range(n):
        t = j * dt
        Z = rng.standard_normal((N, 3)) @ chol.T
        Z = np.concatenate([Z, -Z], axis=0) * sq
        dW, dBd, dBf = Z[:, 0], Z[:, 1], Z[:, 2]

        u = b * lnS                           # ln(S/S0) with S0 = 1
        th = np.tanh(u)
        sig = s0 * (1 - a * th)
        g = -s0 * a * b * (1 - th * th)       # g = S * d sigma / dS

        rd = phi(t, fd, kd, sd) + xd
        rf = phi(t, ff, kf, sf) + xf

        ed, ef = np.exp(kd * t), np.exp(kf * t)
        Dd_det = sd * sd / kd * np.sinh(kd * t)     # tangent of r^d, domestic bump (deterministic)
        Dft_det = sf * sf / kf * np.sinh(kf * t)    # tangent of r^f, foreign bump (deterministic)

        # --- correct tangent systems (Euler) ---------------------------
        L_new = L + L * g * (dW - sig * dt) + (Dd_det - Df + rSd * sig * sd * ed) * dt
        Df_new = Df + (-kf * Df - rSf * sf * g * L + rdf * sf * sd * ed) * dt
        Lt_new = Lt + Lt * g * dW + (Ddt - Dft_det + rSf * sig * sf * ef) * dt
        Ddt_new = Ddt + (-kd * Ddt + rSd * sd * g * Lt + rdf * sd * sf * ef) * dt

        # --- naive variants --------------------------------------------
        L_nc_new = L_nc + L_nc * g * (dW - sig * dt) + (Dd_det - Df_nc + rSd * sig * sd * ed) * dt
        Df_nc_new = Df_nc + (-kf * Df_nc + rdf * sf * sd * ed) * dt
        Lt_nv_new = Lt_nv + Lt_nv * g * (dW - sig * dt) + (Ddt_nv - Dft_det + rSf * sig * sf * ef) * dt
        Ddt_nv_new = Ddt_nv + (-kd * Ddt_nv + rdf * sd * sf * ef) * dt

        L, Df, Lt, Ddt = L_new, Df_new, Lt_new, Ddt_new
        L_nc, Df_nc, Lt_nv, Ddt_nv = L_nc_new, Df_nc_new, Lt_nv_new, Ddt_nv_new

        # --- state -----------------------------------------------------
        intr += rd * dt
        lnS += (rd - rf - 0.5 * sig * sig) * dt + sig * dW
        xd_new = xd - kd * xd * dt + sd * dBd
        xf = xf + (-kf * xf - rSf * sf * sig) * dt + sf * dBf
        xd = xd_new

    S = S0 * np.exp(lnS)
    D = np.exp(-intr)
    cd = xd + sd * sd / (2 * kd * kd) * (1 - np.exp(-kd * T)) ** 2   # r^d_T - f_d(0,T)
    cf = xf + sf * sf / (2 * kf * kf) * (1 - np.exp(-kf * T)) ** 2   # r^f_T - f_f(0,T)
    wd, wf = np.exp(-kd * T), np.exp(-kf * T)
    return dict(S=S, D=D, cd=cd, cf=cf, wd=wd, wf=wf, L=L, Lt=Lt, L_nc=L_nc, Lt_nv=Lt_nv,
                Pd=np.exp(-fd * T), fwd=S0 * np.exp((fd - ff) * T))


def pair_mean(x):
    """mean and standard error using antithetic pairs as the i.i.d. unit"""
    h = len(x) // 2
    p = 0.5 * (x[:h] + x[h:])
    return p.mean(), p.std(ddof=1) / np.sqrt(h)


def report(r, strikes):
    S, D = r["S"], r["D"]
    m, se = pair_mean(D)
    print(f"  E[D]            = {m:.6f}  (P_d(0,T) = {r['Pd']:.6f})")
    m, se = pair_mean(D * S)
    print(f"  E[D S]/P_d      = {m / r['Pd']:.6f}  (forward   = {r['fwd']:.6f})")
    m, se = pair_mean(D * r["cd"])
    print(f"  E[D (r_d-r0_d)]   = {m:+.2e} +- {se:.1e}   (should be 0)")
    m, se = pair_mean(D * S * r["cf"])
    print(f"  E[D S (r_f-r0_f)] = {m:+.2e} +- {se:.1e}   (should be 0)")
    print()
    hdr = f"  {'K':>5} | {'direct':>11} {'tangent':>11} {'z':>6} | {'naive':>11} {'z':>7}"
    for name, direct_w, tan, naive in (
        ("(I)  domestic-rate term", D * r["cd"], r["wd"] * D * S * r["L"], r["wd"] * D * S * r["L_nc"]),
        ("(II) foreign-rate term ", D * S * r["cf"], r["wf"] * D * S * S * r["Lt"], r["wf"] * D * S * S * r["Lt_nv"]),
    ):
        print("  " + name)
        print(hdr)
        for K in strikes:
            a_ = direct_w * np.maximum(S - K, 0.0)
            ind = (S > K)
            b_ = tan * ind
            c_ = naive * ind
            ma, _ = pair_mean(a_)
            mb, _ = pair_mean(b_)
            mc, _ = pair_mean(c_)
            d1, s1 = pair_mean(a_ - b_)
            d2, s2 = pair_mean(a_ - c_)
            print(f"  {K:5.2f} | {ma:+.4e} {mb:+.4e} {d1 / s1:+6.2f} | {mc:+.4e} {d2 / s2:+7.2f}")
        print()


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 400_000
    spy = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    T = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0
    skew = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
    volscale = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    t0 = time.time()
    res = simulate(N, spy, T, seed=7, skew=skew, volscale=volscale)
    print(f"T = {T}, paths = {2 * N} (antithetic), steps/yr = {spy}, skew a = {skew}, rate-vol scale = {volscale}, {time.time() - t0:.0f}s")
    report(res, strikes=[0.7, 0.85, 1.0, 1.15, 1.3, 1.5])
