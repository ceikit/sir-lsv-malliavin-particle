"""
Cross-check: the 'integrated Malliavin' auxiliaries of the second write-up
(Y^{d,S}, Y^{d,f}, Y^{f,S}, Y^{f,d}; Euler, exactly as boxed there) against the
log-tangent processes (L^d, eta^f, L^f, eta^d) of malliavin_fx_fix.

Claim to test:  Y^{d,S}_t = e^{-k_d t} S_t L^d_t   and   Y^{f,S}_t = e^{-k_f t} S_t L^f_t
plus the integration-by-parts identities using the Y's directly.
"""
import sys
import numpy as np


def run(N=300_000, spy=100, T=5.0, seed=11):
    rng = np.random.default_rng(seed)
    kd, sd = 0.05, 0.012
    kf, sf = 0.10, 0.015
    rSd, rSf, rdf = 0.35, -0.30, 0.40
    fd, ff = 0.03, 0.01
    s0, a, b = 0.12, 0.6, 2.0
    chol = np.linalg.cholesky(np.array([[1, rSd, rSf], [rSd, 1, rdf], [rSf, rdf, 1.0]]))
    n = int(round(T * spy)); dt = T / n; sq = np.sqrt(dt)
    phi = lambda t, f0, k, s: f0 + s * s / (2 * k * k) * (1 - np.exp(-k * t)) ** 2

    M = 2 * N
    lnS = np.zeros(M); xd = np.zeros(M); xf = np.zeros(M); intr = np.zeros(M)
    L = np.zeros(M); Df = np.zeros(M); Lt = np.zeros(M); Ddt = np.zeros(M)       # mine
    YdS = np.zeros(M); Ydf = np.zeros(M); YfS = np.zeros(M); Yfd = np.zeros(M)   # theirs

    for j in range(n):
        t = j * dt
        Z = rng.standard_normal((N, 3)) @ chol.T
        Z = np.concatenate([Z, -Z]) * sq
        dW, dBd, dBf = Z[:, 0], Z[:, 1], Z[:, 2]
        S = np.exp(lnS)
        th = np.tanh(b * lnS)
        nu = s0 * (1 - a * th)
        g = -s0 * a * b * (1 - th * th)          # S * dnu/dS
        nuS = g / S
        q = nu + g
        rd = phi(t, fd, kd, sd) + xd
        rf = phi(t, ff, kf, sf) + xf
        ed, ef = np.exp(kd * t), np.exp(kf * t)

        # --- mine -------------------------------------------------------
        L_n = L + L * g * (dW - nu * dt) + (sd * sd / kd * np.sinh(kd * t) - Df + rSd * nu * sd * ed) * dt
        Df_n = Df + (-kf * Df - rSf * sf * g * L + rdf * sf * sd * ed) * dt
        Lt_n = Lt + Lt * g * dW + (Ddt - sf * sf / kf * np.sinh(kf * t) + rSf * nu * sf * ef) * dt
        Ddt_n = Ddt + (-kd * Ddt + rSd * sd * g * Lt + rdf * sd * sf * ef) * dt

        # --- theirs (boxed equations, section 5) -------------------------
        Ydd = sd * sd / (2 * kd) * (1 - np.exp(-2 * kd * t))
        Yff = sf * sf / (2 * kf) * (1 - np.exp(-2 * kf * t))
        YdS_n = YdS + ((rd - rf - kd) * YdS + S * (Ydd - Ydf) + rSd * sd * S * nu) * dt + q * YdS * dW
        Ydf_n = Ydf + (-(kd + kf) * Ydf - rSf * sf * nuS * YdS + rdf * sd * sf) * dt
        YfS_n = YfS + ((rd - rf - kf + S * nu * nuS) * YfS + S * (Yfd - Yff) + rSf * sf * S * nu) * dt + q * YfS * dW
        Yfd_n = Yfd + (-(kd + kf) * Yfd + rSd * sd * nuS * YfS + rdf * sd * sf) * dt

        L, Df, Lt, Ddt = L_n, Df_n, Lt_n, Ddt_n
        YdS, Ydf, YfS, Yfd = YdS_n, Ydf_n, YfS_n, Yfd_n

        intr += rd * dt
        lnS += (rd - rf - 0.5 * nu * nu) * dt + nu * dW
        xd_n = xd - kd * xd * dt + sd * dBd
        xf = xf + (-kf * xf - rSf * sf * nu) * dt + sf * dBf
        xd = xd_n

    S = np.exp(lnS); D = np.exp(-intr)
    cd = xd + sd * sd / (2 * kd * kd) * (1 - np.exp(-kd * T)) ** 2
    cf = xf + sf * sf / (2 * kf * kf) * (1 - np.exp(-kf * T)) ** 2
    mine_d = np.exp(-kd * T) * S * L
    mine_f = np.exp(-kf * T) * S * Lt

    print(f"T={T}, paths={M}, steps/yr={spy}")
    for nm, y, m in (("Y^{d,S} vs e^{-kd t} S L^d", YdS, mine_d), ("Y^{f,S} vs e^{-kf t} S L^f", YfS, mine_f)):
        rel = np.sqrt(np.mean((y - m) ** 2)) / np.sqrt(np.mean(m ** 2))
        print(f"  {nm}:  corr = {np.corrcoef(y, m)[0, 1]:.6f},  relative RMS difference = {rel:.2e}")

    def pm(x):
        h = len(x) // 2; p = 0.5 * (x[:h] + x[h:]); return p.mean(), p.std(ddof=1) / np.sqrt(h)

    print("  IBP identities using the Y's directly (z = difference / s.e.)")
    for K in (0.7, 1.0, 1.3):
        ind = S > K; pay = np.maximum(S - K, 0)
        a1, b1 = D * cd * pay, D * ind * YdS
        a2, b2 = D * S * cf * pay, D * S * ind * YfS
        d1, s1 = pm(a1 - b1); d2, s2 = pm(a2 - b2)
        print(f"   K={K:.2f}  domestic: {pm(a1)[0]:+.4e} vs {pm(b1)[0]:+.4e} (z={d1 / s1:+.2f})"
              f"   foreign: {pm(a2)[0]:+.4e} vs {pm(b2)[0]:+.4e} (z={d2 / s2:+.2f})")


if __name__ == "__main__":
    run(spy=int(sys.argv[1]) if len(sys.argv) > 1 else 100)
