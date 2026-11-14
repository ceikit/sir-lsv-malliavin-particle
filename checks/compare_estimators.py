"""
Fixed-local-volatility comparison of two particle estimators of the stochastic-rates
correction to the local variance, strike by strike:

  A (conventional, single measure, OTM side):
        E^d[ D (K c^d - S c^f) 1_{S>K} ] / ( K^2 C_KK / 2 )        (C_KK taken as known)
  B (integrated Malliavin sensitivities, kernel regression):
        2 * sum_i D_i (l_d - l_f)_i k_h(S_i-K) / sum_i D_i k_h(S_i-K)

Both estimate the same number (the single-measure calibration condition against the main
formula of the paper).  The script reports, over independent batches of particles, the
mean and the standard deviation of each.
This is NOT a calibration benchmark: sigma is fixed and no control variate is used.

usage: python checks/compare_estimators.py [pairs] [steps_per_year] [batch_particles]
"""
import sys

import numpy as np

import _repo  # noqa: F401  (import path)
from sirlsv.model import tanh_local_vol
from sirlsv.prescribed_vol import FD, FF, S0, VOL_B, VOL_LEVEL, simulate

T = 5.0
SKEW = 0.6
QUANTILES = [0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99]
quart = lambda u: np.where(np.abs(u) < 1, 15.0 / 16.0 * (1 - u * u) ** 2, 0.0)


def main(argv):
    N = int(argv[1]) if len(argv) > 1 else 1_000_000
    spy = int(argv[2]) if len(argv) > 2 else 100
    nb = int(argv[3]) if len(argv) > 3 else 20_000
    r = simulate(N, spy, T, SKEW, 1.0, seed=23)
    S, D, cd, cf = r["S"], r["D"], r["cd"], r["cf"]
    lam = r["ld"] - r["lf"]
    fwd = S0 * np.exp((FD - FF) * T)

    qs = QUANTILES
    Ks = np.quantile(S, qs)
    sig = lambda K: tanh_local_vol(np.log(K), VOL_LEVEL, SKEW, VOL_B)[0]

    # reference quantities from the full sample
    h_ref = 0.015
    CKK = np.array([np.mean(D * quart((S - K) / h_ref)) / h_ref for K in Ks])
    B_ref = np.array([2 * np.sum(D * lam * quart((S - K) / h_ref)) / np.sum(D * quart((S - K) / h_ref)) for K in Ks])

    def est_A(idx, K, ckk):
        s, d = S[idx], D[idx]
        w = d * (K * cd[idx] - s * cf[idx])
        num = np.mean(w * (s > K)) if K >= fwd else -np.mean(w * (s < K))
        return num / (0.5 * K * K * ckk)

    def est_B(idx, K, h):
        s, d = S[idx], D[idx]
        k = quart((s - K) / h)
        return 2 * np.sum(d * lam[idx] * k) / np.sum(d * k)

    half = nb // 2
    nbatch = N // half
    h = 1.5 * 0.12 * np.sqrt(T) * nb ** (-0.2)
    A = np.zeros((nbatch, len(Ks)))
    B = np.zeros((nbatch, len(Ks)))
    for b in range(nbatch):
        i = np.arange(b * half, (b + 1) * half)
        idx = np.concatenate([i, i + N])              # antithetic partners stay in the same batch
        for j, K in enumerate(Ks):
            A[b, j] = est_A(idx, K, CKK[j])
            B[b, j] = est_B(idx, K, h)

    print(f"T={T}, {2 * N} paths, {nbatch} batches of {nb} particles, steps/yr={spy}, bandwidth h={h:.3f}")
    print("correction in local-variance units x 1e4; 'vol pts' = std of the implied local vol, in vol points")
    print("  quantile     K  sigma(K) | reference |   A mean    A std  vol pts |   B mean    B std  vol pts | std A / std B")
    for j, K in enumerate(Ks):
        s = sig(K)
        a_m, a_s, b_m, b_s = A[:, j].mean(), A[:, j].std(ddof=1), B[:, j].mean(), B[:, j].std(ddof=1)
        print(f"   {qs[j]:5.2f}  {K:6.3f}  {s:6.3f}  | {B_ref[j] * 1e4:8.3f}  | {a_m * 1e4:8.3f} {a_s * 1e4:8.3f} {a_s / (2 * s) * 100:7.3f}  |"
              f" {b_m * 1e4:8.3f} {b_s * 1e4:8.3f} {b_s / (2 * s) * 100:7.3f}  | {a_s / b_s:8.1f}")
    # Strike-by-strike check of the identity itself, free of any density estimate:
    #     E[ D (K c^d - S c^f) 1_{S>K} ]  =  K^2 E[ D (l_d - l_f) delta(S-K) ]
    # (the right hand side with a narrow kernel; standard error from antithetic pairs).
    h0 = 0.01
    print("identity at each strike (both sides x 1e4), z = difference / standard error")
    for j, K in enumerate(Ks):
        w = D * (K * cd - S * cf)
        lhs = w * (S > K) if K >= fwd else -w * (S < K)
        rhs = K * K * D * lam * quart((S - K) / h0) / h0
        d = lhs - rhs
        p = 0.5 * (d[:N] + d[N:])
        print(f"   q={qs[j]:4.2f}  K={K:6.3f}   {lhs.mean() * 1e4:9.4f}  {rhs.mean() * 1e4:9.4f}   z={p.mean() / (p.std(ddof=1) / np.sqrt(N)):+5.2f}")


if __name__ == "__main__":
    main(sys.argv)
