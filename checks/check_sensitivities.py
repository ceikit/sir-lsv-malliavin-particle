"""
Monte Carlo check of the two identities behind the main formula, in the three-factor
FX model with a prescribed skewed local volatility (see sirlsv/prescribed_vol.py for
the model and sirlsv/sensitivities.py for the auxiliary processes l_d, z_f, l_f, z_d).

Identities tested (integrated in strike; c^i = r^i_T - f^i(0,T), D = exp(-int r^d)):
    (I)   E[ D c^d (S_T-K)^+ ]     = E[ D 1_{S_T>K} S_T   l_d ]
    (II)  E[ D S_T c^f (S_T-K)^+ ] = E[ D 1_{S_T>K} S_T^2 l_f ]

'naive' = same recursions with the coupling terms (g l_d in z_f, g l_f in z_d)
dropped and with g (dW - sigma dt) instead of g dW in l_f.

usage: python checks/check_sensitivities.py [pairs] [steps_per_year] [T] [skew] [rate_vol_scale] [seed]
"""
import sys
import time

import numpy as np

import _repo  # noqa: F401  (import path)
from sirlsv.montecarlo import antithetic_mean as pm
from sirlsv.prescribed_vol import simulate

STRIKES = [0.70, 0.85, 1.00, 1.15, 1.30, 1.50]


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


def main(argv):
    arg = lambda i, d, f: f(argv[i]) if len(argv) > i else d
    N, spy, T, skew, vs = arg(1, 400_000, int), arg(2, 100, int), arg(3, 5.0, float), arg(4, 0.6, float), arg(5, 1.0, float)
    t0 = time.time()
    res = simulate(N, spy, T, skew, vs, seed=arg(6, 7, int))
    print(f"T={T}, paths={2 * N}, steps/yr={spy}, skew={skew}, rate-vol scale={vs}, {time.time() - t0:.0f}s")
    report(res, STRIKES)


if __name__ == "__main__":
    main(sys.argv)
