"""Tables of calibration/RESULTS.md from the five runs of calibrate.py:

    det            --rates none   --estimator none        deterministic rates
    base           --rates base                           base rates, with the correction
    base_nocorr    --rates base   --estimator none        base rates, Dupire local vol
    stress         --rates stress
    stress_nocorr  --rates stress --estimator none

All runs must use the same particles, steps and seed: the repricing runs then share
their random numbers, and the error due to the stochastic rates is the difference
with the deterministic-rates run.

usage: python calibration/summarise.py [--work DIR]
"""
import argparse
from pathlib import Path

import numpy as np

import _repo  # noqa: F401  (import path)
from calibrate import G, RATE_PARAMS
from market import WORK_DIR, Market
from sirlsv.hull_white import variance

RUNS = [("deterministic rates", "det"), ("base, no correction", "base_nocorr"), ("base, with correction", "base"),
        ("stress, no correction", "stress_nocorr"), ("stress, with correction", "stress")]


def flat_vol_correction(p, grid_t, a):
    """u_d - u_f for a spot volatility a(t) that does not depend on the spot (check (b) of the paper):
    the sensitivities are then deterministic and solve two linear ODEs."""
    ud = uf = z = 0.0
    out = np.zeros(len(grid_t) - 1)
    for j in range(len(grid_t) - 1):
        t, dt = grid_t[j], grid_t[j + 1] - grid_t[j]
        out[j] = ud - uf
        vd, vf = variance(t, p.kd, p.sd), variance(t, p.kf, p.sf)
        ud, uf, z = (ud + (-p.kd * ud + vd - z + p.rSd * p.sd * a[j]) * dt,
                     uf + (-p.kf * uf + z - vf + p.rSf * p.sf * a[j]) * dt,
                     z + (-(p.kd + p.kf) * z + p.rdf * p.sd * p.sf) * dt)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", default=str(WORK_DIR))
    work = Path(ap.parse_args().work)
    res = {lab: np.load(work / f"calib_{lab}.npz") for _, lab in RUNS}
    R = {lab: z["R"] for lab, z in res.items()}
    T = R["det"]["T"]
    t_last = float(res["det"]["t_last"]) if "t_last" in res["det"] else 1.1
    mkt_rows, mid, five = T <= t_last + 1e-9, (T > t_last + 1e-9) & (T < 4.0), T > 4.0
    f = lambda x: f"{x:.3f}"

    print("Repricing error (model minus market), vol points\n")
    print("| | " + " | ".join(name for name, _ in RUNS) + " |")
    print("|---|" + "---|" * len(RUNS))
    for title, rows, agg in (("market expiries, mean", mkt_rows, np.mean), ("market expiries, max", mkt_rows, np.max),
                             ("5y (synthetic), max", five, np.max)):
        print(f"| {title} | " + " | ".join(f(agg(np.abs(R[lab]["err"][rows]))) for _, lab in RUNS) + " |")
    se = np.concatenate([R[lab]["se"] for _, lab in RUNS])
    print(f"\nMonte Carlo standard error of a single repricing error: {se.min():.2f} to {se.max():.2f}")
    i = np.argmax(np.abs(R["det"]["err"]))
    print(f"largest error, deterministic rates: {R['det']['err'][i]:+.3f} at T={T[i]:.3f}, call delta {R['det']['delta'][i]:.2f}")

    print("\nError due to the stochastic rates (difference with the deterministic-rates run), vol points\n")
    print("| | market expiries, mean | market expiries, max | 1.5y to 3.5y, max | 5y, max |")
    print("|---|---|---|---|---|")
    for name, lab in (("base, with correction", "base"), ("base, no correction", "base_nocorr"),
                      ("stress, with correction", "stress"), ("stress, no correction", "stress_nocorr")):
        d = np.abs(R[lab]["err"] - R["det"]["err"])
        print(f"| {name} | {f(d[mkt_rows].mean())} | {f(d[mkt_rows].max())} | {f(d[mid].max())} | {f(d[five].max())} |")

    print("\nCalibrated correction at the money against the flat-volatility closed form (variance units x 1e4)\n")
    print("| case | t | calibrated | flat-vol closed form | ratio |")
    print("|---|---|---|---|---|")
    mkt = Market(work)
    for lab in ("base", "stress"):
        grid_t, tab = res[lab]["grid_t"], res[lab]["tab"]
        atm = tab[:, G // 2]
        closed = flat_vol_correction(RATE_PARAMS[lab], grid_t, atm)
        for t in (1.0, 2.0, 3.0, 4.0, 5.0):
            j = min(int(np.searchsorted(grid_t, t)), len(atm) - 1)
            cal = 0.5 * (mkt.dupire_variance(0.0, grid_t[j]) - atm[j] ** 2)
            print(f"| {lab} | {grid_t[j]:.1f} | {cal * 1e4:.3f} | {closed[j] * 1e4:.3f} | {cal / closed[j]:.2f} |")

    print("\nVariance floor (grid points per run where the calibrated variance was floored)")
    for name, lab in RUNS:
        fl = res[lab]["floor_log"] if "floor_log" in res[lab] else None
        if fl is None:
            continue
        inner = fl[fl[:, 1] < 1.0] if len(fl) else fl
        print(f"  {name}: {int(fl[:, 2].sum()) if len(fl) else 0} points at {len(fl)} time steps"
              + (f"; within 1 std dev of the money from t = {inner[:, 0].min():.2f}" if len(inner) else ""))


if __name__ == "__main__":
    main()
