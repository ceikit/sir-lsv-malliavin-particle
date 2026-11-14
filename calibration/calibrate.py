"""
Particle calibration of the three-factor FX local volatility model (local vol +
Hull-White USD and EUR rates) to the CME EUR/USD surface, then independent Monte
Carlo repricing of the listed options.

Domestic = USD, foreign = EUR, S = USD per EUR. Everything is simulated under the
domestic risk-neutral measure Q^d, with discount weights in the kernel regression.

Inputs : futures.csv, parity.csv, svi_fits.csv in the work folder (from build_surface.py)
Outputs: printed repricing table, calib_<label>.npz in the work folder

usage: python calibration/calibrate.py <label> [--rates base|stress|none] [--estimator malliavin|none]
           [--pairs N] [--steps-per-year N] [--seed N] [--smooth W] [--min-neff N] [--work DIR] [--diag]
"""
import argparse
import dataclasses
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.signal import savgol_filter
from scipy.stats import norm

import _repo  # noqa: F401  (import path)
from market import WORK_DIR, Market, black
from sirlsv.hull_white import convexity as conv, variance as var
from sirlsv.model import RateParams
from sirlsv.sensitivities import sensitivity_step

T_MAX = 5.0
T_SYN = [1.5, 2.5, 3.5, 5.0]          # repricing dates beyond the data: synthetic surface
VAR_FLOOR = 0.005 ** 2                # floor on the calibrated local variance
G = 81                                # points of the log-moneyness grid at each time step
SMILE_SD = 0.075                      # rough EUR/USD volatility: sets the width of the grid and of the kernel

# rate parameters (illustrative, not market-calibrated)
RATE_PARAMS = {
    "base":   RateParams(kd=0.03, sd=0.010, kf=0.03, sf=0.008, rSd=-0.30, rSf=0.12, rdf=0.63),
    "stress": RateParams(kd=0.03, sd=0.015, kf=0.03, sf=0.012, rSd=0.30, rSf=-0.30, rdf=0.40),
    "none":   RateParams(kd=0.03, sd=1e-9, kf=0.03, sf=1e-9, rSd=0.0, rSf=0.0, rdf=0.0),
}

quart = lambda u: np.where(np.abs(u) < 1, (1 - u * u) ** 2, 0.0)


def time_grid(mkt, spy):
    """uniform grid with spy steps per year, plus every repricing date"""
    exps = sorted(set(mkt.Ts) | set(T_SYN))
    grid_t = np.unique(np.concatenate([np.linspace(0, T_MAX, int(T_MAX * spy) + 1), exps]))
    return grid_t[grid_t <= T_MAX + 1e-12]


def xgrid(t):
    """grid of x = ln(S / F(t)) at time t: +-5 standard deviations"""
    s = SMILE_SD * np.sqrt(max(t, 1e-3))
    return np.linspace(-5 * s, 5 * s, G)


def conditional_mean(x, y, D, xg, h, min_neff):
    """Kernel estimate of E[D y | x = xg] / E[D | x = xg] on the grid xg, bandwidth h.

    Trusted only where enough particles contribute, on one contiguous region around the
    money; flat beyond it."""
    order = np.argsort(x); xs = x[order]
    num = np.zeros(G); den = np.zeros(G); den2 = np.zeros(G)
    lo = np.searchsorted(xs, xg - h); hi = np.searchsorted(xs, xg + h)
    for g in range(G):
        if hi[g] > lo[g]:
            idx = order[lo[g]:hi[g]]
            kw = D[idx] * quart((x[idx] - xg[g]) / h)
            num[g] = np.sum(kw * y[idx]); den[g] = np.sum(kw); den2[g] = np.sum(kw * kw)
    neff = np.where(den2 > 0, den ** 2 / np.maximum(den2, 1e-300), 0.0)
    ok = (den > 0) & (neff >= max(min_neff, 1))
    ok[: np.argmax(ok)] = False; ok[G - np.argmax(ok[::-1]):] = False
    return np.interp(xg, xg[ok], num[ok] / den[ok])


def simulate(mkt, p, grid_t, n_pairs, seed, calibrate, estimator="malliavin", sig_tab=None, price_at=None,
             smooth=0, min_neff=0, diag=False):
    """Simulate under Q^d with n_pairs antithetic pairs.

    If calibrate, build sigma(t, x) step by step (with the Malliavin correction if
    estimator == "malliavin", Dupire only otherwise); if not, use the table sig_tab.
    Returns the table of sigma on xgrid(t) for each step, the log of the steps where the
    variance floor bound (time, distance to the money in standard deviations, number of
    grid points), and (S, D) at the dates in price_at.
    """
    kd, sd, kf, sf, rSf = p.kd, p.sd, p.kf, p.sf, p.rSf
    lnF, lnPd = mkt.lnF, mkt.lnPd
    chol = np.linalg.cholesky(p.correlation())
    rng = np.random.default_rng(seed)
    M = 2 * n_pairs
    lnS = np.full(M, np.log(mkt.S0)); xd = np.zeros(M); xf = np.zeros(M); intr = np.zeros(M)
    ld = np.zeros(M); zf = np.zeros(M); lf = np.zeros(M); zd = np.zeros(M)
    tab, floor_log, prices = [], [], {}
    g_last = np.zeros(1)
    sig_x = xgrid(0.0); sig_v = np.sqrt(mkt.dupire_variance(sig_x, 1e-3))
    for j in range(len(grid_t) - 1):
        t, t1 = grid_t[j], grid_t[j + 1]; dt = t1 - t
        lnFt = lnF(t)
        if calibrate:
            xg = xgrid(t)
            dup = mkt.dupire_variance(xg, max(t, 1e-3))
            if t > 0 and estimator == "malliavin":
                lam = ld - lf
                h = max(1.2 * SMILE_SD * np.sqrt(t) * (2 * n_pairs) ** (-0.2), 1e-4)
                corr = conditional_mean(lnS - lnFt, lam, np.exp(-intr), xg, h, min_neff)
                if smooth:                                                  # the correction is smooth in K: regularise it
                    corr = savgol_filter(corr, smooth, 2, mode='nearest')
                v = dup - 2 * corr
                if diag and (j % 25 == 0 or (t > 3.7 and j % 5 == 0)):
                    q = np.quantile(np.abs(lam), [0.5, 0.99, 0.9999])
                    print(f'    t={t:5.2f}  |lam| median {q[0]:.2e} p99 {q[1]:.2e} p99.99 {q[2]:.2e} max {np.abs(lam).max():.2e}'
                          f'  corr ATM {corr[G//2]:+.2e}  corr range [{corr.min():+.2e},{corr.max():+.2e}]  g range [{g_last.min():+.3f},{g_last.max():+.3f}]', flush=True)
            else:
                v = dup.copy()
            hit = v < VAR_FLOOR
            if hit.any():
                sd_x = SMILE_SD * np.sqrt(max(t, 1e-3))
                floor_log.append((t, float(np.abs(xg[hit] / sd_x).min()), int(hit.sum())))
            sig_x, sig_v = xg, np.sqrt(np.maximum(v, VAR_FLOOR))
            # light smoothing before differentiating (3-point), needed for g = dsigma/dlnS
            sig_v = np.convolve(np.pad(sig_v, 1, mode="edge"), [0.25, 0.5, 0.25], mode="valid")
            tab.append(sig_v.copy())
        else:
            sig_x, sig_v = xgrid(t), sig_tab[j]
        x = lnS - lnFt
        sig = np.interp(x, sig_x, sig_v)
        dsx = np.gradient(sig_v, sig_x)
        g_last = dsx
        g = np.where((x > sig_x[0]) & (x < sig_x[-1]), np.interp(x, sig_x, dsx), 0.0)

        Z = rng.standard_normal((n_pairs, 3)) @ chol.T
        Z = np.concatenate([Z, -Z]) * np.sqrt(dt)
        dW, dBd, dBf = Z[:, 0], Z[:, 1], Z[:, 2]
        if calibrate and estimator == "malliavin":
            ld, zf, lf, zd = sensitivity_step(ld, zf, lf, zd, g, sig, dW, dt, var(t, kd, sd), var(t, kf, sf), p)
        icd = 0.5 * dt * (conv(t, kd, sd) + conv(t1, kd, sd)); icf = 0.5 * dt * (conv(t, kf, sf) + conv(t1, kf, sf))
        intr += (lnPd(t) - lnPd(t1)) + icd + xd * dt
        lnS += (lnF(t1) - lnF(t)) + icd - icf + (xd - xf) * dt - 0.5 * sig * sig * dt + sig * dW
        xd_n = xd - kd * xd * dt + sd * dBd
        xf = xf + (-kf * xf - rSf * sf * sig) * dt + sf * dBf
        xd = xd_n
        if price_at is not None:
            for Te in price_at:
                if abs(t1 - Te) < 1e-12:
                    prices[Te] = (np.exp(lnS), np.exp(-intr))
    return np.array(tab), floor_log, prices


def reprice(mkt, prices, expiries):
    """Repricing error against the market smile, at 10- and 25-delta calls and puts and at the money.

    Control variates on E[D S] and E[D], both exact in the model. err and se are in vol points."""
    wtot = mkt.total_variance
    rows = []
    for Te in expiries:
        S, D = prices[Te]
        Fm, Pd = np.exp(mkt.lnF(Te)), np.exp(mkt.lnPd(Te))
        h = len(S) // 2
        for delta in (0.10, 0.25, 0.50, 0.75, 0.90):        # call deltas: 10d call ... 10d put
            # strike for this call delta on the market smile
            f = lambda k: norm.cdf((-k + 0.5 * wtot(k, Te)) / np.sqrt(wtot(k, Te))) - delta
            k = brentq(f, -0.6, 0.6); K = Fm * np.exp(k); vol = np.sqrt(wtot(k, Te) / Te)
            cp = "C" if k >= 0 else "P"
            pay = D * (np.maximum(S - K, 0) if cp == "C" else np.maximum(K - S, 0))
            X = np.vstack([D * S - Pd * Fm, D - Pd]).T
            Xc = X - X.mean(0); keep = Xc.std(0) > 1e-12 * max(Xc.std(0).max(), 1e-300)
            keep &= Xc.std(0) > 1e-14                        # with deterministic rates D is constant: drop that column
            beta = np.zeros(2)
            if keep.any():
                beta[keep] = np.linalg.lstsq(Xc[:, keep], pay - pay.mean(), rcond=1e-10)[0]
            adj = pay - X @ beta
            pp = 0.5 * (adj[:h] + adj[h:])
            price, se = pp.mean(), pp.std(ddof=1) / np.sqrt(h)
            vega = Pd * Fm * norm.pdf((-k + 0.5 * vol * vol * Te) / (vol * np.sqrt(Te))) * np.sqrt(Te)
            mkt_price = black(Fm, K, Te, vol, Pd, cp)
            rows.append(dict(T=Te, delta=delta, K=K, vol=100 * vol, err=100 * (price - mkt_price) / vega, se=100 * se / vega))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description="Particle calibration to the CME EUR/USD surface and repricing.")
    ap.add_argument("label", help="name of the run; results go to calib_<label>.npz")
    ap.add_argument("--rates", choices=sorted(RATE_PARAMS), default="base", help="rate parameters (none = deterministic rates)")
    ap.add_argument("--estimator", choices=["malliavin", "none"], default="malliavin", help="none = Dupire local vol, no rates correction")
    ap.add_argument("--pairs", type=int, default=100_000, help="antithetic pairs")
    ap.add_argument("--steps-per-year", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--smooth", type=int, default=0, help="Savitzky-Golay window (odd) on the correction, 0 = off")
    ap.add_argument("--min-neff", type=int, default=0, help="minimum effective particle count for the kernel estimate")
    ap.add_argument("--work", default=str(WORK_DIR), help="folder with the files from build_surface.py (default: calibration/work)")
    ap.add_argument("--diag", action="store_true", help="print traces of the sensitivities during calibration")
    a = ap.parse_args()

    t_start = time.time()
    mkt = Market(a.work)
    p = RATE_PARAMS[a.rates]
    grid_t = time_grid(mkt, a.steps_per_year)
    expiries = list(mkt.Ts) + T_SYN
    print(f"[{a.label}] rates={a.rates} {dataclasses.asdict(p)}  estimator={a.estimator}  particles={2*a.pairs}  "
          f"steps/yr={a.steps_per_year}  smoothing window={a.smooth}  min n_eff={a.min_neff}")
    print(f"  S0={mkt.S0:.5f}  USD zero 1y={-mkt.lnPd(1.0)*100:.2f}%  f_d-f_f at 1y={mkt.carry(1.0)*100:.2f}%  "
          f"calendar-repaired slices={mkt.n_repaired}")
    tab, floor_log, _ = simulate(mkt, p, grid_t, a.pairs, a.seed, True, a.estimator,
                                 smooth=a.smooth, min_neff=a.min_neff, diag=a.diag)
    t_cal = time.time() - t_start
    # independent repricing with the sigma table frozen
    _, _, pr = simulate(mkt, p, grid_t, a.pairs, a.seed + 1000, False, a.estimator, sig_tab=tab, price_at=expiries)
    t_all = time.time() - t_start
    R = reprice(mkt, pr, expiries)

    piv = R.pivot_table(index="T", columns="delta", values="err").round(3)
    piv.columns = ["10dC", "25dC", "ATM", "25dP", "10dP"]
    piv.index = [f"{t:.3f}" + (" (synthetic)" if t > mkt.T_last + 1e-9 else "") for t in piv.index]
    piv["max se"] = R.groupby("T").se.max().round(3).values
    print("\n  repricing error, model minus market, in vol points (Monte Carlo standard error in last column)")
    print(piv.to_string())
    Rm = R[R["T"] <= mkt.T_last + 1e-9]
    floors = sum(n for _, _, n in floor_log)
    print(f"\n  market expiries only: mean |err| {Rm.err.abs().mean():.3f}  max |err| {Rm.err.abs().max():.3f}  mean |err|/se {(Rm.err.abs()/Rm.se).mean():.2f}")
    print(f"  all rows: mean |err| {R.err.abs().mean():.3f}  max |err| {R.err.abs().max():.3f}  "
          f"mean |err|/se {(R.err.abs()/R.se).mean():.2f}  floor hits {floors}  "
          f"time: calibration {t_cal:.0f}s, total {t_all:.0f}s")
    if floor_log:
        FL = np.array(floor_log)
        inner = FL[FL[:, 1] < 3.0]
        print(f"  floor binds at {len(FL)} of {len(grid_t)-1} time steps; closest to the money: "
              f"{FL[:,1].min():.1f} std devs; steps where it binds inside 3 std devs: {len(inner)}"
              + (f" (t from {inner[:,0].min():.2f} to {inner[:,0].max():.2f})" if len(inner) else ""))
        atm1 = FL[FL[:, 1] < 1.0]
        if len(atm1):
            print(f"  floor binds within 1 std dev of the money from t = {atm1[:,0].min():.2f}")
    np.savez(Path(a.work) / f"calib_{a.label}.npz", grid_t=grid_t, tab=tab, R=R.to_records(index=False),
             t_last=mkt.T_last, floor_log=np.array(floor_log).reshape(-1, 3))


if __name__ == "__main__":
    main()
