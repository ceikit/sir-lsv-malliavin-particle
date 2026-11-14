"""Implied-vol surface from CME EUR/USD option settlements (one trade date).

Black-76 on the underlying future, discount factor from put-call parity per expiry.
Uses out-of-the-money options only, with premium >= MIN_TICKS ticks. One SVI smile
per expiry.

Writes to the work folder (default calibration/work, git-ignored):
    implied_vols.csv    one row per option used
    svi_fits.csv        SVI parameters, forward and discount factor per expiry
    parity.csv          put-call parity: discount factor, implied forward, worst residual
    futures.csv         futures settlements with their time to expiry
    surface_quotes.csv  ATM, risk reversals and butterflies per expiry (also printed)

usage: python calibration/build_surface.py <settlements.csv> [--work DIR] [--valuation-time "YYYY-MM-DD HH:MM"]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares
from scipy.stats import norm

from market import WORK_DIR, black

TICK, MIN_TICKS = 5e-5, 3
SETTLEMENT_HOUR_UTC = "21:00"       # valuation time on the trade date, unless given


def iv(p, F, K, T, D, cp):
    intr = D * max(F - K, 0) if cp == "C" else D * max(K - F, 0)
    if p <= intr + 1e-12: return np.nan
    try: return brentq(lambda s: black(F, K, T, s, D, cp) - p, 1e-4, 2.0)
    except ValueError: return np.nan


def svi(k, a, b, rho, m, sig):
    """SVI raw parametrisation of total variance in log-moneyness k"""
    return a + b * (rho * (k - m) + np.sqrt((k - m) ** 2 + sig ** 2))


def strike_from_delta(fit, delta, cp):           # forward delta, premium-unadjusted
    T, F = fit["T"], fit["F"]; par = [fit[c] for c in "a b rho m sig".split()]
    f = lambda k: (norm.cdf((-k + 0.5 * svi(k, *par)) / np.sqrt(svi(k, *par))) - delta) if cp == "C" else \
                  (norm.cdf(-(-k + 0.5 * svi(k, *par)) / np.sqrt(svi(k, *par))) - delta)
    k = brentq(f, -1.5, 1.5); return k, np.sqrt(svi(k, *par) / T)


def gfun(k, a, b, rho, m, sig):
    """Gatheral's g(k); no butterfly arbitrage on the slice if g >= 0"""
    w = svi(k, a, b, rho, m, sig); x = k - m; r = np.sqrt(x*x+sig*sig)
    w1 = b*(rho + x/r); w2 = b*sig*sig/r**3
    return (1 - k*w1/(2*w))**2 - w1**2/4*(1/w + 0.25) + w2/2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("settlements", help="CSV written by data/pull_eurusd_options_databento.py")
    ap.add_argument("--work", default=str(WORK_DIR), help="output folder (default: calibration/work)")
    ap.add_argument("--valuation-time", default=None, help='UTC, e.g. "2026-11-12 21:00" (default: trade date at 21:00 UTC)')
    args = ap.parse_args()
    work = Path(args.work); work.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.settlements, parse_dates=["expiration"])
    if args.valuation_time is None:
        dates = df.trade_date.unique()
        assert len(dates) == 1, f"expected one trade date, found {list(dates)}"
        args.valuation_time = f"{dates[0]} {SETTLEMENT_HOUR_UTC}"
    t0 = pd.Timestamp(args.valuation_time, tz="UTC")
    years = lambda expiration: (expiration - t0).dt.total_seconds() / (365.25 * 86400)

    futures = df[df.instrument_class == "F"].copy()
    futures["T"] = years(futures.expiration)
    futures[["raw_symbol", "expiration", "T", "settlement"]].to_csv(work / "futures.csv", index=False)
    fut = futures.set_index("raw_symbol")["settlement"]
    opt = df[df.instrument_class.isin(["C", "P"]) & (df.settlement > 0)].copy()
    opt["T"] = years(opt.expiration)
    opt["F"] = opt.underlying.map(fut)

    pts, fits, parity = [], [], []
    for exp, g in opt.groupby("expiration"):
        F, T = g.F.iloc[0], g["T"].iloc[0]
        c = g[g.instrument_class == "C"].set_index("strike_price").settlement
        p = g[g.instrument_class == "P"].set_index("strike_price").settlement
        ks = sorted(set(c.index) & set(p.index)); k_ = np.array(ks)
        # put-call parity near the money: C - P = D (F - K), fitted as a_ - D K
        near = [x for x in ks if abs(x - F) < 0.06 and c[x] > 10 * TICK and p[x] > 10 * TICK]
        A = np.vstack([np.ones(len(near)), -np.array(near)]).T
        cmp_ = np.array([c[x] - p[x] for x in near])
        (a_, D), *_ = np.linalg.lstsq(A, cmp_, rcond=None)
        parity.append(dict(expiry=exp, und=g.underlying.iloc[0], T=T, F_fut=F, F_parity=a_ / D, D=D, n_pairs=len(near),
                           max_resid_ticks=np.abs(cmp_ - (a_ - D * np.array(near))).max() / TICK))
        for _, r in g.iterrows():
            K, cp = r.strike_price, r.instrument_class
            if (cp == "C" and K < F) or (cp == "P" and K >= F): continue        # OTM only
            if r.settlement < MIN_TICKS * TICK: continue
            v = iv(r.settlement, F, K, T, D, cp)
            if np.isfinite(v):
                pts.append(dict(expiry=exp, T=T, F=F, D=D, K=K, cp=cp, price=r.settlement, iv=v,
                                k=np.log(K / F), oi=r.open_interest, vega_tick=TICK / (D * F * norm.pdf(np.log(F/K)/(v*np.sqrt(T)) + 0.5*v*np.sqrt(T)) * np.sqrt(T))))
        P = pd.DataFrame([q for q in pts if q["expiry"] == exp])
        w = P.iv.values ** 2 * T; kk = P.k.values
        wt = 1 / np.maximum(P.vega_tick.values * 2 * P.iv.values * T, 1e-7)     # weight by 1/(tick error in total variance)
        atm = np.interp(0, kk[np.argsort(kk)], w[np.argsort(kk)])
        res = least_squares(lambda x: (svi(kk, *x) - w) * wt / wt.max(), [atm * 0.8, 0.05, -0.1, 0.0, 0.05],
                            bounds=([-1, 1e-6, -0.999, -1, 1e-4], [1, 5, 0.999, 1, 2]))
        fits.append(dict(expiry=exp, T=T, F=F, D=D, a=res.x[0], b=res.x[1], rho=res.x[2], m=res.x[3], sig=res.x[4],
                         n=len(P), kmin=kk.min(), kmax=kk.max(),
                         rmse_volpts=100 * np.sqrt(np.mean((np.sqrt(svi(kk, *res.x) / T) - P.iv.values) ** 2))))
    pts = pd.DataFrame(pts); fits = pd.DataFrame(fits); parity = pd.DataFrame(parity)
    pts.to_csv(work / "implied_vols.csv", index=False); fits.to_csv(work / "svi_fits.csv", index=False)
    parity.to_csv(work / "parity.csv", index=False)

    out = []
    for _, f in fits.iterrows():
        par = [f[c] for c in "a b rho m sig".split()]
        atm = np.sqrt(svi(0.0, *par) / f["T"])
        k25c, v25c = strike_from_delta(f, 0.25, "C"); k25p, v25p = strike_from_delta(f, 0.25, "P")
        k10c, v10c = strike_from_delta(f, 0.10, "C"); k10p, v10p = strike_from_delta(f, 0.10, "P")
        out.append(dict(expiry=f.expiry.date(), T=round(f["T"], 3), n=int(f.n), ATM=100*atm,
                        RR25=100*(v25c-v25p), BF25=100*((v25c+v25p)/2-atm), RR10=100*(v10c-v10p), BF10=100*((v10c+v10p)/2-atm),
                        K10P=np.exp(k10p)*f.F, K10C=np.exp(k10c)*f.F, Kmin=np.exp(f.kmin)*f.F, Kmax=np.exp(f.kmax)*f.F,
                        rmse=f.rmse_volpts))
    q = pd.DataFrame(out).round(3); q.to_csv(work / "surface_quotes.csv", index=False)
    print(q.to_string(index=False))
    # calendar arbitrage check: total variance non-decreasing in T at fixed k
    grid = np.linspace(-0.15, 0.15, 61)
    W = np.array([[svi(k, *[f[c] for c in "a b rho m sig".split()]) for k in grid] for _, f in fits.iterrows()])
    dW = np.diff(W, axis=0)
    print("calendar violations (pairs of consecutive expiries, any k in +-15%):", int((dW < -1e-6).any(axis=1).sum()), "of", len(dW))
    # butterfly check on each slice (Gatheral g(k) >= 0)
    print("butterfly violations (slices with g<0 in +-15%):", sum((gfun(grid, *[f[c] for c in "a b rho m sig".split()]) < 0).any() for _, f in fits.iterrows()))
    print(f"put-call parity: worst residual {parity.max_resid_ticks.max():.1f} ticks; "
          f"worst |implied forward - futures settlement| {np.abs(parity.F_parity - parity.F_fut).max():.1e}")


if __name__ == "__main__":
    main()
