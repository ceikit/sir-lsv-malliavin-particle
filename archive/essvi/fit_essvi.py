"""
Arbitrage-free fit of the implied-vol surface: eSSVI (Hendriks-Martini), slice by slice.

  w(k) = theta/2 * (1 + rho*psi/theta*k + sqrt((psi/theta*k + rho)^2 + 1 - rho^2))

Constraints, from one slice to the next (no calendar arbitrage):
  theta_i >= theta_{i-1},  psi_i - psi_{i-1} >= |rho_i psi_i - rho_{i-1} psi_{i-1}|
and on each slice (no butterfly arbitrage, Gatheral-Jacquier):
  psi (1 + |rho|) <= 4,   psi^2 (1 + |rho|) <= 4 theta
Input: implied_vols.csv (build_surface.py). Output: essvi_fits.csv and a fit report.

Archived: this fit was tried for the calibration test and rejected (it misfit the 2027 risk
reversals by up to 1.6 vol points). Run it from the folder that holds implied_vols.csv
(calibration/work by default).
"""
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
from scipy.optimize import brentq

P = pd.read_csv("implied_vols.csv", parse_dates=["expiry"])
def w_essvi(k, th, rho, psi):
    phi = psi / th
    return 0.5 * th * (1 + rho * phi * k + np.sqrt((phi * k + rho) ** 2 + 1 - rho * rho))

fits, prev = [], None
for exp, g in P.groupby("expiry"):
    T = g["T"].iloc[0]; k = g.k.values; vol = g.iv.values
    wt = 1 / np.maximum(g.vega_tick.values, 1e-4)                       # error in vol of one tick
    wt = wt / wt.max()
    obj = lambda x: np.sum((wt * (np.sqrt(w_essvi(k, *x) / T) - vol) * 1e4) ** 2)
    th0 = np.interp(0, np.sort(k), (vol ** 2 * T)[np.argsort(k)])
    cons = []
    for sgn in (1, -1):                                                # |rho| written as two smooth constraints
        cons += [{"type": "ineq", "fun": lambda x, s=sgn: 4 - x[2] * (1 + s * x[1])},
                 {"type": "ineq", "fun": lambda x, s=sgn: 4 * x[0] - x[2] ** 2 * (1 + s * x[1])}]
    if prev is not None:
        tp, rp, pp = prev
        cons += [{"type": "ineq", "fun": lambda x, tp=tp: x[0] - tp - 1e-7}]
        for sgn in (1, -1):
            cons += [{"type": "ineq", "fun": lambda x, rp=rp, pp=pp, s=sgn: (x[2] - pp) - s * (x[1] * x[2] - rp * pp) - 1e-9}]
    feasible = lambda x: all(c["fun"](x) >= -1e-9 for c in cons)
    best = None
    for r0 in (-0.3, -0.1, 0.0):
        for p0 in (0.5 * np.sqrt(th0), np.sqrt(th0), 2 * np.sqrt(th0)):
            x0 = [max(th0, prev[0] + 1e-6) if prev else th0, r0, max(p0, prev[2] + 1e-4) if prev else p0]
            r = minimize(obj, x0, method="SLSQP", constraints=cons,
                         bounds=[(1e-6, 1), (-0.99, 0.99), (1e-5, 4)], options=dict(maxiter=500, ftol=1e-12))
            if feasible(r.x) and (best is None or r.fun < best.fun):
                best = r
    if best is None:
        raise SystemExit(f"no feasible fit for expiry {exp.date()}")
    th, rho, psi = best.x
    fit_vol = np.sqrt(w_essvi(k, th, rho, psi) / T)
    err = 100 * (fit_vol - vol)
    fits.append(dict(expiry=exp, T=T, theta=th, rho=rho, psi=psi, n=len(k),
                     rmse=np.sqrt(np.mean(err ** 2)), maxerr=np.abs(err).max(),
                     wrmse=np.sqrt(np.average(err ** 2, weights=wt ** 2)), tick_vol_atm=100 * g.vega_tick.min()))
    prev = (th, rho, psi)
F = pd.DataFrame(fits); F.to_csv("essvi_fits.csv", index=False)

# quotes from the fit, and the same from the per-slice SVI fits for comparison
svi = pd.read_csv("surface_quotes.csv")
def quotes(th, rho, psi, T):
    w = lambda k: w_essvi(k, th, rho, psi)
    def kd(delta, cp):
        f = (lambda k: norm.cdf((-k + 0.5 * w(k)) / np.sqrt(w(k))) - delta) if cp == "C" else \
            (lambda k: norm.cdf((k - 0.5 * w(k)) / np.sqrt(w(k))) - delta)
        k = brentq(f, -1.5, 1.5); return np.sqrt(w(k) / T)
    atm = np.sqrt(th / T); c25, p25, c10, p10 = kd(.25, "C"), kd(.25, "P"), kd(.10, "C"), kd(.10, "P")
    return 100 * atm, 100 * (c25 - p25), 100 * ((c25 + p25) / 2 - atm), 100 * (c10 - p10), 100 * ((c10 + p10) / 2 - atm)
rows = []
for i, f in F.iterrows():
    a, rr, bf, rr10, bf10 = quotes(f.theta, f.rho, f.psi, f["T"])
    s = svi.iloc[i]
    rows.append(dict(expiry=str(f.expiry.date()), T=round(f["T"], 3), ATM=a, dATM=a - s.ATM, RR25=rr, dRR25=rr - s.RR25,
                     BF25=bf, dBF25=bf - s.BF25, RR10=rr10, BF10=bf10, rmse=f.rmse, maxerr=f.maxerr, tick=f.tick_vol_atm))
Q = pd.DataFrame(rows).round(3); Q.to_csv("essvi_quotes.csv", index=False)
print("eSSVI fit (vol points). d* = difference from the unconstrained per-expiry SVI fit; tick = one tick in vol at the money")
print(Q.to_string(index=False))
K = np.linspace(-0.6, 0.6, 1201)
W = np.array([w_essvi(K, f.theta, f.rho, f.psi) for _, f in F.iterrows()])
print("calendar violations on k in [-0.6, 0.6]:", int((np.diff(W, axis=0) < -1e-12).any(axis=1).sum()))
