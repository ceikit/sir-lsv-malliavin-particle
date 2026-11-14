"""
Test 3: particle calibration of the three-factor FX local volatility model
(local vol + Hull-White USD and EUR rates) to the CME EUR/USD surface, then
independent Monte Carlo repricing of the listed options.

Domestic = USD, foreign = EUR, S = USD per EUR. Everything is simulated under the
domestic risk-neutral measure Q^d with discount weights (second line of (42) in the paper).

Inputs : svi_fits.csv, parity.csv (from build_surface.py) and the settlements CSV
Outputs: printed repricing table, calib_<label>.npz

usage: python calibrate.py <label> <rates: base|stress|none> <estimator: malliavin|none> [pairs] [steps/yr] [seed] [smooth window] [min effective particles]
"""
import sys, time, numpy as np, pandas as pd
from scipy.stats import norm
from scipy.optimize import brentq
from scipy.signal import savgol_filter

label = sys.argv[1] if len(sys.argv) > 1 else "base"
RATES = sys.argv[2] if len(sys.argv) > 2 else "base"
EST = sys.argv[3] if len(sys.argv) > 3 else "malliavin"
NP = int(sys.argv[4]) if len(sys.argv) > 4 else 100_000           # antithetic pairs
SPY = int(sys.argv[5]) if len(sys.argv) > 5 else 100
SEED = int(sys.argv[6]) if len(sys.argv) > 6 else 1
SMOOTH = int(sys.argv[7]) if len(sys.argv) > 7 else 0
MIN_NEFF = int(sys.argv[8]) if len(sys.argv) > 8 else 0                # minimum effective particle count for the kernel estimate                # Savitzky-Golay window (odd) on the correction, 0 = off
T_MAX = 5.0
VAR_FLOOR = 0.005 ** 2
FLOOR_LOG = []
import os
DIAG = bool(os.environ.get('DIAG'))

# ------------------------------------------------------------------ rate parameters (illustrative, not market-calibrated)
PARAMS = {
    "base":   dict(kd=0.03, sd=0.010, kf=0.03, sf=0.008, rSd=-0.30, rSf=0.12, rdf=0.63),
    "stress": dict(kd=0.03, sd=0.015, kf=0.03, sf=0.012, rSd=0.30, rSf=-0.30, rdf=0.40),
    "none":   dict(kd=0.03, sd=1e-9, kf=0.03, sf=1e-9, rSd=0.0, rSf=0.0, rdf=0.0),
}[RATES]
kd, sd, kf, sf = PARAMS["kd"], PARAMS["sd"], PARAMS["kf"], PARAMS["sf"]
rSd, rSf, rdf = PARAMS["rSd"], PARAMS["rSf"], PARAMS["rdf"]

# ------------------------------------------------------------------ market curves
raw = pd.read_csv("eurusd_cme_options_2026-11-12.csv", parse_dates=["expiration"])
t0 = pd.Timestamp("2026-11-12 21:00", tz="UTC")
fut = raw[raw.instrument_class == "F"].copy()
fut["T"] = (fut.expiration - t0).dt.total_seconds() / (365.25 * 86400)
cF = np.polyfit(fut["T"], np.log(fut.settlement), 2)                 # ln F(T) = c2 T^2 + c1 T + c0
lnF = lambda T: np.polyval(cF, T)
S0 = float(np.exp(cF[2]))
par = pd.read_csv("parity.csv")
cD = np.polyfit(par["T"], np.log(par["D"]), 2, w=np.sqrt(par["T"]))    # ln P^d(T), forced through 0 below
cD[2] = 0.0
lnPd = lambda T: np.polyval(cD, T)
fd0 = lambda t: -(2 * cD[0] * t + cD[1])                              # domestic instantaneous forward
mu0 = lambda t: 2 * cF[0] * t + cF[1]                                 # f_d - f_f
ff0 = lambda t: fd0(t) - mu0(t)

# ------------------------------------------------------------------ implied surface: SVI slices, extended to 5Y
svi = pd.read_csv("svi_fits.csv")
MAX_MKT_T = 1.1                                                       # Dec 2026 - Dec 2027; the 2028 expiries are dropped
svi = svi[svi["T"] <= MAX_MKT_T].reset_index(drop=True)                # (no open interest, calendar-inconsistent settlements)
Ts = svi["T"].values
SV = svi[["a", "b", "rho", "m", "sig"]].values
wslice = lambda i, k: SV[i, 0] + SV[i, 1] * (SV[i, 2] * (k - SV[i, 3]) + np.sqrt((k - SV[i, 3]) ** 2 + SV[i, 4] ** 2))
def wslices(k):
    """total variance of every market slice at log-moneyness k, with calendar arbitrage repaired by a running max"""
    k = np.atleast_1d(k)
    W = np.array([wslice(i, k) for i in range(len(Ts))])
    return np.maximum.accumulate(W, axis=0)


KCHK = np.linspace(-0.6, 0.6, 1201)
n_repaired = int((np.array([wslice(i, KCHK) for i in range(len(Ts))]) < wslices(KCHK) - 1e-10).any(axis=1).sum())
T_LAST = Ts[-1]


def wtot(k, T):
    """total implied variance at log-moneyness k (vs forward to T) and maturity T; SVI slices evaluated analytically"""
    scalar = np.ndim(k) == 0
    W = wslices(k)
    if T <= Ts[0]:
        out = W[0] * T / Ts[0]
    elif T >= T_LAST:                                                  # beyond data: same implied-vol smile (synthetic)
        out = W[-1] * T / T_LAST
    else:
        i = np.searchsorted(Ts, T) - 1
        a = (T - Ts[i]) / (Ts[i + 1] - Ts[i])
        out = (1 - a) * W[i] + a * W[i + 1]
    return float(out[0]) if scalar else out


def dupire_var(k, T, eps=1e-4):
    """Dupire local variance (deterministic-rates) in terms of total implied variance"""
    w = wtot(k, T); dk = 1e-3
    wT = (wtot(k, T + eps) - wtot(k, max(T - eps, 1e-6))) / (T + eps - max(T - eps, 1e-6))
    wk = (wtot(k + dk, T) - wtot(k - dk, T)) / (2 * dk)
    wkk = (wtot(k + dk, T) - 2 * w + wtot(k - dk, T)) / dk ** 2
    den = 1 - k / w * wk + 0.25 * (-0.25 - 1 / w + k * k / (w * w)) * wk ** 2 + 0.5 * wkk
    return np.maximum(wT, 1e-8) / np.maximum(den, 1e-3)


# ------------------------------------------------------------------ time grid: uniform + option expiries
T_SYN = [1.5, 2.5, 3.5, 5.0]                                             # beyond the data: synthetic surface
exps = sorted(set(Ts) | set(T_SYN))
grid_t = np.unique(np.concatenate([np.linspace(0, T_MAX, int(T_MAX * SPY) + 1), exps]))
grid_t = grid_t[grid_t <= T_MAX + 1e-12]

G = 81
chol = np.linalg.cholesky(np.array([[1, rSd, rSf], [rSd, 1, rdf], [rSf, rdf, 1.0]]))
conv = lambda t, k, s: s * s / (2 * k * k) * (1 - np.exp(-k * t)) ** 2
var = lambda t, k, s: s * s * (1 - np.exp(-2 * k * t)) / (2 * k)
quart = lambda u: np.where(np.abs(u) < 1, (1 - u * u) ** 2, 0.0)


def xgrid(t):
    s = 0.075 * np.sqrt(max(t, 1e-3))
    return np.linspace(-5 * s, 5 * s, G)


def run(calibrate, sig_tab=None, seed=0, n_pairs=NP, price_at=None):
    """Simulate under Q^d. If calibrate, build sigma(t, x) step by step; otherwise use sig_tab."""
    rng = np.random.default_rng(seed)
    M = 2 * n_pairs
    lnS = np.full(M, np.log(S0)); xd = np.zeros(M); xf = np.zeros(M); intr = np.zeros(M)
    ld = np.zeros(M); zf = np.zeros(M); lf = np.zeros(M); zd = np.zeros(M)
    tab, floor_hits, prices = [], 0, {}
    g_last = np.zeros(1)
    sig_x = xgrid(0.0); sig_v = np.sqrt(dupire_var(sig_x, 1e-3))
    for j in range(len(grid_t) - 1):
        t, t1 = grid_t[j], grid_t[j + 1]; dt = t1 - t
        lnFt = lnF(t)
        if calibrate:
            xg = xgrid(t)
            dup = dupire_var(xg, max(t, 1e-3))
            if t > 0 and EST == "malliavin":
                lam = ld - lf
                D = np.exp(-intr)
                h = max(1.2 * 0.075 * np.sqrt(t) * (2 * n_pairs) ** (-0.2), 1e-4)
                x = lnS - lnFt
                order = np.argsort(x); xs = x[order]
                num = np.zeros(G); den = np.zeros(G); den2 = np.zeros(G)
                lo = np.searchsorted(xs, xg - h); hi = np.searchsorted(xs, xg + h)
                for g in range(G):
                    if hi[g] > lo[g]:
                        idx = order[lo[g]:hi[g]]
                        kw = D[idx] * quart((x[idx] - xg[g]) / h)
                        num[g] = np.sum(kw * lam[idx]); den[g] = np.sum(kw); den2[g] = np.sum(kw * kw)
                neff = np.where(den2 > 0, den ** 2 / np.maximum(den2, 1e-300), 0.0)
                ok = (den > 0) & (neff >= max(MIN_NEFF, 1))                                            # trust the estimate only where enough particles contribute
                ok[: np.argmax(ok)] = False; ok[G - np.argmax(ok[::-1]):] = False    # one contiguous region around the money
                corr = np.interp(xg, xg[ok], num[ok] / den[ok])             # flat beyond particle support
                if SMOOTH:                                                  # the correction is smooth in K: regularise it
                    corr = savgol_filter(corr, SMOOTH, 2, mode='nearest')
                v = dup - 2 * corr
                if DIAG and (j % 25 == 0 or (t > 3.7 and j % 5 == 0)):
                    q = np.quantile(np.abs(lam), [0.5, 0.99, 0.9999])
                    print(f'    t={t:5.2f}  |lam| median {q[0]:.2e} p99 {q[1]:.2e} p99.99 {q[2]:.2e} max {np.abs(lam).max():.2e}'
                          f'  corr ATM {corr[G//2]:+.2e}  corr range [{corr.min():+.2e},{corr.max():+.2e}]  g range [{g_last.min():+.3f},{g_last.max():+.3f}]', flush=True)
            else:
                v = dup.copy()
            hit = v < VAR_FLOOR
            floor_hits += int(hit.sum())
            if hit.any():
                sd_x = 0.075 * np.sqrt(max(t, 1e-3))
                FLOOR_LOG.append((t, float(np.abs(xg[hit] / sd_x).min()), int(hit.sum())))
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
        if calibrate and EST == "malliavin":
            vd, vf = var(t, kd, sd), var(t, kf, sf)
            ld_n = ld * (1 + g * (dW - sig * dt) - kd * dt) + (vd - zf + rSd * sd * sig) * dt
            zf_n = zf + (-(kd + kf) * zf - rSf * sf * g * ld + rdf * sd * sf) * dt
            lf_n = lf * (1 + g * dW - kf * dt) + (zd - vf + rSf * sf * sig) * dt
            zd_n = zd + (-(kd + kf) * zd + rSd * sd * g * lf + rdf * sd * sf) * dt
            ld, zf, lf, zd = ld_n, zf_n, lf_n, zd_n
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
    return np.array(tab), floor_hits, prices


def black(F, K, T, s, D, cp):
    d1 = (np.log(F / K) + 0.5 * s * s * T) / (s * np.sqrt(T)); d2 = d1 - s * np.sqrt(T)
    return D * (F * norm.cdf(d1) - K * norm.cdf(d2)) if cp == "C" else D * (K * norm.cdf(-d2) - F * norm.cdf(-d1))


if __name__ == "__main__":
    t_start = time.time()
    print(f"[{label}] rates={RATES} {PARAMS}  estimator={EST}  particles={2*NP}  steps/yr={SPY}  smoothing window={SMOOTH}  min n_eff={MIN_NEFF}")
    print(f"  S0={S0:.5f}  USD zero 1y={-lnPd(1.0)*100:.2f}%  f_d-f_f at 1y={mu0(1.0)*100:.2f}%  "
          f"calendar-repaired slices={n_repaired}")
    tab, floors, _ = run(True, seed=SEED)
    t_cal = time.time() - t_start
    # independent repricing with frozen sigma table
    _, _, pr = run(False, sig_tab=tab, seed=SEED + 1000, n_pairs=NP, price_at=list(Ts) + T_SYN)
    t_all = time.time() - t_start
    rows = []
    for i, Te in enumerate(list(Ts) + T_SYN):
        S, D = pr[Te]
        Fm, Pd = np.exp(lnF(Te)), np.exp(lnPd(Te))
        h = len(S) // 2
        for delta in (0.10, 0.25, 0.50, 0.75, 0.90):        # call deltas: 10d call ... 10d put
            # strike for this call delta on the market smile
            f = lambda k: norm.cdf((-k + 0.5 * wtot(k, Te)) / np.sqrt(wtot(k, Te))) - delta
            k = brentq(f, -0.6, 0.6); K = Fm * np.exp(k); vol = np.sqrt(wtot(k, Te) / Te)
            cp = "C" if k >= 0 else "P"
            pay = D * (np.maximum(S - K, 0) if cp == "C" else np.maximum(K - S, 0))
            # control variates: E[D S] = Pd*F and E[D] = Pd, both exact in the model
            X = np.vstack([D * S - Pd * Fm, D - Pd]).T
            Xc = X - X.mean(0); keep = Xc.std(0) > 1e-12 * max(Xc.std(0).max(), 1e-300)
            keep &= Xc.std(0) > 1e-14
            beta = np.zeros(2)
            if keep.any():
                beta[keep] = np.linalg.lstsq(Xc[:, keep], pay - pay.mean(), rcond=1e-10)[0]
            adj = pay - X @ beta
            pp = 0.5 * (adj[:h] + adj[h:])
            price, se = pp.mean(), pp.std(ddof=1) / np.sqrt(h)
            vega = Pd * Fm * norm.pdf((-k + 0.5 * vol * vol * Te) / (vol * np.sqrt(Te))) * np.sqrt(Te)
            mkt = black(Fm, K, Te, vol, Pd, cp)
            rows.append(dict(T=Te, delta=delta, K=K, vol=100 * vol, err=100 * (price - mkt) / vega, se=100 * se / vega))
    R = pd.DataFrame(rows)
    piv = R.pivot_table(index="T", columns="delta", values="err").round(3)
    piv.columns = ["10dC", "25dC", "ATM", "25dP", "10dP"]
    piv.index = [f"{t:.3f}" + (" (synthetic)" if t > T_LAST + 1e-9 else "") for t in piv.index]
    piv["max se"] = R.groupby("T").se.max().round(3).values
    print("\n  repricing error, model minus market, in vol points (Monte Carlo standard error in last column)")
    print(piv.to_string())
    Rm = R[R["T"] <= T_LAST + 1e-9]
    print(f"\n  market expiries only: mean |err| {Rm.err.abs().mean():.3f}  max |err| {Rm.err.abs().max():.3f}  mean |err|/se {(Rm.err.abs()/Rm.se).mean():.2f}")
    print(f"  all rows: mean |err| {R.err.abs().mean():.3f}  max |err| {R.err.abs().max():.3f}  "
          f"mean |err|/se {(R.err.abs()/R.se).mean():.2f}  floor hits {floors}  "
          f"time: calibration {t_cal:.0f}s, total {t_all:.0f}s")
    if FLOOR_LOG:
        FL = np.array(FLOOR_LOG)
        inner = FL[FL[:, 1] < 3.0]
        print(f"  floor binds at {len(FL)} of {len(grid_t)-1} time steps; closest to the money: "
              f"{FL[:,1].min():.1f} std devs; steps where it binds inside 3 std devs: {len(inner)}"
              + (f" (t from {inner[:,0].min():.2f} to {inner[:,0].max():.2f})" if len(inner) else ""))
        atm1 = FL[FL[:, 1] < 1.0]
        if len(atm1):
            print(f"  floor binds within 1 std dev of the money from t = {atm1[:,0].min():.2f}")
    np.savez(f"calib_{label}.npz", grid_t=grid_t, tab=tab, R=R.to_records(index=False))
