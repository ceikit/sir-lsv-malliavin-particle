"""Market inputs of the calibration: forward curve, domestic discount curve and the
implied-variance surface, read from the files written by build_surface.py.

    futures.csv    futures settlements        -> forward curve  ln F(T), quadratic in T
    parity.csv     discount factors by expiry -> ln P^d(T), quadratic in T through 0
    svi_fits.csv   one SVI smile per expiry   -> total implied variance w(k, T)

Beyond the last market expiry kept, the surface is extended at constant implied-vol
smile. That extension is synthetic, not market data.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

WORK_DIR = Path(__file__).resolve().parent / "work"     # git-ignored: holds licensed data and files derived from it
MAX_MKT_T = 1.1     # Dec 2026 - Dec 2027; the 2028 expiries are dropped (no open interest, calendar-inconsistent settlements)


def black(F, K, T, s, D, cp):
    """Black-76 price of a call ("C") or put ("P") on a forward F, discount factor D."""
    d1 = (np.log(F / K) + 0.5 * s * s * T) / (s * np.sqrt(T)); d2 = d1 - s * np.sqrt(T)
    return D * (F * norm.cdf(d1) - K * norm.cdf(d2)) if cp == "C" else D * (K * norm.cdf(-d2) - F * norm.cdf(-d1))


class Market:
    def __init__(self, work=WORK_DIR, max_market_T=MAX_MKT_T):
        work = Path(work)
        fut = pd.read_csv(work / "futures.csv")
        self.cF = np.polyfit(fut["T"], np.log(fut.settlement), 2)               # ln F(T) = c2 T^2 + c1 T + c0
        self.S0 = float(np.exp(self.cF[2]))
        par = pd.read_csv(work / "parity.csv")
        cD = np.polyfit(par["T"], np.log(par["D"]), 2, w=np.sqrt(par["T"]))     # ln P^d(T), forced through 0 below
        cD[2] = 0.0
        self.cD = cD
        svi = pd.read_csv(work / "svi_fits.csv")
        svi = svi[svi["T"] <= max_market_T].reset_index(drop=True)
        self.Ts = svi["T"].values                                              # market expiries kept
        self.SV = svi[["a", "b", "rho", "m", "sig"]].values
        self.T_last = self.Ts[-1]
        k = np.linspace(-0.6, 0.6, 1201)
        raw = np.array([self._slice(i, k) for i in range(len(self.Ts))])
        self.n_repaired = int((raw < self.slices(k) - 1e-10).any(axis=1).sum())   # slices touched by the calendar repair

    def lnF(self, T):
        """log of the forward FX rate for maturity T"""
        return np.polyval(self.cF, T)

    def lnPd(self, T):
        """log of the domestic discount factor for maturity T"""
        return np.polyval(self.cD, T)

    def carry(self, t):
        """instantaneous forward rate differential f_d(0,t) - f_f(0,t)"""
        return 2 * self.cF[0] * t + self.cF[1]

    def _slice(self, i, k):
        a, b, rho, m, sig = self.SV[i]
        return a + b * (rho * (k - m) + np.sqrt((k - m) ** 2 + sig ** 2))

    def slices(self, k):
        """total variance of every market slice at log-moneyness k, with calendar arbitrage repaired by a running max"""
        k = np.atleast_1d(k)
        W = np.array([self._slice(i, k) for i in range(len(self.Ts))])
        return np.maximum.accumulate(W, axis=0)

    def total_variance(self, k, T):
        """total implied variance at log-moneyness k (vs forward to T) and maturity T; SVI slices evaluated analytically"""
        scalar = np.ndim(k) == 0
        Ts = self.Ts
        W = self.slices(k)
        if T <= Ts[0]:
            out = W[0] * T / Ts[0]
        elif T >= self.T_last:                                             # beyond data: same implied-vol smile (synthetic)
            out = W[-1] * T / self.T_last
        else:
            i = np.searchsorted(Ts, T) - 1
            a = (T - Ts[i]) / (Ts[i + 1] - Ts[i])
            out = (1 - a) * W[i] + a * W[i + 1]
        return float(out[0]) if scalar else out

    def dupire_variance(self, k, T, eps=1e-4):
        """Dupire local variance (deterministic rates) in terms of total implied variance"""
        wtot = self.total_variance
        w = wtot(k, T); dk = 1e-3
        wT = (wtot(k, T + eps) - wtot(k, max(T - eps, 1e-6))) / (T + eps - max(T - eps, 1e-6))
        wk = (wtot(k + dk, T) - wtot(k - dk, T)) / (2 * dk)
        wkk = (wtot(k + dk, T) - 2 * w + wtot(k - dk, T)) / dk ** 2
        den = 1 - k / w * wk + 0.25 * (-0.25 - 1 / w + k * k / (w * w)) * wk ** 2 + 0.5 * wkk
        return np.maximum(wT, 1e-8) / np.maximum(den, 1e-3)
