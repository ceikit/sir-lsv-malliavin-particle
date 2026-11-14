"""Hull-White quantities with constant mean reversion k and volatility s.

The short rate is r_t = f(0,t) + convexity(t) + x_t, with x the centred factor.
"""
import numpy as np


def convexity(t, k, s):
    """s^2/(2 k^2) (1 - e^{-k t})^2: risk-neutral mean of r_t minus the forward rate f(0,t)."""
    return s * s / (2 * k * k) * (1 - np.exp(-k * t)) ** 2


def variance(t, k, s):
    """v(t) = s^2 (1 - e^{-2 k t})/(2 k): variance of the factor x_t."""
    return s * s * (1 - np.exp(-2 * k * t)) / (2 * k)
