"""Parameters of the three-factor model and the local volatility used in the checks."""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RateParams:
    """Hull-White parameters and correlations. d = domestic, f = foreign, S = spot."""
    kd: float    # mean reversion of the domestic rate
    sd: float    # volatility of the domestic rate
    kf: float    # mean reversion of the foreign rate
    sf: float    # volatility of the foreign rate
    rSd: float   # correlation spot / domestic rate
    rSf: float   # correlation spot / foreign rate
    rdf: float   # correlation domestic rate / foreign rate

    def correlation(self):
        """3x3 correlation matrix of the drivers, ordered (spot, domestic, foreign)."""
        return np.array([[1, self.rSd, self.rSf], [self.rSd, 1, self.rdf], [self.rSf, self.rdf, 1.0]])


def tanh_local_vol(x, s0, a, b):
    """sigma = s0 (1 - a tanh(b x)) at x = ln(S/S0), and g = S dsigma/dS = dsigma/dx."""
    th = np.tanh(b * x)
    return s0 * (1 - a * th), -s0 * a * b * (1 - th * th)
