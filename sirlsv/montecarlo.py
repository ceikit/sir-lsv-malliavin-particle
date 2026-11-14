"""Monte Carlo helpers."""
import numpy as np


def antithetic_mean(x):
    """Mean and standard error of a sample whose second half is the antithetic copy of the first."""
    h = len(x) // 2
    p = 0.5 * (x[:h] + x[h:])
    return p.mean(), p.std(ddof=1) / np.sqrt(h)
