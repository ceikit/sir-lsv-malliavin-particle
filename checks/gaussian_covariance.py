"""Stable Hull-White Gaussian step integrals and covariance square roots.

Positive quadrature evaluates the covariance as a Gram matrix. QR factorization
of its rectangular factor avoids forming a nearly singular covariance before
sampling; no eigenvalue clipping or diagonal jitter is used.
"""
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss


def bond_loading(kappa, remaining):
    """(1-exp(-kappa*remaining))/kappa, including its kappa=0 limit."""
    remaining = np.asarray(remaining, dtype=float)
    if not np.isfinite(kappa) or kappa < 0 or not np.isfinite(remaining).all() or (remaining < 0).any():
        raise ValueError('Require finite nonnegative kappa and remaining time')
    if kappa == 0:
        return remaining
    return -np.expm1(-kappa * remaining) / kappa


@lru_cache(maxsize=32)
def _quadrature(order):
    nodes, weights = leggauss(order)
    return nodes, weights / 2


def _step_kernels(t0, t1, terminal, kd, kf):
    if not np.isfinite([t0, t1, terminal, kd, kf]).all() or not (0 <= t0 < t1 <= terminal) or min(kd, kf) < 0:
        raise ValueError('Require finite 0 <= t0 < t1 <= terminal and nonnegative kappas')
    dt = t1 - t0
    order = max(16, int(np.ceil(max(kd, kf) * dt)) + 16)
    nodes, weights = _quadrature(order)
    # Evaluate relative times directly, avoiding cancellation near the endpoint.
    remaining_step = .5 * dt * (1 - nodes)
    remaining_bond = (terminal - t1) + remaining_step
    ed, ef = np.exp(-kd * remaining_step), np.exp(-kf * remaining_step)
    ad, af = bond_loading(kd, remaining_bond), bond_loading(kf, remaining_bond)
    return weights, np.array([np.ones_like(ed), ed, ad, ef, af])


def step_coefficients(t0, t1, terminal, kd, kf):
    """Step averages, without the dt multiplier, in the manuscript notation."""
    weights, (_, ed, ad, ef, af) = _step_kernels(t0, t1, terminal, kd, kf)
    avg = lambda x: float(weights @ x)
    return dict(J1d=avg(ed), J2d=avg(ed*ed), I1d=avg(ad), I2d=avg(ad*ad), Ld=avg(ed*ad),
                J1f=avg(ef), J2f=avg(ef*ef), I1f=avg(af), I2f=avg(af*af), Lf=avg(ef*af),
                Jdf=avg(ed*ef), L1df=avg(ad*ef), L2df=avg(ed*af), Kdf=avg(ad*af))


def step_covariance_factor(t0, t1, terminal, kd, kf, correlation):
    """Return L with L L.T = Sigma (step-average covariance, excluding dt).

    correlation is the positive-definite 3x3 driver matrix ordered (spot,d,f).
    QR acts on the Gram factor directly, retaining small conditional variances
    that can disappear if Sigma is explicitly formed and then Cholesky factored.
    """
    correlation = np.asarray(correlation, dtype=float)
    if (correlation.shape != (3, 3) or not np.isfinite(correlation).all()
            or not np.allclose(correlation, correlation.T, rtol=0, atol=1e-14)
            or not np.allclose(np.diag(correlation), 1, rtol=0, atol=1e-14)):
        raise ValueError('Require a symmetric 3x3 correlation matrix with unit diagonal')
    driver_factor = np.linalg.cholesky(correlation)
    weights, kernels = _step_kernels(t0, t1, terminal, kd, kf)
    driver_rows = driver_factor[[0, 1, 1, 2, 2]]
    gram_factor = (kernels[:, :, None] * np.sqrt(weights)[None, :, None]
                   * driver_rows[:, None, :]).reshape(5, -1)
    _, triangular = np.linalg.qr(gram_factor.T, mode='reduced')
    factor = triangular.T
    signs = np.where(np.diag(factor) < 0, -1., 1.)
    return factor * signs[None, :]
