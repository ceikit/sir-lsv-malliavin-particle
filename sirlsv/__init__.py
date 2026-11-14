"""Shared code for the particle calibration of the three-factor FX model:
local volatility for the spot, Hull-White domestic and foreign short rates.

    model               rate parameters and the local volatility used in the checks
    hull_white          convexity adjustment and variance of the Hull-White factor
    sensitivities       Euler step of the four auxiliary processes of the paper
    prescribed_vol      simulation with a prescribed local volatility (used by the checks)
    gaussian_covariance stable step integrals and covariance factors for the exact scheme
    montecarlo          small Monte Carlo helpers
"""
