# Theorem extensions and covariance validation

Recorded: 2026-11-14T09:42:29Z (UTC).

This implementation note supplements the earlier timestamped mathematical and
code audit. It records the scope of the new proofs and the checks actually run;
it makes no claim to have solved general-smile McKean--Vlasov existence.

## Mathematical changes

The main manuscript now separates smooth-payoff duality, pointwise calibration,
and existence. The calibration sufficiency proof uses the model-generated strike
density and uniqueness of the deterministic forward Dupire PDE with matching
initial and boundary conditions. The required flow, integrability, martingale,
density and uniqueness assumptions appear with the propositions.

The representation admits local spatial differentiability without a global
bound on the logarithmic volatility derivative. Smooth-payoff duality does not
require a positive density or an invertible correlation matrix. Weighted density
identities can be read almost everywhere at common Lebesgue points; classical
pointwise statements require additional density regularity. These relaxations
do not prove the existence of the assumed flow or cover arbitrary nonsmooth
local volatility.

The included `paper/existence_extensions.tex` adds:

- A conditional extension to deterministic time-dependent Gaussian rate
  parameters and correlations. The covariance implementation still uses constant
  rate parameters.
- Local existence for flat smile input in the original, unregularized equation,
  with uniqueness only within the deterministic-volatility Gaussian subclass.
  A positive variance gap gives continuation; zero spot/rate correlations admit
  an explicit sufficient condition for the full horizon.
- A necessary log-contract compatibility bound applying also to spot-dependent
  volatility, under the stated integrability assumptions. Thus arbitrary target
  surfaces and rate parameters cannot admit an unconditional existence theorem.
- Strong well-posedness of a precisely specified equation with kernel smoothing,
  denominator regularization, bounded features, state cutoffs and a smooth
  positive variance limiter. These modifications are part of that theorem.

The flat-smile construction is a benchmark, not the central research result.
It allows a term structure, but setting the spatial derivative to zero removes
the difficult smile feedback. No proof is supplied for the original general-smile
equation or for removal of the regularizations. A synthetic non-flat recovery
experiment and estimates for the unregularized limit remain outstanding.

## Numerical changes

`checks/gaussian_covariance.py` uses stable `expm1` loadings and positive
Gauss--Legendre quadrature. Sampling uses QR of the rectangular Gram factor,
avoiding formation of a nearly singular covariance matrix before factorization.
No eigenvalue clipping or covariance jitter is introduced. The driver correlation
matrix in this implementation must be positive definite, even though the
analytical representation permits positive semidefinite correlations.

The forward-measure checker also uses stable bond and variance formulas with
the exact zero-mean-reversion limits. Its CLI arguments are now parsed only
when the module is executed directly.

## Validation performed

- Covariance regression passed 54 parameter/time-step cases against independent
  dense integration of all fourteen coefficients, with QR reconstruction checks.
  Additional steps of 1e-6 and 1e-8 years preserve an analytically known small
  conditional variance. Zero-mean-reversion and invalid-input checks passed.
- The Gaussian existence checker integrated the two auxiliary ODEs and compared
  independently integrated forward variance at four maturities. Maximum absolute
  errors were 9.721e-13 with nonzero spot/rate correlations and 9.040e-14 with
  zero spot/rate correlations.
- The log-contract obstruction example gives target variance 0.075000 versus a
  necessary minimum of 0.512311 at thirty years with zero spot/rate correlations.
  Under the proposition's integrability assumptions this excludes spot-dependent
  solutions too; it is not merely failure of the Gaussian construction.
- The forward-measure Monte Carlo check at 100 steps/year and 200,000 paths gave
  sensitivity discrepancies below 0.99 paired standard errors. A separate
  50,000-path zero-mean-reversion run gave sensitivity discrepancies below 1.84
  standard errors. These are prescribed-volatility consistency checks, not a
  full particle calibration or a convergence proof.
- The LaTeX source compiled to 29 pages without undefined-reference or overfull
  box warnings. The added mathematical pages were visually inspected.

## Reproduction

```sh
python checks/check_covariance_stability.py
python checks/check_gaussian_existence.py
python checks/check_forward_measure.py 100 100000
```

The second command requires SciPy. The Monte Carlo command uses antithetic
pairs, so its 100,000 pairs correspond to 200,000 paths.
