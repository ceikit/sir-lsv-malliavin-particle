# Particle calibration of local volatility with stochastic rates: long-dated FX

Working draft by Filippo Perugini (first written in 2016, revised November 2026), with the code that goes with it.

The paper derives the particle (McKean) calibration of local volatility models with stochastic interest rates, following Guyon and Henry-Labordère, and extends it to FX with Hull-White domestic and foreign rates.

## Main result

Under the paper's regularity, integrability and forward-PDE uniqueness assumptions, an existing three-factor FX local volatility model with Hull-White domestic and foreign rates is calibrated to the market smile if and only if

$$\sigma(t,K)^2=\sigma_{Dup}(t,K)^2-2\,\mathbb{E}^{d^t}\big[\ell^d_t-\ell^f_t\,\big|\,S_t=K\big]$$

where $\ell^d$ and $\ell^f$ are two of four auxiliary processes carried by each particle and simulated forward with the state ("Reduction to four scalar processes" in the paper, `sirlsv/sensitivities.py` in the code). The correction contains neither $\partial_K^2 C$ nor tail expectations, and a single system of particles under the domestic forward measure is enough.

## Layout

| Path | What it is |
|---|---|
| `paper/` | The paper. `SIR_LSV_particle_FX.tex` is the main file; the text is in `sections/`, the bibliography in `references.tex`. `notes/` holds two standalone notes: the moment-closure section on its own (with its PDF) and a short note on a cancellation in that system. |
| `sirlsv/` | Shared Python code: the recursion for the four auxiliary processes, Hull-White formulas, the stable Gaussian step covariance, and the prescribed-volatility model used by the checks. |
| `checks/` | Scripts that test the building blocks with a prescribed local volatility (nothing is calibrated). |
| `calibration/` | Particle calibration to the CME EUR/USD surface of 12 November 2026 and repricing. Report: `calibration/RESULTS.md`. |
| `benchmarks/` | Synthetic comparison of three calibration methods (two 3D PDE solvers and the particle method) on a hidden non-flat volatility. Report: `benchmarks/RESULTS.md`. Needs a C++ compiler. |
| `data/` | Script that pulls one day of CME EUR/USD option settlements from Databento (needs your own API key). |
| `docs/reviews/` | Notes by a second agent: audit of the paper and code, comparison with Bloch and Nakashima (2008), novelty assessment. The index gives the status of each finding. |
| `archive/` | Superseded material: the first four-page note with its scripts, and a surface fit that was tried and rejected. |

The checks:

| Script | What it tests |
|---|---|
| `check_sensitivities.py` | The two identities behind the main formula, by Monte Carlo, for a local volatility with skew. |
| `check_forward_measure.py` | The simulation scheme under the domestic T-forward measure, and the particle weights. |
| `check_matrix_form.py` | The matrix form and the scalar form of the auxiliary processes agree path by path. Self-contained on purpose. |
| `compare_estimators.py` | Noise of the conventional estimator against the new one, strike by strike, at fixed local volatility. |
| `check_covariance_stability.py` | Deterministic regression test of the Gaussian step covariance (54 cases, including zero mean reversion and tiny steps). |
| `check_gaussian_existence.py` | The flat-smile construction of the existence section, against independent integration. |
| `check_moment_closure.py` | The first-moment closure of the note: symbolic check of the generator identity and of the two feedback reductions (needs `sympy`), then Monte Carlo check of the weak balances. |
| `check_variance_derivative_cancellation.py` | Symbolic check that the two main moment equations contain no second derivative of the unknown variance (needs `sympy`). |

## Running

Python 3 with the packages in `requirements.txt` (`numpy` 2 or later, `scipy`, `pandas`; `sympy` for two checks, `matplotlib` and `g++` for the benchmark). Run everything from the repository root.

```
make checks-quick      # all checks at a reduced size, about a minute
make checks            # the sizes quoted below, several minutes
make paper             # pdflatex, three passes
make notes             # the two standalone notes in paper/notes
make benchmark         # the synthetic three-method comparison from scratch, several minutes
make calibration SETTLEMENTS=path/to/eurusd_cme_options_2026-11-12.csv
```

or call the scripts directly, for example

```
python checks/check_sensitivities.py 1000000 200 5.0     # pairs, steps per year, maturity
python checks/check_forward_measure.py 200 500000        # steps per year, pairs
python checks/compare_estimators.py 1000000 100 20000    # pairs, steps per year, particles per batch
```

Each script prints what it tests and the difference between the two sides in standard errors. The figure `volatility.jpg` from the 2016 draft is not in the repository; the paper compiles without it and shows a placeholder.

## Status

Checked:

- The formula reduces to the one-factor Hull-White formula of Guyon and Henry-Labordère when there is no foreign rate, and matches the Gaussian closed form for constant spot volatility.
- With a skewed local volatility (2 million paths), both identities hold within one standard error; dropping the measure-change terms puts the foreign term about 2.7% off.
- At fixed local volatility and 20,000 particles, the standard deviation of the new estimator of the rates correction was 14 to 66 times lower than that of the conventional estimator without control variates, across strikes from the 1% to the 99% quantile. This is not a calibration benchmark.
- Full particle calibration to CME EUR/USD option settlements of 12 November 2026 (expiries December 2026 to December 2027), with independent repricing and illustrative rate parameters. With the correction, the stochastic rates add at most 0.005 vol points of repricing error over the market expiries, against 0.19 (base rate parameters) and 0.47 (stressed) without it. In the stressed case the calibrated variance hits its floor near the money from about 3 years, on the synthetic part of the surface. Report: `calibration/RESULTS.md`.
- A deterministic formulation (Section "Exact closure of discounted first moments" of the paper, by the second agent): for a prescribed local volatility, the discounted density and the discounted first moments of the four auxiliary processes solve five linear forward equations in the three state variables, with no higher moments. Dividing the correction by the market density instead of the model density gives an equivalent calibration condition for regular solutions. The generator identity is verified symbolically; in the Monte Carlo test (120,000 paths) the weak balances and compatibility identities hold within 1.6 standard errors.
- A synthetic comparison on spot options (`benchmarks/`, by the second agent): a hidden non-flat volatility is recovered from the same target by a conventional 3D density PDE, by the five-field PDE of the closure, and by the particle method. On their finest settings all three reach a root-mean-square error of 0.002 to 0.003 vol points and reprice within 0.004. The particle method is the cheapest: 20,000 particles give 0.004 in about 3 seconds, against 69 and 189 seconds for the finest PDE grids. One parameter set; at that accuracy the reference target's own error is not negligible, so the methods cannot be ranked.

Not done:

- Market data beyond December 2027 (the surface is extended synthetically to 5 years), and rate parameters calibrated to caps or swaptions.
- A comparison with control-variate particle methods.
- General-smile existence and uniqueness remain unproved; the existence section proves a flat-smile result and a separate regularized-system theorem. The local variance can become negative for some surfaces and rate parameters (see "Scope and caveats" in the paper).
- Prior art. Malliavin integration by parts for the stochastic-rates correction, the cancellation of the density, the domestic-minus-foreign form of the correction and a fixed-point calibration are already in Bloch and Nakashima (2008), who credit Benhamou, Rivoira and Gruz (2008). According to the equation-level comparison in `docs/reviews/`, their implemented two-rate method freezes the drift locally and evaluates the volatility at a fixed strike, and has no equivalent of the four forward equations of this paper. The candidate contribution is therefore narrow: the exact four-process system for the coupled model, with the foreign-measure terms, and the calibration equivalence. Benhamou, Rivoira and Gruz has still not been read in full, and no claim of priority is made.
- An existence theorem for general smiles. A cancellation in the moment system (`paper/notes/variance_derivative_cancellation.tex`) is a possible first step; the estimates are missing.

Known issues:

- The calibration data are options on EUR/USD futures. They are treated as options on the forward to their expiry; futures-forward convexity and the delivery-date mismatch are ignored.
- `data/pull_eurusd_options_databento.py` does not process `update_action`, so a settlement that was published and then deleted would still be used. The script has not been tested against the live API.
- The existence section and the hypotheses added to the propositions come from the second agent's audit (`docs/reviews/`). Its numerical checks pass and the flat-smile and log-contract arguments have been re-derived; the regularized well-posedness proof has been read but not checked line by line. The same holds for the moment-closure section: its equations have been re-derived by hand and match `sirlsv/sensitivities.py`, while the fixed-coefficient uniqueness proof is a sketch of a standard duality argument that has been read, not checked in detail.
- The comparison with Bloch and Nakashima was made by the second agent from the full paper. Its mathematical objections have been checked; its description of what that paper contains has not, because the paper is not in this repository.
- The benchmark has no conventional particle estimator as a baseline, and its intermediate tables (`benchmarks/work/`, about 27 MB) are not committed: `make benchmark` regenerates them.

## Data

Market data pulled with the Databento script is licensed and must not be committed here. `.gitignore` excludes the settlements CSV and the whole of `calibration/work/`, where the calibration scripts write the tables they derive from it.
