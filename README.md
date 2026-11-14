# Particle calibration of local volatility with stochastic rates: long-dated FX

Working draft by Filippo Perugini (first written in 2016, revised November 2026), with the numerical checks that go with it.

The paper derives the particle (McKean) calibration of local volatility models with stochastic interest rates, following Guyon and Henry-Labordère, and extends it to FX with Hull-White domestic and foreign rates.

## Main result

For Hull-White domestic and foreign rates, the three-factor FX local volatility model is calibrated to the market smile if and only if

$$\sigma(t,K)^2=\sigma_{Dup}(t,K)^2-2\,\mathbb{E}^{d^t}\big[\ell^d_t-\ell^f_t\,\big|\,S_t=K\big]$$

where $\ell^d$ and $\ell^f$ are two of four auxiliary processes carried by each particle and simulated forward with the state (equation (40) of the paper). The correction contains neither $\partial_K^2 C$ nor tail expectations, and a single system of particles under the domestic forward measure is enough.

## Contents

| Path | What it is |
|---|---|
| `paper/SIR_LSV_particle_FX.tex` | The paper. Self-contained LaTeX source. |
| `checks/check_sensitivities.py` | Monte Carlo test of the two identities behind the main formula, for a local volatility with skew. |
| `checks/check_forward_measure.py` | Test of the simulation scheme under the domestic T-forward measure and of the particle weights. |
| `checks/check_matrix_form.py` | Shows that the matrix form and the scalar form of the auxiliary processes agree path by path. |
| `checks/compare_estimators.py` | Noise of the conventional estimator against the new one, strike by strike, at fixed local volatility. |
| `data/pull_eurusd_options_databento.py` | Pulls one day of CME EUR/USD option settlements from Databento (needs your own API key). |
| `archive/first_note/` | The first four-page note and its scripts. Superseded by the paper. |
| `archive/paper_build/` | The fragments and the script used to assemble the paper from the 2016 draft. They need the original draft and a results log, which are not in this repository. |

## Running the checks

Requirements: Python 3 with `numpy` (and `pandas`, `databento` for the data script).

```
python checks/check_sensitivities.py 1000000 200 5.0     # pairs, steps per year, maturity
python checks/check_forward_measure.py 200 500000        # steps per year, pairs
python checks/check_matrix_form.py 100
python checks/compare_estimators.py 1000000 100 20000    # pairs, steps per year, particles per batch
```

Each script prints what it tests and the difference between the two sides in standard errors.

## Building the paper

```
cd paper && pdflatex SIR_LSV_particle_FX.tex && pdflatex SIR_LSV_particle_FX.tex
```

The figure `volatility.jpg` from the 2016 draft is not in the repository; the source compiles without it and shows a placeholder.

## Status

Checked so far:

- The formula reduces to the one-factor Hull-White formula of Guyon and Henry-Labordère when there is no foreign rate, and matches the Gaussian closed form for constant spot volatility.
- With a skewed local volatility (2 million paths), both identities hold within one standard error; dropping the measure-change terms puts the foreign term about 2.7% off.
- At fixed local volatility and 20,000 particles, the standard deviation of the new estimator of the rates correction was 14 to 66 times lower than that of the conventional estimator without control variates, across strikes from the 1% to the 99% quantile.

Not done:

- A full calibration against a market smile followed by vanilla repricing. The estimator comparison above is not a calibration benchmark.
- A comparison with control-variate particle methods.
- Existence and uniqueness of the calibrated model are not addressed, and the local variance can become negative for some surfaces and rate parameters (see "Scope and caveats" in the paper).
- No claim of novelty is made; the literature has not been checked systematically.

## Data

Market data pulled with the Databento script is licensed and must not be committed here; `.gitignore` excludes the CSV it writes.
