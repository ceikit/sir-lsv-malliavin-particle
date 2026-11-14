# Calibration to the CME EUR/USD surface and repricing

Trade date 12 November 2026. Data: CME EUR/USD monthly and quarterly options (root EUU) and 6E futures, official settlements pulled from Databento. The data itself is licensed and is not in this repository.

## What was run

1. **Surface.** Implied volatilities from settlements (Black-76 on the underlying future), out-of-the-money options with premium of at least 3 ticks. The discount factor for each expiry comes from put-call parity: it holds within 1.8 ticks at every expiry, and the implied forward is within 0.00002 of the futures settlement. One SVI smile per expiry; fit error 0.01 to 0.13 vol points.
2. **Data cleaning.** The three 2028 expiries were dropped: they have almost no open interest, and their settled smiles are calendar-inconsistent with the 2027 smiles across 27% to 48% of the strike range within three standard deviations. A smaller violation between the September and October 2027 smiles, between 2.7 and 3 standard deviations, is repaired by taking a running maximum of total variance. The market range is therefore December 2026 to December 2027. Beyond that, the surface is extended at constant implied-vol smile to 5 years. **The extension is synthetic, not market data.**
3. **Calibration.** Particle method of the paper, with the integrated Malliavin sensitivities, under the domestic (USD) risk-neutral measure. 200,000 particles, 100 time steps per year plus every option expiry, quartic kernel.
4. **Repricing.** Independent Monte Carlo with the calibrated local volatility frozen: 200,000 paths, a different seed, control variates on E[D] and E[D S]. Error = (model price − market price) / vega, at 10-delta and 25-delta calls and puts and at the money.

The rate parameters are illustrative, not calibrated to caps or swaptions (none were used):

| Case | κ_d | σ_d | κ_f | σ_f | ρ_Sd | ρ_Sf | ρ_df |
|---|---|---|---|---|---|---|---|
| base | 3% | 1.0% | 3% | 0.8% | −0.30 | +0.12 | +0.63 |
| stress | 3% | 1.5% | 3% | 1.2% | +0.30 | −0.30 | +0.40 |

The base correlations are the historical 2012-2014 estimates quoted by Cozma, Mariapragassam and Reisinger (2019).

## Results

Repricing error in vol points, mean of absolute values over the five strikes per expiry, Monte Carlo standard error of a single error 0.02 to 0.07:

| | deterministic rates | base, no correction | base, with correction | stress, no correction | stress, with correction |
|---|---|---|---|---|---|
| market expiries, mean | 0.023 | 0.099 | 0.024 | 0.239 | 0.023 |
| market expiries, max | 0.134 | 0.184 | 0.135 | 0.489 | 0.134 |
| 5y (synthetic), max | 0.102 | 0.613 | 0.133 | 2.892 | 0.537 |

The deterministic-rates column is the error floor of the pipeline itself (time discretisation, Dupire from fitted smiles). The largest single error, −0.13 at the 10-delta put of the 3-week expiry, is time discretisation: it falls to −0.06 with 400 steps per year.

Same random numbers are used in every repricing run, so the error due to the stochastic rates is isolated by differencing against the deterministic-rates run:

| | market expiries, mean | market expiries, max | 1.5y to 3.5y, max | 5y, max |
|---|---|---|---|---|
| base, with correction | 0.001 | 0.005 | 0.017 | 0.030 |
| base, no correction | 0.098 | 0.186 | 0.504 | 0.638 |
| stress, with correction | 0.002 | 0.005 | 0.027 | 0.435 |
| stress, no correction | 0.233 | 0.471 | 1.842 | 2.830 |

As a sanity check, the calibrated correction at the money is compared with the flat-volatility closed form of the paper's check (b), evaluated with the calibrated at-the-money volatility. In the base case the calibrated correction is 93% of the closed form at 1 year, 87% at 2 years, 80% at 3 years and 76% at 4 years. The two are not expected to agree exactly, because the smile is not flat.

## Findings

- **The correction works over the market range.** With it, the stochastic rates add at most 0.005 vol points of error up to December 2027, against 0.19 (base) and 0.47 (stress) without it.
- **The positivity caveat is real.** In the stress case, the rates generate more forward variance than the surface contains, and the calibrated local variance hits the floor within one standard deviation of the money from t ≈ 3.0 years. No local volatility model with those rate parameters reproduces the (synthetic) long end, and the 5-year at-the-money error is 0.54 vol points.
- **Feedback instability seen, but traced to a bug.** A first version became unstable after about 2.6 years: a few particles developed very large sensitivities, which fed back through the estimated local-vol slope. The cause was my surface code (smiles tabulated on a clipped strike grid, giving noisy second derivatives). With the smiles evaluated analytically the instability disappeared without any smoothing. An optional minimum-effective-particle guard is kept in the script; it was not needed afterwards.

## Not done

- Rate parameters calibrated to cap or swaption volatilities.
- Futures-forward convexity and the delivery-date mismatch: options are treated as options on the forward to their expiry.
- Comparison with the conventional estimator and with a control-variate particle method.
- A real long-dated surface beyond December 2027.

## Scripts

From the repository root, with the settlements file from `data/pull_eurusd_options_databento.py`:

```
python calibration/build_surface.py path/to/eurusd_cme_options_2026-11-12.csv   # implied vols, parity, SVI fits, arbitrage checks
python calibration/calibrate.py det --rates none --estimator none                # deterministic rates
python calibration/calibrate.py base --rates base                                # base rates, with correction
python calibration/calibrate.py base_nocorr --rates base --estimator none        # base rates, no correction
python calibration/calibrate.py stress --rates stress
python calibration/calibrate.py stress_nocorr --rates stress --estimator none
python calibration/summarise.py                                                  # the tables above
```

or `make calibration SETTLEMENTS=path/to/eurusd_cme_options_2026-11-12.csv`. The defaults are 100,000 antithetic pairs, 100 steps per year and seed 1; each run takes about half a minute. The 400-steps figure quoted above is `calibrate.py det400 --rates none --estimator none --steps-per-year 400`.

| File | What it does |
|---|---|
| `build_surface.py` | Settlements to implied vols, put-call parity, one SVI smile per expiry. Writes its tables to `calibration/work/`. |
| `market.py` | Forward curve, discount curve, implied-variance surface and Dupire local variance, read from those tables. |
| `calibrate.py` | Particle calibration and independent repricing. Writes `calib_<label>.npz`. |
| `summarise.py` | The tables of this report from the five runs. |

Everything in `calibration/work/` is derived from licensed data and is git-ignored. A joint arbitrage-free fit of the surface (eSSVI) was tried and rejected, because it misfit the 2027 risk reversals by up to 1.6 vol points; the script is kept in `archive/essvi/`.

## Revision note

The numbers above were regenerated on 14 November 2026 after the code was reorganised, and are unchanged to the three decimals shown. The first version of this report was produced with discount factors read from a table rounded to five decimals that the committed scripts did not write; `build_surface.py` now writes that table itself, unrounded.
