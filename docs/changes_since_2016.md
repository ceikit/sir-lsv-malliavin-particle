# Changes with respect to the 2016 draft

This list was an appendix of the paper until the reorganisation of 14 November 2026.

Corrections to the 2016 text:

- One stochastic rate, general model: the centring rate is the forward rate, the expectation of the short rate under the t-forward measure (it was written as an expectation under the risk-neutral measure).
- One stochastic rate, Hull-White: exponent in theta(t); definition of the forward-measure Brownian motion; Xi is defined with the rate volatility; sign of the rate term in the SDE of the forward, and a missing rate-volatility factor in front of the stochastic integral of the discretised forward; initialisation and time indexing of the algorithm.
- Spot SDE of the three-factor model: the diffusion term is sigma(t, S) S dW.
- Calibration condition: unchanged, but also stated under one measure. The second system of particles under the foreign forward measure, and with it the section on exact integration under that measure, are removed.
- Change of measure from domestic to foreign: the spot driver is shifted by sigma(t, S), not by the spot-domestic correlation times sigma(t, S).
- Foreign rate: centred at the market forward rate. The terms in the expectation of sigma disappear, as does a drift term that had been computed with the domestic mean reversion in place of the foreign one.
- Discretised forward: sign of the drift, the domestic variance coefficient, and the cross coefficient.
- General Malliavin formula with two rates: the foreign term has the factor 2/K (it was 2 S0/K^2).
- The lemmas on the Malliavin derivatives of the two rates and of the spot, including the "second approach", are replaced by the propositions on the integrated sensitivities and on the four-process representation.

Added since 2016:

- The equivalent one-process form for one stochastic rate.
- The standing hypotheses and the smooth-payoff duality lemma; the consistency checks; the calibration tests; the existence section; the moment-closure section; scope, caveats and outlook.
- Bibliography: a duplicate entry removed, eleven references added.

Reorganisation of 14 November 2026: introduction first; the Malliavin tools before their first use; the one-rate material as background, with its time discretisation in an appendix; the three-factor material split into model, main result, algorithm, numerical evidence, existence and deterministic formulation, limits.
