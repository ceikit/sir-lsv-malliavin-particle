# Prior art: first-hand reading of the two 2008 papers

14 November 2026. Both papers were read from the PDFs supplied by the author of this repository (SSRN 1107711 and SSRN 1153337). The PDFs are not in the repository.

## Benhamou, Gruz and Rivoira, "Stochastic interest rates for local volatility hybrids models" (March 2008, 10 pages)

- Setting: equity local volatility with jumps and one stochastic rate. No foreign rate.
- Proposition 2.1 is the tail-expectation form of the stochastic-rates correction.
- Appendix, proof of (2.7): Malliavin integration by parts turns the tail term into a term carrying a Dirac mass at the strike; dividing by the density gives a conditional expectation given the terminal spot. Equation (6.18) writes the correction as the time integral of the rate volatility times the correlation times a conditional expectation of the local volatility multiplied by a ratio of first variation processes. This is the correlation part of the process `l_d` of this repository, written as an integral instead of a forward equation.
- In that proof the Malliavin derivative of the log-spot in the direction of the rate driver is taken to be the correlation times its derivative in the direction of the spot driver, on the grounds that the bond price is adapted. The effect of the rate on the drift of the spot is therefore absent, and the correction vanishes at zero correlation. In the Hull-White model it does not (flat-smile check of the paper).
- The formula is not used as it stands: the covariance is approximated with the volatility evaluated at the strike, and the local volatility is obtained by a fixed-point iteration, (2.10) to (2.12), which the paper calls its main result.
- The paper attributes the result, for FX without jumps, to a 2005 conference talk by Balland. Not read.

## Bloch and Nakashima, "Multi-Currency Local Volatility Model" (June 2008, revised August 2008, 37 pages)

Read: Sections 5 to 9 and Appendix B. The description in `2026-11-14_bloch_nakashima_prior_art_comparison.md` is confirmed, with one refinement.

- Section 5 derives the stochastic-rates Dupire formula with domestic and foreign rates, heuristically.
- The paper's own method is not Malliavin calculus: Sections 6 and 7 match the maturity derivative of the variance of the log-spot with and without stochastic rates, then set the terminal spot equal to the strike. Equation (7.10) is the domestic-minus-foreign correction, written with conditional covariances.
- Section 7 works on a locally constant drift. Section 8.1 replaces the volatility along the path by the volatility at the strike and iterates; 8.2 solves a quadratic equation at each maturity. The rate-variance terms are present, as deterministic integrals. Section 9 sets the foreign rate to zero.
- Appendix B reproduces the Malliavin argument of Benhamou, Gruz and Rivoira, with the foreign rate set to zero, including the step discussed above.
- Refinement: the Malliavin integration by parts and the removal of the density should be credited to Benhamou, Gruz and Rivoira (and, in exact form with one rate, to Guyon and Henry-Labordere), not to Bloch and Nakashima, whose contribution is the two-rate setting, the form of the correction and the bootstrap.

## What is not in either paper

No forward system for the conditional quantity, no term from the change to the foreign measure, no proof that the resulting local volatility calibrates the model. The implemented formulas are the Gaussian (deterministic-volatility) correction evaluated with the local volatility frozen at the strike.
