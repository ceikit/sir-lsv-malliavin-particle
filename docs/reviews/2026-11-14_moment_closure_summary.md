# Exact discounted first-moment closure: derivation and checks

> Write-up by the second agent (Codex), reproduced as received. Paths in it refer to that agent's working tree.
> In this repository: the note is `paper/discounted_moment_closure_note.pdf`, the script is `checks/check_moment_closure.py`,
> and the validation record is `docs/reviews/2026-11-14_moment_closure_validation.json` (regenerated here, same numbers).
> The LaTeX source of the note is not in this repository yet, so the main paper does not include the section.


Recorded: 2026-11-14T10:29:01Z. Developed from the proposal in the publication
assessment, following the request to exploit the exact closure.

## Deliverables

- paper/discounted_moment_closure.tex: mathematical section, included by the
main manuscript. Contains the explicit equations and four results with proofs.
- paper/discounted_moment_closure_note.tex: standalone compilation wrapper.
- paper/discounted_moment_closure_note.pdf: standalone readable note.
- checks/check_moment_closure.py: reproducible symbolic and Monte Carlo checks.
- reviews/2026-11-14_moment_closure_validation.json: machine-readable output.

The results concern the original three-factor Gaussian-rate local-volatility
model and its four integrated sensitivities, including genuine nonflat volatility.
They are mathematical derivations, with code for verification. They are not a
new implementation of the full nonlinear calibration algorithm.

## What is proved

1. Weak exact closure. For a prescribed deterministic local volatility,
the discounted physical measure p and signed moment measures m_i solve
five linear forward equations on three physical coordinates. The derivation
retains the shared-noise divergence terms. No second auxiliary moments occur.
This weak result only needs the stated solution/integrability assumptions;
densities and a nonsingular correlation matrix are not required.
2. Identification and uniqueness at fixed coefficient. Under explicit
stronger smoothness, bounded-volatility and nondegeneracy assumptions, the
stochastic moments are the unique solution in an exponentially weighted
measure class. The proof uses the backward evolution of affine auxiliary
payoffs and cutoff duality. This result covers nonflat volatility, but fixes
the coefficient in advance.
3. Full-state compatibility identities. Own-forward Malliavin duality gives
two differential constraints on the moment fields and density. Rate
integration yields J_d = -partial_x n_1 and
J_f = -partial_x n_3 - n_3. The extra foreign term is essential.
These identities are also useful diagnostics for a numerical PDE solver.
4. Nonlinear equivalence, including a known denominator. A closed-system
solution whose induced coefficient belongs to the stated regularity class
is a genuine calibrated model if it satisfies either
  s^2 = q - 2(n_1-n_3)/p_bar
  or
  s^2 = q - 2(n_1-n_3)/p_market.
  Here p_market(t,x) = exp(x) C_KK(t,exp(x)) is the known discounted
log-spot density, not a probability density. In the second case the difference
w = p_bar-p_market satisfies a homogeneous linear forward equation with
coefficient s, zero initial data, and deterministic discount killing.
Uniqueness forces w=0. Thus the substitution does not assume calibration
in advance. The theorem identifies regular solutions; it does not construct
a solution of the nonlinear system.

## Why the closure works

The auxiliary SDE is affine in its four auxiliary variables. Its joint generator
maps functions affine in those variables to functions of the same form. The
quadratic auxiliary diffusion contributes no second moments because the second
auxiliary derivatives of those test functions vanish. Mixed physical/auxiliary
derivatives do survive and produce the divergence terms.

The resulting formulation uses five fields on three physical coordinates,
rather than a joint density on seven coordinates. Discounting is a killing term;
no extra discount-factor coordinate is needed. This does not claim that older
stochastic-rate calibration required a seven-dimensional PDE: the physical
density alone already determines the conventional correction via tail integrals.
The benefit is an exact alternative evolution for the particular sensitivity
moments and a framework for further analysis or deterministic cross-checks.

## Checks performed

The symbolic test constructs the full seven-dimensional drift and covariance
directly from the three correlated Brownian drivers. It then checks the generator
on a generic affine auxiliary test function. The resulting coefficients agree
exactly with the proposed weak system. It separately verifies the marginal
equation and both feedback reductions, including the homogeneous equation for
the market-denominator formulation.

The nonflat simulation uses 120,000 independent paths, 600 Euler steps over
three years, seed 271828, and

sigma(S) = 0.18 (1 - 0.6 tanh(2 log S)).

Rate parameters are kappa_d=0.07, kappa_f=0.13, nu_d=0.04, nu_f=0.05,
with correlations rho_Sd=0.36, rho_Sf=-0.28, rho_df=0.25.
Each weak balance compares a terminal discounted test payoff with its integrated
generator. Standard errors use pathwise paired residuals.

|Identity                         |Paired z-score|
|---------------------------------|-------------:|
|Discounted physical measure      |-0.977        |
|Domestic log sensitivity moment  |-1.286        |
|Domestic-to-foreign cross moment |-0.510        |
|Foreign log sensitivity moment   |+0.858        |
|Foreign-to-domestic cross moment |-1.586        |
|Domestic full-state compatibility|+0.412        |
|Foreign full-state compatibility |-0.810        |

As a diagnostic, omitting the shared-noise flux gives z-scores approximately
-60.4 and +51.5 for the two log-sensitivity moment equations. This strongly
detects that particular incorrect formula in this experiment. It is not a
performance comparison against another calibration method.

The simulation retains Euler time-discretization error; these statistics are
consistency evidence, not a convergence theorem. The symbolic generator and
analytical proofs concern the continuous-time equations.

Both LaTeX entry points were compiled and the new note was visually inspected.
The local TeX installation needs explicit format/font-map paths; this was handled
in the build environment, not by inserting machine-specific paths into the source.

## Scope, usefulness and next mathematical step

The known-market-denominator formulation removes division by an unknown marginal
from the deterministic feedback. It still divides by a possibly tiny market
density and differentiates the unknown coefficient. The first-moment equations
also contain divergences of terms involving partial_x s, so their nonlinear
principal part needs analysis; linear-system well-posedness does not supply it.
The time-zero Dirac law and the positive variance gap are additional obstacles.

The next substantive theoretical task is to prove estimates for this nonlinear
system, for example around an admissible nonflat or Gaussian reference in a
specified weighted function space. No contraction, general-smile existence,
nonlinear uniqueness, propagation of chaos, or numerical superiority is claimed.

Affine moment closure and Feynman–Kac duality are standard mechanisms. The
application, explicit operators, compatibility relations and equivalent
formulation are useful additions to this draft; their novelty against all
previous literature has not been established.
