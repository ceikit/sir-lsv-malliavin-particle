# Bloch–Nakashima (2008): equation-level prior-art comparison

Recorded: 2026-11-14T10:15:07Z. Our manuscript assessed at local commit
`936c7bd`. This note updates the earlier publication assessment after receipt
of the complete Bloch–Nakashima paper.

Source: Daniel Bloch and Yukio Nakashima, *Multi-Currency Local Volatility
Model*, working paper dated 19 June 2008, revised August 2008, 37 pages,
[SSRN 1153337](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1153337).
Page references below are the printed page numbers, which coincide with the
PDF page numbers in the supplied `ssrn-1153337.pdf`. The third-party PDF is
not added to this repository.

## Assessment

This is close and material prior art. The broad idea of using Malliavin
integration by parts to replace the explicit density-normalized stochastic-rate
correction already appears here. The paper also displays a domestic/foreign
correction with essentially the same outward form as ours. Neither idea should
be presented as our invention.

I did not find an equivalent of our four forward auxiliary SDEs, including
the state-dependent foreign-measure tangent correction, in this paper.
Its implemented two-rate method uses a locally frozen drift and fixed-strike
approximations. Its identification of the correction with an ordinary
conditional covariance also has a mathematical problem explained below.
Consequently, this paper does not establish that our exact integrated system
has already been derived. That is a limited comparison result, not a proof of
priority against the full literature.

The publication case survives as a specialized, incremental contribution:
an explicit finite-dimensional realization of established calibration ideas
for the coupled Gaussian-rate FX diffusion, together with a sound equivalence
proof and convincing numerical or theoretical consequences. The current work
still does not establish general nonlinear existence or numerical superiority.

## What is already present

| Item | Location in Bloch–Nakashima | Consequence for our claims |
|---|---|---|
| FX local volatility with stochastic domestic and foreign rates | Introduction; Section 5, p. 13 | The model setting and stochastic-rate Dupire identity are established background. |
| Malliavin integration by parts for the rate correction | Appendix B, pp. 30–32 | The basic Malliavin strategy predates our work. |
| Cancellation of the explicit strike-density factor | Appendix B, p. 32 | Removing the explicit density divisor is not a new conceptual contribution. |
| Conditional expectation involving an integrated first variation | Appendix B, p. 32 | Even the integrated-tangent interpretation has prior art in the one-rate, frozen-drift treatment. |
| Domestic-minus-foreign correction | Section 7.3, equation (7.10), p. 20 | The outward form of our main formula cannot establish novelty. |
| Fixed-point calibration and a quadratic bootstrap | Sections 8.1–8.2, pp. 20–22 | Fast calibration by an iterative or bootstrap correction is not itself new. |
| Explicit four-scalar stochastic system for the fully coupled diffusion | No equivalent located in the paper | This remains the main candidate contribution. |

Appendix B explicitly credits Benhamou, Gruz and Rivoira, *Stochastic interest
rates for local volatility hybrids models* (2008), listed as reference [BGR08]
on p. 36. Their original paper therefore remains another necessary comparison;
Bloch–Nakashima's presentation is not a substitute for reading it.

## Where the implemented method differs

Section 7 is explicitly an analysis of approximated dynamics, invoking the
locally constant drift of Section 2.1. Section 8.1 then replaces the dependence
on earlier random spot values in volatility integrals by values at a fixed
strike. The transition is explicit on p. 21: the integrals are approximated
using `sigma(u, K)` rather than the path-dependent `sigma(u, S_u)` expressions
displayed on p. 20.

Section 8.2 treats volatility as piecewise constant in time at that fixed strike
and obtains a quadratic equation. An exact solution of this quadratic is not
an exact solution of the original continuous-time smile calibration problem.
Section 9's numerical examples set the foreign rate to zero (p. 22), so those
examples do not validate a fully coupled two-stochastic-rate implementation.

This qualification must not be extended indiscriminately to the whole paper:
Section 6 does attempt continuous-dynamics calculations and contains volatility
derivatives. It is nevertheless not the four-process construction in our
manuscript, and its step from unconditional variance expressions to values at a
terminal strike does not supply a payoffwise calibration proof.

Our equations carry the actual pathwise coefficient
`g(t,S) = S partial_S sigma(t,S)`. In particular, they include

\[
\ell_t^d g(t,S_t)(dW_t^S-\sigma(t,S_t)dt),\qquad
\ell_t^f g(t,S_t)dW_t^S,
\]

and the cross-rate feedback terms

\[
-\rho_{Sf}\sigma_f g(t,S_t)\ell_t^d,\qquad
+\rho_{Sd}\sigma_d g(t,S_t)\ell_t^f.
\]

These arise from differentiating the coupled diffusion and the foreign
measure change, including

\[
H=w(\nabla\sigma)^\top,
\qquad w=(S\sigma,\rho_{Sd}\sigma_d,\rho_{Sf}\sigma_f)^\top.
\]

The numerical use of our exact representation will still introduce time-step,
regression, finite-particle and regularization errors. Exactness of the identity
must not be advertised as exactness of its numerical implementation. Also,
our explicit closure relies on Gaussian rates; we are not claiming a more
general short-rate specification than every model discussed in the older paper.

## Mathematical issue 1: the conditional covariance is not the Malliavin kernel

Equation (7.10), p. 20, uses terms of the form

\[
\operatorname{Cov}(r_T^i,Y(T,T)\mid X_T=\log K),
\]

where `X = log S`, `Y(t,T) = log F(t,T)`, and therefore `Y(T,T) = X_T`.
Under the ordinary definition of conditional covariance, for square-integrable
random variables,

\[
\operatorname{Cov}(Z,X\mid X)
=\mathbb E[ZX\mid X]-\mathbb E[Z\mid X]\mathbb E[X\mid X]
=0.
\]

Thus the displayed terms vanish almost surely in the conditioning variable
if read literally. The intended nonzero correction cannot be that ordinary
conditional covariance.

Appendix B, p. 32, makes the issue more than an isolated notation ambiguity:
it drops the conditional mean-product term and applies an unconditional
Malliavin duality identity inside a conditional expectation. Neither operation
is justified by the displayed argument. Conditioning a centered stochastic
integral on the terminal spot does not generally preserve its zero mean.

A simple illustration is `X = Z = W_T`. Then

\[
\operatorname{Cov}(X,Z\mid X)=0,
\qquad
\mathbb E\left[\int_0^T D_sX\,D_sZ\,ds\,\middle|\,X\right]=T.
\]

A charitable interpretation is that the authors intend a conditional
Malliavin covariance kernel. That object can indeed be nonzero and relevant.
But it must be defined through the appropriate duality, rather than identified
with ordinary conditional covariance. This distinction does not erase the
paper's conceptual priority.

## Mathematical issue 2: log-variance matching does not prove smile calibration

Section 7.2, pp. 19–20, writes the variance of the log spot as the expectation
of integrated instantaneous variance, then differentiates and substitutes a
terminal strike. This is not an exact identity for a general continuous-time
local-volatility diffusion.

For example, take zero rates, let `L = log S`, and write
`a(t,x) = sigma(t,exp(x))^2`. Then

\[
dL_t=-\tfrac12a(t,L_t)dt+\sqrt{a(t,L_t)}dW_t.
\]

Itô's formula, with the usual finite-moment assumptions, gives exactly

\[
\frac{d}{dt}\operatorname{Var}(L_t)
=\mathbb E[a(t,L_t)]-\operatorname{Cov}(L_t,a(t,L_t)).
\]

The additional covariance generally does not vanish. A smooth, uniformly
positive bounded example is

\[
a(x)=a_0+\varepsilon\tanh x,
\qquad a_0>\varepsilon>0.
\]

At every positive time the law is nondegenerate, and strict monotonicity gives
`Cov(L_t,a(L_t)) > 0`. Indeed, for an independent copy `L_t'`, this covariance
equals

\[
\tfrac12\mathbb E[(L_t-L_t')(a(L_t)-a(L_t'))]>0.
\]

This is an analytic counterexample to treating the simplified variance identity
as exact for the original diffusion. Section 7 declares an approximation, so
the conclusion is about its scope; it is not a claim that the authors never
mention approximations. Their abstract's broad exact-calibration wording is
nevertheless not supported by this argument.

Furthermore, an unconditional moment equality cannot simply be localized by
setting a terminal random variable equal to a strike. Matching a log variance
does not determine an entire terminal distribution. Our smooth-payoff duality
and forward-PDE uniqueness argument address this missing step, conditional on
existence and the stated regularity assumptions.

## Mathematical issue 3: the rate-to-drift sensitivity matters

Appendix B, p. 31, invokes adaptedness of a bond process to discard its
Malliavin derivative and concludes that the rate-driver sensitivity of the log
spot is proportional to the spot-driver sensitivity through their correlation.
Adaptedness only eliminates sensitivity to *future* noise; it does not eliminate
dependence on past noise. The separate fact that a bond pays one at maturity
does not eliminate the rate dependence of the rest of the spot dynamics.

In the original continuous model with one Hull–White rate, constant spot
volatility and zero spot/rate correlation, the domestic forward-measure
sensitivity is

\[
D_s^{W_r,T}\log S_T
=\sigma_r\frac{1-e^{-\kappa(T-s)}}{\kappa},\qquad s<T.
\]

It is nonzero because the log spot contains the integrated short rate.
Deterministic forward-measure drift shifts do not change this derivative.
At zero mean reversion the expression is `sigma_r (T-s)`. An identity setting
this entire sensitivity to zero at zero correlation can describe a removed or
frozen drift channel, but not the original continuous stochastic-rate model.

## A precise interpretation of our correction

Fix a maturity `T` and let `L = log S_T`. Under each rate's own `T`-forward
measure, write `c_i = r_T^i-r_i^0(T)`. Our smooth-payoff duality implies

\[
\mathbb E^{i^T}[c_i f(L)]
=\mathbb E^{i^T}[\ell_T^i f'(L)],\qquad i=d,f,
\]

for suitable smooth test functions. Consequently
`tau_i(x) = E^{i^T}[ell_T^i | L=x]` is a cross Stein kernel: a function
characterized by this integration-by-parts identity, not a conditional
covariance. This is an interpretation through a standard concept, not a claim
to a new Stein identity.

If `p_i` is the density of `L` in that measure and
`m_i(x) = E^{i^T}[c_i | L=x]`, the distributional identity is

\[
\partial_x(p_i\tau_i)=-p_i m_i.
\]

Under suitable decay at the upper boundary, it gives

\[
\tau_i(x)=\frac{\int_x^\infty m_i(y)p_i(y)\,dy}{p_i(x)}.
\]

This also explains why cancelling an explicit density divisor is not a theorem
eliminating tail estimation difficulties. The same conditional quantity still
contains that information. The foreign/domestic forward change of measure has
a terminal density proportional to `S_T`, so their conditional expectations of
the same integrable `ell_T^f` at a fixed spot agree.

In the jointly Gaussian, deterministic-spot-volatility case the kernel is the
constant unconditional covariance `Cov^{i^T}(c_i,L)`. Thus the ordinary Gaussian
covariance calculation is recovered, while `Cov(c_i,L | L)` is still zero.
This is why flat volatility is a valuable check and a weak novelty claim.
The paper's substantive test must involve a nonconstant spot volatility.

## Revised publication strategy

1. Cite Bloch–Nakashima and its attribution to Benhamou–Gruz–Rivoira prominently.
   Do not claim to originate density cancellation, Malliavin calibration or the
   domestic/foreign correction structure.
2. State the candidate contribution narrowly: four forward auxiliary processes
   realizing the exact coupled Gaussian-rate FX correction, with the foreign
   tangent transformation and a calibration equivalence conditional on existence.
3. Compare a nonflat synthetic recovery problem against the older fixed-strike
   approximation and a modern control-variate particle method. Measure repricing
   error, statistical error, runtime, derivative estimation and tail sample size.
   An exact representation may cost more and need not win every practical test.
4. For a stronger mathematical paper, prioritize a genuine extension beyond the
   representation: nonflat local existence, a justified moment-PDE formulation
   with stability, or a quantitative calibration-error bound with an estimable
   residual. These remain research goals, not consequences proved here.

The remaining priority checks include the original Benhamou–Gruz–Rivoira paper
and the relevant Guyon–Henry-Labordère derivations. This comparison establishes
neither first discovery nor publication readiness. It identifies more precisely
what would have to be new and useful.

## Scope of this review

The comparison inspected the supplied full text, particularly Sections 3 and
5–9 and Appendices A–B, and visually checked the main equations on the rendered
PDF pages. The counterexamples above are analytic derivations; no new Monte
Carlo experiment was run for this note. No claim is made to have re-audited every
formula in the 37-page paper or established novelty against all prior work.
