# Publication, novelty and research assessment

Recorded: 2026-11-14T10:03:00Z. Manuscript and code assessed at local commit
`b6c12c9`. Literature checked on 14 November 2026.

**Follow-up, 14 November 2026:** The complete Bloch–Nakashima (2008) paper has
now been supplied and compared at equation level. The
[dedicated comparison](2026-11-14_bloch_nakashima_prior_art_comparison.md)
supersedes the full-text access limitation below. It establishes substantial
prior art for the Malliavin/density-cancellation idea, while finding no equivalent
of our fully coupled four-process system in that paper. It does not establish
priority against the remaining literature.

## Judgment

The work is worth developing into a focused computational/mathematical-finance
paper. The current evidence does not justify calling it publication-ready, a
solution of general stochastic-rate calibration, or an established improvement
over the best available numerical methods.

The most defensible candidate contribution is the explicit, exact four-scalar
integrated-sensitivity representation for three-factor FX local volatility with
Gaussian domestic and foreign rates, with the foreign change-of-measure tangent
correction correctly included. This is a potentially useful, specialized
extension of existing Malliavin/particle calibration ideas. Its priority is not
yet established by a complete equation-level comparison with the closest older
literature. No acceptance probability or claim of first discovery is warranted.

## Novelty assessment by component

| Component | Assessment |
|---|---|
| Stochastic-rate Dupire calibration identity | Established background; not our novelty. |
| McKean particle calibration and conditional expectation estimation | Established background. |
| Malliavin representation for hybrid calibration | Already part of Guyon–Henry-Labordère's programme. |
| Explicit domestic/foreign integrated tangent system and four-scalar reduction | Main candidate for incremental novelty; requires close prior-art comparison. |
| Rank-one foreign tangent correction | Essential to this derivation, but follows standard differentiation of a Girsanov-adjusted drift; not a new change-of-measure theorem. |
| Complete forward-PDE sufficiency proof | Necessary mathematical completion; standard uniqueness argument, not the main original contribution. |
| Local flow/integrability hypotheses, singular correlations and time-dependent Gaussian parameters | Useful scope improvements; mainly extensions of the representation rather than a new existence theory. |
| Flat-smile existence | Exact benchmark with term structure; modest independent novelty and no general-smile existence conclusion. |
| Fully regularized strong existence | Application of familiar Lipschitz McKean–Vlasov arguments; useful specification, weak standalone novelty. |
| Log-contract compatibility bound | Useful supporting result with broader scope than the Gaussian subclass; related variance-floor ideas already exist. |
| Stable covariance evaluation and QR sampling | Valuable numerical engineering using standard techniques. |
| Lower estimator standard deviation in the existing experiment | Promising numerical evidence; not yet an accuracy-per-runtime or calibration result. |

My earlier additions should not be used to inflate the paper's originality.
Correct assumptions, a complete proof, and robust covariance calculations make
the work credible; they do not by themselves make it new.

## Closest literature and what it changes

1. [Guyon and Henry-Labordère, *The Smile Calibration Problem Solved* (2011)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1885032)
   already describes particle calibration and Malliavin representations for
   hybrid models. The 2013 book is central prior art. The preprint abstract was
   retrieved, but full preprint/book equations were not accessible in this
   search. An equation-by-equation comparison remains a publication prerequisite.

2. [Deelstra and Rayée, *Local Volatility Pricing Models for Long-dated FX Derivatives*](https://arxiv.org/pdf/1204.0633)
   develops the three-factor stochastic-rate FX setting and calibration
   identities. Its discussion also identifies Malliavin integration by parts
   as a way to estimate conditional expectations. Relevant full-text sections
   were inspected. Merely combining FX, two rates and local volatility is not new.

3. [Bloch and Nakashima, *Multi-Currency Local Volatility Model* (2008)](https://www.iupindia.in/1008/IJDM_Multi_Currency_35.html)
   is important missing prior art in our present bibliography. Its publisher
   description discusses stochastic domestic/foreign rates, variance-covariance
   corrections and calibration; the SSRN record also lists Malliavin calculus.
   The full paper was not retrieved. The description alone cannot establish
   whether its formulas coincide with, approximate, or differ from ours.
   [SSRN record](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1153337).

4. [Cozma, Mariapragassam and Reisinger, control-variate particle calibration (2019)](https://arxiv.org/pdf/1701.06001)
   supplies a stronger numerical comparator than a plain tail estimator. Its
   four-factor FX study includes calibration, repricing and variance reduction,
   with a framework adaptable beyond the particular CIR specification. Relevant
   full-text sections were inspected. Our 14–66 standard-deviation ratios cannot
   establish superiority over that method or be compared directly with its
   differently defined variance-reduction figures.

5. [Ögetbil, Ganesan and Hientzsch (2022)](https://arxiv.org/html/2009.14764v2)
   treats FX calibration with Gaussian rates and gives simulation algorithms,
   calibration tests and a Gaussian total-implied-variance lower bound in
   Section 5. That makes generic variance-floor claims unoriginal. Our
   log-contract bound concerns an entire terminal law and minimizes over
   nonnegative, possibly state-dependent volatility; it must not be equated
   with a bound on Black implied variance at each individual strike.

6. [Joseph, Loeper and Obłój (2026)](https://link.springer.com/article/10.1007/s00780-026-00588-6)
   provides an optimal-transport calibration method with stochastic rates and
   discounted-density duality. Its example permits volatility to depend on the
   short rate and calibrates discrete option constraints. This is a different
   admissible model class and optimization problem, not a direct replacement
   for our spot-only conditional representation. Relevant full text was inspected.

7. [Reisinger and Tsianni](https://arxiv.org/abs/2302.00434) and
   [Bayer, Belomestny, Butkovsky and Schoenmakers (2024)](https://link.springer.com/article/10.1007/s00780-024-00541-5)
   already establish substantial regularized LSV theory and numerical analysis.
   The former gives explicit regularization-dependent Euler error estimates;
   the latter develops RKHS regularization and propagation of chaos. A generic
   clipped/smoothed well-posedness proposition is below that existing level of
   contribution. The former abstract and the latter relevant full-text sections
   were inspected.

8. [Friz, Jourdain, Wagenhofer and Zhou, *On the Weak Error for Local Stochastic Volatility Models*](https://arxiv.org/pdf/2506.10817)
   (2025 preprint; [2026 journal publication](https://epubs.siam.org/doi/10.1137/25M1768515))
   shows why general continuous-time existence need not be the first obstacle
   tackled: it analyzes a defined time-discrete calibration construction and
   its weak error. This supplies a useful research strategy, not a theorem
   automatically applicable to our discounted-rate/sensitivity system.

9. [Djete](https://arxiv.org/html/2208.09986v2),
   [Mustapha](https://arxiv.org/abs/2406.14074), and
   [Arnese](https://arxiv.org/abs/2507.22222) clarify the scope of current
   existence/convergence theory. The inspected Djete theorem has a particular
   positive density-dependent structure; Mustapha's abstract assumes a
   finite-valued volatility factor; Arnese's abstract concerns conditional
   expectations in the drift. None supplies the missing general theorem here
   without substantial new verification. These are not interchangeable notions
   of a conditional McKean–Vlasov equation.

Searches found no exact counterpart to our four displayed auxiliary equations
in the accessible material inspected. This is a bounded literature finding,
not proof of novelty. The inaccessible closest sources matter more than a large
number of loosely relevant search hits.

## Practical usefulness and the decisive missing experiment

The representation replaces rate-tail expectations and their explicit market
density divisor with regression of four forward-evolving sensitivity variables.
It avoids storing the entire Malliavin derivative history. It remains nonlinear,
requires conditional estimation, and requires the spatial derivative of the
calibrated coefficient. Existing tail-based particle methods also evolve
forward without an ever-growing derivative history, so this observation alone
does not imply a complexity advantage over those competitors.

The reproduced 14–66 reduction is in standard deviation in one prescribed-volatility
experiment. It does not measure total error or speed. Smoothing bias, unsupported
strikes, regression derivatives, finite-particle dependence and calibration
feedback may change the result. There is no automatic Rao–Blackwell variance
ordering between these two different integration-by-parts estimators.

The derivative is a particularly consequential cost: for ordinary independent
samples and a smooth one-dimensional kernel, the elementary kernel scaling gives
variance of order 1/(N h) for a local mean and 1/(N h^3) for its first derivative,
under regular positive-density conditions. This is a motivation for careful
analysis, not a proven rate for our interacting particle algorithm.

The decisive experiment should use a nonconstant synthetic volatility with
skew and curvature, generate an independent accurate call surface, calibrate
without revealing the generating coefficient, and reprice with independent
paths. Compare a conventional tail estimator, an appropriate control-variate
baseline and our estimator at comparable total runtime. Vary bandwidth, time
step and particle count; record price/IV error, floor activation, derivative
error, conditional effective sample size, and repeated-seed variability.

Surface recovery is diagnostic when comparing against a known generating
solution. It is not an identification theorem: without uniqueness, excellent
call-price recovery need not force convergence to that exact generating surface.
For a desk, lower calibration and Greek instability could be valuable even
without the lowest wall-clock runtime, but that benefit also needs evidence.

## Two concrete mathematical consequences worth developing

These observations follow from our representation and standard stochastic/PDE
calculus. They are research leads, not claims of new general existence theorems.

### A. A quantitative calibration residual bound

For an existing regular coefficient sigma, define the true-law residual

\[
R(t,K)=\sigma^2(t,K)+2E^{d^t}[\ell_t^d-\ell_t^f\mid S_t=K]
       -\sigma_{Dup}^2(t,K).
\]

The sufficiency proof already yields, for e=C_model-C_market,

\[
\partial_t e=\mathcal L_{Dup}e+
  \tfrac12 K^2 R(t,K)\partial_{KK}C_{model}(t,K),
\]

where the zero-order coefficient of the market forward operator is -r_f^0(t).
Under a valid comparison/evolution-operator principle, zero initial and boundary
differences, and integrability of the forcing in the supremum norm, Duhamel's
formula gives

\[
\|e(t,\cdot)\|_\infty\leq\frac12\int_0^t
 e^{-\int_s^t r_f^0(u)du}
 \sup_{K>0}|K^2R(s,K)\partial_{KK}C_{model}(s,K)|\,ds.
\]

This is a direct quantitative extension of sufficiency: it connects a failure
of the calibration identity to vanilla price error. Its PDE mechanism is
standard; the valuable next work is a useful norm, a verifiable estimator of
the true residual, statistical error control and treatment of discrete auxiliary
equations. Substituting the empirical residual into this bound without those
steps is not justified. The current finite-particle auxiliaries are not exact
Malliavin derivatives of the interacting system.

### B. A closed system for discounted first moments

Let Z=(log S,x_d,x_f) be the physical three-dimensional state and
Y=(ell_d,zeta_f,ell_f,zeta_d). For a prescribed deterministic sigma, write

\[
dZ=b(t,Z)dt+\sum_\alpha V_\alpha(t,Z)dW^\alpha,\qquad
dY=(A(t,Z)Y+f(t,Z))dt+\sum_\alpha G_\alpha(t,Z)YdW^\alpha,
\]

using independent Brownian coordinates. The auxiliary equations are affine in Y.
Define discounted measures p=E[D delta_Z] and m_i=E[D Y_i delta_Z].
Applying Ito to D phi(Z) and D Y_i phi(Z) shows, in weak form and assuming the
necessary integrability,

\[
\partial_t p=\mathcal L_Z^*p-r_d p,
\]

\[
\partial_t m_i=\mathcal L_Z^*m_i-r_d m_i+
 \sum_j A_{ij}m_j+f_i p
 -\sum_\alpha\nabla_Z\!\cdot\!\left[V_\alpha(G_\alpha m)_i\right].
\]

The cross-divergence is the quadratic-covariation term. No second moments of Y
are required. The calibration correction is obtained by integrating m_1-m_3
over the two rate coordinates and dividing by the corresponding integral of p.
Thus the relevant information can be represented by five fields over the
three physical coordinates, rather than a full density over all augmented
coordinates. For fixed sigma these are linear equations; imposing calibration
makes their coefficients nonlinear functions of the fields.

This could provide an independent deterministic check, a route to regularity
estimates, and a framework for local non-flat existence. It is not a one-dimensional
closure, and it does not by itself solve positivity, boundary behavior, or
well-posedness of the resulting nonlinear PDE. First-moment closure for affine
auxiliaries is a standard mechanism; novelty would lie in this application and
what can be proved or computed with it.

## Further directions and relative priority

1. **Calibration error and stability:** develop A into a theorem that includes
   regression, auxiliary discretization and particle errors. Most direct link
   between the mathematical representation and practical value.
2. **Non-flat local existence:** choose an admissible Gaussian reference with
   positive variance margin, perturb the target with a smooth genuine smile,
   and prove existence/stability in an explicitly specified function class.
   A contraction or implicit-function argument requires estimates for the law
   map and its derivative; spatial derivative loss, the time-zero Dirac law,
   density tails and positivity are substantive obstacles. This is a potentially
   stronger theorem, not a routine consequence of the flat case.
3. **Finite-dimensional kernel closure:** generalize the integrated directions
   to deterministic rate kernels with a finite-dimensional realization, including
   multiple Gaussian rate factors. Derive the required cross sensitivities and
   dimension rather than asserting two states per factor. The general mechanism
   is familiar; a useful closure/approximation theorem could broaden the paper.
4. **Sharper compatibility results:** investigate equality and attainability of
   the log-contract bound and consistency across maturities. An aggregate
   necessary condition is not sufficient for fitting a whole smile. Standard
   square-completion alone is unlikely to justify a separate paper.
5. **A genuine LSV extension:** adding a stochastic variance factor could increase
   practical relevance, but introduces additional tangent channels and
   regularity issues. Do not claim the current four-process system already
   solves that model.

## Publication strategy

A defensible central claim is: an explicit finite-dimensional Malliavin
representation for Gaussian-rate FX local volatility, with a proved conditional
calibration equivalence and reproducible numerical evaluation. Avoid claiming
the first stochastic-rate calibration algorithm, universal tail stability,
general LSV existence, or an established speedup.

For an applied paper, the non-flat recovery experiment and a strong comparator
are the next priorities after the missing prior-art comparison. An exact
representation plus clear, reproducible practical benefit can be publishable
without solving the full general existence problem. If the numerical advantage
disappears, an independently novel generalization or error theorem would become
more important. If older work contains an equivalent closed system, the present
representation alone would not support a strong novelty claim.

The [Journal of Computational Finance](https://www.risk.net/journal-of-computational-finance)
is a natural scope match for a completed computational study.
[Quantitative Finance](https://www.tandfonline.com/journals/rquf20/about-this-journal)
is another plausible venue if the contribution and applications are substantial.
These are scope judgments, not predictions of acceptance. A stronger
mathematical-finance submission would benefit from a substantive stability,
error, closure, or non-flat existence theorem.

Keep the flat case, routine regularized existence and covariance implementation
as supporting material. Center the manuscript on the four-process representation
and the best verified consequence of that representation.
