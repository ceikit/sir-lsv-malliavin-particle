# Review of the SIR–LV Malliavin particle repository and supplied feedback

Recorded: 2026-11-14T13:22:27+04:00 (Asia/Dubai); 2026-11-14T09:22:27Z (UTC).

Reviewed 14 November 2026. Repository: https://github.com/ceikit/sir-lsv-malliavin-particle

Snapshot: `717f936e063c9e1289d2138df6f5b1d6be513ac7`.

## Assessment

The central three-factor Hull–White/local-volatility representation is mathematically convincing under appropriate assumptions. I found no substantive algebraic error in the integrated sensitivity equations, the foreign Girsanov correction, or the four-scalar reduction. Independent directional perturbations, a Gaussian benchmark, and reruns of the supplied Monte Carlo scripts support that assessment.

The implementation is a collection of research validation scripts, not an implemented calibration engine. The impressive variance reduction is reproducible at prescribed local volatility, but calibration stability, bias, convergence, and performance relative to a strong control-variate baseline remain unestablished.

There are actionable errors and omissions: a sign error and stale algorithm timing in the older one-factor section; unstable floating-point covariance formulas; insufficiently qualified general-rate claims; incomplete handling of market-data statistics; and a missing distinction between the intended spot-option calibration surface and the futures options downloaded by the data script.

The supplied feedback is substantially correct. Some of it is already addressed by the current repository, and it should not be treated as an exhaustive code audit.

## 1. Independent mathematical assessment

### 1.1 Calibration identity and signs

Write

\[
c_t^i=r_t^i-f^i(0,t),\qquad D_t=\exp\left(-\int_0^t r_s^d ds\right).
\]

For the model call surface \(C_m(t,K)=\mathbb E^d[D_t(S_t-K)^+]\), Itô–Tanaka gives

\[
\partial_t C_m+(f^d-f^f)K\partial_K C_m+f^f C_m
=\frac12K^2\sigma^2(t,K)\partial_{KK}C_m
+\mathbb E^d[D_t(Kc_t^d-S_tc_t^f)\mathbf1_{S_t>K}].
\]

The single-measure correction in the paper has the correct sign, discount factor, strike powers, and foreign spot factor.

### 1.2 The foreign-measure correction is correct

For \(X=(S,r^d,r^f)\), let \(b\) be the domestic drift and

\[
w=(S\sigma,\rho_{Sd}\sigma_d,\rho_{Sf}\sigma_f)^\top.
\]

Changing to the foreign risk-neutral measure adds \(\sigma w\) to the drift. Differentiating that drift and converting the tangent noise back to domestic Brownian motions leaves

\[
\nabla(\sigma w)-\sigma\nabla w=w(\nabla\sigma)^\top=H.
\]

Thus the additional matrix in the foreign variational equation is necessary. It is not legitimate simply to identify the Malliavin derivatives under the two stochastic Girsanov transformations.

The logarithmic sensitivities then satisfy exactly the displayed four-process system. In particular, the domestic multiplicative term is \(g(dW-\sigma dt)\), while the foreign term is \(g\,dW\). The foreign rate’s own quanto derivative cancels. The two remaining cross-rate feedback terms also have the correct signs.

Source: [integrated sensitivities and reduction](https://github.com/ceikit/sir-lsv-malliavin-particle/blob/717f936e063c9e1289d2138df6f5b1d6be513ac7/paper/SIR_LSV_particle_FX.tex#L630-L679).

### 1.3 A direct sufficiency proof should be added

Define the model discounted density

\[
m(t,K)=\mathbb E^d[D_t\delta(S_t-K)]=\partial_{KK}C_m(t,K).
\]

Use this model density first, without identifying it with the market density. Malliavin integration by parts gives

\[
\mathbb E^d[D_tc_t^d\mathbf1_{S_t>K}]
=K m(t,K)\mathbb E^{d^t}[\ell_t^d\mid S_t=K],
\]

\[
\mathbb E^d[D_tS_tc_t^f\mathbf1_{S_t>K}]
=K^2m(t,K)\mathbb E^{d^t}[\ell_t^f\mid S_t=K].
\]

The second equality uses the foreign forward measure and

\[
\frac{d\mathbb Q^{f^t}}{d\mathbb Q^{d^t}}\bigg|_{\mathcal F_t}=\frac{S_t}{F_{0t}}.
\]

Conditioning on \(S_t=K\) therefore gives the same conditional expectation under the two forward measures. Substitution yields

\[
\partial_t C_m+(f^d-f^f)K\partial_K C_m+f^f C_m
=\frac12K^2\left(\sigma^2+2\mathbb E^{d^t}[\ell^d-\ell^f\mid S_t=K]\right)\partial_{KK}C_m.
\]

If the proposed calibration identity holds, the coefficient in parentheses is the market Dupire variance. Both market and model calls then solve the same linear forward PDE. State uniqueness in an appropriate solution class, with the common initial condition \((S_0-K)^+\), boundary behavior \(C(t,0)=S_0P^f_{0t}\), decay as \(K\to\infty\), and the required growth/integrability conditions. This establishes sufficiency conditional on existence and regularity of the nonlinear model.

It does not establish existence or uniqueness of that nonlinear model, nor convergence of its particle approximation.

### 1.4 The general-rate formula is broader than its stated hypotheses justify

At [line 553](https://github.com/ceikit/sir-lsv-malliavin-particle/blob/717f936e063c9e1289d2138df6f5b1d6be513ac7/paper/SIR_LSV_particle_FX.tex#L553), the centered rate is represented by one Clark–Ocone integral against its named rate driver. This is valid in the autonomous Hull–White setting, but not for every unspecified adapted drift and diffusion allowed by the initial general SDE.

For example, with independent Brownian drivers, a centered Gaussian variable

\[
R_t=B_t^d+\int_0^t B_u^f du
=\int_0^t dB_s^d+\int_0^t(t-s)dB_s^f
\]

requires both drivers. Its directional derivative in the domestic rate direction alone does not reconstruct the second term. Merely writing an SDE with diffusion \(dB^d\) does not ensure that its adapted drift contains no other noise.

Either restrict the general statement to an autonomous one-factor rate under its own forward measure, or use the full vector martingale representation. This does not invalidate the later Hull–White result.

### 1.5 Finite particles and unknown local volatility

The auxiliary equations differentiate a diffusion with a fixed deterministic coefficient function \(\sigma(t,S)\), including a deterministic function obtained from a limiting law. At finite particle count, \(\sigma_N\) depends on the full simulated population. Exact differentiation of that interacting system would include the empirical estimator’s dependence on other particles and their weights.

The proposed scheme should therefore be described as an approximation to the limiting augmented nonlinear system. Its auxiliary variables are not exact Malliavin derivatives of the finite-N interacting system. Their convergence cannot be justified merely by declaring omitted terms to be small: bandwidth dependence also matters.

### 1.6 Negative estimates do not automatically establish nonexistence

The Gaussian benchmark provides a useful compatibility test in its stated setting. In general, the correction depends on the unknown volatility function, its spatial derivative, and the generated law. A negative trial update does not prove that every candidate solution is impossible.

Similarly, a finite-sample floor can bind because of estimator noise, time discretization, interpolation, or genuine incompatibility. The current instruction that the remedy is necessarily in model inputs is too strong. Report floor locations, frequency and magnitude; diagnose numerical convergence before changing economic inputs.

### 1.7 Smaller manuscript issues

- **One-factor forward sign:** [line 200](https://github.com/ceikit/sir-lsv-malliavin-particle/blob/717f936e063c9e1289d2138df6f5b1d6be513ac7/paper/SIR_LSV_particle_FX.tex#L198-L201) has a minus sign before \(\sigma_rA\,dB^T\). Since bond volatility is \(-\sigma_rA\), the forward has **plus** \(\sigma_rA\,dB^T\). The subsequent exponential discretization and three-factor implementation use the correct sign.
- **One-factor timing:** [lines 295–301](https://github.com/ceikit/sir-lsv-malliavin-particle/blob/717f936e063c9e1289d2138df6f5b1d6be513ac7/paper/SIR_LSV_particle_FX.tex#L295-L301) advance particles to \(t_{j+1}\) but assign the recalibrated coefficient to the preceding interval. Use the corrected three-factor algorithm’s indexing: estimate at \(t_{j+1}\), then apply on the next step.
- **Multi-factor outlook:** “each factor contributes its own pair” is not a general state-count formula. With several domestic and foreign mean-reversion factors, cross sensitivities generally need to be tracked separately; different decay rates prevent arbitrary aggregation into one scalar. The matrix construction extends, but derive the reduction before promising two scalars per factor.
- **Title and scope:** the new developed result is three-factor FX local volatility with two Hull–White rates. A general LSV extension is only sketched. “Exact” should refer to the continuous representation, not finite-N numerical calibration.
- Put differentiability, flow, integrability, true-martingale numeraire changes, density positivity and PDE uniqueness hypotheses before the propositions. The caveats section should not carry the entire theorem specification.

## 2. Code audit and reproduced findings

### 2.1 The supplied dynamics are implemented consistently

The main sensitivity updates use the same Brownian increments as the state. Right-hand sides are evaluated at old states before replacement. The foreign quanto and forward-measure drift adjustments match the mathematics. Antithetic partners are kept together when estimating standard errors and when forming estimator-comparison batches.

The constant-volatility test and independent directional checks below are additional evidence beyond two implementations of the same displayed formulas agreeing with each other.

### 2.2 Reproducible covariance/Cholesky failure

Location: [coefficient formulas](https://github.com/ceikit/sir-lsv-malliavin-particle/blob/717f936e063c9e1289d2138df6f5b1d6be513ac7/checks/check_forward_measure.py#L29-L45) and [Cholesky call](https://github.com/ceikit/sir-lsv-malliavin-particle/blob/717f936e063c9e1289d2138df6f5b1d6be513ac7/checks/check_forward_measure.py#L95-L100).

The analytic formulas are correct but contain severe cancellation: differences of exponentials are subtracted again and divided by small powers of kappa and the time step. The true covariance is a Gram integral and is positive semidefinite; floating-point evaluation can destroy that property.

Using the repository’s correlations, terminal horizon 6 years and step `[0, 0.01]`:

| Mean reversions | Minimum eigenvalue of implemented covariance | Cholesky |
|---|---:|---|
| 0.05, 0.10 | 2.0069964e-7 | succeeds |
| 0.0001, 0.0002 | -2.0336868e-5 | fails |

For the latter case, an independently assembled covariance using positive Gauss–Legendre weights and stable kernel evaluation has minimum eigenvalue **1.3547275e-7** and Cholesky succeeds. This confirms numerical cancellation rather than an invalid correlation input.

With the original mean reversions, failures also appear at some tested times with 5,000 or 10,000 steps per year. Ordinary 100–200-step runs pass.

Reproduction:

```python
import sys
import numpy as np
sys.path.insert(0, "checks")
import check_forward_measure as f
f.kd, f.kf = 1e-4, 2e-4
cov = f.Sigma(f.coeffs(0.0, 0.01))
print(np.linalg.eigvalsh(cov)[0])
np.linalg.cholesky(cov)
```

Use `expm1`, small-argument series for composite formulas, and explicit zero-mean-reversion limits. A covariance constructed as a positive weighted Gram sum is a useful independent reference or fallback. Simply clipping a materially negative eigenvalue would conceal the source of the problem.

### 2.3 Estimator variance improvement reproduced, but bias is unresolved

I reran the advertised configuration: **2,000,000 paths**, 100 time steps/year, 5 years, 100 antithetic batches of 20,000 particles.

| Spot quantile | Conventional SD / Malliavin SD |
|---:|---:|
| 1% | 66.3 |
| 5% | 35.0 |
| 10% | 23.0 |
| 25% | 16.9 |
| 50% | 42.1 |
| 75% | 14.2 |
| 90% | 17.0 |
| 95% | 19.1 |
| 99% | 29.8 |

This independently reproduces the README’s 14–66 range for this benchmark. It is an empirical standard-deviation reduction at a prescribed volatility, not an established speedup or an RMSE comparison.

Important limits of [compare_estimators.py](https://github.com/ceikit/sir-lsv-malliavin-particle/blob/717f936e063c9e1289d2138df6f5b1d6be513ac7/checks/compare_estimators.py#L31-L86):

1. Strikes, the density normalizer and the reference are estimated from the same full sample that supplies the batches. The “known” density is held fixed across batches; it is not an independent exact oracle.
2. The reference is another kernel estimate, with bandwidth 0.015; batch estimates use approximately 0.056. Neither is an exact target.
3. At the 5% quantile, the two Malliavin means differ by 0.218 in the displayed variance-times-1e4 units, versus a per-batch SD of 0.173. This does not measure true bias, because the reference is noisy, but it shows that bandwidth dependence is already comparable to the reported noise.
4. No timing comparison, calibration feedback, derivative estimation from noisy surfaces, or control-variate competitor is tested.
5. The final strike identities replace the Dirac delta by a bandwidth-0.01 kernel. Their z-scores include sampling error, but not a correction for smoothing bias or data-dependent strikes.
6. The compact-support regression has no guard against zero denominator at unsupported strikes.

Use an independent reference, fixed strikes, several bandwidths, particle counts and time grids. Report RMSE and repricing error at equal runtime. Report effective support/ESS and preserve antithetic pairing in uncertainty estimates.

### 2.4 Data handling and instrument mapping

The download script’s `build_table` uses unconditional “last message wins” and does not process `update_action`. Databento defines additions and deletions, and separate flags distinguish final/preliminary and actual/theoretical settlements. A synthetic add followed by delete for the same instrument/date still emerges from `build_table` with a settlement price. This was reproduced without calling the paid API.

Fix: replay the applicable statistics update semantics, then apply an explicit settlement-quality policy. Preserve enough provenance to diagnose corrections and instrument mapping. Source: [Databento statistics documentation](https://databento.com/docs/schemas/statistics), and [repository lines 35–45](https://github.com/ceikit/sir-lsv-malliavin-particle/blob/717f936e063c9e1289d2138df6f5b1d6be513ac7/data/pull_eurusd_options_databento.py#L35-L45).

A separate modeling issue: `EUU.OPT` supplies options on EUR/USD futures. CME’s [2026 product guide](https://www.cmegroup.com/markets/fx/fx-product-guide) states that these options deliver into underlying futures. The paper defines the calibration input as a discounted spot payoff \((S_t-K)^+\). Under stochastic rates, a futures option cannot simply be relabeled as a spot call at the same strike and expiry. Preserve the underlying contract and maturity and address exercise/delivery, premium conventions and futures/forward convexity, or use an appropriately constructed OTC spot-FX option surface.

The repository does not yet perform this mistaken calibration; this is a blocker for using its downloaded CSV as the paper’s input without further work. I did not run the paid download.

### 2.5 Validation infrastructure and missing engine

The checks print diagnostics but do not fail on a statistical/algebraic mismatch. They are useful experiments, not a regression test suite. Add deterministic limit-case assertions, covariance tests across a parameter grid, meaningful stochastic tolerances and recorded seeds. Do not require every stochastic diagnostic to fall within one standard error; that would create a flaky suite.

No implemented module performs the paper’s complete interacting-particle calibration, smooth derivative estimation, maturity-by-maturity surface update and independent vanilla repricing. The synthetic recovery experiment is the most consequential next implementation.

Python byte compilation passed for `checks/` and `data/`. I attempted the documented LaTeX build, but the local TeX installation lacks `pdflatex.fmt`; compilation was not verified. That is an environment limitation, not an identified source error. No repository source was changed or pushed.

## 3. Numerical evidence from this audit

Environment: NumPy 2.3.5, SciPy 1.17.0. BLAS threads limited to one for the simulations.

| Check | Configuration | Observation |
|---|---|---|
| Domestic-measure IBP identities | 200,000 paths; 100 steps/year; T=5; default skew/parameters | max absolute z-score 1.20 over tested identities/strikes |
| Forward-measure identities | 200,000 paths; 100 steps/year; t=4, T=6 | max absolute sensitivity z-score 0.99; tested martingale/curve identities within MC uncertainty |
| Forward Gaussian step coefficients | repository quadrature check | max relative error 3.8e-11 at its test interval |
| Stress IBP identities | 400,000 paths; T=10; skew 0.8; doubled rate vols | max absolute z-score 0.98; naive foreign formula deviates up to 3.44 standard errors |
| Matrix versus scalar forms | 200,000 paths; 100 then 400 steps/year | domestic relative RMS 0.00140 then 0.000648; foreign 0.00139 then 0.000618 |
| Constant-vol ODE versus analytic correction | default rate parameters, spot vol 0.12; T=1,5,10,30 | absolute differences below 1e-13 |
| Independent Cameron–Martin perturbations | 512 paths; T=5; central bump 1e-4; 50/200/800 steps/year | error against paper’s Euler weights decreases approximately linearly with dt |
| Full estimator benchmark | 2 million paths, 20,000-particle batches | SD ratio 14.2–66.3; max absolute strike-identity z-score 1.51 |

For the directional experiment, I perturbed the domestic driver in its own representation, and independently perturbed the foreign-measure driver after reconstructing its Brownian increments from the baseline domestic path. Correlated driver shifts were included. This tests the interpretation of the sensitivities directly, not only agreement between two algebraically equivalent auxiliary systems.

| Steps/year | Domestic relative RMS difference | Foreign relative RMS difference |
|---:|---:|---:|
| 50 | 8.61697e-4 | 1.62171e-3 |
| 200 | 2.18418e-4 | 4.10710e-4 |
| 800 | 5.41167e-5 | 1.01631e-4 |

The finite-step mismatch is expected: bumping the discretized state and Euler-discretizing the continuous integrated sensitivity are not identical operations. These results support consistency; they are not a convergence proof for calibration.

The constant-volatility analytical comparator used here was

\[
2(\ell_t^d-\ell_t^f)=\sigma_d^2A_d(t)^2+\sigma_f^2A_f(t)^2
+2\rho_{Sd}\sigma\sigma_dA_d(t)-2\rho_{Sf}\sigma\sigma_fA_f(t)
-2\rho_{df}\sigma_d\sigma_fA_d(t)A_f(t),
\]

where \(A_i(t)=(1-e^{-\kappa_i t})/\kappa_i\). This is the maturity derivative of the integrated rate contribution for constant spot volatility.

## 4. Evaluation of the supplied feedback

| Feedback item | Assessment against the reviewed repository |
|---|---|
| Main representation, H correction and four-process algebra are sound | Agree, subject to explicit analytic hypotheses. My independent checks support it. |
| A. Add a complete sufficiency proof | Agree. Use the model density before imposing calibration, then invoke forward-PDE uniqueness with matching initial/boundary conditions. |
| B. Wing stability is overstated | Correct warning, but partly already addressed: the current three-factor section explicitly says that Dupire still needs C_KK and regression still needs particles near the strike. Early one-factor language remains stronger. |
| C. Weights are not exact finite-N interacting-system derivatives | Agree. This is a substantive qualification for a future convergence theorem, not an algebraic refutation. |
| D. Restrict the negative-variance nonexistence inference | Agree. Also distinguish genuine incompatibility from finite-sample floor activation. |
| E. Promote hypotheses into theorem statements | Agree. Include the uniqueness assumption needed for the calibration equivalence. |
| Synthetic smile recovery before market calibration | Strongly agree. Reprice independently; do not test only the recovered volatility grid. |
| Scripts unavailable | No longer applies to this repository: all named scripts are present. I reproduced the full two-million-path estimator experiment. |
| Title should emphasize LV and the finite-dimensional representation | Agree. A fully developed stochastic-volatility extension is not supplied. |
| Novelty and numerical superiority remain unproved | Agree. The present review is not a systematic novelty search. |

### 4.1 The feedback’s Gaussian table is numerically credible but incomplete

The feedback omits the constant spot volatility needed to reproduce its table. With its stated parameters, all four values match **sigma = 0.14**:

| Maturity | Analytic value independently recomputed |
|---:|---:|
| 1 | 0.002617005525604 |
| 3 | 0.008378136940401 |
| 7 | 0.019834709486567 |
| 12 | 0.031647702946037 |

Thus the numbers are consistent, and the missing volatility can be inferred. The text also loses the precision thresholds and the strike in its skewed Monte Carlo paragraph. Restore these along with the seed and runnable code. I cannot certify that reviewer’s claimed execution from the pasted text; I can independently corroborate the algebra and Gaussian values and compare its stochastic magnitudes with the repository reruns.

### 4.2 Two refinements to its recommendations

First, ESS is a useful support diagnostic, but \((\sum w)^2/\sum w^2\) is not the exact number of independent samples in an antithetic simulation. Pairwise variance estimation still matters.

Second, recovery of the generating sigma is a valuable synthetic test, but convergence to that particular function requires uniqueness/identifiability assumptions that the manuscript has not established. Independent vanilla repricing is the essential calibration outcome; sigma recovery adds evidence about the inverse problem.

### 4.3 What the feedback misses

It does not identify the older one-factor sign/indexing mistakes, the reproducible covariance numerical failure, the too-general single-driver Clark–Ocone statement, the data-script update semantics, or the futures-versus-spot calibration-input issue. Its positive mathematical assessment is justified; its coverage is not a complete software review.

## 5. Recommended next work

1. Correct the localized manuscript mistakes and add the PDE sufficiency argument with explicit hypotheses. State the finite-N interpretation precisely.
2. Stabilize the Gaussian covariance implementation and add deterministic tests for small kappa, short time steps, zero rate volatility and constant spot volatility.
3. Implement one reusable augmented-state simulator and the actual calibration loop. Use smooth, controlled spatial derivatives; monitor support/ESS, tails, floors and derivative amplification.
4. Run a synthetic nonconstant-volatility recovery experiment. Generate the target with independent high-accuracy data or a suitable PDE calculation. Freeze the calibrated surface and reprice on new paths. Vary particle count, bandwidth, strike grid and time step separately.
5. Compare bias, RMSE and runtime with the conventional estimator using effective control variates, not only its unadjusted version. Published work already provides market-calibrated control-variate and stochastic-rate Monte Carlo baselines.
6. Only then map real quotes to the exact payoff and quote conventions used by the model and evaluate market calibration.

The useful research claim at this stage is a tractable four-scalar integrated sensitivity representation with the foreign measure correction handled consistently. It has strong algebraic and prescribed-volatility numerical support. A robust nonlinear calibrator and a claim of competitive performance require the next set of experiments.

## 6. Reproduction commands and reference links

Run from the reviewed repository root:

```bash
OPENBLAS_NUM_THREADS=1 python checks/check_sensitivities.py 100000 100 5.0
OPENBLAS_NUM_THREADS=1 python checks/check_forward_measure.py 100 100000
OPENBLAS_NUM_THREADS=1 python checks/check_sensitivities.py 200000 100 10.0 0.8 2.0
OPENBLAS_NUM_THREADS=1 python checks/compare_estimators.py 1000000 100 20000
```

The matrix checks called `check_matrix_form.run(N=100000, spy=100)` and then `spy=400`.

External sources consulted for specific context (not a systematic literature search):

- [Ögetbil, Ganesan and Hientzsch, Calibrating Local Volatility Models with Stochastic Drift and Diffusion](https://arxiv.org/abs/2009.14764): Monte Carlo calibration algorithms and numerical convergence/accuracy benchmarks.
- [Cozma, Mariapragassam and Reisinger, Calibration of a Hybrid Local-Stochastic Volatility Stochastic Rates Model with a Control Variate Particle Method](https://arxiv.org/abs/1701.06001): a relevant control-variate, market-calibrated comparator.
- [Databento statistics schema](https://databento.com/docs/schemas/statistics): update actions and settlement quality flags.
- [CME 2026 FX product guide](https://www.cmegroup.com/markets/fx/fx-product-guide): EUR/USD options on futures and exercise/delivery structure.

## Appendix: independent audit implementation

This script runs from the audit workspace; adjust the repository import path when using it elsewhere. It uses NumPy and SciPy.

```python
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0,str(Path(__file__).parents[1]/'sir-lsv-malliavin-particle'/'checks'))
import check_forward_measure as fw

kd,sd,kf,sf=.05,.012,.10,.015
rd,rf,rdf=.35,-.30,.40
R=np.array([[1,rd,rf],[rd,1,rdf],[rf,rdf,1.]])

def directional(spy):
    T=5.;N=512;dt=1/spy;eps=1e-4
    rng=np.random.default_rng(99);C=np.linalg.cholesky(R)
    state=np.zeros((5,3,N)); ld=np.zeros(N);lf=np.zeros(N);zf=np.zeros(N);zd=np.zeros(N)
    for j in range(int(T*spy)):
        t=j*dt;dW,dBd,dBf=(rng.normal(size=(N,3))@C.T*np.sqrt(dt)).T
        y,xd,xf=state[:,0,:],state[:,1,:],state[:,2,:]
        sig=.12*(1-.6*np.tanh(2*y));g=-.144*(1-np.tanh(2*y[0])**2)
        mu=(.03+sd**2/(2*kd**2)*(1-np.exp(-kd*t))**2+xd)-(.01+sf**2/(2*kf**2)*(1-np.exp(-kf*t))**2+xf)
        incW=np.tile(dW,(5,1));incd=np.tile(dBd,(5,1));incf=np.tile(dBf,(5,1))
        incW[3:]-=sig[0]*dt;incd[3:]-=rd*sig[0]*dt;incf[3:]-=rf*sig[0]*dt
        for ix,sign in [(1,1),(2,-1)]:
            kick=sign*eps*sd*np.exp(kd*t)*dt
            incW[ix]+=rd*kick;incd[ix]+=kick;incf[ix]+=rdf*kick
        for ix,sign in [(3,1),(4,-1)]:
            kick=sign*eps*sf*np.exp(kf*t)*dt
            incW[ix]+=rf*kick;incd[ix]+=rdf*kick;incf[ix]+=kick
        new=state.copy()
        new[:3,0]=y[:3]+(mu[:3]-.5*sig[:3]**2)*dt+sig[:3]*incW[:3]
        new[:3,1]=xd[:3]-kd*xd[:3]*dt+sd*incd[:3]
        new[:3,2]=xf[:3]+(-kf*xf[:3]-rf*sf*sig[:3])*dt+sf*incf[:3]
        new[3:,0]=y[3:]+(mu[3:]+.5*sig[3:]**2)*dt+sig[3:]*incW[3:]
        new[3:,1]=xd[3:]+(-kd*xd[3:]+rd*sd*sig[3:])*dt+sd*incd[3:]
        new[3:,2]=xf[3:]-kf*xf[3:]*dt+sf*incf[3:]
        vd=sd**2*(-np.expm1(-2*kd*t))/(2*kd);vf=sf**2*(-np.expm1(-2*kf*t))/(2*kf)
        ld,lf,zf,zd=(ld*(1+g*(dW-sig[0]*dt)-kd*dt)+(vd-zf+rd*sd*sig[0])*dt,
            lf*(1+g*dW-kf*dt)+(zd-vf+rf*sf*sig[0])*dt,
            zf+(-(kd+kf)*zf-rf*sf*g*ld+rdf*sd*sf)*dt,
            zd+(-(kd+kf)*zd+rd*sd*g*lf+rdf*sd*sf)*dt)
        state=new
    fd=(state[1,0]-state[2,0])/(2*eps)*np.exp(-kd*T)
    ff=(state[3,0]-state[4,0])/(2*eps)*np.exp(-kf*T)
    err=lambda x,y:np.linalg.norm(x-y)/np.linalg.norm(x)
    print(f'Directional bumps spy={spy}: domestic relative RMS={err(fd,ld):.8g}, foreign relative RMS={err(ff,lf):.8g}')

def constant_vol():
    sig=.12
    def fun(t,y):
        vd=sd**2*(-np.expm1(-2*kd*t))/(2*kd);vf=sf**2*(-np.expm1(-2*kf*t))/(2*kf)
        d,f,z=y
        return [-kd*d+vd-z+rd*sd*sig,-kf*f+z-vf+rf*sf*sig,-(kd+kf)*z+rdf*sd*sf]
    sol=solve_ivp(fun,[0,30],[0,0,0],rtol=1e-11,atol=1e-14,dense_output=True)
    for t in [1.,5.,10.,30.]:
        d,f,z=sol.sol(t);Ad=-np.expm1(-kd*t)/kd;Af=-np.expm1(-kf*t)/kf
        ref=sd**2*Ad**2+sf**2*Af**2+2*rd*sig*sd*Ad-2*rf*sig*sf*Af-2*rdf*sd*sf*Ad*Af
        print(f'Constant-vol t={t}: correction={2*(d-f):.12g}; closed-form={ref:.12g}; abs error={abs(2*(d-f)-ref):.3g}')

def gram_check():
    fw.kd=.0001;fw.kf=.0002;dt=.01;t0=0.;t1=dt
    bad=fw.Sigma(fw.coeffs(t0,t1))
    u,w=np.polynomial.legendre.leggauss(20);s=.5*dt*(u+1);stable=np.zeros((5,5))
    drivers=np.array([0,1,1,2,2]);corr=R[np.ix_(drivers,drivers)]
    for si,wi in zip(s,w):
        k=np.array([1.,np.exp(-fw.kd*(t1-si)),-np.expm1(-fw.kd*(fw.T-si))/fw.kd,np.exp(-fw.kf*(t1-si)),-np.expm1(-fw.kf*(fw.T-si))/fw.kf])
        stable+=wi*.5*np.outer(k,k)*corr
    print('Small-kappa covariance min eigenvalue: original',np.linalg.eigvalsh(bad)[0], 'stable Gram',np.linalg.eigvalsh(stable)[0])
    np.linalg.cholesky(stable)

constant_vol()
for spy in [50,200,800]:directional(spy)
gram_check()

```

## Appendix: full estimator rerun output

```text
T=5.0, 2000000 paths, 100 batches of 20000 particles, steps/yr=100, bandwidth h=0.056
correction in local-variance units x 1e4; 'vol pts' = std of the implied local vol, in vol points
  quantile     K  sigma(K) | reference |   A mean    A std  vol pts |   B mean    B std  vol pts | std A / std B
    0.01   0.412   0.188  |  181.791  |  179.408   15.196   0.404  |  181.648    0.229   0.006  |     66.3
    0.05   0.574   0.178  |  176.634  |  176.169    6.052   0.170  |  176.416    0.173   0.005  |     35.0
    0.10   0.681   0.167  |  168.453  |  168.561    4.229   0.127  |  168.273    0.184   0.006  |     23.0
    0.25   0.879   0.138  |  145.076  |  146.025    2.452   0.089  |  144.971    0.145   0.005  |     16.9
    0.50   1.110   0.105  |  117.602  |  116.725    2.328   0.111  |  117.660    0.055   0.003  |     42.1
    0.75   1.338   0.082  |   99.457  |  100.116    1.835   0.112  |   99.558    0.129   0.008  |     14.2
    0.90   1.545   0.070  |   89.881  |   90.785    2.270   0.163  |   89.953    0.133   0.010  |     17.0
    0.95   1.672   0.064  |   86.087  |   85.971    2.871   0.223  |   86.150    0.151   0.012  |     19.1
    0.99   1.922   0.058  |   81.328  |   80.657    6.270   0.542  |   81.381    0.211   0.018  |     29.8
identity at each strike (both sides x 1e4), z = difference / standard error
   q=0.01  K= 0.412      1.8460     1.8801   z=-1.08
   q=0.05  K= 0.574      9.8057     9.8645   z=-0.60
   q=0.10  K= 0.681     19.6477    19.6524   z=-0.03
   q=0.25  K= 0.879     45.9224    45.6455   z=+0.95
   q=0.50  K= 1.110     73.2513    74.0277   z=-1.51
   q=0.75  K= 1.338     70.7219    70.1528   z=+1.26
   q=0.90  K= 1.545     45.1247    44.6089   z=+1.33
   q=0.95  K= 1.672     28.6329    28.6287   z=+0.01
   q=0.99  K= 1.922      8.4580     8.6431   z=-0.94

```
