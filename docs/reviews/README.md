# Reviews

Two notes written on 14 November 2026 by a second agent (Codex), which audited the paper and the code as they stood at commit `717f936`. File names, line numbers and commands in them refer to that commit, before the repository was reorganised.

| Finding | Status |
|---|---|
| Sign of the rate term in the one-factor forward SDE | Fixed |
| Time indexing of the one-factor algorithm | Fixed |
| Covariance formulas of the forward-measure scheme lose positive-definiteness for small mean reversion or tiny steps | Fixed: `sirlsv/gaussian_covariance.py`, regression test in `checks/check_covariance_stability.py` |
| General-rate Clark-Ocone statement too broad | Qualified in the paper |
| Finite-particle sensitivities are not exact derivatives of the interacting system | Stated in the paper |
| "The remedy is in the inputs" when the variance floor binds | Rewritten in the paper |
| Multi-factor outlook promised two processes per factor | Rewritten in the paper |
| Data script ignores `update_action` | Open |
| The data are options on futures, the paper's calibration input is options on spot | Open. The calibration test treats them as options on the forward to their expiry; see `calibration/RESULTS.md` |

The existence section of the paper (`paper/sections/3e_fx_existence.tex`) comes from the second note.

## Moment-closure section

`2026-11-14_discounted_moment_closure_results.md` is the second agent's write-up of the section "Exact closure of discounted first moments" (`paper/sections/3f_fx_moment_closure.tex`; standalone version in `paper/notes/`). Paths in it refer to that agent's working tree. `2026-11-14_moment_closure_validation.json` is the output of `checks/check_moment_closure.py` with its default settings.

Review (14 November 2026):

- The auxiliary matrices of the section are those of `sirlsv/sensitivities.py`.
- Theorem 1 (closure), Proposition 3 (compatibility identities and the marginal equation) and both reductions of Theorem 4 were re-derived by hand. No error found.
- The script reproduces the z-scores quoted in the write-up. Its Monte Carlo part uses affine test functions only, so the second-order terms of the generator are covered by the symbolic part, not by the simulation.
- Proposition 2 (uniqueness at fixed coefficient) is a sketch of a standard duality argument. It was read, not checked in detail.
- Theorem 4 is conditional: it assumes a regular solution of the nonlinear system exists.

## Prior art and novelty

- `2026-11-14_bloch_nakashima_prior_art_comparison.md`: equation-level comparison with Bloch and Nakashima (2008), made by the second agent from the full paper. Conclusion: close and material prior art for the idea and for the form of the correction; no equivalent of the exact four-process system was found there. The three mathematical objections in the note (conditional covariance, log-variance matching, rate-to-drift sensitivity) were checked and are correct as stated. What the note says the paper contains was not checked: the paper is not in this repository.
- `2026-11-14_publication_novelty_and_research_assessment.md`: the second agent's assessment of novelty and of possible research directions. It predates the comparison above in part; where they differ, the comparison is the later document.

## Synthetic PDE comparison

Report in `benchmarks/RESULTS.md`. Review (14 November 2026): the operator check and the cancellation check pass; the summary table rebuilds exactly from the saved runs; three of the fifteen calibrations and their repricing were re-run from scratch and agree with the saved results to rounding; the cancellation proposition was re-derived by hand. The archive's intermediate tables are not committed.
