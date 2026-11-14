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

## Moment-closure note

`2026-11-14_moment_closure_summary.md` is the second agent's write-up of `paper/discounted_moment_closure_note.pdf`; `2026-11-14_moment_closure_validation.json` is the output of `checks/check_moment_closure.py` with its default settings.

Review of the note (14 November 2026):

- The auxiliary matrices of the note are those of `sirlsv/sensitivities.py`.
- Theorem 1 (closure), Proposition 3 (compatibility identities and the marginal equation) and both reductions of Theorem 4 were re-derived by hand. No error found.
- The script reproduces the z-scores quoted in the write-up. Its Monte Carlo part uses affine test functions only, so the second-order terms of the generator are covered by the symbolic part, not by the simulation.
- Proposition 2 (uniqueness at fixed coefficient) is a sketch of a standard duality argument. It was read, not checked in detail.
- Theorem 4 is conditional: it assumes a regular solution of the nonlinear system exists.
