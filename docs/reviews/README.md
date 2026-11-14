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
