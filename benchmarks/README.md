# Nonflat synthetic calibration comparison

Completed on 14 November 2026: six PDE calibrations, nine particle calibrations,
and independent Monte Carlo repricing of every frozen calibrated surface.
Read [RESULTS.md](RESULTS.md) before interpreting the measurements.
This is a synthetic spot-option experiment; no CME data enter it.

| File | Purpose |
| --- | --- |
| `forward.cpp` | Native conservative 3D forward operator, one or five fields |
| `synthetic.py` | Model, target generation, calibration, regression, repricing |
| `compare.py` | Experiment stages and report generation |
| `check_solver.py` | Check native operator against the continuous generator |
| `results/comparison.csv` | Recovery, runtime and repricing summary, all 15 runs |
| `results/repricing.csv` | Each option price, paired error and standard error |
| `results/runs.json` | Run parameters and PDE diagnostics |
| `results/diagnostics.json` | Reference-grid sensitivity and diagnostics |
| `results/comparison.png` | Runtime versus recovery error and final profiles |
| `../paper/notes/variance_derivative_cancellation.tex` | Exact cancellation; not an existence proof |

The particle equations are imported from `../sirlsv/`.

## Reproduce

Requires Python 3.10+, NumPy, SciPy, pandas, matplotlib and a C++17 compiler named
`g++`. SymPy is needed for the algebraic check. From the repository root:

```bash
python -m pip install -r requirements.txt matplotlib sympy
python benchmarks/check_solver.py
python checks/check_variance_derivative_cancellation.py
OPENBLAS_NUM_THREADS=1 python -u benchmarks/compare.py reference
OPENBLAS_NUM_THREADS=1 python -u benchmarks/compare.py calibrate
OPENBLAS_NUM_THREADS=1 python -u benchmarks/compare.py validate
OPENBLAS_NUM_THREADS=1 python benchmarks/compare.py report
```

The native operator builds automatically. Reference generation and calibration
reuse existing files in `benchmarks/work/`. Move that directory aside before a
fresh run with changed code or parameters: the cache is filename-based, not
content-addressed. Validation always reruns. A fresh run takes several minutes
and uses substantial memory during the fine solve and multi-model repricing.
Do not run competing CPU-heavy jobs when timing.

The downloadable archive includes reference and reference-check tables, all 15
frozen volatility tables, and metadata in `work/`. With those present, `report`
rebuilds the summary and figure without rerunning calibration; `validate`
independently reprices the stored surfaces. These intermediate binaries are
ignored by git.

The particle implementation supports `bandwidth_scale` and records a direct-kernel
effective-sample-size diagnostic in new runs. The archived study predates that
instrumentation and does not report ESS. This diagnostic approximates the binned
regression weights and does not adjust for antithetic dependence.

The requested nonflat local-existence and stability theorem remains unfinished.
Neither the numerical study nor the algebraic cancellation proves it.
