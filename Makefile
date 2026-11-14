PYTHON ?= python3
SETTLEMENTS ?= eurusd_cme_options_2026-11-12.csv

.PHONY: paper checks checks-quick calibration clean

paper:
	cd paper && pdflatex -interaction=nonstopmode SIR_LSV_particle_FX.tex >/dev/null && pdflatex -interaction=nonstopmode SIR_LSV_particle_FX.tex >/dev/null && pdflatex -interaction=nonstopmode SIR_LSV_particle_FX.tex | grep -E "Output written|Warning" || true

# deterministic checks, then the Monte Carlo checks at a size that runs in about a minute.
# At this size the time-step bias is visible (differences of up to about 2 standard errors);
# use `make checks` for the sizes quoted in the paper.
checks-quick:
	$(PYTHON) checks/check_covariance_stability.py
	$(PYTHON) checks/check_gaussian_existence.py
	$(PYTHON) checks/check_sensitivities.py 100000 50 5.0
	$(PYTHON) checks/check_forward_measure.py 50 50000
	$(PYTHON) checks/check_matrix_form.py 50
	$(PYTHON) checks/compare_estimators.py 100000 50 20000

# the sizes quoted in the paper and in the README (several minutes)
checks:
	$(PYTHON) checks/check_covariance_stability.py
	$(PYTHON) checks/check_gaussian_existence.py
	$(PYTHON) checks/check_sensitivities.py 1000000 200 5.0
	$(PYTHON) checks/check_forward_measure.py 200 500000
	$(PYTHON) checks/check_matrix_form.py 100
	$(PYTHON) checks/compare_estimators.py 1000000 100 20000

# needs the licensed settlements file: make calibration SETTLEMENTS=path/to/file.csv
calibration:
	$(PYTHON) calibration/build_surface.py $(SETTLEMENTS)
	$(PYTHON) calibration/calibrate.py det --rates none --estimator none
	$(PYTHON) calibration/calibrate.py base --rates base
	$(PYTHON) calibration/calibrate.py base_nocorr --rates base --estimator none
	$(PYTHON) calibration/calibrate.py stress --rates stress
	$(PYTHON) calibration/calibrate.py stress_nocorr --rates stress --estimator none
	$(PYTHON) calibration/summarise.py

clean:
	cd paper && rm -f *.aux *.log *.out *.toc *.synctex.gz
	find . -name __pycache__ -type d -prune -exec rm -rf {} +
