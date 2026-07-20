PYTHON ?= .venv/bin/python

.PHONY: check verify verify-core verify-search verify-certificates

check:
	$(PYTHON) -m compileall -q jcsearch scripts
	$(PYTHON) scripts/check_readme_links.py

verify-core:
	$(PYTHON) scripts/verify_counterexample.py
	$(PYTHON) scripts/cubic_model.py
	$(PYTHON) scripts/image_nonproperness.py

verify-search:
	$(PYTHON) scripts/validate_ladder.py
	$(PYTHON) scripts/verify_translated_box_theorem.py
	$(PYTHON) scripts/verify_newton_translation.py
	$(PYTHON) scripts/newton_9_27_regression.py

verify-certificates:
	$(PYTHON) scripts/verify_certificates.py

verify: check verify-core verify-search verify-certificates
