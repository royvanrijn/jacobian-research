PYTHON ?= .venv/bin/python
SYSTEM_PYTHON ?= python3

.PHONY: check verify verify-minimal verify-core verify-geometry verify-quartic verify-normal-forms scan-weighted-seeds

check:
	$(PYTHON) -m compileall -q jcsearch scripts
	$(PYTHON) scripts/check_readme_links.py

verify-minimal:
	$(SYSTEM_PYTHON) scripts/verify_counterexample_independent.py

verify-core: verify-minimal
	$(PYTHON) scripts/verify_counterexample.py
	$(PYTHON) scripts/audit_map_consistency.py
	$(PYTHON) scripts/cubic_model.py

verify-geometry:
	$(PYTHON) scripts/image_nonproperness.py
	$(PYTHON) scripts/verify_exceptional_fibers.py
	$(PYTHON) scripts/verify_image_nonproperness_inclusions.py

verify-quartic:
	$(PYTHON) scripts/verify_weighted_seed_schema.py
	$(PYTHON) scripts/verify_weighted_seed_theorem.py
	$(PYTHON) scripts/verify_canonical_family_image.py
	$(PYTHON) scripts/verify_deformed_seed_boundary.py
	$(PYTHON) scripts/verify_omitted_value_classification.py
	$(PYTHON) scripts/verify_repeated_root_boundary.py
	$(PYTHON) scripts/verify_weighted_chebotarev.py
	$(PYTHON) scripts/verify_quartic_weighted_map.py
	$(PYTHON) scripts/verify_quartic_discriminant.py
	$(PYTHON) scripts/verify_quartic_monodromy.py
	$(PYTHON) scripts/verify_quartic_c0_fibers.py
	$(PYTHON) scripts/verify_quartic_nonproperness_paths.py
	$(PYTHON) scripts/verify_quartic_properness_converse.py
	$(PYTHON) scripts/verify_quartic_singular_locus.py
	$(PYTHON) scripts/verify_quartic_image.py

verify-normal-forms:
	$(PYTHON) scripts/cubic_homogeneous_reduction.py
	$(PYTHON) scripts/verify_cubic_homogeneous_counterexample.py
	$(PYTHON) scripts/cubic_linear_reduction.py
	$(PYTHON) scripts/verify_cubic_linear_counterexample.py

verify: check verify-core verify-geometry verify-quartic

scan-weighted-seeds:
	$(PYTHON) scripts/scan_weighted_seeds.py
