SHELL := /bin/bash

PYTHON ?= .venv/bin/python
SYSTEM_PYTHON ?= python3

.PHONY: check verify verify-logged verify-minimal verify-core verify-geometry \
	verify-theorems verify-regressions verify-derived verify-family \
	verify-quartic verify-normal-forms scan-weighted-seeds

check:
	$(PYTHON) -m compileall -q jcsearch scripts
	$(PYTHON) scripts/check_markdown_links.py
	$(PYTHON) scripts/audit_claim_ledger.py

verify-minimal:
	$(SYSTEM_PYTHON) scripts/verify_counterexample_independent.py

verify-core: verify-minimal
	$(PYTHON) scripts/verify_counterexample.py
	$(PYTHON) scripts/audit_map_consistency.py
	$(PYTHON) scripts/cubic_model.py
	$(PYTHON) scripts/image_nonproperness.py
	$(PYTHON) scripts/verify_exceptional_fibers.py
	$(PYTHON) scripts/verify_image_nonproperness_inclusions.py

verify-geometry: verify-core

verify-theorems:
	$(PYTHON) scripts/verify_weighted_seed_schema.py
	$(PYTHON) scripts/verify_weighted_seed_theorem.py
	$(PYTHON) scripts/verify_universal_discriminant_incidences.py
	$(PYTHON) scripts/verify_contact_partition_strata.py
	$(PYTHON) scripts/verify_uniform_exceptional_seed_theorem.py
	$(PYTHON) scripts/verify_maximal_phi_irreducibility.py
	$(PYTHON) scripts/verify_contact_atom_principle.py
	$(PYTHON) scripts/verify_unique_omitted_value.py
	$(PYTHON) scripts/verify_component_normalization.py
	$(PYTHON) scripts/verify_degree12_branch_intersection.py
	$(PYTHON) scripts/verify_dicritical_divisors.py
	$(PYTHON) scripts/classify_transfer_block_k2.py
	$(PYTHON) scripts/verify_omitted_value_classification.py
	$(PYTHON) scripts/verify_repeated_root_boundary.py

verify-regressions:
	$(PYTHON) scripts/verify_generic_discriminant_geometry.py
	$(PYTHON) scripts/verify_canonical_family_image.py
	$(PYTHON) scripts/verify_deformed_seed_boundary.py
	$(PYTHON) scripts/verify_weighted_chebotarev.py
	$(PYTHON) scripts/verify_quartic_weighted_map.py
	$(PYTHON) scripts/verify_quartic_discriminant.py
	$(PYTHON) scripts/verify_quartic_monodromy.py
	$(PYTHON) scripts/verify_quartic_c0_fibers.py
	$(PYTHON) scripts/verify_quartic_nonproperness_paths.py
	$(PYTHON) scripts/verify_quartic_properness_converse.py
	$(PYTHON) scripts/verify_quartic_singular_locus.py
	$(PYTHON) scripts/verify_quartic_image.py

# Backward-compatible names retained for existing commands and links.
verify-family: verify-theorems verify-regressions
verify-quartic: verify-regressions

verify-normal-forms:
	$(PYTHON) scripts/cubic_homogeneous_reduction.py
	$(PYTHON) scripts/verify_cubic_homogeneous_counterexample.py
	$(PYTHON) scripts/cubic_linear_reduction.py
	$(PYTHON) scripts/verify_cubic_linear_counterexample.py

verify-derived: verify-normal-forms

verify: check verify-core verify-theorems verify-regressions verify-derived

verify-logged:
	mkdir -p artifacts/verification
	$(PYTHON) scripts/record_environment.py | tee artifacts/verification/environment.txt
	set -o pipefail; $(MAKE) verify 2>&1 | tee artifacts/verification/verify.log

scan-weighted-seeds:
	$(PYTHON) scripts/scan_weighted_seeds.py
