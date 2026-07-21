SHELL := /bin/bash

PYTHON ?= .venv/bin/python
SYSTEM_PYTHON ?= python3

.PHONY: check verify verify-logged verify-minimal verify-core verify-geometry \
	verify-theorems verify-regressions verify-derived verify-family \
	verify-master \
	verify-quartic verify-normal-forms verify-formal verify-lean-c01 \
	verify-foundations verify-foundations-formal \
	verify-coincident-root-loci scan-weighted-seeds

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
	$(PYTHON) scripts/verify_marked_root_model.py
	$(PYTHON) scripts/image_nonproperness.py
	$(PYTHON) scripts/verify_exceptional_fibers.py
	$(PYTHON) scripts/verify_image_nonproperness_inclusions.py

verify-geometry: verify-core

verify-theorems:
	$(MAKE) verify-master
	$(PYTHON) scripts/verify_weighted_seed_schema.py
	$(PYTHON) scripts/verify_weighted_seed_theorem.py
	$(PYTHON) scripts/verify_weighted_marked_root_model.py
	$(PYTHON) scripts/verify_universal_discriminant_incidences.py
	$(PYTHON) scripts/verify_contact_partition_strata.py
	$(PYTHON) scripts/verify_uniform_exceptional_seed_theorem.py
	$(PYTHON) scripts/verify_maximal_phi_irreducibility.py
	$(PYTHON) scripts/verify_contact_atom_principle.py
	$(PYTHON) scripts/verify_unique_omitted_value.py
	$(PYTHON) scripts/verify_component_normalization.py
	$(PYTHON) scripts/verify_degree12_branch_intersection.py
	$(PYTHON) scripts/verify_dicritical_divisors.py
	$(PYTHON) scripts/verify_c16_blowup_geometry.py
	$(PYTHON) scripts/classify_transfer_block_k2.py
	$(PYTHON) scripts/classify_transfer_block_k3.py
	$(PYTHON) scripts/classify_transfer_block_k4.py
	$(PYTHON) scripts/verify_all_k_transfer_block.py
	$(SYSTEM_PYTHON) scripts/verify_c22_deformation_audit.py
	$(PYTHON) scripts/verify_global_affine_rigidity.py
	$(PYTHON) scripts/verify_universal_factorization_geometry.py
	$(PYTHON) scripts/verify_allocation_hensel_product.py
	$(PYTHON) scripts/verify_mixed_allocation_equalizer.py
	$(PYTHON) scripts/verify_omitted_value_classification.py
	$(PYTHON) scripts/verify_repeated_root_boundary.py
	$(PYTHON) scripts/verify_effective_chebotarev.py

verify-master:
	$(PYTHON) scripts/verify_master_universal.py
	$(PYTHON) scripts/verify_master_instances.py
	$(PYTHON) scripts/verify_resolvent_ramification_signature.py

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
	$(PYTHON) scripts/audit_c15_independent.py
	$(PYTHON) scripts/generate_c15_consequences.py

verify-derived: verify-normal-forms

# Optional formal replication. This fetches Dean Cureton's separately authored
# Lean project at the audited commit recorded in verified/LEAN_C01.md; it is kept
# out of the default target because it downloads a pinned Lean/mathlib toolchain.
verify-lean-c01:
	bash scripts/verify_lean_c01.sh

verify-formal: verify-lean-c01

verify-foundations: verify-core
	$(PYTHON) scripts/verify_weighted_seed_schema.py
	$(PYTHON) scripts/verify_weighted_seed_theorem.py
	$(PYTHON) scripts/verify_weighted_marked_root_model.py

verify-foundations-formal: verify-foundations verify-lean-c01

# Optional independent bounded-degree comparison with Macaulay2's classical
# CoincidentRootLoci package.  The wrapper uses a pinned Docker image if M2 is
# not installed locally.
verify-coincident-root-loci:
	bash scripts/verify_coincident_root_slices.sh

verify: check verify-core verify-theorems verify-regressions verify-derived

verify-logged:
	mkdir -p artifacts/verification
	$(PYTHON) scripts/record_environment.py | tee artifacts/verification/environment.txt
	set -o pipefail; $(MAKE) verify 2>&1 | tee artifacts/verification/verify.log

scan-weighted-seeds:
	$(PYTHON) scripts/scan_weighted_seeds.py
