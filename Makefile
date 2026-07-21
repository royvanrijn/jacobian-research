SHELL := /bin/bash

PYTHON ?= .venv/bin/python
SYSTEM_PYTHON ?= python3

.PHONY: check verify verify-logged verify-minimal verify-core verify-geometry \
	verify-theorems verify-regressions verify-derived verify-family \
	verify-master \
	verify-quartic verify-normal-forms verify-formal verify-lean-foundational \
	verify-foundations verify-foundations-formal \
	verify-coincident-root-loci verify-papers clean-papers

check:
	$(PYTHON) -m compileall -q jcsearch scripts
	$(PYTHON) scripts/check_markdown_links.py
	$(PYTHON) scripts/audit_status.py

verify-minimal:
	$(SYSTEM_PYTHON) scripts/verify_counterexample_independent.py

verify-core: verify-minimal
	$(PYTHON) scripts/verify_counterexample.py
	$(PYTHON) scripts/audit_map_consistency.py
	$(PYTHON) scripts/cubic_model.py
	$(PYTHON) scripts/audit_foundational_invariance_regression.py
	$(PYTHON) scripts/verify_symplectic_weyl_lift.py
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
	$(SYSTEM_PYTHON) scripts/audit_weighted_independent.py
	$(PYTHON) scripts/verify_universal_discriminant_incidences.py
	$(PYTHON) scripts/verify_contact_partition_strata.py
	$(PYTHON) scripts/verify_uniform_exceptional_seed_theorem.py
	$(PYTHON) scripts/verify_maximal_phi_irreducibility.py
	$(PYTHON) scripts/verify_contact_atom_principle.py
	$(PYTHON) scripts/verify_unique_omitted_value.py
	$(PYTHON) scripts/verify_component_normalization.py
	$(PYTHON) scripts/verify_degree12_branch_intersection.py
	$(SYSTEM_PYTHON) scripts/audit_degree_twelve_independent.py
	$(PYTHON) scripts/verify_dicritical_divisors.py
	$(PYTHON) scripts/verify_dicritical_blowup_geometry.py
	$(PYTHON) scripts/verify_omitted_value_classification.py
	$(PYTHON) scripts/verify_repeated_root_boundary.py
	$(PYTHON) scripts/verify_effective_chebotarev.py

verify-master:
	$(SYSTEM_PYTHON) scripts/audit_boundary_exhaustion_independent.py
	$(SYSTEM_PYTHON) scripts/audit_thick_intersection_local.py
	$(PYTHON) scripts/verify_master_universal.py
	$(PYTHON) scripts/verify_master_instances.py
	$(PYTHON) scripts/verify_resolvent_ramification_signature.py
	$(PYTHON) scripts/verify_target_fixed_parameter_rigidity.py
	$(PYTHON) scripts/verify_boundary_intersection_obstruction.py
	$(PYTHON) scripts/verify_scheme_boundary_all_parameters.py
	$(PYTHON) scripts/verify_parameter_irreducibility.py
	$(PYTHON) scripts/verify_parameter_discriminant.py
	$(PYTHON) scripts/verify_parameter_galois_groups.py
	$(PYTHON) scripts/verify_parameter_galois_jordan.py
	$(PYTHON) scripts/verify_generalized_cancellation.py
	$(PYTHON) scripts/verify_three_weight_cancellation.py
	$(PYTHON) scripts/verify_two_factor_resolvent.py
	$(PYTHON) scripts/verify_target_dependent_resolvent.py

verify-regressions:
	$(SYSTEM_PYTHON) scripts/audit_quartic_independent.py
	$(PYTHON) scripts/verify_generic_discriminant_geometry.py
	$(PYTHON) scripts/verify_canonical_family_image.py
	$(PYTHON) scripts/verify_deformed_seed_boundary.py
	$(PYTHON) scripts/verify_weighted_chebotarev.py
	$(PYTHON) scripts/verify_quartic_weighted_map.py
	$(PYTHON) scripts/verify_quartic_discriminant.py
	$(PYTHON) scripts/verify_quartic_monodromy.py
	$(PYTHON) scripts/verify_external_quartic_islands.py
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
	$(PYTHON) scripts/audit_stable_normal_form_independent.py
	$(PYTHON) scripts/generate_stable_normal_form_consequences.py

verify-derived: verify-normal-forms

# Optional formal replication. This fetches Dean Cureton's separately authored
# Lean project at the audited commit recorded in verified/LEAN_FOUNDATIONAL_MAP.md; it is kept
# out of the default target because it downloads a pinned Lean/mathlib toolchain.
verify-lean-foundational:
	bash scripts/verify_lean_foundational_map.sh

verify-formal: verify-lean-foundational

verify-foundations: verify-core
	$(PYTHON) scripts/verify_weighted_seed_schema.py
	$(PYTHON) scripts/verify_weighted_seed_theorem.py
	$(PYTHON) scripts/verify_weighted_marked_root_model.py
	$(SYSTEM_PYTHON) scripts/audit_weighted_independent.py

verify-foundations-formal: verify-foundations verify-lean-foundational

# Optional independent bounded-degree comparison with Macaulay2's classical
# CoincidentRootLoci package.  The wrapper uses a pinned Docker image if M2 is
# not installed locally.
verify-coincident-root-loci:
	bash scripts/verify_coincident_root_slices.sh

verify-papers:
	latexmk -cd -pdf -interaction=nonstopmode -halt-on-error papers/core-counterexample/main.tex
	latexmk -cd -pdf -interaction=nonstopmode -halt-on-error papers/discriminant-pencils/main.tex
	latexmk -cd -pdf -interaction=nonstopmode -halt-on-error papers/marked-root-multiplicity/main.tex

clean-papers:
	latexmk -cd -C papers/core-counterexample/main.tex
	latexmk -cd -C papers/discriminant-pencils/main.tex
	latexmk -cd -C papers/marked-root-multiplicity/main.tex
	$(RM) papers/core-counterexample/main.bbl

verify: check verify-core verify-theorems verify-regressions verify-derived

verify-logged:
	mkdir -p artifacts/verification
	$(PYTHON) scripts/record_environment.py | tee artifacts/verification/environment.txt
	set -o pipefail; $(MAKE) verify 2>&1 | tee artifacts/verification/verify.log
