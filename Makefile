SHELL := /bin/bash

PYTHON ?= .venv/bin/python
SYSTEM_PYTHON ?= python3

.PHONY: check verify verify-logged verify-minimal verify-core verify-geometry \
	verify-theorems verify-regressions verify-derived verify-family \
	verify-external-consequences verify-plane-jc verify-weighted-boundary \
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

verify-plane-jc:
	$(SYSTEM_PYTHON) plane-jc/cas/frontier_125_150.py
	$(PYTHON) plane-jc/cas/boundary_lattice_prefilter.py

verify-weighted-boundary:
	Singular -q scripts/verify_foundational_constant_c_boundary.sing
	Singular -q scripts/verify_foundational_reduced_gluing.sing

verify-core: verify-minimal
	$(PYTHON) scripts/verify_counterexample.py
	$(PYTHON) scripts/audit_map_consistency.py
	$(PYTHON) scripts/verify_normalized_factorization_slice.py
	$(PYTHON) scripts/verify_quadratic_cubic_factorization_invariants.py
	$(PYTHON) scripts/verify_quadratic_cubic_modification_topology.py
	$(PYTHON) scripts/verify_quadratic_cubic_additive_actions.py
	$(PYTHON) scripts/verify_weighted_invariant_jacobian_reduction.py
	$(PYTHON) scripts/verify_weighted_tangent_suspension.py
	$(PYTHON) scripts/verify_foundational_weighted_coefficient_scheme.py
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
	$(PYTHON) scripts/verify_full_boundary_diagram.py
	$(PYTHON) scripts/verify_contact_resultant_endpoint_reduction.py
	$(PYTHON) scripts/verify_counterexample_ladder.py
	$(PYTHON) scripts/verify_parameter_irreducibility.py
	$(PYTHON) scripts/verify_parameter_discriminant.py
	$(PYTHON) scripts/verify_fixed_r_newton_ramification.py
	$(PYTHON) scripts/verify_parameter_galois_groups.py
	$(PYTHON) scripts/verify_parameter_galois_jordan.py
	$(PYTHON) scripts/verify_generalized_cancellation.py
	$(PYTHON) scripts/verify_log_geometry_of_suspensions.py
	$(PYTHON) scripts/verify_three_weight_cancellation.py
	$(PYTHON) scripts/verify_two_factor_resolvent.py
	$(PYTHON) scripts/verify_target_dependent_resolvent.py

verify-external-consequences:
	$(SYSTEM_PYTHON) scripts/verify_long_gaussian_moments.py
	$(SYSTEM_PYTHON) scripts/verify_long_xz_mathieu.py
	$(PYTHON) scripts/verify_long_su2_haar.py
	$(PYTHON) scripts/verify_long_foundational_normalization.py
	$(PYTHON) scripts/verify_rank_two_poisson_preaudit.py
	$(PYTHON) scripts/verify_rank_two_poisson_completion.py
	$(SYSTEM_PYTHON) scripts/audit_rank_two_poisson_completion_independent.py
	$(PYTHON) scripts/verify_long_bcw_79_route.py
	$(SYSTEM_PYTHON) scripts/audit_long_bcw_79_independent.py
	$(PYTHON) scripts/verify_shared_bcw_33_route.py
	$(SYSTEM_PYTHON) scripts/audit_shared_bcw_33_independent.py
	$(PYTHON) scripts/verify_rank_compressed_bcw_24_route.py
	$(SYSTEM_PYTHON) scripts/audit_rank_compressed_bcw_24_independent.py
	$(PYTHON) scripts/verify_constant_kernel_bcw_22_route.py
	$(SYSTEM_PYTHON) scripts/audit_constant_kernel_bcw_22_independent.py
	$(PYTHON) scripts/verify_essential_bcw_profile.py
	$(PYTHON) scripts/verify_essential_bcw_candidate.py
	$(PYTHON) scripts/verify_essential_bcw_21_route.py
	$(SYSTEM_PYTHON) scripts/audit_essential_bcw_21_independent.py
	$(PYTHON) scripts/generate_image_vanishing_counterexamples.py
	$(PYTHON) scripts/generate_identity_slice_counterexamples.py
	$(SYSTEM_PYTHON) scripts/audit_identity_slice_counterexamples_independent.py
	$(PYTHON) scripts/verify_inverse_coordinate_recurrence.py
	$(PYTHON) scripts/audit_bcw_21_linear_quotients.py
	$(PYTHON) scripts/verify_two_parameter_bcw_obstruction.py
	$(SYSTEM_PYTHON) scripts/verify_fixed_gmc_sic_bridge.py
	$(PYTHON) scripts/verify_formal_gaussian_lagrange.py
	$(PYTHON) scripts/verify_weighted_gaussian_bridge.py
	$(PYTHON) scripts/verify_gaussian_moment_fingerprint.py
	$(SYSTEM_PYTHON) scripts/audit_weighted_gaussian_bridge_independent.py

verify-regressions: verify-external-consequences
	$(PYTHON) scripts/verify_degree_five_stable_moduli.py
	$(PYTHON) scripts/verify_degree_five_rank_two_descent.py
	$(SYSTEM_PYTHON) scripts/audit_quartic_independent.py
	$(PYTHON) scripts/verify_generic_discriminant_geometry.py
	$(PYTHON) scripts/verify_canonical_family_image.py
	$(PYTHON) scripts/verify_deformed_seed_boundary.py
	$(PYTHON) scripts/verify_weighted_chebotarev.py
	$(PYTHON) scripts/verify_quartic_weighted_map.py
	$(PYTHON) scripts/verify_quartic_discriminant.py
	$(PYTHON) scripts/verify_quartic_monodromy.py
	$(PYTHON) scripts/verify_external_quartic_islands.py
	$(PYTHON) scripts/verify_decorated_normalization.py
	$(PYTHON) scripts/verify_affine_branch_mark_audit.py
	$(PYTHON) scripts/verify_stable_generator_rigidity.py
	$(PYTHON) scripts/verify_multicluster_ll_comparison.py
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
# CoincidentRootLoci package.  The wrapper uses a pinned Docker image if M2
# is not installed locally.
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

verify: check verify-plane-jc verify-core verify-theorems verify-regressions verify-derived

verify-logged:
	mkdir -p artifacts/verification
	$(PYTHON) scripts/record_environment.py | tee artifacts/verification/environment.txt
	set -o pipefail; $(MAKE) verify 2>&1 | tee artifacts/verification/verify.log
