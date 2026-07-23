SHELL := /bin/bash

PYTHON ?= .venv/bin/python
SYSTEM_PYTHON ?= python3

.PHONY: check verify verify-logged verify-minimal verify-core verify-geometry \
	verify-theorems verify-regressions verify-derived verify-family \
	verify-external-consequences verify-restricted-minima verify-two-real-gmc verify-counterexample-scoreboard verify-plane-jc verify-plane-case2-residue-strata verify-plane-case2-j1-endpoint verify-plane-case2-maximal-gcd verify-plane-case2-gcd6 verify-plane-poisson-radical verify-plane-poisson-primary-charts verify-plane-poisson-separators verify-plane-poisson-primary-filtration verify-plane-poisson-filtered-modules verify-weighted-boundary \
	verify-master \
	verify-quartic verify-normal-forms verify-formal verify-lean-foundational \
	verify-foundations verify-foundations-formal \
	verify-coincident-root-loci verify-papers verify-ritt-boundary \
	verify-ritt-2-complex verify-hessian-synchronization \
	verify-contact-r6 verify-contact-branch-schema verify-contact-r7-asymptotic \
	verify-parameter-dusart-frontier verify-parameter-sharp-dusart-frontier \
	verify-parameter-adaptive-dusart-frontier \
	verify-minimal-boundary \
	render-status clean-papers

check:
	$(PYTHON) -m compileall -q jcsearch scripts
	$(PYTHON) scripts/check_markdown_links.py
	$(PYTHON) scripts/audit_status.py

render-status:
	$(SYSTEM_PYTHON) scripts/render_status.py

verify-minimal:
	$(SYSTEM_PYTHON) scripts/verify_counterexample_independent.py

verify-plane-jc:
	$(SYSTEM_PYTHON) plane-jc/cas/frontier_125_150.py
	$(PYTHON) plane-jc/cas/boundary_lattice_prefilter.py
	$(PYTHON) plane-jc/cas/test_intrinsic_a2_boundary.py
	$(PYTHON) plane-jc/cas/test_log_boundary_compiler.py
	$(PYTHON) plane-jc/cas/test_poisson_square_rigidity.py
	$(PYTHON) plane-jc/cas/test_poisson_square_filtered_modules.py

verify-plane-case2-residue-strata:
	$(PYTHON) plane-jc/cas/audit_case2_residue_strata.py

verify-plane-case2-j1-endpoint:
	$(PYTHON) plane-jc/cas/case2_infinity_resolution.py

verify-plane-case2-maximal-gcd:
	$(PYTHON) plane-jc/cas/audit_case2_maximal_gcd.py

verify-plane-case2-gcd6:
	$(PYTHON) plane-jc/cas/audit_case2_gcd6.py

verify-plane-poisson-radical:
	Singular -q plane-jc/cas/poisson_square_radical.sing

verify-plane-poisson-primary-charts:
	Singular -q plane-jc/cas/poisson_square_primary_charts.sing

verify-plane-poisson-separators:
	Singular -q plane-jc/cas/poisson_square_separator_primary.sing

verify-plane-poisson-primary-filtration:
	Singular -q plane-jc/cas/poisson_square_normalized_defect.sing

verify-plane-poisson-filtered-modules:
	$(PYTHON) plane-jc/cas/test_poisson_square_filtered_modules.py

verify-weighted-boundary:
	Singular -q scripts/verify_foundational_constant_c_boundary.sing
	Singular -q scripts/verify_foundational_reduced_gluing.sing

verify-ritt-boundary:
	bash scripts/verify_degree_six_ritt_boundary_atlas.sh

verify-ritt-2-complex:
	$(PYTHON) scripts/verify_degree30_ritt_2_complex.py

verify-hessian-synchronization:
	$(PYTHON) scripts/verify_hessian_synchronization_lifts.py

audit-degree30-hessian-synchronization-pairs:
	$(PYTHON) scripts/audit_degree30_hessian_synchronization_pairs.py

verify-minimal-boundary:
	$(PYTHON) scripts/verify_minimal_boundary_cubic.py
	$(PYTHON) scripts/verify_cubic_marking_frontier.py
	$(PYTHON) scripts/verify_cubic_normalization_frontend.py
	Singular -q scripts/verify_cubic_double_saturation.sing
	$(PYTHON) scripts/verify_cubic_gauge_straightening.py
	$(PYTHON) scripts/verify_cubic_gauge_first_obstruction.py

verify-contact-r6:
	$(PYTHON) scripts/verify_contact_resultant_r6_effective.py

verify-parameter-dusart-frontier:
	$(PYTHON) scripts/verify_parameter_irreducibility_dusart_frontier.py

verify-parameter-sharp-dusart-frontier:
	$(PYTHON) scripts/verify_parameter_irreducibility_sharp_dusart_frontier.py

verify-parameter-adaptive-dusart-frontier:
	$(PYTHON) scripts/verify_parameter_irreducibility_adaptive_dusart_frontier.py

verify-contact-branch-schema:
	$(PYTHON) scripts/verify_contact_resultant_fixed_r_branch_schema.py

verify-contact-r7-asymptotic:
	$(PYTHON) scripts/verify_contact_resultant_r7_asymptotic.py

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
	$(MAKE) verify-minimal-boundary
	$(PYTHON) scripts/verify_weighted_seed_schema.py
	$(PYTHON) scripts/verify_weighted_seed_theorem.py
	$(PYTHON) scripts/verify_all_degree_rational_fibers.py
	$(PYTHON) scripts/verify_real_fiber_spectrum.py
	$(PYTHON) scripts/verify_hasse_keller_fiber.py
	$(PYTHON) scripts/verify_weighted_marked_root_model.py
	$(SYSTEM_PYTHON) scripts/audit_weighted_independent.py
	$(PYTHON) scripts/verify_universal_discriminant_incidences.py
	$(PYTHON) scripts/verify_contact_partition_strata.py
	$(PYTHON) scripts/verify_uniform_exceptional_seed_theorem.py
	$(PYTHON) scripts/verify_maximal_phi_irreducibility.py
	$(PYTHON) scripts/verify_contact_atom_principle.py
	$(PYTHON) scripts/verify_unique_omitted_value.py
	$(PYTHON) scripts/verify_component_normalization.py
	$(PYTHON) scripts/verify_nonsurjective_enumerative_geometry.py
	$(PYTHON) scripts/verify_degree12_branch_intersection.py
	$(SYSTEM_PYTHON) scripts/audit_degree_twelve_independent.py
	$(PYTHON) scripts/verify_exceptional_partition_lattice.py
	$(PYTHON) scripts/verify_degree18_triple_intersection.py
	$(PYTHON) scripts/verify_omitted_intersection_algebra.py
	$(PYTHON) scripts/verify_ferrand_norm_transfer_blocks.py
	$(PYTHON) scripts/verify_maximally_collided_transfer_cones.py --compare-current-power-ideal
	$(PYTHON) scripts/verify_maximally_collided_transfer_cones.py --method raw --min-degree 12 --max-degree 12
	$(PYTHON) scripts/verify_maximally_collided_transfer_cones.py --method modular --basis-engine slimgb --prime 3 --min-degree 11 --max-degree 16
	$(PYTHON) scripts/verify_maximally_collided_transfer_cones.py --method modular --basis-engine slimgb --prime 3 --min-degree 2 --max-degree 10 --compare-characteristic-three-model
	$(PYTHON) scripts/verify_maximally_collided_transfer_cones.py --method divided-three --basis-engine slimgb --min-degree 17 --max-degree 17 --timeout 240
	$(PYTHON) scripts/verify_dicritical_divisors.py
	$(PYTHON) scripts/verify_dicritical_blowup_geometry.py
	$(PYTHON) scripts/verify_omitted_value_classification.py
	$(PYTHON) scripts/verify_repeated_root_boundary.py
	$(PYTHON) scripts/verify_effective_chebotarev.py
	$(PYTHON) scripts/verify_global_sunada_keller.py
	$(PYTHON) scripts/verify_davenport_cox_boundary.py
	$(PYTHON) scripts/verify_davenport_tangent_mark_curve.py
	$(PYTHON) scripts/verify_davenport_proportional_tangent_sections.py
	$(PYTHON) scripts/verify_davenport_weighted_glue_obstruction.py
	$(PYTHON) scripts/verify_davenport_derivative_center_mismatch.py
	$(PYTHON) scripts/verify_davenport_boundary_involution.py
	$(PYTHON) scripts/verify_davenport_node_separation.py
	$(PYTHON) scripts/verify_davenport_post_coordinate_attacks.py
	$(PYTHON) scripts/verify_stratified_adelic_engineering.py

verify-master:
	$(SYSTEM_PYTHON) scripts/audit_boundary_exhaustion_independent.py
	$(SYSTEM_PYTHON) scripts/audit_thick_intersection_local.py
	$(SYSTEM_PYTHON) scripts/verify_degreewise_multiplicity_count.py
	$(PYTHON) scripts/verify_master_universal.py
	$(PYTHON) scripts/verify_root_engineered_quadratic_gauge.py
	$(PYTHON) scripts/verify_master_instances.py
	$(PYTHON) scripts/verify_resolvent_ramification_signature.py
	$(PYTHON) scripts/verify_target_fixed_parameter_rigidity.py
	$(PYTHON) scripts/verify_cancellation_parameter_faithfulness.py
	$(PYTHON) scripts/verify_boundary_intersection_obstruction.py
	$(PYTHON) scripts/verify_scheme_boundary_all_parameters.py
	$(PYTHON) scripts/verify_full_boundary_diagram.py
	$(PYTHON) scripts/verify_contact_resultant_endpoint_reduction.py
	$(PYTHON) scripts/verify_contact_resultant_r4.py
	$(PYTHON) scripts/verify_contact_resultant_irreducible_ranges.py
	$(PYTHON) scripts/verify_counterexample_ladder.py
	$(PYTHON) scripts/verify_parameter_irreducibility.py
	$(PYTHON) scripts/verify_parameter_discriminant.py
	$(PYTHON) scripts/verify_fixed_r_newton_ramification.py
	$(PYTHON) scripts/verify_parameter_galois_groups.py
	$(PYTHON) scripts/verify_parameter_galois_jordan.py
	$(PYTHON) scripts/verify_generalized_cancellation.py
	$(PYTHON) scripts/verify_log_geometry_of_suspensions.py
	$(PYTHON) scripts/verify_reciprocal_link_classifier.py
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
	$(MAKE) verify-restricted-minima
	$(PYTHON) scripts/generate_image_vanishing_counterexamples.py
	$(PYTHON) scripts/generate_identity_slice_counterexamples.py
	$(SYSTEM_PYTHON) scripts/audit_identity_slice_counterexamples_independent.py
	$(PYTHON) scripts/verify_inverse_coordinate_recurrence.py
	$(PYTHON) scripts/audit_bcw_21_linear_quotients.py
	$(SYSTEM_PYTHON) scripts/audit_bcw_21_affine_vector_symmetries.py
	$(PYTHON) scripts/verify_two_parameter_bcw_obstruction.py
	$(SYSTEM_PYTHON) scripts/verify_fixed_gmc_sic_bridge.py
	$(PYTHON) scripts/verify_formal_gaussian_lagrange.py
	$(PYTHON) scripts/verify_weighted_gaussian_bridge.py
	$(PYTHON) scripts/verify_gaussian_moment_fingerprint.py
	$(SYSTEM_PYTHON) scripts/audit_weighted_gaussian_bridge_independent.py

verify-restricted-minima:
	$(PYTHON) scripts/verify_index_reduced_bcw_22_route.py
	$(SYSTEM_PYTHON) scripts/audit_index_reduced_bcw_22_independent.py
	$(PYTHON) scripts/verify_rank_reduced_bcw_24_route.py
	$(SYSTEM_PYTHON) scripts/audit_rank_reduced_bcw_24_independent.py
	$(PYTHON) scripts/verify_hessian_rank_reduced_bcw_22_route.py
	$(SYSTEM_PYTHON) scripts/audit_hessian_rank_reduced_bcw_22_independent.py
	$(PYTHON) scripts/analyze_cotangent_kernel_excess.py
	$(PYTHON) scripts/verify_index_three_inverse_model.py
	$(PYTHON) scripts/derive_index_three_tree_obstruction.py
	$(PYTHON) scripts/verify_index_three_rank_normal_form.py
	$(PYTHON) scripts/verify_restricted_minima_frontier.py

verify-two-real-gmc:
	$(PYTHON) scripts/verify_two_real_gmc_frontier.py
	$(PYTHON) scripts/verify_two_real_gmc_symmetric_chart.py
	$(PYTHON) scripts/verify_two_real_gmc_remaining_four_weight.py

verify-counterexample-scoreboard: verify-two-real-gmc
	$(PYTHON) scripts/verify_minimal_counterexample_scoreboard.py

verify-regressions: verify-external-consequences
	$(PYTHON) scripts/verify_degree_five_stable_moduli.py
	$(PYTHON) scripts/verify_degree_five_rank_two_descent.py
	$(PYTHON) scripts/verify_degree_five_torus_module.py
	$(PYTHON) scripts/search_rees_torsion_witnesses.py --max-target-degree 0
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
	$(PYTHON) scripts/verify_generic_affine_mark_faithfulness.py
	$(PYTHON) scripts/verify_hasse_typical_seed_recovery.py
	$(PYTHON) scripts/verify_degree_six_gaussian_moment_geometry.py
	$(PYTHON) scripts/verify_moment_prony_determinantal_geometry.py
	$(PYTHON) scripts/verify_hessian_ritt_degree_six.py
	$(PYTHON) scripts/verify_degree_six_ritt_atlas.py
	$(PYTHON) scripts/verify_degree_six_ritt_boundary_atlas.py
	$(PYTHON) scripts/verify_hessian_ritt_degrees_eight_twelve.py
	$(PYTHON) scripts/verify_hessian_synchronization_lifts.py
	$(PYTHON) scripts/verify_degree30_ritt_2_complex.py
	$(PYTHON) scripts/verify_stable_generator_rigidity.py
	$(PYTHON) scripts/verify_multicluster_ll_comparison.py
	$(PYTHON) scripts/verify_labelled_node_saturation.py
	$(PYTHON) scripts/verify_branch_wonderful_pullback.py
	$(PYTHON) scripts/verify_source_vertex_rigidity.py
	$(PYTHON) scripts/verify_general_radial_source_atlas.py
	$(PYTHON) scripts/verify_polynomial_monodromy_forests.py
	$(PYTHON) scripts/verify_monodromy_inertia_characters.py
	$(PYTHON) scripts/verify_recursive_resonance_atlas.py
	$(PYTHON) scripts/verify_branch_scale_fan.py
	$(PYTHON) scripts/verify_degree_six_branch_target_graph.py
	$(PYTHON) scripts/verify_degree_six_admissible_equal_scale.py
	$(PYTHON) scripts/verify_degree_six_admissible_radial_atlas.py
	$(PYTHON) scripts/verify_degree_six_admissible_maxwell_atlas.py
	$(PYTHON) scripts/verify_degree_six_central_hurwitz_selection.py
	$(PYTHON) scripts/verify_degree_six_stack_inertia.py
	$(PYTHON) scripts/verify_degree_six_stacky_fan_descent.py
	$(PYTHON) scripts/verify_rerooting_groupoid_boundary.py
	$(PYTHON) scripts/verify_coarse_affine_mark_descent.py
	$(PYTHON) scripts/verify_restricted_ll_degree.py
	$(PYTHON) scripts/verify_caustic_maxwell_boundary.py
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
	$(PYTHON) scripts/verify_all_degree_rational_fibers.py
	$(PYTHON) scripts/verify_real_fiber_spectrum.py
	$(PYTHON) scripts/verify_weighted_marked_root_model.py
	$(SYSTEM_PYTHON) scripts/audit_weighted_independent.py

verify-foundations-formal: verify-foundations verify-lean-foundational

# Optional independent bounded-degree comparison with Macaulay2's classical
# CoincidentRootLoci package.  The wrapper uses a pinned Docker image if M2
# is not installed locally.
verify-coincident-root-loci:
	bash scripts/verify_coincident_root_slices.sh

verify-papers:
	@set -e; for paper in papers/*/main.tex; do \
		latexmk -cd -pdf -interaction=nonstopmode -halt-on-error "$$paper"; \
	done

clean-papers:
	@set -e; for paper in papers/*/main.tex; do \
		latexmk -cd -C "$$paper"; \
	done
	$(RM) papers/core-counterexample/main.bbl

verify: check verify-plane-jc verify-core verify-theorems verify-regressions verify-derived

verify-logged:
	mkdir -p artifacts/verification
	$(PYTHON) scripts/record_environment.py | tee artifacts/verification/environment.txt
	set -o pipefail; $(MAKE) verify 2>&1 | tee artifacts/verification/verify.log
