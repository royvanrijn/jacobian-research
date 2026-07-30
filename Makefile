SHELL := /bin/bash

PYTHON ?= .venv/bin/python
SYSTEM_PYTHON ?= python3

FINALIZED_PAPERS := \
	papers/gaussian-moments-two-variables \
	papers/sparse-minimality-gaussian-moments-dimension-three
ACTIVE_PAPERS := \
	papers/common-arithmetic-fibers \
	papers/fixed-map-hasse-failures
PARKED_PAPERS := \
	papers/exact-real-chamber-spectra \
	papers/discriminant-pencils
COMPANION_PAPERS := papers/quadratic-gauge-nonproperness
VERIFIED_PAPERS := $(FINALIZED_PAPERS) $(ACTIVE_PAPERS)
ALL_PAPERS := $(VERIFIED_PAPERS) $(PARKED_PAPERS) $(COMPANION_PAPERS)

.PHONY: check verify verify-logged verify-minimal verify-core verify-geometry \
	verify-theorems verify-regressions verify-derived verify-family \
	verify-external-consequences verify-restricted-minima verify-two-real-gmc verify-sic2c4 verify-factorial-moments verify-factorial-frontier verify-counterexample-scoreboard verify-plane-jc verify-plane-case2-residue-strata verify-plane-case2-j1-endpoint verify-plane-case2-maximal-gcd verify-plane-case2-gcd6 verify-plane-poisson-radical verify-plane-poisson-primary-charts verify-plane-poisson-separators verify-plane-poisson-primary-filtration verify-plane-poisson-filtered-modules verify-weighted-boundary verify-quartic-degree-drop-quantization \
	verify-lr-rees-sagbi \
	verify-plane-sparse-supports verify-plane-support-bridge \
	verify-linear-torus-free verify-algebraic-torus-free \
	verify-master \
	verify-quartic verify-normal-forms verify-formal verify-lean-local \
	verify-lean-foundational \
	verify-foundations verify-foundations-formal \
	verify-gq2-local-fibers \
	verify-coincident-root-loci verify-papers verify-ritt-boundary \
	verify-ritt-2-complex verify-ll-ritt-reduction verify-ritt-deformation-complex \
	verify-unified-deformation-complex \
	verify-boundary-obstruction-theory \
	verify-universal-cubic-cotangent-saturation \
	verify-cubic-formal-gauge-cokernel-atlas \
	verify-nodal-cubic-formal-slice \
	verify-universal-cubic-filtered-syzygy-frontier \
	verify-marked-root-ore-bridge \
	verify-degree-five-relative-quantization-family \
	verify-degree-five-cubic-h7-unit \
	verify-relative-fiber-connection-complex \
	verify-degree42-ritt-relative-cone \
	verify-hessian-synchronization \
	verify-common-right-factor-synchronization \
	verify-degree42-hessian-normal-jets \
	verify-degree42-conormal-rees-synchronization \
	verify-degree42-divisor-rees-reduction \
	verify-degree42-kuranishi-branches \
	verify-degree42-discriminant-quartics \
	verify-degree42-depth-reduction \
	verify-degree42-ab-residual-quartics \
	verify-degree42-ab-residual-factors \
	verify-degree42-ab-residual-quintics \
	verify-degree42-higher-gcd-strata \
	verify-degree42-kuranishi-cutoff-chain \
	verify-backward-cubic-reduction \
	verify-macfarlane-f12 \
	verify-k12-coordinate-pair-frontier \
	verify-k12-parameterized-completion \
	verify-k12-z8-cubic-completion \
	verify-k12-single-defect-quartic-completion \
	verify-hvc38-cross-frontier \
	verify-hvc38-gap-closure \
	verify-hvc38-maximal-block-closure \
	verify-support-saturation-compiler refresh-support-saturation-cases \
	verify-degree30-hessian-pairs refresh-degree30-hessian-pairs \
	verify-contact-r6 verify-contact-branch-schema verify-contact-r7-asymptotic \
	verify-contact-r8-asymptotic verify-contact-r8 \
	verify-parameter-dusart-frontier verify-parameter-sharp-dusart-frontier \
	verify-parameter-adaptive-dusart-frontier \
	verify-minimal-boundary verify-minimal-boundary-pipeline \
	render-status prepare-arxiv-uploads clean-papers

.PHONY: verify-normal-covering-hasse verify-arithmetic-compilation \
	refresh-arithmetic-compilation \
	verify-common-arithmetic-fibers-correspondence \
	refresh-common-arithmetic-fibers-example \
	verify-fixed-map-hasse-failures \
	refresh-multiplicative-hasse-artifact
.PHONY: verify-hilbert14-invariants
.PHONY: verify-rank-three-collision-descent
.PHONY: verify-rank-four-collision-cross-ratio
.PHONY: verify-all-rank-collision-projective-descent
.PHONY: verify-rank-four-nonprojective-keller-lift
.PHONY: verify-rank-four-degree-eighteen-target-obstruction

check:
	$(PYTHON) -m compileall -q jcsearch scripts
	$(PYTHON) scripts/check_markdown_links.py
	$(PYTHON) scripts/audit_status.py

verify-rank-three-collision-descent:
	$(PYTHON) scripts/verify_rank_three_collision_descent.py

verify-rank-four-collision-cross-ratio:
	$(PYTHON) scripts/verify_rank_four_collision_cross_ratio.py

verify-all-rank-collision-projective-descent:
	$(PYTHON) scripts/verify_all_rank_collision_projective_descent.py

verify-rank-four-nonprojective-keller-lift:
	$(PYTHON) scripts/verify_rank_four_nonprojective_keller_lift.py

verify-rank-four-degree-eighteen-target-obstruction:
	$(PYTHON) scripts/verify_rank_four_degree_eighteen_target_obstruction.py

verify-backward-cubic-reduction:
	$(PYTHON) scripts/verify_backward_cubic_suite.py

verify-macfarlane-f12:
	$(PYTHON) scripts/verify_macfarlane_f12_suite.py

verify-k12-coordinate-pair-frontier:
	$(PYTHON) scripts/audit_k12_coordinate_pair_frontier.py

verify-k12-parameterized-completion:
	$(PYTHON) scripts/audit_k12_parameterized_completion.py

verify-k12-z8-cubic-completion:
	$(PYTHON) scripts/audit_k12_z8_cubic_completion.py

verify-k12-single-defect-quartic-completion:
	$(PYTHON) scripts/audit_k12_single_defect_quartic_completion.py

verify-hvc38-cross-frontier:
	$(PYTHON) scripts/audit_hvc38_cross_construction_frontier.py

verify-hvc38-gap-closure:
	$(PYTHON) scripts/audit_hvc38_gap_closure.py

verify-hvc38-maximal-block-closure:
	$(PYTHON) scripts/audit_hvc38_maximal_block_closure.py

verify-support-saturation-compiler:
	$(PYTHON) scripts/verify_support_saturation_compiler.py

refresh-support-saturation-cases:
	$(PYTHON) scripts/compile_support_saturation_cases.py

verify-gq2-local-fibers:
	$(SYSTEM_PYTHON) scripts/verify_gq2_permutation_action.py arithmetic/certificates/gq2_s3_x3_minus_2.json
	$(SYSTEM_PYTHON) scripts/verify_gq2_permutation_action.py arithmetic/certificates/gq2_s4_mixed_action.json
	$(SYSTEM_PYTHON) scripts/verify_gq2_permutation_action.py arithmetic/certificates/gq2_common_quintic_stable_pair.json
	$(PYTHON) scripts/verify_gq2_action_first_keller.py
	$(PYTHON) scripts/verify_gq2_s4_quartic_keller.py
	$(PYTHON) scripts/verify_marked_q2_stable_separation.py
	gp -q scripts/verify_gq2_s4_local_models.gp
	gp -q scripts/verify_gq2_local_decompositions.gp

verify-normal-covering-hasse:
	$(PYTHON) scripts/verify_normal_covering_certificates.py
	$(SYSTEM_PYTHON) scripts/verify_banks_degree_5_10_candidates.py
	$(PYTHON) scripts/verify_degree_six_normal_cover_keller.py

verify-arithmetic-compilation:
	$(SYSTEM_PYTHON) scripts/verify_arithmetic_keller_certificate.py
	gp -q -f scripts/verify_arithmetic_keller_certificate.gp
	$(SYSTEM_PYTHON) scripts/verify_arithmetic_keller_certificate.py artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json
	ARITHMETIC_CERTIFICATE=artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json gp -q -f scripts/verify_arithmetic_keller_certificate.gp
	$(SYSTEM_PYTHON) scripts/verify_arithmetic_keller_certificate.py artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json
	ARITHMETIC_CERTIFICATE=artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json gp -q -f scripts/verify_arithmetic_keller_certificate.gp
	cd formal/finite-etale-keller && lake env lean FiniteEtaleKeller/GeneratedArithmeticQuintic.lean
	cd formal/finite-etale-keller && lake env lean FiniteEtaleKeller/GeneratedArithmeticQuinticStableM2.lean
	cd formal/finite-etale-keller && lake env lean FiniteEtaleKeller/GeneratedArithmeticCubicStableN7.lean

refresh-arithmetic-compilation:
	$(PYTHON) scripts/compile_arithmetic_keller_certificate.py
	$(PYTHON) scripts/compile_arithmetic_keller_certificate.py --stable-parameter 2 --certificate artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json --lean-module FiniteEtaleKeller.GeneratedArithmeticQuinticStableM2 --lean formal/finite-etale-keller/FiniteEtaleKeller/GeneratedArithmeticQuinticStableM2.lean
	$(PYTHON) scripts/compile_arithmetic_keller_certificate.py --spec arithmetic/specifications/connected_cubic_stable_n7.json --certificate artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json --lean formal/finite-etale-keller/FiniteEtaleKeller/GeneratedArithmeticCubicStableN7.lean

verify-common-arithmetic-fibers-correspondence:
	$(PYTHON) scripts/verify_common_arithmetic_fibers_correspondence.py

refresh-common-arithmetic-fibers-example:
	$(PYTHON) scripts/compile_common_arithmetic_fibers_example.py

verify-fixed-map-hasse-failures:
	$(PYTHON) scripts/verify_infinite_hasse_keller_fibers.py
	$(PYTHON) scripts/verify_multiplicative_hasse_artifact.py

refresh-multiplicative-hasse-artifact:
	$(PYTHON) scripts/count_multiplicative_hasse_parameters.py \
		--bound 1000000 \
		--output artifacts/generated-results/multiplicative_hasse_parameters_1000000.json

render-status:
	$(SYSTEM_PYTHON) scripts/render_status.py

verify-minimal:
	$(SYSTEM_PYTHON) scripts/verify_counterexample_independent.py

verify-plane-jc:
	$(SYSTEM_PYTHON) plane-jc/cas/frontier_125_150.py
	$(PYTHON) plane-jc/cas/boundary_lattice_prefilter.py
	$(PYTHON) plane-jc/cas/test_intrinsic_a2_boundary.py
	$(PYTHON) plane-jc/cas/test_plane_boundary_exclusion.py
	$(PYTHON) plane-jc/cas/test_degree_zero_endpoint_pairing.py
	$(PYTHON) plane-jc/cas/test_finite_normalization_signatures.py
	$(PYTHON) plane-jc/cas/test_target_conductor_atlas.py
	$(PYTHON) plane-jc/cas/verify_unibranch_spectator_models.py
	$(PYTHON) plane-jc/cas/test_log_boundary_compiler.py
	$(PYTHON) plane-jc/cas/test_poisson_square_rigidity.py
	$(PYTHON) plane-jc/cas/test_poisson_square_filtered_modules.py

verify-plane-sparse-supports:
	$(PYTHON) plane-jc/cas/verify_sparse_support_exclusions.py

verify-plane-support-bridge:
	$(PYTHON) plane-jc/cas/verify_affine_support_newton_bridge.py
	$(PYTHON) plane-jc/cas/classify_f2_75_125_layers.py
	$(PYTHON) plane-jc/cas/audit_f2_75_125_boundary_handoff.py
	$(PYTHON) plane-jc/cas/test_f2_75_125_frontend.py

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

verify-linear-torus-free:
	$(PYTHON) scripts/verify_no_linear_torus_counterexample.py
	$(SYSTEM_PYTHON) scripts/audit_no_linear_torus_counterexample_independent.py

verify-algebraic-torus-free:
	$(PYTHON) scripts/verify_quartic_algebraic_torus_obstruction.py

verify-ritt-boundary:
	bash scripts/verify_degree_six_ritt_boundary_atlas.sh

verify-ritt-2-complex:
	$(PYTHON) scripts/verify_degree30_ritt_2_complex.py

verify-ll-ritt-reduction:
	$(PYTHON) scripts/verify_ll_ritt_reduction.py

verify-ritt-deformation-complex:
	$(PYTHON) scripts/verify_hessian_ritt_deformation_complex.py

verify-unified-deformation-complex:
	$(PYTHON) scripts/verify_unified_deformation_complex.py

verify-boundary-obstruction-theory:
	$(PYTHON) scripts/verify_boundary_obstruction_theory.py

verify-universal-cubic-cotangent-saturation:
	PYTHONPATH=scripts $(PYTHON) scripts/verify_universal_cubic_cotangent_saturation.py

verify-cubic-formal-gauge-cokernel-atlas:
	PYTHONPATH=scripts $(PYTHON) scripts/verify_cubic_formal_gauge_cokernel_atlas.py

verify-nodal-cubic-formal-slice:
	PYTHONPATH=scripts $(PYTHON) scripts/verify_nodal_cubic_formal_slice.py

verify-universal-cubic-filtered-syzygy-frontier:
	$(PYTHON) scripts/verify_universal_cubic_filtered_syzygy_frontier.py

verify-marked-root-ore-bridge:
	$(PYTHON) scripts/verify_marked_root_ore_bridge.py

verify-degree-five-relative-quantization-family:
	PYTHONPATH=scripts $(PYTHON) scripts/verify_degree_five_relative_quantization_family.py

verify-degree-five-cubic-h7-unit:
	$(PYTHON) scripts/verify_degree_five_cubic_h7_unit_certificate.py

verify-relative-fiber-connection-complex:
	$(PYTHON) scripts/verify_relative_fiber_connection_complex.py

verify-degree42-ritt-relative-cone:
	$(PYTHON) scripts/verify_degree42_ritt_relative_cotangent_cone.py

verify-hessian-synchronization:
	$(PYTHON) scripts/verify_hessian_synchronization_lifts.py

verify-common-right-factor-synchronization:
	$(PYTHON) scripts/verify_common_right_factor_synchronization.py

verify-degree42-hessian-normal-jets:
	$(PYTHON) scripts/verify_degree42_transported_27_normal_jets.py

verify-degree42-conormal-rees-synchronization:
	$(PYTHON) scripts/verify_degree42_conormal_rees_synchronization.py

verify-degree42-divisor-rees-reduction:
	$(PYTHON) scripts/verify_degree42_divisor_rees_reduction.py

verify-degree42-kuranishi-branches:
	$(PYTHON) scripts/verify_degree42_kuranishi_branches.py

verify-degree42-discriminant-quartics:
	$(PYTHON) scripts/verify_degree42_discriminant_quartics.py

verify-degree42-depth-reduction:
	$(PYTHON) scripts/verify_degree42_depth_certificate.py --method pivots

verify-degree42-ab-residual-quartics:
	$(PYTHON) scripts/verify_degree42_ab_residual_quartics.py

verify-degree42-ab-residual-factors:
	$(PYTHON) scripts/verify_degree42_ab_residual_factors.py

verify-degree42-ab-residual-quintics:
	$(PYTHON) scripts/verify_degree42_ab_residual_quintics.py

verify-degree42-higher-gcd-strata:
	$(PYTHON) scripts/verify_degree42_higher_gcd_strata.py

verify-degree42-kuranishi-cutoff-chain:
	PYTHON=$(PYTHON) scripts/verify_degree42_kuranishi_cutoff_chain.sh

verify-degree30-hessian-pairs:
	$(PYTHON) scripts/verify_cubic_remainder_synchronization.py
	$(PYTHON) scripts/verify_degree30_transported_23_synchronization.py
	$(PYTHON) scripts/verify_degree30_transported_25_synchronization.py
	$(PYTHON) scripts/verify_degree30_transported_35_synchronization.py

refresh-degree30-hessian-pairs:
	$(PYTHON) scripts/verify_cubic_remainder_synchronization.py
	$(PYTHON) scripts/verify_degree30_transported_23_synchronization.py --refresh
	$(PYTHON) scripts/verify_degree30_transported_25_synchronization.py --refresh
	$(PYTHON) scripts/verify_degree30_transported_35_synchronization.py --refresh

audit-degree30-hessian-synchronization-pairs:
	$(PYTHON) scripts/audit_degree30_hessian_synchronization_pairs.py

verify-minimal-boundary-pipeline:
	$(PYTHON) scripts/verify_minimal_boundary_pipeline.py

verify-minimal-boundary:
	$(PYTHON) scripts/verify_minimal_boundary_cubic.py
	$(PYTHON) scripts/verify_cubic_marking_frontier.py
	$(PYTHON) scripts/verify_cubic_normalization_frontend.py
	$(PYTHON) scripts/verify_cubic_symbol_double_saturation.py
	$(PYTHON) scripts/verify_cubic_symbol_deformation_saturation.py
	Singular -q scripts/verify_cubic_double_saturation.sing
	$(PYTHON) scripts/verify_cubic_gauge_straightening.py
	$(PYTHON) scripts/verify_cubic_gauge_first_obstruction.py

verify-cubic-quartic-tangent-saturation:
	$(PYTHON) scripts/verify_cubic_symbol_quartic_tangent_saturation.py

verify-smooth-cubic-quartic-plane-saturation:
	$(PYTHON) scripts/verify_smooth_cubic_quartic_plane_saturation.py

verify-singular-cubic-quartic-plane-saturation:
	$(PYTHON) scripts/verify_singular_cubic_quartic_plane_saturation.py

verify-smooth-cubic-quartic-three-space-saturation:
	$(PYTHON) scripts/verify_smooth_cubic_quartic_three_space_saturation.py

verify-cubic-dense-quartic-plane-saturation:
	$(PYTHON) scripts/verify_cubic_symbol_dense_quartic_plane_saturation.py

verify-cubic-affine-dense-quartic-plane-saturation:
	$(PYTHON) scripts/verify_cubic_symbol_affine_dense_quartic_plane_saturation.py

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

verify-contact-r8-asymptotic:
	$(PYTHON) scripts/verify_contact_resultant_r8_asymptotic.py

verify-contact-r8:
	$(PYTHON) scripts/verify_contact_resultant_r8_effective.py

verify-core: verify-minimal
	$(PYTHON) scripts/verify_counterexample.py
	$(PYTHON) scripts/audit_map_consistency.py
	$(PYTHON) scripts/verify_normalized_factorization_slice.py
	$(PYTHON) scripts/verify_quadratic_cubic_factorization_invariants.py
	$(PYTHON) scripts/verify_quadratic_cubic_modification_topology.py
	$(PYTHON) scripts/verify_quadratic_cubic_additive_actions.py
	$(PYTHON) scripts/verify_quadratic_cubic_saturated_lnd.py
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

verify-hilbert14-invariants:
	$(PYTHON) scripts/verify_quadratic_cubic_saturated_lnd.py
	$(PYTHON) scripts/verify_quadratic_quartic_additive_invariants.py
	$(PYTHON) scripts/verify_hilbert14_saturation_ladders.py
	$(PYTHON) scripts/verify_multiboundary_hilbert14_antichain.py

verify-theorems:
	$(PYTHON) scripts/verify_controlled_boundary_suspensions.py
	$(PYTHON) scripts/verify_puncture_rank_frontier.py
	$(PYTHON) scripts/verify_three_puncture_nonlinear_frontier.py
	$(MAKE) verify-master
	$(MAKE) verify-minimal-boundary
	$(PYTHON) scripts/verify_weighted_seed_schema.py
	$(PYTHON) scripts/verify_weighted_seed_theorem.py
	$(PYTHON) scripts/verify_all_degree_rational_fibers.py
	$(PYTHON) scripts/verify_common_arithmetic_fibers.py
	$(PYTHON) scripts/verify_universal_quartic_fiber_multiplicity.py
	$(PYTHON) scripts/verify_universal_quartic_gauge_multiplicity.py
	$(PYTHON) scripts/verify_universal_cubic_gauge_multiplicity.py
	$(PYTHON) scripts/verify_universal_power_shifted_gauge_multiplicity.py
	$(PYTHON) scripts/verify_whole_plane_stable_multiplicity.py
	$(PYTHON) scripts/verify_universal_quintic_fiber_multiplicity.py
	$(PYTHON) scripts/verify_universal_higher_degree_fiber_multiplicity.py
	$(PYTHON) scripts/verify_universal_multiplicity_witness_cards.py
	$(PYTHON) scripts/verify_low_rank_multiplicity_boundaries.py
	$(PYTHON) scripts/verify_real_fiber_spectrum.py
	$(PYTHON) scripts/verify_hasse_keller_fiber.py
	$(PYTHON) scripts/verify_infinite_hasse_keller_fibers.py
	$(PYTHON) scripts/verify_multiplicative_hasse_artifact.py
	$(PYTHON) scripts/verify_weighted_marked_root_model.py
	$(PYTHON) scripts/verify_intrinsic_selector_attack.py
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
	$(MAKE) verify-linear-torus-free
	$(MAKE) verify-algebraic-torus-free
	$(PYTHON) scripts/verify_quadratic_cancellation_intersection.py
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
	$(PYTHON) scripts/verify_index_three_degree_bound_counterexample.py
	$(PYTHON) scripts/derive_index_three_tree_obstruction.py
	$(PYTHON) scripts/verify_index_three_rank_normal_form.py
	$(PYTHON) scripts/verify_restricted_minima_frontier.py

verify-two-real-gmc:
	$(PYTHON) scripts/verify_two_real_gmc_frontier.py
	$(PYTHON) scripts/verify_two_real_gmc_symmetric_chart.py
	$(PYTHON) scripts/verify_two_real_gmc_remaining_four_weight.py
	$(PYTHON) scripts/verify_two_real_gmc_five_weight_frontier.py
	$(PYTHON) scripts/verify_cubic_gaussian_null_cone_closure.py
	$(PYTHON) scripts/audit_prime_endpoint_rigidity_independent.py
	$(PYTHON) scripts/audit_two_real_gmc_lower_face.py
	$(PYTHON) scripts/verify_two_real_gmc_unit_star_rigidity.py
	$(PYTHON) scripts/verify_two_real_gmc_three_level_rigidity.py
	$(PYTHON) scripts/verify_two_real_gmc_first_cycle_rigidity.py
	$(PYTHON) scripts/verify_two_real_gmc_three_weight_low_degree.py
	$(PYTHON) scripts/verify_two_real_gmc_resolvent_system.py

verify-sic2c4:
	$(PYTHON) scripts/verify_two_pair_image_mathieu_counterexample.py
	$(SYSTEM_PYTHON) scripts/audit_two_pair_image_mathieu_coefficient_extraction.py
	$(PYTHON) scripts/verify_two_pair_sic_characteristic_p.py
	cd formal/finite-etale-keller && lake build FiniteEtaleKeller.SIC2C4FiniteSum

verify-counterexample-scoreboard: verify-two-real-gmc verify-sic2c4
	$(PYTHON) scripts/verify_three_real_gmc_rank_one_classification.py
	$(PYTHON) scripts/audit_dvorsky_gvc5_counterexample.py
	$(PYTHON) scripts/verify_separable_gvc_escape_obstructions.py
	$(PYTHON) scripts/verify_binary_heat_quadratic_gvc.py
	$(PYTHON) scripts/verify_minimal_counterexample_scoreboard.py

verify-factorial-moments:
	$(SYSTEM_PYTHON) scripts/verify_factorial_moment_witnesses.py

verify-factorial-frontier:
	$(PYTHON) scripts/verify_sparse_factorial_moment_frontier.py

verify-lr-rees-sagbi:
	$(PYTHON) scripts/compute_lr_rees_sagbi_modules.py
	$(SYSTEM_PYTHON) scripts/audit_lr_rees_sagbi_module_certificate.py

verify-regressions: verify-external-consequences verify-factorial-moments verify-factorial-frontier
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
	$(PYTHON) scripts/verify_ll_ritt_reduction.py
	$(PYTHON) scripts/verify_hessian_ritt_deformation_complex.py
	$(PYTHON) scripts/verify_unified_deformation_complex.py
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

verify-quartic-degree-drop-quantization:
	PYTHONPATH=scripts $(PYTHON) scripts/verify_quartic_degree_drop_quantization.py

verify-normal-forms:
	$(PYTHON) scripts/cubic_homogeneous_reduction.py
	$(PYTHON) scripts/verify_cubic_homogeneous_counterexample.py
	$(PYTHON) scripts/cubic_linear_reduction.py
	$(PYTHON) scripts/verify_cubic_linear_counterexample.py
	$(PYTHON) scripts/audit_stable_normal_form_independent.py
	$(PYTHON) scripts/generate_stable_normal_form_consequences.py

verify-derived: verify-normal-forms verify-boundary-obstruction-theory

# Optional formal replication. This fetches Dean Cureton's separately authored
# Lean project at the audited commit recorded in verified/LEAN_FOUNDATIONAL_MAP.md; it is kept
# out of the default target because it downloads a pinned Lean/mathlib toolchain.
verify-lean-local:
	$(SYSTEM_PYTHON) scripts/check_lean_placeholders.py
	$(SYSTEM_PYTHON) scripts/check_paper_certificate_imports.py
	cd formal/discriminant-pencils && lake exe cache get
	cd formal/discriminant-pencils && lake build
	cd formal/finite-etale-keller && lake exe cache get
	cd formal/finite-etale-keller && lake build
	cd formal/gmc2 && lake exe cache get
	cd formal/gmc2 && lake build

verify-lean-foundational:
	bash scripts/verify_lean_foundational_map.sh

verify-formal: verify-lean-local verify-lean-foundational

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
	@set -e; for paper_dir in $(VERIFIED_PAPERS); do \
		latexmk -cd -pdf -interaction=nonstopmode -halt-on-error "$$paper_dir/main.tex"; \
	done
	mkdir -p output/pdf
	cp papers/gaussian-moments-two-variables/main.pdf \
		output/pdf/gaussian-moments-two-variables.pdf
	cp papers/sparse-minimality-gaussian-moments-dimension-three/main.pdf \
		output/pdf/sparse-minimality-gaussian-moments-dimension-three.pdf
	cp papers/common-arithmetic-fibers/main.pdf \
		output/pdf/common-arithmetic-fibers.pdf
	cp papers/fixed-map-hasse-failures/main.pdf \
		output/pdf/fixed-map-hasse-failures.pdf

prepare-arxiv-uploads:
	bash scripts/prepare_arxiv_uploads.sh

clean-papers:
	@set -e; for paper_dir in $(ALL_PAPERS); do \
		latexmk -cd -C "$$paper_dir/main.tex"; \
	done
	$(RM) papers/core-counterexample/main.bbl

verify: check verify-plane-jc verify-core verify-theorems verify-regressions verify-derived

verify-logged:
	mkdir -p artifacts/verification
	$(PYTHON) scripts/record_environment.py | tee artifacts/verification/environment.txt
	set -o pipefail; $(MAKE) verify 2>&1 | tee artifacts/verification/verify.log
