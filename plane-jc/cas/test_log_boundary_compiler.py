#!/usr/bin/env python3
"""Regression tests for the certified Newton/log-boundary front end."""

import sympy as sp

from intrinsic_a2_boundary import (
    infer_finite_model_dicritical_projection_budget,
)
from log_boundary_compiler import (
    BranchScale,
    NewtonBoundaryCertificate,
    PullbackMonomial,
    ToroidalCluster,
    branch_scale_regression_certificate,
    compile_log_boundary,
    corner_direction,
    frontier_72_108_case1_boundary_package_audit,
    frontier_72_108_case2_gcd_six_audit,
    frontier_72_108_case2_maximal_gcd_audit,
    frontier_72_108_case2_j1_endpoint_audit,
    frontier_72_108_case2_lower_jet_audit,
    frontier_72_108_case2_boundary_package_audit,
    frontier_72_108_case2_residue_stratum_audit,
    frontier_72_108_common_graph_pole_audit,
    frontier_72_108_dicritical_residue_degree_audits,
    frontier_72_108_exceptional_pole_audit,
    frontier_72_108_forced_first_block_cluster_audit,
    frontier_72_108_incomplete_certificate,
    frontier_72_108_local_report,
    frontier_72_108_minimal_dicritical_extension_audit,
    frontier_72_108_plane_return_edge_audit,
    frontier_72_108_plane_return_partition_audits,
    frontier_72_108_poisson_ramification_audits,
    frontier_72_108_residue_case_tree,
    frontier_72_108_translation_records,
    frontier_72_108_two_step_dicritical_witnesses,
    frontier_72_108_unselected_factor_audit,
    hirzebruch_transition_audit,
    fill_temporary_boundary,
    laurent_translation_base_ideal_audit,
    laurent_translation_branch_certificate,
    laurent_translation_composition_audit,
    laurent_translation_graph_certificate,
)


compiled = compile_log_boundary(branch_scale_regression_certificate())
assert compiled.passes_prefilter
assert compiled.intrinsic_a2_audit is not None
assert compiled.intrinsic_a2_audit.passes
assert compiled.as_dict()["intrinsic_a2_boundary"]["passes"]
assert compiled.boundary.names == ("L", "E1", "E2", "E3")
assert tuple(step.center for step in compiled.blowups) == (
    ("L",),
    ("E1",),
    ("E1", "E2"),
)
assert tuple(step.ray for step in compiled.blowups) == (
    (1, 1),
    (1, 2),
    (2, 3),
)
assert compiled.boundary.class_matrix == sp.Matrix(
    [
        [1, 0, 0, 0],
        [-1, 1, 0, 0],
        [0, -1, 1, 0],
        [0, -1, -1, 1],
    ]
)
assert abs(int(compiled.intersection_matrix.det())) == 1
scale = compiled.clusters[0].scales[0]
assert scale.equality_ray == (2, 3)
assert scale.dicritical_divisor == "E3"
assert scale.common_base_order == 6
assert scale.residue_degree == 1
assert scale.projection_differents == (1, 2)
assert scale.semigroup_conductor == 2
valuation = {
    entry.divisor: entry for entry in compiled.clusters[0].valuations
}["E3"]
assert valuation.pullback_orders == (("target_u", 2), ("target_v", 3))
assert valuation.tame_differents == (("target_u", 1), ("target_v", 2))


# Two scales at one infinitely-near point share their fan prefix rather than
# creating two disconnected copies of the same preliminary blowups.
merged = compile_log_boundary(
    NewtonBoundaryCertificate(
        name="merged scales",
        chart="A2",
        clusters=(
            ToroidalCluster(
                name="p",
                root_boundary="L",
                scales=(
                    BranchScale("left", 3, 2),   # ray (2,3)
                    BranchScale("right", 2, 3),  # ray (3,2)
                ),
                pullbacks=(PullbackMonomial("boundary", 1, 1),),
            ),
        ),
        transformations=("one common toroidal chart",),
        theorem_source="regression fixture",
        exhaustive=True,
    )
)
assert tuple(step.ray for step in merged.blowups) == (
    (1, 1),
    (2, 1),
    (3, 2),
    (1, 2),
    (2, 3),
)
assert len(merged.boundary.names) == 6
assert merged.localization.smith_diagonal == (1, 1, 1, 1, 1, 1)


# A later resonance chart may be rooted on an exceptional prime created by an
# earlier radial chart.  Its coordinate orders are local to that nested chart.
nested = compile_log_boundary(
    NewtonBoundaryCertificate(
        name="nested resonance scale",
        chart="A2",
        clusters=(
            ToroidalCluster(
                name="radial",
                root_boundary="L",
                scales=(BranchScale("first", 1, 1),),
            ),
            ToroidalCluster(
                name="resonance",
                root_boundary="E1",
                scales=(BranchScale("second", 2, 1),),
                local_parameters=("epsilon", "cross_ratio"),
            ),
        ),
        transformations=("radial chart", "second-scale resonance chart"),
        theorem_source="regression fixture",
        exhaustive=True,
    )
)
assert tuple(step.center for step in nested.blowups) == (
    ("L",),
    ("E1",),
    ("E2",),
)
assert nested.clusters[1].local_parameters == ("epsilon", "cross_ratio")
assert nested.clusters[1].scales[0].dicritical_divisor == "E3"


# A bounded Farey grid exercises the general Stern--Brocot insertion rather
# than only the two hand-written rays above.
farey_scales = tuple(
    BranchScale(f"s_{a}_{b}", a, b)
    for a in range(1, 6)
    for b in range(1, 6)
    if sp.gcd(a, b) == 1
)
farey = compile_log_boundary(
    NewtonBoundaryCertificate(
        name="bounded Farey fan",
        chart="A2",
        clusters=(
            ToroidalCluster(
                name="grid",
                root_boundary="L",
                scales=farey_scales,
            ),
        ),
        transformations=("one regular toroidal chart",),
        theorem_source="bounded algorithm regression",
        exhaustive=True,
    )
)
fan_rays = tuple(ray for ray, _ in farey.clusters[0].rays)
assert all(
    left[0] * right[1] - left[1] * right[0] == 1
    for left, right in zip(fan_rays, fan_rays[1:])
)
assert {scale.equality_ray for scale in farey_scales} <= set(fan_rays)
assert abs(int(farey.boundary.class_matrix.det())) == 1


# Boundary blowups preserve the rank-one unit lattice of the Laurent chart.
laurent = compile_log_boundary(
    NewtonBoundaryCertificate(
        name="Laurent chart scale",
        chart="GmA1",
        clusters=(
            ToroidalCluster(
                name="x_zero",
                root_boundary="X0",
                scales=(BranchScale("scale", 2, 1),),
            ),
        ),
        transformations=("local branch-scale chart",),
        theorem_source="regression fixture",
        exhaustive=True,
    )
)
assert laurent.expected_unit_rank == 1
assert laurent.localization.unit_rank == 1
assert laurent.passes_prefilter
assert laurent.intrinsic_a2_audit is None


# Laurent translations y -> y+lambda*x^-q are resolved at the SNC crossing
# x=0, w=1/y=0 by the scale [w:x^q].  Both fan endpoints are boundary primes.
crossing = compile_log_boundary(
    NewtonBoundaryCertificate(
        name="order-three Laurent translation",
        chart="GmA1",
        clusters=(
            ToroidalCluster(
                name="x0_yinf",
                root_boundary="Yinf",
                transverse_boundary="X0",
                local_parameters=("w", "x"),
                scales=(BranchScale("w_over_x3", 1, 3),),
            ),
        ),
        transformations=("y -> y + lambda*x^-3",),
        theorem_source="GGHV 2022 Proposition 4.3 local translation",
        exhaustive=True,
    )
)
assert tuple(step.center for step in crossing.blowups) == (
    ("Yinf", "X0"),
    ("Yinf", "E1"),
    ("Yinf", "E2"),
)
assert crossing.clusters[0].scales[0].equality_ray == (3, 1)
assert crossing.clusters[0].scales[0].dicritical_divisor == "E3"
assert crossing.localization.unit_rank == 1


# The primary Case (8,28) proof supplies local Laurent translations of
# orders 2, 3, and 4.  In (w,x) parameters these are [w:x^q], with equality
# rays (q,1).  This is an audited local skeleton, not a claimed global graph.
frontier_records = frontier_72_108_translation_records()
assert tuple(record.order for record in frontier_records) == (2, 3, 4)
assert tuple(
    record.branch_scale.equality_ray for record in frontier_records
) == ((2, 1), (3, 1), (4, 1))
assert tuple(
    record.branch_scale.common_base_order for record in frontier_records
) == (2, 3, 4)
assert tuple(
    record.branch_scale.projection_differents for record in frontier_records
) == ((1, 0), (2, 0), (3, 0))
assert all(
    record.branch_scale.branch_count == 1
    and record.branch_scale.semigroup_conductor == 0
    for record in frontier_records
)
assert tuple(
    record.common_after_branching for record in frontier_records
) == (False, False, True)
frontier_local_compilations = tuple(
    compile_log_boundary(laurent_translation_branch_certificate(record))
    for record in frontier_records
)
assert tuple(
    len(compilation.blowups) for compilation in frontier_local_compilations
) == (2, 3, 4)
assert all(
    compilation.blowups[0].center == ("Yinf", "X0")
    and compilation.blowups[-1].center
    == ("Yinf", f"E{len(compilation.blowups) - 1}")
    for compilation in frontier_local_compilations
)
assert all(
    compilation.passes_prefilter
    for compilation in frontier_local_compilations
)
base_ideal_audits = tuple(
    laurent_translation_base_ideal_audit(record)
    for record in frontier_records
)
assert tuple(
    audit.branch_fan_length for audit in base_ideal_audits
) == (2, 3, 4)
assert tuple(
    audit.base_ideal_length for audit in base_ideal_audits
) == (4, 6, 8)
assert tuple(
    audit.adapted_generators for audit in base_ideal_audits
) == (("t", "x^4"), ("t", "x^6"), ("t", "x^8"))
assert all(audit.verify_identity() for audit in base_ideal_audits)
assert tuple(
    audit.target_yinf_orders for audit in base_ideal_audits
) == (
    (2, 2, 1, 0),
    (3, 3, 3, 2, 1, 0),
    (4, 4, 4, 4, 3, 2, 1, 0),
)
assert all(
    audit.target_x_orders == (1,) * (2 * audit.order)
    for audit in base_ideal_audits
)
composition_audits = tuple(
    laurent_translation_composition_audit(orders)
    for orders in ((2, 4), (3, 4), (2, 3, 4))
)
assert all(audit.verify_identity() for audit in composition_audits)
assert tuple(
    audit.adapted_generators for audit in composition_audits
) == (("t", "x^8"),) * 3
assert tuple(audit.graph_length for audit in composition_audits) == (8, 8, 8)
unselected_factor = frontier_72_108_unselected_factor_audit()
assert unselected_factor.verify_identities()
assert unselected_factor.case_presence == (
    ("a", False),
    ("b", False),
    ("c", True),
)
assert not unselected_factor.creates_additional_basepoint
assert not unselected_factor.meets_filled_x_zero
exceptional_poles = frontier_72_108_exceptional_pole_audit()
assert exceptional_poles.verify_orders()
assert exceptional_poles.p_pole_orders == (1, 1, 1, 8, 6, 4, 2, 1)
assert exceptional_poles.q_pole_orders == (0, 0, 1, 12, 9, 6, 3, 2)
assert exceptional_poles.target_infinity_orders == (
    1,
    1,
    1,
    12,
    9,
    6,
    3,
    2,
)
assert not exceptional_poles.complete_global_pullback
common_poles = frontier_72_108_common_graph_pole_audit()
assert common_poles.pole_vector == (
    1,
    24,
    1,
    1,
    1,
    12,
    9,
    6,
    3,
    2,
)
assert common_poles.geometric_degree == 427
assert common_poles.hyperplane_intersections == (
    0,
    12,
    0,
    0,
    11,
    10,
    0,
    0,
    2,
    1,
)
assert not common_poles.dicritical_candidates
assert not common_poles.passes
assert len(common_poles.failures) == 1
first_block = frontier_72_108_forced_first_block_cluster_audit()
assert first_block.crossing == ("E3", "E4")
assert first_block.local_parameters == ("s=Y", "t=1/(X*Y)")
assert first_block.laurent_coordinates == ("T=s/t", "z=1/s")
assert (first_block.u_degree, first_block.v_degree) == (7, 10)
assert first_block.root_count == 10
assert (
    first_block.first_base_order,
    first_block.first_exceptional_pole,
) == (10, 3)
assert (
    first_block.child_base_order,
    first_block.child_exceptional_pole,
    first_block.child_target_degree,
) == (1, 2, 1)
assert first_block.smooth_e3_basepoints == 0
assert first_block.boundary_audit.passes
assert first_block.pole_audit.geometric_degree == 317
assert first_block.pole_audit.hyperplane_intersections[-11:] == (
    0,
    *(1 for _ in range(10)),
)
assert not first_block.pole_audit.dicritical_candidates
case2_package = frontier_72_108_case2_boundary_package_audit()
assert case2_package.first_center == ("Yinf",)
assert case2_package.local_parameters == ("X", "z=1/Y")
assert (
    case2_package.first_base_order,
    case2_package.first_exceptional_pole,
    case2_package.second_base_order,
    case2_package.dicritical_pole,
) == (12, 12, 12, 0)
assert case2_package.dicritical_residue == "[1:C(r):G(r)]"
assert case2_package.dicritical_degree == 12
assert case2_package.dicritical_ramification_intersection == 35
assert case2_package.boundary_audit.passes
assert case2_package.pole_audit.passes
assert case2_package.pole_audit.geometric_degree == 29
assert case2_package.pole_audit.dicritical_candidates == ("T2",)
case1_package = frontier_72_108_case1_boundary_package_audit()
assert case1_package.selected_partition == (4,)
assert case1_package.final_chart_relation == "Y2=Y1-beta/X+delta"
assert case1_package.residual_point == "u=beta"
assert case1_package.adapted_parameter == "v=u-beta+delta*e*u"
assert (
    case1_package.first_base_order,
    case1_package.first_exceptional_pole,
    case1_package.second_base_order,
    case1_package.dicritical_pole,
) == (12, 12, 12, 0)
assert case1_package.dicritical_degree == 12
assert case1_package.dicritical_ramification_intersection == 35
assert case1_package.boundary_audit.passes
assert case1_package.pole_audit.passes
assert case1_package.pole_audit.geometric_degree == 29
assert case1_package.pole_audit.dicritical_candidates == ("R1_1",)

# Contracting the nine H-null boundary curves in either terminal graph gives
# the same formal finite Stein-model row.  These use the Keller class
# K+3H; the later transformed-Poisson correction [P,Q]=X^2 is deliberately
# a separate, nonconstant-Jacobian audit.
case1_finite_budgets = tuple(
    infer_finite_model_dicritical_projection_budget(
        case1_package.pole_audit, "R1_1", image_degree
    )
    for image_degree in (3, 6, 12)
)
case2_finite_budget = infer_finite_model_dicritical_projection_budget(
    case2_package.pole_audit, "T2", 12
)
assert all(budget.budgets_match for budget in case1_finite_budgets)
assert case2_finite_budget.budgets_match
assert all(
    len(budget.contraction.contracted_names) == 9
    for budget in (*case1_finite_budgets, case2_finite_budget)
)
assert all(
    budget.finite_self_intersection == sp.Rational(33, 8)
    for budget in (*case1_finite_budgets, case2_finite_budget)
)
assert tuple(
    (
        budget.target_curve_degree,
        budget.residue_degree,
        budget.target_normalization_correction,
    )
    for budget in case1_finite_budgets
) == ((3, 4, 8), (6, 2, 40), (12, 1, 110))
assert (
    case2_finite_budget.residue_degree,
    case2_finite_budget.target_normalization_correction,
) == (1, 110)

poisson_ramification = {
    audit.terminal_case: audit
    for audit in frontier_72_108_poisson_ramification_audits()
}
assert poisson_ramification["Case 1"].jacobian == "X^2"
assert poisson_ramification["Case 1"].coordinate_jacobian == (
    "1/(1+delta*e)"
)
assert poisson_ramification["Case 1"].pulled_back_jacobian == (
    "e^2*u(e,r)^2/(1+delta*e)"
)
assert (
    poisson_ramification["Case 1"].dicritical_divisor,
    poisson_ramification["Case 1"].dicritical_ramification_coefficient,
    poisson_ramification["Case 1"].dicritical_generic_ramification_index,
    poisson_ramification[
        "Case 1"
    ].dicritical_boundary_ramification_intersection,
    poisson_ramification[
        "Case 1"
    ].dicritical_total_ramification_intersection,
) == ("R1_1", 2, 3, 35, 35)
assert tuple(
    (name, value)
    for name, value in zip(
        poisson_ramification["Case 1"].boundary_names,
        poisson_ramification["Case 1"].affine_boundary_intersections,
    )
    if value
) == (("D0", 1),)
assert poisson_ramification["Case 2"].pulled_back_jacobian == "e^4*r^2"
assert (
    poisson_ramification["Case 2"].dicritical_divisor,
    poisson_ramification["Case 2"].dicritical_ramification_coefficient,
    poisson_ramification["Case 2"].dicritical_generic_ramification_index,
    poisson_ramification[
        "Case 2"
    ].dicritical_boundary_ramification_intersection,
    poisson_ramification[
        "Case 2"
    ].dicritical_total_ramification_intersection,
) == ("T2", 4, 5, 33, 35)
assert tuple(
    (name, value)
    for name, value in zip(
        poisson_ramification["Case 2"].boundary_names,
        poisson_ramification["Case 2"].affine_boundary_intersections,
    )
    if value
) == (("T2", 1),)
assert not any(
    audit.constant_keller_ramification_applicable
    for audit in poisson_ramification.values()
)
residue_degrees = {
    audit.terminal_case: audit
    for audit in frontier_72_108_dicritical_residue_degree_audits()
}
assert residue_degrees["Case 1"].possible_normalization_cover_degrees == (
    1,
    2,
    4,
)
assert residue_degrees[
    "Case 1"
].a_priori_normalization_cover_degrees == (1, 2, 4)
assert not residue_degrees[
    "Case 1"
].excluded_normalization_cover_degrees
assert residue_degrees["Case 1"].possible_image_degrees == (12, 6, 3)
assert residue_degrees["Case 1"].possible_normalization_differents == (
    0,
    2,
    6,
)
assert residue_degrees["Case 1"].homogeneous_initial_cover_degree == 4
assert residue_degrees["Case 1"].forced_normalization_cover_degree is None
assert residue_degrees[
    "Case 1"
].valuation_contributions_to_degree_29 == (3, 6, 12)
assert residue_degrees["Case 1"].remaining_valuation_degrees == (
    26,
    23,
    17,
)
assert residue_degrees[
    "Case 2"
].a_priori_normalization_cover_degrees == (1, 2, 4)
assert residue_degrees[
    "Case 2"
].excluded_normalization_cover_degrees == (2, 4)
assert residue_degrees[
    "Case 2"
].possible_normalization_cover_degrees == (1,)
assert residue_degrees["Case 2"].possible_image_degrees == (12,)
assert residue_degrees["Case 2"].possible_normalization_differents == (0,)
assert residue_degrees["Case 2"].forced_normalization_cover_degree == 1
assert residue_degrees["Case 2"].forced_image_degree == 12
assert residue_degrees[
    "Case 2"
].valuation_contributions_to_degree_29 == (5,)
assert residue_degrees["Case 2"].remaining_valuation_degrees == (
    24,
)
case2_residue_strata = frontier_72_108_case2_residue_stratum_audit()
assert case2_residue_strata.excluded_cover_degrees == (2, 4)
assert case2_residue_strata.surviving_cover_degrees == (1,)
assert case2_residue_strata.constraint_counts == ((2, 9), (4, 12))
assert case2_residue_strata.singular_input_sha256[0][1] == (
    "8bd56f701554cfe909f5c5dbdc78139e"
    "68583043dc1783ce1bffafd48da0297a"
)
assert "no J1 compatibility" in case2_residue_strata.equation_scope
assert "no J0" in case2_residue_strata.equation_scope
case2_j1_endpoint = frontier_72_108_case2_j1_endpoint_audit()
assert case2_j1_endpoint.forced_endpoint == "G_12"
assert case2_j1_endpoint.compatibility_equation_count == 7
assert case2_j1_endpoint.compatibility_term_counts == (8,) * 7
assert case2_j1_endpoint.generic_infinity_characteristic == (4, 13)
assert case2_j1_endpoint.generic_infinity_self_intersections == (
    -2,
    -2,
    -5,
    -1,
    -2,
    -2,
    -2,
)
assert case2_j1_endpoint.unit_ideal
assert "no J0" in case2_j1_endpoint.equation_scope
case2_maximal_gcd = frontier_72_108_case2_maximal_gcd_audit()
assert case2_maximal_gcd.gcd_degree == 7
assert case2_maximal_gcd.selected_remainder_degrees == (0, 1, 2)
assert case2_maximal_gcd.j0_coefficient_degree == 19
assert case2_maximal_gcd.selected_constraint_term_counts == (
    155,
    155,
    155,
    9,
)
assert case2_maximal_gcd.selected_constraint_parameter_degrees == (
    13,
    13,
    13,
    4,
)
assert case2_maximal_gcd.excluded_exact_gcd_degrees == (7,)
assert case2_maximal_gcd.surviving_exact_gcd_degrees == (1, 2, 3, 4, 5, 6)
assert case2_maximal_gcd.singular_input_sha256 == (
    "1ac0b4db7ddd0b50fcbef6c93d49c28f"
    "7f80cb4133a73b5f2158af6c78f3b069"
)
case2_gcd_six = frontier_72_108_case2_gcd_six_audit()
assert case2_gcd_six.gcd_degree == 6
assert case2_gcd_six.selected_factor_conditions == ("C'(0)=0", "H(0)=0")
assert case2_gcd_six.selected_g_remainder_degrees == (4, 5)
assert case2_gcd_six.selected_constraint_term_counts == (
    5,
    30,
    656,
    352,
    9,
)
assert case2_gcd_six.selected_constraint_parameter_degrees == (
    2,
    7,
    16,
    15,
    4,
)
assert case2_gcd_six.excluded_exact_gcd_degrees == (6,)
assert case2_gcd_six.surviving_precompatibility_gcd_degrees == (
    1,
    2,
    3,
    4,
    5,
)
case2_lower_jet = frontier_72_108_case2_lower_jet_audit()
assert case2_lower_jet.j0_conclusion == ("B=K*c", "F=K*g")
assert case2_lower_jet.multiplier_degree_bound == "deg(K)<=deg(H)+1"
assert case2_lower_jet.multiplier_origin_constraint == (
    "K=0 or t divides K"
)
assert case2_lower_jet.wronskian_gcd_factorization == (
    "C'*G''-C''*G'=H^2*(c*g'-c'*g)"
)
assert case2_lower_jet.reduced_j1_wronskian_identity == (
    "K^2*(C'*G''-C''*G')=-2*H^3*(A*g-c*E)"
)
assert case2_lower_jet.forced_origin_jets[-1] == (
    "t divides H=gcd(C',G')"
)
assert case2_lower_jet.cover_degree_gcd_bounds == (
    (1, 1),
    (2, 1),
    (4, 3),
)
assert case2_lower_jet.degree_one_gcd_orders == (
    "ord(H)=ord(C')=1",
    "ord(K)=ord(B)=1",
    "ord(E)=2",
    "ord(g)=2",
    "ord(G')=ord(F)=3",
)
assert case2_lower_jet.degree_one_gcd_excluded_order_pairs == ((2, 1),)
assert case2_lower_jet.linear_series_degree == 12
assert case2_lower_jet.degree_one_vanishing_sequences == (
    ("forced affine point t=0", (0, 2, 4)),
    ("point at infinity", (0, 4, 12)),
)
assert case2_lower_jet.degree_one_plucker_weights == (
    ("forced affine point t=0", 3),
    ("point at infinity", 13),
)
assert case2_lower_jet.degree_one_affine_wronskian == (
    "C'*G''-C''*G'=t^3*W_14(t), deg(W_14)=14"
)
assert case2_lower_jet.degree_one_remaining_plucker_weight == 14
assert case2_lower_jet.gcd_degree_residual_wronskian_degrees == (
    (1, 15),
    (2, 13),
    (3, 11),
    (4, 9),
    (5, 7),
    (6, 5),
    (7, 3),
)
plane_return_edge = frontier_72_108_plane_return_edge_audit()
assert plane_return_edge.common_factor_degree == 4
assert plane_return_edge.forced_factorization == (
    "A=a*r^2",
    "C=c*r^3",
)
assert plane_return_edge.case2_root_partition == (4,)
assert plane_return_edge.case1_root_partitions == (
    (4,),
    (3, 1),
    (2, 2),
    (2, 1, 1),
    (1, 1, 1, 1),
)
plane_return_partitions = frontier_72_108_plane_return_partition_audits()
assert tuple(audit.partition for audit in plane_return_partitions) == (
    (4,),
    (3, 1),
    (2, 2),
    (2, 1, 1),
    (1, 1, 1, 1),
)
assert all(
    audit.boundary_audit.passes
    and audit.pole_audit.passes
    and audit.pole_audit.geometric_degree == 29
    and abs(int(audit.boundary.class_matrix.det())) == 1
    for audit in plane_return_partitions
)
assert tuple(
    tuple(chain.hyperplane_degree for chain in audit.root_chains)
    for audit in plane_return_partitions
) == ((12,), (3, 3), (6, 6), (6, 3, 3), (3, 3, 3, 3))
assert tuple(
    audit.pole_audit.dicritical_candidates
    for audit in plane_return_partitions
) == (
    ("R1_1",),
    ("R1_4", "R2_4"),
    ("R1_2", "R2_2"),
    ("R1_2", "R2_4", "R3_4"),
    ("R1_4", "R2_4", "R3_4", "R4_4"),
)
root_chain_shapes = {
    chain.root_multiplicity: (
        chain.blowup_rays,
        chain.target_infinity_poles,
        chain.dicritical_residue,
        chain.normalization_cover_degree,
    )
    for audit in plane_return_partitions
    for chain in audit.root_chains
}
assert root_chain_shapes == {
    1: (((1, 1), (1, 2), (1, 3), (1, 4)), (9, 6, 3, 0), (2, 3), 1),
    2: (((1, 1), (1, 2)), (6, 0), (4, 6), 2),
    3: (
        ((1, 1), (1, 2), (2, 3), (3, 4)),
        (3, 0, 0, 0),
        (2, 3),
        1,
    ),
    4: (((1, 1),), (0,), (8, 12), 4),
}
minimal_dicritical = frontier_72_108_minimal_dicritical_extension_audit()
assert len(minimal_dicritical.tested_centers) == 19
assert minimal_dicritical.admissible_centers == (("E3",),)
assert minimal_dicritical.unique_center == ("E3",)
assert minimal_dicritical.pole_vector[-1] == 0
assert minimal_dicritical.hyperplane_intersections[-1] == 1
assert minimal_dicritical.canonical_coefficients[-1] == 0
assert minimal_dicritical.geometric_degree == 426
assert minimal_dicritical.ordinary_ramification[-1] == 0
assert minimal_dicritical.log_ramification[-1] == 1
assert minimal_dicritical.dicritical_candidates == ("D1",)
two_step_witnesses = frontier_72_108_two_step_dicritical_witnesses()
assert tuple(
    (
        witness.first_center,
        witness.first_pole,
        witness.geometric_degree,
        witness.dicritical_degree,
        witness.dicritical_canonical_coefficient,
    )
    for witness in two_step_witnesses
) == (
    ("Yinf", 12, 139, 12, 0),
    ("E4", 2, 323, 2, 0),
    ("E7", 1, 422, 1, 3),
    ("E8", 1, 425, 1, 4),
)
frontier_graph_compilations = tuple(
    compile_log_boundary(laurent_translation_graph_certificate(record))
    for record in frontier_records
)
assert tuple(
    len(compilation.blowups) for compilation in frontier_graph_compilations
) == (4, 6, 8)
for q, compilation in zip((2, 3, 4), frontier_graph_compilations):
    assert tuple(
        step.center for step in compilation.blowups[:q]
    ) == (("Yinf", "X0"),) + tuple(
        ("Yinf", f"E{index}") for index in range(1, q)
    )
    assert tuple(
        step.center for step in compilation.blowups[q:]
    ) == tuple(
        (f"E{index}",) for index in range(q, 2 * q)
    )
    assert compilation.clusters[1].root_boundary == f"E{q}"
    assert compilation.clusters[1].scales[0].equality_ray == (1, q)
    assert compilation.passes_prefilter


# Printed corners determine a signed Newton normal, not a positive toroidal
# cluster.  The published pair therefore remains a typed incomplete input.
assert corner_direction(
    (sp.Rational(8), sp.Rational(28)),
    (sp.Rational(11, 4), sp.Rational(7)),
) == (4, -1)
incomplete = frontier_72_108_incomplete_certificate()
assert tuple(record.order for record in frontier_records) == (2, 3, 4)
assert len(incomplete.transformations) == 5
assert len(incomplete.missing_data) == 1
frontier_report = frontier_72_108_local_report()
assert tuple(
    item["order"] for item in frontier_report["translations"]
) == (2, 3, 4)
assert tuple(
    item["base_ideal"]["graph_length"]
    for item in frontier_report["translations"]
) == (4, 6, 8)
assert frontier_report["boundary_architecture_missing_data"] == []
assert len(frontier_report["generic_ir_limitations"]) == 1
assert (
    frontier_report["residue_case_tree"]
    == frontier_72_108_residue_case_tree()
)
assert (
    frontier_report["residue_case_tree"][
        "pred_1_minus_3_two_factors"
    ]["constraint"]
    == "alpha_1 != alpha_2"
)
assert frontier_report["unselected_factor_continuation"] == {
    "case_presence": {"a": False, "b": False, "c": True},
    "residue_difference": "beta=alpha_2-alpha_1 != 0",
    "after_selected_order_three": "x^3*y-beta",
    "after_common_order_four_numerator": (
        "x^4*y+alpha-beta*x"
    ),
    "order_four_center_value": "alpha != 0",
    "after_hirzebruch_transition": "X*(Y+alpha)-beta",
    "creates_additional_basepoint": False,
    "meets_filled_x_zero": False,
}
assert frontier_report["exceptional_target_poles"] == {
    "divisors": [f"E{index}" for index in range(1, 9)],
    "P": [1, 1, 1, 8, 6, 4, 2, 1],
    "Q": [0, 0, 1, 12, 9, 6, 3, 2],
    "target_infinity": [1, 1, 1, 12, 9, 6, 3, 2],
    "first_fan_R_orders": [17, 10, 3, -4],
    "residual_R_orders": [-3, -2, -1, 0],
    "complete_global_pullback": False,
}
assert frontier_report["common_graph_keller_pole_audit"]["pole_vector"] == [
    1,
    24,
    1,
    1,
    1,
    12,
    9,
    6,
    3,
    2,
]
assert (
    frontier_report["common_graph_keller_pole_audit"]["geometric_degree"]
    == 427
)
assert not frontier_report["common_graph_keller_pole_audit"][
    "dicritical_candidates"
]
assert frontier_report["forced_first_block_cluster"][
    "crossing"
] == ["E3", "E4"]
assert frontier_report["forced_first_block_cluster"]["root_count"] == 10
assert frontier_report["forced_first_block_cluster"][
    "smooth_E3_basepoints"
] == 0
assert frontier_report["forced_first_block_cluster"]["pole_audit"][
    "geometric_degree"
] == 317
assert frontier_report["terminal_case2_boundary_package"][
    "dicritical_degree"
] == 12
assert frontier_report["terminal_case2_boundary_package"][
    "dicritical_ramification_intersection"
] == 35
assert frontier_report["terminal_case2_boundary_package"]["pole_audit"][
    "geometric_degree"
] == 29
assert frontier_report["terminal_case2_boundary_package"]["pole_audit"][
    "dicritical_candidates"
] == ["T2"]
assert frontier_report["plane_return_edge"][
    "forced_factorization"
] == ["A=a*r^2", "C=c*r^3"]
assert frontier_report["plane_return_edge"][
    "edge_only_root_partitions"
] == [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
assert frontier_report["terminal_case1_boundary_package"][
    "selected_partition"
] == [4]
assert frontier_report["terminal_case1_boundary_package"][
    "adapted_parameter"
] == "v=u-beta+delta*e*u"
assert frontier_report["terminal_case1_boundary_package"][
    "dicritical_ramification_intersection"
] == 35
assert frontier_report["terminal_case1_boundary_package"]["pole_audit"][
    "geometric_degree"
] == 29
assert frontier_report["terminal_case1_boundary_package"]["pole_audit"][
    "dicritical_candidates"
] == ["R1_1"]
report_ramification = {
    item["terminal_case"]: item
    for item in frontier_report["transformed_poisson_ramification"]
}
assert report_ramification["Case 1"][
    "dicritical_generic_ramification_index"
] == 3
assert report_ramification["Case 2"][
    "dicritical_generic_ramification_index"
] == 5
assert report_ramification["Case 2"][
    "dicritical_boundary_ramification_intersection"
] == 33
report_residue_degrees = {
    item["terminal_case"]: item
    for item in frontier_report["dicritical_residue_degree_splits"]
}
assert report_residue_degrees["Case 1"][
    "possible_normalization_cover_degrees"
] == [1, 2, 4]
assert report_residue_degrees["Case 2"][
    "possible_normalization_cover_degrees"
] == [1]
assert report_residue_degrees["Case 2"][
    "remaining_valuation_degrees"
] == [24]
assert frontier_report["terminal_case2_J1_endpoint_exclusion"][
    "unit_ideal"
]
assert frontier_report["terminal_case2_J1_endpoint_exclusion"][
    "generic_infinity_characteristic"
] == [4, 13]
assert "before J0" in frontier_report[
    "terminal_case2_J1_endpoint_exclusion"
]["conclusion"]
assert frontier_report["terminal_case2_lower_jet_reduction"][
    "J0_conclusion"
] == ["B=K*c", "F=K*g"]
assert frontier_report["terminal_case2_maximal_gcd_exclusion"][
    "excluded_exact_gcd_degrees"
] == [7]
assert frontier_report["terminal_case2_maximal_gcd_exclusion"][
    "surviving_exact_gcd_degrees"
] == [1, 2, 3, 4, 5, 6]
assert frontier_report["terminal_case2_gcd_six_exclusion"][
    "excluded_exact_gcd_degrees"
] == [6]
assert frontier_report["terminal_case2_gcd_six_exclusion"][
    "surviving_precompatibility_gcd_degrees"
] == [1, 2, 3, 4, 5]
return_boundaries = frontier_report["plane_return_partition_boundaries"]
assert [item["partition"] for item in return_boundaries] == [
    [4],
    [3, 1],
    [2, 2],
    [2, 1, 1],
    [1, 1, 1, 1],
]
assert all(
    item["pole_audit"]["passes"]
    and item["pole_audit"]["geometric_degree"] == 29
    and len(item["boundary_matrix"]) == len(item["boundary_names"])
    and len(item["intersection_matrix"]) == len(item["boundary_names"])
    for item in return_boundaries
)
assert return_boundaries[1]["root_chains"][0][
    "target_infinity_poles"
] == [3, 0, 0, 0]
assert return_boundaries[1]["root_chains"][0][
    "dicritical_divisor"
] == "R1_4"
assert frontier_report["minimal_dicritical_extension"][
    "tested_center_count"
] == 19
assert frontier_report["minimal_dicritical_extension"][
    "admissible_centers"
] == [["E3"]]
assert frontier_report["minimal_dicritical_extension"][
    "geometric_degree"
] == 426
assert tuple(
    item["first_center"]
    for item in frontier_report["two_step_dicritical_witnesses"]
) == ("Yinf", "E4", "E7", "E8")
assert frontier_report["common_translation_graph"]["chart"] == "GmA1_F4"
assert len(
    frontier_report["common_translation_graph"]["blowups"]
) == 8
assert frontier_report["final_monomial_transition"] == {
    "chart": "GmA1_F4",
    "valuation_matrix": [[-1, 0], [4, 1]],
    "boundary_correspondence": [
        ["pre:X0", "post:Xinf"],
        ["pre:Xinf", "post:X0"],
        ["pre:Yinf", "post:Yinf"],
    ],
    "determinant": -1,
    "involutive": True,
    "additional_blowups": 0,
}
assert frontier_report["affine_plane_fill"]["filled_component"] == "Xinf"
assert (
    frontier_report["affine_plane_fill"]["posttransition_filled_component"]
    == "X0"
)
assert frontier_report["affine_plane_fill"][
    "posttransition_boundary_names"
] == ["Xinf", "Yinf", *(f"E{index}" for index in range(1, 9))]
assert frontier_report["affine_plane_fill"]["passes_prefilter"]
assert frontier_report["affine_plane_fill"]["unit_rank"] == 0
assert frontier_report["affine_plane_fill"]["picard_free_rank"] == 0
assert frontier_report["affine_plane_fill"]["picard_torsion"] == []
assert frontier_report["affine_plane_fill"]["smith_diagonal"] == [1] * 10
assert frontier_report["affine_plane_fill"]["intrinsic_a2_boundary"]["passes"]
assert (
    frontier_report["affine_plane_fill"]["intrinsic_a2_boundary"][
        "canonical_square"
    ]
    == "0"
)
transition = hirzebruch_transition_audit(4)
assert transition.involutive and transition.determinant == -1
filled_common = fill_temporary_boundary(
    frontier_graph_compilations[-1],
    component="Xinf",
    expected_unit_rank=0,
)
assert filled_common.passes_prefilter
try:
    compile_log_boundary(incomplete)
except ValueError as error:
    assert "generic NewtonBoundaryCertificate serialization" in str(error)
else:
    raise AssertionError("the case-unspecialized certificate must not compile")

print("PASS: [u^3:v^2] compiles to its regular toroidal proximity chain")
print("PASS: shared scale prefixes merge and produce a complete unimodular boundary")
print("PASS: certified nested resonance charts attach to earlier exceptional primes")
print("PASS: the bounded Farey grid remains a regular fan")
print("PASS: A2 and GmA1 localization invariants remain chart-aware")
print("PASS: Laurent translations compile at two-boundary SNC crossings")
print("PASS: the (72,108) local orders 2,3,4 compile to regular branch fans")
print("PASS: Laurent map base ideals compile to nested chains of lengths 4,6,8")
print("PASS: the unselected q=3 factor is separated from the return corner")
print("PASS: target-infinity poles are exact on all eight exceptionals")
print("PASS: the F4 transition swaps base divisors before the affine fill")
print("PASS: the common graph has no dicritical and cannot be globally exhaustive")
print("PASS: J4 forces the E3-E4 cluster and excludes a smooth E3 basepoint")
print("PASS: terminal Case 2 compiles to a complete intrinsic boundary package")
print("PASS: the alternate q=3 chart selects and resolves Case-1 partition [4]")
print("PASS: terminal graphs contract to corrected finite-model conductor budgets")
print("PASS: div(X^2) corrects the terminal ramification indices to 3 and 5")
print("PASS: J1 compatibility excludes the forced Case-2 degree-twelve endpoint")
print("PASS: compact J0 certificates exclude precompatibility gcd degrees six,seven")
print("PASS: the degree-one gcd row leaves fourteen Pluecker weights")
print("PASS: all five edge-only quartic models compile intrinsically")
print("PASS: E3 is only the numerical one-blowup zero-pole candidate")
print("PASS: four exact two-step witnesses delimit the minimality claim")
