#!/usr/bin/env python3
"""Regression tests for band-limited Poisson-square rigidity."""

import sympy as sp

from poisson_square_rigidity import (
    c_zero_curvilinear_family,
    c_zero_dual_number_family,
    c_zero_local_slice_audit,
    component_tangent_space_audit,
    degenerate_reduced_families,
    intersection_tangent_space_audit,
    intersection_tangent_kernel_slice_audit,
    lower_layer_degree_audit,
    normalized_residual_system,
    primary_chart_scaling_audit,
    reduced_intersection_strata,
    reduced_three_band_family,
    standard_three_band_problem,
    three_band_layers,
    verify_reduced_classification,
)


problem = standard_three_band_problem()
assert tuple(problem.layers) == (4, 3, 2)
assert len(problem.P_bands) == 2
assert len(problem.Q_bands) == 2
assert len(problem.coefficient_equations) == 17
assert primary_chart_scaling_audit() == {
    "coefficient_weights": (-2, -1, 0, 1),
    "layer_weights": (-2, -1, 0),
    "normalizing_parameter": "lambda=d0^-1",
    "normalized_constant": 1,
}

verify_reduced_classification()
family = reduced_three_band_family()
assert sp.expand(family.bracket - family.X**2) == 0
assert tuple(
    sp.expand(layer)
    for layer in three_band_layers(
        family.A, family.B, family.C, family.D, family.t
    )
) == (0, 0, 1)

# Recover the foundational weighted-tangent bands after h=t-1:
# a=1/2, d=2, e=-4.
specialization = {
    family.root: 1,
    family.cubic_scale: sp.Rational(1, 2),
    family.gamma_scale: 2,
    family.tangent_linear: -4,
}
assert sp.factor(family.A.subs(specialization)) == (family.t - 1) ** 3 / 2
assert sp.factor(family.B.subs(specialization)) == (
    (family.t - 2) * (family.t - 1) / 2
)
assert sp.factor(family.C.subs(specialization)) == -3 * (family.t - 1) ** 2
assert sp.expand(
    family.D.subs(specialization) + 2 * (2 * family.t - 3)
) == 0

residual = normalized_residual_system()
assert tuple(degree for degree, _ in residual.equations_by_degree) == (
    6,
    5,
    4,
    3,
    1,
    0,
)
assert "d2^2" in residual.nilpotent_warning
degenerate = degenerate_reduced_families()
assert tuple(family.name for family in degenerate) == (
    "C_zero_D_constant",
    "C_zero_D_linear",
    "A_zero_D_constant",
    "A_zero_D_linear",
)
assert all(
    sp.expand(family.bracket - family.X**2) == 0 for family in degenerate
)
degree_audit = lower_layer_degree_audit()
assert len(degree_audit["degree_2"]) == 5
assert len(degree_audit["degree_3"]) == 7
c_zero_dual, epsilon = c_zero_dual_number_family()
assert sp.rem(
    sp.Poly(sp.expand(c_zero_dual.bracket - c_zero_dual.X**2), epsilon),
    sp.Poly(epsilon**2, epsilon),
).as_expr() == 0
c_zero_length_three, epsilon = c_zero_curvilinear_family()
length_three_residual = sp.Poly(
    sp.expand(
        c_zero_length_three.bracket
        - c_zero_length_three.X**2
    ),
    epsilon,
)
assert sp.rem(
    length_three_residual,
    sp.Poly(epsilon**3, epsilon),
).as_expr() == 0
assert sp.rem(
    length_three_residual,
    sp.Poly(epsilon**4, epsilon),
).as_expr() != 0
local_slice = c_zero_local_slice_audit()
assert sp.factor(
    local_slice["constant_after_middle"]
    - local_slice["expected_constant"]
) == 0
assert tuple(
    (report.component, report.jacobian_rank, report.tangent_dimension)
    for report in component_tangent_space_audit()
) == (
    ("tangent_closure", 11, 5),
    ("C_zero", 11, 5),
    ("A_zero", 12, 4),
)
intersection_strata = reduced_intersection_strata()
assert tuple(stratum.name for stratum in intersection_strata) == (
    "triple_lower_wronskian_core",
    "tangent_C_zero_boundary",
    "tangent_A_zero_boundary",
)
assert tuple(
    stratum.component_membership for stratum in intersection_strata
) == (
    ("tangent", "C_zero", "A_zero"),
    ("tangent", "C_zero"),
    ("tangent", "A_zero"),
)
assert all(
    sp.expand(stratum.family.bracket - stratum.family.X**2) == 0
    for stratum in intersection_strata
)
assert tuple(
    (report.component, report.jacobian_rank, report.tangent_dimension)
    for report in intersection_tangent_space_audit()
) == (
    ("triple_lower_wronskian_core", 8, 8),
    ("tangent_C_zero_boundary", 9, 7),
    ("tangent_A_zero_boundary", 10, 6),
)
kernel_slices = intersection_tangent_kernel_slice_audit()
assert tuple(
    (
        report.stratum,
        report.tangent_kernel_dimension,
        report.length,
        report.hilbert_vector,
        report.socle_dimension,
    )
    for report in kernel_slices
) == (
    ("tangent_C_zero_boundary", 4, 5, (1, 4), 4),
    ("tangent_A_zero_boundary", 3, 5, (1, 3, 1), 1),
)

print("PASS: generic bounded supports compile to exactly three Wronskian layers")
print("PASS: the reduced nondegenerate locus is the tangent-pencil family")
print("PASS: four explicit charts classify every degenerate reduced point")
print("PASS: the foundational weighted pair is its normalized cubic member")
print("PASS: generic multiplicities are 2 on tangent, 3 on C=0, 1 on A=0")
print("PASS: dense charts cover the three reduced intersection branches")
print("PASS: intersection tangent dimensions are 8,7,6 on the dense charts")
print("PASS: two extra-branch kernel slices have distinct length-five algebras")
print("PASS: the d0 principal chart has an exact Gm normalization")
