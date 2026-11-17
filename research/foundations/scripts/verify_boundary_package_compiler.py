#!/usr/bin/env python3
"""Exact regression for the stage-one boundary-package compiler."""

from __future__ import annotations

import json
from dataclasses import replace

from boundary_package_compiler import (
    BoundaryOutputExpressionDatum,
    ConductorBranchJetDatum,
    ConductorBranchSensitivityDatum,
    ContactExpression,
    NormalJetInputDatum,
    PackageStatus,
    StageTwoRealizationCertificate,
    a4_cone_branch_obstruction_package,
    a4_three_puncture_package,
    affine_ramification_obstruction_package,
    all_benchmark_packages,
    compile_boundary_package,
    determinant_ledger_obstruction_package,
    elliptic_selected_boundary_package,
    gl3f2_triangle_package,
    nodal_rational_package,
    presented_core_extension_obstruction_package,
    presented_core_invariants,
    retained_root_euler_obstruction_package,
    semigroup_hole_obstruction_package,
    torus_class_group_obstruction_package,
    core_class_order,
    torus_core_invariants,
)


def invariant(compilation, name):
    return dict(compilation.invariants)[name]


def obstruction_codes(compilation):
    return {
        diagnostic.code
        for diagnostic in compilation.diagnostics
        if diagnostic.obstruction
    }


a4 = compile_boundary_package(a4_three_puncture_package())
assert a4.status is PackageStatus.UNKNOWN
assert invariant(a4, "generated_group_order") == 12
assert invariant(a4, "transitive") is True
assert invariant(a4, "branch_product_identity") is True
assert invariant(a4, "selected_ramified_prime_count") == 2
assert invariant(a4, "negative_principal_part_slots") == 3
assert invariant(a4, "riemann_hurwitz") == {
    "lhs": -2,
    "rhs": -2,
    "ramification": 6,
}
a4_curve = invariant(a4, "curve:three_punctured_rational")
assert a4_curve["unit_rank"] == 2
assert a4_curve["arithmetic_genus"] == 0
assert a4_curve["adjunction_lhs"] == -2
assert invariant(a4, "retained_root_euler_gate")["status"] == "not_applicable"
assert invariant(a4, "conductor_jet_truncation")["status"] == "not_declared"
assert invariant(a4, "conductor_jet_sensitivity")["status"] == "not_declared"
assert not obstruction_codes(a4)
print("PASS: A4 double-transposition package passes all stage-one filters")

gl3f2 = compile_boundary_package(gl3f2_triangle_package())
assert gl3f2.status is PackageStatus.UNKNOWN
assert invariant(gl3f2, "generated_group_order") == 168
assert invariant(gl3f2, "transitive") is True
assert invariant(gl3f2, "branch_product_identity") is True
assert invariant(gl3f2, "selected_ramified_prime_count") == 2
assert invariant(gl3f2, "riemann_hurwitz") == {
    "lhs": -2,
    "rhs": -2,
    "ramification": 12,
}
assert not obstruction_codes(gl3f2)
print("PASS: GL3(F2) (2,3,7) Fano package passes all stage-one filters")

nodal = compile_boundary_package(nodal_rational_package(preserve_conductor=True))
assert nodal.status is PackageStatus.UNKNOWN
nodal_curve = invariant(nodal, "curve:nodal_rational_gm")
assert nodal_curve["normalization_genus"] == 0
assert nodal_curve["arithmetic_genus"] == 1
assert nodal_curve["unit_rank"] == 1
assert nodal_curve["unit_congruence_index"] == 2
assert nodal_curve["adjunction_lhs"] == 0
assert not obstruction_codes(nodal)
print("PASS: nodal rational conductor and finite-index unit lattice compile")

elliptic = compile_boundary_package(elliptic_selected_boundary_package())
assert elliptic.status is PackageStatus.UNKNOWN
assert invariant(elliptic, "generated_group_order") == 6
assert invariant(elliptic, "riemann_hurwitz") == {
    "lhs": 0,
    "rhs": 0,
    "ramification": 6,
}
elliptic_curve = invariant(elliptic, "curve:elliptic_one_puncture")
assert elliptic_curve["normalization_genus"] == 1
assert elliptic_curve["unit_rank"] == 0
assert elliptic_curve["adjunction_lhs"] == 0
assert not obstruction_codes(elliptic)
print("PASS: elliptic selected-boundary package passes exact genus and unit gates")

bad_conductor = compile_boundary_package(
    nodal_rational_package(preserve_conductor=False)
)
assert bad_conductor.status is PackageStatus.OBSTRUCTED
assert "curve.conductor_not_preserved" in obstruction_codes(bad_conductor)
print("PASS: non-equivariant node pairing gives an intrinsic obstruction")

bad_affine = compile_boundary_package(affine_ramification_obstruction_package())
assert bad_affine.status is PackageStatus.OBSTRUCTED
assert "keller.affine_ramification" in obstruction_codes(bad_affine)
print("PASS: ramification inside the Keller source is rejected")

cone_branch = compile_boundary_package(a4_cone_branch_obstruction_package())
assert cone_branch.status is PackageStatus.OBSTRUCTED
assert invariant(cone_branch, "selected_ramified_prime_count") == 1
assert invariant(cone_branch, "riemann_hurwitz") == {
    "lhs": -2,
    "rhs": -2,
    "ramification": 6,
}
assert obstruction_codes(cone_branch) == {"keller.affine_ramification"}
print("PASS: the exact A4 cone (e,f)=(2,2) profile fails only affine coloring")

bad_ledger = compile_boundary_package(determinant_ledger_obstruction_package())
assert bad_ledger.status is PackageStatus.OBSTRUCTED
assert "ledger.determinant_balance" in obstruction_codes(bad_ledger)
print("PASS: an unbalanced target divisor ledger is rejected")

for degree in (2, 3, 5, 6, 7):
    profile = torus_core_invariants(
        ((1, 0), (degree - 2, degree - 1))
    )
    assert profile["unit_rank"] == 0
    assert profile["class_group_free_rank"] == 0
    expected_torsion = () if degree == 2 else (degree - 1,)
    assert profile["class_group_torsion"] == expected_torsion
    assert core_class_order(
        ((1, 0), (degree - 2, degree - 1)), (0, 1)
    ) == degree - 1
    assert core_class_order(
        ((1, 0), (degree - 2, degree - 1)), (1, degree - 2)
    ) == 1
assert core_class_order(((1,), (0,)), (0, 1)) is None
balanced_torus = compile_boundary_package(
    torus_class_group_obstruction_package(5)
)
assert balanced_torus.status is PackageStatus.OBSTRUCTED
assert invariant(balanced_torus, "factorial_core") == {
    "boundary_count": 3,
    "character_rank": 3,
    "matrix_rank": 3,
    "unit_rank": 0,
    "class_group_free_rank": 0,
    "class_group_torsion": (4,),
    "smith_diagonal": (1, 1, 4),
}
assert invariant(balanced_torus, "core_class_orders") == {
    "L1": 4,
    "div(x)": 1,
}
assert "affine_source.core_class_group" in obstruction_codes(balanced_torus)
assert "affine_source.core_required_class" in obstruction_codes(balanced_torus)
print("PASS: the factorial-core Smith gate detects group and class torsion")

split_core = presented_core_invariants(((2,),), ((2,),), ((0,),))
nonsplit_core = presented_core_invariants(((2,),), ((2,),), ((1,),))
assert split_core["smith_diagonal"] == (2, 2)
assert split_core["class_group_torsion"] == (2, 2)
assert nonsplit_core["smith_diagonal"] == (1, 4)
assert nonsplit_core["class_group_torsion"] == (4,)
free_core = presented_core_invariants(((),), ((),), ((),))
assert free_core["unit_rank"] == 0
assert free_core["class_group_free_rank"] == 2
assert free_core["smith_diagonal"] == ()
assert core_class_order(free_core["presentation_matrix"], (0, 0)) == 1
assert core_class_order(free_core["presentation_matrix"], (1, 0)) is None
constant_unit_torsion_core = presented_core_invariants(
    ((),), ((3,),), ((1,),)
)
assert constant_unit_torsion_core["unit_rank"] == 0
assert constant_unit_torsion_core["class_group_torsion"] == ()
assert constant_unit_torsion_core["class_group_free_rank"] == 1
presented = compile_boundary_package(
    presented_core_extension_obstruction_package()
)
assert presented.status is PackageStatus.OBSTRUCTED
assert invariant(presented, "presented_core")["presentation_matrix"] == (
    (2, 1),
    (0, 2),
)
assert invariant(presented, "presented_core_class_orders") == {"G": 4}
assert "affine_source.presented_core_class_group" in obstruction_codes(presented)
assert "affine_source.presented_core_required_class" in obstruction_codes(presented)
print("PASS: lifted core relations distinguish split Z/2+Z/2 from nonsplit Z/4")

bad_semigroup = compile_boundary_package(semigroup_hole_obstruction_package())
assert bad_semigroup.status is PackageStatus.OBSTRUCTED
assert "semigroup.membership" in obstruction_codes(bad_semigroup)
print("PASS: an actual-semigroup hole is not silently filled by saturation")

retained_quartic = compile_boundary_package(
    retained_root_euler_obstruction_package(4)
)
assert retained_quartic.status is PackageStatus.OBSTRUCTED
retained_quartic_gate = invariant(
    retained_quartic, "retained_root_euler_gate"
)
assert retained_quartic_gate["status"] == "obstructed"
assert retained_quartic_gate["geometric_euler_characteristic"] == 4
assert retained_quartic_gate["geometric_open_class"] == "L^2+3*L"
assert retained_quartic_gate["finite_field_open_count"] == (
    "q^2+(n_q(A)-1)*q"
)
assert "affine_source.retained_root_euler" in obstruction_codes(
    retained_quartic
)

retained_linear = compile_boundary_package(
    retained_root_euler_obstruction_package(1)
)
assert retained_linear.status is PackageStatus.UNKNOWN
assert invariant(retained_linear, "retained_root_euler_gate")["status"] == (
    "passes"
)

uncertified_datum = replace(
    retained_root_euler_obstruction_package(4).retained_root_euler,
    different_support_certificate="",
)
uncertified_retained = compile_boundary_package(
    replace(
        a4_three_puncture_package(),
        name="uncertified retained-root input",
        retained_root_euler=uncertified_datum,
    )
)
assert uncertified_retained.status is PackageStatus.UNKNOWN
assert invariant(
    uncertified_retained, "retained_root_euler_gate"
)["status"] == "uncertified"
assert not obstruction_codes(uncertified_retained)
print("PASS: retained-root Euler gate rejects only certified nonlinear rows")


def conductor_branch(
    available_jet_order,
    *,
    certificates=True,
):
    certificate = "exact regression certificate" if certificates else ""
    return ConductorBranchJetDatum(
        name="cusp-branch",
        conductor_exponent=2,
        differential_order=1,
        pole_order=2,
        additional_contact_loss=1,
        available_jet_order=available_jet_order,
        conductor_certificate=certificate,
        expression_certificate=certificate,
        valuation_certificate=certificate,
    )


jet_pass = compile_boundary_package(
    replace(
        a4_three_puncture_package(),
        name="conductor jet pass",
        conductor_jet_branches=(conductor_branch(6),),
    )
)
assert invariant(jet_pass, "conductor_jet_truncation")["status"] == "passes"
assert invariant(jet_pass, "conductor_jet_truncation")["branches"][0][
    "required_jet_order"
] == 6
assert jet_pass.status is PackageStatus.UNKNOWN

jet_short = compile_boundary_package(
    replace(
        a4_three_puncture_package(),
        name="conductor jet short",
        conductor_jet_branches=(conductor_branch(5),),
    )
)
assert invariant(jet_short, "conductor_jet_truncation")["status"] == (
    "insufficient"
)
assert jet_short.status is PackageStatus.UNKNOWN
assert "boundary_module.conductor_jet_truncation" not in obstruction_codes(
    jet_short
)

jet_uncertified = compile_boundary_package(
    replace(
        a4_three_puncture_package(),
        name="conductor jet uncertified",
        conductor_jet_branches=(
            conductor_branch(None, certificates=False),
        ),
    )
)
assert invariant(jet_uncertified, "conductor_jet_truncation")["status"] == (
    "uncertified"
)
assert jet_uncertified.status is PackageStatus.UNKNOWN
print("PASS: conductor/contact-loss truncation is proof-bearing and non-obstructing")


detailed_expression = ContactExpression.combine(
    "add",
    ContactExpression.shift("derivative", ContactExpression.input("P"), 1),
    ContactExpression.shift(
        "pole",
        ContactExpression.shift(
            "derivative", ContactExpression.input("Q"), 2
        ),
        2,
    ),
)
detailed_branch = ConductorBranchSensitivityDatum(
    name="asymmetric-cusp",
    conductor_exponent=2,
    inputs=(
        NormalJetInputDatum("P", 3, "P valuation"),
        NormalJetInputDatum("Q", 6, "Q valuation"),
    ),
    outputs=(
        BoundaryOutputExpressionDatum(
            "rho", detailed_expression, "complete expression tree"
        ),
    ),
    conductor_certificate="cusp conductor",
    dependency_completeness_certificate="all P/Q paths listed",
)
detailed_pass = compile_boundary_package(
    replace(
        a4_three_puncture_package(),
        name="detailed conductor pass",
        conductor_jet_sensitivity=(detailed_branch,),
    )
)
detailed_invariant = invariant(detailed_pass, "conductor_jet_sensitivity")
assert detailed_invariant["status"] == "passes"
assert detailed_invariant["maximum_deficit"] == 0
assert [
    (row["input_name"], row["required_jet_order"])
    for row in detailed_invariant["requirements"]
] == [("P", 3), ("Q", 6)]
assert detailed_pass.status is PackageStatus.UNKNOWN

detailed_short = compile_boundary_package(
    replace(
        a4_three_puncture_package(),
        name="detailed conductor short",
        conductor_jet_sensitivity=(
            replace(
                detailed_branch,
                inputs=(
                    detailed_branch.inputs[0],
                    NormalJetInputDatum("Q", 5, "short Q valuation"),
                ),
            ),
        ),
    )
)
assert invariant(detailed_short, "conductor_jet_sensitivity")["status"] == (
    "insufficient"
)
assert invariant(detailed_short, "conductor_jet_sensitivity")[
    "maximum_deficit"
] == 1
assert detailed_short.status is PackageStatus.UNKNOWN

try:
    compile_boundary_package(
        replace(
            a4_three_puncture_package(),
            conductor_jet_branches=(conductor_branch(6),),
            conductor_jet_sensitivity=(detailed_branch,),
        )
    )
except ValueError:
    pass
else:
    raise AssertionError("conflicting scalar and detailed ledgers were accepted")
print("PASS: detailed conductor sensitivities preserve asymmetric P/Q jet budgets")

unreplayed_stage_two = replace(
    a4_three_puncture_package(),
    stage_two_realization=StageTwoRealizationCertificate(
        root_equation="candidate/root-equation.cert",
        local_factorization_certificate="candidate/local-factorization.cert",
        reconstruction_identities="candidate/reconstruction.cert",
        polynomial_ring_isomorphism="candidate/affine-source.cert",
        constant_jacobian_certificate="candidate/jacobian.cert",
        monodromy_certificate="candidate/monodromy.cert",
    ),
)
unreplayed = compile_boundary_package(unreplayed_stage_two)
assert unreplayed.status is PackageStatus.UNKNOWN
assert any(
    item.code == "stage_two.unverified" for item in unreplayed.diagnostics
)
print("PASS: unreplayed stage-two references cannot claim realization")

reports = [
    compile_boundary_package(package).to_dict()
    for package in all_benchmark_packages()
]
assert [report["status"] for report in reports] == [
    "unknown",
    "unknown",
    "unknown",
    "unknown",
    "obstructed",
    "obstructed",
    "obstructed",
    "obstructed",
    "obstructed",
    "obstructed",
    "obstructed",
    "obstructed",
]
print(json.dumps(reports, indent=2, sort_keys=True))
print("PASS: stage-one compiler preserves the realized/obstructed/unknown boundary")
