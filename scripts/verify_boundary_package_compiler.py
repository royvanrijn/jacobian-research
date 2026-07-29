#!/usr/bin/env python3
"""Exact regression for the stage-one boundary-package compiler."""

from __future__ import annotations

import json
from dataclasses import replace

from boundary_package_compiler import (
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
    semigroup_hole_obstruction_package,
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

bad_semigroup = compile_boundary_package(semigroup_hole_obstruction_package())
assert bad_semigroup.status is PackageStatus.OBSTRUCTED
assert "semigroup.membership" in obstruction_codes(bad_semigroup)
print("PASS: an actual-semigroup hole is not silently filled by saturation")

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
]
print(json.dumps(reports, indent=2, sort_keys=True))
print("PASS: stage-one compiler preserves the realized/obstructed/unknown boundary")
