#!/usr/bin/env python3
"""Regression tests for the plane residue/conductor obstruction."""

from plane_boundary_exclusion import (
    OneDicriticalNormalizationCertificate,
    audit_one_dicritical_normalization,
    conductor_collision_budget,
    first_free_depth_package,
    one_puncture_budget,
    two_puncture_budgets,
)


# An immersive finite A1 -> A1 map has no affine ramification.  The exact
# budget therefore forces degree one.
for degree in range(1, 9):
    budget = one_puncture_budget(degree)
    assert budget.forced_affine_ramification == degree - 1
    assert (budget.forced_affine_ramification == 0) == (degree == 1)


# Every finite Gm -> A1 map has positive ramification on Gm, independently
# of how its degree is split between the two poles.
for degree in range(2, 9):
    budgets = two_puncture_budgets(degree)
    assert len(budgets) == degree - 1
    for budget in budgets:
        assert budget.forced_affine_ramification == degree
        assert budget.forced_affine_ramification > 0


# In a primitive minimal link, d=e+1.  A conductor identification costs two
# boundary packets, so it requires 2e<=d.  This fails directly for e>=2;
# e=1 leaves a forbidden connected degree-two finite etale cover of A2.
for transverse_index in range(1, 8):
    budget = conductor_collision_budget(transverse_index)
    assert budget.affine_degree == 1
    assert budget.generic_degree == transverse_index + 1
    assert budget.minimum_collision_length == 2 * transverse_index
    assert budget.length_deficit == transverse_index - 1
    if transverse_index == 1:
        assert "finite etale" in budget.verdict
    else:
        assert budget.minimum_collision_length > budget.generic_degree

# More generally a collision is numerically possible only if the total
# affine degree is at least the transverse boundary index.
for transverse_index in range(2, 8):
    for affine_degree in range(1, 8):
        budget = conductor_collision_budget(
            transverse_index, affine_degree
        )
        assert (budget.length_deficit <= 0) == (
            affine_degree >= transverse_index
        )


# The existing intersection gate accepts its first numerical package in
# degree six.  The new residue gate sees target-line degree one and excludes
# the resulting affine-line Jelonek component.
package = first_free_depth_package()
assert package["passes_intrinsic_gates"]
assert package["canonical_coefficients"] == (-3, -2, -1, 0)
assert package["pole_vector"] == (3, 2, 1, 0)
assert package["hyperplane_intersections"] == (2, 0, 0, 1)
assert package["geometric_degree"] == 6
assert package["dicritical_candidates"] == ("E3",)
assert package["dicritical_target_line_degree"] == 1


minimal_link = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="minimal link",
        generic_degree=3,
        transverse_index=2,
        residue_degree=1,
        affine_degree=1,
        punctures=1,
        single_normalization_boundary=True,
        log_pure=True,
        exhaustive_pullback=True,
        target_transfer_certified=True,
    )
)
assert minimal_link.status == "excluded"
assert minimal_link.sheet_deficient
assert minimal_link.conductor_length_deficit == 1

two_puncture = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="two punctures",
        generic_degree=7,
        transverse_index=2,
        residue_degree=2,
        affine_degree=3,
        punctures=2,
        single_normalization_boundary=True,
        log_pure=True,
        exhaustive_pullback=True,
        target_transfer_certified=True,
    )
)
assert two_puncture.status == "excluded"
assert two_puncture.hurwitz_forced_affine_ramification == 2

nontrivial_one_puncture_cover = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="nontrivial one-puncture residue cover",
        generic_degree=7,
        transverse_index=2,
        residue_degree=2,
        affine_degree=3,
        punctures=1,
        single_normalization_boundary=True,
        log_pure=True,
        exhaustive_pullback=True,
        target_transfer_certified=True,
    )
)
assert nontrivial_one_puncture_cover.status == "excluded"
assert nontrivial_one_puncture_cover.hurwitz_forced_affine_ramification == 1

index_one = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="index-one boundary",
        generic_degree=2,
        transverse_index=1,
        residue_degree=1,
        affine_degree=1,
        punctures=1,
        single_normalization_boundary=True,
        log_pure=True,
        exhaustive_pullback=True,
        target_transfer_certified=True,
    )
)
assert index_one.status == "excluded"
assert "etale" in index_one.reasons[0]

case1_untransferred = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="72,108 Case 1 source package",
        generic_degree=29,
        transverse_index=3,
        residue_degree=1,
        affine_degree=26,
        punctures=1,
        single_normalization_boundary=False,
        log_pure=True,
        exhaustive_pullback=True,
        target_transfer_certified=False,
    )
)
assert case1_untransferred.status == "incomplete"
assert not case1_untransferred.residue_immersion_certified
assert "finite normalization" in case1_untransferred.reasons[0]

# This is a numerical preview, not a claim that the missing target transfer
# has been proved.  Even after supplying it, the large affine remainder pays
# for a two-point conductor collision.
case1_preview = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="72,108 Case 1 transferred preview",
        generic_degree=29,
        transverse_index=3,
        residue_degree=1,
        affine_degree=26,
        punctures=1,
        single_normalization_boundary=True,
        log_pure=True,
        exhaustive_pullback=True,
        target_transfer_certified=True,
    )
)
assert case1_preview.status == "survives_conductor_budget"
assert not case1_preview.sheet_deficient
assert case1_preview.conductor_length_deficit == -23

case2_preview = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="72,108 Case 2 source package",
        generic_degree=29,
        transverse_index=5,
        residue_degree=1,
        affine_degree=24,
        punctures=1,
        single_normalization_boundary=False,
        log_pure=False,
        exhaustive_pullback=True,
        target_transfer_certified=False,
    )
)
assert case2_preview.status == "incomplete"
assert "ramification" in case2_preview.reasons[-1]


print("PASS: one-puncture residue immersion forces degree one")
print("PASS: two-puncture residue immersion contradicts Riemann--Hurwitz")
print("PASS: minimal-sheet fiber length excludes conductor gluing")
print("PASS: the first numerical degree-six package is residue-excluded")
print("PASS: typed normalization gate refuses source-only (72,108) data")
print("PASS: the Case-1 sheet budget permits conductor gluing")
