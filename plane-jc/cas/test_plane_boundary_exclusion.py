#!/usr/bin/env python3
"""Regression tests for the plane residue/conductor obstruction."""

from math import gcd

import sympy as sp

from plane_boundary_exclusion import (
    AffinePrimeLedgerEntry,
    BoundaryPrimeLedgerEntry,
    ConductorPacketPoint,
    OneDicriticalNormalizationCertificate,
    TargetComponentLedger,
    TargetFiberPacket,
    audit_one_dicritical_normalization,
    audit_target_component_ledger,
    conductor_collision_budget,
    conductor_packet_budget,
    first_free_depth_package,
    one_puncture_budget,
    puncture_profile_budgets,
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

# The rational-boundary theorem makes the arbitrary-puncture profile
# available.  For every ordered positive pole partition, Riemann--Hurwitz
# forces f+s-2 affine ramification.  An immersive residue map can therefore
# occur only for the one-puncture degree-one profile.
for degree in range(1, 9):
    for punctures in range(1, degree + 1):
        budgets = puncture_profile_budgets(degree, punctures)
        assert len(budgets) > 0
        for budget in budgets:
            assert (
                budget.forced_affine_ramification
                == degree + punctures - 2
            )
            assert (
                budget.forced_affine_ramification == 0
            ) == (degree == 1 and punctures == 1)

assert puncture_profile_budgets(3, 4) == ()

three_puncture_gate = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="three-puncture rational boundary",
        generic_degree=7,
        transverse_index=2,
        residue_degree=3,
        affine_degree=1,
        punctures=3,
        single_normalization_boundary=True,
        log_pure=True,
        exhaustive_pullback=True,
        target_transfer_certified=True,
    )
)
assert three_puncture_gate.status == "excluded"
assert three_puncture_gate.hurwitz_forced_affine_ramification == 4

impossible_profile_gate = audit_one_dicritical_normalization(
    OneDicriticalNormalizationCertificate(
        name="too many punctures for residue degree",
        generic_degree=4,
        transverse_index=1,
        residue_degree=2,
        affine_degree=2,
        punctures=3,
        single_normalization_boundary=True,
        log_pure=True,
        exhaustive_pullback=True,
        target_transfer_certified=True,
    )
)
assert impossible_profile_gate.status == "excluded"
assert "no positive pole profile" in impossible_profile_gate.reasons[0]


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

# A packet is local to one closed target fiber.  Its transverse lengths add,
# while its comparison degree is the full generic boundary contribution plus
# the affine contribution over the same target component.
mixed_packet = conductor_packet_budget(
    transverse_indices=(2, 3, 4),
    generic_boundary_degree=5,
    affine_degree=3,
)
assert mixed_packet.minimum_packet_length == 9
assert mixed_packet.generic_degree == 8
assert mixed_packet.minimum_affine_degree == 4
assert mixed_packet.length_deficit == 1
assert "excluded" in mixed_packet.verdict

paid_mixed_packet = conductor_packet_budget(
    transverse_indices=(2, 3, 4),
    generic_boundary_degree=5,
    affine_degree=4,
)
assert paid_mixed_packet.generic_degree == 9
assert paid_mixed_packet.length_deficit == 0

# For one boundary prime of residue degree one and constant index e, an
# r-point packet requires a >= (r-1)e.
for transverse_index in range(1, 6):
    for packet_size in range(2, 6):
        required_affine_degree = (packet_size - 1) * transverse_index
        packet = conductor_packet_budget(
            transverse_indices=(transverse_index,) * packet_size,
            generic_boundary_degree=transverse_index,
            affine_degree=required_affine_degree,
        )
        assert packet.minimum_affine_degree == required_affine_degree
        assert packet.length_deficit == 0

# The target-transfer ledger keeps generic divisor degrees separate from
# special-fiber packets.  This is the primitive d=3, e=2, a=1 exclusion.
primitive_ledger = TargetComponentLedger(
    name="primitive two-point conductor packet",
    generic_degree=3,
    boundary_primes=(
        BoundaryPrimeLedgerEntry("D", transverse_index=2, residue_degree=1),
    ),
    affine_primes=(AffinePrimeLedgerEntry("E", residue_degree=1),),
    packets=(
        TargetFiberPacket(
            name="node fiber",
            points=(
                ConductorPacketPoint("p", "D", 2, residue_immersive=True),
                ConductorPacketPoint("q", "D", 2, residue_immersive=True),
            ),
            same_target_fiber_certified=True,
        ),
    ),
    finite_flat_certified=True,
    target_transfer_certified=True,
    exhaustive_generic_pullback=True,
)
primitive_ledger_audit = audit_target_component_ledger(primitive_ledger)
assert primitive_ledger_audit.generic_boundary_degree == 2
assert primitive_ledger_audit.affine_degree == 1
assert primitive_ledger_audit.degree_identity_holds
assert primitive_ledger_audit.status == "excluded"
assert primitive_ledger_audit.packet_audits[0].length_deficit == 1

# Multiple normalization-boundary primes and affine sheets can pay exactly
# for a mixed packet.
mixed_ledger = TargetComponentLedger(
    name="paid mixed packet",
    generic_degree=9,
    boundary_primes=(
        BoundaryPrimeLedgerEntry("D2", transverse_index=2, residue_degree=1),
        BoundaryPrimeLedgerEntry("D3", transverse_index=3, residue_degree=1),
    ),
    affine_primes=(
        AffinePrimeLedgerEntry("A1", residue_degree=1),
        AffinePrimeLedgerEntry("A2", residue_degree=3),
    ),
    packets=(
        TargetFiberPacket(
            name="three-point fiber",
            points=(
                ConductorPacketPoint("p2", "D2", 2, residue_immersive=True),
                ConductorPacketPoint("p3a", "D3", 3, residue_immersive=True),
                ConductorPacketPoint("p3b", "D3", 3, residue_immersive=True),
            ),
            same_target_fiber_certified=True,
        ),
    ),
    finite_flat_certified=True,
    target_transfer_certified=True,
    exhaustive_generic_pullback=True,
)
mixed_ledger_audit = audit_target_component_ledger(mixed_ledger)
assert mixed_ledger_audit.generic_boundary_degree == 5
assert mixed_ledger_audit.affine_degree == 4
assert mixed_ledger_audit.status == "survives_packet_budget"
assert mixed_ledger_audit.packet_audits[0].minimum_packet_length == 8
assert mixed_ledger_audit.packet_audits[0].length_deficit == -1

# A numerical packet is deliberately refused when target grouping or
# pointwise residue immersion has not been certified.
untransferred_ledger = TargetComponentLedger(
    name="source-only packet preview",
    generic_degree=3,
    boundary_primes=(
        BoundaryPrimeLedgerEntry("D", transverse_index=2, residue_degree=1),
    ),
    affine_primes=(AffinePrimeLedgerEntry("E", residue_degree=1),),
    packets=(
        TargetFiberPacket(
            name="uncertified collision",
            points=(
                ConductorPacketPoint("p", "D", 2, residue_immersive=True),
                ConductorPacketPoint("q", "D", 2, residue_immersive=False),
            ),
            same_target_fiber_certified=False,
        ),
    ),
    finite_flat_certified=True,
    target_transfer_certified=False,
    exhaustive_generic_pullback=True,
)
untransferred_ledger_audit = audit_target_component_ledger(untransferred_ledger)
assert untransferred_ledger_audit.status == "incomplete"
assert not untransferred_ledger_audit.packet_audits[0].applicable
assert "one target fiber" in untransferred_ledger_audit.packet_audits[0].reasons[0]

# Large curve conductor does not force a large point packet.  These immersed
# one-place parametrizations have two points over the origin and arbitrary
# branch contact m.
t, c = sp.symbols("t c")
for contact_order in range(2, 8):
    x = t**2 - 1
    y = c * x + t * x**contact_order
    implicit_pullback = sp.expand(
        (y - c * x) ** 2 - x ** (2 * contact_order) * (x + 1)
    )
    assert implicit_pullback == 0
    assert (x.subs(t, 1), y.subs(t, 1)) == (0, 0)
    assert (x.subs(t, -1), y.subs(t, -1)) == (0, 0)
    assert sp.diff(y, t).subs(t, 0) == (-1) ** contact_order
    assert gcd(2, 2 * contact_order + 1) == 1


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
print("PASS: arbitrary puncture profiles force affine ramification f+s-2")
print("PASS: minimal-sheet fiber length excludes conductor gluing")
print("PASS: mixed and constant-index conductor packets obey the additive bound")
print("PASS: target-transfer ledgers separate generic and special-fiber data")
print("PASS: arbitrary conductor tangency can remain a two-point packet")
print("PASS: the first numerical degree-six package is residue-excluded")
print("PASS: typed normalization gate refuses source-only (72,108) data")
print("PASS: the Case-1 sheet budget permits conductor gluing")
