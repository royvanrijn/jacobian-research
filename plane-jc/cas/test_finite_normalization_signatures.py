#!/usr/bin/env python3
"""Regression tests for the bounded finite-normalization signature atlas."""

from finite_normalization_signatures import (
    enumerate_target_signatures,
    pareto_minimal_signatures,
)


assert enumerate_target_signatures(1) == ()

expected_counts = {
    2: (1, 0, 0),
    3: (6, 1, 1),
    4: (19, 4, 4),
    5: (54, 15, 11),
    6: (132, 42, 25),
    7: (312, 115, 51),
    8: (691, 282, 97),
}

for degree in range(2, 9):
    signatures = enumerate_target_signatures(degree)
    ramified = enumerate_target_signatures(
        degree,
        require_transverse_ramification=True,
    )
    ramified_immersive = enumerate_target_signatures(
        degree,
        require_transverse_ramification=True,
        require_residue_immersion=True,
    )
    assert (
        len(signatures),
        len(ramified),
        len(ramified_immersive),
    ) == expected_counts[degree]
    assert signatures
    for signature in signatures:
        assert signature.total_degree == degree
        assert signature.affine_degree >= 1
        assert signature.boundary_degree <= degree - 1
        for row in signature.boundary_rows:
            assert 1 <= row.punctures <= row.residue_degree
            assert (
                row.forced_affine_residue_ramification
                == row.residue_degree + row.punctures - 2
            )
            assert (
                row.puncture_ramification
                + row.forced_affine_residue_ramification
                == row.total_residue_ramification
            )

# A connected counterexample needs some transverse branch row.  The affine
# companion sheet then excludes geometric degree two.
assert enumerate_target_signatures(
    2,
    require_transverse_ramification=True,
) == ()

# In geometric degree three, residue immersion and transverse ramification
# force exactly the boundary row (e,f,s)=(2,1,1) and one affine degree-one
# row.
degree_three = enumerate_target_signatures(
    3,
    require_transverse_ramification=True,
    require_residue_immersion=True,
)
assert len(degree_three) == 1
degree_three_signature = degree_three[0]
assert tuple(
    (
        row.ramification_index,
        row.residue_degree,
        row.punctures,
    )
    for row in degree_three_signature.boundary_rows
) == ((2, 1, 1),)
assert degree_three_signature.affine_residue_degrees == (1,)

# Residue immersion uniformly collapses every boundary row to f=s=1.
for degree in range(3, 9):
    immersive = enumerate_target_signatures(
        degree,
        require_residue_immersion=True,
    )
    assert immersive
    for signature in immersive:
        assert signature.residual_ramification_cost == 0
        assert all(
            row.residue_degree == row.punctures == 1
            for row in signature.boundary_rows
        )

# Pareto extraction is intrinsic to the chosen complexity coordinates and
# leaves no signature dominated by another retained signature.
for degree in range(2, 9):
    signatures = enumerate_target_signatures(degree)
    pareto = pareto_minimal_signatures(signatures)
    assert pareto
    complexities = [signature.pareto_complexity[1:] for signature in pareto]
    for index, left in enumerate(complexities):
        for right in complexities[index + 1 :]:
            assert not (
                all(x <= y for x, y in zip(left, right))
                and any(x < y for x, y in zip(left, right))
            )
            assert not (
                all(y <= x for x, y in zip(left, right))
                and any(y < x for x, y in zip(left, right))
            )

print("PASS: finite normalization signatures exhaust degrees two through eight")
print("PASS: affine-sheet positivity removes every ramified degree-two row")
print("PASS: the immersive ramified degree-three atlas is the forced 2+1 ledger")
print("PASS: residue immersion collapses every boundary row to f=s=1")
print("PASS: Pareto extraction leaves an antichain of normalization signatures")
