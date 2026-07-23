#!/usr/bin/env python3
"""Exact contact-atom certificates with optional legacy Mason regressions."""

import sys
from math import lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    affine_difference_mason_defect,
    collision_precedes,
    contact_atom_refinement,
    contact_semigroup_atoms,
    maximal_threshold_partitions,
    multiplicity_excess,
    threshold_component_count,
    threshold_dominant_length,
    threshold_mason_margin,
)


def threshold_partitions(total, threshold, maximum=None):
    """Partitions of total with every part at least the threshold."""
    if total == 0:
        yield ()
        return
    maximum = min(total, maximum or total)
    for first in range(maximum, threshold - 1, -1):
        for tail in threshold_partitions(total - first, threshold, first):
            yield (first,) + tail


# For S_r={m>=r}, the indecomposable elements are exactly r,...,2r-1.
for threshold in range(2, 9):
    atoms = set(contact_semigroup_atoms(threshold))
    assert atoms == set(range(threshold, 2 * threshold))
    for multiplicity in range(threshold, 8 * threshold + 1):
        decomposable = any(
            left >= threshold and multiplicity - left >= threshold
            for left in range(threshold, multiplicity)
        )
        assert (multiplicity in atoms) == (not decomposable)

# Every allowed partition refines to an atomic partition, and collision merges
# the returned atoms back to the original multiplicity type.
for threshold in range(2, 7):
    atoms = set(contact_semigroup_atoms(threshold))
    for degree in range(threshold, 6 * threshold + 1):
        for partition in threshold_partitions(degree, threshold):
            refinement = contact_atom_refinement(partition, threshold)
            assert set(refinement) <= atoms
            assert sum(refinement) == degree
            assert collision_precedes(refinement, partition)

# Atomic partitions are exactly the maximal threshold types.  Their counting
# series is product_{m=r}^{2r-1}(1-t^m)^(-1).
for threshold in range(2, 8):
    for degree in range(threshold, 8 * threshold + 1):
        atomic_partitions = maximal_threshold_partitions(degree, threshold)
        assert len(atomic_partitions) == threshold_component_count(
            degree, threshold
        )
        assert all(
            set(partition) <= set(range(threshold, 2 * threshold))
            for partition in atomic_partitions
        )

# Optional legacy comparison for r=2.  The atoms are {2,3}, and the Mason
# defect for any two degree-n types is exactly half their total excess above
# double contact.  This is a regression for the old distinct-type proof, not
# the uniqueness input; verify_unique_omitted_value.py covers all pairs and
# closes the all-double equality case.
assert contact_semigroup_atoms(2) == (2, 3)
for degree in range(3, 25):
    partitions = tuple(threshold_partitions(degree, 2))
    for left_index, left in enumerate(partitions):
        for right in partitions[left_index + 1 :]:
            defect = affine_difference_mason_defect(left, right)
            excess_sum = multiplicity_excess(left) + multiplicity_excess(right)
            assert 2 * defect == excess_sum
            assert defect >= 1

# For a hypothetical threshold r>=3, abc separation is no longer borderline:
# each radical support has size at most n/r, so two supports total < n.
for threshold in range(3, 8):
    for degree in range(threshold, 8 * threshold + 1):
        partitions = tuple(threshold_partitions(degree, threshold))
        for left in partitions:
            for right in partitions:
                support_sum = len(left) + len(right)
                assert threshold * support_sum <= 2 * degree
                assert support_sum < degree
        assert threshold_mason_margin(degree, threshold, 1) >= 1

# If the difference has degree at most d, the strict automatic range is
# n-2*floor(n/r)-d+1>0.  This records the version needed by higher-dimensional
# polynomial pencils rather than only affine pencils.
for threshold in range(3, 8):
    for difference_degree in range(0, 6):
        for degree in range(threshold, 12 * threshold + 1):
            margin = threshold_mason_margin(
                degree, threshold, difference_degree
            )
            assert margin == (
                degree
                - 2 * (degree // threshold)
                - difference_degree
                + 1
            )

# Maximizing the atomic-stratum dimension a+b-1 over 2a+3b=n gives the
# omitted-value phase diagram.  The normalized seed space has dimension n-3.
for degree in range(4, 101):
    atomic_types = [
        (a, b)
        for a in range(degree // 2 + 1)
        for b in range(degree // 3 + 1)
        if 2 * a + 3 * b == degree
    ]
    maximal_dimension = max(a + b - 1 for a, b in atomic_types)
    expected_dimension = degree // 2 - 1
    expected_codimension = (degree - 3) // 2
    assert maximal_dimension == expected_dimension
    assert degree - 3 - maximal_dimension == expected_codimension

assert (4 - 3) // 2 == 0
assert all((degree - 3) // 2 >= 1 for degree in range(5, 101))

# In the standard length-minus-delta dimension model, the dominant dimension
# is floor(n/r)-delta.  Hence its codimension in any affine-linear ambient
# dimension is a degree-one quasipolynomial of period dividing r.
for threshold in range(2, 9):
    for degree in range(threshold, 100):
        atomic_types = maximal_threshold_partitions(degree, threshold)
        assert max(map(len, atomic_types)) == threshold_dominant_length(
            degree, threshold
        )
    for degree in range(threshold, 100 - threshold):
        codimension = degree - (degree // threshold - 1)
        shifted_codimension = degree + threshold - (
            (degree + threshold) // threshold - 1
        )
        assert shifted_codimension - codimension == threshold - 1

# The number of primitive components is itself a quasipolynomial.  For fixed
# r its degree is r-1 and a period divides lcm(r,...,2r-1).  Vanishing r-th
# finite differences on every residue class gives an exact finite audit for
# thresholds through five.
for threshold in range(2, 6):
    atoms = tuple(range(threshold, 2 * threshold))
    period = lcm(*atoms)
    limit = (threshold + 1) * period
    counts = [1] + [0] * (limit - 1)
    for atom in atoms:
        for total in range(atom, limit):
            counts[total] += counts[total - atom]
    for residue in range(period):
        values = counts[residue::period]
        for _ in range(threshold):
            values = [right - left for left, right in zip(values, values[1:])]
        assert values and set(values) == {0}

print("PASS: the atoms of {m>=r} are exactly r,...,2r-1")
print("PASS: every allowed contact partition is a collision of atom partitions")
print("PASS: atomic partitions count the primitive threshold components")
print("PASS: the r=2 atoms are precisely the double and triple contacts")
print("PASS: optional Mason regression separates distinct r=2 contact types")
print("PASS: thresholds r>=3 satisfy the stronger bounded-difference abc bound")
print("PASS: component counts and dominant codimensions are quasipolynomial")
print("PASS: nonsurjective dimension floor(n/2)-1 and phase transition at n=4")
