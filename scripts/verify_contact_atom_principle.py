#!/usr/bin/env python3
"""Exact combinatorial and Mason certificates for primitive contact orders."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    affine_difference_mason_defect,
    collision_precedes,
    contact_atom_refinement,
    contact_semigroup_atoms,
    multiplicity_excess,
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

# In the actual reconstruction problem r=2.  The atoms are {2,3}, and the
# Mason defect for any two degree-n types is exactly half their total excess
# above double contact.  It is positive for every distinct pair because the
# unique zero-excess partition is the all-double type.
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

print("PASS: the atoms of {m>=r} are exactly r,...,2r-1")
print("PASS: every allowed contact partition is a collision of atom partitions")
print("PASS: the r=2 atoms are precisely the double and triple contacts")
print("PASS: Mason separates every pair of distinct full-contact types")
print("PASS: thresholds r>=3 satisfy the stronger automatic abc bound")
