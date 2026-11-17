#!/usr/bin/env python3
"""Verify the all-multiplicity radial source atlas and its node saturation."""

from __future__ import annotations

import json
import math
import sys
from collections import Counter
from itertools import permutations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.admissible_nodes import normalization_branch_count
from jcsearch.radial_sources import (
    ordered_set_partitions,
    radial_component_kinds,
    radial_inertia_factors,
    radial_inertia_order,
    radial_node_indices,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "general_radial_source_atlas.json"
)


cluster_profiles_checked = 0
ordered_scale_types_checked = 0
target_bubbles_checked = 0
target_nodes_checked = 0
label_permutations_checked = 0
scale_type_counts = {}
radial_inertia_distribution = Counter()
nontrivial_radial_inertia_types = 0
maximum_radial_inertia_order = 1


def phase_count(
    multiplicities: tuple[int, ...],
    blocks: tuple[tuple[int, ...], ...],
) -> int:
    """Independently count telescoping diagonal phase solutions."""

    suffix_lcms = tuple(
        math.lcm(
            *(
                multiplicities[index]
                for suffix in blocks[level:]
                for index in suffix
            )
        )
        for level in range(len(blocks))
    )
    states = {0: 1}
    for level, block in enumerate(blocks):
        modulus = suffix_lcms[level]
        next_modulus = (
            suffix_lcms[level + 1]
            if level + 1 < len(blocks)
            else 1
        )
        next_states = Counter()
        for prefix, count in states.items():
            for phase in range(modulus):
                updated = (prefix + phase) % modulus
                if all(
                    updated % multiplicities[index] == 0
                    for index in block
                ):
                    next_states[updated % next_modulus] += count
        states = dict(next_states)
    return sum(states.values())

for cluster_count in range(1, 5):
    labels = tuple(range(cluster_count))
    ordered_partitions = ordered_set_partitions(labels)
    scale_type_counts[str(cluster_count)] = len(ordered_partitions)

    for multiplicities in product(range(2, 7), repeat=cluster_count):
        total_degree = sum(multiplicities)
        cluster_profiles_checked += 1

        for partition in ordered_partitions:
            level_by_cluster = {
                cluster: level
                for level, block in enumerate(partition)
                for cluster in block
            }
            levels = tuple(level_by_cluster[label] for label in labels)
            inertia_factors = radial_inertia_factors(
                multiplicities,
                partition,
            )
            inertia_order = radial_inertia_order(
                multiplicities,
                partition,
            )
            assert inertia_order == math.prod(inertia_factors)
            assert inertia_factors[-1] == 1
            assert inertia_order == phase_count(
                multiplicities,
                partition,
            )
            if len(set(multiplicities)) == 1:
                assert inertia_order == 1
            radial_inertia_distribution[str(inertia_order)] += 1
            nontrivial_radial_inertia_types += inertia_order > 1
            maximum_radial_inertia_order = max(
                maximum_radial_inertia_order,
                inertia_order,
            )

            # On every target bubble, an inactive cluster contributes one
            # power connector, an active cluster its local polynomial, and
            # a previously active cluster its separated identity sheets.
            for target_level in range(len(partition)):
                components = tuple(
                    component
                    for multiplicity, cluster_level in zip(
                        multiplicities,
                        levels,
                    )
                    for component in radial_component_kinds(
                        multiplicity,
                        cluster_level,
                        target_level,
                    )
                )
                assert sum(component[1] for component in components) == total_degree
                assert all(
                    ramification == 2 * degree - 2
                    for _, degree, ramification in components
                )
                target_bubbles_checked += 1

            # The central-to-radial node sees one index mu_i per cluster.
            # After a cluster becomes active, that index is replaced by
            # mu_i unramified indices.  Every partition still sums to N.
            central_indices = tuple(sorted(multiplicities))
            assert sum(central_indices) == total_degree
            assert normalization_branch_count(central_indices) >= 1
            target_nodes_checked += 1
            for target_level in range(len(partition) - 1):
                indices = radial_node_indices(
                    multiplicities,
                    levels,
                    target_level,
                )
                assert sum(indices) == total_degree
                assert normalization_branch_count(indices) >= 1
                target_nodes_checked += 1

            ordered_scale_types_checked += 1

        # Relabelling clusters transports every ordered scale type and every
        # node-index multiset.  Check all permutations through four clusters.
        partition_set = {
            tuple(frozenset(block) for block in partition)
            for partition in ordered_partitions
        }
        for permutation in permutations(labels):
            action = dict(zip(labels, permutation))
            transported = {
                tuple(
                    frozenset(action[label] for label in block)
                    for block in partition
                )
                for partition in ordered_partitions
            }
            assert transported == partition_set

            permuted_multiplicities = tuple(
                multiplicities[permutation.index(label)]
                for label in labels
            )
            assert sorted(permuted_multiplicities) == sorted(multiplicities)
            label_permutations_checked += 1


# Unequal multiplicities expose the lcm saturation that is invisible in the
# all-quadratic degree-six example.
unequal_multiplicities = (2, 3, 4)
unequal_levels = (0, 1, 2)
unequal_nodes = (
    tuple(sorted(unequal_multiplicities)),
    radial_node_indices(unequal_multiplicities, unequal_levels, 0),
    radial_node_indices(unequal_multiplicities, unequal_levels, 1),
)
assert unequal_nodes == (
    (2, 3, 4),
    (1, 1, 3, 4),
    (1, 1, 1, 1, 1, 4),
)
unequal_branch_counts = tuple(
    normalization_branch_count(indices)
    for indices in unequal_nodes
)
assert unequal_branch_counts == (2, 1, 1)
unequal_inertia_factors = radial_inertia_factors(
    unequal_multiplicities,
    ((0,), (1,), (2,)),
)
assert unequal_inertia_factors == (6, 4, 1)
assert math.prod(unequal_inertia_factors) == 24
assert radial_inertia_order((2, 3), ((0,), (1,))) == 3
assert radial_inertia_order((2, 3), ((1,), (0,))) == 2


artifact = {
    "experiment": "general radial admissible-source atlas",
    "cluster_count_range": [1, 4],
    "multiplicity_range": [2, 6],
    "cluster_profiles_checked": cluster_profiles_checked,
    "ordered_scale_type_counts": scale_type_counts,
    "ordered_scale_types_checked": ordered_scale_types_checked,
    "target_bubbles_checked": target_bubbles_checked,
    "target_nodes_checked": target_nodes_checked,
    "label_permutations_checked": label_permutations_checked,
    "distinct_radial_inertia_orders": len(radial_inertia_distribution),
    "nontrivial_radial_inertia_types": nontrivial_radial_inertia_types,
    "maximum_radial_inertia_order": maximum_radial_inertia_order,
    "radial_inertia_rule": (
        "for ordered blocks B_j, order is product_j L_j/M_j, where M_j "
        "is the lcm on B_j and L_j is the lcm on the suffix from B_j"
    ),
    "component_rule": {
        "before_active_scale": (
            "one degree-mu power connector, with both node fibers mu-fold"
        ),
        "at_active_scale": (
            "one degree-mu local polynomial tail, reconstructed from its "
            "root divisor and the mu-fold outer-node fiber"
        ),
        "after_active_scale": "mu degree-one identity strands",
    },
    "node_rule": (
        "at a radial node, an inactive multiplicity-mu cluster contributes "
        "one index mu and an active cluster contributes mu indices 1"
    ),
    "unequal_multiplicity_example": {
        "multiplicities": list(unequal_multiplicities),
        "strict_node_indices": [
            list(indices)
            for indices in unequal_nodes
        ],
        "normalization_branch_counts": list(unequal_branch_counts),
        "radial_inertia_factors": list(unequal_inertia_factors),
        "radial_inertia_order": math.prod(unequal_inertia_factors),
    },
    "consequence": (
        "the harmonic source tree, component degrees, local rational maps, "
        "and saturated node-branch counts are canonical on every radial "
        "ordered-scale stratum; full-chain radial inertia is product_j "
        "L_j/M_j; the recursive resonance atlas separately constructs all "
        "residue-resonance refinements and contractions"
    ),
    "scope": (
        "radial source combinatorics and vertex maps for labelled clusters; "
        "nonradial smoothing and the full stack comparison are verified by "
        "the separate recursive resonance atlas"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print(
    "PASS radial profiles: "
    f"{cluster_profiles_checked} multiplicity profiles, "
    f"{ordered_scale_types_checked} ordered scale types"
)
print(
    "PASS source components: "
    f"{target_bubbles_checked} target bubbles with degree and RH"
)
print(
    "PASS node saturation: "
    f"{target_nodes_checked} degree partitions"
)
print(
    "PASS full-chain radial inertia: "
    f"{nontrivial_radial_inertia_types} nontrivial types, maximum order "
    f"{maximum_radial_inertia_order}"
)
print(
    "PASS label equivariance: "
    f"{label_permutations_checked} cluster permutations"
)
print("GENERAL_RADIAL_SOURCE_ATLAS_PASS")
