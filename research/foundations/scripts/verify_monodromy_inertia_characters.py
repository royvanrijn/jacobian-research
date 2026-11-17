#!/usr/bin/env python3
"""Verify deck centralizers and labelled resonance-node inertia."""

from __future__ import annotations

import json
import math
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.monodromy_forests import (
    centralizer,
    fixes_labels,
    forest_components,
    reduced_cycle_factorizations,
    standard_cycle,
    transposition,
)
from jcsearch.admissible_nodes import (
    compose_character_tables,
    simultaneous_inertia_character_count,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "monodromy_inertia_characters.json"
)


def simple_resonance_inertia_order(
    node_indices: tuple[int, ...],
    zero_flag_on_bubble: bool,
) -> int:
    """Compute diagonal inertia from unanchored quadratic components."""

    common_index = math.lcm(*node_indices)
    constrained_indices = tuple(
        index
        for index in node_indices
        if index != 2 or zero_flag_on_bubble
    )
    constrained_lcm = (
        math.lcm(*constrained_indices)
        if constrained_indices
        else 1
    )
    assert common_index % constrained_lcm == 0
    return common_index // constrained_lcm


# A connected simple-branch tree generates S_k.  Its deck group is the
# centralizer: mu_2 in degree two and trivial from degree three onward.
tree_centralizer_checks = 0
for degree in range(2, 7):
    for edge_sequence in reduced_cycle_factorizations(degree):
        generators = tuple(
            transposition(degree, *edge)
            for edge in edge_sequence
        )
        deck_group = centralizer(generators)
        expected_order = 2 if degree == 2 else 1
        assert len(deck_group) == expected_order
        if degree == 2:
            assert sum(fixes_labels(deck, (0,)) for deck in deck_group) == 1
        tree_centralizer_checks += 1


# A connector with only its two node flags has cyclic monodromy and cyclic
# deck group.  A label on that connector kills it directly.  A label on a
# different downstream component requires the full-chain character test in
# verify_recursive_resonance_atlas.py.
connector_checks = {}
for degree in range(2, 9):
    connector_deck_group = centralizer((standard_cycle(degree),))
    anchored_group = tuple(
        deck
        for deck in connector_deck_group
        if fixes_labels(deck, (0,))
    )
    assert len(connector_deck_group) == degree
    assert len(anchored_group) == 1
    connector_checks[str(degree)] = {
        "unanchored_deck_order": degree,
        "one_sheet_anchored_order": 1,
    }


# Exhaust every collision subforest through degree six.  Away from the
# target zero flag, each size-two component supplies a flip.  The inertia on
# one normalized lift is the intersection with the diagonal node phase.  If
# the zero flag lies on the bubble, its labelled regular fiber kills all
# those flips.
partition_occurrences = Counter()
unanchored_inertia_occurrences = Counter()
anchored_inertia_occurrences = Counter()
collision_nodes_checked = 0
for degree in range(2, 7):
    for edge_sequence in reduced_cycle_factorizations(degree):
        for mask in range(1 << (degree - 1)):
            selected_edges = tuple(
                edge
                for position, edge in enumerate(edge_sequence)
                if mask & (1 << position)
            )
            node_indices = tuple(
                sorted(
                    len(component)
                    for component in forest_components(
                        degree,
                        selected_edges,
                    )
                )
            )
            unanchored_order = simple_resonance_inertia_order(
                node_indices,
                zero_flag_on_bubble=False,
            )
            anchored_order = simple_resonance_inertia_order(
                node_indices,
                zero_flag_on_bubble=True,
            )
            assert unanchored_order in (1, 2)
            assert anchored_order == 1
            if 2 not in node_indices:
                assert unanchored_order == 1

            key = "+".join(map(str, node_indices))
            partition_occurrences[key] += 1
            unanchored_inertia_occurrences[str(unanchored_order)] += 1
            anchored_inertia_occurrences[str(anchored_order)] += 1
            collision_nodes_checked += 1


# The same partitions have different inertia depending only on whether the
# labelled zero fiber anchors their quadratic components.
pairwise_maxwell = (1, 1, 2, 2)
triple_maxwell = (2, 2, 2)
two_edge_caustic = (1, 1, 1, 3)
assert simple_resonance_inertia_order(pairwise_maxwell, False) == 2
assert simple_resonance_inertia_order(pairwise_maxwell, True) == 1
assert simple_resonance_inertia_order(triple_maxwell, False) == 2
assert simple_resonance_inertia_order(triple_maxwell, True) == 1
assert simple_resonance_inertia_order(two_edge_caustic, False) == 1


# Higher-codimension inertia is a simultaneous character intersection, not a
# product of generic node orders.  The first example recovers independent
# pair/triple Maxwell roots.  The second has one connector rotation acting
# at two nodes and therefore retains one mu_2, not their naive product.
pair_triple_profiles = (pairwise_maxwell, triple_maxwell)
pair_triple_moduli = tuple(
    index
    for profile in pair_triple_profiles
    for index in profile
)
pair_flip = (0, 0, 1, 1, 0, 0, 0)
triple_flip = (0, 0, 0, 0, 1, 1, 1)
pair_triple_characters = compose_character_tables(
    pair_triple_moduli,
    (
        ((0,) * 7, pair_flip),
        ((0,) * 7, triple_flip),
    ),
)
assert len(pair_triple_characters) == 4
assert simultaneous_inertia_character_count(
    pair_triple_profiles,
    pair_triple_characters,
) == 4

coupled_profiles = ((2, 2), (2, 2))
coupled_characters = (
    (0, 0, 0, 0),
    (1, 1, 1, 1),
)
assert simultaneous_inertia_character_count(
    coupled_profiles,
    coupled_characters,
) == 2

off_diagonal_characters = (
    (0, 0, 0, 0),
    (1, 0, 1, 0),
)
assert simultaneous_inertia_character_count(
    coupled_profiles,
    off_diagonal_characters,
) == 1


artifact = {
    "experiment": "monodromy centralizers and labelled node inertia",
    "degree_range": [2, 8],
    "tree_centralizer_checks": tree_centralizer_checks,
    "connector_checks": connector_checks,
    "collision_nodes_checked": collision_nodes_checked,
    "distinct_node_partitions": len(partition_occurrences),
    "node_partition_occurrences": dict(sorted(partition_occurrences.items())),
    "unanchored_inertia_occurrences": dict(
        sorted(unanchored_inertia_occurrences.items())
    ),
    "anchored_inertia_occurrences": dict(
        sorted(anchored_inertia_occurrences.items())
    ),
    "rules": {
        "simple_branch_component": (
            "a connected k-sheet transposition tree has deck centralizer "
            "mu_2 for k=2 and trivial centralizer for k>=3"
        ),
        "power_connector": (
            "an unanchored k-cycle connector has cyclic deck group mu_k; "
            "one labelled sheet on that connector kills it directly"
        ),
        "resonance_node": (
            "only unanchored size-two forest components contribute phases; "
            "intersect their flips with the diagonal mu_lcm"
        ),
        "zero_flag_anchor": (
            "if the target zero flag lies on the bubble, the fully labelled "
            "regular zero fiber kills every quadratic flip"
        ),
    },
    "degree_six": {
        "pairwise_maxwell_unanchored": 2,
        "pairwise_same_partition_anchored": 1,
        "triple_maxwell_unanchored": 2,
        "triple_same_partition_anchored": 1,
        "two_edge_caustic": 1,
    },
    "simultaneous_inertia": {
        "independent_pair_triple_maxwell": 4,
        "one_flip_coupled_across_two_nodes": 2,
        "off_diagonal_flip_rejected": 1,
        "rule": (
            "compose local character tables with multiplicity, then require "
            "the phase segment at every target node to lie in its diagonal "
            "mu_lcm"
        ),
    },
    "consequence": (
        "local radial and simple-resonance character tables are determined "
        "by monodromy centralizers; the recursive resonance checker performs "
        "the required full-chain matching and verifies radial order "
        "product_j L_j/M_j, which can be nontrivial for unequal "
        "multiplicities"
    ),
    "scope": (
        "deck and node inertia for radial connectors and simple finite "
        "branch resonance components; higher non-simple branch profiles "
        "use the same centralizer-character algorithm but are not "
        "exhausted here"
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
    "PASS tree deck centralizers: "
    f"{tree_centralizer_checks} polynomial monodromy trees"
)
print("PASS cyclic connectors: degrees two through eight")
print(
    "PASS labelled resonance inertia: "
    f"{collision_nodes_checked} collision nodes"
)
print("PASS degree-six radial/Maxwell/caustic distinction")
print("PASS simultaneous higher-codimension character intersection")
print("MONODROMY_INERTIA_CHARACTERS_PASS")
