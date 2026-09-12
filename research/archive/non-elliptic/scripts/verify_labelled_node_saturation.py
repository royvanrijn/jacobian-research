#!/usr/bin/env python3
"""Verify the general labelled-node normalization and inertia rule."""

from __future__ import annotations

import json
import math
import sys
from itertools import permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.admissible_nodes import (
    diagonal_phase_group,
    inertia_character_count,
    node_lcm,
    normalization_branch_count,
    phase_group,
)

ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "labelled_node_saturation.json"
)


# Exhaust the local arithmetic well beyond the indices needed by the sextic
# experiment.  The diagonal q -> zeta_L q has order L and its quotient in
# prod_j mu_{e_j} is the set of normalized branches.
profiles_checked = 0
for length in range(1, 5):
    for indices in product(range(1, 7), repeat=length):
        ambient = phase_group(indices)
        diagonal = diagonal_phase_group(indices)
        common_index = node_lcm(indices)
        assert len(ambient) == math.prod(indices)
        assert len(diagonal) == common_index
        assert normalization_branch_count(indices) == (
            len(ambient) // len(diagonal)
        )

        # Reordering labelled source nodes transports the phase quotient.
        for permutation in permutations(range(length)):
            permuted_indices = tuple(indices[index] for index in permutation)
            transported_diagonal = {
                tuple(value[index] for index in permutation)
                for value in diagonal
            }
            assert transported_diagonal == diagonal_phase_group(
                permuted_indices
            )
            assert normalization_branch_count(permuted_indices) == (
                normalization_branch_count(indices)
            )
        profiles_checked += 1


# If every cyclic component deck transformation is allowed, its character
# image is the full phase group and the inertia of one normalized lift is
# exactly the diagonal mu_L.  If one labelled natural tail anchors a phase,
# the intersection with the diagonal is trivial for equal indices.
maxwell_examples = {}
anchored_examples = {}
for index in range(2, 7):
    for source_node_count in range(2, 6):
        indices = (index,) * source_node_count
        full_characters = phase_group(indices)
        maxwell_inertia = inertia_character_count(
            indices,
            full_characters,
        )
        assert maxwell_inertia == index

        anchored_characters = (
            (0,) + character
            for character in product(
                range(index),
                repeat=source_node_count - 1,
            )
        )
        anchored_inertia = inertia_character_count(
            indices,
            anchored_characters,
        )
        assert anchored_inertia == 1

        key = f"e={index},r={source_node_count}"
        maxwell_examples[key] = {
            "normalization_branches": index ** (source_node_count - 1),
            "labelled_inertia_order": maxwell_inertia,
        }
        anchored_examples[key] = {
            "normalization_branches": index ** (source_node_count - 1),
            "labelled_inertia_order": anchored_inertia,
        }


# With unequal indices, one constrained phase leaves L/e_anchor diagonal
# elements rather than automatically killing the whole diagonal.
unequal_anchor_checks = 0
unequal_anchor_examples = {}
for source_node_count in range(2, 5):
    for indices in product(range(1, 7), repeat=source_node_count):
        common_index = node_lcm(indices)
        for anchor in range(source_node_count):
            constrained_characters = (
                character
                for character in phase_group(indices)
                if character[anchor] == 0
            )
            inertia = inertia_character_count(
                indices,
                constrained_characters,
            )
            assert inertia == common_index // indices[anchor]
            unequal_anchor_checks += 1
            if indices == (2, 3):
                unequal_anchor_examples[f"anchor={anchor}"] = inertia
assert unequal_anchor_examples == {
    "anchor=0": 3,
    "anchor=1": 2,
}


# The stack quotient underlying H2 is independent of the compactification:
# for S_{n-1} fixing one label inside S_n, the left-coset set has n elements.
quotient_degrees = {}
for label_count in range(2, 9):
    symmetric_group = tuple(permutations(range(label_count)))
    stabilizer = tuple(
        permutation
        for permutation in symmetric_group
        if permutation[0] == 0
    )
    assert len(symmetric_group) == math.factorial(label_count)
    assert len(stabilizer) == math.factorial(label_count - 1)
    assert len(symmetric_group) // len(stabilizer) == label_count
    quotient_degrees[str(label_count)] = label_count


artifact = {
    "experiment": "general labelled admissible-node saturation",
    "profiles_checked": profiles_checked,
    "profile_range": {
        "source_node_count": [1, 4],
        "indices": [1, 6],
    },
    "normalization_theorem": {
        "common_index": "L=lcm(e_j)",
        "phase_group": "product_j mu_(e_j)",
        "diagonal_reparametrization": (
            "mu_L -> product_j mu_(e_j), "
            "zeta maps to (zeta^(L/e_j))_j"
        ),
        "branch_count": "product_j(e_j)/L",
    },
    "inertia_theorem": {
        "character": (
            "chi: Aut_labelled(cover) -> product_j mu_(e_j)"
        ),
        "inertia_on_a_normalized_lift": "chi^(-1)(diagonal mu_L)",
        "warning": (
            "the diagonal phase subgroup is a reparametrization group; "
            "it is inertia only to the extent that labelled cover "
            "automorphisms realize it"
        ),
    },
    "families": {
        "unanchored_cyclic_tails": maxwell_examples,
        "one_labelled_anchor": anchored_examples,
        "unequal_single_constrained_phase": unequal_anchor_examples,
    },
    "unequal_anchor_checks": unequal_anchor_checks,
    "equivariance": (
        "the phase quotient and inertia test commute with every permutation "
        "of labelled source nodes"
    ),
    "corrected_graph_quotient": {
        "label_count_to_degree": quotient_degrees,
        "rule": (
            "for any S_n-equivariant labelled normalized graph Gamma, "
            "[Gamma/S_(n-1)] -> [Gamma/S_n] is finite representable etale "
            "of degree n"
        ),
    },
    "scope": (
        "local node normalization, label-preserving inertia, permutation "
        "equivariance, and the formal subgroup quotient; construction of "
        "the global corrected graph and its explicit recursive fan remain "
        "separate"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print(f"PASS node saturation: {profiles_checked} index profiles")
print("PASS label equivariance: every source-node permutation")
print("PASS inertia rule: unanchored diagonal versus labelled anchor")
print(
    "PASS unequal constrained phases: "
    f"{unequal_anchor_checks} anchored profiles"
)
print("PASS corrected graph quotient: degrees 2 through 8")
print("LABELLED_NODE_SATURATION_PASS")
