#!/usr/bin/env python3
"""Verify the forest rule for every reduced polynomial factorization to d=6."""

from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.admissible_nodes import normalization_branch_count
from jcsearch.monodromy_forests import (
    cycle_lengths,
    forest_components,
    is_tree,
    permutation_product,
    reduced_cycle_factorizations,
    standard_cycle,
    transposition,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "polynomial_monodromy_forests.json"
)


factorization_counts = {}
collision_subsets_checked = 0
nested_pairs_checked = 0
riemann_hurwitz_components_checked = 0
representative_degree_six = None

for degree in range(2, 7):
    edges = tuple(combinations(range(degree), 2))
    edge_permutations = {
        edge: transposition(degree, *edge)
        for edge in edges
    }
    target_cycle = standard_cycle(degree)
    factorizations = reduced_cycle_factorizations(degree)

    for edge_sequence in factorizations:
        permutations = tuple(edge_permutations[edge] for edge in edge_sequence)
        assert permutation_product(permutations, degree) == target_cycle
        assert len(set(edge_sequence)) == degree - 1
        assert is_tree(degree, edge_sequence)

    expected_count = degree ** (degree - 2)
    assert len(factorizations) == expected_count
    factorization_counts[str(degree)] = len(factorizations)

    for edge_sequence in factorizations:
        subset_data = {}
        for mask in range(1 << (degree - 1)):
            selected_edges = tuple(
                edge
                for index, edge in enumerate(edge_sequence)
                if mask & (1 << index)
            )
            selected_permutations = tuple(
                edge_permutations[edge]
                for edge in selected_edges
            )
            boundary_monodromy = permutation_product(
                selected_permutations,
                degree,
            )
            components = forest_components(degree, selected_edges)
            component_sizes = tuple(
                sorted(len(component) for component in components)
            )
            assert cycle_lengths(boundary_monodromy) == component_sizes

            # A component with k sheets contains k-1 selected simple branch
            # transpositions and one k-cycle at the attaching node.
            for component in components:
                component_set = set(component)
                internal_edges = sum(
                    set(edge) <= component_set
                    for edge in selected_edges
                )
                assert internal_edges == len(component) - 1
                finite_ramification = internal_edges
                node_ramification = len(component) - 1
                assert (
                    finite_ramification + node_ramification
                    == 2 * len(component) - 2
                )
                riemann_hurwitz_components_checked += 1

            branch_count = normalization_branch_count(component_sizes)
            assert branch_count >= 1
            subset_data[mask] = components
            collision_subsets_checked += 1

        # Nested collision subsets give nested forests: every source
        # component for the smaller subset lies in one for the larger.
        for smaller_mask in subset_data:
            for larger_mask in subset_data:
                if smaller_mask & ~larger_mask:
                    continue
                smaller_components = subset_data[smaller_mask]
                larger_components = subset_data[larger_mask]
                assert all(
                    any(
                        set(smaller) <= set(larger)
                        for larger in larger_components
                    )
                    for smaller in smaller_components
                )
                nested_pairs_checked += 1

        if degree == 6 and representative_degree_six is None:
            matching_positions = next(
                (
                    positions
                    for positions in combinations(range(5), 3)
                    if len(
                        {
                            vertex
                            for position in positions
                            for vertex in edge_sequence[position]
                        }
                    )
                    == 6
                ),
                None,
            )
            if matching_positions is not None:
                pair_positions = matching_positions[:2]
                triple_edges = tuple(
                    edge_sequence[position]
                    for position in matching_positions
                )
                pair_edges = triple_edges[:2]
                triple_indices = tuple(
                    sorted(
                        len(component)
                        for component in forest_components(6, triple_edges)
                    )
                )
                pair_indices = tuple(
                    sorted(
                        len(component)
                        for component in forest_components(6, pair_edges)
                    )
                )
                assert triple_indices == (2, 2, 2)
                assert pair_indices == (1, 1, 2, 2)
                assert normalization_branch_count(triple_indices) == 4
                assert normalization_branch_count(pair_indices) == 2

                adjacent_positions = next(
                    positions
                    for positions in combinations(range(5), 2)
                    if len(
                        set(edge_sequence[positions[0]])
                        | set(edge_sequence[positions[1]])
                    )
                    == 3
                )
                adjacent_edges = tuple(
                    edge_sequence[position]
                    for position in adjacent_positions
                )
                adjacent_indices = tuple(
                    sorted(
                        len(component)
                        for component in forest_components(6, adjacent_edges)
                    )
                )
                assert adjacent_indices == (1, 1, 1, 3)

                representative_degree_six = {
                    "factorization_edges": [
                        list(edge)
                        for edge in edge_sequence
                    ],
                    "pairwise_maxwell": {
                        "selected_edge_positions": list(pair_positions),
                        "node_indices": list(pair_indices),
                        "normalization_branches": 2,
                    },
                    "triple_maxwell": {
                        "selected_edge_positions": list(matching_positions),
                        "node_indices": list(triple_indices),
                        "normalization_branches": 4,
                    },
                    "two-edge_caustic": {
                        "selected_edge_positions": list(adjacent_positions),
                        "node_indices": list(adjacent_indices),
                        "normalization_branches": 1,
                    },
                }


assert representative_degree_six is not None
assert sum(factorization_counts.values()) == 1441

artifact = {
    "experiment": "polynomial monodromy forest source rule",
    "degree_range": [2, 6],
    "reduced_factorization_counts": factorization_counts,
    "expected_count_formula": "degree^(degree-2)",
    "factorizations_checked": sum(factorization_counts.values()),
    "collision_subsets_checked": collision_subsets_checked,
    "nested_subset_pairs_checked": nested_pairs_checked,
    "riemann_hurwitz_components_checked": (
        riemann_hurwitz_components_checked
    ),
    "theorem": {
        "global_monodromy": (
            "the d-1 simple finite branch transpositions of a degree-d "
            "polynomial factor a d-cycle and form the edges of a labelled "
            "tree"
        ),
        "collision_subset": (
            "source components over its target bubble are the connected "
            "components of the selected edge subforest"
        ),
        "node_indices": (
            "the source-node expansion indices are the subforest component "
            "sizes, equivalently the cycle lengths of the boundary "
            "monodromy product"
        ),
        "nested_compatibility": (
            "subforest components refine under inclusion of nested branch "
            "subsets"
        ),
    },
    "degree_six_recovery": representative_degree_six,
    "consequence": (
        "the source dual tree, component degrees, node indices, and "
        "normalization branch counts on every nested simple-branch "
        "resonance stratum are determined by the labelled polynomial "
        "monodromy tree"
    ),
    "scope": (
        "combinatorial source enhancement for simple finite branch cycles; "
        "the separate recursive resonance atlas verifies algebraic "
        "target-flag positions, residue equations, and the global stack"
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
    "PASS polynomial factorizations: "
    f"{sum(factorization_counts.values())} through degree six"
)
print(
    "PASS collision forests: "
    f"{collision_subsets_checked} branch subsets"
)
print(
    "PASS nested compatibility: "
    f"{nested_pairs_checked} subset pairs"
)
print("PASS degree-six Maxwell and caustic node partitions")
print("POLYNOMIAL_MONODROMY_FORESTS_PASS")
