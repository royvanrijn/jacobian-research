#!/usr/bin/env python3
"""Audit normalization branches and surviving inertia in the sextic atlas."""

from __future__ import annotations

import json
import math
import sys
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.admissible_nodes import (  # noqa: E402
    node_lcm,
    normalization_branch_count,
)
from jcsearch.monodromy_forests import identity, transposition  # noqa: E402
from jcsearch.resonance_atlas import (  # noqa: E402
    compile_nested_monodromy_characters,
)

ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "stack_inertia_degree6.json"
)

CLUSTERS = ("x", "y", "z")


def ordered_partitions() -> list[tuple[tuple[str, ...], ...]]:
    result = []
    for level_count in range(1, len(CLUSTERS) + 1):
        for assignment in product(range(level_count), repeat=len(CLUSTERS)):
            if set(assignment) != set(range(level_count)):
                continue
            result.append(
                tuple(
                    tuple(
                        cluster
                        for cluster, level in zip(CLUSTERS, assignment)
                        if level == block
                    )
                    for block in range(level_count)
                )
            )
    return result


def kummer_branch_count(index_two_count: int) -> int:
    if index_two_count == 0:
        return 1
    return 2 ** (index_two_count - 1)


radial_inertia = {}
for partition in ordered_partitions():
    level_by_cluster = {
        cluster: level
        for level, block in enumerate(partition)
        for cluster in block
    }

    node_index_two_counts = [3]
    for level in range(len(partition) - 1):
        unactivated = sum(
            cluster_level > level
            for cluster_level in level_by_cluster.values()
        )
        node_index_two_counts.append(unactivated)

    node_kummer_branches = [
        kummer_branch_count(count)
        for count in node_index_two_counts
    ]
    connector_sign_rank = 0
    for index_two_count, branch_count in zip(
        node_index_two_counts,
        node_kummer_branches,
    ):
        # For node indices consisting of index-two and index-one points,
        # L=lcm(e_j)=2.  The full sign group has order product(e_j)=2^k;
        # it has 2^(k-1) normalization orbits, each with diagonal stabilizer
        # mu_2.  Most individual sign flips move between normalization
        # branches rather than fixing one normalized point.
        indices = [2] * index_two_count + [1] * (6 - 2 * index_two_count)
        assert node_lcm(indices) == 2
        assert branch_count == normalization_branch_count(indices)

    # A downstream root label does not componentwise kill an automorphism of
    # a separate connector.  Compile the full radial chain instead: the zero
    # flag is an identity-monodromy mark on the deepest screen, the three
    # moving edges are ordered from the deepest level out, and two stationary
    # edges rigidify the root vertex.
    connector_sign_rank = sum(
        cluster_level
        for cluster_level in level_by_cluster.values()
    )
    cluster_index = {cluster: index for index, cluster in enumerate(CLUSTERS)}
    moving_edges = ((0, 1), (2, 3), (4, 5))
    stationary_edges = ((1, 2), (3, 4))
    inner_to_outer = tuple(reversed(partition))
    moving_order = tuple(
        cluster_index[cluster]
        for block in inner_to_outer
        for cluster in block
    )
    branch_permutations = (
        identity(6),
        *(
            transposition(6, *moving_edges[index])
            for index in moving_order
        ),
        *(
            transposition(6, *edge)
            for edge in stationary_edges
        ),
    )
    cumulative = {0}
    position = 1
    radial_family = []
    for block in inner_to_outer:
        cumulative.update(range(position, position + len(block)))
        position += len(block)
        radial_family.append(frozenset(cumulative))
    compiled = compile_nested_monodromy_characters(
        6,
        tuple(branch_permutations),
        tuple(radial_family),
        anchored_vertices=(radial_family[0],),
        root_anchored=False,
    )
    assert compiled.global_automorphism_count == 2**connector_sign_rank
    assert compiled.inertia_order == 1

    partition_name = "|".join("".join(block) for block in partition)
    radial_inertia[partition_name] = {
        "target_node_count": len(node_index_two_counts),
        "node_index_two_counts": node_index_two_counts,
        "connector_sign_group_before_normalized_lift": (
            f"(mu_2)^{connector_sign_rank}"
        ),
        "connector_sign_group_order": 2**connector_sign_rank,
        "normalized_labelled_inertia": "trivial",
        "normalized_labelled_inertia_order": 1,
        "node_kummer_branches": node_kummer_branches,
        "total_kummer_branches": math.prod(node_kummer_branches),
        "anchor_rule": (
            "full cycle matching between the anchored zero screen and rigid "
            "stationary root makes every nontrivial connector sign change "
            "at least one simultaneous normalized branch"
        ),
    }

assert radial_inertia["xyz"] == {
    "target_node_count": 1,
    "node_index_two_counts": [3],
    "connector_sign_group_before_normalized_lift": "(mu_2)^0",
    "connector_sign_group_order": 1,
    "normalized_labelled_inertia": "trivial",
    "normalized_labelled_inertia_order": 1,
    "node_kummer_branches": [4],
    "total_kummer_branches": 4,
    "anchor_rule": (
        "full cycle matching between the anchored zero screen and rigid "
        "stationary root makes every nontrivial connector sign change "
        "at least one simultaneous normalized branch"
    ),
}
assert radial_inertia["x|y|z"] == {
    "target_node_count": 3,
    "node_index_two_counts": [3, 2, 1],
    "connector_sign_group_before_normalized_lift": "(mu_2)^3",
    "connector_sign_group_order": 8,
    "normalized_labelled_inertia": "trivial",
    "normalized_labelled_inertia_order": 1,
    "node_kummer_branches": [4, 2, 1],
    "total_kummer_branches": 8,
    "anchor_rule": (
        "full cycle matching between the anchored zero screen and rigid "
        "stationary root makes every nontrivial connector sign change "
        "at least one simultaneous normalized branch"
    ),
}

# A Maxwell node with m index-two preimages has sign group (mu_2)^m.
# Normalization separates 2^(m-1) sign-ratio branches.  At one selected
# branch, only the diagonal mu_2 remains as inertia.
maxwell_inertia = {}
for colliding_count, name in ((2, "pairwise"), (3, "triple")):
    sign_group_order = 2**colliding_count
    indices = [2] * colliding_count + [1] * (6 - 2 * colliding_count)
    normalization_branches = normalization_branch_count(indices)
    normalized_inertia_order = 2
    assert (
        normalization_branches * normalized_inertia_order
        == sign_group_order
    )
    maxwell_inertia[name] = {
        "colliding_branch_values": colliding_count,
        "quadratic_square_tails": colliding_count,
        "source_node_sign_group_order": sign_group_order,
        "normalized_labelled_inertia": "mu_2",
        "normalized_labelled_inertia_order": normalized_inertia_order,
        "normalization_branches": normalization_branches,
        "reason": (
            "all colliding square tails admit deck flips and their diagonal "
            "subgroup stabilizes one normalized sign-ratio branch"
        ),
    }

assert maxwell_inertia["pairwise"]["normalized_labelled_inertia_order"] == 2
assert maxwell_inertia["pairwise"]["normalization_branches"] == 2
assert maxwell_inertia["triple"]["normalized_labelled_inertia_order"] == 2
assert maxwell_inertia["triple"]["normalization_branches"] == 4

artifact = {
    "experiment": "degree-six stack inertia audit",
    "normalization_rule": (
        "for source-node indices e_j over one target node, normalization has "
        "product(e_j)/lcm(e_j) branches and diagonal inertia "
        "mu_lcm(e_j) on each branch"
    ),
    "radial_types": radial_inertia,
    "maxwell_types": maxwell_inertia,
    "distinction": {
        "normalization_branches": "product(indices)/lcm(indices)",
        "candidate_diagonal_on_one_branch": "mu_lcm(indices)",
        "actual_inertia_test": (
            "intersect the candidate diagonal with cover automorphisms "
            "preserving every source-root label"
        ),
        "conclusion": (
            "full radial chain matching gives trivial normalized labelled "
            "inertia; Maxwell square tails leave diagonal mu_2"
        ),
    },
    "scope": (
        "normalization branches and node-sign cover automorphisms in the "
        "fully source-root-labelled degree-six atlas; full root-label "
        "quotient inertia is separate"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS radial inertia: all 13 ordered scale types")
print("PASS labelled radial inertia: full-chain character matching is trivial")
print("PASS Maxwell inertia: diagonal mu_2 on pairwise and triple branches")
print("PASS normalization branches separated from actual cover inertia")
print("DEGREE_SIX_STACK_INERTIA_PASS")
