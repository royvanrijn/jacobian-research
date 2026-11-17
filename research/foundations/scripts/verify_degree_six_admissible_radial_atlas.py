#!/usr/bin/env python3
"""Verify every radial admissible source type for three quadratic clusters."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.admissible_nodes import normalization_branch_count
from jcsearch.radial_sources import (
    ordered_set_partitions,
    radial_node_indices,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "admissible_radial_atlas_degree6.json"
)

CLUSTERS = ("x", "y", "z")


partitions = ordered_set_partitions(CLUSTERS)
assert len(partitions) == 13
assert sum(len(partition) == 1 for partition in partitions) == 1
assert sum(len(partition) == 2 for partition in partitions) == 6
assert sum(len(partition) == 3 for partition in partitions) == 6


def component_description(
    cluster: str,
    cluster_level: int,
    target_level: int,
) -> list[dict[str, object]]:
    """Components contributed by one quadratic cluster on one target bubble."""

    if cluster_level > target_level:
        return [
            {
                "cluster": cluster,
                "kind": "connector",
                "degree": 2,
                "outer_node_index": 2,
                "inner_node_indices": [2],
                "ramification_contribution": 2,
                "model": "S=V^2",
            }
        ]
    if cluster_level == target_level:
        return [
            {
                "cluster": cluster,
                "kind": "quadratic_tail",
                "degree": 2,
                "outer_node_index": 2,
                "inner_node_indices": [1, 1],
                "ramification_contribution": 2,
                "model": "S=u_i*X*(X-1)",
            }
        ]
    return [
        {
            "cluster": cluster,
            "kind": "identity_strand",
            "sheet": sheet,
            "degree": 1,
            "outer_node_index": 1,
            "inner_node_indices": [1],
            "ramification_contribution": 0,
            "model": "S=V",
        }
        for sheet in (1, 2)
    ]


def node_indices(
    level_by_cluster: dict[str, int],
    after_level: int,
) -> list[int]:
    """Ramification indices over the node after one target level."""

    return list(
        radial_node_indices(
            (2, 2, 2),
            tuple(level_by_cluster[cluster] for cluster in CLUSTERS),
            after_level,
        )
    )


def kummer_branch_count(indices: list[int]) -> int:
    """Normalization branches for the common smoothing relation."""

    return normalization_branch_count(indices)


atlas = []
for partition in partitions:
    level_by_cluster = {
        cluster: level
        for level, block in enumerate(partition)
        for cluster in block
    }

    bubbles = []
    for target_level, active_block in enumerate(partition):
        components = []
        for cluster in CLUSTERS:
            components.extend(
                component_description(
                    cluster,
                    level_by_cluster[cluster],
                    target_level,
                )
            )
        assert sum(component["degree"] for component in components) == 6
        assert all(
            component["ramification_contribution"]
            == 2 * component["degree"] - 2
            for component in components
        )
        bubbles.append(
            {
                "level": target_level,
                "active_clusters": list(active_block),
                "components": components,
                "total_degree": 6,
            }
        )

    nodes = []
    main_indices = [2, 2, 2]
    assert sum(main_indices) == 6
    nodes.append(
        {
            "between": ["main", 0],
            "indices": main_indices,
            "kummer_normalization_branches": kummer_branch_count(main_indices),
        }
    )
    for target_level in range(len(partition) - 1):
        indices = node_indices(level_by_cluster, target_level)
        assert sum(indices) == 6
        nodes.append(
            {
                "between": [target_level, target_level + 1],
                "indices": indices,
                "kummer_normalization_branches": kummer_branch_count(indices),
            }
        )

    # Every node partition has matching degree six on its two branches.
    assert all(sum(node["indices"]) == 6 for node in nodes)

    atlas.append(
        {
            "ordered_partition": [list(block) for block in partition],
            "target_bubble_count": len(partition),
            "bubbles": bubbles,
            "nodes": nodes,
            "total_kummer_branches": (
                math.prod(
                    node["kummer_normalization_branches"] for node in nodes
                )
            ),
        }
    )

# Representative checks for the three possible depths.
equal_scale = next(
    item
    for item in atlas
    if item["ordered_partition"] == [["x", "y", "z"]]
)
assert [node["indices"] for node in equal_scale["nodes"]] == [[2, 2, 2]]
assert equal_scale["total_kummer_branches"] == 4

strict_xyz = next(
    item for item in atlas
    if item["ordered_partition"] == [["x"], ["y"], ["z"]]
)
assert [node["indices"] for node in strict_xyz["nodes"]] == [
    [2, 2, 2],
    [1, 1, 2, 2],
    [1, 1, 1, 1, 2],
]
assert strict_xyz["total_kummer_branches"] == 8

two_then_one = next(
    item for item in atlas
    if item["ordered_partition"] == [["x", "y"], ["z"]]
)
assert [node["indices"] for node in two_then_one["nodes"]] == [
    [2, 2, 2],
    [1, 1, 1, 1, 2],
]
assert two_then_one["total_kummer_branches"] == 4

one_then_two = next(
    item for item in atlas
    if item["ordered_partition"] == [["x"], ["y", "z"]]
)
assert [node["indices"] for node in one_then_two["nodes"]] == [
    [2, 2, 2],
    [1, 1, 2, 2],
]
assert one_then_two["total_kummer_branches"] == 8

atlas_summary = {}
for item in atlas:
    bubble_summaries = []
    for bubble in item["bubbles"]:
        kinds = [component["kind"] for component in bubble["components"]]
        bubble_summaries.append(
            "active={}; connectors={}; identity_strands={}; "
            "quadratic_tails={}; degree={}".format(
                "".join(bubble["active_clusters"]),
                kinds.count("connector"),
                kinds.count("identity_strand"),
                kinds.count("quadratic_tail"),
                bubble["total_degree"],
            )
        )
    partition_name = "|".join(
        "".join(block) for block in item["ordered_partition"]
    )
    atlas_summary[partition_name] = {
        "bubbles": bubble_summaries,
        "node_indices": [
            "+".join(str(index) for index in node["indices"])
            for node in item["nodes"]
        ],
        "node_kummer_branches": [
            node["kummer_normalization_branches"] for node in item["nodes"]
        ],
        "total_kummer_branches": item["total_kummer_branches"],
    }

artifact = {
    "experiment": "degree-six radial admissible-cover atlas",
    "cluster_profile": [2, 2, 2],
    "ordered_partition_count": len(partitions),
    "component_rule": {
        "before_cluster_scale": "one degree-two connector S=V^2",
        "at_cluster_scale": "one degree-two tail S=u_i*X*(X-1)",
        "after_cluster_scale": "two degree-one identity strands",
    },
    "node_rule": {
        "unactivated_cluster": "one node of index two",
        "activated_cluster": "two nodes of index one",
        "normalization_branches": "2^(number_of_index_two_nodes-1)",
    },
    "atlas": atlas_summary,
    "scope": (
        "all radial cones and equality faces for the labelled (2,2,2) "
        "profile; the extension morphism, component selection, and label "
        "descent remain separate"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS radial atlas: all 13 ordered set partitions")
print("PASS every target bubble: source degree six")
print("PASS every source component: Riemann-Hurwitz")
print("PASS every target node: matching index partition of six")
print("PASS Kummer saturation counts on all radial cones and faces")
print("DEGREE_SIX_ADMISSIBLE_RADIAL_ATLAS_PASS")
