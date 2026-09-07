#!/usr/bin/env python3
"""Validate persistent K3 surface/fibration/character object graphs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GRAPH = (
    ROOT
    / "elkies-k3/data/object-graph/determinant-948-surface-graph-v1.json"
)

CERTIFIED = {"CERTIFIED_EXACT", "CERTIFIED_PARTIAL"}
KINDS = {
    "surface",
    "transcendental_lattice",
    "neron_severi_lattice",
    "clifford_order",
    "marking_group",
    "period_point",
    "fibration",
    "frame_lattice",
    "carrier",
    "character",
    "field",
    "surface_cover",
    "mw_block",
    "divisor",
    "utility_profile",
}
STATUSES = CERTIFIED | {"DEFINITIONAL", "UNKNOWN", "DIAGNOSTIC_ONLY"}

ALLOWED_ENDPOINTS = {
    "HAS_TRANSCENDENTAL": {("surface", "transcendental_lattice")},
    "HAS_NERON_SEVERI": {("surface", "neron_severi_lattice")},
    "HAS_CLIFFORD_ORDER": {("surface", "clifford_order")},
    "HAS_STABLE_MARKING_GROUP": {("surface", "marking_group")},
    "HAS_PERIOD_POINT": {("surface", "period_point")},
    "HAS_FIBRATION": {
        ("surface", "fibration"),
        ("surface_cover", "fibration"),
    },
    "HAS_FRAME": {("fibration", "frame_lattice")},
    "HAS_UTILITY_PROFILE": {("surface", "utility_profile")},
    "HOP_TO": {("fibration", "fibration")},
    "HAS_BASE_FIELD": {("surface", "field"), ("surface_cover", "field")},
    "HAS_BASE_FUNCTION_FIELD": {("fibration", "field")},
    "SUBFIELD_OF": {("field", "field")},
    "HAS_AMBIENT_FIELD": {("character", "field")},
    "PULLS_BACK_TO": {("character", "character")},
    "DESCENDS_TO_BASE": {("character", "fibration")},
    "SAME_SURFACE_CHARACTER": {("character", "character")},
    "SUPPORTS_CARRIER": {("fibration", "carrier")},
    "DEFINES_CHARACTER": {("carrier", "character")},
    "BASE_CHANGES_TO": {
        ("carrier", "surface_cover"),
        ("fibration", "fibration"),
    },
    "COVERS_SURFACE": {("surface_cover", "surface")},
    "HAS_MW_BLOCK": {("fibration", "mw_block")},
    "EIGENSPACE_FOR": {("mw_block", "character")},
    "REPRESENTED_BY": {("mw_block", "divisor")},
    "PRESENTED_BY": {("divisor", "fibration")},
    "DEFINED_OVER": {("divisor", "field")},
    "SAME_DIVISOR": {("divisor", "divisor")},
}


def fail(message: str) -> None:
    raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_pointer(document: Any, pointer: str) -> Any:
    if not pointer.startswith("/"):
        fail(f"JSON pointer must begin with '/': {pointer!r}")
    current = document
    for raw_token in pointer[1:].split("/"):
        token = raw_token.replace("~1", "/").replace("~0", "~")
        if isinstance(current, list):
            try:
                current = current[int(token)]
            except (ValueError, IndexError) as error:
                raise ValueError(f"invalid list token {token!r} in {pointer!r}") from error
        elif isinstance(current, dict):
            if token not in current:
                fail(f"missing object token {token!r} in {pointer!r}")
            current = current[token]
        else:
            fail(f"pointer {pointer!r} descends through a scalar")
    return current


def require_keys(record: dict[str, Any], keys: set[str], context: str) -> None:
    missing = keys - set(record)
    if missing:
        fail(f"{context} is missing keys {sorted(missing)}")


def validate_graph(graph_path: Path) -> dict[str, int]:
    graph = json.loads(graph_path.read_text())
    require_keys(
        graph,
        {"schema", "graph_id", "status", "sources", "nodes", "edges", "claim_boundary"},
        "graph",
    )
    if graph["schema"] != "elkies-k3.surface-fibration-character-graph.v1":
        fail(f"unexpected schema {graph['schema']!r}")

    source_records: dict[str, dict[str, Any]] = {}
    source_documents: dict[str, Any] = {}
    for source in graph["sources"]:
        require_keys(source, {"id", "path", "sha256"}, "source")
        source_id = source["id"]
        if source_id in source_records:
            fail(f"duplicate source id {source_id!r}")
        relative = Path(source["path"])
        if relative.is_absolute() or ".." in relative.parts:
            fail(f"source path must be repository-relative: {relative}")
        absolute = ROOT / relative
        if not absolute.is_file():
            fail(f"missing evidence source {relative}")
        actual_hash = sha256(absolute)
        if actual_hash != source["sha256"]:
            fail(
                f"evidence hash mismatch for {relative}: "
                f"expected {source['sha256']}, got {actual_hash}"
            )
        source_records[source_id] = source
        if absolute.suffix == ".json":
            source_documents[source_id] = json.loads(absolute.read_text())

    nodes: dict[str, dict[str, Any]] = {}
    for node in graph["nodes"]:
        require_keys(node, {"id", "kind", "label", "status", "data", "evidence"}, "node")
        node_id = node["id"]
        if node_id in nodes:
            fail(f"duplicate node id {node_id!r}")
        if node["kind"] not in KINDS:
            fail(f"unknown node kind {node['kind']!r} for {node_id}")
        if node["status"] not in STATUSES:
            fail(f"unknown node status {node['status']!r} for {node_id}")
        if node["status"] in CERTIFIED and not node["evidence"]:
            fail(f"certified node {node_id} has no evidence")
        nodes[node_id] = node

    edges: dict[str, dict[str, Any]] = {}
    for edge in graph["edges"]:
        require_keys(
            edge,
            {"id", "relation", "source", "target", "status", "data", "evidence"},
            "edge",
        )
        edge_id = edge["id"]
        if edge_id in edges:
            fail(f"duplicate edge id {edge_id!r}")
        if edge["source"] not in nodes or edge["target"] not in nodes:
            fail(f"edge {edge_id} has a missing endpoint")
        if edge["relation"] not in ALLOWED_ENDPOINTS:
            fail(f"unknown relation {edge['relation']!r} for {edge_id}")
        endpoint_kinds = (nodes[edge["source"]]["kind"], nodes[edge["target"]]["kind"])
        if endpoint_kinds not in ALLOWED_ENDPOINTS[edge["relation"]]:
            fail(
                f"edge {edge_id} has endpoint kinds {endpoint_kinds}, "
                f"not allowed for {edge['relation']}"
            )
        if edge["status"] not in STATUSES:
            fail(f"unknown edge status {edge['status']!r} for {edge_id}")
        if edge["status"] in CERTIFIED and not edge["evidence"]:
            fail(f"certified edge {edge_id} has no evidence")
        edges[edge_id] = edge

    for owner in [*nodes.values(), *edges.values()]:
        for evidence in owner["evidence"]:
            require_keys(evidence, {"source_id"}, f"evidence on {owner['id']}")
            source_id = evidence["source_id"]
            if source_id not in source_records:
                fail(f"{owner['id']} cites unknown source {source_id!r}")
            assertions = evidence.get("assertions", [])
            if assertions and source_id not in source_documents:
                fail(f"{owner['id']} puts JSON assertions on non-JSON source {source_id}")
            for assertion in assertions:
                require_keys(assertion, {"pointer", "equals"}, f"assertion on {owner['id']}")
                actual = json_pointer(source_documents[source_id], assertion["pointer"])
                if actual != assertion["equals"]:
                    fail(
                        f"assertion failed on {owner['id']} from {source_id} at "
                        f"{assertion['pointer']}: expected {assertion['equals']!r}, got {actual!r}"
                    )

    exact_has_fibration = [
        edge
        for edge in edges.values()
        if edge["relation"] == "HAS_FIBRATION" and edge["status"] == "CERTIFIED_EXACT"
    ]
    fibration_owners: dict[str, set[str]] = {}
    for edge in exact_has_fibration:
        fibration_owners.setdefault(edge["target"], set()).add(edge["source"])

    for node in nodes.values():
        if node["kind"] != "fibration" or node["status"] != "CERTIFIED_EXACT":
            continue
        owners = fibration_owners.get(node["id"], set())
        if len(owners) != 1:
            fail(f"exact fibration {node['id']} must have exactly one exact ambient owner")

    exact_or_definitional = CERTIFIED | {"DEFINITIONAL"}
    frame_edges = [
        edge
        for edge in edges.values()
        if edge["relation"] == "HAS_FRAME" and edge["status"] in exact_or_definitional
    ]
    base_function_field_edges = [
        edge
        for edge in edges.values()
        if edge["relation"] == "HAS_BASE_FUNCTION_FIELD"
        and edge["status"] in exact_or_definitional
    ]
    for fibration_id, owners in fibration_owners.items():
        base_fields = {
            edge["target"] for edge in base_function_field_edges if edge["source"] == fibration_id
        }
        if len(base_fields) != 1:
            fail(f"exact fibration {fibration_id} must name exactly one base function field")
        owner_id = next(iter(owners))
        if nodes[owner_id]["kind"] == "surface":
            frames = {edge["target"] for edge in frame_edges if edge["source"] == fibration_id}
            if len(frames) != 1:
                fail(f"same-K3 fibration {fibration_id} must name exactly one frame")

    for edge in edges.values():
        if edge["relation"] != "HOP_TO" or edge["status"] != "CERTIFIED_EXACT":
            continue
        source_owners = fibration_owners.get(edge["source"], set())
        target_owners = fibration_owners.get(edge["target"], set())
        if len(source_owners) != 1 or source_owners != target_owners:
            fail(f"exact fibration hop {edge['id']} does not stay on one ambient surface")
        owner_id = next(iter(source_owners))
        if nodes[owner_id]["kind"] != "surface":
            fail(f"exact HOP_TO edge {edge['id']} is on a surface cover, not the fixed K3")

    exact_ambient_edges = [
        edge
        for edge in edges.values()
        if edge["relation"] == "HAS_AMBIENT_FIELD"
        and edge["status"] in exact_or_definitional
    ]
    character_fields: dict[str, set[str]] = {}
    for edge in exact_ambient_edges:
        character_fields.setdefault(edge["source"], set()).add(edge["target"])
    for node in nodes.values():
        if node["kind"] == "character" and node["status"] != "UNKNOWN":
            fields = character_fields.get(node["id"], set())
            if len(fields) != 1:
                fail(f"character {node['id']} must name exactly one ambient field")

    subfield_edges = {
        (edge["source"], edge["target"])
        for edge in edges.values()
        if edge["relation"] == "SUBFIELD_OF" and edge["status"] in exact_or_definitional
    }
    for edge in edges.values():
        if edge["relation"] != "PULLS_BACK_TO" or edge["status"] != "CERTIFIED_EXACT":
            continue
        source_fields = character_fields.get(edge["source"], set())
        target_fields = character_fields.get(edge["target"], set())
        if len(source_fields) != 1 or len(target_fields) != 1:
            fail(f"exact character pullback {edge['id']} lacks unique ambient fields")
        if (next(iter(source_fields)), next(iter(target_fields))) not in subfield_edges:
            fail(f"exact character pullback {edge['id']} lacks its ambient field embedding")

    for edge in edges.values():
        if edge["status"] != "CERTIFIED_EXACT":
            continue
        if edge["relation"] == "DESCENDS_TO_BASE":
            base_character_id = edge["data"].get("base_character_id")
            if base_character_id not in nodes or nodes[base_character_id]["kind"] != "character":
                fail(f"exact descent {edge['id']} must name data.base_character_id")
            if not any(
                candidate["relation"] == "PULLS_BACK_TO"
                and candidate["source"] == base_character_id
                and candidate["target"] == edge["source"]
                and candidate["status"] == "CERTIFIED_EXACT"
                for candidate in edges.values()
            ):
                fail(f"exact descent {edge['id']} lacks the matching pullback identity")
            target_base_fields = {
                candidate["target"]
                for candidate in base_function_field_edges
                if candidate["source"] == edge["target"]
            }
            if character_fields.get(base_character_id, set()) != target_base_fields:
                fail(f"exact descent {edge['id']} uses a character on the wrong base field")
        if edge["relation"] == "SAME_SURFACE_CHARACTER":
            comparison_field = edge["data"].get("comparison_field_id")
            if comparison_field not in nodes or nodes[comparison_field]["kind"] != "field":
                fail(f"exact character identity {edge['id']} lacks a comparison field")
        if edge["relation"] == "SAME_DIVISOR":
            common_lattice = edge["data"].get("common_lattice_id")
            if common_lattice not in nodes or nodes[common_lattice]["kind"] not in {
                "neron_severi_lattice",
                "frame_lattice",
            }:
                fail(f"exact divisor identity {edge['id']} lacks a common marked lattice")

    utility_nodes = [node for node in nodes.values() if node["kind"] == "utility_profile"]
    for utility in utility_nodes:
        owners = {
            edge["source"]
            for edge in edges.values()
            if edge["relation"] == "HAS_UTILITY_PROFILE"
            and edge["target"] == utility["id"]
            and edge["status"] in CERTIFIED
        }
        if len(owners) != 1:
            fail(f"utility profile {utility['id']} must belong to exactly one surface")
        surface_id = next(iter(owners))
        surface_fibrations = {
            edge["target"]
            for edge in exact_has_fibration
            if edge["source"] == surface_id
        }
        rootless = {
            fibration_id
            for fibration_id in surface_fibrations
            if nodes[fibration_id]["data"].get("root_rank") == 0
        }
        rootless_classes = {
            nodes[fibration_id]["data"].get("frame_class") for fibration_id in rootless
        }
        hop_pairs = {
            frozenset((edge["source"], edge["target"]))
            for edge in edges.values()
            if edge["relation"] == "HOP_TO"
            and edge["status"] == "CERTIFIED_EXACT"
            and edge["source"] in surface_fibrations
            and edge["target"] in surface_fibrations
        }
        surface_carriers = {
            edge["target"]
            for edge in edges.values()
            if edge["relation"] == "SUPPORTS_CARRIER"
            and edge["status"] == "CERTIFIED_EXACT"
            and edge["source"] in surface_fibrations
        }
        exact_descents = sum(
            edge["relation"] == "DESCENDS_TO_BASE"
            and edge["status"] == "CERTIFIED_EXACT"
            and edge["target"] in surface_fibrations
            for edge in edges.values()
        )
        expected = {
            "certified_same_surface_fibrations_lower_bound": len(surface_fibrations),
            "certified_rootless_fibrations_lower_bound": len(rootless),
            "certified_rootless_frame_isometry_classes_lower_bound": len(rootless_classes),
            "certified_undirected_hop_edges_lower_bound": len(hop_pairs),
            "certified_surface_carriers_lower_bound": len(surface_carriers),
            "certified_cross_fibration_character_descents": exact_descents,
        }
        for key, value in expected.items():
            if utility["data"].get(key) != value:
                fail(
                    f"utility profile {utility['id']} has {key}={utility['data'].get(key)!r}; "
                    f"graph implies {value}"
                )

    return {
        "sources": len(source_records),
        "nodes": len(nodes),
        "edges": len(edges),
        "exact_hops": sum(
            edge["relation"] == "HOP_TO" and edge["status"] == "CERTIFIED_EXACT"
            for edge in edges.values()
        ),
        "unknown_descents": sum(
            edge["relation"] == "DESCENDS_TO_BASE" and edge["status"] == "UNKNOWN"
            for edge in edges.values()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "graphs",
        nargs="*",
        type=Path,
        default=[DEFAULT_GRAPH],
        help="graph JSON files; defaults to the determinant-948 seed",
    )
    arguments = parser.parse_args()
    for supplied_path in arguments.graphs:
        graph_path = supplied_path if supplied_path.is_absolute() else ROOT / supplied_path
        counts = validate_graph(graph_path)
        relative = graph_path.relative_to(ROOT) if graph_path.is_relative_to(ROOT) else graph_path
        print(
            f"PASS {relative}: {counts['sources']} sources, {counts['nodes']} nodes, "
            f"{counts['edges']} edges, {counts['exact_hops']} exact hops, "
            f"{counts['unknown_descents']} unknown character descents"
        )


if __name__ == "__main__":
    main()
