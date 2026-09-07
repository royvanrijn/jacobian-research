#!/usr/bin/env sage-python
"""Extract the selected NS0031 source-to-F017 physical corridor.

The first four edges are selected from a retained bounded rank-first search
work directory.  The last edge is reconstructed independently from the
complete minimum-norm degree-two shell of the preferred rootless frame F017.
The resulting JSON is self-contained and is replayed by
``certify_ns0031_f017_physical_corridor.sage``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json"
)
FOUNDRY = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-ns0031-same-ns-compiler-routes-f017-witness-v1.json"
)
SOURCE_ID = "NS0031-S001"
TARGET_ID = "NS0031-F017"
U = matrix(ZZ, ((0, 1), (1, 0)))

route_path = ROOT / "elkies-k3/scripts/search_lattice_foundry_same_ns_compiler_routes.sage"
route = {"__name__": "route", "__file__": str(route_path)}
exec(compile(route_path.read_text(), str(route_path), "exec"), route)
engine_path = ROOT / "elkies-k3/scripts/exact_neighbor_engine.sage"
engine = {"__name__": "engine", "__file__": str(engine_path)}
exec(compile(engine_path.read_text(), str(engine_path), "exec"), engine)


def load_frame(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search-workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    work = arguments.search_workdir.resolve()
    sources = json.loads(SOURCES.read_text())
    source = next(
        row["source"] for row in sources["sources"] if row["source_id"] == SOURCE_ID
    )
    foundry = json.loads(FOUNDRY.read_text())
    ns = next(row for row in foundry["ns_classes"] if row["ns_id"] == "NS0031")
    target = next(row for row in ns["frames"] if row["frame_id"] == TARGET_ID)

    current = matrix(ZZ, source["root_adapted_gram"])
    source_ns = block_diagonal_matrix(U, -current)
    composed = identity_matrix(ZZ, 19)
    edges = []
    # Expansion indices in the retained rank-first beam.  This is a four-edge
    # route; graph distance was not the beam objective.
    selected = ((1, 3), (3, 19), (19, 22), (22, 31))
    for edge_index, (parent_index, child_index) in enumerate(selected, start=1):
        child = load_frame(work / f"frame-{child_index:05d}.txt")
        payload = json.loads((work / f"neighbors-{parent_index:05d}.json").read_text())
        record = next(
            row
            for row in payload["neighbors"]
            if row.get("child_root_adapted_frame") is not None
            and matrix(ZZ, row["child_root_adapted_frame"]) == child
        )
        adaptation = matrix(ZZ, record["child_root_adapted_basis"])
        neighbor_basis = matrix(ZZ, record["neighbor_basis"])
        transport = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * neighbor_basis
        context = route["edge_context"](current, int(payload["input_root_data"][0]))
        fibre = vector(ZZ, record["fiber"])
        nef = route["physical_nef_profile"](fibre, context)
        pole, norm, closest = route["section_pole_order"](
            context, record["mw_projection"]
        )
        walls = route["negative_horizontal_walls"](fibre, current)
        assert nef["passes"] and not walls
        current_ns = block_diagonal_matrix(U, -current)
        child_ns = block_diagonal_matrix(U, -child)
        assert abs(transport.det()) == 1
        assert transport * current_ns * transport.transpose() == child_ns
        composed = transport * composed
        assert composed * source_ns * composed.transpose() == child_ns
        edges.append(
            {
                "edge_index": edge_index,
                "q": int(record["q"]),
                "factor_order": list(map(int, record["factor_order"])),
                "old_fibre_degree": int(record["old_fiber_degree"]),
                "orbit_index": int(record["orbit_index"]),
                "horizontal_P_dot_O": pole,
                "horizontal_minimum_frame_norm": norm,
                "horizontal_mw_quotient": list(map(int, record["mw_projection"])),
                "horizontal_closest_root_coordinates": closest,
                "physical_weyl_repairs": 0,
                "physical_weyl_repair_sequence": [],
                "source_root_type": payload["input_ade"],
                "source_root_rank": int(payload["input_root_data"][0]),
                "source_mw_rank": int(payload["input_mw_rank"]),
                "target_root_type": record["child_ade"],
                "target_root_rank": int(record["child_root_data"][0]),
                "target_mw_rank": int(record["child_mw_rank"]),
                "fibre": list(map(int, fibre)),
                "physical_nef_profile": nef,
                "negative_horizontal_walls": walls,
                "edge_transport_child_to_parent": rows(transport),
            }
        )
        current = child

    assert edges[-1]["target_root_type"] == "A1"
    target_gram = matrix(ZZ, target["gram"])
    reduction = target_gram.LLL_gram().transpose()
    reduced_target = reduction * target_gram * reduction.transpose()
    assert abs(reduction.det()) == 1
    target_ns = block_diagonal_matrix(U, -reduced_target)
    target_w = vector(ZZ, [1, 1, 1, -1, -1, -1, 1, -1, 1, -1, 0, -1, 0, 1, -1, 1, -1])
    target_fibre = vector(ZZ, [2, 2] + list(target_w))
    split = engine["primitive_hyperbolic_split"](target_ns, target_fibre)
    isometry = route["integral_isometry"](split["child_frame"], current)
    assert isometry is not None
    bridge = block_diagonal_matrix(identity_matrix(ZZ, 2), isometry) * split["transport"]
    current_ns = block_diagonal_matrix(U, -current)
    assert bridge * target_ns * bridge.transpose() == current_ns
    raw_fibre = vector(ZZ, vector(ZZ, [1, 0] + [0] * 17) * bridge.inverse())

    simple_root = vector(ZZ, [0] * 19)
    simple_root[2] = -1
    reflection = identity_matrix(ZZ, 19) + (current_ns * simple_root.column()) * simple_root.row()
    repaired_fibre = vector(ZZ, raw_fibre * reflection)
    transport = bridge.inverse() * reflection
    assert vector(ZZ, transport[0]) == repaired_fibre
    assert abs(transport.det()) == 1
    assert transport * current_ns * transport.transpose() == target_ns
    context = route["edge_context"](current, 1)
    nef = route["physical_nef_profile"](repaired_fibre, context)
    pole, norm, closest = route["section_pole_order"](context, list(repaired_fibre[3:]))
    walls = route["negative_horizontal_walls"](repaired_fibre, current)
    assert nef["passes"] and pole == 6 and not walls
    composed = transport * composed
    assert composed * source_ns * composed.transpose() == target_ns
    edges.append(
        {
            "edge_index": 5,
            "q": int(repaired_fibre[0] * repaired_fibre[1]),
            "factor_order": [int(repaired_fibre[0]), int(repaired_fibre[1])],
            "old_fibre_degree": int(repaired_fibre[1]),
            "orbit_index": -1,
            "reverse_target_shell": {
                "target_frame_id": TARGET_ID,
                "target_side_q": 4,
                "target_side_minimum_norm": 8,
                "target_side_shell_classes": 63902,
                "target_side_matches_to_retained_low_root_frontier": 6,
                "selected_shell_index": 33100,
                "target_side_representative_in_reduced_basis": list(map(int, target_w)),
            },
            "horizontal_P_dot_O": pole,
            "horizontal_minimum_frame_norm": norm,
            "horizontal_mw_quotient": list(map(int, repaired_fibre[3:])),
            "horizontal_closest_root_coordinates": closest,
            "physical_weyl_repairs": 1,
            "physical_weyl_repair_sequence": ["A1 finite simple root"],
            "pre_repair_fibre": list(map(int, raw_fibre)),
            "source_root_type": "A1",
            "source_root_rank": 1,
            "source_mw_rank": 16,
            "target_root_type": "",
            "target_root_rank": 0,
            "target_mw_rank": 17,
            "fibre": list(map(int, repaired_fibre)),
            "physical_nef_profile": nef,
            "negative_horizontal_walls": walls,
            "edge_transport_child_to_parent": rows(transport),
        }
    )

    cost = [
        max(edge["horizontal_P_dot_O"] for edge in edges),
        max(edge["q"] for edge in edges),
        max(edge["old_fibre_degree"] for edge in edges),
        len(edges),
        sum(edge["physical_weyl_repairs"] for edge in edges),
    ]
    assert cost == [6, 8, 2, 5, 1]
    route_record = {
        "target_frame_id": TARGET_ID,
        "cost": cost,
        "edges": edges,
        "terminal_frame": rows(reduced_target),
        "terminal_to_target_frame_isometry": rows(reduction.inverse().change_ring(ZZ)),
        "composed_transport_terminal_to_source": rows(composed),
    }
    payload = {
        "schema": "elkies-k3.same-ns-compiler-route-search.v1",
        "status": "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_HIT",
        "results": [
            {
                "case": "ns0031",
                "source": {
                    "source_id": SOURCE_ID,
                    "root_type": source["root_type"],
                    "root_rank": int(source["root_rank"]),
                    "mw_rank": int(source["mw_rank_for_rho_19"]),
                    "frame_gram_sha256": source["gram_sha256"],
                },
                "requested_targets": [TARGET_ID],
                "search": {
                    "forward_beam": "bounded rank-first q=4,6,8, degree=2, P.O.<=8",
                    "reverse_target_shell": "complete F017 minimum-norm-eight degree-two shell",
                    "cost_key_order": [
                        "maximum_horizontal_P_dot_O",
                        "maximum_q",
                        "maximum_old_fibre_degree",
                        "edge_count",
                        "total_physical_weyl_repairs",
                    ],
                },
                "best_routes_by_target": {TARGET_ID: route_record},
            }
        ],
        "proof_boundary": {
            "proved": "Selected exact marking-level physical corridor with unimodular transport.",
            "not_proved": "No global route optimality, equation lift, or QQ descent.",
        },
        "reproduce": "Replay with certify_ns0031_f017_physical_corridor.sage --check.",
    }
    arguments.output.resolve().write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("NS0031F017WITNESS|edges=5|max_PO=6|max_q=8|repairs=1|status=PASS")


if __name__ == "__main__":
    main()
