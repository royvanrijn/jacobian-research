#!/usr/bin/env sage-python
"""Extract the selected six-edge Golay-720 corridor witness.

This discovery helper consumes a retained compiler-search work directory.  Its
output is self-contained and is independently replayed by
``certify_golay_det720_physical_corridor.sage``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-golay-det720-same-ns-compiler-routes-pole4-witness-v1.json"
U = matrix(ZZ, ((0, 1), (1, 0)))

route_path = ROOT / "elkies-k3/scripts/search_lattice_foundry_same_ns_compiler_routes.sage"
route = {"__name__": "route", "__file__": str(route_path)}
exec(compile(route_path.read_text(), str(route_path), "exec"), route)
engine_path = ROOT / "elkies-k3/scripts/exact_neighbor_engine.sage"
engine = {"__name__": "engine", "__file__": str(engine_path)}
exec(compile(engine_path.read_text(), str(engine_path), "exec"), engine)


def load_frame(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")])


def rows(M):
    return [[int(x) for x in row] for row in M.rows()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search-workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    work = args.search_workdir.resolve()
    sources = json.loads((ROOT / "artifacts/generated-results/elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json").read_text())
    source = next(x["source"] for x in sources["sources"] if x["source_id"] == "G720-S0128")
    current = matrix(ZZ, source["root_adapted_gram"])
    edges = []
    for depth, (parent_index, child_index) in enumerate(((1, 7), (7, 16), (16, 18), (18, 29), (29, 34)), start=1):
        child = load_frame(work / f"frame-{child_index:05d}.txt")
        payload = json.loads((work / f"neighbors-{parent_index:05d}.json").read_text())
        record = next(x for x in payload["neighbors"] if x.get("child_root_adapted_frame") is not None and matrix(ZZ, x["child_root_adapted_frame"]) == child)
        transport = block_diagonal_matrix(matrix(ZZ, 2, 2, 1), matrix(ZZ, record["child_root_adapted_basis"])) * matrix(ZZ, record["neighbor_basis"])
        context = route["edge_context"](current, int(payload["input_root_data"][0]))
        fibre = vector(ZZ, record["fiber"])
        nef = route["physical_nef_profile"](fibre, context)
        pole, norm, closest = route["section_pole_order"](context, record["mw_projection"])
        walls = route["negative_horizontal_walls"](fibre, current)
        assert nef["passes"] and pole <= 2 and not walls
        edges.append({"edge_index": depth, "q": int(record["q"]), "factor_order": list(map(int, record["factor_order"])), "old_fibre_degree": int(record["old_fiber_degree"]), "orbit_index": int(record["orbit_index"]), "horizontal_P_dot_O": pole, "horizontal_minimum_frame_norm": norm, "horizontal_mw_quotient": list(map(int, record["mw_projection"])), "horizontal_closest_root_coordinates": closest, "physical_weyl_repairs": 0, "source_root_type": source["root_type"] if depth == 1 else edges[-1]["target_root_type"], "source_root_rank": int(payload["input_root_data"][0]), "source_mw_rank": 17-int(payload["input_root_data"][0]), "target_root_type": record["child_ade"], "target_root_rank": int(record["child_root_data"][0]), "target_mw_rank": int(record["child_mw_rank"]), "fibre": list(map(int, fibre)), "physical_nef_profile": nef, "negative_horizontal_walls": walls, "edge_transport_child_to_parent": rows(transport)})
        current = child

    direct = json.loads((ROOT / "artifacts/generated-results/elkies-k3-golay-det720-degree2-direct-3a5-corridor-v1.json").read_text())
    reduced = matrix(ZZ, direct["target"]["reduced_target_frame"])
    reduction = matrix(ZZ, direct["target"]["target_basis_to_reduced_basis"])
    target_ns = block_diagonal_matrix(U, -reduced)
    target_fibre = vector(ZZ, [2,2,2,1,2,1,0,2,0,-1,1,1,-1,0,0,0,0,0,0])
    split = engine["primitive_hyperbolic_split"](target_ns, target_fibre)
    iso = route["integral_isometry"](split["child_frame"], current)
    assert iso is not None
    C = block_diagonal_matrix(matrix(ZZ, 2, 2, 1), iso) * split["transport"]
    current_ns = block_diagonal_matrix(U, -current)
    D = vector(ZZ, vector(ZZ, [1,0]+[0]*17) * C.inverse())
    W = identity_matrix(ZZ, 19)
    for name in ("a0", "s1"):
        r = vector(ZZ, [0]*19)
        if name == "a0": r[0], r[2] = 1, 1
        else: r[3] = -1
        reflection = identity_matrix(ZZ, 19) + (current_ns * r.column()) * r.row()
        D = D * reflection
        W = W * reflection
    transport = C.inverse() * W
    assert vector(ZZ, transport[0]) == D
    assert transport * current_ns * transport.transpose() == target_ns
    context = route["edge_context"](current, 2)
    nef = route["physical_nef_profile"](D, context)
    pole, norm, closest = route["section_pole_order"](context, list(D[4:]))
    walls = route["negative_horizontal_walls"](D, current)
    assert nef["passes"] and pole == 4 and not walls
    edges.append({"edge_index": 6, "q": 4, "factor_order": [2,2], "old_fibre_degree": 2, "orbit_index": -1, "horizontal_P_dot_O": pole, "horizontal_minimum_frame_norm": norm, "horizontal_mw_quotient": list(map(int,D[4:])), "horizontal_closest_root_coordinates": closest, "physical_weyl_repairs": 2, "source_root_type": edges[-1]["target_root_type"], "source_root_rank": 2, "source_mw_rank": 15, "target_root_type": "", "target_root_rank": 0, "target_mw_rank": 17, "fibre": list(map(int,D)), "physical_nef_profile": nef, "negative_horizontal_walls": walls, "edge_transport_child_to_parent": rows(transport)})
    composed = identity_matrix(ZZ,19)
    ns = block_diagonal_matrix(U, -matrix(ZZ,source["root_adapted_gram"]))
    for edge in edges: composed = matrix(ZZ,edge["edge_transport_child_to_parent"])*composed
    selected = {"target_frame_id":"G720-F001", "cost":[4,4,2,6,2], "edges":edges, "terminal_frame":rows(reduced), "terminal_to_target_frame_isometry":rows(reduction.inverse()), "composed_transport_terminal_to_source":rows(composed)}
    output = {"schema":"elkies-k3.same-ns-compiler-route-search.v1", "status":"PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_HIT", "results":[{"case":"golay720", "source":{"source_id":"G720-S0128","root_type":source["root_type"],"root_rank":15,"mw_rank":2,"frame_gram_sha256":source["gram_sha256"]}, "requested_targets":["G720-F001"], "search":{"witness_extraction":True}, "accounting":[], "best_routes_by_target":{"G720-F001":selected}, "inputs":{}}], "proof_boundary":{"proved":"Selected exact witness only.","not_proved":"No optimality or equation lift."}, "reproduce":"See build_golay_det720_physical_corridor_witness.sage."}
    args.output.resolve().write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print("GOLAY720CORRIDORWITNESS|edges=6|max_PO=4|status=PASS")


if __name__ == "__main__": main()
