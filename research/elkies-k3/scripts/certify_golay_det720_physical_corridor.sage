#!/usr/bin/env sage-python
"""Replay a low-pole physical corridor from G720-S0128 to G720-F001.

The input route is discovery output from the bounded same-NS compiler search.
This checker does not repeat the beam search.  It independently replays every
selected fibre, physical-nef gate, horizontal-wall gate, unimodular transport,
and the terminal integral isometry.  The result is a certificate for the
selected route, not a completeness or optimality theorem for the graph.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, gcd, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
ROUTE = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-golay-det720-same-ns-compiler-routes-pole4-witness-v1.json"
)
SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
)
TARGET = ROOT / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
REJECTION = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-golay-det720-3a5-saturation-rejection-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-golay-det720-3a5-to-mw17-physical-corridor-v1.json"
)
U = matrix(ZZ, ((0, 1), (1, 0)))

route_source = ROOT / "elkies-k3/scripts/search_lattice_foundry_same_ns_compiler_routes.sage"
route_engine = {"__file__": str(route_source), "__name__": "corridor_route_engine"}
exec(compile(route_source.read_text(), str(route_source), "exec"), route_engine)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rows(value) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in value.rows()]


def root_signature(frame) -> tuple[int, int, int]:
    minimum = pari(frame).qfminim(2)
    signed_roots = int(minimum[0])
    if signed_roots == 0:
        return (0, 0, 1)
    roots = matrix(ZZ, minimum[2].sage()).transpose().rows()
    roots = tuple(roots) + tuple(-root for root in roots)
    root_lattice = matrix(ZZ, roots).row_module(ZZ)
    basis = root_lattice.basis_matrix()
    return (basis.nrows(), signed_roots, abs(int((basis * frame * basis.transpose()).det())))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--route", type=Path, default=ROUTE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    route_path = arguments.route.resolve()
    output_path = arguments.output.resolve()
    route_payload = json.loads(route_path.read_text())
    source_payload = json.loads(SOURCES.read_text())
    target_payload = json.loads(TARGET.read_text())
    rejection_payload = json.loads(REJECTION.read_text())

    assert route_payload["status"] == "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_HIT"
    result = next(row for row in route_payload["results"] if row["case"] == "golay720")
    selected = result["best_routes_by_target"]["G720-F001"]
    assert selected["cost"][0] == 4
    assert result["source"]["source_id"] == "G720-S0128"
    source_row = next(
        row for row in source_payload["sources"] if row["source_id"] == "G720-S0128"
    )
    source = source_row["source"]
    current_frame = matrix(ZZ, source["root_adapted_gram"])
    source_ns = block_diagonal_matrix(U, -current_frame)
    current_ns = source_ns
    composed = matrix(ZZ, 19, 19, 1)
    chain = [source["root_type"]]
    replayed_edges = []

    for index, edge in enumerate(selected["edges"], start=1):
        assert edge["edge_index"] == index
        fibre = vector(ZZ, edge["fibre"])
        assert gcd(tuple(current_ns * fibre)) == 1
        context = route_engine["edge_context"](current_frame, edge["source_root_rank"])
        nef = route_engine["physical_nef_profile"](fibre, context)
        assert nef == edge["physical_nef_profile"] and nef["passes"]
        pole, section_norm, closest_root = route_engine["section_pole_order"](
            context, edge["horizontal_mw_quotient"]
        )
        assert pole == edge["horizontal_P_dot_O"] <= selected["cost"][0]
        assert section_norm == edge["horizontal_minimum_frame_norm"]
        assert closest_root == edge["horizontal_closest_root_coordinates"]
        walls = route_engine["negative_horizontal_walls"](fibre, current_frame)
        assert walls == edge["negative_horizontal_walls"] == []

        transport = matrix(ZZ, edge["edge_transport_child_to_parent"])
        assert abs(transport.det()) == 1
        next_ns = transport * current_ns * transport.transpose()
        assert next_ns[:2, :2] == U
        assert next_ns[:2, 2:] == 0 and next_ns[2:, :2] == 0
        next_frame = -next_ns[2:, 2:]
        assert next_frame.is_positive_definite() and next_frame.det() == 720
        signature = root_signature(next_frame)
        assert signature[0] == edge["target_root_rank"]
        composed = transport * composed
        assert composed * source_ns * composed.transpose() == next_ns
        chain.append(edge["target_root_type"] or "rootless")
        replayed_edges.append(
            {
                "edge_index": index,
                "source_root_rank": edge["source_root_rank"],
                "target_root_rank": edge["target_root_rank"],
                "target_root_type": edge["target_root_type"],
                "q": edge["q"],
                "old_fibre_degree": edge["old_fibre_degree"],
                "horizontal_P_dot_O": pole,
                "physical_weyl_repairs": edge["physical_weyl_repairs"],
                "primitive_fibre": True,
                "physical_nef": True,
                "all_section_gate": True,
                "negative_horizontal_walls": 0,
                "unimodular_transport": True,
            }
        )
        current_frame = next_frame
        current_ns = next_ns

    assert rows(current_frame) == selected["terminal_frame"]
    assert rows(composed) == selected["composed_transport_terminal_to_source"]
    target_frame = matrix(ZZ, target_payload["frame"]["gram"])
    terminal_isometry = matrix(ZZ, selected["terminal_to_target_frame_isometry"])
    assert abs(terminal_isometry.det()) == 1
    assert terminal_isometry * current_frame * terminal_isometry.transpose() == target_frame
    assert root_signature(current_frame) == (0, 0, 1)
    cost = [
        max(edge["horizontal_P_dot_O"] for edge in selected["edges"]),
        max(edge["q"] for edge in selected["edges"]),
        max(edge["old_fibre_degree"] for edge in selected["edges"]),
        len(selected["edges"]),
        sum(edge["physical_weyl_repairs"] for edge in selected["edges"]),
    ]
    assert cost == selected["cost"]
    assert max(edge["horizontal_P_dot_O"] for edge in selected["edges"][:-1]) <= 2
    assert selected["edges"][-1]["horizontal_P_dot_O"] == 4
    assert rejection_payload["status"] == (
        "PASS_EXACT_RATIONAL_POINT_REJECTED_NS_DET20_TORSION3_HALF_SECTION"
    )

    payload = {
        "schema": "elkies-k3.golay-det720-physical-corridor.v1",
        "status": "PASS_EXACT_MARKING_LEVEL_PHYSICAL_CORRIDOR_FINAL_PO_4",
        "inputs": {
            relative(route_path): digest(route_path),
            relative(SOURCES): digest(SOURCES),
            relative(TARGET): digest(TARGET),
            relative(REJECTION): digest(REJECTION),
            relative(route_source): digest(route_source),
        },
        "source": {
            "source_id": "G720-S0128",
            "root_type": source["root_type"],
            "mw_rank": source["mw_rank_for_rho_19"],
            "interpretation": "prescribed-root marking in the determinant-720 NS",
        },
        "target": {"frame_id": "G720-F001", "root_type": "rootless", "mw_rank": 17},
        "route": {
            "root_type_chain": chain,
            "cost": cost,
            "edges": replayed_edges,
            "terminal_integrally_isometric_to_target": True,
            "composed_transport_unimodular": True,
        },
        "proof_boundary": {
            "proved": (
                "The selected marking-level route lies in one determinant-720 NS. "
                "Every edge has primitive fibre, passes the finite/affine and all-section "
                "nef gates, has no negative finite-degree horizontal wall, and carries "
                "an exact unimodular NS transport. The first five edges have P.O. at "
                "most two; the final edge has P.O. four after two component-Weyl repairs. The terminal "
                "frame is integrally isometric to the rootless MW17 target."
            ),
            "not_proved": (
                "This is not an equation-level corridor over QQ. The stored rational "
                "3I6 specialization has full NS determinant 20 and is explicitly not "
                "the determinant-720 source marking. No rational determinant-720 start "
                "or algebraic neighbour equations are supplied here."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/certify_golay_det720_physical_corridor.sage --check"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "GOLAY720CORRIDOR|edges={}|max_PO={}|max_q={}|target=MW17|status=PASS".format(
            len(selected["edges"]), cost[0], cost[1]
        )
    )


if __name__ == "__main__":
    main()
