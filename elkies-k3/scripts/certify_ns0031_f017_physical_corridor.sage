#!/usr/bin/env sage-python
"""Replay the selected physical NS0031-S001 to NS0031-F017 corridor."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, gcd, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
ROUTE = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-ns0031-same-ns-compiler-routes-f017-witness-v1.json"
)
SOURCES = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json"
)
FOUNDRY = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
PILOT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0031-rootless-mw17-multisection-pilot-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-ns0031-a1-2a7-to-f017-physical-corridor-v1.json"
)
U = matrix(ZZ, ((0, 1), (1, 0)))

route_source = ROOT / "elkies-k3/scripts/search_lattice_foundry_same_ns_compiler_routes.sage"
route_engine = {"__file__": str(route_source), "__name__": "corridor_route_engine"}
exec(compile(route_source.read_text(), str(route_source), "exec"), route_engine)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def root_signature(frame):
    minimum = pari(frame).qfminim(2)
    signed_roots = int(minimum[0])
    if signed_roots == 0:
        return (0, 0, 1)
    roots = matrix(ZZ, minimum[2].sage()).transpose().rows()
    roots = tuple(roots) + tuple(-root for root in roots)
    basis = matrix(ZZ, roots).row_module(ZZ).basis_matrix()
    return (
        basis.nrows(),
        signed_roots,
        abs(int((basis * frame * basis.transpose()).det())),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--route", type=Path, default=ROUTE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    route_path = arguments.route.resolve()
    route_payload = json.loads(route_path.read_text())
    sources = json.loads(SOURCES.read_text())
    foundry = json.loads(FOUNDRY.read_text())
    pilot = json.loads(PILOT.read_text())
    result = next(row for row in route_payload["results"] if row["case"] == "ns0031")
    selected = result["best_routes_by_target"]["NS0031-F017"]
    source = next(
        row["source"] for row in sources["sources"] if row["source_id"] == "NS0031-S001"
    )
    ns = next(row for row in foundry["ns_classes"] if row["ns_id"] == "NS0031")
    target = next(row for row in ns["frames"] if row["frame_id"] == "NS0031-F017")
    richness = next(row for row in pilot["targets"] if row["frame_id"] == "NS0031-F017")

    current_frame = matrix(ZZ, source["root_adapted_gram"])
    source_ns = block_diagonal_matrix(U, -current_frame)
    current_ns = source_ns
    composed = identity_matrix(ZZ, 19)
    chain = [source["root_type"]]
    replayed = []
    for edge_index, edge in enumerate(selected["edges"], start=1):
        assert edge["edge_index"] == edge_index
        fibre = vector(ZZ, edge["fibre"])
        assert fibre * current_ns * fibre == 0
        assert gcd(tuple(current_ns * fibre)) == 1
        context = route_engine["edge_context"](current_frame, edge["source_root_rank"])
        nef = route_engine["physical_nef_profile"](fibre, context)
        assert nef == edge["physical_nef_profile"] and nef["passes"]
        pole, norm, closest = route_engine["section_pole_order"](
            context, edge["horizontal_mw_quotient"]
        )
        assert pole == edge["horizontal_P_dot_O"]
        assert norm == edge["horizontal_minimum_frame_norm"]
        assert closest == edge["horizontal_closest_root_coordinates"]
        walls = route_engine["negative_horizontal_walls"](fibre, current_frame)
        assert walls == edge["negative_horizontal_walls"] == []

        if edge["physical_weyl_repairs"]:
            assert edge_index == len(selected["edges"])
            assert edge["physical_weyl_repairs"] == 1
            assert edge["physical_weyl_repair_sequence"] == ["A1 finite simple root"]
            simple_root = vector(ZZ, [0] * 19)
            simple_root[2] = -1
            reflection = identity_matrix(ZZ, 19) + (
                current_ns * simple_root.column()
            ) * simple_root.row()
            assert vector(ZZ, vector(ZZ, edge["pre_repair_fibre"]) * reflection) == fibre
        else:
            assert edge["physical_weyl_repair_sequence"] == []

        transport = matrix(ZZ, edge["edge_transport_child_to_parent"])
        assert abs(transport.det()) == 1
        next_ns = transport * current_ns * transport.transpose()
        assert next_ns[:2, :2] == U
        assert next_ns[:2, 2:] == 0 and next_ns[2:, :2] == 0
        next_frame = -next_ns[2:, 2:]
        assert next_frame.is_positive_definite() and abs(next_frame.det()) == 1184
        assert root_signature(next_frame)[0] == edge["target_root_rank"]
        composed = transport * composed
        assert composed * source_ns * composed.transpose() == next_ns
        chain.append(edge["target_root_type"] or "rootless")
        replayed.append(
            {
                "edge_index": edge_index,
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
    terminal_isometry = matrix(ZZ, selected["terminal_to_target_frame_isometry"])
    target_frame = matrix(ZZ, target["gram"])
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
    assert cost == selected["cost"] == [6, 8, 2, 5, 1]

    payload = {
        "schema": "elkies-k3.ns0031-f017-physical-corridor.v1",
        "status": "PASS_EXACT_MARKING_LEVEL_PHYSICAL_CORRIDOR_TO_F017",
        "inputs": {
            relative(route_path): digest(route_path),
            relative(SOURCES): digest(SOURCES),
            relative(FOUNDRY): digest(FOUNDRY),
            relative(PILOT): digest(PILOT),
            relative(route_source): digest(route_source),
        },
        "source": {
            "source_id": "NS0031-S001",
            "root_type": source["root_type"],
            "mw_rank": source["mw_rank_for_rho_19"],
            "interpretation": "prescribed-root marking in the determinant-1184 NS",
        },
        "target": {
            "frame_id": "NS0031-F017",
            "root_type": "rootless",
            "mw_rank": 17,
            "richness_coordinates": richness["richness_coordinates"],
        },
        "route": {
            "root_type_chain": chain,
            "cost": cost,
            "edges": replayed,
            "terminal_integrally_isometric_to_target": True,
            "composed_transport_unimodular": True,
        },
        "search_scope": {
            "forward": "bounded rank-first beam with q=4,6,8, degree two and P.O.<=8",
            "reverse": (
                "complete 63,902-class minimum-norm-eight degree-two shell on F017; "
                "six matches to the retained low-root forward frontier"
            ),
            "optimality": "selected compiler-cost leader inside the declared bounded meeting search",
        },
        "proof_boundary": {
            "proved": (
                "The selected marking-level route stays in the determinant-1184 NS. "
                "Every edge is primitive and physically nef, passes the complete all-section "
                "and finite horizontal-wall gates, and has exact unimodular transport. "
                "The terminal frame is integrally isometric to the trisection-first F017 target."
            ),
            "not_proved": (
                "This is not an equation-level corridor over QQ and is not a global route-"
                "optimality theorem. The NS0031 source has an exact GF(7) marking and finite "
                "7-adic lift, but no rational equation or formal characteristic-zero family "
                "is supplied by this certificate."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/certify_ns0031_f017_physical_corridor.sage --check"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print("NS0031F017CORRIDOR|edges=5|max_PO=6|max_q=8|repairs=1|status=PASS")


if __name__ == "__main__":
    main()
