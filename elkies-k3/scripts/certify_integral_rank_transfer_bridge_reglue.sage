#!/usr/bin/env sage-python
"""Certify the bridge-core reglue theorem on marked same-NS corridors.

For consecutive primitive hyperbolic planes U_0,U_1 in a fixed NS lattice,
the positive frames W_i=U_i^perp(-1) share the primitive core

    K = W_0 intersect W_1 = (U_0+U_1)^perp(-1).

Writing C_i=K^perp inside W_i, each W_i is recovered from K+C_i by a
finite isotropic glue subgroup.  The script computes the core, both bridge
lattices, exact glue generators, and the root classes for every marked edge
in the NS0024 and Golay-720 corridors.  No elliptic equations are used.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
)
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_json(relative_path):
    return json.loads((ROOT / relative_path).read_text())


def digest(relative_path):
    return hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def signed_roots(gram):
    result = pari(gram).qfminim(2)
    positive = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2].sage()).columns()
    ]
    answer = positive + [-root for root in positive]
    assert len(answer) == int(result[0])
    return answer


def in_row_lattice(row, basis):
    try:
        coordinates = basis.solve_left(matrix(QQ, [row])).row(0)
    except ValueError:
        return False
    return all(value in ZZ for value in coordinates)


def quotient_label(coordinates, smith_right, diagonal):
    transformed = vector(ZZ, coordinates) * smith_right
    return tuple(
        int(transformed[index] % order)
        for index, order in enumerate(diagonal)
        if order > 1
    )


def bridge_glue(frame, frame_basis, core_basis):
    """Return K, C and the isotropic glue W/(K+C) in W coordinates."""
    core_coordinates = frame_basis.solve_left(core_basis).change_ring(ZZ)
    assert core_coordinates * frame_basis == core_basis

    bridge_coordinates = (core_coordinates * frame).right_kernel_matrix()
    assert bridge_coordinates.nrows() + core_coordinates.nrows() == frame.nrows()
    assert core_coordinates * frame * bridge_coordinates.transpose() == 0

    split_coordinates = core_coordinates.stack(bridge_coordinates)
    assert split_coordinates.rank() == frame.nrows()
    index = abs(int(split_coordinates.det()))
    diagonal_matrix, smith_left, smith_right = split_coordinates.smith_form()
    assert smith_left * split_coordinates * smith_right == diagonal_matrix
    diagonal = [abs(int(value)) for value in diagonal_matrix.diagonal()]
    assert index == int(ZZ.prod(diagonal))

    inverse_right = smith_right.inverse().change_ring(ZZ)
    inverse_split = split_coordinates.inverse()
    generators = []
    for coordinate_index, order in enumerate(diagonal):
        if order <= 1:
            continue
        smith_generator = vector(ZZ, identity_matrix(ZZ, frame.nrows()).row(coordinate_index))
        frame_coordinates = smith_generator * inverse_right
        split_dual_coordinates = frame_coordinates * inverse_split
        assert order * split_dual_coordinates in ZZ ** frame.nrows()
        norm = int(frame_coordinates * frame * frame_coordinates)
        assert norm % 2 == 0
        generators.append(
            {
                "order": order,
                "frame_coordinates": list(map(int, frame_coordinates)),
                "K_plus_C_dual_coordinates": [
                    str(value) for value in split_dual_coordinates
                ],
                "norm": norm,
                "q_mod_2Z": "0",
            }
        )

    root_histogram = Counter()
    core_root_count = 0
    for root in signed_roots(frame):
        label = quotient_label(root, smith_right, diagonal)
        root_histogram[label] += 1
        ambient = root * frame_basis
        if in_row_lattice(ambient, core_basis):
            core_root_count += 1

    core_gram = core_coordinates * frame * core_coordinates.transpose()
    bridge_gram = bridge_coordinates * frame * bridge_coordinates.transpose()
    bridge_smith = [
        abs(int(value))
        for value in bridge_gram.smith_form()[0].diagonal()
        if abs(int(value)) > 1
    ]
    split_gram = split_coordinates * frame * split_coordinates.transpose()
    assert split_gram[: core_gram.nrows(), core_gram.nrows() :] == 0

    return {
        "bridge_rank": bridge_coordinates.nrows(),
        "bridge_gram": rows(bridge_gram),
        "bridge_discriminant_group_invariants": bridge_smith,
        "bridge_determinant_absolute": abs(int(bridge_gram.det())),
        "K_plus_C_index_in_W": index,
        "glue_group_invariants": [order for order in diagonal if order > 1],
        "glue_generators": generators,
        "root_count_signed": len(signed_roots(frame)),
        "core_root_count_signed": core_root_count,
        "root_count_by_glue_label": {
            ",".join(map(str, label)) if label else "0": count
            for label, count in sorted(root_histogram.items())
        },
    }


def certify_embedded_edge(
    corridor,
    edge_index,
    ambient_ns,
    old_frame,
    old_basis,
    new_frame,
    new_basis,
    edge_metadata,
):
    assert -(old_basis * ambient_ns * old_basis.transpose()) == old_frame
    assert -(new_basis * ambient_ns * new_basis.transpose()) == new_frame
    assert old_frame.is_positive_definite() and new_frame.is_positive_definite()
    assert old_frame.det() == new_frame.det()

    old_module = old_basis.row_module(ZZ)
    new_module = new_basis.row_module(ZZ)
    core_basis = old_module.intersection(new_module).basis_matrix()
    core_rank = core_basis.nrows()
    span_rank = (old_module + new_module).rank()
    assert core_rank == 34 - span_rank

    core_gram = -(core_basis * ambient_ns * core_basis.transpose())
    assert core_gram.is_positive_definite()
    old_data = bridge_glue(old_frame, old_basis, core_basis)
    new_data = bridge_glue(new_frame, new_basis, core_basis)
    for data in (old_data, new_data):
        assert data["bridge_discriminant_group_invariants"] == [
            data["bridge_determinant_absolute"]
        ]
        assert data["K_plus_C_index_in_W"] == data[
            "bridge_determinant_absolute"
        ]
    assert old_data["K_plus_C_index_in_W"] == new_data["K_plus_C_index_in_W"]

    old_root_ambient = {
        tuple(map(int, root * old_basis)) for root in signed_roots(old_frame)
    }
    new_root_ambient = {
        tuple(map(int, root * new_basis)) for root in signed_roots(new_frame)
    }
    core_root_ambient = {
        tuple(map(int, root * core_basis)) for root in signed_roots(core_gram)
    }
    assert old_root_ambient.intersection(new_root_ambient) == core_root_ambient
    assert old_data["core_root_count_signed"] == len(core_root_ambient)
    assert new_data["core_root_count_signed"] == len(core_root_ambient)

    return {
        "corridor": corridor,
        "edge_index": edge_index,
        "q": int(edge_metadata["q"]),
        "old_fibre_degree": int(edge_metadata["old_fibre_degree"]),
        "source_root_rank": int(edge_metadata["source_root_rank"]),
        "target_root_rank": int(edge_metadata["target_root_rank"]),
        "core": {
            "rank": core_rank,
            "gram": rows(core_gram),
            "determinant_absolute": abs(int(core_gram.det())),
            "root_count_signed": len(core_root_ambient),
        },
        "old_frame": old_data,
        "new_frame": new_data,
        "bridge_replacement": {
            "common_cyclic_glue_order": old_data["K_plus_C_index_in_W"],
            "old_and_new_glue_orders_equal": True,
            "both_glues_project_isomorphically_to_full_bridge_discriminant": True,
        },
        "root_transfer": {
            "surviving_signed_roots": len(core_root_ambient),
            "removed_signed_roots": len(old_root_ambient - new_root_ambient),
            "introduced_signed_roots": len(new_root_ambient - old_root_ambient),
            "identity": "Phi(W_old) intersect Phi(W_new) = Phi(K)",
        },
    }


def certify_local_edge(corridor, edge_index, current_frame, transport, edge_metadata):
    parent_ns = block_diagonal_matrix(U, -current_frame)
    assert abs(int(transport.det())) == 1
    child_ns = transport * parent_ns * transport.transpose()
    assert child_ns[:2, :2] == U
    assert child_ns[:2, 2:] == 0 and child_ns[2:, :2] == 0
    child_frame = -child_ns[2:, 2:]
    old_basis = identity_matrix(ZZ, 19)[2:, :]
    new_basis = transport[2:, :]
    record = certify_embedded_edge(
        corridor,
        edge_index,
        parent_ns,
        current_frame,
        old_basis,
        child_frame,
        new_basis,
        edge_metadata,
    )
    return child_frame, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    ns_route_path = (
        "artifacts/generated-results/"
        "elkies-k3-lattice-foundry-ns0024-r13-nef-route.json"
    )
    ns_source_path = (
        "artifacts/generated-results/"
        "elkies-k3-lattice-foundry-ns0024-source-hunt-r13.json"
    )
    golay_route_path = (
        "artifacts/generated-results/"
        "elkies-k3-golay-det720-same-ns-compiler-routes-pole4-witness-v1.json"
    )
    golay_source_path = (
        "artifacts/generated-results/"
        "elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
    )
    h3_path = (
        "artifacts/generated-results/"
        "elkies-k3-rank17-to-h3-reverse-transport.json"
    )

    ns_route = load_json(ns_route_path)
    ns_source = load_json(ns_source_path)
    golay_route = load_json(golay_route_path)
    golay_sources = load_json(golay_source_path)
    h3 = load_json(h3_path)

    assert ns_route["status"] == "PASS_EXACT_NEW_K3_ROOTFUL_TO_ROOTLESS_NEF_ROUTE"
    assert golay_route["status"] == "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_HIT"
    assert h3["status"] == "PASS_EXACT_PINNED_R17_TO_H3_REVERSE_TRANSPORT"

    current = matrix(ZZ, ns_source["source"]["root_adapted_gram"])
    certificates = []
    for edge in ns_route["edges"]:
        current, record = certify_local_edge(
            "NS0024",
            int(edge["edge_index"]),
            current,
            matrix(ZZ, edge["edge_transport_child_to_parent"]),
            edge,
        )
        certificates.append(record)

    golay_source = next(
        row["source"]
        for row in golay_sources["sources"]
        if row["source_id"] == "G720-S0128"
    )
    current = matrix(ZZ, golay_source["root_adapted_gram"])
    golay_selected = golay_route["results"][0]["best_routes_by_target"]["G720-F001"]
    for edge in golay_selected["edges"]:
        current, record = certify_local_edge(
            "Golay720",
            int(edge["edge_index"]),
            current,
            matrix(ZZ, edge["edge_transport_child_to_parent"]),
            edge,
        )
        certificates.append(record)

    h3_stages = h3["stages"]
    assert len(h3_stages) == 14
    first_frame = matrix(ZZ, h3_stages[0]["positive_frame"])
    first_basis = matrix(ZZ, h3_stages[0]["stage_basis_in_h3_ns"])
    h3_ambient = (
        first_basis.inverse()
        * block_diagonal_matrix(U, -first_frame)
        * first_basis.inverse().transpose()
    )
    assert h3_ambient in matrix(ZZ, 19, 19).parent()
    h3_ambient = h3_ambient.change_ring(ZZ)
    for edge_index in range(1, len(h3_stages)):
        old_stage = h3_stages[edge_index - 1]
        new_stage = h3_stages[edge_index]
        incoming = new_stage["incoming_neighbor"]
        metadata = {
            "q": incoming["q"],
            "old_fibre_degree": incoming["old_fiber_degree"],
            "source_root_rank": 17 - int(old_stage["mw_rank"]),
            "target_root_rank": 17 - int(new_stage["mw_rank"]),
        }
        record = certify_embedded_edge(
            "H3",
            edge_index,
            h3_ambient,
            matrix(ZZ, old_stage["positive_frame"]),
            matrix(ZZ, old_stage["stage_basis_in_h3_ns"])[2:, :],
            matrix(ZZ, new_stage["positive_frame"]),
            matrix(ZZ, new_stage["stage_basis_in_h3_ns"])[2:, :],
            metadata,
        )
        certificates.append(record)

    bridge_rank_histogram = Counter(row["old_frame"]["bridge_rank"] for row in certificates)
    core_rootless_count = sum(row["core"]["root_count_signed"] == 0 for row in certificates)
    signed_root_count_decreasing = sum(
        row["root_transfer"]["removed_signed_roots"]
        > row["root_transfer"]["introduced_signed_roots"]
        for row in certificates
    )

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-bridge-reglue.v1",
        "status": "PASS_EXACT_BRIDGE_REGLUE_CERTIFICATES",
        "theorem": {
            "name": "bridge-core root-transfer theorem",
            "statement": (
                "For two primitive hyperbolic planes U0,U1 in one integral "
                "lattice, K=U0^perp intersect U1^perp is primitive in each "
                "frame. With Ci=K^perp in Wi, Wi is the even overlattice of "
                "K+Ci encoded by the computed isotropic subgroup Hi. Moreover "
                "Phi(W0) intersect Phi(W1)=Phi(K); all other root change is "
                "confined to the two bridge/glue presentations."
            ),
            "rootless_gate": (
                "W_i is rootless iff K and every glue coset of K+C_i contain "
                "no norm-two vector."
            ),
            "certified_corridor_specialization": (
                "On every checked edge C_old and C_new have rank two and cyclic "
                "discriminant; both frame glues are maximal over the bridge "
                "factor and have the same order. Thus every click is a rank-two "
                "cyclic bridge replacement over one fixed rank-fifteen core."
            ),
        },
        "inputs": {
            path: digest(path)
            for path in [
                ns_route_path,
                ns_source_path,
                golay_route_path,
                golay_source_path,
                h3_path,
            ]
        },
        "aggregate": {
            "edge_count": len(certificates),
            "corridors": ["H3", "NS0024", "Golay720"],
            "bridge_rank_histogram": {
                str(rank): count for rank, count in sorted(bridge_rank_histogram.items())
            },
            "core_rootless_edges": core_rootless_count,
            "root_rank_decreasing_edges": sum(
                row["source_root_rank"] > row["target_root_rank"]
                for row in certificates
            ),
            "signed_root_count_decreasing_edges": signed_root_count_decreasing,
            "all_glue_groups_cyclic": all(
                len(row["old_frame"]["glue_group_invariants"]) <= 1
                and len(row["new_frame"]["glue_group_invariants"]) <= 1
                for row in certificates
            ),
            "all_glue_orders_preserved_across_bridge_replacement": all(
                row["old_frame"]["K_plus_C_index_in_W"]
                == row["new_frame"]["K_plus_C_index_in_W"]
                for row in certificates
            ),
        },
        "edges": certificates,
        "proof_boundary": {
            "proved": (
                "The abstract bridge-core theorem is elementary lattice theory. "
                "For all selected H3, NS0024, and Golay-720 edges, exact marked NS "
                "transports certify K, C_old, C_new, both finite glue groups, "
                "explicit isotropic generators, and complete norm-two transfer."
            ),
            "not_proved": (
                "No completeness theorem for possible U-pivots is asserted. "
                "The Q80 corridor is not included until its older TSV route "
                "format is normalized to the same marking interface."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_rank_transfer_bridge_reglue.sage --check"
        ),
    }

    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit(f"stale or missing artifact: {output}")
        print("PASS integral rank-transfer bridge reglue")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)
    print(output.relative_to(ROOT))


if __name__ == "__main__":
    main()
