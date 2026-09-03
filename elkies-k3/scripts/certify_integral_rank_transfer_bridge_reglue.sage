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
import csv
import hashlib
import json
from pathlib import Path

from sage.all import (
    Genus,
    QQ,
    ZZ,
    block_diagonal_matrix,
    gcd,
    identity_matrix,
    matrix,
    pari,
    vector,
    xgcd,
)


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
)
RELATIVE_U_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-relative-u-bridge-lifting-regression-v1.json"
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


def discriminant_form_key(gram):
    normal = Genus(gram).discriminant_form().normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rational_rows(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def load_matrix(relative_path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in (ROOT / relative_path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        if not pairing:
            continue
        common, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = common
    assert abs(current) == 1
    return vector(
        ZZ, coefficients if current == 1 else [-value for value in coefficients]
    )


def split_neighbor(parent, q, a, b, coordinates):
    ns = block_diagonal_matrix(U, -parent)
    coordinates = vector(ZZ, coordinates)
    fibre = vector(ZZ, [a, b] + list(coordinates))
    assert a * b == q
    assert coordinates * parent * coordinates == 2 * q
    assert fibre * ns * fibre == 0
    assert gcd([abs(ZZ(value)) for value in ns * fibre]) == 1
    mate = bezout_vector(list(ns * fibre))
    mate -= ZZ(mate * ns * mate) // 2 * fibre
    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    transport = matrix(ZZ, [list(fibre), list(mate)] + complement.rows())
    assert abs(transport.det()) == 1
    assert transport * ns * transport.transpose() == block_diagonal_matrix(U, -child)
    return child, transport


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
                "primary_components": [
                    {
                        "prime": int(prime),
                        "order": int(prime**exponent),
                        "frame_coordinates": [
                            int((order // (prime**exponent)) * value)
                            for value in frame_coordinates
                        ],
                    }
                    for prime, exponent in ZZ(order).factor()
                ],
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
        "bridge_discriminant_form": discriminant_form_key(bridge_gram),
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


def relative_u_projection(
    ambient_ns,
    source_u_basis,
    source_frame,
    source_frame_basis,
    target_u_basis,
    core_basis,
    stored_bridge_determinant,
):
    """Certify the ordered relative-U projection and its saturated bridge."""
    assert source_u_basis * ambient_ns * source_u_basis.transpose() == U
    assert target_u_basis * ambient_ns * target_u_basis.transpose() == U
    assert source_u_basis * ambient_ns * source_frame_basis.transpose() == 0

    cross_pairing = source_u_basis * ambient_ns * target_u_basis.transpose()
    source_projection = cross_pairing.transpose() * U * source_u_basis
    residual = target_u_basis - source_projection
    assert residual * ambient_ns * source_u_basis.transpose() == 0

    positive_gram = cross_pairing.transpose() * U * cross_pairing - U
    assert -(residual * ambient_ns * residual.transpose()) == positive_gram
    assert positive_gram.is_positive_definite()

    frame_coordinates = source_frame_basis.solve_left(residual)
    assert frame_coordinates in matrix(ZZ, 2, source_frame.nrows()).parent()
    frame_coordinates = frame_coordinates.change_ring(ZZ)
    assert frame_coordinates * source_frame_basis == residual
    assert frame_coordinates * source_frame * frame_coordinates.transpose() == positive_gram

    core_coordinates = source_frame_basis.solve_left(core_basis).change_ring(ZZ)
    saturated_bridge = (core_coordinates * source_frame).right_kernel_matrix()
    raw_module = frame_coordinates.row_module(ZZ)
    saturated_module = raw_module.saturation()
    stored_module = saturated_bridge.row_module(ZZ)
    assert saturated_module == stored_module
    saturation_index = abs(int(raw_module.index_in(saturated_module)))
    assert abs(int(positive_gram.det())) == (
        saturation_index**2 * stored_bridge_determinant
    )

    d = int(cross_pairing[0, 0])
    s = int(cross_pairing[0, 1] - d)
    t = int(cross_pairing[1, 0] - d)
    z = int(cross_pairing[1, 1] - d - s - t)
    assert cross_pairing == matrix(
        ZZ,
        [[d, d + s], [d + t, d + s + t + z]],
    )
    return {
        "cross_pairing_A": rows(cross_pairing),
        "intersection_coordinates": {
            "F_dot_F_prime": d,
            "F_dot_O_prime": s,
            "O_dot_F_prime": t,
            "O_dot_O_prime": z,
        },
        "positive_projection_gram_G_A": rows(positive_gram),
        "projected_vectors_in_source_frame": rows(frame_coordinates),
        "raw_bridge_determinant": abs(int(positive_gram.det())),
        "saturation_index": saturation_index,
        "saturated_bridge_gram": rows(
            saturated_bridge * source_frame * saturated_bridge.transpose()
        ),
        "saturated_bridge_determinant": stored_bridge_determinant,
        "determinant_square_index_identity": (
            "det(G_A) = saturation_index^2 * det(C)"
        ),
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
    old_u_basis=None,
    new_u_basis=None,
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

    relative_u = None
    if old_u_basis is not None or new_u_basis is not None:
        assert old_u_basis is not None and new_u_basis is not None
        old_to_new = relative_u_projection(
            ambient_ns,
            old_u_basis,
            old_frame,
            old_basis,
            new_u_basis,
            core_basis,
            old_data["bridge_determinant_absolute"],
        )
        new_to_old = relative_u_projection(
            ambient_ns,
            new_u_basis,
            new_frame,
            new_basis,
            old_u_basis,
            core_basis,
            new_data["bridge_determinant_absolute"],
        )
        assert old_to_new["cross_pairing_A"] == rows(
            matrix(ZZ, new_to_old["cross_pairing_A"]).transpose()
        )
        assert (
            old_to_new["intersection_coordinates"]["F_dot_F_prime"]
            == int(edge_metadata["old_fibre_degree"])
        )
        relative_u = {
            "old_to_new": old_to_new,
            "new_to_old": new_to_old,
            "old_fibre_degree_matches_edge_record": True,
        }

    record = {
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
            "discriminant_form": discriminant_form_key(core_gram),
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
    if relative_u is not None:
        record["relative_u"] = relative_u
    return record


def certify_local_edge(
    corridor,
    edge_index,
    current_frame,
    transport,
    edge_metadata,
    include_relative_u=False,
):
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
        old_u_basis=identity_matrix(ZZ, 19)[:2, :] if include_relative_u else None,
        new_u_basis=transport[:2, :] if include_relative_u else None,
    )
    return child_frame, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--relative-u-output",
        type=Path,
        help="write the opt-in relative-U regression without changing the pinned v1 artifact",
    )
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
    q80_start_path = "elkies-k3/data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt"
    q80_prefix_path = "elkies-k3/data/fibrations/kumar_q80_lowq_alternate_prefix.tsv"
    q80_suffix_path = "elkies-k3/data/fibrations/kumar_q80_new_lowq_rootless_path.tsv"

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
            include_relative_u=arguments.relative_u_output is not None,
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
            include_relative_u=arguments.relative_u_output is not None,
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
            old_u_basis=(
                matrix(ZZ, old_stage["stage_basis_in_h3_ns"])[:2, :]
                if arguments.relative_u_output is not None
                else None
            ),
            new_u_basis=(
                matrix(ZZ, new_stage["stage_basis_in_h3_ns"])[:2, :]
                if arguments.relative_u_output is not None
                else None
            ),
        )
        certificates.append(record)

    q80_steps = []
    for path, count in ((q80_prefix_path, 2), (q80_suffix_path, 8)):
        with (ROOT / path).open() as handle:
            selected_rows = list(csv.DictReader(handle, delimiter="\t"))[:count]
        assert len(selected_rows) == count
        q80_steps.extend(selected_rows)
    assert len(q80_steps) == 10
    current_frame = load_matrix(q80_start_path)
    common_ns = block_diagonal_matrix(U, -current_frame)
    current_basis = identity_matrix(ZZ, 19)
    current_root_rank = 14
    for edge_index, step in enumerate(q80_steps, 1):
        q = ZZ(step["q"])
        a = ZZ(step["a"])
        b = ZZ(step["b"])
        coordinates = vector(ZZ, map(ZZ, step["v"].split(",")))
        child_frame, local_transport = split_neighbor(
            current_frame, q, a, b, coordinates
        )
        new_basis = local_transport * current_basis
        target_root_rank = int(step["root_rank"])
        metadata = {
            "q": q,
            "old_fibre_degree": b,
            "source_root_rank": current_root_rank,
            "target_root_rank": target_root_rank,
        }
        record = certify_embedded_edge(
            "Q80",
            edge_index,
            common_ns,
            current_frame,
            current_basis[2:, :],
            child_frame,
            new_basis[2:, :],
            metadata,
            old_u_basis=(
                current_basis[:2, :]
                if arguments.relative_u_output is not None
                else None
            ),
            new_u_basis=(
                new_basis[:2, :]
                if arguments.relative_u_output is not None
                else None
            ),
        )
        certificates.append(record)
        current_frame = child_frame
        current_basis = new_basis
        current_root_rank = target_root_rank
    assert current_root_rank == 0 and not signed_roots(current_frame)

    bridge_rank_histogram = Counter(row["old_frame"]["bridge_rank"] for row in certificates)
    corridor_edge_counts = Counter(row["corridor"] for row in certificates)
    assert corridor_edge_counts == {
        "H3": 13,
        "Q80": 10,
        "NS0024": 13,
        "Golay720": 6,
    }
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
                q80_start_path,
                q80_prefix_path,
                q80_suffix_path,
            ]
        },
        "aggregate": {
            "edge_count": len(certificates),
            "corridors": ["H3", "Q80", "NS0024", "Golay720"],
            "corridor_edge_counts": {
                corridor: count
                for corridor, count in sorted(corridor_edge_counts.items())
            },
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
                "For all selected H3, Q80, NS0024, and Golay-720 edges, exact marked NS "
                "transports certify K, C_old, C_new, both finite glue groups, "
                "explicit isotropic generators, and complete norm-two transfer."
            ),
            "not_proved": (
                "No completeness theorem for possible U-pivots is asserted, "
                "and no equation-level realization follows from this lattice replay."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_rank_transfer_bridge_reglue.sage --check"
        ),
    }

    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.relative_u_output is not None:
        relative_edges = [
            {
                "corridor": row["corridor"],
                "edge_index": row["edge_index"],
                "q": row["q"],
                "old_fibre_degree": row["old_fibre_degree"],
                "source_root_rank": row["source_root_rank"],
                "target_root_rank": row["target_root_rank"],
                "stored_bridges": {
                    "old_gram": row["old_frame"]["bridge_gram"],
                    "old_determinant": row["old_frame"][
                        "bridge_determinant_absolute"
                    ],
                    "new_gram": row["new_frame"]["bridge_gram"],
                    "new_determinant": row["new_frame"][
                        "bridge_determinant_absolute"
                    ],
                },
                "relative_u": row["relative_u"],
            }
            for row in certificates
        ]
        relative_payload = {
            "schema": "elkies-k3.relative-u-bridge-lifting-regression.v1",
            "status": "PASS_EXACT_RELATIVE_U_BRIDGE_LIFTING_REGRESSION",
            "theorem": {
                "projection_formula": "G_A = A^T J A - J",
                "lift_formula": "u'_j = (J A)_{bullet j} + w_j",
                "saturation_formula": "det(G_A) = [C:<w_1,w_2>]^2 det(C)",
            },
            "aggregate": {
                "edge_count": len(certificates),
                "all_old_fibre_degrees_recovered": all(
                    row["relative_u"]["old_fibre_degree_matches_edge_record"]
                    for row in certificates
                ),
                "saturation_index_histogram_old_to_new": {
                    str(index): count
                    for index, count in sorted(
                        Counter(
                            row["relative_u"]["old_to_new"]["saturation_index"]
                            for row in certificates
                        ).items()
                    )
                },
                "saturation_index_histogram_new_to_old": {
                    str(index): count
                    for index, count in sorted(
                        Counter(
                            row["relative_u"]["new_to_old"]["saturation_index"]
                            for row in certificates
                        ).items()
                    )
                },
            },
            "inputs": payload["inputs"],
            "edges": relative_edges,
            "proof_boundary": {
                "proved": (
                    "For all 42 stored ordered U-pairs, exact integral projection gives "
                    "G_A=A^T J A-J; saturating the projected pair recovers the stored "
                    "rank-two bridge on both orientations, with the determinant corrected "
                    "by the square of the recorded saturation index. The (1,1) entry of A "
                    "recovers every stored old-fibre degree."
                ),
                "not_proved": (
                    "The transported second isotropic basis vector need not be the "
                    "equation-effective zero marking, and this regression does not enumerate "
                    "new relative U embeddings or construct equations."
                ),
            },
            "reproduce": (
                "sage -python elkies-k3/scripts/"
                "certify_integral_rank_transfer_bridge_reglue.sage "
                "--relative-u-output artifacts/generated-results/"
                "elkies-k3-relative-u-bridge-lifting-regression-v1.json"
            ),
        }
        relative_serialized = json.dumps(relative_payload, indent=2, sort_keys=True) + "\n"
        relative_output = (
            arguments.relative_u_output
            if arguments.relative_u_output.is_absolute()
            else ROOT / arguments.relative_u_output
        )
        relative_output.parent.mkdir(parents=True, exist_ok=True)
        relative_output.write_text(relative_serialized)
        try:
            print(relative_output.relative_to(ROOT))
        except ValueError:
            print(relative_output)
        return
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit(f"stale or missing artifact: {output}")
        print("PASS integral rank-transfer bridge reglue")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)
    try:
        print(output.relative_to(ROOT))
    except ValueError:
        print(output)


if __name__ == "__main__":
    main()
