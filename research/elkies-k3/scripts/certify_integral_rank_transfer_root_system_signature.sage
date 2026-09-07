#!/usr/bin/env sage-python
"""Certify metric physical-witness signatures for completed root systems.

The bridge controls reconstruct every norm-two vector as an exact ``k+c``
split vector and retain all pairwise inner products.  The Q80 controls check
that the same metric classifier recognizes the requested 4A1 and A1 targets.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import (
    GF,
    QQ,
    ZZ,
    block_diagonal_matrix,
    gcd,
    identity_matrix,
    lcm,
    matrix,
    pari,
    vector,
    xgcd,
)
from sage.coding.golay_code import GolayCode


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
MASKED = GENERATED / "elkies-k3-integral-rank-transfer-masked-core-controls-v1.json"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
NS_ROUTE_SCRIPT = ROOT / "elkies-k3/scripts/certify_ns0024_new_rootless_source_route.sage"
BASE_SCRIPT = ROOT / "elkies-k3/scripts/generate_integral_rank_transfer_masked_core_neighbors.sage"
SEARCH_SCRIPT = ROOT / "elkies-k3/scripts/search_integral_rank_transfer_masked_core_controls.sage"
CORE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_core_generation.sage"
REVERSE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_reverse_theta_masks.sage"
Q80_FRAME = ROOT / "elkies-k3/data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt"
Q80_PATH = ROOT / "elkies-k3/data/fibrations/kumar_q80_to_rootless_path.tsv"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-root-system-signature-v1.json"
U = matrix(ZZ, [[0, 1], [1, 0]])


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_vector(value):
    return [str(entry) for entry in value]


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def signed_roots(gram):
    result = pari(gram).qfminim(2)
    signed_count = int(result[0])
    if signed_count == 0:
        return []
    representatives = matrix(ZZ, result[2].sage()).columns()
    answer = []
    for representative in representatives:
        root = vector(ZZ, representative)
        answer.extend((root, -root))
    assert len(answer) == signed_count
    assert all(root * gram * root == 2 for root in answer)
    return answer


def canonical_root_lines(gram):
    lines = set()
    for root in signed_roots(gram):
        positive = tuple(map(int, root))
        negative = tuple(-entry for entry in positive)
        lines.add(min(positive, negative))
    return [vector(ZZ, root) for root in sorted(lines)]


def graph_components(adjacency):
    unseen = set(range(len(adjacency)))
    answer = []
    while unseen:
        seed = min(unseen)
        stack = [seed]
        component = []
        unseen.remove(seed)
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        answer.append(sorted(component))
    return answer


def irreducible_ade_type(rank, signed_count):
    if signed_count == rank * (rank + 1):
        return f"A{rank}"
    if rank >= 4 and signed_count == 2 * rank * (rank - 1):
        return f"D{rank}"
    exceptional = {(6, 72): "E6", (7, 126): "E7", (8, 240): "E8"}
    if (rank, signed_count) in exceptional:
        return exceptional[(rank, signed_count)]
    raise ArithmeticError(
        f"unclassified simply-laced component rank={rank}, roots={signed_count}"
    )


def aggregate_ade(labels):
    counts = Counter(labels)
    family_order = {"A": 0, "D": 1, "E": 2}
    ordered = sorted(counts, key=lambda label: (family_order[label[0]], int(label[1:])))
    return "+".join(
        (str(counts[label]) if counts[label] > 1 else "") + label
        for label in ordered
    ) or "rootless"


def discriminant_invariants(gram):
    if gram.nrows() == 0:
        return []
    return [
        abs(int(value))
        for value in gram.elementary_divisors()
        if abs(int(value)) > 1
    ]


def metric_root_signature(frame, split_gram=None, split_basis=None, glue=None):
    """Return the full root-line metric plus its marked embedding invariants."""

    roots = canonical_root_lines(frame)
    root_matrix = matrix(ZZ, roots) if roots else matrix(ZZ, 0, frame.nrows())
    pairings = root_matrix * frame * root_matrix.transpose()
    adjacency = [set() for _ in roots]
    edges = []
    for left in range(len(roots)):
        assert pairings[left, left] == 2
        for right in range(left + 1, len(roots)):
            pairing = int(pairings[left, right])
            assert pairing in (-1, 0, 1)
            if pairing:
                adjacency[left].add(right)
                adjacency[right].add(left)
                edges.append([left, right, pairing])

    component_rows = []
    labels = []
    component_discriminants = []
    for indices in graph_components(adjacency):
        vectors = matrix(ZZ, [roots[index] for index in indices])
        rank = int(vectors.rank())
        signed_count = 2 * len(indices)
        label = irreducible_ade_type(rank, signed_count)
        labels.append(label)
        component_root_module = vectors.row_module(ZZ)
        component_basis = component_root_module.basis_matrix()
        component_gram = component_basis * frame * component_basis.transpose()
        component_discriminant = abs(int(component_gram.det()))
        component_discriminants.append(component_discriminant)
        component_rows.append(
            {
                "type": label,
                "root_line_indices": indices,
                "rank": rank,
                "signed_root_count": signed_count,
                "root_lattice_discriminant": component_discriminant,
            }
        )

    root_rank = int(root_matrix.rank())
    if roots:
        root_module = root_matrix.row_module(ZZ)
        root_basis = root_module.basis_matrix()
        root_gram = root_basis * frame * root_basis.transpose()
        saturation = root_module.saturation()
        saturation_basis = saturation.basis_matrix()
        saturation_gram = saturation_basis * frame * saturation_basis.transpose()
        saturation_index = int(root_module.index_in(saturation))
        torsion_invariants = list(map(int, saturation.quotient(root_module).invariants()))
        root_determinant = abs(int(root_gram.det()))
        saturated_determinant = abs(int(saturation_gram.det()))
        assert root_determinant == saturation_index**2 * saturated_determinant
    else:
        root_basis = matrix(ZZ, 0, frame.nrows())
        saturation_basis = root_basis
        root_determinant = 1
        saturated_determinant = 1
        saturation_index = 1
        torsion_invariants = []

    ade = aggregate_ade(labels)
    component_product = 1
    for value in component_discriminants:
        component_product *= value
    assert component_product == root_determinant

    physical = None
    if split_gram is not None:
        assert split_basis is not None and glue is not None
        assert split_basis * split_gram * split_basis.transpose() == frame
        split_roots = root_matrix * split_basis
        split_pairings = split_roots * split_gram * split_roots.transpose()
        assert split_pairings == pairings
        glue_order = lcm(entry.denominator() for entry in glue)
        physical_rows = []
        for frame_root, split_root in zip(roots, split_roots.rows()):
            labels_for_root = [
                label
                for label in range(glue_order)
                if split_root - label * glue in ZZ ** split_gram.nrows()
            ]
            assert len(labels_for_root) == 1
            physical_rows.append(
                {
                    "frame_coordinates": list(map(int, frame_root)),
                    "core_dual_coordinates": rational_vector(split_root[:-2]),
                    "bridge_dual_coordinates": rational_vector(split_root[-2:]),
                    "graph_glue_label": labels_for_root[0],
                }
            )
        physical = {
            "split_rank": [split_gram.nrows() - 2, 2],
            "cyclic_glue_order": int(glue_order),
            "root_lines": physical_rows,
            "pairing_identity": "<k+c,k'+c'>=<k,k'>+<c,c'>",
            "pairing_identity_verified_for_every_pair": True,
        }

    triangular_pairings = [
        [int(pairings[left, right]) for right in range(left, len(roots))]
        for left in range(len(roots))
    ]
    metric_digest = hashlib.sha256(
        json.dumps(triangular_pairings, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "ade_type": ade,
        "root_rank": root_rank,
        "signed_root_count": 2 * len(roots),
        "root_line_count": len(roots),
        "root_lines_frame_coordinates": [list(map(int, root)) for root in roots],
        "pairwise_inner_products_upper_triangular": triangular_pairings,
        "pairwise_inner_products_sha256": metric_digest,
        "nonorthogonal_root_line_edges": edges,
        "components": component_rows,
        "root_lattice_basis": rows(root_basis),
        "root_lattice_discriminant": root_determinant,
        "root_lattice_discriminant_group_invariants": discriminant_invariants(
            root_basis * frame * root_basis.transpose()
        ),
        "primitive_closure_basis": rows(saturation_basis),
        "primitive_closure_index": saturation_index,
        "primitive_closure_discriminant": saturated_determinant,
        "primitive_closure_quotient_invariants": torsion_invariants,
        "mw_torsion_invariants_for_frame": torsion_invariants,
        "physical_bridge_witnesses": physical,
    }


def overlattice_with_basis(core, bridge, glue):
    split = block_diagonal_matrix(core, bridge)
    denominator = lcm(value.denominator() for value in glue)
    generators = (
        denominator * identity_matrix(QQ, split.nrows())
    ).stack(matrix(QQ, [denominator * glue])).change_ring(ZZ)
    basis = generators.row_module(ZZ).basis_matrix().change_ring(QQ) / denominator
    frame = basis * split * basis.transpose()
    assert all(value in ZZ for value in frame.list())
    frame = frame.change_ring(ZZ)
    assert not any(value % 2 for value in frame.diagonal())
    return frame, split, basis


def ns0024_controls():
    masked = json.loads(MASKED.read_text())
    bridges = json.loads(BRIDGES.read_text())
    theta = json.loads(THETA.read_text())
    ns_route = runpy.run_path(str(NS_ROUTE_SCRIPT))
    base = runpy.run_path(str(BASE_SCRIPT))
    search = runpy.run_path(str(SEARCH_SCRIPT))
    core_tools = runpy.run_path(str(CORE_SCRIPT))
    reverse = runpy.run_path(str(REVERSE_SCRIPT))

    prepared = search["prepare_corridor"](
        "NS0024", bridges, theta, base, core_tools, reverse
    )
    search["configure_order"](base, prepared["order"])
    masked_row = next(row for row in masked["corridors"] if row["corridor"] == "NS0024")
    bridge_index = masked_row["completion"]["bridge_class_index"]
    terminal_multiplier = masked_row["completion"]["glue_multiplier"]
    bridge = next(
        row for row in prepared["viable_bridges"]
        if row["bridge_class_index"] == bridge_index
    )

    quadratic = base["quadratic_form"](prepared["seed"])
    cores = [base["lll_reduce"](quadratic.Hessian_matrix())]
    for prime, raw_vector in ns_route["PATH"]:
        witness = vector(ZZ, raw_vector)
        assert quadratic(witness) % prime == 0
        quadratic = quadratic.find_p_neighbor_from_vec(prime, witness)
        cores.append(base["lll_reduce"](quadratic.Hessian_matrix()))

    expected = ((13, 280, 4), (5, 12, 24), (5, 12, 24), (0, 0, 1))
    rows_out = []
    for stage, core_gram in enumerate(cores):
        choices = ns_route["completed_frames"](
            core_gram, bridge, prepared["order"], base, core_tools
        )
        matches = [
            (multiplier, frame)
            for multiplier, frame in choices
            if ns_route["root_data"](frame) == expected[stage]
        ]
        assert matches
        if stage == len(cores) - 1:
            multiplier, expected_frame = next(
                row for row in matches if row[0] == terminal_multiplier
            )
        else:
            multiplier, expected_frame = min(matches, key=lambda row: row[0])
        core_generator = base["primary_generator"](core_gram, prepared["order"])
        glue = vector(
            QQ,
            list(multiplier * core_generator) + list(bridge["generator"]),
        )
        frame, split, basis = overlattice_with_basis(core_gram, bridge["gram"], glue)
        assert frame == expected_frame
        signature = metric_root_signature(frame, split, basis, glue)
        rows_out.append(
            {
                "stage": stage,
                "incoming_core_neighbor_prime": (
                    None if stage == 0 else ns_route["PATH"][stage - 1][0]
                ),
                "glue_multiplier": int(multiplier),
                "signature": signature,
            }
        )

    assert [row["signature"]["ade_type"] for row in rows_out] == [
        "D5+E8", "3A1+A2", "3A1+A2", "rootless"
    ]
    assert all(
        row["signature"]["primitive_closure_index"] == 1 for row in rows_out
    )
    assert {
        witness["graph_glue_label"]
        for witness in rows_out[2]["signature"]["physical_bridge_witnesses"]["root_lines"]
    } != {0}
    return rows_out


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    return vector(ZZ, coefficients if current == 1 else [-value for value in coefficients])


def elliptic_neighbor(parent, qnorm, left, right, coordinates):
    ns = block_diagonal_matrix(U, -parent)
    fiber = vector(ZZ, [left, right] + list(coordinates))
    assert left * right == qnorm
    assert coordinates * parent * coordinates == 2 * qnorm
    assert gcd([abs(ZZ(value)) for value in ns * fiber]) == 1
    mate = bezout_vector(list(ns * fiber))
    mate -= ZZ(mate * ns * mate) // 2 * fiber
    complement = matrix(ZZ, [list(fiber * ns), list(mate * ns)]).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    transport = matrix(ZZ, [list(fiber), list(mate)] + complement.rows())
    assert abs(transport.det()) == 1
    return child


def q80_controls():
    frame = load_matrix(Q80_FRAME)
    answer = []
    with Q80_PATH.open() as handle:
        steps = list(csv.DictReader(handle, delimiter="\t"))
    for row in steps:
        frame = elliptic_neighbor(
            frame,
            ZZ(row["q"]),
            ZZ(row["a"]),
            ZZ(row["b"]),
            vector(ZZ, map(ZZ, row["v"].split(","))),
        )
        if row["ADE"] not in ("4A1", "A1"):
            continue
        signature = metric_root_signature(frame)
        assert signature["ade_type"] == row["ADE"]
        assert signature["root_rank"] == int(row["root_rank"])
        assert signature["signed_root_count"] == int(row["roots"])
        assert signature["root_lattice_discriminant"] == int(row["rootdet"])
        assert signature["primitive_closure_index"] == 1
        answer.append({"q80_path_step": int(row["step"]), "signature": signature})
    assert [row["signature"]["ade_type"] for row in answer] == ["4A1", "A1"]
    return answer


def marked_embedding_counterexample():
    """Verify that the abstract 24A1 metric does not determine saturation."""

    split_signature = metric_root_signature(2 * identity_matrix(ZZ, 24))
    code = GolayCode(GF(2), extended=True)
    generator = code.generator_matrix().echelon_form()
    assert generator[:, :12] == identity_matrix(GF(2), 12)
    physical_basis = matrix(QQ, 24, 24)
    for index in range(12):
        physical_basis[index] = vector(
            QQ, [QQ(int(entry)) / 2 for entry in generator[index]]
        )
        physical_basis[12 + index, 12 + index] = 1
    niemeier = physical_basis * (2 * identity_matrix(QQ, 24)) * physical_basis.transpose()
    assert niemeier in matrix(ZZ, 24, 24).parent()
    niemeier = niemeier.change_ring(ZZ)
    assert niemeier.det() == 1
    niemeier_signature = metric_root_signature(niemeier)
    assert code.dimension() == 12 and code.minimum_distance() == 8
    assert split_signature["ade_type"] == niemeier_signature["ade_type"] == "24A1"
    assert (
        split_signature["pairwise_inner_products_upper_triangular"]
        == niemeier_signature["pairwise_inner_products_upper_triangular"]
    )
    assert split_signature["primitive_closure_index"] == 1
    assert split_signature["primitive_closure_quotient_invariants"] == []
    assert niemeier_signature["primitive_closure_index"] == 2**12
    assert niemeier_signature["primitive_closure_quotient_invariants"] == [2] * 12
    return {
        "common_metric_ade_type": "24A1",
        "common_pairwise_inner_products_sha256": split_signature[
            "pairwise_inner_products_sha256"
        ],
        "split_A1_24": {
            "ambient_determinant": 2**24,
            "root_primitive_closure_index": 1,
            "root_primitive_closure_quotient_invariants": [],
        },
        "niemeier_N_24A1": {
            "ambient_determinant": int(niemeier.det()),
            "root_primitive_closure_index": 2**12,
            "root_primitive_closure_quotient_invariants": [2] * 12,
            "quotient_identification": "extended binary Golay code",
        },
        "conclusion": (
            "The common pairwise root metric determines 24A1 but not its "
            "primitive closure in the marked ambient lattice."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    ns0024 = ns0024_controls()
    q80 = q80_controls()
    counterexample = marked_embedding_counterexample()
    input_paths = (
        MASKED,
        BRIDGES,
        THETA,
        NS_ROUTE_SCRIPT,
        BASE_SCRIPT,
        SEARCH_SCRIPT,
        CORE_SCRIPT,
        REVERSE_SCRIPT,
        Q80_FRAME,
        Q80_PATH,
    )
    payload = {
        "schema": "elkies-k3.integral-rank-transfer-root-system-signature.v1",
        "status": "PASS_EXACT_METRIC_PHYSICAL_WITNESS_SIGNATURE",
        "theorem_boundary": {
            "metric_signature_determines": [
                "complete signed norm-two root system",
                "nonorthogonal root-line graph",
                "ADE decomposition",
                "root rank",
                "root-lattice and component discriminants",
            ],
            "marked_embedding_additionally_determines": [
                "primitive closure of the root lattice in the frame",
                "primitive-closure quotient",
                "exact Mordell-Weil torsion contribution",
            ],
            "pairwise_metric_alone_does_not_determine": [
                "primitive closure in an ambient frame",
                "exact Mordell-Weil torsion",
            ],
        },
        "ns0024_physical_bridge_controls": ns0024,
        "q80_metric_target_controls": q80,
        "marked_embedding_counterexample": counterexample,
        "target_types_verified": ["D5+E8", "3A1+A2", "4A1", "A1"],
        "inputs": {relative(path): digest(path) for path in input_paths},
        "reproduce": [
            "sage -python elkies-k3/scripts/certify_integral_rank_transfer_root_system_signature.sage",
            "sage -python elkies-k3/scripts/certify_integral_rank_transfer_root_system_signature.sage --check",
        ],
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"missing or stale artifact: {output}")
        print("PASS metric physical-witness root-system signature")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
