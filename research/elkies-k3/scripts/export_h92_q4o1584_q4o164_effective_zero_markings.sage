#!/usr/bin/env sage -python
"""Export all equation-effective zero markings after q4/o1584 -> q4/o164."""

import hashlib
import json
from pathlib import Path

from sage.all import (
    Polyhedron, QQ, ZZ, block_diagonal_matrix, matrix, pari, vector,
)


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE = GENERATED / "elkies-k3-h3-q4o208-physical-q4o1584-second_i6_affine_component-marking.json"
CERT = GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-a1a1a3a3-certificate.json"
SCORES = GENERATED / "elkies-k3-h3-q4o208-q4o1584-second_i6_affine_component-q4d2-equation-cost.json"
OUTPUT = GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-effective-zero-markings.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))
FIBRE = vector(ZZ, [2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 0, 0, 0])


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compatible_simple_roots(complement, known_positive):
    """Return a simple system whose positive cone contains known_positive."""
    half = [vector(ZZ, column) for column in matrix(ZZ, pari(complement).qfminim(2)[2]).columns()]
    roots = half + [-root for root in half]
    inequalities = [[-1] + entries(root * complement) for root in known_positive]
    point = vector(QQ, Polyhedron(ieqs=inequalities).representative_point())
    # Move off any residual root hyperplane without leaving the known-positive
    # cone.  A sufficiently large multiple preserves every strict inequality.
    functional = None
    for scale in range(1, 1000):
        base = scale * point
        for tweak in [vector(QQ, [0] * complement.nrows())] + [
            vector(QQ, [1 if i == j else 0 for i in range(complement.nrows())])
            for j in range(complement.nrows())
        ]:
            trial = base + tweak
            if all(trial * complement * root > 0 for root in known_positive) and all(
                trial * complement * root != 0 for root in roots
            ):
                functional = trial
                break
        if functional is not None:
            break
    assert functional is not None
    positive = [root for root in roots if functional * complement * root > 0]
    positive_set = {tuple(root) for root in positive}
    simple = [
        root for root in positive
        if not any(tuple(root - left) in positive_set for left in positive)
    ]
    simple_matrix = matrix(ZZ, [list(root) for root in simple])
    root_rank = matrix(ZZ, [list(root) for root in roots]).rank()
    assert simple_matrix.nrows() == simple_matrix.rank() == root_rank
    for root in known_positive:
        coefficients = simple_matrix.solve_left(root)
        assert all(value in ZZ and value >= 0 for value in coefficients)
    return simple_matrix


source = json.loads(SOURCE.read_text())
cert = json.loads(CERT.read_text())
scores = json.loads(SCORES.read_text())
assert source["status"] == "PASS_EXACT_Q4O1584_PHYSICAL_EFFECTIVE_ZERO_MARKING"
assert cert["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert cert["candidate_id"]["label"] == "q4o164_after_q4o1584"
best = scores["best_candidate"]
assert best["candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 164}

source_frame = load_matrix(ROOT / source["frame_output"])
source_gram = block_diagonal_matrix(U2, -source_frame)
explicit = {
    name: vector(ZZ, value)
    for name, value in source["equation_explicit_curves_in_child"].items()
}
assert FIBRE == vector(ZZ, best["fibre"])
assert vector(ZZ, best["horizontal"]["section"]) == explicit["old_A11_component_0"]
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
assert old_zero == explicit["second_I6_affine_component"]
assert FIBRE - old_zero - explicit["old_A11_component_0"] == vector(
    ZZ, [2, 0] + [0, 0, 0, 0, 0, 0, 0, 1, 0, 1] + [0] * 7
)

zero_names = sorted(
    name for name, curve in explicit.items()
    if curve * source_gram * FIBRE == 1
)
outputs = {}
for zero_name in zero_names:
    zero = explicit[zero_name]
    u_rows = matrix(ZZ, [list(FIBRE), list(zero + FIBRE)])
    kernel = matrix(ZZ, [list(FIBRE * source_gram), list((zero + FIBRE) * source_gram)]).right_kernel_matrix()
    preliminary = u_rows.stack(kernel)
    assert abs(preliminary.det()) == 1
    prelim_inverse = preliminary.inverse().change_ring(ZZ)
    prelim_gram = preliminary * source_gram * preliminary.transpose()
    complement = -prelim_gram[2:, 2:]
    root_result = pari(complement).qfminim(2)
    assert ZZ(root_result[0]) == 28

    known_positive_source = []
    for curve in {tuple(value): value for value in explicit.values()}.values():
        if curve * source_gram * FIBRE != 0:
            continue
        incidence = ZZ(zero * source_gram * curve)
        assert incidence in (0, 1)
        positive = curve if incidence == 0 else FIBRE - curve
        assert positive * source_gram * FIBRE == 0
        assert positive * source_gram * zero == 0
        assert positive * source_gram * positive == -2
        if tuple(positive) not in {tuple(item) for item in known_positive_source}:
            known_positive_source.append(positive)
    known_positive = [
        vector(ZZ, (curve * prelim_inverse)[2:])
        for curve in known_positive_source
    ]
    simple_tail = compatible_simple_roots(complement, known_positive)
    simple_source = [vector(ZZ, [0, 0] + list(root)) * preliminary for root in simple_tail.rows()]
    cartan = simple_tail * complement * simple_tail.transpose()
    assert abs(cartan.det()) == 64 and ZZ(pari(cartan).qfminim(2)[0]) == 28

    # The search convention has basis row -C for each effective simple curve.
    partial = u_rows.stack(matrix(ZZ, [list(-curve) for curve in simple_source]))
    smith, left, right = partial.smith_form()
    assert smith[:, :10] == matrix.identity(ZZ, 10)
    assert not any(smith[:, 10:].list())
    completion = right.inverse().change_ring(ZZ)
    change = block_diagonal_matrix(left.inverse().change_ring(ZZ), matrix.identity(ZZ, 9))
    basis_in_source = change * completion
    assert basis_in_source[:10] == partial
    for index in range(10, 19):
        row = basis_in_source.row(index)
        row -= ZZ(row * source_gram * (zero + FIBRE)) * FIBRE
        row -= ZZ(row * source_gram * FIBRE) * (zero + FIBRE)
        basis_in_source[index] = row
    assert abs(basis_in_source.det()) == 1
    source_in_basis = basis_in_source.inverse().change_ring(ZZ)
    gram = basis_in_source * source_gram * basis_in_source.transpose()
    assert gram[:2, :2] == U2 and not any(gram[:2, 2:].list())
    frame = -gram[2:, 2:]
    assert frame[:8, :8] == cartan

    slug = zero_name.lower()
    frame_path = GENERATED / f"elkies-k3-h3-q4o208-q4o1584-q4o164-{slug}-frame.txt"
    marking_path = GENERATED / f"elkies-k3-h3-q4o208-q4o1584-q4o164-{slug}-marking.json"
    frame_path.write_text(
        f"# q4/o164 physical 2A1+2A3 frame; zero={zero_name}\n"
        + "\n".join(" ".join(map(str, row)) for row in frame.rows()) + "\n"
    )
    explicit_child = {name: entries(curve * source_in_basis) for name, curve in explicit.items()}
    targets = {
        name: entries(vector(ZZ, value) * source_in_basis)
        for name, value in source["target_fibres_in_root_adapted_hub"].items()
    }
    equation_to_source = matrix(ZZ, source["equation_A11_to_root_adapted_hub_basis"])
    equation_to_child = basis_in_source * equation_to_source
    payload = {
        "schema": "elkies-k3.h3-q4o1584-q4o164-effective-zero-marking.v1",
        "status": "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING",
        "hub": f"q4o208_q4o1584_q4o164_physical_2A1_2A3_{zero_name}_zero",
        "zero": zero_name,
        "root_data": [8, 28, 64],
        "ade": "2A1+2A3",
        "frame_output": str(frame_path.relative_to(ROOT)),
        "frame_sha256": sha256(frame_path),
        "physical_simple_components_in_source": [entries(curve) for curve in simple_source],
        "basis_in_source": rows(basis_in_source),
        "source_in_basis": rows(source_in_basis),
        "equation_A11_to_root_adapted_hub_basis": rows(equation_to_child),
        "target_fibres_in_root_adapted_hub": targets,
        "equation_explicit_curves_in_child": explicit_child,
        "prefix_operational_score": (source.get("prefix_operational_score") or 0) - 1620,
        "literal_divisor": {
            "fibre_in_source": entries(FIBRE),
            "zero": "second_I6_affine_component",
            "horizontal": "old_A11_component_0",
            "effective_vertical": "(F-r7)+(F-r9)",
            "P_dot_O": 0,
            "expected_RR_ambient": 4,
        },
        "proof_boundary": (
            "Exact physical effective zero and component chamber compatible with every "
            "inherited vertical curve, literal attached horizontal, marked targets, and "
            "bidirectional unimodular NS transport."
        ),
        "inputs": {
            "paths": [str(path.relative_to(ROOT)) for path in (SOURCE, CERT, SCORES)],
            "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (SOURCE, CERT, SCORES)},
        },
    }
    marking_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    outputs[zero_name] = {
        "frame": str(frame_path.relative_to(ROOT)),
        "marking": str(marking_path.relative_to(ROOT)),
        "frame_max_abs": max(abs(int(value)) for value in basis_in_source.list()),
    }

summary = {
    "schema": "elkies-k3.h3-q4o1584-q4o164-effective-zero-markings.v1",
    "status": "PASS_EXACT_Q4O164_ALL_PHYSICAL_EFFECTIVE_ZERO_MARKINGS",
    "edge": {
        "q": 4,
        "orbit": 164,
        "ade": "2A1+2A3",
        "mw_rank": 9,
        "literal_divisor": "O(second_I6_affine)+C0+(F-r7)+(F-r9)",
    },
    "outputs": outputs,
    "proof_boundary": "Each output is a complete physical marked state; no suffix is inferred from ADE/MW.",
}
OUTPUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(f"Q4O164ZEROS|count={len(outputs)}|status={summary['status']}|output={OUTPUT}")
