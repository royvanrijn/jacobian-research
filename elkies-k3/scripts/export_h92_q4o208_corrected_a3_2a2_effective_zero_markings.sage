#!/usr/bin/env sage -python
"""Export every physical effective-zero marking of corrected q4/o323."""

import hashlib
import json
from pathlib import Path

from sage.all import Polyhedron, QQ, ZZ, block_diagonal_matrix, matrix, pari, vector

ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
CERT = GENERATED / "elkies-k3-h3-q4o208-physical-q4o323-corrected-a3-2a2-certificate.json"
OUTPUT = GENERATED / "elkies-k3-h3-q4o208-corrected-a3-2a2-effective-zero-markings.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def entries(value):
    return [int(x) for x in vector(ZZ, value)]


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_simple_system(complement, known_positive):
    half = [vector(ZZ, col) for col in matrix(ZZ, pari(complement).qfminim(2)[2]).columns()]
    roots = half + [-root for root in half]
    cone = Polyhedron(ieqs=[[-1] + entries(root * complement) for root in known_positive])
    point = vector(QQ, cone.representative_point())
    functional = None
    tweaks = [vector(QQ, [0] * complement.nrows())] + [
        vector(QQ, [int(i == j) for i in range(complement.nrows())])
        for j in range(complement.nrows())
    ]
    for scale in range(1, 1000):
        for tweak in tweaks:
            trial = scale * point + tweak
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
    simple = matrix(ZZ, [list(root) for root in positive
                         if not any(tuple(root-left) in positive_set for left in positive)])
    root_rank = matrix(ZZ, [list(root) for root in roots]).rank()
    assert simple.nrows() == simple.rank() == root_rank
    for root in known_positive:
        coefficients = simple.solve_left(root)
        assert all(x in ZZ and x >= 0 for x in coefficients)
    return simple


source = json.loads(SOURCE.read_text())
cert = json.loads(CERT.read_text())
assert source["status"] == "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING"
assert cert["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
source_frame = load_matrix(ROOT / source["frame_output"])
source_gram = block_diagonal_matrix(U2, -source_frame)
fibre = vector(ZZ, cert["source_to_child_basis"][0])
assert fibre[1] == 2 and fibre * source_gram * fibre == 0
explicit = {name: vector(ZZ, value) for name, value in source["equation_explicit_curves_in_child"].items()}
zero_names = sorted(name for name, curve in explicit.items() if curve * source_gram * fibre == 1)
assert zero_names == sorted(["old_A11_component_1", "old_A11_component_10", "old_A11_component_2",
                             "old_A11_component_3", "old_A11_component_4", "old_A11_component_9", "old_zero"])

outputs = {}
for zero_name in zero_names:
    zero = explicit[zero_name]
    u_rows = matrix(ZZ, [list(fibre), list(zero + fibre)])
    kernel = matrix(ZZ, [list(fibre * source_gram), list((zero + fibre) * source_gram)]).right_kernel_matrix()
    preliminary = u_rows.stack(kernel)
    assert abs(preliminary.det()) == 1
    preliminary_inverse = preliminary.inverse().change_ring(ZZ)
    prelim_gram = preliminary * source_gram * preliminary.transpose()
    complement = -prelim_gram[2:, 2:]
    assert ZZ(pari(complement).qfminim(2)[0]) == 24

    known_source = []
    for curve in {tuple(value): value for value in explicit.values()}.values():
        if curve * source_gram * fibre != 0:
            continue
        incidence = ZZ(zero * source_gram * curve)
        assert incidence in (0, 1)
        positive = curve if incidence == 0 else fibre - curve
        assert positive * source_gram * positive == -2 and positive * source_gram * zero == 0
        if tuple(positive) not in {tuple(item) for item in known_source}:
            known_source.append(positive)
    known_tail = [vector(ZZ, (curve * preliminary_inverse)[2:]) for curve in known_source]
    simple_tail = physical_simple_system(complement, known_tail)
    simple_source = [vector(ZZ, [0, 0] + list(root)) * preliminary for root in simple_tail.rows()]
    cartan = simple_tail * complement * simple_tail.transpose()
    assert abs(cartan.det()) == 36 and ZZ(pari(cartan).qfminim(2)[0]) == 24

    partial = u_rows.stack(matrix(ZZ, [list(-curve) for curve in simple_source]))
    smith, left, right = partial.smith_form()
    assert smith[:, :9] == matrix.identity(ZZ, 9) and not any(smith[:, 9:].list())
    basis = block_diagonal_matrix(left.inverse().change_ring(ZZ), matrix.identity(ZZ, 10)) * right.inverse().change_ring(ZZ)
    assert basis[:9] == partial
    for index in range(9, 19):
        row = basis.row(index)
        row -= ZZ(row * source_gram * (zero + fibre)) * fibre
        row -= ZZ(row * source_gram * fibre) * (zero + fibre)
        basis[index] = row
    assert abs(basis.det()) == 1
    inverse = basis.inverse().change_ring(ZZ)
    gram = basis * source_gram * basis.transpose()
    assert gram[:2, :2] == U2 and not any(gram[:2, 2:].list())
    frame = -gram[2:, 2:]
    assert frame[:7, :7] == cartan

    slug = zero_name.lower()
    frame_path = GENERATED / f"elkies-k3-h3-q4o208-corrected-a3-2a2-{slug}-frame.txt"
    marking_path = GENERATED / f"elkies-k3-h3-q4o208-corrected-a3-2a2-{slug}-marking.json"
    frame_path.write_text(f"# corrected q4/o323 physical A3+2A2; zero={zero_name}\n" +
                          "\n".join(" ".join(map(str, row)) for row in frame.rows()) + "\n")
    equation_to_source = matrix(ZZ, source["equation_A11_to_root_adapted_hub_basis"])
    payload = {
        "schema": "elkies-k3.h3-q4o208-corrected-a3-2a2-effective-zero-marking.v1",
        "status": "PASS_EXACT_CORRECTED_A3_2A2_PHYSICAL_EFFECTIVE_ZERO_MARKING",
        "hub": f"q4o208_corrected_A3_2A2_{zero_name}_zero",
        "zero": zero_name,
        "root_data": [7, 24, 36],
        "ade": "A3+2A2",
        "frame_output": str(frame_path.relative_to(ROOT)),
        "frame_sha256": sha256(frame_path),
        "physical_simple_components_in_source": [entries(curve) for curve in simple_source],
        "basis_in_source": rows(basis),
        "source_in_basis": rows(inverse),
        "equation_A11_to_root_adapted_hub_basis": rows(basis * equation_to_source),
        "target_fibres_in_root_adapted_hub": {
            name: entries(vector(ZZ, value) * inverse)
            for name, value in source["target_fibres_in_root_adapted_hub"].items()
        },
        "equation_explicit_curves_in_child": {name: entries(curve * inverse) for name, curve in explicit.items()},
        "prefix_operational_score": source.get("prefix_operational_score"),
        "proof_boundary": (
            "Exact physical effective zero and component chamber for the one-wall-corrected "
            "q4/o323 fibre, with all inherited curves, targets, and unimodular transports."
        ),
        "inputs": {"paths": [str(SOURCE.relative_to(ROOT)), str(CERT.relative_to(ROOT))],
                   "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (SOURCE, CERT)}},
    }
    marking_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    outputs[zero_name] = {"frame": str(frame_path.relative_to(ROOT)),
                          "marking": str(marking_path.relative_to(ROOT)),
                          "frame_max_abs": max(abs(int(x)) for x in basis.list())}

summary = {
    "schema": "elkies-k3.h3-corrected-a3-2a2-effective-zero-markings.v1",
    "status": "PASS_EXACT_CORRECTED_A3_2A2_ALL_PHYSICAL_EFFECTIVE_ZERO_MARKINGS",
    "outputs": outputs,
    "proof_boundary": "No ADE-only or pseudo-zero presentation is included.",
}
OUTPUT.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
print(f"CORRA3A2ZEROS|count={len(outputs)}|status={summary['status']}|output={OUTPUT}")
