#!/usr/bin/env sage-python
"""Re-zero the certified A11/q8-o1991/q6-o84 child by its explicit section.

The q6/o84 fibre has degree one on ``old_A11_component_10``.  The historical
certificate uses an abstract Bezout mate instead, which makes physical curves
look artificially fibre-twisted.  This replay installs that already-explicit
curve as the zero, root-adapts the complement, and exports determinant-one
transports and all equation-A11 target fibres for continued EC search.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CERTIFICATE = GENERATED / "elkies-k3-h3-o1991-q6-orbit84-lattice-certificate.json"
A11_MARKING = GENERATED / "elkies-k3-h3-equation-a11-marking.json"
DEFAULT_FRAME = GENERATED / "elkies-k3-h3-o1991-q6o84-explicit-zero-frame.txt"
DEFAULT_OUTPUT = GENERATED / "elkies-k3-h3-o1991-q6o84-explicit-zero-marking.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--frame-output", type=Path, default=DEFAULT_FRAME)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
FRAME_OUTPUT = args.frame_output.resolve()
OUTPUT = args.output.resolve()
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-root for root in half])
    root_basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (root_basis.rank(), count, abs(ZZ(root_gram.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    positive = [root for root in roots if next(value for value in root if value != 0) > 0]
    positive_set = {tuple(root) for root in positive}
    simple = [root for root in positive if not any(tuple(root - left) in positive_set for left in positive)]
    result = matrix(ZZ, [list(root) for root in simple])
    assert result.nrows() == result.rank() == data[0]
    return result


def root_adaptation(gram):
    unused, root_basis, invariants = roots_and_data(gram)
    root_rank = int(invariants[0])
    smith, smith_left, smith_right = root_basis.smith_form()
    assert smith == smith_left * root_basis * smith_right
    assert tuple(abs(smith[index, index]) for index in range(root_rank)) == (1,) * root_rank
    simple = deterministic_simple_roots(gram)
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[root_rank:])
    adapted = adapted_basis * gram * adapted_basis.transpose()
    cartan = adapted[:root_rank, :root_rank]
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, root_rank),
        matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram()).transpose(),
    )
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * gram * adapted_basis.transpose()
    assert abs(adapted_basis.det()) == 1
    return adapted, adapted_basis, invariants


certificate = json.loads(CERTIFICATE.read_text())
a11 = json.loads(A11_MARKING.read_text())
assert certificate["status"] == "PASS_EXACT_O1991_Q6_ORBIT84_LATTICE_CERTIFICATE"
assert a11["status"] == "PASS_EXACT_A11_EQUATION_MARKING"

raw_frame_path = ROOT / certificate["child"]["frame_output"]
a11_frame_path = ROOT / a11["frame_output"]
raw_frame = load_matrix(raw_frame_path)
a11_frame = load_matrix(a11_frame_path)
g_raw = block_diagonal_matrix(U2, -raw_frame)
g_a11 = block_diagonal_matrix(U2, -a11_frame)

# Row-coordinate convention verified against both Gram matrices.
a11_to_raw = matrix(ZZ, certificate["transport"]["child_to_equation_A11_basis"])
raw_to_a11 = matrix(ZZ, certificate["transport"]["equation_A11_to_child_basis"])
assert a11_to_raw * g_raw * a11_to_raw.transpose() == g_a11
assert a11_to_raw * raw_to_a11 == identity_matrix(ZZ, 19)

F_a11 = vector(ZZ, [1, 0] + [0] * 17)
O_a11 = vector(ZZ, [-1, 1] + [0] * 17)
a11_simple = [
    vector(ZZ, [0, 0] + [-1 if index == node else 0 for index in range(17)])
    for node in range(11)
]
cartan = a11_frame[:11, :11]
half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
roots = tuple(half) + tuple(-root for root in half)
highest = max((root for root in roots if all(value >= 0 for value in root)), key=lambda root: sum(root))
affine_a11 = F_a11 + vector(ZZ, [0, 0] + list(highest) + [0] * 6)
explicit_a11 = {"old_A11_zero": O_a11, "old_A11_affine": affine_a11}
explicit_a11.update({f"old_A11_component_{index}": root for index, root in enumerate(a11_simple)})
explicit_raw = {name: curve * a11_to_raw for name, curve in explicit_a11.items()}
assert all(curve * g_raw * curve == -2 for curve in explicit_raw.values())

fibre_raw = vector(ZZ, [1, 0] + [0] * 17)
zero_raw = explicit_raw["old_A11_component_10"]
assert zero_raw * g_raw * fibre_raw == 1
mate = zero_raw + fibre_raw
complement = matrix(ZZ, [list(fibre_raw * g_raw), list(mate * g_raw)]).right_kernel_matrix()
split = matrix(ZZ, [list(fibre_raw), list(mate)] + list(complement.rows()))
assert abs(split.det()) == 1
rezero_frame = -(complement * g_raw * complement.transpose())
g_rezero = block_diagonal_matrix(U2, -rezero_frame)
assert split * g_raw * split.transpose() == g_rezero
split_inverse = split.inverse().change_ring(ZZ)

frame, adaptation, root_data = root_adaptation(rezero_frame)
full_adaptation = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation)
full_inverse = full_adaptation.inverse().change_ring(ZZ)
g_final = block_diagonal_matrix(U2, -frame)
assert full_adaptation * g_rezero * full_adaptation.transpose() == g_final

a11_to_final = a11_to_raw * split_inverse * full_inverse
final_to_a11 = full_adaptation * split * raw_to_a11
assert abs(a11_to_final.det()) == 1
assert a11_to_final * final_to_a11 == identity_matrix(ZZ, 19)
assert a11_to_final * g_final * a11_to_final.transpose() == g_a11

targets_final = {
    name: vector(ZZ, value) * a11_to_final
    for name, value in a11["target_fibres_in_root_adapted_hub"].items()
}
explicit_final = {name: curve * a11_to_final for name, curve in explicit_a11.items()}
assert explicit_final["old_A11_component_10"] == vector(ZZ, [-1, 1] + [0] * 17)
assert all(curve * g_final * curve == -2 for curve in explicit_final.values())
assert all(value * g_final * value == 0 for value in targets_final.values())

FRAME_OUTPUT.write_text(
    "# certified q6/o84 child re-zeroed by old_A11_component_10\n"
    + "\n".join(" ".join(map(str, row)) for row in frame.rows())
    + "\n"
)

source_degrees = {name: int(curve[1]) for name, curve in explicit_final.items()}
payload = {
    "schema": "elkies-k3.h3-o1991-q6o84-explicit-zero-marking.v1",
    "status": "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
    "source_hub": "a11_q8o1991_q6o84_explicit_zero",
    "root_data": list(map(int, root_data)),
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "equation_A11_to_child_basis": rows(a11_to_final),
    "child_to_equation_A11_basis": rows(final_to_a11),
    "target_fibres_in_child": {name: entries(value) for name, value in targets_final.items()},
    "equation_explicit_curves_in_child": {name: entries(value) for name, value in explicit_final.items()},
    "explicit_curve_source_degrees": source_degrees,
    "selected_zero": "old_A11_component_10",
    "prefix_operational_score": 2128,
    "prefix_edges": [
        {"q": 8, "orbit_index": 1991, "operational_score": 1449},
        {"q": 6, "orbit_index": 84, "operational_score": 679},
    ],
    "marked_U": {
        "fibre_in_raw_child": entries(fibre_raw),
        "zero_in_raw_child": entries(zero_raw),
        "mate_in_raw_child": entries(mate),
        "split_determinant": int(split.det()),
    },
    "proof_boundary": (
        "Exact explicit-zero U, primitive root complement, determinant-one transports in both "
        "directions, root data, equation-explicit curves, and all marked target fibres. New "
        "outgoing neighbours require exact component and Proposition-C2 nef replay. The prefix "
        "cost is an equation-planning estimate; no q6/o84 equation lift is claimed."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (CERTIFICATE, A11_MARKING, raw_frame_path, a11_frame_path)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (CERTIFICATE, A11_MARKING, raw_frame_path, a11_frame_path)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "O1991Q6O84ZERO|root={}|zero={}|vertical_explicit={}|sections={}|orbit12_degree={}|prefix=2128|det={}|status={}".format(
        ",".join(map(str, root_data)), payload["selected_zero"],
        sum(value == 0 for value in source_degrees.values()),
        sum(value == 1 for value in source_degrees.values()), targets_final["orbit12"][1],
        int(a11_to_final.det()), payload["status"],
    ), flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
