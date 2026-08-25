#!/usr/bin/env sage -python
"""Export a root-adapted equation-A11 marking of the pinned semistable hub."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3/data/fibrations"
GENERATED = ROOT / "artifacts/generated-results"
FRAME = DATA / "mw2_a5_a4_a3a3_frame.txt"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
O12 = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
SUFFIX = GENERATED / "elkies-k3-h3-semistable-mw2-reverse-suffix-nef.json"
FRAME_OUTPUT = GENERATED / "elkies-k3-h3-semistable-mw2-root-adapted-frame.txt"
OUTPUT = GENERATED / "elkies-k3-h3-semistable-mw2-equation-marking.json"
INPUTS = (FRAME, CROSSOVERS, O12, SUFFIX)
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
    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == data[0]
    return simple, simple * gram * simple.transpose()


def root_adaptation(gram):
    unused, root_basis, invariants = roots_and_data(gram)
    root_rank = invariants[0]
    smith, smith_left, smith_right = root_basis.smith_form()
    assert smith == smith_left * root_basis * smith_right
    assert tuple(abs(smith[index, index]) for index in range(root_rank)) == (1,) * root_rank
    simple, cartan = deterministic_simple_roots(gram)
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[root_rank:])
    adapted = adapted_basis * gram * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    height = adapted[root_rank:, root_rank:] - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, root_rank), matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram()).transpose()
    )
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * gram * adapted_basis.transpose()
    assert abs(adapted_basis.det()) == 1
    return adapted, adapted_basis, invariants


frame = load_matrix(FRAME)
adapted, adapted_basis, root_data = root_adaptation(frame)
assert root_data == (15, 74, 480)
root_rank = root_data[0]
full_adaptation = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis)
g_frame = block_diagonal_matrix(U2, -frame)
g_adapted = block_diagonal_matrix(U2, -adapted)
assert full_adaptation * g_frame * full_adaptation.transpose() == g_adapted

crossovers = json.loads(CROSSOVERS.read_text())
o12 = json.loads(O12.read_text())
suffix = json.loads(SUFFIX.read_text())
pinned_in_equation = matrix(ZZ, crossovers["pinned_R17_basis_in_equation_A11"])
semistable_in_pinned = matrix(ZZ, suffix["composite_pinned_R17_to_semistable"])
semistable_in_equation = semistable_in_pinned * pinned_in_equation
adapted_in_equation = full_adaptation * semistable_in_equation
equation_in_adapted = adapted_in_equation.inverse().change_ring(ZZ)
assert abs(adapted_in_equation.det()) == 1

equation_fibre = vector(ZZ, [1, 0] + [0] * 17)
orbit12_in_equation = matrix(ZZ, o12["selected"]["equation_A11_to_explicit_zero_basis"])
targets_equation = {
    "equation_A11_fibre": equation_fibre,
    "orbit12_fibre": vector(ZZ, orbit12_in_equation.row(0)),
    "pinned_R17_fibre": vector(ZZ, pinned_in_equation.row(0)),
}
targets_adapted = {name: value * equation_in_adapted for name, value in targets_equation.items()}
assert all(value * g_adapted * value == 0 for value in targets_adapted.values())

FRAME_OUTPUT.write_text(
    "# root-adapted A5+A4+2A3/MW2 frame, equation-A11 marked\n"
    + "\n".join(" ".join(map(str, row)) for row in adapted.rows())
    + "\n"
)
payload = {
    "schema": "elkies-k3.h3-semistable-mw2-equation-marking.v1",
    "status": "PASS_EXACT_SEMISTABLE_MW2_EQUATION_MARKING",
    "root_data": list(map(int, root_data)),
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "standard_to_root_adapted_frame_basis": rows(adapted_basis),
    "root_adapted_to_standard_frame_basis": rows(adapted_basis.inverse().change_ring(ZZ)),
    "equation_A11_to_root_adapted_semistable_basis": rows(adapted_in_equation),
    "root_adapted_semistable_to_equation_A11_basis": rows(equation_in_adapted),
    "target_fibres_in_root_adapted_semistable": {name: entries(value) for name, value in targets_adapted.items()},
    "proof_boundary": "Exact root adaptation and full determinant-one equation-A11 marking. This artifact introduces no new neighbour or equation model.",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("SEMMARK|root=15,74,480|frame_max={}|basis_max={}|orbit12_degree={}|status={}".format(
    max(abs(value) for value in adapted.list()), max(abs(value) for value in adapted_in_equation.list()),
    targets_adapted["orbit12_fibre"][1], payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT}", flush=True)
