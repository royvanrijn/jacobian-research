#!/usr/bin/env sage -python
"""Root-adapt and equation-mark a certified reverse R17 hub."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
HUBS = {
    "q25_mw7": ROOT / "elkies-k3/data/fibrations/q25_mw7_frame.txt",
    "q25_mw4": ROOT / "elkies-k3/data/fibrations/q25_mw4_frame.txt",
    "mw3_a5_d4_2a2_a1": ROOT / "elkies-k3/data/fibrations/mw3_a5_d4_a2a2_a1_frame.txt",
}
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--hub", choices=sorted(HUBS), required=True)
parser.add_argument("--frame-output", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
FRAME = HUBS[args.hub]
FRAME_OUTPUT = args.frame_output.resolve()
OUTPUT = args.output.resolve()
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
O12 = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
INPUTS = (FRAME, CROSSOVERS, O12)
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
    return simple


def root_adaptation(gram):
    unused, root_basis, invariants = roots_and_data(gram)
    root_rank = invariants[0]
    smith, smith_left, smith_right = root_basis.smith_form()
    assert smith == smith_left * root_basis * smith_right
    assert tuple(abs(smith[index, index]) for index in range(root_rank)) == (1,) * root_rank
    simple = deterministic_simple_roots(gram)
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[root_rank:])
    adapted = adapted_basis * gram * adapted_basis.transpose()
    cartan = adapted[:root_rank, :root_rank]
    coupling = adapted[:root_rank, root_rank:]
    height = adapted[root_rank:, root_rank:] - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, root_rank),
        matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram()).transpose(),
    )
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * gram * adapted_basis.transpose()
    assert abs(adapted_basis.det()) == 1
    return adapted, adapted_basis, invariants


frame = load_matrix(FRAME)
adapted, adapted_basis, root_data = root_adaptation(frame)
full_adaptation = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis)
g_frame = block_diagonal_matrix(U2, -frame)
g_adapted = block_diagonal_matrix(U2, -adapted)
assert full_adaptation * g_frame * full_adaptation.transpose() == g_adapted

crossovers = json.loads(CROSSOVERS.read_text())
o12 = json.loads(O12.read_text())
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"
assert o12["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"
pinned_in_equation = matrix(ZZ, crossovers["pinned_R17_basis_in_equation_A11"])
hub_in_pinned = matrix(ZZ, crossovers["reverse_target_bases_in_pinned_R17"][args.hub])
hub_in_equation = hub_in_pinned * pinned_in_equation
adapted_in_equation = full_adaptation * hub_in_equation
equation_in_adapted = adapted_in_equation.inverse().change_ring(ZZ)
assert abs(adapted_in_equation.det()) == 1

o12_in_equation = matrix(ZZ, o12["selected"]["equation_A11_to_explicit_zero_basis"])
target_bases_equation = {
    "equation_A11": identity_matrix(ZZ, 19),
    "orbit12": o12_in_equation,
    "pinned_R17": pinned_in_equation,
}
for name, basis in crossovers["reverse_target_bases_in_pinned_R17"].items():
    target_bases_equation[name] = matrix(ZZ, basis) * pinned_in_equation
targets = {
    name: vector(ZZ, basis.row(0)) * equation_in_adapted
    for name, basis in target_bases_equation.items()
}
assert all(value in ZZ**19 and value * g_adapted * value == 0 for value in targets.values())

FRAME_OUTPUT.write_text(
    f"# root-adapted certified reverse hub {args.hub}, equation-A11 marked\n"
    + "\n".join(" ".join(map(str, row)) for row in adapted.rows())
    + "\n"
)
payload = {
    "schema": "elkies-k3.h3-reverse-hub-equation-marking.v1",
    "status": "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING",
    "hub": args.hub,
    "root_data": list(map(int, root_data)),
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "standard_to_root_adapted_frame_basis": rows(adapted_basis),
    "root_adapted_to_standard_frame_basis": rows(adapted_basis.inverse().change_ring(ZZ)),
    "equation_A11_to_root_adapted_hub_basis": rows(adapted_in_equation),
    "root_adapted_hub_to_equation_A11_basis": rows(equation_in_adapted),
    "target_fibres_in_root_adapted_hub": {name: entries(value) for name, value in targets.items()},
    "proof_boundary": (
        "Exact root adaptation and determinant-one equation-A11 marking of a hub "
        "already certified on the reverse suffix. New outgoing neighbours require "
        "their own nef and full transport certificates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "REVHUBMARK|hub={}|root={}|orbit12_degree={}|det={}|status={}".format(
        args.hub, ",".join(map(str, root_data)), targets["orbit12"][1],
        int(adapted_in_equation.det()), payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
