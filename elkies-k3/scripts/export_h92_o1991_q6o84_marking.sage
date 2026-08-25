#!/usr/bin/env sage-python
"""Export the certified A11/q8-o1991/q6-o84 child as a marked beam state.

The q6/o84 certificate already contains determinant-one transports between
the root-adapted ``A1+A2+D10/MW4`` child and equation-A11 coordinates.  This
replay turns those transports into the standard marking contract consumed by
``rank_h92_marked_root_adapted_frontier.sage`` and also transports the
equation-explicit A11 fibre components for later equation-cost scoring.
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
DEFAULT_OUTPUT = GENERATED / "elkies-k3-h3-o1991-q6o84-equation-marking.json"
DEFAULT_FRAME_OUTPUT = GENERATED / "elkies-k3-h3-o1991-q6o84-physical-frame.txt"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--frame-output", type=Path, default=DEFAULT_FRAME_OUTPUT)
args = parser.parse_args()
OUTPUT = args.output.resolve()
FRAME_OUTPUT = args.frame_output.resolve()

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


certificate = json.loads(CERTIFICATE.read_text())
a11 = json.loads(A11_MARKING.read_text())
assert certificate["status"] == "PASS_EXACT_O1991_Q6_ORBIT84_LATTICE_CERTIFICATE"
assert a11["status"] == "PASS_EXACT_A11_EQUATION_MARKING"

frame_path = ROOT / certificate["child"]["frame_output"]
a11_frame_path = ROOT / a11["frame_output"]
raw_frame = load_matrix(frame_path)
a11_frame = load_matrix(a11_frame_path)
g_raw_child = block_diagonal_matrix(U2, -raw_frame)
g_a11 = block_diagonal_matrix(U2, -a11_frame)

# The certificate field names describe the basis identity, while its matrices
# act on row coordinates in the opposite direction: rows of ``child_to_*``
# are the equation-A11 basis expressed in child coordinates.
a11_to_child = matrix(ZZ, certificate["transport"]["child_to_equation_A11_basis"])
child_to_a11 = matrix(ZZ, certificate["transport"]["equation_A11_to_child_basis"])
assert child_to_a11 * a11_to_child == identity_matrix(ZZ, 19)
assert a11_to_child * child_to_a11 == identity_matrix(ZZ, 19)
assert abs(child_to_a11.det()) == 1
assert a11_to_child * g_raw_child * a11_to_child.transpose() == g_a11

targets_a11 = {
    name: vector(ZZ, value)
    for name, value in a11["target_fibres_in_root_adapted_hub"].items()
}
targets_raw_child = {
    name: value * a11_to_child
    for name, value in targets_a11.items()
}
assert all(value in ZZ**19 and value * g_raw_child * value == 0 for value in targets_raw_child.values())

# In the equation-A11 root-adapted frame the first eleven positive coordinates
# are the non-identity A11 components.  The affine component is F plus the
# A11 highest root, and O is the usual (-1,1,0) class for U=<F,O+F>.
F = vector(ZZ, [1, 0] + [0] * 17)
O = vector(ZZ, [-1, 1] + [0] * 17)
simple = [vector(ZZ, [0, 0] + [-1 if i == j else 0 for i in range(17)]) for j in range(11)]
cartan = a11_frame[:11, :11]
half_roots = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
roots = tuple(half_roots) + tuple(-root for root in half_roots)
highest = max((root for root in roots if all(value >= 0 for value in root)), key=lambda root: sum(root))
affine = F + vector(ZZ, [0, 0] + list(highest) + [0] * 6)
explicit_a11 = {"old_A11_zero": O, "old_A11_affine": affine}
explicit_a11.update({f"old_A11_component_{index}": root for index, root in enumerate(simple)})
explicit_raw_child = {name: value * a11_to_child for name, value in explicit_a11.items()}
assert all(value in ZZ**19 and value * g_raw_child * value == -2 for value in explicit_raw_child.values())

# The certificate's deterministic root basis is not the physical chamber: in
# it every enumerated dominant neighbour is negative on two transported old
# A11 components.  Reorient the finite A1+A2+D10 root system to a chamber in
# which every equation-explicit vertical curve is effective.  Positive roots
# in the positive-definite frame are the negatives of actual vertical curve
# classes because NS uses ``-frame``.
root_rank = int(certificate["child"]["root_data"][0])
raw_root = raw_frame[:root_rank, :root_rank]
physical_positive = [
    -vector(ZZ, value[2:2 + root_rank])
    for value in explicit_raw_child.values()
    if value[1] == 0
]
assert physical_positive and all(root * raw_root * root == 2 for root in physical_positive)
cone = Polyhedron(
    ieqs=[[-1] + list(raw_root * root) for root in physical_positive],
    base_ring=QQ,
)
assert not cone.is_empty()
chamber_vector = vector(QQ, cone.an_element())
half = [vector(ZZ, column) for column in matrix(ZZ, pari(raw_root).qfminim(2)[2]).columns()]
all_roots = tuple(half + [-root for root in half])
perturbation = vector(QQ, range(1, root_rank + 1))
for denominator in range(1001, 1000000, 1000):
    trial = chamber_vector + perturbation / denominator
    if (
        min(trial * raw_root * root for root in physical_positive) > 0
        and all(trial * raw_root * root != 0 for root in all_roots)
    ):
        chamber_vector = trial
        break
else:
    raise ArithmeticError("failed to choose a regular physical root chamber")
positive_roots = [root for root in all_roots if chamber_vector * raw_root * root > 0]
positive_set = {tuple(root) for root in positive_roots}
physical_simple = [
    root for root in positive_roots
    if not any(tuple(root - other) in positive_set for other in positive_roots)
]
root_change = matrix(ZZ, [list(root) for root in physical_simple])
assert root_change.nrows() == root_change.rank() == root_rank
assert abs(root_change.det()) == 1
assert all(
    all(value >= 0 for value in root * root_change.inverse())
    for root in physical_positive
)
positive_change = block_diagonal_matrix(root_change, identity_matrix(ZZ, 17 - root_rank))
full_change = block_diagonal_matrix(identity_matrix(ZZ, 2), positive_change)
full_inverse = full_change.inverse().change_ring(ZZ)
frame = positive_change * raw_frame * positive_change.transpose()
g_child = block_diagonal_matrix(U2, -frame)
targets_child = {name: value * full_inverse for name, value in targets_raw_child.items()}
explicit_child = {name: value * full_inverse for name, value in explicit_raw_child.items()}
a11_to_child = a11_to_child * full_inverse
child_to_a11 = full_change * child_to_a11
assert a11_to_child * g_child * a11_to_child.transpose() == g_a11
assert child_to_a11 * a11_to_child == identity_matrix(ZZ, 19)
assert all(value * g_child * value == 0 for value in targets_child.values())
assert all(value * g_child * value == -2 for value in explicit_child.values())
FRAME_OUTPUT.write_text(
    "# physical root chamber for certified A11/q8-o1991/q6-o84 child\n"
    + "\n".join(" ".join(map(str, row)) for row in frame.rows())
    + "\n"
)

payload = {
    "schema": "elkies-k3.h3-o1991-q6o84-equation-marking.v1",
    "status": "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
    "source_hub": "a11_q8o1991_q6o84",
    "root_data": certificate["child"]["root_data"],
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "raw_certificate_frame": str(frame_path.relative_to(ROOT)),
    "physical_chamber_basis_in_raw_child": rows(full_change),
    "raw_child_basis_in_physical_chamber": rows(full_inverse),
    "child_to_equation_A11_basis": rows(child_to_a11),
    "equation_A11_to_child_basis": rows(a11_to_child),
    "target_fibres_in_child": {name: entries(value) for name, value in targets_child.items()},
    "equation_explicit_curves_in_child": {name: entries(value) for name, value in explicit_child.items()},
    "prefix_operational_score": 2128,
    "prefix_edges": [
        {"q": 8, "orbit_index": 1991, "operational_score": 1449},
        {"q": 6, "orbit_index": 84, "operational_score": 679},
    ],
    "proof_boundary": (
        "Exact marked U, full determinant-one equation-A11 transports, target fibres, "
        "and transported equation-explicit A11 curves for the certified q8/o1991 then "
        "q6/o84 branch. New outgoing neighbours require their own exact nef and route "
        "certificates; the 2128 prefix is an equation-cost estimate, not an equation lift."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (CERTIFICATE, A11_MARKING, frame_path, a11_frame_path)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (CERTIFICATE, A11_MARKING, frame_path, a11_frame_path)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "O1991Q6O84MARK|root={}|orbit12_degree={}|pinned_degree={}|prefix=2128|det={}|status={}".format(
        ",".join(map(str, payload["root_data"])), targets_child["orbit12"][1],
        targets_child["pinned_R17"][1], int(child_to_a11.det()), payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
