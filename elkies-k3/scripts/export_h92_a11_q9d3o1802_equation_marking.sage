#!/usr/bin/env sage-python
"""Export the certified compiler-cheap q9/degree-3/orbit1802 A11 child.

The existing lattice certificate already selects old A11 component 7 as its
zero and carries full determinant-one equation-A11 transports.  This replay
adds the equation-explicit physical A11 curves to the generic EC-beam marking
contract so its outgoing exact frontier can be priced by actual availability.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CERTIFICATE = GENERATED / "elkies-k3-h3-a11-q9d3o1802-lattice-certificate.json"
A11_MARKING = GENERATED / "elkies-k3-h3-equation-a11-marking.json"
DEFAULT_OUTPUT = GENERATED / "elkies-k3-h3-a11-q9d3o1802-equation-marking.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
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


certificate = json.loads(CERTIFICATE.read_text())
a11 = json.loads(A11_MARKING.read_text())
assert certificate["status"] == "PASS_EXACT_A11_DEGREE3_CANDIDATE_LATTICE_CERTIFICATE"
assert a11["status"] == "PASS_EXACT_A11_EQUATION_MARKING"

frame_path = ROOT / certificate["frame_output"]
a11_frame_path = ROOT / a11["frame_output"]
frame = load_matrix(frame_path)
a11_frame = load_matrix(a11_frame_path)
g_child = block_diagonal_matrix(U2, -frame)
g_a11 = block_diagonal_matrix(U2, -a11_frame)

# As verified by the Gram identity, rows of child_to_equation_A11_basis are
# equation-A11 basis vectors expressed in child coordinates.
a11_to_child = matrix(ZZ, certificate["transport"]["child_to_equation_A11_basis"])
child_to_a11 = matrix(ZZ, certificate["transport"]["equation_A11_to_child_basis"])
assert a11_to_child * g_child * a11_to_child.transpose() == g_a11
assert a11_to_child * child_to_a11 == identity_matrix(ZZ, 19)
assert abs(a11_to_child.det()) == 1

F = vector(ZZ, [1, 0] + [0] * 17)
O = vector(ZZ, [-1, 1] + [0] * 17)
simple = [
    vector(ZZ, [0, 0] + [-1 if index == node else 0 for index in range(17)])
    for node in range(11)
]
cartan = a11_frame[:11, :11]
half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
roots = tuple(half) + tuple(-root for root in half)
highest = max((root for root in roots if all(value >= 0 for value in root)), key=lambda root: sum(root))
affine = F + vector(ZZ, [0, 0] + list(highest) + [0] * 6)
explicit_a11 = {"old_A11_zero": O, "old_A11_affine": affine}
explicit_a11.update({f"old_A11_component_{index}": root for index, root in enumerate(simple)})
explicit_child = {name: curve * a11_to_child for name, curve in explicit_a11.items()}
assert all(curve * g_child * curve == -2 for curve in explicit_child.values())
assert explicit_child["old_A11_component_7"] == vector(ZZ, [-1, 1] + [0] * 17)

targets_child = {
    name: vector(ZZ, value) * a11_to_child
    for name, value in a11["target_fibres_in_root_adapted_hub"].items()
}
assert all(value * g_child * value == 0 for value in targets_child.values())
for name, degree in certificate["marked_target_degrees"].items():
    assert int(targets_child[name][1]) == int(degree)

source_degrees = {name: int(curve[1]) for name, curve in explicit_child.items()}
payload = {
    "schema": "elkies-k3.h3-a11-q9d3o1802-equation-marking.v1",
    "status": "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
    "source_hub": "a11_q9d3o1802_explicit_zero",
    "root_data": certificate["child"]["root_data"],
    "frame_output": str(frame_path.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest(),
    "equation_A11_to_child_basis": rows(a11_to_child),
    "child_to_equation_A11_basis": rows(child_to_a11),
    "target_fibres_in_child": {name: entries(value) for name, value in targets_child.items()},
    "equation_explicit_curves_in_child": {name: entries(value) for name, value in explicit_child.items()},
    "explicit_curve_source_degrees": source_degrees,
    "selected_zero": "old_A11_component_7",
    "prefix_operational_score": 914,
    "prefix_edges": [{"q": 9, "old_fibre_degree": 3, "orbit_index": 1802, "operational_score": 914}],
    "proof_boundary": (
        "Exact marked U with equation-explicit zero, root data, target fibres, and full "
        "determinant-one equation-A11 transports. Outgoing edge costs remain planning "
        "estimates and a selected continuation requires its own composed route certificate."
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
    "Q9D3O1802MARK|root={}|zero={}|deg0={}|deg1={}|orbit12_degree={}|prefix=914|det={}|status={}".format(
        ",".join(map(str, payload["root_data"])), payload["selected_zero"],
        sum(value == 0 for value in source_degrees.values()),
        sum(value == 1 for value in source_degrees.values()), targets_child["orbit12"][1],
        int(a11_to_child.det()), payload["status"],
    ), flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
