#!/usr/bin/env sage -python
"""Carry equation-explicit curves through the physical q4/o323 correction."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
CERT = GENERATED / "elkies-k3-h3-q4o208-physical-q4o323-corrected-a3-2a2-certificate.json"
OUTPUT = GENERATED / "elkies-k3-h3-q4o208-physical-q4o323-corrected-a3-2a2-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


source = json.loads(SOURCE.read_text())
cert = json.loads(CERT.read_text())
assert source["status"] == "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING"
assert cert["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert cert["candidate_id"]["label"] == "q4o323-physical-wall-corrected"
assert cert["first_edge_exact_horizontal_nef_gate"]

source_frame = load_matrix(ROOT / source["frame_output"])
child_frame = load_matrix(ROOT / cert["frame_output"])
source_gram = block_diagonal_matrix(U2, -source_frame)
child_gram = block_diagonal_matrix(U2, -child_frame)
child_in_source = matrix(ZZ, cert["source_to_child_basis"])
source_in_child = matrix(ZZ, cert["child_to_source_basis"])
assert child_in_source * source_gram * child_in_source.transpose() == child_gram
assert child_in_source * source_in_child == matrix.identity(ZZ, 19)

explicit = {
    name: vector(ZZ, value)
    for name, value in source["equation_explicit_curves_in_child"].items()
}
explicit_child = {
    name: entries(curve * source_in_child)
    for name, curve in explicit.items()
}
assert all(
    vector(ZZ, value) * child_gram * vector(ZZ, value) == -2
    for value in explicit_child.values()
)

inputs = (SOURCE, CERT)
payload = {
    "schema": "elkies-k3.h3-q4o208-corrected-a3-2a2-marking.v1",
    "status": "PASS_EXACT_A5A5_CANDIDATE_SUFFIX_MARKING",
    "hub": "q4o208_physical_q4o323_corrected_A3_2A2",
    "root_data": cert["child"]["root_data"],
    "frame_output": cert["frame_output"],
    "frame_sha256": cert["frame_sha256"],
    "target_fibres_in_child": cert["target_fibres_in_child"],
    "equation_A11_to_child_basis": cert["equation_A11_to_child_basis"],
    "equation_explicit_curves_in_child": explicit_child,
    "source_to_child_basis": cert["source_to_child_basis"],
    "child_to_source_basis": cert["child_to_source_basis"],
    "prefix_operational_score": None,
    "proof_boundary": (
        "Exact inherited-curve and marked-target transport through the fully certified "
        "physical q4/o323 wall-corrected edge.  No successor is promoted without its own "
        "physical curve, component, all-section and finite-horizontal-wall gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"Q4O323PHYSMARK|root={payload['root_data']}|curves={len(explicit_child)}|output={OUTPUT}")
