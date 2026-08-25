#!/usr/bin/env sage -python
"""Export the exact equation-A11 frame with pinned reverse-hub markings."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
NEIGHBORS = ROOT / "artifacts/local/elkies-k3/q24-a11-orbit64-q8-all.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
ORBIT12 = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
OUTPUT = GENERATED / "elkies-k3-h3-equation-a11-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


neighbors = json.loads(NEIGHBORS.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
orbit12 = json.loads(ORBIT12.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"
assert orbit12["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
frame_path = ROOT / neighbors["frame"]
frame = load_matrix(frame_path)
g = block_diagonal_matrix(U2, -frame)

targets = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"] if item["state"] == "equation_A11"
}
targets["equation_A11"] = vector(ZZ, [1, 0] + [0] * 17)
targets["orbit12"] = vector(ZZ, orbit12["edge"]["primitive_nef_isotropic_fibre"])
assert all(value * g * value == 0 for value in targets.values())
assert pari(frame).qfminim(2)[0] == 132

inputs = (NEIGHBORS, CROSSOVERS, ORBIT12, frame_path)
payload = {
    "schema": "elkies-k3.h3-equation-a11-marking.v1",
    "status": "PASS_EXACT_A11_EQUATION_MARKING",
    "hub": "equation_A11",
    "root_data": [11, 132, 12],
    "frame_output": str(frame_path.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest(),
    "equation_A11_to_root_adapted_hub_basis": [
        [int(row == column) for column in range(19)] for row in range(19)
    ],
    "root_adapted_hub_to_equation_A11_basis": [
        [int(row == column) for column in range(19)] for row in range(19)
    ],
    "target_fibres_in_root_adapted_hub": {
        name: list(map(int, value)) for name, value in targets.items()
    },
    "proof_boundary": (
        "Exact equation-A11 root-adapted frame and full marked fibres transported "
        "from the certified pinned crossover audit. This export proposes no new edge."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11MARK|targets={}|pinned_degree={}|orbit12_degree={}|status={}|output={}".format(
        len(targets), targets["pinned_R17"][1], targets["orbit12"][1],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
