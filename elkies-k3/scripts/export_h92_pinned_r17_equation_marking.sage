#!/usr/bin/env sage -python
"""Export pinned R17 with exact equation-A11 and reverse-hub fibre markings."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
ORBIT12 = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
OUTPUT = GENERATED / "elkies-k3-h3-pinned-r17-equation-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


crossovers = json.loads(CROSSOVERS.read_text())
orbit12 = json.loads(ORBIT12.read_text())
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"
assert orbit12["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"

frame = load_matrix(PINNED_FRAME)
g = block_diagonal_matrix(U2, -frame)
assert frame.dimensions() == (17, 17) and frame.det() == 948
assert pari(frame).qfminim(2)[0] == 0

# Rows of P are the pinned basis expressed in equation-A11 coordinates.
pinned_in_equation = matrix(ZZ, crossovers["pinned_R17_basis_in_equation_A11"])
equation_in_pinned = pinned_in_equation.inverse().change_ring(ZZ)
assert abs(pinned_in_equation.det()) == 1
targets = {
    "pinned_R17": vector(ZZ, [1, 0] + [0] * 17),
    "equation_A11": vector(ZZ, equation_in_pinned.row(0)),
}
orbit12_in_equation = vector(ZZ, orbit12["edge"]["primitive_nef_isotropic_fibre"])
targets["orbit12"] = orbit12_in_equation * equation_in_pinned
for name, basis_rows in crossovers["reverse_target_bases_in_pinned_R17"].items():
    targets[name] = vector(ZZ, matrix(ZZ, basis_rows).row(0))
assert all(value in ZZ**19 and value * g * value == 0 for value in targets.values())
assert targets["pinned_R17"] == vector(ZZ, [1, 0] + [0] * 17)

payload = {
    "schema": "elkies-k3.h3-pinned-r17-equation-marking.v1",
    "status": "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING",
    "hub": "pinned_R17",
    "root_data": [0, 0, 1],
    "frame_output": str(PINNED_FRAME.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(PINNED_FRAME.read_bytes()).hexdigest(),
    "equation_A11_to_root_adapted_hub_basis": rows(pinned_in_equation),
    "root_adapted_hub_to_equation_A11_basis": rows(equation_in_pinned),
    "target_fibres_in_root_adapted_hub": {
        name: list(map(int, value)) for name, value in targets.items()
    },
    "proof_boundary": (
        "Exact pinned-rootless marking export. All target fibres are transported by "
        "the certified full unimodular pinned-R17/equation-A11 basis map. This export "
        "does not assert that a new reverse neighbor exists or is nef."
    ),
    "inputs": {
        "paths": [str(CROSSOVERS.relative_to(ROOT)), str(ORBIT12.relative_to(ROOT)),
                  str(PINNED_FRAME.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (CROSSOVERS, ORBIT12, PINNED_FRAME)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "R17MARK|targets={}|equation_A11_degree={}|orbit12_degree={}|det={}|status={}|output={}".format(
        len(targets), targets["equation_A11"][1], targets["orbit12"][1],
        pinned_in_equation.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
