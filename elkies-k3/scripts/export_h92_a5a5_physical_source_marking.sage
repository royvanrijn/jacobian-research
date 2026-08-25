#!/usr/bin/env sage -python
"""Export the equation-effective component-9-zero 2A5 source marking."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
ZERO = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-physical-source-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


zero = json.loads(ZERO.read_text())
manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
assert zero["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"
assert zero["selected_zero_curve"] == "old_A11_component_9"
assert manifest["status"] == "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"

selected = zero["selected"]
frame_path = ROOT / zero["selected_frame_output"]
frame = load_matrix(frame_path)
assert frame == matrix(ZZ, selected["frame"])
g_source = block_diagonal_matrix(U2, -frame)
equation_to_source = matrix(ZZ, selected["equation_A11_to_explicit_zero_basis"])
source_to_equation = equation_to_source.inverse().change_ring(ZZ)
assert abs(equation_to_source.det()) == 1

historical_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
stage_key_by_label = {
    "2A5/MW7": "current_A5A5",
    "3A3/MW8": "current_3A3",
    "A3+2A2/MW10": "current_A3_2A2",
    "5A1/MW12": "current_5A1",
    "4A1/MW13": "current_4A1",
    "3A1/MW14": "current_3A1",
    "2A1/MW15": "current_2A1",
    "A1/MW16": "current_A1",
    "rootless/MW17": "pinned_R17",
}
cumulative = identity_matrix(ZZ, 19)
suffix_bases_equation = {"current_A11": historical_in_equation}
for index, step in enumerate(manifest["forward_steps"]):
    if index < 2:
        continue
    cumulative = matrix(ZZ, step["transition"]) * cumulative
    suffix_bases_equation[stage_key_by_label[step["child"]]] = cumulative * historical_in_equation
suffix_bases_source = {
    name: basis * source_to_equation for name, basis in suffix_bases_equation.items()
}
assert all(abs(basis.det()) == 1 for basis in suffix_bases_source.values())

targets_equation = {
    record["target"]: vector(ZZ, record["target_fibre_in_state"])
    for record in crossovers["records"]
    if record["state"] == "equation_A11"
}
targets_source = {
    name: value * source_to_equation for name, value in targets_equation.items()
}
for name, basis in suffix_bases_source.items():
    targets_source[name] = vector(ZZ, basis.row(0))
assert all(value * g_source * value == 0 for value in targets_source.values())
assert targets_source["current_A5A5"] == vector(ZZ, [1, 0] + [0] * 17)

inputs = (ZERO, MANIFEST, FINGERPRINT, CROSSOVERS, frame_path)
payload = {
    "schema": "elkies-k3.h3-a5a5-physical-source-marking.v1",
    "status": "PASS_EXACT_A5A5_PHYSICAL_SOURCE_MARKING",
    "hub": "equation_effective_component9_zero_2A5",
    "root_data": selected["root_data"],
    "frame_output": str(frame_path.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest(),
    "equation_A11_to_root_adapted_hub_basis": rows(equation_to_source),
    "root_adapted_hub_to_equation_A11_basis": rows(source_to_equation),
    "target_fibres_in_root_adapted_hub": {
        name: list(map(int, value)) for name, value in targets_source.items()
    },
    "current_suffix_stage_bases_in_root_adapted_hub": {
        name: rows(basis) for name, basis in suffix_bases_source.items()
    },
    "proof_boundary": (
        "Exact equation-effective component-9-zero 2A5 frame, full equation-A11 "
        "transport, named reverse targets, and current-suffix bases. No outgoing "
        "neighbour is asserted."
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
    "A5PHYSICAL|targets={}|current3A3_degree={}|pinned_degree={}|status={}|output={}".format(
        len(targets_source), targets_source["current_3A3"][1],
        targets_source["pinned_R17"][1], payload["status"], OUTPUT,
    ),
    flush=True,
)
