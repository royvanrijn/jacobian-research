#!/usr/bin/env sage -python
"""Export an explicit-zero A5A5-exit child with exact marked route targets.

The lattice certificate and the explicit equation zero use different U
splittings on the same child. This exporter deliberately uses the latter: it
composes the pinned suffix from the historical manifest and reverse-hub targets
from equation-A11 coordinates through the selected explicit-zero basis.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
SUFFIX_CROSSOVERS = GENERATED / "elkies-k3-h3-candidate-current-suffix-crossovers.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--certificate", type=Path, required=True)
parser.add_argument("--explicit-zero-frame", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
CERTIFICATE = args.certificate.resolve()
ZERO_FRAME = args.explicit_zero_frame.resolve()
OUTPUT = args.output.resolve()

certificate = json.loads(CERTIFICATE.read_text())
zero_payload = json.loads(ZERO_FRAME.read_text())
manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
suffix_crossovers = json.loads(SUFFIX_CROSSOVERS.read_text())
assert certificate["status"] in {
    "PASS_EXACT_A5A5_EXPLICIT_ZERO_CANDIDATE_LATTICE_CERTIFICATE",
    "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE",
    "PASS_EXACT_Q6O1307_PHYSICAL_WEYL_REPAIR_REJECT_C10_ZERO",
}
physical_nef_repair = (
    certificate["status"]
    == "PASS_EXACT_Q6O1307_PHYSICAL_WEYL_REPAIR_REJECT_C10_ZERO"
)
assert zero_payload["status"] == "PASS_EXACT_CANDIDATE_A11_COMPONENT_EXPLICIT_ZERO_FRAMES"
assert manifest["status"] == "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"
assert suffix_crossovers["status"] == "PASS_EXACT_CANDIDATE_CURRENT_SUFFIX_CROSSOVER_AUDIT"

selected = zero_payload["selected"]
candidate = certificate.get("selection", {}).get("candidate_id")
if candidate is None:
    candidate = certificate["candidate_id"]
candidate_name = (
    "q{}_orbit{}".format(candidate["q"], candidate["orbit_index"])
    if "orbit_index" in candidate
    else candidate["label"]
)
assert selected["root_data"] == certificate["child"]["root_data"]
frame = matrix(ZZ, selected["frame"])
g_candidate = block_diagonal_matrix(U2, -frame)
equation_to_candidate = matrix(ZZ, selected["equation_A11_to_explicit_zero_basis"])
candidate_to_equation = equation_to_candidate.inverse().change_ring(ZZ)
assert abs(equation_to_candidate.det()) == 1

# Full suffix bases in the equation-A11 coordinates used by the explicit-zero
# construction. The fingerprint's historical basis is already the A11 parent,
# so the pre-A11 manifest steps must not enter this composition.
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
suffix_bases_candidate = {
    name: basis * candidate_to_equation
    for name, basis in suffix_bases_equation.items()
}
assert all(abs(basis.det()) == 1 for basis in suffix_bases_candidate.values())

# Named reverse targets were independently certified in this equation-A11 basis.
targets_equation = {
    record["target"]: vector(ZZ, record["target_fibre_in_state"])
    for record in crossovers["records"]
    if record["state"] == "equation_A11"
}
targets_candidate = {
    name: value * candidate_to_equation for name, value in targets_equation.items()
}
explicit_curves_equation = {
    **{
        "old_A11_component_{}".format(index):
        vector(ZZ, [0, 0] + [-ZZ(other == index) for other in range(17)])
        for index in range(11)
    },
    "old_A11_zero": vector(ZZ, [-1, 1] + [0] * 17),
}
explicit_curves_candidate = {
    name: value * candidate_to_equation
    for name, value in explicit_curves_equation.items()
}
assert all(curve * g_candidate * curve == -2 for curve in explicit_curves_candidate.values())
for name, basis in suffix_bases_candidate.items():
    targets_candidate[name] = vector(ZZ, basis.row(0))

# Regression against the separate direct-crossover audit detects mixing the
# abstract certificate zero with this equation-explicit zero.
if "orbit_index" in candidate and not physical_nef_repair:
    for record in suffix_crossovers["records"]:
        if record["candidate"] != candidate_name:
            continue
        key = stage_key_by_label[record["target_stage"]]
        assert targets_candidate[key] == vector(ZZ, record["target_fibre_in_candidate"])
assert targets_candidate["pinned_R17"] == vector(ZZ, suffix_bases_candidate["pinned_R17"].row(0))
assert all(value * g_candidate * value == 0 for value in targets_candidate.values())

frame_path = ROOT / zero_payload["selected_frame_output"]
assert frame == load_matrix(frame_path)
inputs = (CERTIFICATE, ZERO_FRAME, MANIFEST, FINGERPRINT, CROSSOVERS, SUFFIX_CROSSOVERS, frame_path)
payload = {
    "schema": "elkies-k3.h3-a5a5-candidate-marking.v2",
    "status": "PASS_EXACT_A5A5_CANDIDATE_SUFFIX_MARKING",
    "hub": "a5a5_{}_explicit_zero".format(candidate_name),
    "candidate_id": candidate,
    "root_data": selected["root_data"],
    "frame_output": str(frame_path.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest(),
    "equation_A11_to_root_adapted_hub_basis": rows(equation_to_candidate),
    "root_adapted_hub_to_equation_A11_basis": rows(candidate_to_equation),
    "target_fibres_in_root_adapted_hub": {
        name: list(map(int, value)) for name, value in targets_candidate.items()
    },
    "equation_explicit_curves_in_child": {
        name: list(map(int, value)) for name, value in explicit_curves_candidate.items()
    },
    "current_suffix_stage_bases_in_root_adapted_hub": {
        name: rows(basis) for name, basis in suffix_bases_candidate.items()
    },
    "proof_boundary": (
        "Exact explicit-equation-zero marking, the old A11 components and zero, all "
        "named reverse target fibres, and full determinant-one current-suffix bases. "
        "The two old-I12 affine components are not included in the explicit-curve list; "
        "no continuation edge is asserted."
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
    "A5CANDMARK|candidate={}|targets={}|pinned_degree={}|status={}|output={}".format(
        candidate_name, len(targets_candidate), targets_candidate["pinned_R17"][1],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
