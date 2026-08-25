#!/usr/bin/env sage -python
"""Export every current-route suffix frame and fibre in pinned-R17 coordinates."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/local/elkies-k3/q24-equation-d13-to-pinned-r17.json"
PINNED_MARKING = ROOT / "artifacts/generated-results/elkies-k3-h3-pinned-r17-equation-marking.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-pinned-r17-current-suffix-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


manifest = json.loads(MANIFEST.read_text())
pinned_marking = json.loads(PINNED_MARKING.read_text())
assert manifest["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert pinned_marking["status"] == "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING"
pinned = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned)

# The manifest's cumulative equation-D13/pinned matrix has rows of the pinned
# basis in equation-D13 coordinates.  Each step transition has rows of the
# child basis in its parent coordinates.
pinned_in_equation_d13 = matrix(
    ZZ, manifest["equation_d13_to_pinned_r17_transition"]
)
equation_d13_in_pinned = pinned_in_equation_d13.inverse().change_ring(ZZ)
cumulative = identity_matrix(ZZ, 19)
stage_names = (
    "current_D12", "current_A11", "current_A5A5", "current_3A3",
    "current_A3_2A2", "current_5A1", "current_4A1", "current_3A1",
    "current_2A1", "current_A1", "current_rootless",
)
stage_bases = {}
stage_data = {}
for name, step in zip(stage_names, manifest["steps"]):
    cumulative = matrix(ZZ, step["transition"]) * cumulative
    basis = cumulative * equation_d13_in_pinned
    assert abs(basis.det()) == 1
    gram = basis * g_pinned * basis.transpose()
    assert gram[:2, :2] == U2 and gram[:2, 2:] == 0
    stage_bases[name] = basis
    stage_data[name] = {
        "stage": step["stage"],
        "root_data": step["root_data"],
        "mw_rank": step["mw_rank"],
        "basis_in_pinned_R17": rows(basis),
        "pinned_R17_basis_in_stage": rows(basis.inverse().change_ring(ZZ)),
    }

targets = {
    name: vector(ZZ, value)
    for name, value in pinned_marking["target_fibres_in_root_adapted_hub"].items()
}
for name, basis in stage_bases.items():
    targets[name] = vector(ZZ, basis.row(0))
assert targets["current_A11"] == targets["equation_A11"]
assert targets["current_A5A5"] == targets["orbit12"]
assert targets["current_rootless"] == targets["pinned_R17"]
assert all(value * g_pinned * value == 0 for value in targets.values())

inputs = (MANIFEST, PINNED_MARKING, PINNED_FRAME)
payload = {
    "schema": "elkies-k3.h3-pinned-r17-current-suffix-marking.v1",
    "status": "PASS_EXACT_PINNED_R17_CURRENT_SUFFIX_MARKING",
    "hub": "pinned_R17",
    "root_data": [0, 0, 1],
    "frame_output": str(PINNED_FRAME.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(PINNED_FRAME.read_bytes()).hexdigest(),
    "equation_A11_to_root_adapted_hub_basis": pinned_marking[
        "equation_A11_to_root_adapted_hub_basis"
    ],
    "root_adapted_hub_to_equation_A11_basis": pinned_marking[
        "root_adapted_hub_to_equation_A11_basis"
    ],
    "target_fibres_in_root_adapted_hub": {
        name: list(map(int, value)) for name, value in targets.items()
    },
    "current_suffix_stages": stage_data,
    "exact_overlap_checks": {
        "current_A11_equals_equation_A11": True,
        "current_A5A5_equals_orbit12": True,
        "current_rootless_equals_pinned_R17": True,
    },
    "proof_boundary": (
        "Exact full determinant-one transport of every current-route suffix frame "
        "into pinned R17. This is a marking export and asserts no alternative edge."
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
    "R17SUFFIXMARK|stages={}|targets={}|A11={}|A5A5={}|rootless={}|status={}|output={}".format(
        len(stage_data), len(targets), int(targets["current_A11"] == targets["equation_A11"]),
        int(targets["current_A5A5"] == targets["orbit12"]),
        int(targets["current_rootless"] == targets["pinned_R17"]),
        payload["status"], OUTPUT,
    ),
    flush=True,
)
