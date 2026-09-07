#!/usr/bin/env sage -python
"""Identify the pinned target-directed q8 child with the current-route A1 frame."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "artifacts/local/elkies-k3/q24-equation-d13-to-pinned-r17.json"
CERTIFICATE = ROOT / "artifacts/generated-results/elkies-k3-h3-pinned-r17-q8-orbit12-cvp-lattice-certificate.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-pinned-r17-q8-current-route-a1-identification.json"
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
certificate = json.loads(CERTIFICATE.read_text())
assert manifest["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert certificate["status"] == "PASS_EXACT_PINNED_R17_TARGETED_CANDIDATE_CERTIFICATE"
pinned = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned)
a1 = matrix(ZZ, manifest["final_a1_frame"])
g_a1 = block_diagonal_matrix(U2, -a1)

# The stored manifest matrix has rows of the pinned basis in the A1 basis;
# invert it to put both A1 bases in pinned coordinates.
current_a1_in_pinned = matrix(
    ZZ, manifest["final_a1_to_pinned_r17_transition"]
).inverse().change_ring(ZZ)
selected_a1_in_pinned = matrix(ZZ, certificate["source_to_child_basis"])
assert current_a1_in_pinned * g_pinned * current_a1_in_pinned.transpose() == g_a1
selected_frame = load_matrix(ROOT / certificate["frame_output"])
g_selected = block_diagonal_matrix(U2, -selected_frame)
assert selected_a1_in_pinned * g_pinned * selected_a1_in_pinned.transpose() == g_selected
assert current_a1_in_pinned.row(0) == selected_a1_in_pinned.row(0)
assert vector(ZZ, [1, 0] + [0] * 17) == vector(
    ZZ, certificate["target_fibres_in_child"]["pinned_R17"]
) * selected_a1_in_pinned

current_to_selected = current_a1_in_pinned * selected_a1_in_pinned.inverse()
assert current_to_selected in MatrixSpace(ZZ, 19)
assert abs(current_to_selected.det()) == 1
assert current_to_selected * g_selected * current_to_selected.transpose() == g_a1
assert current_to_selected.row(0) == vector(ZZ, [1, 0] + [0] * 17)
simple_component = vector(ZZ, [0, 0, 1] + [0] * 16)
assert simple_component * current_to_selected == simple_component

inputs = (MANIFEST, CERTIFICATE, PINNED_FRAME)
payload = {
    "schema": "elkies-k3.h3-pinned-r17-q8-current-route-a1-identification.v1",
    "status": "PASS_EXACT_PINNED_Q8_EQUALS_CURRENT_ROUTE_A1_FIBRE",
    "pinned_q8_fibre": list(map(int, selected_a1_in_pinned.row(0))),
    "current_route_stage": "A1/MW16 before q6 orbit2247 rootless",
    "root_data": [1, 2, 2],
    "current_a1_basis_in_pinned_R17": rows(current_a1_in_pinned),
    "selected_a1_basis_in_pinned_R17": rows(selected_a1_in_pinned),
    "current_a1_basis_in_selected_a1": rows(current_to_selected),
    "selected_a1_basis_in_current_a1": rows(
        current_to_selected.inverse().change_ring(ZZ)
    ),
    "full_unimodular_marking_identification": True,
    "simple_and_affine_A1_component_chamber_aligned": True,
    "edge_profiles": {
        "pinned_R17_to_A1": {
            "q": certificate["candidate_id"]["q"],
            "old_fibre_degree": certificate["candidate_id"]["old_fibre_degree"],
            "P_dot_O": certificate["candidate_id"]["q"] // 2 - 2,
        },
        "A1_to_pinned_R17": {
            "q": certificate["target_profiles"]["pinned_R17"]["q"],
            "old_fibre_degree": certificate["target_profiles"]["pinned_R17"]["old_fibre_degree"],
            "P_dot_O": certificate["target_profiles"]["pinned_R17"]["P_dot_O"],
            "nef_audit": certificate["target_profiles"]["pinned_R17"]["nef_audit"],
        },
    },
    "conclusion": (
        "The target-directed pinned q8 child is the existing current-route A1 fibre, "
        "with a full integral marking identification. Subsequent target-directed states "
        "are lateral detours from the current suffix, not a new pinned approach."
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
    "PINNEDQ8A1|same_fibre=1|det={}|status={}|output={}".format(
        current_to_selected.det(), payload["status"], OUTPUT
    ),
    flush=True,
)
