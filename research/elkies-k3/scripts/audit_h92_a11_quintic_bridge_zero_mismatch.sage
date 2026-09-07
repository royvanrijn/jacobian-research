#!/usr/bin/env sage -python
"""Audit the D12-zero compatibility of the proposed A11 quintic AJ bridge.

The equation-side orbit64 A11 shell is generated from the R3-zero D12 frame.
The proposed quintic bridge instead applies that transition after coordinates
from the A0-zero D12 frame.  This checker replays both unimodular changes and
records the resulting exact mismatch without modifying the active lift files.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json",
)
args = parser.parse_args()

SPIN = LOCAL / "q24-orbit42-spinor-zero-profiles.json"
Q6 = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
CURVES = LOCAL / "q24-downstream-lift/explicit-curves-a11-span-p100003.json"
BRIDGE = LOCAL / "q24-a11-target-coset-bridge.json"
FRAME = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
INPUTS = (SPIN, Q6, CURVES, BRIDGE, FRAME)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


spin = json.loads(SPIN.read_text())
q6 = json.loads(Q6.read_text())
curves = json.loads(CURVES.read_text())
bridge = json.loads(BRIDGE.read_text())
assert spin["status"] == "PASS_Q24_ORBIT42_EXACT_SPINOR_ZERO_PROFILES"
assert q6["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert bridge["status"] == "PASS_EXACT_A11_TARGET_COSET_BRIDGE"

profiles = {row["zero"]: row for row in spin["profiles"]}
assert set(profiles) == {"A0", "R3"}
selected_parent = load_matrix(FRAME)
assert matrix(ZZ, profiles["R3"]["frame"]) == selected_parent
assert matrix(ZZ, profiles["A0"]["frame"]) != selected_parent

record = next(row for row in q6["neighbors"] if int(row["orbit_index"]) == 64)
transition = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, record["child_root_adapted_basis"])
) * matrix(ZZ, record["neighbor_basis"])
child = matrix(ZZ, record["child_root_adapted_frame"])
g_child = block_diagonal_matrix(U2, -child)
assert abs(transition.det()) == 1
assert transition * block_diagonal_matrix(U2, -selected_parent) * transition.transpose() == g_child

by_name = {row["name"]: row for row in curves["explicit_curve_records"]}
source_curves = {
    name: vector(ZZ, by_name[name]["class"])
    for name in ("close_P24", "oldI9_A0")
}
assert all(int(by_name[name]["square"]) == -2 for name in source_curves)


def transported(profile_name, source):
    parent_basis = matrix(ZZ, profiles[profile_name]["parent_to_child_basis"])
    assert abs(parent_basis.det()) == 1
    return source * parent_basis.inverse().change_ring(ZZ) * transition.inverse().change_ring(ZZ)


wrong = {name: transported("A0", source) for name, source in source_curves.items()}
correct = {name: transported("R3", source) for name, source in source_curves.items()}
assert wrong["close_P24"] * g_child * wrong["close_P24"] == -3210
assert wrong["oldI9_A0"] * g_child * wrong["oldI9_A0"] == -2
assert all(value * g_child * value == -2 for value in correct.values())
assert (correct["close_P24"][1], correct["oldI9_A0"][1]) == (46, 4)
assert vector(ZZ, correct["close_P24"][-6:]) == vector(ZZ, (33, -77, 31, -38, 7, 1))
assert vector(ZZ, correct["oldI9_A0"][-6:]) == vector(ZZ, (2, -6, 3, -3, 1, 0))

shell = [
    vector(ZZ, value)
    for value in bridge["exact_identity_shell"]["MW_vectors_in_equation_order"]
]
coefficients = vector(
    ZZ,
    (-23, 21, -53, 0, 31, 0, 0, 0, -9, 0, 0, 0, 0, 0, 0, 0, 0, -28),
)
shell_sum = sum(
    (coefficient * point for coefficient, point in zip(coefficients, shell)),
    vector(ZZ, [0] * 6),
)
claimed_word_in_selected_marking = (
    -vector(ZZ, correct["close_P24"][-6:])
    + 3 * vector(ZZ, correct["oldI9_A0"][-6:])
    + shell_sum
)
target = vector(ZZ, bridge["selected_bridge"]["mw"])
assert claimed_word_in_selected_marking == vector(ZZ, (-23, 73, -40, 34, -4, -1))
assert claimed_word_in_selected_marking != target

payload = {
    "schema": "elkies-k3.h3-a11-quintic-bridge-zero-mismatch.v1",
    "status": "REJECT_A11_QUINTIC_BRIDGE_ZERO_MISMATCH",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "selected_equation_parent": {
        "zero": "R3",
        "frame": str(FRAME.relative_to(ROOT)),
        "q6_orbit": 64,
    },
    "mismatched_composition": {
        "zero_used_before_orbit64_transition": "A0",
        "close_P24_child_coordinates": entries(wrong["close_P24"]),
        "close_P24_square_in_selected_A11_NS": int(
            wrong["close_P24"] * g_child * wrong["close_P24"]
        ),
        "rejected_degree_five_coordinates_replayed": True,
    },
    "correct_selected_R3_transport": {
        name: {
            "child_coordinates": entries(value),
            "square": int(value * g_child * value),
            "A11_degree": int(value[1]),
            "A11_MW_Abel_Jacobi": entries(value[-6:]),
        }
        for name, value in correct.items()
    },
    "claimed_word_replay_in_selected_marking": {
        "claimed_target": entries(target),
        "actual_sum": entries(claimed_word_in_selected_marking),
        "matches": False,
    },
    "conclusion": (
        "The quintic bridge composes the orbit64 transition for the R3-zero D12 "
        "frame with coordinates taken in the distinct A0-zero D12 frame. In the "
        "selected equation marking close_P24 has degree 46, not 5, and the stored "
        "group word does not equal M. Do not use the quintic route until an exact "
        "A0-to-R3 frame isometry is inserted and the final selected marking is replayed."
    ),
    "proof_boundary": (
        "Exact integral NS replay of the two stored zero choices and the selected "
        "orbit64 transition. This rejects the present bridge artifact; it does not "
        "prove that no other low-degree Abel-Jacobi curve or corrected composition exists."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11QUINTICZERO|stored_square=-3210|selected_degree=46|"
    "selected_A0_degree=4|word_match=0|"
    f"status={payload['status']}",
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
