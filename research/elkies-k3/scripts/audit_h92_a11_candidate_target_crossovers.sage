#!/usr/bin/env sage -python
"""Audit marked direct fibre degrees from A11 route candidates to R17 hubs."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
U2 = matrix(ZZ, ((0, 1), (1, 0)))
INPUTS = (
    LOCAL / "q24-equation-d13-to-pinned-r17.json",
    LOCAL / "h3-r17-backward-exact-lift-manifest.json",
    LOCAL / "q24-a11-q8-construction-fingerprint.json",
    GENERATED / "elkies-k3-h3-a11-q8-orbit1991-lattice-certificate.json",
    GENERATED / "elkies-k3-h3-a11-q8-orbit1991-explicit-zero-frames.json",
    GENERATED / "elkies-k3-h3-o1991-q6-orbit84-lattice-certificate.json",
    GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json",
    GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4-orbit32-lattice-certificate.json",
    GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q6-orbit3372-lattice-certificate.json",
    GENERATED / "elkies-k3-h3-semistable-mw2-pinned-transport.json",
    GENERATED / "elkies-k3-h3-semistable-mw2-reverse-suffix-nef.json",
)
OUTPUT = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    if current == -1:
        coefficients = [-value for value in coefficients]
    return vector(ZZ, coefficients)


def reconstruct_step(parent, a, b, witness):
    ns_parent = block_diagonal_matrix(U2, -parent)
    fibre = vector(ZZ, [a, b] + list(witness))
    assert fibre * ns_parent * fibre == 0
    mate = bezout_vector(ns_parent * fibre)
    mate -= ZZ(mate * ns_parent * mate // 2) * fibre
    assert fibre * ns_parent * mate == 1 and mate * ns_parent * mate == 0
    orthogonal = matrix(
        ZZ, [list(fibre * ns_parent), list(mate * ns_parent)]
    ).right_kernel_matrix()
    child = -(orthogonal * ns_parent * orthogonal.transpose())
    transition = matrix(ZZ, [list(fibre), list(mate)] + orthogonal.rows())
    assert abs(transition.det()) == 1
    return child, transition


closeout = json.loads(INPUTS[0].read_text())
manifest = json.loads(INPUTS[1].read_text())
fingerprint = json.loads(INPUTS[2].read_text())
o1991 = json.loads(INPUTS[3].read_text())
explicit = json.loads(INPUTS[4].read_text())
o84 = json.loads(INPUTS[5].read_text())
o12_explicit = json.loads(INPUTS[6].read_text())
o32 = json.loads(INPUTS[7].read_text())
o3372 = json.loads(INPUTS[8].read_text())
semistable = json.loads(INPUTS[9].read_text())
semistable_nef = json.loads(INPUTS[10].read_text())

equation_d13_to_d12 = matrix(
    ZZ, closeout["q24"]["equation_d13_to_d12_transition"]
)
equation_d13_to_pinned = matrix(
    ZZ, closeout["equation_d13_to_pinned_r17_transition"]
)
pinned_in_equation_d12 = equation_d13_to_pinned * equation_d13_to_d12.inverse()
historical_d12_to_a11 = matrix(ZZ, manifest["forward_steps"][1]["transition"])
pinned_in_historical_a11 = pinned_in_equation_d12 * historical_d12_to_a11.inverse()
historical_in_equation_frame = matrix(
    ZZ,
    fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"],
)
historical_in_equation_ns = block_diagonal_matrix(identity_matrix(ZZ, 2), historical_in_equation_frame)
pinned_in_equation_a11 = pinned_in_historical_a11 * historical_in_equation_ns
assert abs(pinned_in_equation_a11.det()) == 1

rank17_frame = load_matrix(ROOT / "elkies-k3/data/lattice/rank17_gram.txt")
pinned_frame = load_matrix(ROOT / closeout["pinned_rank17_frame"])
assert rank17_frame == pinned_frame

reverse_specs = (
    ("q25_mw7", ROOT / "elkies-k3/data/fibrations/q25_mw7_frame.txt", 5, 5,
     (-1, 0, -4, 3, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, -3, 0, 0)),
    ("q25_mw4", ROOT / "elkies-k3/data/fibrations/q25_mw4_frame.txt", 2, 2,
     (-1, -2, 1, 0, 1, 1, 2, -3, 0, -2, 0, 1, 0, 0, -1, 0, 0)),
    ("mw3_a5_d4_2a2_a1", ROOT / "elkies-k3/data/fibrations/mw3_a5_d4_a2a2_a1_frame.txt", 2, 2,
     (-1, 0, 0, 2, 0, 2, -1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0)),
    ("mw2_e6_d4_2a2_a1", ROOT / "elkies-k3/data/fibrations/mw2_e6_d4_a2a2_a1_frame.txt", 2, 2,
     (0, -2, -2, 0, 1, 0, 3, 2, 4, 0, 3, 0, -1, -3, -1, -4, 0)),
)

target_bases_in_pinned = {"pinned_R17": identity_matrix(ZZ, 19)}
parent_frame = rank17_frame
cumulative = identity_matrix(ZZ, 19)
for name, expected_path, a, b, witness in reverse_specs:
    child_frame, transition = reconstruct_step(parent_frame, a, b, witness)
    assert child_frame == load_matrix(expected_path)
    cumulative = transition * cumulative
    target_bases_in_pinned[name] = cumulative
    parent_frame = child_frame
target_bases_in_pinned["mw2_a5_a4_2a3_semistable"] = matrix(
    ZZ, semistable_nef["composite_pinned_R17_to_semistable"]
)

states = {
    "equation_A11": identity_matrix(ZZ, 19),
    "q8_orbit1991_abstract_zero": matrix(
        ZZ, o1991["transport"]["parent_to_child_basis"]
    ),
    "q8_orbit1991_explicit_zero": matrix(
        ZZ, explicit["selected"]["equation_A11_to_explicit_zero_basis"]
    ),
    "q6_orbit84": matrix(
        ZZ, o84["transport"]["equation_A11_to_child_basis"]
    ),
    "q8_orbit12_explicit_zero": matrix(
        ZZ, o12_explicit["selected"]["equation_A11_to_explicit_zero_basis"]
    ),
    "q4_orbit32": matrix(
        ZZ, o32["transport"]["equation_A11_to_child_basis"]
    ),
    "q6_orbit3372": matrix(
        ZZ, o3372["transport"]["equation_A11_to_child_basis"]
    ),
}

records = []
for state_name, state_in_equation in states.items():
    assert abs(state_in_equation.det()) == 1
    equation_in_state = state_in_equation.inverse().change_ring(ZZ)
    for target_name, target_in_pinned in target_bases_in_pinned.items():
        target_in_equation = target_in_pinned * pinned_in_equation_a11
        target_in_state = target_in_equation * equation_in_state
        fibre = vector(ZZ, target_in_state.row(0))
        records.append({
            "state": state_name,
            "target": target_name,
            "target_fibre_in_state": entries(fibre),
            "direct_old_fibre_degree": int(fibre[1]),
            "fibre_max_abs_coordinate": int(max(abs(value) for value in fibre)),
            "basis_max_abs_coordinate": int(max(abs(value) for value in target_in_state.list())),
        })

payload = {
    "schema": "elkies-k3.h3-a11-candidate-target-crossovers.v1",
    "status": "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "pinned_R17_basis_in_equation_A11": rows(pinned_in_equation_a11),
    "reverse_target_bases_in_pinned_R17": {
        name: rows(value) for name, value in target_bases_in_pinned.items()
    },
    "records": records,
    "conclusion": (
        "These are exact marked direct fibre degrees, not certified one-edge nef "
        "neighbors. They are a search heuristic and do not promote a route."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for item in records:
    print(
        "A11CROSSOVER|state={state}|target={target}|degree={direct_old_fibre_degree}|"
        "fibre_max={fibre_max_abs_coordinate}|basis_max={basis_max_abs_coordinate}".format(**item),
        flush=True,
    )
print(f"A11CROSSOVER|status={payload['status']}|output={OUTPUT.resolve()}", flush=True)
