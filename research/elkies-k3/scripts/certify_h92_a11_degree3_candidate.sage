#!/usr/bin/env sage -python
"""Certify and export a selected full-nef equation-side A11 degree-3 exit."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--q", type=int, default=9)
parser.add_argument("--orbit", type=int, default=1802)
parser.add_argument(
    "--frame-output", type=Path,
    default=GENERATED / "elkies-k3-h3-a11-q9d3o1802-root-adapted-frame.txt",
)
parser.add_argument(
    "--output", type=Path,
    default=GENERATED / "elkies-k3-h3-a11-q9d3o1802-lattice-certificate.json",
)
args = parser.parse_args()
NEIGHBORS = GENERATED / "elkies-k3-h3-a11-q6q9q12-degree3-all.json"
GATE = GENERATED / "elkies-k3-h3-a11-q6q9q12-degree3-full-nef-equation-cost.json"
FRAME_OUTPUT = args.frame_output.resolve()
OUTPUT = args.output.resolve()
INPUTS = (NEIGHBORS, GATE)
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


neighbors = json.loads(NEIGHBORS.read_text())
gate = json.loads(GATE.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert gate["status"] == "PASS_EXACT_A11_DEGREE3_FULL_NEF_EQUATION_COST_GATE"
candidate_pool = (
    gate["ranked_candidates_top_200"]
    + gate["closest_to_orbit12_top_100"]
    + gate["closest_to_pinned_top_100"]
    + gate["orbit12_degree_equation_cost_pareto_front"]
)
matches = [
    item for item in candidate_pool
    if item["candidate_id"]["q"] == args.q and item["candidate_id"]["orbit_index"] == args.orbit
]
assert matches
selected = matches[0]
assert all(item == selected for item in matches)
assert selected["nef_audit"]["nef"]
raw = selected["source_neighbor_record"]

parent_frame_path = ROOT / neighbors["frame"]
parent_frame = load_matrix(parent_frame_path)
g_parent = block_diagonal_matrix(U2, -parent_frame)
neighbor_basis = matrix(ZZ, raw["neighbor_basis"])
child_adaptation = matrix(ZZ, raw["child_root_adapted_basis"])
child_frame = matrix(ZZ, raw["child_root_adapted_frame"])
transition = block_diagonal_matrix(identity_matrix(ZZ, 2), child_adaptation) * neighbor_basis
transition_inverse = transition.inverse().change_ring(ZZ)
g_child = block_diagonal_matrix(U2, -child_frame)
assert abs(neighbor_basis.det()) == abs(child_adaptation.det()) == abs(transition.det()) == 1
assert transition * g_parent * transition.transpose() == g_child
assert transition_inverse * g_child * transition_inverse.transpose() == g_parent
fibre = vector(ZZ, transition.row(0))
zero = vector(ZZ, transition.row(1))
assert fibre == vector(ZZ, raw["fiber"])
assert fibre * g_parent * fibre == zero * g_parent * zero == 0
assert fibre * g_parent * zero == 1

# Recover target vectors from the gate's prerequisites using the raw target
# fields stored in the A11 crossover and orbit12 artifacts.
crossovers_path = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
o12_path = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
crossovers = json.loads(crossovers_path.read_text())
o12 = json.loads(o12_path.read_text())
targets_parent = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"] if item["state"] == "equation_A11"
}
targets_parent["orbit12"] = vector(
    ZZ, matrix(ZZ, o12["selected"]["equation_A11_to_explicit_zero_basis"]).row(0)
)
assert {
    name: int(fibre * g_parent * target) for name, target in targets_parent.items()
} == selected["marked_target_degrees"]
targets_child = {name: value * transition_inverse for name, value in targets_parent.items()}
assert all(value in ZZ**19 and value * g_child * value == 0 for value in targets_child.values())

FRAME_OUTPUT.write_text(
    "# root-adapted child of equation A11 q{} degree3 orbit{}\n".format(args.q, args.orbit)
    + "\n".join(" ".join(map(str, row)) for row in child_frame.rows())
    + "\n"
)
payload = {
    "schema": "elkies-k3.h3-a11-degree3-candidate-lattice-certificate.v1",
    "status": "PASS_EXACT_A11_DEGREE3_CANDIDATE_LATTICE_CERTIFICATE",
    "candidate_id": selected["candidate_id"],
    "child": selected["child"],
    "equation_cost_score": selected["equation_cost_score"],
    "horizontal": selected["horizontal"],
    "expected_RR_ambient": selected["expected_RR_ambient"],
    "explicit_curve_degrees": selected["explicit_curve_degrees"],
    "target_coset_mod_exact_sections": selected["target_coset_mod_exact_sections"],
    "nef_audit": selected["nef_audit"],
    "marked_target_degrees": selected["marked_target_degrees"],
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "transport": {
        "equation_A11_to_child_basis": rows(transition),
        "child_to_equation_A11_basis": rows(transition_inverse),
        "forward_determinant": int(transition.det()),
        "inverse_determinant": int(transition_inverse.det()),
        "marked_U_in_equation_A11": [entries(fibre), entries(zero)],
    },
    "target_fibres_in_child": {name: entries(value) for name, value in targets_child.items()},
    "proof_boundary": (
        "Exact marked U, primitive nef isotropic class, full child root data, and "
        "unimodular NS transports in both directions. The child is not yet connected "
        "to pinned R17 by a certified continuation."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS]
        + [str(parent_frame_path.relative_to(ROOT)), str(crossovers_path.relative_to(ROOT)), str(o12_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (parent_frame_path, crossovers_path, o12_path)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11D3CERT|q={}|orbit={}|child={}/MW{}|det={}|orbit12={}|pinned={}|cost={}|PO={}|RR={}|status={}".format(
        args.q, args.orbit, selected["child"]["ade"], selected["child"]["mw_rank"],
        int(transition.det()), selected["marked_target_degrees"]["orbit12"],
        selected["marked_target_degrees"]["pinned_R17"], selected["equation_cost_score"],
        selected["horizontal"]["P_dot_O"], selected["expected_RR_ambient"], payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
