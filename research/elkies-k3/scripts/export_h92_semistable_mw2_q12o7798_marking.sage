#!/usr/bin/env sage -python
"""Export the exact equation marking at the closest q12 semistable child.

This is a search checkpoint, not a promoted route certificate.  The selected
first hop has already passed the exact affine and all-section nef gates in the
semistable frontier scorer; here we compose its two unimodular basis changes
and retain the full equation-A11 marking for a second bidirectional layer.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--q", type=int, default=12)
parser.add_argument("--orbit", type=int, default=7798)
parser.add_argument(
    "--frame-output", type=Path,
    default=GENERATED / "elkies-k3-h3-semistable-mw2-q12o7798-root-adapted-frame.txt",
)
parser.add_argument(
    "--output", type=Path,
    default=GENERATED / "elkies-k3-h3-semistable-mw2-q12o7798-equation-marking.json",
)
args = parser.parse_args()
SEM_MARKING = GENERATED / "elkies-k3-h3-semistable-mw2-equation-marking.json"
SEM_SCORE = GENERATED / "elkies-k3-h3-semistable-mw2-q8q10q12-equation-cost.json"
FRAME_OUTPUT = args.frame_output.resolve()
OUTPUT = args.output.resolve()
INPUTS = (SEM_MARKING, SEM_SCORE)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


sem_marking = json.loads(SEM_MARKING.read_text())
sem_score = json.loads(SEM_SCORE.read_text())
assert sem_marking["status"] == "PASS_EXACT_SEMISTABLE_MW2_EQUATION_MARKING"
assert sem_score["status"] == "PASS_EXACT_SEMISTABLE_MW2_FRONTIER_COST_SCORING"

matches = [
    item for item in sem_score["ranked_candidates"]
    if item["candidate_id"]["q"] == args.q and item["candidate_id"]["orbit_index"] == args.orbit
]
assert len(matches) == 1
selected = matches[0]
raw = selected["source_neighbor_record"]

sem_frame_path = ROOT / sem_marking["frame_output"]
sem_frame = matrix(
    ZZ,
    [[ZZ(value) for value in line.split()] for line in sem_frame_path.read_text().splitlines()
     if line.strip() and not line.lstrip().startswith("#")],
)
g_sem = block_diagonal_matrix(U2, -sem_frame)
neighbor_basis = matrix(ZZ, raw["neighbor_basis"])
child_adaptation = matrix(ZZ, raw["child_root_adapted_basis"])
child_frame = matrix(ZZ, raw["child_root_adapted_frame"])
transition = block_diagonal_matrix(identity_matrix(ZZ, 2), child_adaptation) * neighbor_basis
g_child = block_diagonal_matrix(U2, -child_frame)
assert abs(neighbor_basis.det()) == abs(child_adaptation.det()) == abs(transition.det()) == 1
assert transition * g_sem * transition.transpose() == g_child

sem_in_equation = matrix(ZZ, sem_marking["equation_A11_to_root_adapted_semistable_basis"])
child_in_equation = transition * sem_in_equation
equation_in_child = child_in_equation.inverse().change_ring(ZZ)
assert abs(child_in_equation.det()) == 1
targets_sem = {
    name: vector(ZZ, value)
    for name, value in sem_marking["target_fibres_in_root_adapted_semistable"].items()
}
targets_child = {name: value * transition.inverse() for name, value in targets_sem.items()}
assert all(value in ZZ**19 and value * g_child * value == 0 for value in targets_child.values())
new_fibre_sem = vector(ZZ, transition.row(0))
assert new_fibre_sem == vector(ZZ, raw["fiber"])
assert {
    name: int(new_fibre_sem * g_sem * value) for name, value in targets_sem.items()
} == selected["marked_target_degrees"]

FRAME_OUTPUT.write_text(
    "# root-adapted A1+A2+A5+D5/MW4 frame at semistable q12 orbit7798\n"
    + "\n".join(" ".join(map(str, row)) for row in child_frame.rows())
    + "\n"
)
payload = {
    "schema": "elkies-k3.h3-semistable-mw2-frontier-equation-marking.v1",
    "status": "PASS_EXACT_SEMISTABLE_FRONTIER_EQUATION_MARKING_SEARCH_CHECKPOINT",
    "candidate_id": selected["candidate_id"],
    "child": selected["child"],
    "first_hop_equation_cost_score": selected["equation_cost_score"],
    "first_hop_marked_target_degrees": selected["marked_target_degrees"],
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "root_adapted_semistable_to_child_basis": rows(transition),
    "child_to_root_adapted_semistable_basis": rows(transition.inverse().change_ring(ZZ)),
    "equation_A11_to_child_basis": rows(child_in_equation),
    "child_to_equation_A11_basis": rows(equation_in_child),
    "target_fibres_in_child": {name: entries(value) for name, value in targets_child.items()},
    "proof_boundary": (
        "Exact determinant-one marking composition for a two-step search checkpoint. "
        "The first hop's nef gates are inherited from the scored frontier artifact; "
        "this file alone is not a promoted complete-route certificate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(sem_frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (sem_frame_path,)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "SEM2MARK|q={}|orbit={}|root={}|orbit12={}|det={}|status={}".format(
        args.q,
        args.orbit,
        ",".join(map(str, selected["child"]["root_data"])),
        selected["marked_target_degrees"]["orbit12_fibre"],
        int(transition.det()),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
