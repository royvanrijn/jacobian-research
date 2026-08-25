#!/usr/bin/env sage -python
"""Rank arbitrary-degree A11 neighbour shells by exact marked hub distances."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--neighbors", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
NEIGHBORS = args.neighbors.resolve()
OUTPUT = args.output.resolve()
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
INPUTS = (NEIGHBORS, CROSSOVERS)
U2 = matrix(ZZ, ((0, 1), (1, 0)))
TARGETS = (
    "pinned_R17",
    "q25_mw7",
    "q25_mw4",
    "mw3_a5_d4_2a2_a1",
    "mw2_e6_d4_2a2_a1",
    "mw2_a5_a4_2a3_semistable",
)


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


neighbors = json.loads(NEIGHBORS.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"
frame_path = ROOT / neighbors["frame"]
frame = load_matrix(frame_path)
g = block_diagonal_matrix(U2, -frame)
targets = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"] if item["state"] == "equation_A11"
}
assert set(targets) == set(TARGETS)

records = []
identity_component_rejected = 0
for raw in neighbors["neighbors"]:
    degree = int(raw["old_fiber_degree"])
    labels = list(map(int, raw["dominant_labels"]))
    # The old A11 affine component has highest-root coefficients all one.
    identity_degree = degree - sum(labels)
    if identity_degree < 0:
        identity_component_rejected += 1
        continue
    fibre = vector(ZZ, raw["fiber"])
    marked = {name: int(fibre * g * target) for name, target in targets.items()}
    assert min(marked.values()) >= 0
    records.append({
        "candidate_id": {
            "q": int(raw["q"]),
            "old_fibre_degree": degree,
            "orbit_index": int(raw["orbit_index"]),
        },
        "child": {
            "ade": raw["child_ade"],
            "mw_rank": int(raw["child_mw_rank"]),
            "root_data": raw["child_root_data"],
        },
        "component_degrees": labels + [identity_degree],
        "marked_target_degrees": marked,
    })

rankings = {
    target: sorted(
        records,
        key=lambda item: (
            item["marked_target_degrees"][target],
            item["candidate_id"]["q"],
            item["candidate_id"]["orbit_index"],
        ),
    )[:100]
    for target in TARGETS
}
payload = {
    "schema": "elkies-k3.h3-a11-general-degree-marked-target-ranking.v1",
    "status": "PASS_EXACT_A11_GENERAL_DEGREE_MARKED_TARGET_RANKING",
    "input_candidate_count": len(neighbors["neighbors"]),
    "affine_nef_candidate_count": len(records),
    "identity_component_rejected_count": identity_component_rejected,
    "rankings_top_100": rankings,
    "proof_boundary": (
        "Exact marked degrees after the full old-A11 component gate. This fast ranking "
        "does not include the all-section nef gate, equation-cost profile, or a full "
        "transport certificate for any selected child."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (frame_path,)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for target in TARGETS:
    best = rankings[target][0]
    print(
        "A11GEN|target={}|q={}|d={}|orbit={}|degree={}|child={}/MW{}".format(
            target, best["candidate_id"]["q"], best["candidate_id"]["old_fibre_degree"],
            best["candidate_id"]["orbit_index"], best["marked_target_degrees"][target],
            best["child"]["ade"], best["child"]["mw_rank"],
        ),
        flush=True,
    )
print(
    "A11GEN|inputs={}|affine_nef={}|status={}|output={}".format(
        len(neighbors["neighbors"]), len(records), payload["status"], OUTPUT
    ),
    flush=True,
)
