#!/usr/bin/env python3
"""Compact an exact neighbour shell by already-explicit curve nefness.

This inexpensive gate is intended to run before the Sage all-section CVP gate
on very large marked shells.  It uses only integral matrix arithmetic and does
not claim that a surviving fibre is nef: it proves only that none of the
declared, already-explicit (-2)-curves has negative degree.
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--neighbors", type=Path, required=True)
parser.add_argument("--marking", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()

neighbors_path = args.neighbors.resolve()
marking_path = args.marking.resolve()
output_path = args.output.resolve()
neighbors = json.loads(neighbors_path.read_text())
marking = json.loads(marking_path.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert marking["status"] in {
    "PASS_EXACT_A11_EQUATION_MARKING",
    "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
}

frame_path = ROOT / marking["frame_output"]
frame = [
    [int(value) for value in line.split()]
    for line in frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
]
assert len(frame) == 17 and all(len(row) == 17 for row in frame)
curves = {
    name: [int(value) for value in curve]
    for name, curve in marking["equation_explicit_curves_in_child"].items()
}


def pairing(left, right):
    value = left[0] * right[1] + left[1] * right[0]
    value -= sum(
        left[i + 2] * frame[i][j] * right[j + 2]
        for i in range(17) for j in range(17)
    )
    return value


assert all(pairing(curve, curve) == -2 for curve in curves.values())
survivors = []
rejected_by_q = Counter()
input_by_q = Counter()
for raw in neighbors["neighbors"]:
    q = int(raw["q"])
    input_by_q[q] += 1
    fibre = [int(value) for value in raw["fiber"]]
    degrees = {name: pairing(curve, fibre) for name, curve in curves.items()}
    negative = sorted(name for name, degree in degrees.items() if degree < 0)
    if negative:
        rejected_by_q[q] += 1
        continue
    survivors.append({
        "candidate_id": {
            "q": q,
            "old_fibre_degree": int(raw["old_fiber_degree"]),
            "orbit_index": int(raw["orbit_index"]),
        },
        "explicit_degree_zero_curves": sorted(
            name for name, degree in degrees.items() if degree == 0
        ),
        "explicit_degree_one_curves": sorted(
            name for name, degree in degrees.items() if degree == 1
        ),
    })


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


inputs = (neighbors_path, marking_path, frame_path)
payload = {
    "schema": "elkies-k3.h3-marked-explicit-curve-gate.v1",
    "status": "PASS_EXACT_MARKED_EXPLICIT_CURVE_GATE",
    "source_hub": marking.get("hub", marking.get("source_hub", "marked_hub")),
    "input_candidate_count": len(neighbors["neighbors"]),
    "survivor_count": len(survivors),
    "rejected_count": len(neighbors["neighbors"]) - len(survivors),
    "counts_by_q": {
        str(q): {
            "input": input_by_q[q],
            "rejected": rejected_by_q[q],
            "survived": input_by_q[q] - rejected_by_q[q],
        }
        for q in sorted(input_by_q)
    },
    "declared_explicit_curve_count": len(curves),
    "survivors": survivors,
    "proof_boundary": (
        "Exact integral nonnegative-intersection gate against every declared "
        "already-explicit (-2)-curve. Component/affine, all-section, root-adaptation, "
        "and endpoint-identification gates remain separate requirements."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
}
output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "EXPLICITCURVEGATE|hub={}|inputs={}|survivors={}|rejected={}|status={}|output={}".format(
        payload["source_hub"], payload["input_candidate_count"],
        payload["survivor_count"], payload["rejected_count"], payload["status"],
        output_path,
    ),
    flush=True,
)
