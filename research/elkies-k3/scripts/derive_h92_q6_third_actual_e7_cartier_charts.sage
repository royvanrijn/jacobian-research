#!/usr/bin/env sage -python
"""Attach the third divisor's integral E7 factors to actual H92 edge charts."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, vector

ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-local-target.json"
RESOLUTION = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-resolution-full.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-actual-e7-cartier-charts.json"
TARGET_SHA256 = "a0699e4ec75930cc93a9706ddf96f4ffc744954809e53a24029bd8c6668843f7"
RESOLUTION_SHA256 = "14378f4718d3fbe781d5b351ba4943a962a430d9827c0ce285cd1125a9e8c500"

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--resolution", type=Path, default=RESOLUTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

if args.target == TARGET:
    assert digest(args.target) == TARGET_SHA256
if args.resolution == RESOLUTION:
    assert digest(args.resolution) == RESOLUTION_SHA256
target = json.loads(args.target.read_text())
resolution = json.loads(args.resolution.read_text())
assert target["status"] == "PASS_EXACT_Q6_THIRD_E7_LATTICE_TARGET"
assert resolution["status"] == "PASS_EXACT_H92_E7_BLOWUP_TREE"
orders = -vector(ZZ, target["resolved_exceptional_coefficients"])
assert orders == vector(ZZ, (22, 44, 66, 44, 33, 33, 55))

edge_data = (
    ("E7_2--E7_5", "E7_5", "Z", "E7_2", "Y"),
    ("E7_1--E7_4", "E7_4", "U", "E7_1", "Y"),
    ("E7_4--E7_3", "E7_4", "Z", "E7_3", "Y"),
    ("E7_3--E7_7", "E7_7", "U", "E7_3", "Y"),
    ("E7_7--E7_2", "E7_7", "Z", "E7_2", "Y"),
    ("E7_3--E7_6", "E7_6", "Z", "E7_3", "Y"),
)
assert set(resolution["edge_charts"]) == {entry[0] for entry in edge_data}
charts = []
for name, exceptional, exceptional_variable, old, old_variable in edge_data:
    exceptional_order = int(orders[int(exceptional[-1]) - 1])
    old_order = int(orders[int(old[-1]) - 1])
    factor = "{}^{}*{}^{}".format(exceptional_variable, exceptional_order, old_variable, old_order)
    charts.append({
        "name": name,
        "actual_h92_surface_equation": resolution["edge_charts"][name],
        "exceptional_component": exceptional,
        "exceptional_equation": exceptional_variable,
        "strict_component": old,
        "strict_equation": old_variable,
        "cartier_factor": factor,
        "membership_condition": "regular representative belongs to ({})".format(factor),
    })

payload = {
    "schema": "elkies-k3.h92-q6-third-actual-e7-cartier-charts.v1",
    "status": "PASS_EXACT_H92_Q6_THIRD_E7_CARTIER_CHARTS",
    "inputs": {
        "lattice_target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)},
        "actual_resolution": {"path": str(args.resolution.relative_to(ROOT)), "sha256": digest(args.resolution)},
    },
    "resolved_vanishing_orders": [int(value) for value in orders],
    "charts": charts,
    "compiler_instruction": "Trivialize O(V_E7) by the inverse displayed factor in each actual chart, then evaluate the marked-chord DAG and reduce in bounded local quotients.",
    "boundary": "This is the actual H92 E7 integral-vertical input; it does not yet evaluate the marked chord or choose finite quotient jets.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92ACTUALE7CARTIER|charts=6|orders=22,44,66,44,33,33,55|status=PASS_EXACT_H92_Q6_THIRD_E7_CARTIER_CHARTS", flush=True)
