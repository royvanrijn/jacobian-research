#!/usr/bin/env sage -python
"""Certify the marked E7 chart for a nested source-q8 endpoint enlargement.

The actual marked E7 frame permits a coefficient u^i/h^k on x^a*m^b exactly
when i<=4k+d(a,b), where d is the certified marked-frame denominator power.
This script checks that the r=4 smooth-candidate ambient retains that bound.
It deliberately says nothing about the five other E7 edge charts.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ENDPOINT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
MARKED = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-actual-e7-marked-frame.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-extra4-marked-e7-cover.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--extra-h-power", type=int, default=4)
parser.add_argument("--extra-e7-pole", type=int, default=0)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.extra_h_power < 0 or args.extra_e7_pole < 0:
    raise ValueError("enlargement parameters must be nonnegative")

endpoint = json.loads(ENDPOINT.read_text())
marked = json.loads(MARKED.read_text())
assert endpoint["status"] == "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT"
assert marked["status"] == "PASS_EXACT_Q8_ACTUAL_E7_MARKED_FRAME"

marked_poles = {
    (entry["basis"]["kind"], int(entry["basis"]["x_power"]), int(entry["basis"]["m_power"])):
    int(entry["t_denominator_power"])
    for entry in marked["normalized_q8_generators"]
}
ambient = []
for family in endpoint["families"]:
    key = (family["generic_basis"]["kind"], int(family["generic_basis"]["x_power"]), int(family["generic_basis"]["m_power"]))
    k = int(family["h_power"]) + args.extra_h_power
    upper = 4*k + marked_poles[key] + args.extra_e7_pole
    for i in range(int(family["e8_minimal_u_power"]), upper + 1):
        ambient.append({"basis": key, "u_power": i, "h_power": k, "marked_e7_upper": 4*k + marked_poles[key]})

if args.extra_e7_pole:
    conclusion = "extra E7 slack was requested, so the marked E7 bound is not certified for all displayed generators"
    certified = False
else:
    assert all(item["u_power"] <= item["marked_e7_upper"] for item in ambient)
    conclusion = "every enlarged endpoint generator satisfies the exact actual marked-E7 frame inequality"
    certified = True

payload = {
    "schema": "elkies-k3.h92-q8-enlarged-endpoint-marked-e7-cover.v1",
    "status": "PASS_EXACT_Q8_ENLARGED_MARKED_E7_COVER" if certified else "DIAGNOSTIC_E7_SLACK_NOT_COVERED",
    "inputs": {"endpoint": {"path": str(ENDPOINT.relative_to(ROOT)), "sha256": digest(ENDPOINT)}, "marked_frame": {"path": str(MARKED.relative_to(ROOT)), "sha256": digest(MARKED)}},
    "parameters": {"extra_h_power": args.extra_h_power, "extra_e7_pole": args.extra_e7_pole},
    "ambient_dimension": len(ambient),
    "marked_chart": {"chart": marked["actual_chart"]["name"], "condition": "u_power<=4*h_power+marked_t_denominator_power", "certified": certified},
    "boundary": "Only the marked E7 chart is covered. Other E7 edges, E8 compatibility of any enlarged ambient, characteristic-zero smooth kernel reconstruction, and all later equation/bisection work remain open.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q8EXTRAENDPOINTE7|extra_h={}|ambient={}|marked_chart={}|status={}".format(args.extra_h_power, len(ambient), int(certified), payload["status"]), flush=True)
