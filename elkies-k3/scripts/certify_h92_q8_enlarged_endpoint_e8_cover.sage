#!/usr/bin/env sage -python
"""Certify E8 membership for a nested source-q8 endpoint enlargement.

Increasing an h-denominator does not alter E8 order because h(0) is a unit.
This checks that every r=4 candidate generator retains its certified E8
minimal u-power.  It is only the E8 end of the local cover.
"""

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENDPOINT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-endpoint-rr-ambient.json"
E8 = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e8-ambient-weights.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-extra4-e8-cover.json"

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--extra-h-power", type=int, default=4)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
if args.extra_h_power < 0:
    raise ValueError("extra-h-power must be nonnegative")

endpoint = json.loads(ENDPOINT.read_text())
e8 = json.loads(E8.read_text())
assert endpoint["status"] == "PASS_EXACT_Q8_ENDPOINT_RR_AMBIENT"
assert e8["status"] == "PASS_EXACT_Q8_E8_AMBIENT_WEIGHTS"
assert endpoint["collision_polynomial"]["h_at_E8"] == "1"
floors = {(row["basis"]["kind"], int(row["basis"]["x_power"]), int(row["basis"]["m_power"])): int(row["minimal_u_power"]) for row in e8["basis_weight_floors"]}
ambient = []
for family in endpoint["families"]:
    key = (family["generic_basis"]["kind"], int(family["generic_basis"]["x_power"]), int(family["generic_basis"]["m_power"]))
    k = int(family["h_power"]) + args.extra_h_power
    e = int(family["e8_minimal_u_power"])
    assert e == floors[key]
    for i in range(e, 4*k + int(family["e7_allowed_t_denominator_power"]) + 1):
        ambient.append({"basis": key, "u_power": i, "h_power": k})
assert all(row["u_power"] >= floors[row["basis"]] for row in ambient)
payload = {
    "schema": "elkies-k3.h92-q8-enlarged-endpoint-e8-cover.v1",
    "status": "PASS_EXACT_Q8_ENLARGED_E8_COVER",
    "inputs": {"endpoint": {"path": str(ENDPOINT.relative_to(ROOT)), "sha256": digest(ENDPOINT)}, "e8_weights": {"path": str(E8.relative_to(ROOT)), "sha256": digest(E8)}},
    "extra_h_power": args.extra_h_power,
    "ambient_dimension": len(ambient),
    "reason": "h(0)=1 and every coefficient retains u_power>=its certified E8 floor",
    "boundary": "Only E8 membership is covered. Five non-marked E7 edges, characteristic-zero smooth kernel reconstruction, and all equation/bisection work remain open.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q8EXTRAENDPOINTE8|extra_h={}|ambient={}|status=PASS_EXACT_Q8_ENLARGED_E8_COVER".format(args.extra_h_power, len(ambient)), flush=True)
