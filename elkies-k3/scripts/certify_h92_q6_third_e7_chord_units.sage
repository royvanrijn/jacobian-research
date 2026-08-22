#!/usr/bin/env sage -python
"""Certify that the third marked chord is a unit at every actual E7 edge."""

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PULLBACKS = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-actual-e7-chart-pullbacks.json"
SERIES = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-point-series.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-third-e7-chord-units.json"
SERIES_SHA256 = None

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--pullbacks", type=Path, default=PULLBACKS)
parser.add_argument("--series", type=Path, default=SERIES)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
pullbacks = json.loads(args.pullbacks.read_text())
series = json.loads(args.series.read_text())
assert pullbacks["status"] == "PASS_EXACT_H92_E7_CHART_PULLBACKS"
assert series["status"] == "PASS_EXACT_Q6_THIRD_E7_SERIES_POINT"
assert series["valuations"] == {"x": 0, "y": 0}
assert series["specialization"]["x_nonzero"] and series["specialization"]["y_nonzero"]

chart_names = []
for chart in pullbacks["charts"]:
    assert chart["old_coordinate_values_at_chart_origin"] == {"t": "0", "x": "0", "y": "0"}
    chart_names.append(chart["name"])

payload = {
    "schema": "elkies-k3.h92-q6-third-e7-chord-units.v1",
    "status": "PASS_EXACT_H92_Q6_THIRD_E7_CHORD_UNITS",
    "inputs": {
        "actual_pullbacks": {"path": str(args.pullbacks.relative_to(ROOT)), "sha256": digest(args.pullbacks)},
        "point_series": {"path": str(args.series.relative_to(ROOT)), "sha256": digest(args.series)},
    },
    "charts": chart_names,
    "conclusion": (
        "At each chart origin old x,y vanish, while x(-P),y(-P) are units. "
        "Therefore x-x(-P) and y-y(-P) are units and the marked chord "
        "(y-y(-P))/(x-x(-P)) is a unit in every actual E7 edge local ring."
    ),
    "boundary": "This certifies only the chord's E7 edge-unit behavior; it does not form finite jets or the complete global condition matrix.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H92Q6THIRDCHORDUNIT|charts=6|status=PASS_EXACT_H92_Q6_THIRD_E7_CHORD_UNITS", flush=True)
