#!/usr/bin/env sage -python
"""Aggregate actual H92 q8 E7-node modular resolved-chart obstructions.

This is deliberately a regression aggregator rather than a proof upgrade.  It
checks that each supplied local-normal-form image or one-way Artinian-corner
obstruction was computed from the pinned six-node clearing atlas, over one
prime and one ambient, and records their finite-image ranks.  It makes no
characteristic-zero, overlap, pencil, or child-model assertion.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CLEARINGS = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-node-principal-clearings.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-e7-node-local-normal-form-mod-43.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--clearings", type=Path, default=CLEARINGS)
parser.add_argument("--node-image", type=Path, action="append", required=True)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
args.clearings = args.clearings.resolve()
args.node_image = [path.resolve() for path in args.node_image]
args.output = args.output.resolve()

clearings = json.loads(args.clearings.read_text())
assert clearings["status"] == "PASS_EXACT_Q8_E7_NODE_PRINCIPAL_CLEARINGS"
expected_charts = {entry["chart"] for entry in clearings["nodes"]}
assert len(expected_charts) == 6
assert len(args.node_image) == len(expected_charts)

records = []
prime = ambient_dimension = None
for path in args.node_image:
    payload = json.loads(path.read_text())
    assert payload["status"] in {
        "EXPERIMENTAL_MODULAR_Q8_E7_NODE_LOCAL_NORMAL_FORM_BLOCK",
        "EXPERIMENTAL_MODULAR_Q8_E7_NODE_FINITE_CORNER_OBSTRUCTION",
    }
    chart = payload["local_ring"]["chart"]
    image = payload["finite_ambient_image"]
    if payload["status"] == "EXPERIMENTAL_MODULAR_Q8_E7_NODE_LOCAL_NORMAL_FORM_BLOCK":
        assert payload["local_ring"]["order"] == "Singular ds local degree order at (Z,U,Y)"
        condition_mode = "exact_local_normal_form"
    else:
        assert payload["local_ring"]["order"] == "Singular dp degree order for the Artinian corner quotient"
        assert payload["local_ring"]["principal_ideal"].startswith("(surface,Z^")
        condition_mode = "one_way_artinian_corner_obstruction"
    assert payload["inputs"]["node_clearings"]["sha256"] == digest(args.clearings)
    assert payload["inputs"]["node_clearings"]["path"] == str(args.clearings.relative_to(ROOT))
    if prime is None:
        prime = int(payload["prime"])
        ambient_dimension = int(image["ambient_dimension"])
    assert int(payload["prime"]) == prime
    assert int(image["ambient_dimension"]) == ambient_dimension
    assert int(image["rank"]) == ambient_dimension
    assert int(image["kernel_dimension"]) == 0
    records.append({
        "chart": chart,
        "path": str(path.relative_to(ROOT)),
        "sha256": digest(path),
        "normal_form_coordinates": int(image["rows"]),
        "rank": int(image["rank"]),
        "kernel_dimension": int(image["kernel_dimension"]),
        "condition_mode": condition_mode,
    })
assert {entry["chart"] for entry in records} == expected_charts

payload = {
    "schema": "elkies-k3.h92-q8-all-e7-node-resolved-modular-obstruction.v1",
    "status": "EXPERIMENTAL_MODULAR_Q8_ALL_E7_NODE_RESOLVED_OBSTRUCTION",
    "prime": prime,
    "inputs": {
        "node_clearings": {"path": str(args.clearings.relative_to(ROOT)), "sha256": digest(args.clearings)},
        "node_images": records,
    },
    "summary": {
        "nodes": len(records),
        "ambient_dimension": ambient_dimension,
        "each_resolved_chart_obstruction_has_full_column_rank": True,
        "condition_modes": {record["chart"]: record["condition_mode"] for record in records},
    },
    "boundary": (
        "This aggregates six good-prime finite-ambient resolved-chart "
        "obstructions. A chart recorded as an Artinian corner is one-way: it "
        "rules out true local solutions but is not the local quotient. This "
        "does not give characteristic-zero local matrices, chart-overlap "
        "compatibility, a common q8 kernel, h0(D), a pencil, or a child model."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "H92Q8ALLE7NODERESOLVED|prime={}|nodes={}|ambient={}|"
    "status=EXPERIMENTAL_MODULAR_Q8_ALL_E7_NODE_RESOLVED_OBSTRUCTION".format(
        prime, len(records), ambient_dimension,
    ),
    flush=True,
)
