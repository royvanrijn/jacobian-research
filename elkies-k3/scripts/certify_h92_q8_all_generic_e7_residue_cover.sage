#!/usr/bin/env sage -python
"""Certify the complete all-component *generic* E7 q=8 residue cover."""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
Y_BRANCH = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-1-3-generic-residue-rows.json"
CONICS = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-5-6-generic-residue-rows.json"
SIMPLE = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residue-rows.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-all-generic-e7-residue-cover.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--y-branch", type=Path, default=Y_BRANCH)
parser.add_argument("--conics", type=Path, default=CONICS)
parser.add_argument("--simple", type=Path, default=SIMPLE)
parser.add_argument(
    "--allow-enlarged", action="store_true",
    help="accept a common enlarged endpoint ambient instead of the pinned 54-column seed",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

payloads = [json.loads(path.read_text()) for path in (args.y_branch, args.conics, args.simple)]
assert [payload["status"] for payload in payloads] == [
    "PASS_EXACT_Q8_E7_1_3_GENERIC_RESIDUE_ROWS",
    "PASS_EXACT_Q8_E7_5_6_GENERIC_RESIDUE_ROWS",
    "PASS_EXACT_Q8_E7_4_7_GENERIC_RESIDUE_ROWS",
]
assert len({payload["ambient_basis_sha256"] for payload in payloads}) == 1
components = [component for payload in payloads for component in payload["components"]]
assert [component["component"] for component in sorted(components, key=lambda value: value["component"])] == [
    "E7_1", "E7_2", "E7_3", "E7_4", "E7_5", "E7_6", "E7_7",
]
counts = {
    component["component"]: len(component["non_singleton_residue_rows"])
    for component in components
}
if not args.allow_enlarged:
    assert counts == {"E7_1": 42, "E7_2": 189, "E7_3": 391, "E7_4": 42,
                      "E7_5": 228, "E7_6": 57, "E7_7": 34}

result = {
    "schema": "elkies-k3.h92-q8-all-generic-e7-residue-cover.v1",
    "status": "PASS_EXACT_Q8_ALL_GENERIC_E7_RESIDUE_COVER",
    "inputs": {
        "y_branch_rows": {"path": path_label(args.y_branch), "sha256": digest(args.y_branch)},
        "conic_rows": {"path": path_label(args.conics), "sha256": digest(args.conics)},
        "simple_rows": {"path": path_label(args.simple), "sha256": digest(args.simple)},
    },
    "ambient_basis_sha256": payloads[0]["ambient_basis_sha256"],
    "enlarged_ambient": bool(args.allow_enlarged),
    "component_row_counts": counts,
    "total_non_singleton_residue_rows": sum(counts.values()),
    "conclusion": (
        "Every actual resolved E7 component now has an exact generic-point "
        "leading-residue evaluator for every non-singleton negative-order "
        "source-q8 seed group."
    ),
    "boundary": (
        "This is the complete generic-component residue cover only. It does "
        "not include edge nodes, marked branch jets, chart overlaps, E8/smooth "
        "assembly, a characteristic-zero global kernel, or a q8 pencil."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
if args.allow_enlarged:
    print(
        "H92Q8ALLGENERICRESIDUES|components=7|rows={}|enlarged=1|status="
        "PASS_EXACT_Q8_ALL_GENERIC_E7_RESIDUE_COVER".format(
            result["total_non_singleton_residue_rows"]
        ),
        flush=True,
    )
else:
    print(
        "H92Q8ALLGENERICRESIDUES|components=7|rows={}|status="
        "PASS_EXACT_Q8_ALL_GENERIC_E7_RESIDUE_COVER".format(
            result["total_non_singleton_residue_rows"]
        ),
        flush=True,
    )
