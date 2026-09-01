#!/usr/bin/env python3
"""Convert a direct finite-field NS0024 MW4 seed to compact point format.

The converter makes no geometric claim of its own.  It binds the declared
abstract basis marking and copies the direct-search P4 coordinates losslessly;
the marked-point adapter independently replays the surface, four sections,
absolute component profiles, and full intersection Gram before compilation.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASIS = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-mw4-minimum-basis.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def coefficients(value):
    return [int(item) for item in value.split(",")]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--seed", type=Path, required=True)
parser.add_argument("--seed-index", type=int, default=0)
parser.add_argument(
    "--basis-marking",
    choices=("original", "resolved_component_depth_recommendation"),
    default="original",
)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

seed_path = args.seed.resolve()
output_path = args.output.resolve()
records = [line for line in seed_path.read_text().splitlines() if line.strip()]
if args.seed_index < 0 or args.seed_index >= len(records):
    raise SystemExit("--seed-index is outside the seed file")
fields = {}
for item in records[args.seed_index].strip().split("|")[1:]:
    key, value = item.split("=", 1)
    fields[key] = value
required = {"p", "P1X", "P1Y", "P2X", "P2Y", "P3X", "P3Y", "P4X", "P4Y", "P4H"}
if not required.issubset(fields):
    raise ValueError("selected record is not a complete direct MW4 seed")

basis = json.loads(BASIS.read_text())
if basis.get("status") != "PASS_EXACT_MINIMUM_POLE_FOUR_SECTION_BASIS":
    raise ValueError("pinned NS0024 MW4 basis certificate is unavailable")
if args.basis_marking == "original":
    profiles = [item["components_I7_I5_I4"] for item in basis["basis"]]
    gram = basis["section_intersection_gram"]
else:
    recommendation = basis["enumeration"]["resolved_component_depth_recommendation"]
    profiles = recommendation["profiles_I7_I5_I4"]
    gram = recommendation["section_intersection_gram"]

payload = {
    "schema": "elkies-k3.lattice-foundry-ns0024-mw4-point-modp.v1",
    "status": "PASS_EXACT_MW4_MARKED_POINT_OVER_FINITE_FIELD",
    "prime": int(fields["p"]),
    "basis_marking": args.basis_marking,
    "mw3_seed": display_path(seed_path),
    "mw3_seed_index": args.seed_index,
    "P4": {
        "X_coefficients_low_to_high": coefficients(fields["P4X"]),
        "Y_coefficients_low_to_high": coefficients(fields["P4Y"]),
        "H_coefficients_low_to_high": coefficients(fields["P4H"]),
    },
    "component_profiles_I7_I5_I4": profiles,
    "section_intersection_gram": gram,
    "inputs": {
        "seed_sha256": digest(seed_path),
        "minimum_basis_sha256": digest(BASIS),
    },
    "proof_boundary": (
        "This lossless conversion delegates every geometric source-marking "
        "claim to adapt_lattice_foundry_ns0024_mw4_point_for_edge1.sage."
    ),
    "reproduce": (
        "python3 elkies-k3/scripts/convert_lattice_foundry_ns0024_mw4_seed_to_point.py "
        "--seed {} --seed-index {} --basis-marking {} --output {}"
    ).format(
        display_path(seed_path),
        args.seed_index,
        args.basis_marking,
        display_path(output_path),
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.check:
    if output_path.read_text() != serialized:
        raise SystemExit("converted NS0024 MW4 point artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "NS0024MW4SEEDPOINT|p={}|basis={}|seed_index={}|status=PASS".format(
        fields["p"], args.basis_marking, args.seed_index
    )
)
