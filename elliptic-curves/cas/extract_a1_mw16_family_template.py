#!/usr/bin/env python3
"""Extract the target-free priority-16875 A1/MW16 family template."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
OUTPUT = ROOT / "elliptic-curves/data/a1_mw16_family_template_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def build(source: Path) -> dict:
    payload = json.loads(source.read_text())
    matches = [
        row for row in payload["parents"]
        if row["parent_id"] == "curve398-p16875"
    ]
    if len(matches) != 1:
        raise ArithmeticError("canonical A1/MW16 source presentation changed")
    row = matches[0]
    template = {
        "schema": "elliptic-curves.a1-mw16-family-template.v1",
        "status": "PASS_TARGET_FREE_A1_MW16_FAMILY_TEMPLATE",
        "family_id": "a1-mw16-11952-norm8-canonical",
        "generic_rank": 16,
        "generic_height_gram": row["generic_height_gram"],
        "pencil": row["pencil"],
        "source_marking": {
            key: row["source_marking"][key]
            for key in (
                "trace_section_basis_w",
                "new_zero_source_section_basis_coordinates",
                "generic_source_section_basis_coordinates",
                "generic_coordinates_in_compiled_mw_basis",
            )
        },
        "source_hash": digest(source),
        "claim_boundary": (
            "This is a sanitized equation-and-section template. It contains no "
            "target curve, target parameter, public point, public rank, target "
            "j-invariant, displayed jump, or target isomorphism."
        ),
    }
    forbidden = (
        "curve_id", "target_parameter", "target_short_model",
        "specialized_generic_points", "target_isomorphism",
    )
    serialized = json.dumps(template, sort_keys=True)
    if any(f'"{key}"' in serialized for key in forbidden):
        raise ArithmeticError("target field leaked into family template")
    return template


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build(args.source)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise ArithmeticError("stored A1/MW16 family template differs")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print("A1MW16TEMPLATE|status=PASS_TARGET_FREE_A1_MW16_FAMILY_TEMPLATE")


if __name__ == "__main__":
    main()
