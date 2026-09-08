#!/usr/bin/env python3
"""Extract anonymous target-free A1/MW16 family presentations."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
AUDIT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_presentation_audit_v1.json"
OUTPUT = ROOT / "elliptic-curves/data/a1_mw16_family_template_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def build(source: Path, audit_path: Path) -> dict:
    payload = json.loads(source.read_text())
    audit = json.loads(audit_path.read_text())
    if len(payload.get("parents", [])) != 9 or audit.get("exact_fibration_class_count") != 5:
        raise ArithmeticError("expected nine presentations of five fibrations")
    presentation_ids = {
        row["parent_id"]: f"a1-presentation-{index:02d}"
        for index, row in enumerate(payload["parents"], 1)
    }
    fibration_ids = {}
    for index, cluster in enumerate(audit["clusters"], 1):
        for source_id in cluster["presentation_ids"]:
            fibration_ids[source_id] = f"a1-fibration-{index:02d}"
    if set(presentation_ids) != set(fibration_ids):
        raise ArithmeticError("presentation clustering is incomplete")
    presentations = []
    for row in payload["parents"]:
        source_id = row["parent_id"]
        presentations.append({
            "presentation_id": presentation_ids[source_id],
            "fibration_id": fibration_ids[source_id],
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
                    "base_maps_lambda_of_old_t",
                )
            },
        })
    template = {
        "schema": "elliptic-curves.a1-mw16-family-template.v1",
        "status": "PASS_TARGET_FREE_A1_MW16_FAMILY_PRESENTATIONS",
        "fibration_count": 5,
        "presentation_count": 9,
        "presentations": presentations,
        "extraction_provenance_sha256": {
            "equation_and_marking_source": digest(source),
            "presentation_equivalence_source": digest(audit_path),
        },
        "claim_boundary": (
            "This sanitized equation-and-section template contains five exact "
            "A1/MW16 fibrations in nine anonymous coordinate presentations. It "
            "contains no target curve, target parameter, public point, public "
            "rank, target j-invariant, displayed jump, or target isomorphism."
        ),
    }
    forbidden = (
        '"parent_id"',
        '"curve_id"',
        '"priority_rank"',
        '"target_parameter"',
        '"target_short_model"',
        '"specialized_generic_points"',
        '"target_isomorphism"',
        '"known_target_fibre_local_control"',
        "curve398",
        "curve400",
        "curve401",
        "curve542",
        "curve548",
    )
    serialized = json.dumps(template, sort_keys=True).lower()
    if any(key in serialized for key in forbidden):
        raise ArithmeticError("target field leaked into family template")
    return template


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--audit", type=Path, default=AUDIT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build(args.source, args.audit)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file():
            raise ArithmeticError("stored A1/MW16 family template is absent")
        stored = json.loads(args.output.read_text())
        stored_provenance = stored.pop("extraction_provenance_sha256", None)
        fresh_provenance = payload.pop("extraction_provenance_sha256", None)
        if stored != payload:
            raise ArithmeticError("stored A1/MW16 family template differs")
        if (
            not isinstance(stored_provenance, dict)
            or stored_provenance.get("equation_and_marking_source")
            != fresh_provenance["equation_and_marking_source"]
        ):
            raise ArithmeticError("sealed equation-and-marking provenance changed")
        if stored_provenance != fresh_provenance:
            print(
                "A1MW16TEMPLATE|audit_bytes_changed=true|"
                "anonymous_partition_and_equations_unchanged=true"
            )
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print("A1MW16TEMPLATE|status=PASS_TARGET_FREE_A1_MW16_FAMILY_PRESENTATIONS")


if __name__ == "__main__":
    main()
