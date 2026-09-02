#!/usr/bin/env sage-python
"""Classify the determinant-384 A2+A5+A8/MW2 pole-[0,1] source cut."""

from __future__ import annotations

import argparse
import hashlib
import json
import runpy
from pathlib import Path

from sage.all import ZZ, identity_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
SOURCES = ROOT / "artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-semistable-mw0-2-sources-large-a-partner1-v1.json"
HELPERS_PATH = ROOT / "elkies-k3/scripts/classify_k3_6ce_a5_a10_mw2_sources.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-a2-a5-a8-mw2-source-isometries-v1.json"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    helpers = runpy.run_path(str(HELPERS_PATH), run_name="k36ce_a5_a10_helpers")
    integral_isometry = helpers["integral_isometry"]
    marking_profile = helpers["marking_profile"]
    matrix_rows = helpers["rows"]

    source_payload = json.loads(SOURCES.read_text())
    candidates = [
        row
        for row in source_payload["sources"]
        if row["source"]["root_type"] == "A2+A5+A8"
        and row["source"]["mw_rank_for_rho_19"] == 2
        and row["source"]["support_count"] == 3
        and row["source"]["torsion"] == 1
        and row["source"]["pole_audit"]["basis"] is not None
        and [
            section["pole_order"]
            for section in row["source"]["pole_audit"]["basis"]
        ]
        == [0, 1]
    ]
    candidates.sort(key=lambda row: row["source_id"])
    if len(candidates) != 17:
        raise ArithmeticError("A2+A5+A8 pole-[0,1] MW2 cut changed")

    classes = []
    for row in candidates:
        gram = matrix(ZZ, row["source"]["gram"])
        selected_class = None
        selected_isometry = None
        for class_row in classes:
            isometry = integral_isometry(gram, class_row["representative_gram"])
            if isometry is not None:
                selected_class = class_row
                selected_isometry = isometry
                break
        if selected_class is None:
            selected_class = {
                "representative_source_id": row["source_id"],
                "representative_gram": gram,
                "members": [],
            }
            classes.append(selected_class)
            selected_isometry = identity_matrix(ZZ, 17)
        selected_class["members"].append(
            {
                "source_id": row["source_id"],
                "source_gram_sha256": row["source"]["gram_sha256"],
                "source_to_representative_isometry": matrix_rows(selected_isometry),
                "isometry_determinant": int(selected_isometry.det()),
                "marking_profile": marking_profile(row["source"], (2, 5, 8)),
            }
        )

    output_classes = []
    profile_count = 0
    for class_index, class_row in enumerate(classes, start=1):
        groups = {}
        for member in class_row["members"]:
            key = json.dumps(member["marking_profile"], sort_keys=True)
            groups.setdefault(key, []).append(member["source_id"])
        profiles = [
            {"marking_profile": json.loads(key), "source_ids": source_ids}
            for key, source_ids in sorted(groups.items())
        ]
        profile_count += len(profiles)
        output_classes.append(
            {
                "class_index": class_index,
                "representative_source_id": class_row["representative_source_id"],
                "member_count": len(class_row["members"]),
                "marked_basis_profiles": profiles,
                "members": class_row["members"],
            }
        )

    payload = {
        "schema": "elkies-k3.k3-6ce-a2-a5-a8-mw2-source-isometries.v1",
        "status": "PASS_EXACT_INTEGRAL_ISOMETRY_AND_MARKING_PROFILE_CLASSIFICATION",
        "ideal_cut": {
            "surface_id": "K3-6ce16abb9de3c7c5",
            "target_mw_rank": 15,
            "root_type": "A2+A5+A8",
            "mw_rank": 2,
            "semistable": True,
            "support_count": 3,
            "torsion_order": 1,
            "complete_basis_pole_profile": [0, 1],
        },
        "accounting": {
            "reduced_gram_rows": len(candidates),
            "integral_isometry_classes": len(classes),
            "class_sizes": sorted(len(row["members"]) for row in classes),
            "selected_physical_basis_profiles": profile_count,
        },
        "classes": output_classes,
        "proof_boundary": {
            "proved": (
                "Every selected reduced Gram is placed in an exact integral-isometry "
                "class, and every displayed basis profile has exact component depths, "
                "height Gram, pole orders, and section intersection."
            ),
            "not_proved": (
                "No rational marking, equation, characteristic-zero lift, arithmetic "
                "moduli curve, or neighbour route is asserted."
            ),
        },
        "inputs": {
            relative(SOURCES): digest(SOURCES),
            relative(HELPERS_PATH): digest(HELPERS_PATH),
        },
        "reproduce": "sage -python elkies-k3/scripts/classify_k3_6ce_a2_a5_a8_mw2_sources.sage",
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "K36CEA2A5A8|"
        f"rows={len(candidates)}|classes={len(classes)}|profiles={profile_count}|status=PASS"
    )


if __name__ == "__main__":
    main()
