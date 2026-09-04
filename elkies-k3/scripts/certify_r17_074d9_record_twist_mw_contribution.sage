#!/usr/bin/env sage-python
"""Combine four twist bounds and exact record-fibre quotient images."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from sage.all import matrix, QQ


ROOT = Path(__file__).resolve().parents[2]
BOUNDS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-twist-good-reduction-bounds-v1.json"
)
SECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-record-twist-sections-v1.json"
)
DESCENT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-twist-2descent-audit-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-record-twist-mw-contribution-v1.json"
)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def build_payload():
    bounds = json.loads(BOUNDS.read_text())
    sections = json.loads(SECTIONS.read_text())
    descent = json.loads(DESCENT.read_text())
    if bounds.get("status") != "PASS_EXACT_FOUR_TWIST_GOOD_REDUCTION_RANK_BOUNDS":
        raise ValueError("unexpected twist-bound status")
    if sections.get("status") != "PASS_EXACT_FOUR_RECORD_TWIST_SECTIONS_AND_SPECIALIZATIONS":
        raise ValueError("unexpected twist-section status")
    if descent.get("status") != "INCOMPLETE_GOOD_REDUCTION_2DESCENTS_TIMED_OUT":
        raise ValueError("unexpected twist-descent status")

    bound_by_label = {row["label"]: row for row in bounds["twists"]}
    descent_by_label = {row["label"]: row for row in descent["twists"]}
    twist_rows = []
    for section in sections["records"]:
        label = section["label"]
        bound = int(bound_by_label[label]["best_geometric_MW_rank_upper_bound"])
        if (
            not section["exact_identities_verified"]
            or section["height_on_twist_surface"] != 6
            or section["P_dot_O"] != 0
        ):
            raise ArithmeticError("known twist section lost its exact lower-bound gates")
        if descent_by_label[label]["characteristic_zero_function_field_rank_status"] != "UNKNOWN":
            raise ArithmeticError("descent status changed without updating this certificate")
        image = list(map(int, section["twist_section_image_quotient_class"]))
        branch = list(map(int, section["bisection_branch_quotient_class"]))
        if image != [2 * value for value in branch]:
            raise ArithmeticError("twist-section quotient image is not twice the branch image")
        twist_rows.append(
            {
                "label": label,
                "curve_id": int(section["curve_id"]),
                "QQ_u_MW_rank_lower_bound": 1,
                "geometric_MW_rank_upper_bound": bound,
                "QQ_u_MW_rank_status": "UNKNOWN",
                "known_section": {
                    "P_dot_O": 0,
                    "height": 6,
                    "record_fibre_quotient_image": image,
                    "specialization_identity": section["specialization_identity"],
                },
            }
        )

    fibre_rows = []
    for curve_id in (356, 385):
        rows = [row for row in twist_rows if row["curve_id"] == curve_id]
        if len(rows) != 2:
            raise ArithmeticError("record fibre did not retain two twists")
        images = [row["known_section"]["record_fibre_quotient_image"] for row in rows]
        known_image_rank = int(matrix(QQ, images).rank())
        if known_image_rank != 2:
            raise ArithmeticError("known record-twist images lost independence")
        combined_upper = sum(row["geometric_MW_rank_upper_bound"] for row in rows)
        quotient_rank = len(images[0])
        fibre_rows.append(
            {
                "curve_id": curve_id,
                "exceptional_quotient_rank": quotient_rank,
                "twist_labels": [row["label"] for row in rows],
                "certified_known_specialization_image_rank": known_image_rank,
                "combined_specialization_image_rank_upper_bound": combined_upper,
                "quotient_equality": False,
                "rank_obstruction": f"{combined_upper} < {quotient_rank}",
                "minimum_fibre_specific_quotient_directions": quotient_rank - combined_upper,
                "all_four_twists_rank_one_scenario_status": "UNRESOLVED",
            }
        )
    return {
        "schema": "elkies-k3.r17-074d9-record-twist-mw-contribution.v1",
        "status": "PROVED_NO_FULL_RECORD_QUOTIENT_LIFT_RANKS_REMAIN_UNKNOWN",
        "claim": (
            "Good-reduction geometric bounds disprove both proposed full exceptional-"
            "quotient specialization equalities. One exact height-six section per "
            "twist supplies the current rank-one lower bound and exact quotient image; "
            "the four full QQ(u)-ranks remain unknown."
        ),
        "twists": twist_rows,
        "fibres": fibre_rows,
        "proof_boundary": (
            "The geometric upper bounds are sufficient for the two negative equality "
            "tests. They do not determine the individual QQ(u)-ranks. The scenario in "
            "which all four ranks equal one is neither proved nor disproved."
        ),
        "inputs": {
            relative(path): digest(path) for path in (BOUNDS, SECTIONS, DESCENT)
        },
        "generation": {
            "command": (
                "sage -python elkies-k3/scripts/"
                "certify_r17_074d9_record_twist_mw_contribution.sage"
            )
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("stale record-twist MW contribution certificate")
        terminal = "PASS"
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        terminal = "WROTE"
    print(
        "R17074D9MWCONTRIBUTION|"
        + ",".join(
            f"{row['curve_id']}:{row['combined_specialization_image_rank_upper_bound']}"
            f"<{row['exceptional_quotient_rank']}"
            for row in payload["fibres"]
        )
        + f"|status={terminal}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
