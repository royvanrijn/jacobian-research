#!/usr/bin/env sage-python
"""Freeze generic-MW17-only inputs for the refreshed R17 jump ladder.

This boundary builder may read the complete refresh quotient certificate and
public point projection.  Its output contains only each curve, the seventeen
exactly transported generic sections, and the exact generic height form.  The
displayed complement coordinates, displayed rank, and jump label are omitted.

Curve 499 is excluded before any blind search because the exact refresh audit
proves that its displayed subgroup does not contain the specialized generic
subgroup.  Hence its displayed-subgroup quotient by MW17 is not defined.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from sage.all import EllipticCurve, Matrix, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
OVERVIEW = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "icarm_curve_refresh_475_573_overview_v1.json"
)
QUOTIENTS = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-r17-norm12-refresh-priority-quotients-v1.json"
)
OUTPUT = ROOT / "elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json"

DIRECT_BY_CHART = {
    "norm12-orbit-07ca9": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit07ca9-direct-fibration-v1.json",
    "norm12-orbit-08234": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08234-direct-fibration-v1.json",
    "norm12-orbit-11952": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json",
    "norm12-orbit-103b2": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit103b2-direct-fibration-v1.json",
    "norm12-orbit-08f72": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08f72-direct-fibration-v1.json",
}
WGXLI = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"

# The 074d9 direct records are the first seventeen marked sections, whereas
# the stored saturated generic Gram uses these unimodular words in them.
GENERIC_WORDS = (
    ((2, 1),), ((3, 1),), ((4, 1),), ((5, 1),), ((8, 1),),
    ((11, 1),), ((13, 1),), ((15, 1),), ((16, 1),), ((17, 1),),
    ((1, 1), (2, -1)), ((1, 1), (6, -1)), ((1, 1), (7, -1)),
    ((1, 1), (9, -1)), ((1, 1), (10, -1)), ((1, 1), (12, -1)),
    ((1, 1), (14, -1)),
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    return str(value)


def short_curve_and_points(public_record):
    a1, a2, a3, a4, a6 = map(QQ, public_record["ainvs"])
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    curve = EllipticCurve(QQ, [-c4 / 48, -c6 / 864])
    points = tuple(
        curve(
            QQ(x_value) + b2 / 12,
            QQ(y_value) + (a1 * QQ(x_value) + a3) / 2,
        )
        for x_value, y_value in public_record["points"]
    )
    return curve, points


def generic_gram(native_chart: str):
    if native_chart != "norm12-orbit-074d9":
        direct = json.loads(DIRECT_BY_CHART[native_chart].read_text())
        gram = Matrix(ZZ, direct["sections"]["height_gram"])
    else:
        lineage = json.loads(WGXLI.read_text())
        word_gram = Matrix(ZZ, lineage["generic_basis"]["height_gram"])
        word_matrix = Matrix(
            ZZ,
            17,
            17,
            lambda row, column: next(
                (
                    coefficient
                    for index, coefficient in GENERIC_WORDS[row]
                    if index - 1 == column
                ),
                0,
            ),
        )
        if abs(word_matrix.det()) != 1:
            raise ArithmeticError("the 074d9 generic word change is not unimodular")
        inverse = word_matrix.inverse()
        gram = inverse * word_gram * inverse.transpose()
    if gram.dimensions() != (17, 17) or gram.det() != 948:
        raise ArithmeticError(f"{native_chart}: generic height form changed")
    if any(value not in ZZ for value in gram.list()):
        raise ArithmeticError(f"{native_chart}: generic height form is not integral")
    return Matrix(ZZ, gram)


def point_record(point):
    return {"x": rational_text(point[0]), "y": rational_text(point[1])}


def build():
    overview = json.loads(OVERVIEW.read_text())
    quotient = json.loads(QUOTIENTS.read_text())
    if overview.get("status") != "PASS_EXACT_OVERVIEW_OF_ICARM_CURVES_475_THROUGH_573":
        raise ArithmeticError("the refresh overview is not exact")
    if quotient.get("status") != "PASS_EXACT_REFRESH_ATLAS_HIT_SPECIALIZATION_AUDIT":
        raise ArithmeticError("the refreshed quotient audit is not exact")
    if quotient["summary"]["new_atlas_hit_count"] != 17:
        raise ArithmeticError("the refresh no longer contains exactly seventeen atlas hits")
    if quotient["summary"]["noninclusive_displayed_subgroup_curve_ids"] != [499]:
        raise ArithmeticError("the predeclared quotient-ineligible case changed")

    public = {int(row["id"]): row for row in overview["snapshot"]["records"]}
    cases = []
    for fibre in sorted(quotient["fibres"], key=lambda row: int(row["curve_id"])):
        curve_id = int(fibre["curve_id"])
        curve, public_points = short_curve_and_points(public[curve_id])
        coordinates = Matrix(
            ZZ,
            fibre["specialized_generic_subgroup"][
                "coordinate_matrix_rows_in_ordered_public_points"
            ],
        )
        if coordinates.dimensions() != (len(public_points), 17):
            raise ArithmeticError(f"curve {curve_id}: generic coordinate shape changed")
        generic_points = []
        for column in range(17):
            point = curve(0)
            for coefficient, public_point in zip(coordinates.column(column), public_points):
                point += int(coefficient) * public_point
            if point == curve(0):
                raise ArithmeticError(f"curve {curve_id}: a generic point became torsion")
            generic_points.append(point)
        gram = generic_gram(fibre["native_chart"])
        cases.append(
            {
                "curve_id": curve_id,
                "representative_class": fibre["representative_class"],
                "native_chart": fibre["native_chart"],
                "short_model": [
                    rational_text(0),
                    rational_text(0),
                    rational_text(0),
                    rational_text(curve.a4()),
                    rational_text(curve.a6()),
                ],
                "generic_points": [point_record(point) for point in generic_points],
                "generic_height_gram": [
                    [int(value) for value in row] for row in gram.rows()
                ],
                "generic_height_gram_determinant": int(gram.det()),
            }
        )

    expected_ids = [478, 498, 531, 532, 534, 535, 536, 537, 538, 539, 540, 541, 543, 544, 545, 546]
    if [row["curve_id"] for row in cases] != expected_ids:
        raise ArithmeticError("the blinded quotient-eligible case inventory changed")
    source_paths = [OVERVIEW, QUOTIENTS, WGXLI, *DIRECT_BY_CHART.values()]
    return {
        "schema": "elliptic-curves.r17-refresh-jump-ladder-blind-inputs.v1",
        "status": "FROZEN_MW17_ONLY_NO_PUBLIC_COMPLEMENT",
        "case_count": len(cases),
        "cases": cases,
        "redaction": {
            "builder_read_public_complement": True,
            "blind_runner_may_read_public_complement": False,
            "contains_displayed_complement_coordinates": False,
            "contains_displayed_rank_or_jump": False,
            "contains_public_points_other_than_reconstructed_generic_MW17": False,
            "case_order": "increasing curve id, independent of displayed rank or jump",
        },
        "pre_search_exclusion": {
            "curve_id": 499,
            "reason": "displayed quotient by specialized MW17 is not defined because MW17 is not contained in the displayed subgroup",
            "decision_used_no_blind_search_outcome": True,
        },
        "source_hashes_available_only_to_boundary_builder": {
            relative(path): digest(path) for path in source_paths
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
        },
        "claim_boundary": (
            "This file certifies an information boundary, not a search outcome. "
            "The blind runner receives only the curve, specialized generic MW17, "
            "exact generic height form, curve id, and atlas class."
        ),
    }


def main():
    payload = build()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"R17JUMPLADDERINPUT|status=FROZEN|cases={payload['case_count']}|"
        f"sha256={digest(OUTPUT)}|output={relative(OUTPUT)}"
    )


if __name__ == "__main__":
    main()
