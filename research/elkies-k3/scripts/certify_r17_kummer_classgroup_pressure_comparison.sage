#!/usr/bin/env sage-python
"""Certify the six-fibre Kummer/class-group pressure comparison.

Status: ACTIVE_PROOF
Claim: exact finite comparison of point-forced cubic 2-class directions.
Inputs: the five-fibre pressure certificate and the native alternate-Q80
curve-12 quotient certificate.
Output: artifacts/generated-results/
elkies-k3-r17-kummer-classgroup-pressure-comparison-v1.json.

For the five published-R17 controls this replays the existing half-ideal and
bad-valuation calculation byte-for-byte.  For alternate-Q80 curve 12 it first
changes the 29 public points by the certified unimodular matrix to the ordered
basis

    G1,...,G17,Q1,...,Q12,

where the G columns are the specialized generic MW17 basis and the Q columns
are the certified quotient basis.  It then invokes exactly the same cubic,
valuation, unit-ambiguity, and point-half-ideal audit.

The reported residual pressure is a lower bound for the quotient of the
known point-forced 2-class image by its generic-MW17 contribution.  It is not
an exact class-group rank, S-class rank, Selmer upper bound, or prospective
rank predictor: the exceptional point directions are inputs to the invariant.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
from typing import Any

from sage.all import EllipticCurve, GF, Matrix, QQ, ZZ, pari
from sage.version import version as sage_version


ROOT = Path(__file__).resolve().parents[2]
BASE_SCRIPT = ROOT / "elkies-k3/scripts/certify_r17_kummer_classgroup_pressure.sage"
BASE_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-kummer-classgroup-pressure-v1.json"
)
PUBLIC = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
CURVE12 = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json"
)
RANK29 = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_elkies_klagsbrun_rank29_certificate.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-kummer-classgroup-pressure-comparison-v1.json"
)
SOURCE = Path(__file__).resolve()

SCHEMA = "elkies-k3.r17-kummer-classgroup-pressure-comparison.v1"
STATUS = "PROVED_SIX_FIBRE_KUMMER_CLASSGROUP_PRESSURE_COMPARISON"
PROTOCOL = "R17KUMMERCL2COMPARE"
TARGET_IDS = (351, 356, 376, 377, 385, 12)
GENERIC_RANK = 17
CURVE12_QUOTIENT_PUBLIC_INDICES = (2, 11, 4, 3, 6, 8, 17, 10, 28, 24, 19, 15)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_base_module() -> dict[str, Any]:
    return runpy.run_path(str(BASE_SCRIPT), run_name="r17_kummer_pressure_base")


def transformed_signature_certificate(rank29: dict[str, Any], change: Matrix):
    signatures = []
    for signature in rank29["finite_reduction_certificate"]["signatures"]:
        old_rows = Matrix(GF(2), signature["rows"])
        new_rows = old_rows * change.change_ring(GF(2))
        signatures.append(
            {
                **signature,
                "rows": [
                    [int(value) for value in row] for row in new_rows.rows()
                ],
            }
        )
    return {
        "proved_displayed_subgroup_rank": 29,
        "mod2_reduction_signatures": signatures,
    }


def relabel_curve12_record(record: dict[str, Any]) -> dict[str, Any]:
    labels = [f"G{index}" for index in range(1, 18)] + [
        f"Q{index}=P{public_index}"
        for index, public_index in enumerate(CURVE12_QUOTIENT_PUBLIC_INDICES, start=1)
    ]
    label_map = {f"P{index}": label for index, label in enumerate(labels, start=1)}
    for row in record["known_point_bad_valuation_parity_rows"]:
        row["label"] = label_map[row["label"]]
    for row in record["point_half_ideals"]:
        row["label"] = label_map[row["label"]]
    for row in record["residual_adjustments_by_mw17"]:
        row["residual_label"] = label_map[row["residual_label"]]
        row["mw17_labels_to_add"] = [
            label_map[label] for label in row["mw17_labels_to_add"]
        ]
    kernel = record["everywhere_even_known_kummer_kernel"]
    kernel["ordered_basis_labels"] = labels
    kernel["basis_rows_in_ordered_basis_coordinates"] = kernel.pop(
        "basis_rows_in_P1_through_Pn_coordinates"
    )
    record["ordered_half_ideal_basis_labels"] = labels
    record["generic_mw17_basis_labels"] = labels[:GENERIC_RANK]
    record["displayed_quotient_basis_labels"] = labels[GENERIC_RANK:]
    record["native_frame"] = "alternate-Q80"
    record["native_chart"] = "norm12-orbit-11952"
    return record


def audit_curve12(base: dict[str, Any]) -> dict[str, Any]:
    public = json.loads(PUBLIC.read_text())
    curve12 = json.loads(CURVE12.read_text())
    rank29 = json.loads(RANK29.read_text())
    if public.get("status") != base["PUBLIC_STATUS"]:
        raise ArithmeticError("the public-fibre certificate is not passing")
    if curve12.get("status") != "PROVED_CURVE12_NATIVE_ALTERNATE_Q80_AND_DISPLAYED_QUOTIENT":
        raise ArithmeticError("the native alternate-Q80 quotient certificate changed")
    if rank29.get("status") != "exact_unconditional_rank_at_least_29":
        raise ArithmeticError("the curve-12 rank lower-bound certificate changed")
    if rank29["finite_reduction_certificate"]["combined_exact_rank_over_F2"] != 29:
        raise ArithmeticError("the curve-12 public points lost mod-two independence")

    public_record = next(row for row in public["records"] if int(row["id"]) == 12)
    if public_record["ainvs"] != rank29["general_weierstrass_coefficients"]:
        raise ArithmeticError("the two curve-12 model certificates disagree")
    if public_record["points"] != [
        [row["x"], row["y"]] for row in rank29["published_points"]
    ]:
        raise ArithmeticError("the two curve-12 public point inventories disagree")

    generic_coordinates = Matrix(
        ZZ,
        curve12["specialized_generic_subgroup"][
            "coordinate_matrix_rows_in_ordered_29_public_points"
        ],
    )
    quotient_coordinates = Matrix(
        ZZ,
        29,
        12,
        lambda row, column: int(
            row + 1 == CURVE12_QUOTIENT_PUBLIC_INDICES[column]
        ),
    )
    change = generic_coordinates.augment(quotient_coordinates)
    if abs(change.det()) != 1:
        raise ArithmeticError("the curve-12 generic-plus-quotient basis is not unimodular")
    expected_labels = [
        f"P{index}" for index in CURVE12_QUOTIENT_PUBLIC_INDICES
    ]
    if (
        curve12["displayed_exceptional_quotient"]
        ["free_basis_modulo_specialized_generic"]
        != expected_labels
    ):
        raise ArithmeticError("the curve-12 quotient basis changed")

    curve = EllipticCurve(QQ, [QQ(value) for value in public_record["ainvs"]])
    public_points = [
        curve(QQ(point[0]), QQ(point[1])) for point in public_record["points"]
    ]
    changed_points = []
    for column in change.columns():
        point = curve(0)
        for coefficient, public_point in zip(column, public_points):
            point += ZZ(coefficient) * public_point
        if point.is_zero():
            raise ArithmeticError("the curve-12 change of basis produced the zero point")
        changed_points.append([str(point[0]), str(point[1])])

    transformed = {
        **public_record,
        "points": changed_points,
        "snapshot_rank_lower_bound": 29,
    }
    independence = transformed_signature_certificate(rank29, change)
    record = base["audit_curve"](transformed, independence)
    if record["residual_gain_over_mw17"] != 12:
        raise ArithmeticError("the alternate-Q80 residual gain changed")
    return relabel_curve12_record(record)


def comparison_row(record: dict[str, Any]) -> dict[str, Any]:
    jump = int(record["residual_gain_over_mw17"])
    valuation_obstruction = int(record["residual_bad_valuation_rank_modulo_mw17"])
    unit_ambiguity = int(record["norm_square_unit_squareclass_dimension"])
    lower = int(record["proved_adjusted_residual_class_group_image_dimension_lower_bound"])
    if lower != max(0, jump - valuation_obstruction - unit_ambiguity):
        raise ArithmeticError("the residual pressure formula changed")
    return {
        "curve_id": int(record["curve_id"]),
        "frame": record.get("native_frame", "published-R17"),
        "known_rank_lower_bound": int(record["point_count"]),
        "known_jump_over_generic_mw17": jump,
        "bad_valuation_rank_generic_mw17": int(
            record["generic_mw17_bad_valuation_rank"]
        ),
        "bad_valuation_rank_all_known_points": int(
            record["known_point_bad_valuation_rank"]
        ),
        "residual_bad_valuation_obstruction_dimension": valuation_obstruction,
        "norm_positive_unit_squareclass_ambiguity_dimension": unit_ambiguity,
        "known_kummer_half_ideal_span_modulo_generic_mw17_dimension_lower_bound": lower,
        "full_cubic_class_group_2rank_lower_bound": int(
            record["proved_class_group_2rank_lower_bound"]
        ),
        "formula": f"max(0,{jump}-{valuation_obstruction}-{unit_ambiguity})={lower}",
    }


def build() -> dict[str, Any]:
    base = load_base_module()
    replayed_base = base["build"]()
    stored_base = json.loads(BASE_OUTPUT.read_text())
    if replayed_base != stored_base:
        raise ArithmeticError("the stored five-fibre pressure certificate does not replay")
    if replayed_base.get("status") != base["STATUS"]:
        raise ArithmeticError("the five-fibre pressure status changed")

    curves = list(replayed_base["curves"])
    for record in curves:
        record["native_frame"] = "published-R17"
        record["native_chart"] = "norm12-orbit-074d9"
        record["ordered_half_ideal_basis_labels"] = [
            f"P{index}" for index in range(1, int(record["point_count"]) + 1)
        ]
        record["generic_mw17_basis_labels"] = [f"P{index}" for index in range(1, 18)]
        record["displayed_quotient_basis_labels"] = [
            f"P{index}" for index in range(18, int(record["point_count"]) + 1)
        ]
    curves.append(audit_curve12(base))
    if tuple(int(record["curve_id"]) for record in curves) != TARGET_IDS:
        raise ArithmeticError("the six-fibre comparison inventory changed")

    rows = [comparison_row(record) for record in curves]
    ordered = sorted(rows, key=lambda row: row["known_jump_over_generic_mw17"])
    distinct_jump_strata = []
    for jump in sorted({row["known_jump_over_generic_mw17"] for row in rows}):
        bounds = [
            row[
                "known_kummer_half_ideal_span_modulo_generic_mw17_dimension_lower_bound"
            ]
            for row in rows
            if row["known_jump_over_generic_mw17"] == jump
        ]
        distinct_jump_strata.append(
            {
                "known_jump_over_generic_mw17": jump,
                "curve_ids": [
                    row["curve_id"]
                    for row in rows
                    if row["known_jump_over_generic_mw17"] == jump
                ],
                "pressure_lower_bound_minimum": min(bounds),
                "pressure_lower_bound_maximum": max(bounds),
            }
        )
    strict_separation = all(
        left["pressure_lower_bound_maximum"]
        < right["pressure_lower_bound_minimum"]
        for left, right in zip(distinct_jump_strata, distinct_jump_strata[1:])
    )
    if not strict_separation:
        raise ArithmeticError("the observed jump strata are no longer strictly separated")
    if any(row["residual_bad_valuation_obstruction_dimension"] for row in rows):
        raise ArithmeticError("a residual direction now adds bad-valuation rank")

    return {
        "schema": SCHEMA,
        "status": STATUS,
        "summary": {
            "curve_ids": list(TARGET_IDS),
            "frames": ["published-R17", "alternate-Q80"],
            "comparison_rows": rows,
            "jump_strata": distinct_jump_strata,
            "all_exceptional_blocks_add_zero_bad_valuation_rank_modulo_generic_mw17": True,
            "strict_pressure_lower_bound_separation_between_observed_jump_strata": True,
            "finite_comparison_conclusion": (
                "On these six exact fibres the residual pressure lower-bound strata "
                "increase strictly with the known rank-jump strata: +5 maps to 3, "
                "+6 to 5, +8 to 6, and every +12 fibre lies in the final certified "
                "interval recorded above."
            ),
            "interpretation": (
                "The comparison explains why full cubic class-group computation is "
                "pressured by known exceptional points.  It is not a prospective "
                "rank predictor because those point directions are inputs."
            ),
        },
        "curves": curves,
        "theorem": {
            "spaces": (
                "Let B be the F2-space of certified known Kummer point classes, G "
                "the specialized generic MW17 subspace, v the parity map at all "
                "prime ideals above 2*Delta, and c:ker(v)->Cl(K)[2] the ideal "
                "square-root map."
            ),
            "residual_quotient": (
                "The measured object is the image induced by adjusted known classes "
                "in Cl(K)[2]/c(ker(v) intersect G)."
            ),
            "lower_bound": (
                "dim residual image >= dim(B/G) - (rank(v(B))-rank(v(G))) "
                "- (r1+r2-1)."
            ),
            "reason": (
                "After killing the residual valuation obstruction with generic "
                "classes, any remaining kernel is represented by norm-positive unit "
                "squareclasses, whose dimension is at most r1+r2-1 because K has "
                "odd degree and Norm(-1)=-1."
            ),
            "specialization_to_this_dataset": (
                "rank(v(B))-rank(v(G))=0 on all six fibres, so the lower bound is "
                "known jump minus the norm-positive unit ambiguity."
            ),
        },
        "method": {
            **replayed_base["method"],
            "curve_12_basis_change": (
                "the certified unimodular 29-by-29 matrix [specialized MW17 | "
                "P2,P11,P4,P3,P6,P8,P17,P10,P28,P24,P19,P15]"
            ),
            "comparability": (
                "all six rows use the identical completed-square cubic, all-bad-prime "
                "valuation-parity, half-ideal, and norm-positive-unit calculation"
            ),
        },
        "claim_boundary": [
            "The residual values are unconditional lower bounds for a known point-forced quotient of the full cubic ideal class-group 2-torsion; they are not exact dimensions.",
            "They are not S-class-group lower bounds after bad-prime ideals are inverted.",
            "They do not compute a Selmer upper bound, a complete residual Selmer group, or an exact elliptic-curve rank.",
            "Strict ordering is a theorem only for this six-fibre finite dataset; no population-level or out-of-sample correlation is claimed.",
            "The invariant uses the already known exceptional points and therefore explains class-group pressure but is not by itself a prospective structural predictor of a rank jump.",
        ],
        "inputs": {
            relative(path): digest(path)
            for path in (BASE_SCRIPT, BASE_OUTPUT, PUBLIC, CURVE12, RANK29, SOURCE)
        },
        "software_assumptions": {
            "sage": str(sage_version),
            "pari": ".".join(str(part) for part in pari.version()),
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_kummer_classgroup_pressure_comparison.sage --check"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    document = build()
    serialized = json.dumps(document, indent=2, sort_keys=True) + "\n"
    output = args.output.resolve()
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ArithmeticError("stored six-fibre pressure comparison differs")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    bounds = ",".join(
        f"{row['curve_id']}:{row['known_kummer_half_ideal_span_modulo_generic_mw17_dimension_lower_bound']}"
        for row in document["summary"]["comparison_rows"]
    )
    print(
        f"{PROTOCOL}|curves=6|bounds={bounds}|status={document['status']}|"
        f"output={relative(output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
