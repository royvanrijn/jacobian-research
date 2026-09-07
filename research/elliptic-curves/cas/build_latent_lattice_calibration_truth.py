#!/usr/bin/env python3
"""Build withheld exact truth for latent-lattice control calibration.

The output is not an inference result.  It supplies exact embeddings that a
blind selector is scored against after selection:

* the published rank-17 subgroup in each rank-25--28 public Elkies fibre;
* the known rank-12 Fermigier--Mestre subgroup in ICARM 245.

Coordinates are proposed by a high-precision height-dual solve and accepted
only after exact elliptic group-law replay.  Smith factors state whether the
resulting embedding is primitive in the displayed subgroup.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import importlib
import json
from pathlib import Path
import platform
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
)
from ecsearch.fermigier_rank import specialize_fermigier_rank_sections  # noqa: E402
from elliptic_candidate_record import source_point_to_target  # noqa: E402
from icarm_curve245 import (  # noqa: E402
    GENERAL_WEIERSTRASS_COEFFICIENTS as CURVE245_MODEL,
    POINTS as CURVE245_POINTS,
)
from icarm_curve245_mestre import (  # noqa: E402
    ANCHOR_SHORT_TO_PUBLIC_CHANGE,
    CANONICAL_PARAMETER,
    CONSTRUCTION,
    PUBLIC_MODEL,
    extra_quartic_point,
    primitive_short_model,
)
from nagao_1994 import (  # noqa: E402
    primitive_visible_points,
    quartic_point_to_short_jacobian,
)
from elliptic_candidate_record import (  # noqa: E402
    WeierstrassChange,
    change_weierstrass_model,
)
from latent_lattice import (  # noqa: E402
    EllipticCurve,
    height_gram,
    primitive_column_closure,
    recover_exact_embedding,
)


MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
PINNED_R17 = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
PUBLISHED_R17_CERTIFICATE = (
    ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "latent_lattice_calibration_truth_v1.json"
)
CURVE282_SOURCE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_7fff_zip_public_source_281_282_285_286.json"
)
FERMIGIER_RANK20_SOURCE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)
POSITIVE_PARAMETERS = (
    (25, -2, 377),
    (26, -308, 251),
    (27, 2456, 135),
    (28, -9529, 5471),
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def float_gram(curve: EllipticCurve, points, digits: int) -> np.ndarray:
    return np.array(height_gram(curve, points, digits=digits), dtype=float)


def gram_fit(left: np.ndarray, right: np.ndarray) -> dict[str, str]:
    scale = float(np.sum(left * right) / np.sum(left * left))
    residual = float(np.linalg.norm(right - scale * left) / np.linalg.norm(right))
    correlation = float(np.corrcoef(left.ravel(), right.ravel())[0, 1])
    return {
        "right_over_left_scale": f"{scale:.17g}",
        "relative_frobenius_residual": f"{residual:.17g}",
        "entry_correlation": f"{correlation:.17g}",
    }


def positive_controls(digits: int) -> tuple[list[dict[str, object]], list[np.ndarray]]:
    data = load_q12o5867_data(MODEL, SECTIONS)
    records = []
    grams = []
    for rank, numerator, denominator in POSITIVE_PARAMETERS:
        module = importlib.import_module(f"elkies_rank{rank}")
        displayed_points = tuple(module.POINTS)
        specialization = evaluate_projective_specialization(
            data, numerator, denominator
        )
        minimal_model, minimal_change, _metadata = global_minimal_model_with_change(
            specialization.model
        )
        if tuple(minimal_model) != tuple(module.GENERAL_WEIERSTRASS_COEFFICIENTS):
            raise ArithmeticError("a positive control no longer matches its public model")
        generic_points = tuple(
            source_point_to_target(point, minimal_change)
            for point in specialization.points
        )
        curve = EllipticCurve(minimal_model)
        embedding = recover_exact_embedding(
            curve,
            displayed_points,
            generic_points,
            digits=max(120, digits),
        )
        gram = float_gram(curve, generic_points, digits)
        grams.append(gram)
        smith = embedding.smith_invariant_factors()
        records.append(
            {
                "label": f"rank_at_least_{rank}",
                "parameter": f"{numerator}/{denominator}",
                "displayed_rank": rank,
                "truth_rank": 17,
                "embedding_matrix_rows": [list(row) for row in embedding.rows()],
                "embedding_matrix_columns": [
                    list(column) for column in embedding.columns
                ],
                "smith_invariant_factors": list(smith),
                "primitive_in_displayed_subgroup": all(value == 1 for value in smith),
                "max_abs_coordinate": embedding.max_abs_coordinate,
                "nonzero_coordinate_count": embedding.nonzero_coordinates,
                "height_dual_numerical_residual_max": embedding.numerical_residual_max,
                "exact_group_law_replay": True,
                "canonical_height_gram": [list(row) for row in height_gram(curve, generic_points, digits=digits)],
            }
        )
    return records, grams


def negative_control(digits: int) -> dict[str, object]:
    parameter = CANONICAL_PARAMETER
    change = WeierstrassChange.from_values(ANCHOR_SHORT_TO_PUBLIC_CHANGE)
    if change_weierstrass_model(primitive_short_model(parameter), change) != PUBLIC_MODEL:
        raise ArithmeticError("the curve-245 family transport changed")
    quartic_points = list(primitive_visible_points(CONSTRUCTION, parameter))
    quartic_points.append(extra_quartic_point(parameter))
    short_points = tuple(
        quartic_point_to_short_jacobian(CONSTRUCTION, parameter, point)
        for point in quartic_points
    )
    generic_points_with_relation = tuple(
        source_point_to_target(point, change) for point in short_points
    )
    # The first twelve visible points sum to zero.  Retain eleven plus the
    # independent Fermigier section in position thirteen.
    generic_basis = generic_points_with_relation[:11] + (
        generic_points_with_relation[12],
    )
    curve = EllipticCurve(CURVE245_MODEL)
    embedding = recover_exact_embedding(
        curve,
        CURVE245_POINTS,
        generic_basis,
        digits=max(150, digits),
    )
    smith = embedding.smith_invariant_factors()
    primitive_closure = primitive_column_closure(embedding.rows())
    closure_columns = tuple(zip(*primitive_closure))
    closure_points = tuple(
        curve.linear_combination(CURVE245_POINTS, column)
        for column in closure_columns
    )
    if any(point is None or not curve.is_on_curve(point) for point in closure_points):
        raise ArithmeticError("a primitive-closure generator failed exact replay")
    if curve.linear_combination(generic_points_with_relation[:12], (1,) * 12) is not None:
        raise ArithmeticError("the Fermigier visible-section relation changed")
    return {
        "label": "ICARM_245_Fermigier_negative_control",
        "displayed_rank": 20,
        "truth_rank": 12,
        "generic_point_selection_one_based": list(range(1, 12)) + [13],
        "visible_point_relation": [1] * 12 + [0],
        "embedding_matrix_rows": [list(row) for row in embedding.rows()],
        "embedding_matrix_columns": [list(column) for column in embedding.columns],
        "primitive_closure_embedding_matrix_rows": [
            list(row) for row in primitive_closure
        ],
        "primitive_closure_points": [
            [str(point[0]), str(point[1])] for point in closure_points
        ],
        "generic_subgroup_index_in_primitive_closure": int(
            np.prod(np.array(smith, dtype=object))
        ),
        "smith_invariant_factors": list(smith),
        "primitive_in_displayed_subgroup": all(value == 1 for value in smith),
        "max_abs_coordinate": embedding.max_abs_coordinate,
        "nonzero_coordinate_count": embedding.nonzero_coordinates,
        "height_dual_numerical_residual_max": embedding.numerical_residual_max,
        "exact_group_law_replay": True,
        "canonical_height_gram": [
            list(row) for row in height_gram(curve, generic_basis, digits=digits)
        ],
    }


def _point_from_xy(record: dict[str, str]):
    return Fraction(record["x"]), Fraction(record["y"])


def curve282_control(digits: int) -> dict[str, object]:
    document = json.loads(CURVE282_SOURCE.read_text())
    source = next(record for record in document["curves"] if int(record["id"]) == 282)
    model = tuple(Fraction(value) for value in source["ainvs"])
    displayed = tuple(
        (Fraction(x_value), Fraction(y_value)) for x_value, y_value in source["points"]
    )
    parameter = Fraction(11671, 42)
    specialization = specialize_fermigier_rank_sections(parameter)
    canonical = specialization.canonical_model
    u_value = Fraction(1, 882)
    s_value = Fraction(-1, 2)
    r_value = (u_value * u_value - canonical[1] - Fraction(1, 4)) / 3
    change = WeierstrassChange(
        u_value, r_value, s_value, -(1 + r_value) / 2
    )
    if change_weierstrass_model(canonical, change) != model:
        raise ArithmeticError("the curve-282 Fermigier transport changed")
    generic = tuple(
        source_point_to_target(point, change)
        for point in specialization.section_differences
    )
    curve = EllipticCurve(model)
    embedding = recover_exact_embedding(
        curve, displayed, generic, digits=max(120, digits)
    )
    smith = embedding.smith_invariant_factors()
    return {
        "label": "ICARM_282_Fermigier_sibling",
        "parameter": "11671/42",
        "displayed_rank": len(displayed),
        "truth_rank": 12,
        "embedding_matrix_rows": [list(row) for row in embedding.rows()],
        "embedding_matrix_columns": [list(column) for column in embedding.columns],
        "smith_invariant_factors": list(smith),
        "primitive_in_displayed_subgroup": all(value == 1 for value in smith),
        "max_abs_coordinate": embedding.max_abs_coordinate,
        "nonzero_coordinate_count": embedding.nonzero_coordinates,
        "height_dual_numerical_residual_max": embedding.numerical_residual_max,
        "exact_group_law_replay": True,
        "canonical_height_gram": [
            list(row) for row in height_gram(curve, generic, digits=digits)
        ],
    }


def fermigier_rank20_control(digits: int) -> dict[str, object]:
    document = json.loads(FERMIGIER_RANK20_SOURCE.read_text())
    model = tuple(
        Fraction(value)
        for value in document["models"]["global_minimal"]["coefficients"]
    )
    selected = document["imported_selected_twenty_basis"]["basis"]
    displayed = tuple(
        _point_from_xy(record["points"]["global_minimal"]) for record in selected
    )
    generic = displayed[:12]
    curve = EllipticCurve(model)
    embedding = recover_exact_embedding(
        curve, displayed, generic, digits=max(120, digits)
    )
    smith = embedding.smith_invariant_factors()
    return {
        "label": "Fermigier_u_28917_over_20_sibling",
        "parameter": "28917/20",
        "displayed_rank": len(displayed),
        "truth_rank": 12,
        "embedding_matrix_rows": [list(row) for row in embedding.rows()],
        "embedding_matrix_columns": [list(column) for column in embedding.columns],
        "smith_invariant_factors": list(smith),
        "primitive_in_displayed_subgroup": all(value == 1 for value in smith),
        "max_abs_coordinate": embedding.max_abs_coordinate,
        "nonzero_coordinate_count": embedding.nonzero_coordinates,
        "height_dual_numerical_residual_max": embedding.numerical_residual_max,
        "exact_group_law_replay": True,
        "canonical_height_gram": [
            list(row) for row in height_gram(curve, generic, digits=digits)
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=100)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.digits < 60:
        raise SystemExit("--digits must be at least 60")
    positive, grams = positive_controls(args.digits)
    pairwise = {}
    for left in range(len(positive)):
        for right in range(left + 1, len(positive)):
            key = f"{positive[left]['label']}/{positive[right]['label']}"
            pairwise[key] = gram_fit(grams[left], grams[right])
    negative = negative_control(args.digits)
    fermigier_controls = [
        negative,
        curve282_control(args.digits),
        fermigier_rank20_control(args.digits),
    ]
    fermigier_grams = [
        np.array(record["canonical_height_gram"], dtype=float)
        for record in fermigier_controls
    ]
    fermigier_pairwise = {}
    for left in range(len(fermigier_controls)):
        for right in range(left + 1, len(fermigier_controls)):
            key = f"{fermigier_controls[left]['label']}/{fermigier_controls[right]['label']}"
            fermigier_pairwise[key] = gram_fit(
                fermigier_grams[left], fermigier_grams[right]
            )
    source_paths = [
        Path(__file__).resolve(),
        MODEL,
        SECTIONS,
        PINNED_R17,
        PUBLISHED_R17_CERTIFICATE,
        *(ELLIPTIC / "cas" / f"elkies_rank{rank}.py" for rank in range(25, 29)),
        ELLIPTIC / "cas/icarm_curve245.py",
        ELLIPTIC / "cas/icarm_curve245_mestre.py",
        CURVE282_SOURCE,
        FERMIGIER_RANK20_SOURCE,
    ]
    payload = {
        "schema": "elliptic-curves.latent-lattice-calibration-truth.v1",
        "status": "PASS_EXACT_WITHHELD_CALIBRATION_TRUTH",
        "role": (
            "Ground truth used only after blind selection. Supplying these points "
            "or matrices to the selector would invalidate calibration."
        ),
        "decimal_precision_digits": args.digits,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "inputs": {display(path): digest(path) for path in source_paths},
        "positive_controls": positive,
        "positive_truth_pairwise_height_fits": pairwise,
        "negative_control": negative,
        "fermigier_family_controls": fermigier_controls,
        "fermigier_truth_pairwise_height_fits": fermigier_pairwise,
        "proof_boundary": (
            "Point identities, coordinate embeddings, Smith factors, and ranks of the "
            "displayed coordinate matrices are exact. Canonical-height Grams and their "
            "cross-fibre fits are numerical at the declared precision. Primitivity is "
            "only relative to each displayed independent subgroup."
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != rendered:
            raise SystemExit(f"FAIL: {args.output} differs from recomputation")
        print(f"PASS|{args.output}|sha256={sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        "LATENTTRUTH|positive_ranks=17,17,17,17|negative_rank=12|"
        f"output={args.output}|sha256={sha256(rendered.encode()).hexdigest()}"
    )


if __name__ == "__main__":
    main()
