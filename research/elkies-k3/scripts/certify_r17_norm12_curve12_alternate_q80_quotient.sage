#!/usr/bin/env sage-python
"""Identify the specialized alternate-Q80 subgroup inside ICARM curve 12.

The complete ICARM sweep recognizes the Elkies--Klagsbrun rank-at-least-29
curve (ICARM id 12) as an untwisted fibre of norm12-orbit-11952.  This exact
follow-up specializes the already certified saturated generic MW17 basis,
places it in the independently certified 29-point public subgroup, verifies
all seventeen group-law relations, and computes the Smith quotient.

The displayed 29-point subgroup modulo the specialized generic subgroup is
free of rank 12.  No assertion that the displayed subgroup is the full
Mordell--Weil group, and hence no exact-rank assertion, is made.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, Matrix, PolynomialRing, QQ, ZZ, matrix
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[2]
ATLAS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
)
SWEEP = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v1.json"
)
DIRECT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
RANK29 = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/elliptic_elkies_klagsbrun_rank29_certificate.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-curve12-alternate-q80-quotient-v1.json"
)
CURVE_ID = 12
CHART = "norm12-orbit-11952"

# Columns are the seventeen specialized generic sections in the order stored
# by the direct-fibration artifact.  Rows are the ordered 29 public points.
# The matrix was discovered from a 100-decimal Neron--Tate solve whose largest
# distance to an integer was < 5e-87; every resulting relation is rechecked
# below by exact elliptic-curve group law, so the numerical solve is not part
# of the proof.
GENERIC_IN_PUBLIC_COORDINATES = (
    (1, 2, 0, 2, 1, 2, -1, 3, 1, 2, 3, 2, 2, 0, 3, 5, 2),
    (-1, -4, 0, -3, -2, -3, 1, -6, -2, -5, -5, -3, -4, -1, -3, -5, -4),
    (1, 2, 1, 1, 2, 1, 1, 3, 1, 2, 2, 1, 2, 1, 1, 1, 2),
    (0, -1, -1, -1, -3, -1, -1, -3, -1, -2, -2, -1, -2, -1, 0, -1, -2),
    (1, 1, 0, 2, 1, 2, -1, 2, 1, 1, 1, 1, 1, 0, 0, 1, 1),
    (2, 6, 1, 4, 5, 4, 1, 10, 4, 6, 6, 3, 6, 1, 4, 4, 6),
    (0, 0, 0, 1, 0, 1, -1, -1, 0, 0, 0, 1, 0, 0, 0, 2, 0),
    (1, 3, 1, 2, 3, 2, 1, 5, 2, 3, 3, 2, 3, 1, 2, 2, 3),
    (-2, -4, 0, -2, -2, -3, 0, -5, -2, -3, -4, -2, -3, 0, -3, -4, -3),
    (-2, -4, -1, -3, -3, -3, 0, -6, -2, -4, -4, -3, -4, -1, -3, -4, -4),
    (0, 1, 0, 1, 2, 1, 0, 3, 1, 2, 2, 1, 2, 0, 1, 2, 2),
    (1, 2, 1, 1, 2, 1, 1, 3, 1, 2, 2, 1, 2, 0, 2, 1, 1),
    (0, 1, 0, 1, 1, 1, 0, 2, 1, 1, 1, 1, 1, 0, 0, 0, 1),
    (-1, -2, 0, -1, -1, -1, 0, -3, -1, -2, -2, -1, -2, 0, -2, -2, -2),
    (-1, -2, 0, -1, -1, -1, 0, -3, -1, -2, -2, -1, -2, 0, -2, -2, -2),
    (4, 9, 0, 7, 5, 8, -2, 13, 5, 9, 10, 7, 9, 1, 8, 12, 9),
    (1, -1, -1, -2, -1, -1, 0, -2, -1, -2, -1, -1, -1, -1, 1, 0, -1),
    (-1, -2, 1, -2, 0, -2, 2, -3, -1, -2, -2, -2, -2, 1, -3, -4, -2),
    (1, 2, 0, 1, 1, 1, 0, 3, 1, 2, 2, 1, 2, 0, 2, 2, 2),
    (-1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, -2, -2, 1),
    (-1, -5, -1, -3, -4, -3, -1, -8, -3, -6, -6, -3, -5, -2, -3, -4, -5),
    (0, -1, -1, -1, -2, -1, -1, -2, -1, -1, -1, -1, -1, -1, 0, 0, -1),
    (-2, -4, 0, -3, -2, -3, 1, -6, -2, -4, -4, -3, -4, 0, -4, -5, -4),
    (-2, -4, 0, -3, -2, -3, 1, -6, -2, -4, -4, -3, -4, 0, -4, -5, -4),
    (2, 5, 1, 3, 3, 4, 0, 7, 3, 5, 5, 3, 5, 1, 4, 5, 5),
    (2, 5, 1, 4, 4, 4, 0, 8, 3, 5, 5, 3, 5, 1, 3, 4, 5),
    (0, -4, -1, -4, -3, -3, 0, -7, -3, -5, -4, -3, -4, -1, -2, -3, -4),
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (1, 2, 0, 1, 0, 1, 0, 2, 1, 1, 1, 1, 1, 0, 2, 1, 1),
)
QUOTIENT_BASIS_PUBLIC_POINT_INDICES = (2, 11, 4, 3, 6, 8, 17, 10, 28, 24, 19, 15)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial(ring, values):
    return ring([QQ(value) for value in values])


def rational_function(ring, record):
    return ring.fraction_field()(
        polynomial(ring, record["numerator_coefficients_low_to_high"])
        / polynomial(ring, record["denominator_coefficients_low_to_high"])
    )


def short_invariants(ainvs):
    a1, a2, a3, a4, a6 = map(QQ, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    return a1, a3, b2, -c4 / 48, -c6 / 864


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    args.output = args.output.resolve()

    atlas = json.loads(ATLAS.read_text())
    sweep = json.loads(SWEEP.read_text())
    direct = json.loads(DIRECT.read_text())
    rank29 = json.loads(RANK29.read_text())
    if rank29["status"] != "exact_unconditional_rank_at_least_29":
        raise ArithmeticError("the curve-12 rank lower-bound input is not exact")
    if rank29["finite_reduction_certificate"]["combined_exact_rank_over_F2"] != 29:
        raise ArithmeticError("the 29 public points lost their independence certificate")
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError("the alternate-Q80 generic basis is no longer saturated")

    hit = next(
        record
        for record in sweep["rational_j_hits_and_twists"]
        if record["curve_id"] == CURVE_ID and record["representative"] == CHART
    )
    if hit["snapshot_rank_lower_bound"] != 29:
        raise ArithmeticError("curve 12 lost its pinned rank-at-least-29 status")
    if any(
        record["twist"]["status"] != "QQ_ISOMORPHIC_UNTWISTED"
        for record in hit["native_chart_twists"]
    ):
        raise ArithmeticError("curve 12 acquired a nontrivial native-chart twist")
    numerator = ZZ(hit["representative_parameter"]["numerator"])
    denominator = ZZ(hit["representative_parameter"]["denominator"])
    parameter = QQ(numerator) / denominator

    ring = PolynomialRing(QQ, "u")
    direct_A = polynomial(
        ring, direct["weierstrass_model"]["A_coefficients_low_to_high"]
    )
    direct_B = polynomial(
        ring, direct["weierstrass_model"]["B_coefficients_low_to_high"]
    )
    atlas_chart = next(
        record for record in atlas["atlas"]["charts"] if record["label"] == CHART
    )
    atlas_A = polynomial(
        ring, atlas_chart["weierstrass_model"]["A_coefficients_low_to_high"]
    )
    atlas_B = polynomial(
        ring, atlas_chart["weierstrass_model"]["B_coefficients_low_to_high"]
    )
    if direct_A != atlas_A or direct_B != atlas_B:
        raise ArithmeticError("the direct and six-class atlas 11952 models differ")

    ainvs = rank29["general_weierstrass_coefficients"]
    a1, a3, b2, target_A, target_B = short_invariants(ainvs)
    fibre_A = direct_A(parameter)
    fibre_B = direct_B(parameter)
    q_value = target_B * fibre_A / (fibre_B * target_A)
    if not q_value.is_square():
        raise ArithmeticError("the representative curve-12 twist is nontrivial")
    s_value = q_value.sqrt()
    if target_A != q_value**2 * fibre_A or target_B != q_value**3 * fibre_B:
        raise ArithmeticError("the curve-12 short-model isomorphism failed")

    target_curve = EllipticCurve(QQ, [target_A, target_B])
    public_points = []
    for record in rank29["published_points"]:
        x_value = QQ(record["x"])
        y_value = QQ(record["y"])
        public_points.append(
            target_curve(
                x_value + b2 / 12,
                y_value + (a1 * x_value + a3) / 2,
            )
        )
    if len(public_points) != 29:
        raise ArithmeticError("the public point inventory changed")

    generic_points = []
    for record in direct["sections"]["records"]:
        x_coordinate = rational_function(ring, record["X"])(parameter)
        y_coordinate = rational_function(ring, record["Y"])(parameter)
        generic_points.append(
            target_curve(
                q_value * x_coordinate,
                s_value**3 * y_coordinate,
            )
        )
    if len(generic_points) != 17 or len(set(generic_points)) != 17:
        raise ArithmeticError("the specialized generic point inventory changed")

    coordinates = matrix(ZZ, GENERIC_IN_PUBLIC_COORDINATES)
    if coordinates.dimensions() != (29, 17) or coordinates.rank() != 17:
        raise ArithmeticError("the generic-in-public coordinate matrix lost rank")
    for column, generic_point in enumerate(generic_points):
        reconstructed = target_curve(0)
        for coefficient, public_point in zip(
            coordinates.column(column), public_points
        ):
            reconstructed += int(coefficient) * public_point
        if reconstructed != generic_point:
            raise ArithmeticError(f"generic section {column + 1} relation failed")

    smith = coordinates.elementary_divisors()
    if smith != [ZZ(1)] * 17 + [ZZ(0)] * 12:
        raise ArithmeticError("the displayed exceptional quotient is not free rank 12")
    complement = matrix(
        ZZ,
        29,
        12,
        lambda row, column: int(
            row + 1 == QUOTIENT_BASIS_PUBLIC_POINT_INDICES[column]
        ),
    )
    augmented = coordinates.augment(complement)
    if abs(augmented.det()) != 1:
        raise ArithmeticError("the displayed quotient basis is not unimodular")

    payload = {
        "schema": "elkies-k3.r17-norm12-curve12-alternate-q80-quotient.v1",
        "status": "PROVED_CURVE12_NATIVE_ALTERNATE_Q80_AND_DISPLAYED_QUOTIENT",
        "curve": {
            "icarm_id": CURVE_ID,
            "unconditional_rank_lower_bound": 29,
            "displayed_independent_point_count": 29,
            "exact_rank_not_asserted": True,
        },
        "native_fibre": {
            "chart": CHART,
            "frame_class": "alternate-Q80",
            "parameter": rational_text(parameter),
            "quadratic_twist": "trivial",
            "short_model_scale_q": rational_text(q_value),
            "short_model_scale_s_with_s_squared_q": rational_text(s_value),
            "isomorphism": "x_target=q*x_fibre; y_target=s^3*y_fibre",
        },
        "specialized_generic_subgroup": {
            "rank": 17,
            "source": "the saturated generic basis of norm12-orbit-11952",
            "coordinate_matrix_rows_in_ordered_29_public_points": [
                list(map(int, row)) for row in coordinates.rows()
            ],
            "all_seventeen_relations_verified_by_exact_group_law": True,
        },
        "displayed_exceptional_quotient": {
            "quotient": "Z^12",
            "smith_nonzero_invariant_factors": [1] * 17,
            "smith_zero_count": 12,
            "free_basis_modulo_specialized_generic": [
                f"P{index}" for index in QUOTIENT_BASIS_PUBLIC_POINT_INDICES
            ],
            "augmented_basis_determinant": int(augmented.det()),
            "scope": (
                "quotient of the subgroup generated by the 29 displayed public points; "
                "not a quotient of a proved full Mordell-Weil group"
            ),
        },
        "discovery_note": {
            "method": "100-decimal Neron--Tate linear solve followed by exact recovery",
            "maximum_distance_to_recovered_integer": "<5e-87",
            "proof_uses_numerics": False,
        },
        "claim_boundary": {
            "proved": [
                "ICARM curve 12 is a QQ-isomorphic untwisted fibre of the native norm12-orbit-11952 alternate-Q80 family",
                "the seventeen saturated generic sections specialize injectively into the displayed 29-point subgroup",
                "the displayed subgroup modulo the specialized generic subgroup is free of rank twelve with the stated basis",
            ],
            "not_proved": [
                "that the displayed 29-point subgroup is the full Mordell-Weil group",
                "an unconditional rank upper bound or exact rank 29",
            ],
        },
        "inputs": {
            relative(ATLAS): digest(ATLAS),
            relative(SWEEP): digest(SWEEP),
            relative(DIRECT): digest(DIRECT),
            relative(RANK29): digest(RANK29),
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact QQ elliptic-curve group law",
                "exact Smith normal form",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "certify_r17_norm12_curve12_alternate_q80_quotient.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored curve-12 quotient differs from exact replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17CURVE12Q80|curve=12|rank_lower_bound=29|chart=11952|twist=trivial|"
        "generic_rank=17|displayed_rank=29|quotient=Z^12|status=PROVED|output={}".format(
            relative(args.output)
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
