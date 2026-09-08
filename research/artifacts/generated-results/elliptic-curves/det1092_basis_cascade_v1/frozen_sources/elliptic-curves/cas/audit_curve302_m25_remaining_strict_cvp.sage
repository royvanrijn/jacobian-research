#!/usr/bin/env sage-python
"""Retrospective M25 re-translation of the six remaining strict 302 charts.

This is a narrow diagnostic, not an exhaustive M25/2M25 scan.  For each of
the six residual strict directions, it takes the previously declared
minimum-coordinate chart among that direction's fixed 32 M17-vetted charts.
The same M17 centre is padded into the certified rank-25 subgroup, and a
cross-precision CVP chooses its best M25 translate in the rounded 384-bit
canonical-height metric.  The resulting pointed quartic coordinate is then
checked exactly.

Known target points enter only this retrospective diagnostic.  Nothing here
is an execution selector or a point search.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import json
from decimal import Decimal, getcontext
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path

from sage.all import EllipticCurve, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
VISIBILITY = ART / "curve302_residual_visibility_geometry_v1.json"
M25 = ART / "curve302_m24_bisection_orbit117420_mod2_v1.json"
PARENT = ART / "curve302_recovered_mw17_parent_v1.json"
MAPPER = CAS / "factor_free_pari_mapping.sage"
CVP_SOURCE = CAS / "audit_curve302_residual_visibility_geometry.sage"
OUTPUT = ART / "curve302_m25_remaining_strict_cvp_v1.json"
SCALE = 1_000_000
M17 = 17
M25_DIMENSION = 25
PRECISION = 192


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load(path: Path):
    return json.loads(path.read_text())


def qtext(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def as_ints(values) -> list[int]:
    return [int(value) for value in values]


def primitive(a, b):
    a, b = Fraction(a), Fraction(b)
    denominator = a.denominator * b.denominator // gcd(a.denominator, b.denominator)
    left, right = int(a * denominator), int(b * denominator)
    common = gcd(abs(left), abs(right))
    if not common:
        raise ArithmeticError("zero projective coordinate")
    left //= common
    right //= common
    if right < 0:
        left, right = -left, -right
    return left, right


def square_root(value: Fraction) -> Fraction:
    if value < 0:
        raise ArithmeticError("negative quartic value")
    numerator, denominator = isqrt(value.numerator), isqrt(value.denominator)
    if numerator * numerator != value.numerator or denominator * denominator != value.denominator:
        raise ArithmeticError("quartic value is not a square")
    return Fraction(numerator, denominator)


def rounded_gram(curve, points):
    getcontext().prec = 110
    pari.default("realprecision", 110)
    pari.allocatemem(256_000_000, 1_073_741_824, silent=True)
    height = pari(curve).ellheightmatrix([list(point) for point in points], precision=384)
    decimal = [[Decimal(str(height[i, j])) for j in range(len(points))] for i in range(len(points))]
    asymmetry = max(abs(decimal[i][j] - decimal[j][i]) for i in range(len(points)) for j in range(len(points)))
    if asymmetry > Decimal("1e-90"):
        raise ArithmeticError("PARI height matrix is asymmetric")
    answer = matrix(ZZ, len(points), len(points), [int((entry * Decimal(SCALE)).to_integral_value()) for row in decimal for entry in row])
    if not answer.is_positive_definite():
        raise ArithmeticError("rounded M25 target height matrix is not positive definite")
    return answer, str(asymmetry)


def coordinate(mapper, model, basis, centre_word, target):
    """Return an exact reduced coordinate and a checked quartic identity."""

    from alternate_quartic_covers import alternate_cover
    from half_lattice_pointed_sieve import linear_combination

    centre = linear_combination(model, basis, centre_word)
    if centre is None:
        raise ArithmeticError("M25 centre is at infinity")
    mapped = mapper.mapping(model, basis, {"representative": as_ints(centre_word)})
    raw_t, raw_w = alternate_cover(model, centre).curve_point_to_cover(target)
    a, b, c, d = [Fraction(value) for value in mapped["matrix"]]
    numerator, denominator = primitive(d * raw_t - b, -c * raw_t + a)
    coefficients = [Fraction(value) for value in mapped["discriminant_quartic"]]
    value = sum(coefficient * numerator**index * denominator ** (4 - index) for index, coefficient in enumerate(coefficients))
    transformed_denominator = c * numerator + d * denominator
    ratio = Fraction(mapped["square_ratio"])
    root = abs(raw_w * transformed_denominator**2 / square_root(ratio))
    if value != root * root:
        raise ArithmeticError("exact reduced quartic identity failed")
    return {
        "centre_point_short_model": [qtext(value) for value in centre],
        "mapping": mapped,
        "reduced_coordinate": [str(numerator), str(denominator)],
        "reduced_coordinate_height": str(max(abs(numerator), denominator)),
        "reduced_coordinate_decimal_digits": len(str(max(abs(numerator), denominator))),
        "reduced_quartic_square_root_absolute": str(root),
    }


def m25_translate(cvp, gram, change, centre, target_cross):
    """Choose R+m in the fixed centre chart by rounded-metric CVP."""

    inverse = change.inverse()
    reduced = change * gram * change.transpose()
    reduced_centre = vector(ZZ, (matrix(ZZ, 1, M25_DIMENSION, centre) * inverse).row(0))
    parity = vector(ZZ, [int(value) % 2 for value in reduced_centre])
    u = vector(ZZ, [(reduced_centre[index] - parity[index]) // 2 for index in range(M25_DIMENSION)])
    original_projection = gram.solve_right(target_cross)
    projection = vector(QQ, (matrix(QQ, 1, M25_DIMENSION, original_projection).row(0) * inverse))
    dd = cvp.RoundedMetricCVP(reduced, "dd", PRECISION)
    mpfr = cvp.RoundedMetricCVP(reduced, "mpfr", PRECISION)
    z, distance = cvp.cross_precision_closest(dd, mpfr, vector(QQ, [(projection[i] - parity[i]) / 2 for i in range(M25_DIMENSION)]))
    translated = vector(ZZ, (matrix(ZZ, 1, M25_DIMENSION, u - z) * change).row(0))
    nearest_double = vector(ZZ, (matrix(ZZ, 1, M25_DIMENSION, parity + 2 * z) * change).row(0))
    return {
        "centre_M25_word": as_ints(centre),
        "M25_parity_mask": sum((int(value) & 1) << index for index, value in enumerate(centre)),
        "closest_to_double_target_M25_word": as_ints(nearest_double),
        "target_translation_M25_word": as_ints(translated),
        "CVP_distance_in_reduced_metric": str(distance),
    }


def build() -> dict:
    visibility, m25, parent = load(VISIBILITY), load(M25), load(PARENT)
    if visibility["status"] != "PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC":
        raise ArithmeticError("M17 residual visibility diagnostic is unavailable")
    if m25["status"] != "COMPLETE_DECLARED_FINITE_AUDIT" or m25["rank_lower_bound"] != M25_DIMENSION:
        raise ArithmeticError("rank-25 M24-to-M25 cloud is unavailable")
    basis = tuple(tuple(QQ(value) for value in point) for point in m25["independent_points"])
    short_model = tuple(Fraction(value) for value in m25["curve"])
    map_basis = tuple(tuple(Fraction(value) for value in point) for point in m25["independent_points"])
    if len(basis) != M25_DIMENSION:
        raise ArithmeticError("M25 basis input changed")
    curve = EllipticCurve(QQ, [QQ(value) for value in m25["curve"]])
    if any(point not in curve for point in basis):
        raise ArithmeticError("M25 basis contains a point off the short curve")
    mapper = importlib.machinery.SourceFileLoader("curve302_m25_strict_mapper", str(MAPPER)).load_module()
    cvp = importlib.machinery.SourceFileLoader("curve302_m25_strict_cvp", str(CVP_SOURCE)).load_module()
    targets = [row for row in visibility["directions"] if row["cohort"] == "residual_strict" and row["id"] != "residual-strict-06"]
    if len(targets) != 6:
        raise ArithmeticError("the already recovered residual direction changed")
    target_points = []
    chosen = []
    for row in targets:
        parity = row["minimum_finite_coordinate_parity_in_vetted_set"]
        chart = next(item for item in row["exact_pointed_quartic_rows"] if item["parity"] == parity)
        target = tuple(QQ(value) for value in chart["exact_target_point_short_model"])
        if target not in curve:
            raise ArithmeticError("retrospective target is off the M25 curve")
        target_points.append(target)
        chosen.append((row, chart, target, tuple(Fraction(value) for value in chart["exact_target_point_short_model"])))
    height, asymmetry = rounded_gram(curve, (*basis, *target_points))
    gram = height[:M25_DIMENSION, :M25_DIMENSION]
    if not gram.is_positive_definite():
        raise ArithmeticError("M25 rounded Gram is not positive definite")
    change = matrix(ZZ, pari(gram).qflllgram()).transpose()
    if abs(change.det()) != 1:
        raise ArithmeticError("M25 LLL change is not unimodular")
    rows = []
    from alternate_quartic_covers import short_add
    from half_lattice_pointed_sieve import linear_combination
    for index, (entry, chart, target, target_for_map) in enumerate(chosen):
        centre = vector(ZZ, chart["centre_m17_word"] + [0] * (M25_DIMENSION - M17))
        target_cross = vector(ZZ, [height[row, M25_DIMENSION + index] * 2 for row in range(M25_DIMENSION)])
        translation = m25_translate(cvp, gram, change, centre, target_cross)
        added = linear_combination(short_model, map_basis, translation["target_translation_M25_word"])
        moved_target = target_for_map if added is None else short_add(short_model, target_for_map, added)
        if moved_target is None:
            raise ArithmeticError("M25 translation sent a target to infinity")
        exact = coordinate(mapper, short_model, map_basis, centre, moved_target)
        old_height = int(entry["minimum_finite_reduced_coordinate_height_in_vetted_set"])
        old_center = chart["centre_m17_word"]
        rows.append({
            "id": entry["id"],
            "old_M17_vetted_minimum": {
                "coordinate_height": str(old_height),
                "coordinate_decimal_digits": len(str(old_height)),
                "parity": chart["parity"],
                "centre_M17_word": old_center,
                "centre_l1": sum(abs(value) for value in old_center),
                "centre_linf": max(abs(value) for value in old_center),
                "degree_two_translation_orbit": chart["degree_two_translation_orbit"],
                "old_initial_SHA_sample_contains_parity": False,
                "old_initial_selected_centres_contain_parity": False,
            },
            "M25_fixed_chart_CVP": {
                **translation,
                "translated_target_short_model": [qtext(value) for value in moved_target],
                **exact,
                "improves_declared_M17_vetted_height": int(exact["reduced_coordinate_height"]) < old_height,
                "inside_historical_height_125000": int(exact["reduced_coordinate_height"]) <= 125000,
            },
        })
    return {
        "schema": "elliptic-curves.curve302-m25-remaining-strict-cvp.v1",
        "status": "PASS_RETROSPECTIVE_FIXED_CHART_M25_CVP",
        "inputs": {rel(path): sha(path) for path in (VISIBILITY, M25, PARENT, MAPPER, CVP_SOURCE, Path(__file__))},
        "protocol": {
            "target_count": 6,
            "chart_rule": "for each remaining direction, reuse the old declared minimum-coordinate chart from its vetted set of 32 and only retranslate its target through M25",
            "M25_CVP": "384-bit PARI height pairings rounded at 10^6; DD/MPFR-agreeing CVP in the rank-25 rounded lattice",
            "boundary": "No all-parity M25 scan, global coordinate minimum, point search, or prospective selector claim.",
        },
        "M25_rounded_height_gram": [as_ints(row) for row in gram.rows()],
        "M25_LLL_change": [as_ints(row) for row in change.rows()],
        "maximum_height_matrix_asymmetry": asymmetry,
        "directions": rows,
        "reproducing_command": "sage -python elliptic-curves/cas/audit_curve302_m25_remaining_strict_cvp.sage --check",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.build == arguments.check:
        parser.error("choose exactly one of --build or --check")
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.build:
        if OUTPUT.exists():
            raise FileExistsError("preserve immutable M25 residual CVP diagnostic")
        OUTPUT.write_text(rendered)
    elif OUTPUT.read_text() != rendered:
        raise ArithmeticError("stored M25 residual CVP diagnostic did not replay")
    print("CURVE302M25STRICTCVP|directions=6|status=PASS_RETROSPECTIVE_FIXED_CHART_M25_CVP", flush=True)


if __name__ == "__main__":
    main()
