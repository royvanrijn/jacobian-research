#!/usr/bin/env python3
"""Certify a real double-shear eight-moment point and its mu_9 obstruction.

Work on the direct internal-gauge quotient with pivot rows 0 and 4,

    U_0 = e_0 + a e_1,       U_1 = e_4 + d e_3,

and the dense off-diagonal factor torus normalized by B[0,4]=B[1,0]=1.
The first moment eliminates d=-a*b01/b13.  Exact rational Krawczyk
arithmetic isolates a real zero of mu_2,...,mu_8 in a displayed rational
box.  Every torus coordinate stays nonzero, the pivot minor is -1, and
interval evaluation proves mu_9>0 throughout the box.

This certifies one finite-prefix exact-rank-two point.  It is not a global
exclusion of the double-shear chart or of the other double-shear orbits.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, factorial
from pathlib import Path

import sympy as sp

from verify_two_pair_sic_bidegree33_rank_two_finite_prefix import (
    Interval,
    as_interval,
    derivative,
    evaluate,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_rank_two_double_shear_real_prefix.json"
)
CENTER_STRINGS = (
    "-0.3524378371261593",
    "-1.3307525332074560",
    "-1.6245630868373802",
    "-0.7928610509028164",
    "0.3680936242580407",
    "0.9407113877550523",
    "1.2756140806891972",
)
RADIUS = Fraction(1, 10_000_000_000)


def multiply(left: dict[int, sp.Expr], right: dict[int, sp.Expr]) -> dict[int, sp.Expr]:
    answer: dict[int, sp.Expr] = {}
    for left_degree, left_value in left.items():
        for right_degree, right_value in right.items():
            degree = left_degree + right_degree
            answer[degree] = answer.get(degree, 0) + left_value * right_value
    return {degree: sp.expand(value) for degree, value in answer.items()}


def power(polynomial: dict[int, sp.Expr], exponent: int) -> dict[int, sp.Expr]:
    answer = {0: sp.Integer(1)}
    for _ in range(exponent):
        answer = multiply(answer, polynomial)
    return answer


def raw_moment(
    left_0: dict[int, sp.Expr],
    left_1: dict[int, sp.Expr],
    right_0: dict[int, sp.Expr],
    right_1: dict[int, sp.Expr],
    order: int,
) -> sp.Expr:
    answer = 0
    for channel_count in range(order + 1):
        left = multiply(
            power(left_0, channel_count),
            power(left_1, order - channel_count),
        )
        right = multiply(
            power(right_0, channel_count),
            power(right_1, order - channel_count),
        )
        answer += comb(order, channel_count) * sum(
            factorial(degree)
            * factorial(4 * order - degree)
            * coefficient
            * right.get(degree, 0)
            for degree, coefficient in left.items()
        )
    return sp.together(answer)


def integer_terms(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> dict[tuple[int, ...], int]:
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    denominator = 1
    for coefficient in polynomial.coeffs():
        denominator = sp.ilcm(denominator, coefficient.q)
    return {
        exponents: int(coefficient * denominator)
        for exponents, coefficient in polynomial.terms()
    }


def interval_payload(interval: Interval) -> list[str]:
    return [str(interval.lo), str(interval.hi)]


def reversal(item: tuple[int, int, tuple[tuple[int, int], ...]]):
    r, s, support = item
    return (4 - s, 4 - r, tuple(sorted((4 - row, 1 - column) for row, column in support)))


def double_shear_census() -> dict[str, object]:
    charts = []
    for r, s in combinations(range(5), 2):
        positions = [
            (row, column)
            for row in range(5)
            if row not in {r, s}
            for column in range(2)
        ]
        for support in combinations(positions, 2):
            charts.append((r, s, tuple(support)))
    representatives = {min(chart, reversal(chart)) for chart in charts}
    profiles = {
        "same_row_different_columns": 0,
        "different_rows_same_column": 0,
        "different_rows_different_columns": 0,
    }
    for _, _, support in representatives:
        same_row = support[0][0] == support[1][0]
        same_column = support[0][1] == support[1][1]
        if same_row:
            profiles["same_row_different_columns"] += 1
        elif same_column:
            profiles["different_rows_same_column"] += 1
        else:
            profiles["different_rows_different_columns"] += 1
    assert len(charts) == 150
    assert len(representatives) == 78
    assert profiles == {
        "same_row_different_columns": 16,
        "different_rows_same_column": 30,
        "different_rows_different_columns": 32,
    }
    return {
        "labeled_chart_count": len(charts),
        "reversal_orbit_count": len(representatives),
        "orbit_profiles": profiles,
    }


def main() -> None:
    a, b01, b02, b03, b11, b12, b13 = sp.symbols(
        "a b01 b02 b03 b11 b12 b13"
    )
    variables = (a, b01, b02, b03, b11, b12, b13)
    eliminated_shear = -a * b01 / b13
    left_0 = {0: sp.Integer(1), 1: a}
    left_1 = {4: sp.Integer(1), 3: eliminated_shear}
    right_0 = {1: b01, 2: b02, 3: b03, 4: sp.Integer(1)}
    right_1 = {0: sp.Integer(1), 1: b11, 2: b12, 3: b13}

    rational_moments = {
        order: raw_moment(left_0, left_1, right_0, right_1, order)
        for order in range(1, 10)
    }
    assert rational_moments[1] == 0
    moment_terms = {}
    moment_denominators = {}
    for order in range(2, 10):
        numerator, denominator = sp.fraction(rational_moments[order])
        moment_terms[order] = integer_terms(numerator, variables)
        moment_denominators[order] = integer_terms(denominator, variables)

    system = [moment_terms[order] for order in range(2, 9)]
    center = [Fraction(value) for value in CENTER_STRINGS]
    box = [Interval(value - RADIUS, value + RADIUS) for value in center]
    assert all(not (interval.lo <= 0 <= interval.hi) for interval in box)
    eliminated_interval = -(box[0] * box[1] / box[6])
    assert not (eliminated_interval.lo <= 0 <= eliminated_interval.hi)

    jacobian = [
        [derivative(polynomial, column) for column in range(7)]
        for polynomial in system
    ]
    jacobian_center = sp.Matrix(
        [
            [
                sp.Rational(
                    evaluate(entry, center).numerator,
                    evaluate(entry, center).denominator,
                )
                for entry in row
            ]
            for row in jacobian
        ]
    )
    assert jacobian_center.det() != 0
    inverse = [
        [Fraction(int(value.p), int(value.q)) for value in row]
        for row in jacobian_center.inv().tolist()
    ]
    residual = [evaluate(polynomial, center) for polynomial in system]
    newton_center = [
        center[i] - sum(inverse[i][j] * residual[j] for j in range(7))
        for i in range(7)
    ]
    jacobian_box = [
        [evaluate(entry, box) for entry in row]
        for row in jacobian
    ]
    defect = []
    for i in range(7):
        row = []
        for j in range(7):
            value = as_interval(1 if i == j else 0)
            value -= sum(
                (inverse[i][k] * jacobian_box[k][j] for k in range(7)),
                start=as_interval(0),
            )
            row.append(value)
        defect.append(row)
    delta = Interval(-RADIUS, RADIUS)
    krawczyk = [
        as_interval(newton_center[i])
        + sum(
            (defect[i][j] * delta for j in range(7)),
            start=as_interval(0),
        )
        for i in range(7)
    ]
    assert all(
        image.strict_subset(domain)
        for image, domain in zip(krawczyk, box, strict=True)
    )

    ninth_numerator = evaluate(moment_terms[9], box)
    ninth_denominator = evaluate(moment_denominators[9], box)
    ninth = ninth_numerator / ninth_denominator / factorial(37)
    assert ninth.lo > 0

    census = double_shear_census()
    artifact = {
        "format": "two-pair-sic-bidegree44-rank-two-double-shear-real-prefix-v1",
        "field": "real characteristic zero",
        "direct_factor_chart": {
            "pivot_rows": [0, 4],
            "U_columns": ["e0+a*e1", "e4+d*e3"],
            "B_normalization": ["b04=1", "b10=1"],
            "mu1_elimination": "d=-a*b01/b13",
            "pivot_minor_columns_0_4": "-1",
            "variables": [str(variable) for variable in variables],
            "coefficient_torus": True,
        },
        "double_shear_census": census,
        "moment_orders_vanishing": [1, 8],
        "moment_term_counts": {
            str(order): len(moment_terms[order]) for order in range(2, 10)
        },
        "center": list(CENTER_STRINGS),
        "radius": str(RADIUS),
        "box": [interval_payload(interval) for interval in box],
        "krawczyk_image": [interval_payload(interval) for interval in krawczyk],
        "eliminated_shear_interval": interval_payload(eliminated_interval),
        "jacobian_center_det_sha256": sha256(
            str(jacobian_center.det()).encode()
        ).hexdigest(),
        "normalized_mu9_interval": interval_payload(ninth),
        "conclusion": (
            "the box contains a unique real exact-rank-two coefficient-torus "
            "point with mu_1 through mu_8 zero, and mu_9 is positive there"
        ),
        "relative_period": "mu_m/(4*m+1)!=CT_u integral_0^1 P(u,t)^m dt",
        "recurrence_stage": (
            "not entered: the exact ninth moment already obstructs this point"
        ),
        "scope": (
            "one different-row/different-column double-shear orbit and one "
            "isolated real component; global complex orbit exclusion and the "
            "other 77 reversal orbits remain open"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS 150 double-shear charts reduce to 78 reversal orbits")
    print("PASS exact rational Krawczyk box isolates a real mu_1,...,mu_8 zero")
    print("PASS coefficient rank is exactly two and every torus coordinate is nonzero")
    print("PASS normalized mu_9 is positive throughout the isolating box")


if __name__ == "__main__":
    main()
