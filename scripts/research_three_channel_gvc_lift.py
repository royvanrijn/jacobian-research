#!/usr/bin/env python3
"""Exact audits for the weighted three-channel GVC(3) lift.

The positive grading 2*z+t+y=4 contains ordinary degrees 2, 3, and 4.
It therefore keeps every surviving pure contraction scalar while allowing
the channel-count convolution that is absent from the two-channel tag.

This script checks five exact characteristic-zero statements.

1. All rank-three four-point parallelograms on the weighted quartic plane
   fail the pure premise by moment four.
2. The exceptional parallelogram plus the middle point z*t*y has pure
   radical d*e=0 and h^2+4=0.
3. In isotropic coordinates, the complete quartic polynomial jet for
   Lambda=d_z*d_L*d_M has radical a2=a3=a4=0 through moment eight.
4. A sparse activation of the quadratic and quartic operator endpoints has
   a terminal radical through moment ten.
5. The complete odd-quartic operator/polynomial jet has radical
   (A, S, R*U) through moment six.

The accompanying note gives the all-order support cutoffs on the surviving
radicals.  No GVC(3) counterexample is claimed.
"""

from __future__ import annotations

import json
from itertools import combinations
from math import factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "three_channel_gvc_lift.json"
)
Exponent = tuple[int, int, int]
SparsePolynomial = dict[Exponent, object]


def exponent_add(left: Exponent, right: Exponent) -> Exponent:
    return tuple(left[index] + right[index] for index in range(3))


def factorial_weight(exponent: Exponent) -> int:
    result = 1
    for value in exponent:
        result *= factorial(value)
    return result


def multiply(
    left: SparsePolynomial,
    right: SparsePolynomial,
) -> SparsePolynomial:
    result: SparsePolynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = exponent_add(left_exponent, right_exponent)
            result[exponent] = sp.expand(
                result.get(exponent, 0)
                + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient != 0
    }


def moment_prefix(
    operator: SparsePolynomial,
    polynomial: SparsePolynomial,
    cutoff: int,
) -> list[object]:
    operator_power: SparsePolynomial = {(0, 0, 0): 1}
    polynomial_power: SparsePolynomial = {(0, 0, 0): 1}
    moments: list[object] = []
    for _order in range(1, cutoff + 1):
        operator_power = multiply(operator_power, operator)
        polynomial_power = multiply(polynomial_power, polynomial)
        moment = sum(
            operator_coefficient
            * polynomial_power.get(exponent, 0)
            * factorial_weight(exponent)
            for exponent, operator_coefficient in operator_power.items()
        )
        moments.append(sp.factor(moment))
    return moments


def is_unit_basis(basis: sp.GroebnerBasis) -> bool:
    return (
        len(basis.polys) == 1
        and basis.polys[0].as_expr() == 1
    )


def weighted_points() -> list[Exponent]:
    return [
        (z_exponent, t_exponent, 4 - 2 * z_exponent - t_exponent)
        for z_exponent in range(3)
        for t_exponent in range(5 - 2 * z_exponent)
    ]


def rank_three_parallelograms():
    points = weighted_points()
    seen = set()
    answer = []
    for four_points in combinations(points, 4):
        for endpoints in combinations(four_points, 2):
            middle = [
                point
                for point in four_points
                if point not in endpoints
            ]
            if exponent_add(*endpoints) != exponent_add(*middle):
                continue
            matrix = sp.Matrix.hstack(
                *[
                    sp.Matrix(point)
                    for point in (endpoints[0], middle[0], middle[1])
                ]
            )
            if matrix.det() == 0:
                continue
            key = (tuple(sorted(endpoints)), tuple(sorted(middle)))
            if key in seen:
                continue
            seen.add(key)
            answer.append((endpoints, middle))
    return answer


def parallelogram_audit() -> dict[str, object]:
    d, e = sp.symbols("d e")
    histogram: dict[int, int] = {}
    delayed = []
    parallelograms = rank_three_parallelograms()
    for endpoints, middle in parallelograms:
        for q1, q4 in (endpoints, endpoints[::-1]):
            terminal_coefficient = -sp.Rational(
                factorial_weight(q1),
                factorial_weight(q4),
            )
            operator = {q1: 1, q4: 1}
            polynomial = {
                q1: 1,
                middle[0]: d,
                middle[1]: e,
                q4: terminal_coefficient,
            }
            moments = moment_prefix(operator, polynomial, 4)
            first_unit = None
            bases = []
            for order in range(1, 5):
                basis = sp.groebner(
                    moments[:order],
                    d,
                    e,
                    order="grevlex",
                )
                bases.append(basis)
                if is_unit_basis(basis):
                    first_unit = order
                    break
            assert first_unit is not None
            histogram[first_unit] = histogram.get(first_unit, 0) + 1
            if first_unit == 4:
                delayed.append(
                    {
                        "operator_endpoints": [list(q1), list(q4)],
                        "polynomial_middle": [list(x) for x in middle],
                        "moments": [str(moment) for moment in moments],
                        "moment_3_basis": [
                            str(sp.factor(item.as_expr()))
                            for item in bases[2].polys
                        ],
                    }
                )
    assert len(parallelograms) == 28
    assert histogram == {3: 54, 4: 2}
    return {
        "weighted_points": [list(point) for point in weighted_points()],
        "unoriented_parallelograms": len(parallelograms),
        "oriented_normalizations": 56,
        "first_unit_moment_histogram": {
            str(key): value
            for key, value in sorted(histogram.items())
        },
        "delayed_orientations": delayed,
    }


def five_term_persistent_branch() -> dict[str, object]:
    d, e, h = sp.symbols("d e h")
    operator = {(1, 0, 2): 1, (1, 2, 0): 1}
    polynomial = {
        (1, 0, 2): 1,
        (1, 2, 0): -1,
        (0, 2, 2): d,
        (2, 0, 0): e,
        (1, 1, 1): h,
    }
    moments = moment_prefix(operator, polynomial, 10)
    basis = sp.groebner(moments, d, e, h, order="grevlex")
    basis_expressions = [
        sp.factor(item.as_expr())
        for item in basis.polys
    ]
    assert basis_expressions == [d * e, h**2 + 4]
    return {
        "operator": "z*(t^2+y^2)",
        "polynomial": (
            "z*(y^2-t^2+h*t*y)+d*t^2*y^2+e*z^2"
        ),
        "orders": list(range(1, 11)),
        "groebner_basis": [str(item) for item in basis_expressions],
        "all_order_terminal_reason": (
            "on d*e=0 and h^2=-4 the cubic core is z times an "
            "isotropic square; a fixed multiplier bridges only bounded "
            "z-tag depth, after which the required transverse derivative "
            "degree exceeds its fixed supply"
        ),
    }


def complete_polynomial_quartic() -> dict[str, object]:
    coefficients = sp.symbols("a0:5")
    operator = {(1, 1, 1): 1}
    polynomial = {
        (1, 2, 0): 1,
        (2, 0, 0): 1,
        **{
            (0, 4 - index, index): coefficient
            for index, coefficient in enumerate(coefficients)
        },
    }
    moments = moment_prefix(operator, polynomial, 8)
    basis_7 = sp.groebner(
        moments[:7],
        *coefficients,
        order="grevlex",
    )
    basis_8 = sp.groebner(
        moments,
        *coefficients,
        order="grevlex",
    )
    expected_7 = [
        coefficients[4] * (coefficients[0] + 1),
        coefficients[1] * coefficients[4],
        coefficients[2],
        coefficients[3],
    ]
    expected_8 = [
        coefficients[4] * (coefficients[0] + 1),
        coefficients[1] * coefficients[4],
        coefficients[4] ** 2,
        coefficients[2],
        coefficients[3],
    ]
    assert [
        sp.factor(item.as_expr()) for item in basis_7.polys
    ] == expected_7
    assert [
        sp.factor(item.as_expr()) for item in basis_8.polys
    ] == expected_8
    return {
        "operator": "d_z*d_L*d_M",
        "polynomial": (
            "z*L^2+z^2+a0*L^4+a1*L^3*M+a2*L^2*M^2"
            "+a3*L*M^3+a4*M^4"
        ),
        "moment_7_basis": [str(item) for item in expected_7],
        "moment_8_basis": [str(item) for item in expected_8],
        "radical": "(a2,a3,a4)",
        "status": "exact characteristic-zero terminal radical",
    }


def activated_endpoint_model() -> dict[str, object]:
    a2, b, c, f, g = sp.symbols("A2 B C F G")
    variables = (a2, b, c, f, g)
    operator = {
        (1, 1, 1): 1,
        (2, 0, 0): a2,
        (0, 4, 0): b,
        (0, 0, 4): c,
    }
    polynomial = {
        (1, 2, 0): 1,
        (2, 0, 0): 1,
        (0, 4, 0): f,
        (0, 0, 4): g,
    }
    moments = moment_prefix(operator, polynomial, 10)
    basis = sp.groebner(moments, *variables, order="grevlex")
    actual = [
        sp.factor(item.as_expr())
        for item in basis.polys
    ]
    expected = [
        a2 * b**2,
        a2 * (3 * a2 + b),
        a2 + 12 * b * f + 12 * c * g,
        a2 * g,
        b * g,
        g * (f + 1),
        g**2,
    ]
    assert actual == expected
    return {
        "operator": "z*L*M+A2*z^2+B*L^4+C*M^4",
        "polynomial": "z*L^2+z^2+F*L^4+G*M^4",
        "orders": list(range(1, 11)),
        "groebner_basis": [str(item) for item in actual],
        "radical": "(A2,G,B*F)",
        "status": "exact characteristic-zero terminal radical",
    }


def complete_odd_quartic_jet() -> dict[str, object]:
    a2, u, v, r, s = sp.symbols("A2 U V R S")
    variables = (a2, u, v, r, s)
    operator = {
        (1, 1, 1): 1,
        (2, 0, 0): a2,
        (0, 3, 1): u,
        (0, 1, 3): v,
    }
    polynomial = {
        (1, 2, 0): 1,
        (2, 0, 0): 1,
        (0, 3, 1): r,
        (0, 1, 3): s,
    }
    moments = moment_prefix(operator, polynomial, 6)
    basis = sp.groebner(moments, *variables, order="grevlex")
    actual = [
        sp.factor(item.as_expr())
        for item in basis.polys
    ]
    expected = [
        a2**2,
        a2 + 3 * r * u + 3 * s * v,
        a2 * s,
        s * (10 * u + 1),
        r * s,
        s**2,
    ]
    assert actual == expected
    return {
        "operator": "d_z*d_L*d_M+A2*d_z^2+U*d_L^3*d_M+V*d_L*d_M^3",
        "polynomial": "z*L^2+z^2+R*L^3*M+S*L*M^3",
        "orders": list(range(1, 7)),
        "groebner_basis": [str(item) for item in actual],
        "radical": "(A2,S,R*U)",
        "components": [
            "(A2,S,R)",
            "(A2,S,U)",
        ],
        "all_order_terminal_reason": (
            "on R=0 every operator factor has positive M-degree while "
            "P is M-free; on U=0, choosing k copies of d_L*d_M^3 "
            "forces k at most deg_M(Q)/2 and then z-degree balance "
            "forces m at most 2*deg_M(Q)+deg_z(Q)"
        ),
        "status": "exact characteristic-zero terminal radical",
    }


def main() -> None:
    artifact = {
        "format": "three-channel-gvc-lift-v2",
        "field": "Q",
        "grading": "2*z+t+y=4",
        "parallelogram_audit": parallelogram_audit(),
        "five_term_persistent_branch": five_term_persistent_branch(),
        "complete_polynomial_quartic": complete_polynomial_quartic(),
        "activated_endpoint_model": activated_endpoint_model(),
        "complete_odd_quartic_jet": complete_odd_quartic_jet(),
        "conclusion": {
            "gvc3_counterexample_found": False,
            "closed": (
                "all weighted-quartic rank-three parallelograms and three "
                "minimal endpoint repairs, including the complete "
                "odd-quartic operator/polynomial jet"
            ),
            "remaining": (
                "the simultaneous complete even-and-odd quartic jet, a "
                "different positive grading, or a different purification"
            ),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS 56 oriented weighted-quartic parallelogram charts")
    print("PASS exact five-term persistent radical and mixed cutoff")
    print("PASS complete polynomial-quartic moment-eight terminal radical")
    print("PASS activated endpoint moment-ten terminal radical")
    print("PASS complete odd-quartic moment-six terminal radical")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
