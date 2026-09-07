#!/usr/bin/env python3
"""Derive the algebraic identities behind Kihara's fifteen sections.

The computation has three deliberately separate layers:

* a universal degree-five obstruction for six paired centers;
* Kihara's two-parameter center relation and generic ``P13`` square;
* the two further squares ``P14,P15`` after Kihara's one-parameter base
  change.

Everything is checked in a rational polynomial or rational-function field.
This proves point-on-quartic identities only; it is not a new rank claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import reduce
from math import gcd

import sympy as sp


def square_approximant(product: sp.Poly, variable: sp.Symbol) -> sp.Poly:
    """Return the monic sextic matching a monic degree-12 square through degree 6."""

    if product.degree() != 12 or product.LC() != 1:
        raise ValueError("a monic degree-twelve polynomial is required")
    coefficients = [sp.Integer(0)] * 7
    coefficients[6] = sp.Integer(1)
    for index in range(5, -1, -1):
        approximant = sp.Poly(
            sum(coefficients[degree] * variable**degree for degree in range(7)),
            variable,
        )
        degree = 6 + index
        coefficients[index] = sp.cancel(
            (product.nth(degree) - (approximant * approximant).nth(degree)) / 2
        )
    return sp.Poly(
        sum(coefficients[degree] * variable**degree for degree in range(7)),
        variable,
    )


def paired_remainder(
    centers: tuple[sp.Expr, ...], shift: sp.Expr, variable: sp.Symbol
) -> sp.Poly:
    product = sp.Poly(1, variable)
    for center in centers:
        product *= sp.Poly(variable - center - shift, variable)
        product *= sp.Poly(variable - center + shift, variable)
    approximant = square_approximant(product, variable)
    return approximant * approximant - product


def kihara_centers(p: sp.Symbol, q: sp.Symbol) -> tuple[sp.Expr, ...]:
    return (
        sp.Integer(0),
        (2 * p**2 + p * q + 2 * q**2) ** 2,
        2 * (p + q) ** 2 * (2 * p**2 + p * q + q**2),
        q**2 * (4 * p**2 - p * q + 4 * q**2),
        p * (2 * p - q) * (2 * p**2 + 4 * p * q + 5 * q**2),
        4 * p**4
        + 8 * p**3 * q
        + 9 * p**2 * q**2
        - 2 * p * q**3
        + 2 * q**4,
    )


def center_relation() -> dict[str, object]:
    x, u = sp.symbols("x u")
    c0, c1, c2, c3, c4, c5 = sp.symbols("c0 c1 c2 c3 c4 c5")
    center_polynomial = (
        x**6 + c5 * x**5 + c4 * x**4 + c3 * x**3 + c2 * x**2 + c1 * x + c0
    )
    product = sp.Poly(
        sp.expand(
            center_polynomial.subs(x, x - u)
            * center_polynomial.subs(x, x + u)
        ),
        x,
    )
    remainder = square_approximant(product, x) ** 2 - product
    obstruction = sp.factor(remainder.nth(5))
    expected = -u**2 * (
        24 * c1
        - 8 * c2 * c5
        - 12 * c3 * c4
        + 7 * c3 * c5**2
        + 8 * c4**2 * c5
        - 6 * c4 * c5**3
        + c5**5
    )
    if sp.expand(obstruction - expected) != 0:
        raise AssertionError("universal degree-five obstruction changed")
    centered = sp.factor(obstruction.subs(c5, 0))
    if centered != -12 * u**2 * (2 * c1 - c3 * c4):
        raise AssertionError("centered obstruction changed")

    p, q, z = sp.symbols("p q z")
    centers = kihara_centers(p, q)
    mean = sp.cancel(sum(centers) / 6)
    polynomial = sp.Poly(sp.prod(z - (center - mean) for center in centers), z)
    alpha = (
        16 * p**8
        + 56 * p**7 * q
        + 116 * p**6 * q**2
        + 110 * p**5 * q**3
        + 86 * p**4 * q**4
        + 41 * p**3 * q**5
        + 82 * p**2 * q**6
        + 24 * p * q**7
        + 12 * q**8
    )
    beta = (
        p
        * (p - q)
        * (p + 2 * q)
        * (2 * p - q)
        * (2 * p + 5 * q)
        * (4 * p + q)
        * (p**2 + 2 * q**2)
        * (2 * p**2 + 2 * p * q + 3 * q**2)
        * (4 * p**2 + 4 * p * q + 3 * q**2)
    )
    expected_coefficients = {
        4: -2 * alpha / 3,
        3: 2 * beta / 27,
        1: -2 * alpha * beta / 81,
    }
    for degree, value in expected_coefficients.items():
        if sp.cancel(polynomial.nth(degree) - value) != 0:
            raise AssertionError(f"Kihara centered coefficient z^{degree} changed")
    if sp.cancel(2 * polynomial.nth(1) - polynomial.nth(3) * polynomial.nth(4)) != 0:
        raise AssertionError("Kihara's intrinsic center relation failed")
    return {
        "universal_remainder_x5_coefficient": str(obstruction),
        "mean_zero_remainder_x5_coefficient": str(centered),
        "mean_zero_elementary_relation": "2*e5=e2*e3",
        "kihara_centered_c4": str(sp.factor(polynomial.nth(4))),
        "kihara_centered_c3": str(sp.factor(polynomial.nth(3))),
        "kihara_centered_c1": str(sp.factor(polynomial.nth(1))),
        "kihara_relation_exact": True,
    }


def generic_p13_identity() -> tuple[sp.Poly, dict[str, object]]:
    x, p, q, u = sp.symbols("x p q u")
    remainder = paired_remainder(kihara_centers(p, q), u, x)
    if remainder.degree() != 4:
        raise AssertionError("Kihara's generic remainder is not quartic")
    denominator = 2 * p**2 + 2 * p * q + 3 * q**2
    linear_u = 2 * p**2 + 4 * p * q + 5 * q**2
    degree_six = (
        8 * p**6
        + 28 * p**5 * q
        + 58 * p**4 * q**2
        + 69 * p**3 * q**3
        + 76 * p**2 * q**4
        + 40 * p * q**5
        + 22 * q**6
    )
    x13 = sp.cancel((linear_u * u + degree_six) / denominator)
    h0 = (
        q
        * (p + q)
        * (4 * p**2 + 5 * q**2)
        * (2 * p**2 + p * q + 8 * q**2)
        * denominator**2
        * (2 * p**2 + 3 * p * q + 2 * q**2)
        * linear_u
        * (3 * p**2 + 2 * p * q + 2 * q**2)
    )
    h1 = 2 * (
        64 * p**12
        + 576 * p**11 * q
        + 2656 * p**10 * q**2
        + 8192 * p**9 * q**3
        + 19120 * p**8 * q**4
        + 35264 * p**7 * q**5
        + 52339 * p**6 * q**6
        + 62738 * p**5 * q**7
        + 60769 * p**4 * q**8
        + 45996 * p**3 * q**9
        + 26455 * p**2 * q**10
        + 10260 * p * q**11
        + 2352 * q**12
    )
    h2 = (
        4
        * linear_u
        * (
            8 * p**6
            + 40 * p**5 * q
            + 84 * p**4 * q**2
            + 125 * p**3 * q**3
            + 157 * p**2 * q**4
            + 105 * p * q**5
            + 48 * q**6
        )
    )
    h3 = 24 * q * (p + q) * (2 * p**2 + 3 * p * q + 4 * q**2)
    h = h0 + h1 * u + h2 * u**2 + h3 * u**3
    y13 = sp.cancel(
        2
        * p
        * u
        * q**2
        * (2 * p - q)
        * (p + q) ** 2
        * (2 * p**2 + p * q + 2 * q**2)
        * h
        / denominator**2
    )
    if sp.cancel(remainder.eval(x13) - y13**2) != 0:
        raise AssertionError("the generic P13 square identity failed")
    return remainder, {
        "field": "QQ(p,q,u)",
        "x13": str(x13),
        "y13_factorized": str(sp.factor(y13)),
        "h_coefficients_low_to_high_in_u": [
            str(sp.factor(value)) for value in (h0, h1, h2, h3)
        ],
        "identity_exact": True,
    }


def _rational_square_root(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    if value < 0:
        raise ArithmeticError("negative rational constant")
    numerator, numerator_exact = sp.integer_nthroot(int(value.p), 2)
    denominator, denominator_exact = sp.integer_nthroot(int(value.q), 2)
    if not numerator_exact or not denominator_exact:
        raise ArithmeticError("rational constant is not a square")
    return sp.Rational(numerator, denominator)


def _primitive_coefficients(polynomial: sp.Poly) -> tuple[int, ...]:
    _, integral = polynomial.clear_denoms(convert=True)
    coefficients = [int(integral.nth(index)) for index in range(integral.degree() + 1)]
    content = reduce(gcd, (abs(value) for value in coefficients if value))
    coefficients = [value // content for value in coefficients]
    if coefficients[-1] < 0:
        coefficients = [-value for value in coefficients]
    return tuple(coefficients)


def _coefficient_digest(polynomial: sp.Poly) -> str:
    coefficients = _primitive_coefficients(polynomial)
    return hashlib.sha256(
        ("\n".join(map(str, coefficients)) + "\n").encode("ascii")
    ).hexdigest()


def factored_square_root(
    value: sp.Expr, parameter: sp.Symbol, *, include_polynomials: bool
) -> tuple[sp.Expr, dict[str, object]]:
    numerator, denominator = sp.fraction(sp.cancel(value))
    numerator_constant, numerator_factors = sp.factor_list(numerator, parameter)
    denominator_constant, denominator_factors = sp.factor_list(denominator, parameter)
    constant = sp.Rational(numerator_constant) / sp.Rational(denominator_constant)
    root = _rational_square_root(constant)
    for factor, exponent in numerator_factors:
        if exponent % 2:
            raise ArithmeticError("odd numerator factor exponent")
        root *= factor ** (exponent // 2)
    for factor, exponent in denominator_factors:
        if exponent % 2:
            raise ArithmeticError("odd denominator factor exponent")
        root /= factor ** (exponent // 2)
    root = sp.cancel(root)
    if sp.cancel(root**2 - value) != 0:
        raise AssertionError("reconstructed rational-function square root failed")

    def records(factors: list[tuple[sp.Expr, int]]) -> list[dict[str, object]]:
        answer = []
        for factor, exponent in factors:
            polynomial = sp.Poly(factor, parameter, domain=sp.QQ)
            coefficients = _primitive_coefficients(polynomial)
            record: dict[str, object] = {
                "degree": int(polynomial.degree()),
                "exponent_in_square": int(exponent),
                "exponent_in_root": int(exponent // 2),
                "primitive_coefficients_sha256": _coefficient_digest(polynomial),
                "first_three_ascending_coefficients": list(coefficients[:3]),
                "last_three_ascending_coefficients": list(coefficients[-3:]),
            }
            if include_polynomials or polynomial.degree() <= 8:
                record["factor"] = str(factor)
            answer.append(record)
        return answer

    return root, {
        "square_constant": str(constant),
        "square_root_constant": str(_rational_square_root(constant)),
        "numerator_factors": records(numerator_factors),
        "denominator_factors": records(denominator_factors),
    }


def specialized_p14_p15_identities(*, include_polynomials: bool) -> dict[str, object]:
    x, t = sp.symbols("x t")
    p = t**2 * (8 + 3 * t**2)
    q = -6 * (2 + t**2) * (4 + t**2)
    first = 2304 + 2400 * t**2 + 928 * t**4 + 150 * t**6 + 9 * t**8
    common = 1152 + 1632 * t**2 + 860 * t**4 + 201 * t**6 + 18 * t**8
    u = 4 * (2 + t**2) * first * common / t
    remainder = paired_remainder(kihara_centers(p, q), u, x)
    if remainder.degree() != 4:
        raise AssertionError("the specialized Kihara remainder is not quartic")
    polynomial14 = (
        10616832
        - 18579456 * t
        + 33619968 * t**2
        - 51535872 * t**3
        + 45895680 * t**4
        - 61848576 * t**5
        + 35397888 * t**6
        - 41945856 * t**7
        + 16968640 * t**8
        - 17591104 * t**9
        + 5232272 * t**10
        - 4675248 * t**11
        + 1035180 * t**12
        - 769824 * t**13
        + 126252 * t**14
        - 71874 * t**15
        + 8559 * t**16
        - 2916 * t**17
        + 243 * t**18
    )
    denominator14 = 2304 + 3168 * t**2 + 1580 * t**4 + 339 * t**6 + 27 * t**8
    x14 = sp.cancel(-4 * common * polynomial14 / (t * denominator14))
    left15 = -48 + 24 * t - 34 * t**2 + 16 * t**3 - 6 * t**4 + 3 * t**5
    right15 = 96 + 80 * t**2 + 4 * t**3 + 18 * t**4 + 3 * t**5
    x15 = sp.cancel(4 * left15 * right15 * common / t)
    answer: dict[str, object] = {}
    expected_high_factors = {
        "P14": (58, "b6b5bb73a584cacddcbe1e45d0fb7e839f1489ffaffa3f2db14e606d0784326a"),
        "P15": (44, "bedb18dc27516b36336fc7690319a05e2b32b1052ca4068c5cd0412860fcd585"),
    }
    for name, abscissa in (("P14", x14), ("P15", x15)):
        value = sp.cancel(remainder.eval(abscissa))
        ordinate, factorization = factored_square_root(
            value, t, include_polynomials=include_polynomials
        )
        high_factors = [
            record
            for record in factorization["numerator_factors"]
            if int(record["degree"]) > 8
        ]
        if len(high_factors) != 1:
            raise AssertionError(f"{name} did not have one high-degree factor")
        expected_degree, expected_digest = expected_high_factors[name]
        if (
            high_factors[0]["degree"] != expected_degree
            or high_factors[0]["primitive_coefficients_sha256"] != expected_digest
        ):
            raise AssertionError(f"{name} high-degree factor changed")
        answer[name] = {
            "field": "QQ(t)",
            "abscissa": str(abscissa),
            "ordinate": str(ordinate) if include_polynomials else None,
            "factorization": factorization,
            "identity_exact": True,
        }
    return answer


def derive(*, include_polynomials: bool = False) -> dict[str, object]:
    center = center_relation()
    _, p13 = generic_p13_identity()
    p14_p15 = specialized_p14_p15_identities(
        include_polynomials=include_polynomials
    )
    return {
        "schema": "scratch.kihara-rank14-algebraic-identities.v1",
        "claim_level": "exact symbolic point identities; no new rank claim",
        "center_geometry": center,
        "extra_sections": {"P13": p13, **p14_p15},
        "all_identities_exact": True,
        "interpretation": (
            "the centered relation is necessary and sufficient for the paired "
            "degree-twelve square remainder to have degree at most four when u is "
            "nonzero; the extra-section formulas are construction-specific sufficient "
            "identities"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-polynomials",
        action="store_true",
        help="include the large P14/P15 square-root factors and ordinates",
    )
    arguments = parser.parse_args()
    print(
        json.dumps(
            derive(include_polynomials=arguments.full_polynomials),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
