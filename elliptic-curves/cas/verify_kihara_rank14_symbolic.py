#!/usr/bin/env python3
"""Symbolically verify the generic point identities in Kihara's 2001 family.

This is intentionally a narrow verifier.  It checks the printed construction
in ``QQ(t)[x]`` and proves that the three additional printed abscissae give
squares in ``QQ(t)``.  It does not recompute the paper's numerical height
pairing or independently certify generic rank 14.
"""

from __future__ import annotations

import argparse
import json

import sympy as sym


def _constant_is_rational_square(value: sym.Expr) -> bool:
    value = sym.Rational(value)
    if value < 0:
        return False
    numerator_root, numerator_exact = sym.integer_nthroot(int(value.p), 2)
    denominator_root, denominator_exact = sym.integer_nthroot(int(value.q), 2)
    return bool(
        numerator_exact
        and denominator_exact
        and numerator_root**2 == value.p
        and denominator_root**2 == value.q
    )


def _square_signature(
    value: sym.Expr, parameter_t: sym.Symbol
) -> dict[str, object]:
    numerator, denominator = sym.fraction(sym.cancel(value))
    numerator_constant, numerator_factors = sym.factor_list(numerator, parameter_t)
    denominator_constant, denominator_factors = sym.factor_list(
        denominator, parameter_t
    )
    factors = tuple(numerator_factors) + tuple(denominator_factors)
    constant = sym.Rational(numerator_constant) / sym.Rational(denominator_constant)
    return {
        "is_square": _constant_is_rational_square(constant)
        and all(exponent % 2 == 0 for _, exponent in factors),
        "numerator_degree": int(sym.degree(numerator, parameter_t)),
        "denominator_degree": int(sym.degree(denominator, parameter_t)),
        "numerator_factor_degree_exponents": [
            [int(sym.degree(factor, parameter_t)), int(exponent)]
            for factor, exponent in numerator_factors
        ],
        "denominator_factor_degree_exponents": [
            [int(sym.degree(factor, parameter_t)), int(exponent)]
            for factor, exponent in denominator_factors
        ],
    }


def symbolic_verification() -> dict[str, object]:
    x, t = sym.symbols("x t")
    function_field = sym.QQ.frac_field(t)

    p = t**2 * (8 + 3 * t**2)
    q = -6 * (2 + t**2) * (4 + t**2)
    u = (
        4
        * (2 + t**2)
        * (2304 + 2400 * t**2 + 928 * t**4 + 150 * t**6 + 9 * t**8)
        * (1152 + 1632 * t**2 + 860 * t**4 + 201 * t**6 + 18 * t**8)
        / t
    )
    a_values = (
        0,
        (2 * p**2 + p * q + 2 * q**2) ** 2,
        2 * (p + q) ** 2 * (2 * p**2 + p * q + q**2),
        q**2 * (4 * p**2 - p * q + 4 * q**2),
        p * (2 * p - q) * (2 * p**2 + 4 * p * q + 5 * q**2),
        4 * p**4 + 8 * p**3 * q + 9 * p**2 * q**2 - 2 * p * q**3 + 2 * q**4,
    )
    b_values = tuple(u + value for value in a_values) + tuple(
        -u + value for value in a_values
    )

    product = sym.Poly(1, x, domain=function_field)
    for root in b_values:
        product *= sym.Poly(x - root, x, domain=function_field)

    approximant_coefficients = [sym.S(0)] * 7
    approximant_coefficients[6] = sym.S(1)
    for index in range(5, -1, -1):
        approximant = sym.Poly(
            sum(
                approximant_coefficients[degree] * x**degree
                for degree in range(7)
            ),
            x,
            domain=function_field,
        )
        square = approximant * approximant
        approximant_coefficients[index] = sym.cancel(
            (product.nth(6 + index) - square.nth(6 + index)) / 2
        )
    approximant = sym.Poly(
        sum(
            approximant_coefficients[degree] * x**degree for degree in range(7)
        ),
        x,
        domain=function_field,
    )
    remainder = approximant * approximant - product
    if remainder.degree() != 4:
        raise AssertionError("Kihara's symbolic remainder is not quartic")

    visible_identities = [
        sym.cancel(remainder.eval(root) - approximant.eval(root) ** 2) == 0
        for root in b_values
    ]

    denominator_13 = 2 * p**2 + 2 * p * q + 3 * q**2
    numerator_13 = 2 * p**2 + 4 * p * q + 5 * q**2
    degree_six = (
        8 * p**6
        + 28 * p**5 * q
        + 58 * p**4 * q**2
        + 69 * p**3 * q**3
        + 76 * p**2 * q**4
        + 40 * p * q**5
        + 22 * q**6
    )
    x_13 = sym.cancel(
        numerator_13 * u / denominator_13 + degree_six / denominator_13
    )

    common = 1152 + 1632 * t**2 + 860 * t**4 + 201 * t**6 + 18 * t**8
    polynomial_14 = (
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
    x_14 = sym.cancel(
        -4
        * common
        * polynomial_14
        / (t * (2304 + 3168 * t**2 + 1580 * t**4 + 339 * t**6 + 27 * t**8))
    )
    x_15 = sym.cancel(
        4
        * (-48 + 24 * t - 34 * t**2 + 16 * t**3 - 6 * t**4 + 3 * t**5)
        * (96 + 80 * t**2 + 4 * t**3 + 18 * t**4 + 3 * t**5)
        * common
        / t
    )
    extra_signatures = {
        name: _square_signature(remainder.eval(abscissa), t)
        for name, abscissa in (("P13", x_13), ("P14", x_14), ("P15", x_15))
    }

    # A single exact smooth specialization proves that the generic binary
    # quartic discriminant is not the zero rational function.  Forming the
    # full expanded discriminant in QQ(t) is unnecessary and much larger.
    e, d, c, b, a = (
        sym.cancel(remainder.nth(index).subs(t, 2)) for index in range(5)
    )
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    discriminant_at_two = (4 * invariant_i**3 - invariant_j**2) / 27

    return {
        "coefficient_domain": "QQ(t)[x]",
        "product_degree": product.degree(),
        "remainder_degree": remainder.degree(),
        "visible_section_count": len(visible_identities),
        "visible_section_identities_exact": all(visible_identities),
        "extra_section_square_signatures": extra_signatures,
        "all_fifteen_sections_exact": all(visible_identities)
        and all(
            bool(signature["is_square"])
            for signature in extra_signatures.values()
        ),
        "generic_quartic_discriminant_nonzero": discriminant_at_two != 0,
        "discriminant_nonzero_witness_t": 2,
        "scope": (
            "point-on-curve identities only; the published numerical height "
            "determinant and independence are not recomputed"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(symbolic_verification(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
