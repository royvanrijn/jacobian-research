#!/usr/bin/env python3
"""Bounded exact search for a four-variable Dvorsky--Long descendant.

The searched family is

    P = (r*t + u*a + v*b + w*d) * (a*d + b*t),
    Lambda = d/dt * R(d/da, d/db, d/dd),

where (r,u,v,w) is induced by identifying c with a small linear form and
R is the general ternary quadratic form.  All arithmetic is integral and
the coefficient boxes are exhaustive as declared below.
"""

from __future__ import annotations

import json
from functools import reduce
from itertools import product
from math import factorial, gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "dvorsky_gvc4_bounded_search.json"
)
VARIABLES = ("t", "a", "b", "d")
Exponent = tuple[int, int, int, int]
Polynomial = dict[Exponent, int]

LINEAR_FORM_BOX = (-1, 0, 1)
QUADRATIC_FORM_BOX = (-2, -1, 0, 1, 2)
PURE_CUTOFF = 12
MIXED_WINDOW = tuple(range(5, PURE_CUTOFF + 1))

# a^2, ab, ad, b^2, bd, d^2, with one additional t derivative.
R_EXPONENTS: tuple[Exponent, ...] = (
    (0, 2, 0, 0),
    (0, 1, 1, 0),
    (0, 1, 0, 1),
    (0, 0, 2, 0),
    (0, 0, 1, 1),
    (0, 0, 0, 2),
)
LAMBDA_EXPONENTS = tuple(
    (exponent[0] + 1, *exponent[1:]) for exponent in R_EXPONENTS
)


def add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, 0) + coefficient
            if result[exponent] == 0:
                del result[exponent]
    return result


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(len(VARIABLES))
            )
            result[exponent] = (
                result.get(exponent, 0) + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def powers(polynomial: Polynomial, cutoff: int) -> list[Polynomial]:
    result: list[Polynomial] = [{(0, 0, 0, 0): 1}]
    for _ in range(cutoff):
        result.append(multiply(result[-1], polynomial))
    return result


def apolar_scalar(symbol: Polynomial, polynomial: Polynomial) -> int:
    return sum(
        symbol_coefficient
        * polynomial.get(exponent, 0)
        * factorial_product(exponent)
        for exponent, symbol_coefficient in symbol.items()
    )


def factorial_product(exponent: Exponent) -> int:
    result = 1
    for entry in exponent:
        result *= factorial(entry)
    return result


def linear_output(
    symbol: Polynomial, polynomial: Polynomial
) -> tuple[int, int, int, int]:
    """Return coefficients of symbol(partial) polynomial when degree is one."""
    result = []
    for output_variable in range(len(VARIABLES)):
        coefficient = 0
        for exponent, symbol_coefficient in symbol.items():
            source_exponent = tuple(
                entry + (index == output_variable)
                for index, entry in enumerate(exponent)
            )
            coefficient += (
                symbol_coefficient
                * polynomial.get(source_exponent, 0)
                * factorial_product(source_exponent)
            )
        result.append(coefficient)
    return tuple(result)  # type: ignore[return-value]


def monomial(variable: int) -> Polynomial:
    return {
        tuple(index == variable for index in range(len(VARIABLES))): 1
    }


def primitive_normalized_vectors() -> list[tuple[int, ...]]:
    result = []
    for coefficients in product(QUADRATIC_FORM_BOX, repeat=6):
        nonzero = [coefficient for coefficient in coefficients if coefficient]
        if not nonzero or nonzero[0] < 0:
            continue
        if reduce(gcd, (abs(coefficient) for coefficient in nonzero)) != 1:
            continue
        result.append(coefficients)
    return result


def make_p(linear_coefficients: tuple[int, int, int, int]) -> Polynomial:
    linear = {
        tuple(index == variable for index in range(4)): coefficient
        for variable, coefficient in enumerate(linear_coefficients)
        if coefficient
    }
    quadratic = {
        (0, 1, 0, 1): 1,  # ad
        (1, 0, 1, 0): 1,  # bt
    }
    return multiply(linear, quadratic)


def make_lambda(coefficients: tuple[int, ...]) -> Polynomial:
    return {
        exponent: coefficient
        for exponent, coefficient in zip(LAMBDA_EXPONENTS, coefficients)
        if coefficient
    }


def main() -> None:
    r_vectors = primitive_normalized_vectors()
    linear_vectors = [
        coefficients
        for coefficients in product(LINEAR_FORM_BOX, repeat=4)
        if any(coefficients)
        and next(value for value in coefficients if value) > 0
    ]

    tested_pairs = 0
    pure_four_survivors = 0
    pure_cutoff_survivors = 0
    persistent_linear_mixed_survivors = 0
    first_failure_histogram: dict[int, int] = {}
    delayed_examples: list[dict[str, object]] = []

    linear_data = [
        (
            linear_coefficients,
            powers(make_p(linear_coefficients), PURE_CUTOFF),
        )
        for linear_coefficients in linear_vectors
    ]

    for r_coefficients in r_vectors:
        lambda_polynomial = make_lambda(r_coefficients)
        lambda_powers = powers(lambda_polynomial, PURE_CUTOFF)

        for linear_coefficients, p_powers in linear_data:
            tested_pairs += 1
            pure_values: list[int] = []
            first_failure: int | None = None

            for order in range(1, PURE_CUTOFF + 1):
                value = apolar_scalar(
                    lambda_powers[order], p_powers[order]
                )
                pure_values.append(value)
                if value and first_failure is None:
                    first_failure = order
                    break

            if first_failure is not None and first_failure <= 4:
                continue

            pure_four_survivors += 1
            if first_failure is not None:
                first_failure_histogram[first_failure] = (
                    first_failure_histogram.get(first_failure, 0) + 1
                )
                if len(delayed_examples) < 8:
                    delayed_examples.append(
                        {
                            "linear_coefficients_r_u_v_w": linear_coefficients,
                            "R_coefficients_a2_ab_ad_b2_bd_d2": r_coefficients,
                            "pure_values_through_first_failure": pure_values,
                            "first_failure": first_failure,
                        }
                    )
                continue

            pure_cutoff_survivors += 1
            mixed_orders_by_coordinate: list[set[int]] = []
            for multiplier_variable in range(len(VARIABLES)):
                mixed_nonzero_orders: set[int] = set()
                for order in MIXED_WINDOW:
                    mixed_input = multiply(
                        monomial(multiplier_variable), p_powers[order]
                    )
                    output = linear_output(
                        lambda_powers[order], mixed_input
                    )
                    if any(output):
                        mixed_nonzero_orders.add(order)
                mixed_orders_by_coordinate.append(mixed_nonzero_orders)
            # Over characteristic zero, finitely many nonzero linear maps
            # have a common input outside the union of their kernels.  Thus
            # a fixed linear multiplier exists on the whole window exactly
            # when the mixed map is nonzero at each individual order.
            covered_orders = set().union(*mixed_orders_by_coordinate)
            if covered_orders == set(MIXED_WINDOW):
                persistent_linear_mixed_survivors += 1

    artifact = {
        "format": "dvorsky-gvc4-bounded-search-v1",
        "field": "integers (exact arithmetic)",
        "family": {
            "P": "(r*t+u*a+v*b+w*d)(a*d+b*t)",
            "Lambda": "d/dt R(d/da,d/db,d/dd)",
            "R": "A*a^2+B*a*b+C*a*d+D*b^2+E*b*d+F*d^2",
            "interpretation": (
                "c is identified with "
                "(r-1)t+u*a+v*b+w*d in the Dvorsky--Long polynomial"
            ),
        },
        "boxes": {
            "linear_coefficients_r_u_v_w": list(LINEAR_FORM_BOX),
            "quadratic_coefficients_A_through_F": list(
                QUADRATIC_FORM_BOX
            ),
            "normalization": (
                "both coefficient vectors are sign-normalized; R is primitive"
            ),
        },
        "cutoffs": {
            "pure_moments": [1, PURE_CUTOFF],
            "required_initial_scheme": [1, 4],
            "coordinate_mixed_window": [
                MIXED_WINDOW[0],
                MIXED_WINDOW[-1],
            ],
        },
        "counts": {
            "linear_forms": len(linear_vectors),
            "primitive_ternary_quadratic_forms": len(r_vectors),
            "pairs_tested": tested_pairs,
            "pure_m1_through_m4_survivors": pure_four_survivors,
            "pure_m1_through_m12_survivors": pure_cutoff_survivors,
            "m1_through_m12_survivors_admitting_a_fixed_linear_multiplier_nonzero_at_every_m5_through_m12": (
                persistent_linear_mixed_survivors
            ),
        },
        "first_pure_failure_after_m4": {
            str(order): count
            for order, count in sorted(first_failure_histogram.items())
        },
        "delayed_resonance_examples": delayed_examples,
        "conclusion": (
            "No four-variable witness occurs in this declared lattice "
            "slice through the stated cutoffs. This is a bounded negative "
            "search, not an all-order nonexistence theorem."
        ),
        "provenance": (
            "Search descends Long's SU(2) seed via Dvorsky's "
            "homogenization/lift."
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    assert tested_pairs == len(linear_vectors) * len(r_vectors)
    assert persistent_linear_mixed_survivors == 0
    print(
        "PASS bounded GVC(4) search:",
        tested_pairs,
        "exact normalized coefficient pairs",
    )
    print(
        "PASS pure survivors m<=4 / m<=12:",
        pure_four_survivors,
        "/",
        pure_cutoff_survivors,
    )
    print(
        "PASS fixed-linear-multiplier survivors on m=5..12:",
        persistent_linear_mixed_survivors,
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
