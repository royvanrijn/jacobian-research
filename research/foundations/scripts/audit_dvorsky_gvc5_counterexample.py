#!/usr/bin/env python3
"""Dependency-free audit of the Dvorsky--Long GVC(5)/SIC(5) witness."""

from __future__ import annotations

import json
from math import factorial
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "dvorsky_gvc5_counterexample.json"
VARIABLES = ("t", "a", "b", "c", "d")
Exponent = tuple[int, int, int, int, int]
Polynomial = dict[Exponent, int]


def add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, 0) + coefficient
            if result[exponent] == 0:
                del result[exponent]
    return result


def scale(polynomial: Polynomial, scalar: int) -> Polynomial:
    return {
        exponent: scalar * coefficient
        for exponent, coefficient in polynomial.items()
        if scalar * coefficient
    }


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


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    result: Polynomial = {(0, 0, 0, 0, 0): 1}
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        remaining //= 2
    return result


def derivative(polynomial: Polynomial, variable: int) -> Polynomial:
    result: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        if exponent[variable] == 0:
            continue
        new_exponent = list(exponent)
        new_exponent[variable] -= 1
        result[tuple(new_exponent)] = coefficient * exponent[variable]
    return result


def mixed_derivative(polynomial: Polynomial, *variables: int) -> Polynomial:
    result = polynomial
    for variable in variables:
        result = derivative(result, variable)
    return result


def d_operator(polynomial: Polynomial) -> Polynomial:
    return add(
        mixed_derivative(polynomial, 1, 4),
        scale(mixed_derivative(polynomial, 2, 3), -1),
    )


def lambda_operator(polynomial: Polynomial) -> Polynomial:
    return derivative(d_operator(polynomial), 0)


def iterate(operator, polynomial: Polynomial, count: int) -> Polynomial:
    result = polynomial
    for _ in range(count):
        result = operator(result)
    return result


def monomial(coefficient: int, **powers: int) -> Polynomial:
    exponent = tuple(powers.get(variable, 0) for variable in VARIABLES)
    return {exponent: coefficient}


def main() -> None:
    # P=(t+c)(ad+bt)=tad+cad+bt^2+bct and Q=-c.
    p = add(
        monomial(1, t=1, a=1, d=1),
        monomial(1, c=1, a=1, d=1),
        monomial(1, b=1, t=2),
        monomial(1, b=1, c=1, t=1),
    )
    q = monomial(-1, c=1)

    checked_orders = []
    for m in range(1, 9):
        p_to_m = power(p, m)
        qp_to_m = multiply(q, p_to_m)

        expected_d_p = monomial(factorial(m) ** 2, c=m)
        expected_d_qp = add(
            monomial(
                factorial(m) ** 2 * (-1) ** (m + 1),
                t=m + 1,
            ),
            monomial(-factorial(m) ** 2, c=m + 1),
            monomial(factorial(m) ** 2, c=m, t=1),
        )
        assert iterate(d_operator, p_to_m, m) == expected_d_p
        assert iterate(d_operator, qp_to_m, m) == expected_d_qp
        assert iterate(lambda_operator, p_to_m, m) == {}

        if m == 1:
            expected_mixed = add(monomial(1, c=1), monomial(2, t=1))
        else:
            expected_mixed = monomial(
                (-1) ** (m + 1) * factorial(m + 1) * factorial(m) ** 2,
                t=1,
            )
        assert iterate(lambda_operator, qp_to_m, m) == expected_mixed
        checked_orders.append(m)

    artifact = {
        "format": "dvorsky-gvc5-counterexample-v1",
        "field": "characteristic zero",
        "variables": list(VARIABLES),
        "P": "(t+c)(ad+bt)",
        "Lambda": "d/dt (d/da d/dd - d/db d/dc)",
        "Q": "-c",
        "all_order_identities": {
            "D^m(P^m)": "(m!)^2 c^m",
            "D^m(Q P^m)": (
                "(m!)^2 [(-1)^(m+1) t^(m+1) - c^m(c-t)]"
            ),
            "Lambda^m(P^m)": "0 for every m>=1",
            "Lambda^m(Q P^m)": (
                "(-1)^(m+1) (m+1)! (m!)^2 t for every m>=2"
            ),
        },
        "bounded_sparse_polynomial_replay": {
            "orders": checked_orders,
            "implementation": "Python standard library sparse exact arithmetic",
        },
        "consequences": {
            "unrestricted_constant_coefficient_GVC": "fails in 5 variables",
            "SIC": "fails in 5 contraction pairs",
            "ordinary_Laplacian_GVC": "not affected; Lambda has order 3",
            "homogeneous_quartic_HN_VC": "not affected",
        },
        "provenance": {
            "SU2_seed": {
                "author": "Christopher D. Long",
                "url": "https://arxiv.org/abs/2607.19012",
            },
            "five_variable_GVC_SIC_lift": {
                "author": "Alexander Dvorsky",
                "date": "2026-07-23",
                "url": (
                    "https://sbseminar.wordpress.com/2026/07/20/"
                    "the-new-counterexample-to-the-jacobian-conjecture/"
                ),
            },
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS Dvorsky--Long GVC(5): exact sparse replay through m=8")
    print("PASS all-order formulas documented for GVC(5) and SIC(5)")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
