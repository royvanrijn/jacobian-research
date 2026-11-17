#!/usr/bin/env python3
"""Exact bounded replay for a homogeneous three-variable GVC counterexample.

The all-order proof is mathematical (Hopf/sphere coefficient extraction plus
Wick's formula).  This script independently verifies the polynomial identity,
homogeneity, and the first six pure and mixed differential identities using
integer arithmetic only.
"""

from __future__ import annotations

import json
from math import factorial, gcd
from pathlib import Path
from typing import Dict, Tuple

Exponent = Tuple[int, int, int]  # x, y, t
Polynomial = Dict[Exponent, int]

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "gvc3_homogeneous_counterexample.json"

ZERO: Exponent = (0, 0, 0)
X: Polynomial = {(1, 0, 0): 1}
Y: Polynomial = {(0, 1, 0): 1}
T: Polynomial = {(0, 0, 1): 1}


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
                for index in range(3)
            )
            result[exponent] = (
                result.get(exponent, 0)
                + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    result: Polynomial = {ZERO: 1}
    base = polynomial
    value = exponent
    while value:
        if value & 1:
            result = multiply(result, base)
        value >>= 1
        if value:
            base = multiply(base, base)
    return result


def monomial(exponent: Exponent, coefficient: int = 1) -> Polynomial:
    return {exponent: coefficient} if coefficient else {}


def factorial_product(exponent: Exponent) -> int:
    result = 1
    for entry in exponent:
        result *= factorial(entry)
    return result


def apolar_scalar(operator: Polynomial, polynomial: Polynomial) -> int:
    return sum(
        operator_coefficient
        * polynomial.get(exponent, 0)
        * factorial_product(exponent)
        for exponent, operator_coefficient in operator.items()
    )


def apply_operator(operator: Polynomial, polynomial: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for derivative, operator_coefficient in operator.items():
        for exponent, polynomial_coefficient in polynomial.items():
            if any(
                exponent[index] < derivative[index]
                for index in range(3)
            ):
                continue
            output_exponent = tuple(
                exponent[index] - derivative[index]
                for index in range(3)
            )
            falling_factorial = 1
            for index in range(3):
                falling_factorial *= (
                    factorial(exponent[index])
                    // factorial(output_exponent[index])
                )
            result[output_exponent] = (
                result.get(output_exponent, 0)
                + operator_coefficient
                * polynomial_coefficient
                * falling_factorial
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def homogeneous_degree(polynomial: Polynomial) -> int:
    degrees = {sum(exponent) for exponent in polynomial}
    assert len(degrees) == 1
    return next(iter(degrees))


def double_factorial(value: int) -> int:
    result = 1
    for entry in range(value, 0, -2):
        result *= entry
    return result


def gaussian_expectation(polynomial: Polynomial) -> int:
    """Expectation under x=X+iY, y=X-iY, t=T with X,Y,T iid N(0,1)."""
    result = 0
    for (x_degree, y_degree, t_degree), coefficient in polynomial.items():
        if x_degree != y_degree or t_degree % 2:
            continue
        circular = (2 ** x_degree) * factorial(x_degree)
        real = double_factorial(t_degree - 1) if t_degree else 1
        result += coefficient * circular * real
    return result


def expected_mixed_gaussian(order: int) -> int:
    numerator = (
        2 ** (2 * order)
        * factorial(2 * order)
        * double_factorial(12 * order + 3)
    )
    denominator = double_factorial(4 * order + 1)
    assert numerator % denominator == 0
    return numerator // denominator


def expected_mixed_detector(order: int) -> int:
    # Delta^(6m+1)(x^2 P^m)
    # = 2^(6m+1)(6m+1)! E[x^2 P(G)^m]
    # and E[x^2 P(G)^m]
    # = (12m+3)!! * integral_0^1 (1-v^2)^(2m) dv.
    numerator = (
        2 ** (8 * order + 1)
        * factorial(6 * order + 1)
        * factorial(2 * order)
        * double_factorial(12 * order + 3)
    )
    denominator = double_factorial(4 * order + 1)
    assert numerator % denominator == 0
    return numerator // denominator


def polynomial_terms(polynomial: Polynomial) -> list[dict[str, object]]:
    return [
        {
            "exponent_x_y_t": list(exponent),
            "coefficient": coefficient,
        }
        for exponent, coefficient in sorted(polynomial.items())
    ]


def main() -> None:
    # q = t^2 + xy, A = q + x^2,
    # C = y q^2 - 2 x t^2 q - x^3 t^2,
    # P = A C^2.
    q = add(monomial((0, 0, 2)), monomial((1, 1, 0)))
    x_squared = monomial((2, 0, 0))
    a = add(q, x_squared)
    c = add(
        multiply(Y, power(q, 2)),
        scale(multiply(multiply(X, power(T, 2)), q), -2),
        scale(monomial((3, 0, 2)), -1),
    )
    p = multiply(a, power(c, 2))

    # The polynomiality identity behind the homogeneous lift:
    # x C = q^3 - t^2 (q+x^2)^2.
    assert multiply(X, c) == add(
        power(q, 3),
        scale(multiply(power(T, 2), power(a, 2)), -1),
    )

    assert homogeneous_degree(c) == 5
    assert homogeneous_degree(p) == 12
    assert len(c) == 6
    assert len(p) == 23
    assert gcd(*(abs(value) for value in p.values())) == 1

    # In x=X+iY, y=X-iY coordinates, the ordinary real Laplacian is
    # Delta = 4 d_x d_y + d_t^2.  Lambda = Delta^6.
    delta_symbol = add(
        monomial((1, 1, 0), 4),
        monomial((0, 0, 2), 1),
    )

    checks: list[dict[str, object]] = []
    p_power: Polynomial = {ZERO: 1}
    for order in range(1, 7):
        p_power = multiply(p_power, p)
        delta_6m = power(delta_symbol, 6 * order)
        pure = apolar_scalar(delta_6m, p_power)
        assert pure == 0
        assert gaussian_expectation(p_power) == 0

        mixed_input = multiply(x_squared, p_power)
        gaussian_mixed = gaussian_expectation(mixed_input)
        assert gaussian_mixed == expected_mixed_gaussian(order)
        assert gaussian_mixed != 0
        mixed_output = apply_operator(delta_6m, mixed_input)
        assert mixed_output
        assert homogeneous_degree(mixed_output) == 2

        detector = apolar_scalar(
            power(delta_symbol, 6 * order + 1),
            mixed_input,
        )
        expected = expected_mixed_detector(order)
        assert detector == expected
        assert detector != 0
        assert apply_operator(delta_symbol, mixed_output) == {ZERO: expected}

        checks.append(
            {
                "m": order,
                "pure_Delta_6m_Pm": pure,
                "pure_gaussian_moment": 0,
                "mixed_gaussian_moment": gaussian_mixed,
                "mixed_output_terms": polynomial_terms(mixed_output),
                "Delta_6m_plus_1_x2_Pm": detector,
            }
        )

    artifact = {
        "format": "gvc3-homogeneous-counterexample-v1",
        "field": "characteristic zero",
        "variables": ["x", "y", "t"],
        "quadratic_form_q": "t^2+x*y",
        "polynomial": (
            "P=(q+x^2)*(y*q^2-2*x*t^2*q-x^3*t^2)^2"
        ),
        "polynomial_degree": 12,
        "polynomial_term_count": len(p),
        "operator": "Lambda=(4*d_x*d_y+d_t^2)^6",
        "multiplier": "Q=x^2",
        "all_order_claim": {
            "pure": "Lambda^m(P^m)=0 for every m>=1",
            "mixed": "Lambda^m(x^2*P^m)!=0 for every m>=1",
            "detector": (
                "Delta^(6m+1)(x^2*P^m)="
                "2^(8m+1)*(6m+1)!*(2m)!*(12m+3)!!/(4m+1)!!"
            ),
        },
        "bounded_exact_replay": checks,
        "expanded_P": polynomial_terms(p),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS polynomiality identity x*C=q^3-t^2*(q+x^2)^2")
    print("PASS P is primitive homogeneous degree 12 with 23 terms")
    print("PASS pure identities through m=6")
    print("PASS independent Gaussian moments through m=6")
    print("PASS mixed identities and closed detector through m=6")
    print(f"PASS wrote {OUTPUT}")


if __name__ == "__main__":
    main()
