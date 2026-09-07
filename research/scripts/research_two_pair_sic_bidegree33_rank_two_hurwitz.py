#!/usr/bin/env python3
"""Explore the generic rank-two Hurwitz chart over a function field.

For an exact-rank-two cubic coefficient matrix write

    F(x,y) = (1+x) B(y) + x^2 (lambda+x) D(y),

where

    B = 1+a1*y+a2*y^2+a3*y^3,
    D = b0+b1*y+b2*y^2+d3*y^3.

The first moment eliminates

    d3 = -1-a1/3-lambda*b2/3.

This script constructs the remaining contraction moments directly from

    mu_m = sum_n n! (3m-n)! [x^n y^n] F(x,y)^m

over GF(p), exports them to Singular with lambda in the coefficient field
GF(p)(lambda), and optionally computes a standard basis.  Treating lambda
as a coefficient avoids the large extra elimination variable used by the
earlier reconnaissance.

This is an exploratory chart calculation, not a global rank-two theorem:
the binary-cubic-pencil boundary, the projective B_0=0 chart, and every
channel-minor open not explicitly requested remain separate.
"""

from __future__ import annotations

import argparse
from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial, gcd, lcm
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
PARAMETER_NAMES = ("lambda", "a1", "a2", "a3", "b0", "b1", "b2")
RING_VARIABLES = ("h", "k", "a1", "a2", "a3", "b0", "b1", "b2")
GENERIC_MSOLVE_VARIABLES = (
    "h", "k", "a1", "a2", "a3", "b0", "b1", "b2", "lambda"
)
COMBINED_FIXED_VARIABLES = ("r", "a1", "a2", "a3", "b0", "b1", "b2")
COMBINED_GENERIC_VARIABLES = (
    "r", "a1", "a2", "a3", "b0", "b1", "b2", "lambda"
)
MU2_ELIMINATION_VARIABLES = (
    "r", "a1", "a2", "b0", "b1", "b2", "lambda"
)
MU2_BOUNDARY_REDUCED_VARIABLES = (
    "r", "a1", "a3", "b0", "b2", "lambda"
)
MU2_BOUNDARY_SECONDARY_VARIABLES = (
    "r", "a1", "a2", "a3", "b0", "lambda"
)
MU2_BOUNDARY_TERTIARY_VARIABLES = (
    "r", "a1", "b0", "b2", "lambda"
)
MU2_GENERIC_REDUCED_VARIABLES = (
    "r", "a1", "b0", "b1", "b2", "lambda"
)
EXACT_LIFT_PRIME = (1 << 521) - 1
Exponent = tuple[int, int, int, int, int, int, int]
ParameterPolynomial = dict[Exponent, int]
RationalParameterPolynomial = dict[Exponent, Fraction]
YPolynomial = dict[int, ParameterPolynomial]
FixedExponent = tuple[int, int, int, int, int, int]
FixedPolynomial = dict[FixedExponent, int]
RingExponent = tuple[int, int, int, int, int, int, int, int]
RingPolynomial = dict[RingExponent, int]
GenericRingExponent = tuple[int, int, int, int, int, int, int, int, int]
GenericRingPolynomial = dict[GenericRingExponent, int]


def zero_exponent() -> Exponent:
    return (0, 0, 0, 0, 0, 0, 0)


def monomial(variable: int, coefficient: int, prime: int) -> ParameterPolynomial:
    exponent = [0] * len(PARAMETER_NAMES)
    exponent[variable] = 1
    return {tuple(exponent): coefficient % prime}


def add_parameter_polynomial(
    target: ParameterPolynomial,
    source: ParameterPolynomial,
    scale: int,
    prime: int,
) -> None:
    scale %= prime
    for exponent, coefficient in source.items():
        value = (target.get(exponent, 0) + scale * coefficient) % prime
        if value:
            target[exponent] = value
        else:
            target.pop(exponent, None)


def multiply_parameter_polynomials(
    left: ParameterPolynomial,
    right: ParameterPolynomial,
    prime: int,
) -> ParameterPolynomial:
    answer: ParameterPolynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(len(PARAMETER_NAMES))
            )
            value = (
                answer.get(exponent, 0)
                + left_coefficient * right_coefficient
            ) % prime
            if value:
                answer[exponent] = value
            else:
                answer.pop(exponent, None)
    return answer


def sum_parameter_polynomials(
    prime: int,
    *summands: tuple[int, ParameterPolynomial],
) -> ParameterPolynomial:
    answer: ParameterPolynomial = {}
    for scale, polynomial in summands:
        add_parameter_polynomial(answer, polynomial, scale, prime)
    return answer


def constant_parameter_polynomial(
    coefficient: int, prime: int
) -> ParameterPolynomial:
    coefficient %= prime
    return {zero_exponent(): coefficient} if coefficient else {}


def multiply_y_polynomials(
    left: YPolynomial,
    right: YPolynomial,
    prime: int,
) -> YPolynomial:
    answer: YPolynomial = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            product = multiply_parameter_polynomials(
                left_coefficient, right_coefficient, prime
            )
            bucket = answer.setdefault(degree, {})
            add_parameter_polynomial(bucket, product, 1, prime)
            if not bucket:
                del answer[degree]
    return answer


def powers(base: YPolynomial, maximum: int, prime: int) -> list[YPolynomial]:
    answer: list[YPolynomial] = [{0: {zero_exponent(): 1}}]
    for _ in range(maximum):
        answer.append(multiply_y_polynomials(answer[-1], base, prime))
    return answer


def base_polynomials(prime: int) -> tuple[YPolynomial, YPolynomial]:
    one = {zero_exponent(): 1}
    b_polynomial: YPolynomial = {
        0: one,
        1: monomial(1, 1, prime),
        2: monomial(2, 1, prime),
        3: monomial(3, 1, prime),
    }

    inverse_three = pow(3, -1, prime)
    d3: ParameterPolynomial = {zero_exponent(): -1 % prime}
    add_parameter_polynomial(
        d3, monomial(1, 1, prime), -inverse_three, prime
    )
    lambda_b2_exponent = [0] * len(PARAMETER_NAMES)
    lambda_b2_exponent[0] = 1
    lambda_b2_exponent[6] = 1
    add_parameter_polynomial(
        d3,
        {tuple(lambda_b2_exponent): 1},
        -inverse_three,
        prime,
    )
    d_polynomial: YPolynomial = {
        0: monomial(4, 1, prime),
        1: monomial(5, 1, prime),
        2: monomial(6, 1, prime),
        3: d3,
    }
    return b_polynomial, d_polynomial


def moment(
    order: int,
    b_powers: list[YPolynomial],
    d_powers: list[YPolynomial],
    prime: int,
) -> ParameterPolynomial:
    """Return the reduced mu_order in the seven parameter symbols."""

    answer: ParameterPolynomial = {}
    factorials = [factorial(index) % prime for index in range(3 * order + 1)]
    for b_power in range(order + 1):
        d_power = order - b_power
        y_polynomial = multiply_y_polynomials(
            b_powers[b_power], d_powers[d_power], prime
        )
        outer_scale = comb(order, b_power) % prime

        # (1+x)^b_power * x^(2*d_power) * (lambda+x)^d_power.
        x_terms: dict[int, ParameterPolynomial] = {}
        for first_degree in range(b_power + 1):
            first_coefficient = comb(b_power, first_degree)
            for second_degree in range(d_power + 1):
                degree = first_degree + 2 * d_power + second_degree
                lambda_exponent = d_power - second_degree
                exponent = [0] * len(PARAMETER_NAMES)
                exponent[0] = lambda_exponent
                coefficient = (
                    first_coefficient * comb(d_power, second_degree)
                ) % prime
                bucket = x_terms.setdefault(degree, {})
                add_parameter_polynomial(
                    bucket, {tuple(exponent): coefficient}, 1, prime
                )

        for degree, x_coefficient in x_terms.items():
            y_coefficient = y_polynomial.get(degree)
            if not y_coefficient:
                continue
            diagonal = multiply_parameter_polynomials(
                x_coefficient, y_coefficient, prime
            )
            weight = (
                outer_scale
                * factorials[degree]
                * factorials[3 * order - degree]
            ) % prime
            add_parameter_polynomial(answer, diagonal, weight, prime)
    return answer


def format_power(name: str, exponent: int) -> str:
    if exponent == 0:
        return ""
    if exponent == 1:
        return name
    return f"{name}^{exponent}"


def singular_polynomial(polynomial: ParameterPolynomial, prime: int) -> str:
    """Serialize a parameter polynomial with lambda as a coefficient."""

    if not polynomial:
        return "0"
    terms: list[str] = []
    for exponent in sorted(polynomial, reverse=True):
        coefficient = polynomial[exponent] % prime
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        lambda_factor = format_power("lambda", exponent[0])
        if lambda_factor:
            factors.append(lambda_factor)
        for name, power in zip(PARAMETER_NAMES[1:], exponent[1:]):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        terms.append("*".join(factors) if factors else "1")
    return "+".join(terms)


def evaluate_lambda(
    polynomial: ParameterPolynomial,
    value: int,
    prime: int,
) -> FixedPolynomial:
    answer: FixedPolynomial = {}
    value %= prime
    for exponent, coefficient in polynomial.items():
        fixed_exponent = exponent[1:]
        fixed_coefficient = coefficient * pow(value, exponent[0], prime)
        updated = (
            answer.get(fixed_exponent, 0) + fixed_coefficient
        ) % prime
        if updated:
            answer[fixed_exponent] = updated
        else:
            answer.pop(fixed_exponent, None)
    return answer


def symmetric_residue(value: int, prime: int) -> int:
    value %= prime
    return value if value <= prime // 2 else value - prime


def exact_moment_polynomials(
    orders: Iterable[int],
) -> dict[int, ParameterPolynomial]:
    """Recover the integer moments of 3F from one provably large prime.

    The coefficient l1 norm of

        3F = 3(1+x)B + x^2(lambda+x)(3D)

    is at most 52.  Hence every coefficient of its order-m contraction is
    bounded by (3m)!*52^m.  The 521-bit Mersenne prime is larger than twice
    this bound for every order accepted below, so symmetric representatives
    are the exact integer coefficients rather than a heuristic modular
    reconstruction.
    """

    maximum = max(orders)
    prime = EXACT_LIFT_PRIME
    if factorial(3 * maximum) * 52**maximum >= prime // 2:
        raise ValueError("the exact-lift coefficient bound is too small")
    b_base, d_base = base_polynomials(prime)
    b_powers = powers(b_base, maximum, prime)
    d_powers = powers(d_base, maximum, prime)
    answer: dict[int, ParameterPolynomial] = {}
    for order in orders:
        modular = moment(order, b_powers, d_powers, prime)
        scale = pow(3, order, prime)
        answer[order] = {
            exponent: symmetric_residue(coefficient * scale, prime)
            for exponent, coefficient in modular.items()
            if symmetric_residue(coefficient * scale, prime)
        }
        assert max(
            (abs(coefficient) for coefficient in answer[order].values()),
            default=0,
        ) <= factorial(3 * order) * 52**order
    return answer


def exact_parameter_polynomial(
    modular: ParameterPolynomial,
    prime: int = EXACT_LIFT_PRIME,
    scale: int = 1,
) -> ParameterPolynomial:
    return {
        exponent: symmetric_residue(coefficient * scale, prime)
        for exponent, coefficient in modular.items()
        if symmetric_residue(coefficient * scale, prime)
    }


def evaluate_lambda_exact(
    polynomial: ParameterPolynomial,
    value: int,
) -> FixedPolynomial:
    answer: FixedPolynomial = {}
    for exponent, coefficient in polynomial.items():
        fixed_exponent = exponent[1:]
        updated = (
            answer.get(fixed_exponent, 0)
            + coefficient * value ** exponent[0]
        )
        if updated:
            answer[fixed_exponent] = updated
        else:
            answer.pop(fixed_exponent, None)
    return answer


def exact_polynomial_string(
    polynomial: ParameterPolynomial,
) -> str:
    if not polynomial:
        return "0"
    terms: list[str] = []
    for exponent in sorted(polynomial, reverse=True):
        coefficient = polynomial[exponent]
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        lambda_factor = format_power("lambda", exponent[0])
        if lambda_factor:
            factors.append(lambda_factor)
        for name, power in zip(PARAMETER_NAMES[1:], exponent[1:]):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        terms.append("*".join(factors) if factors else "1")
    return "+".join(terms).replace("+-", "-")


def exact_fixed_ring_string(
    polynomial: FixedPolynomial,
    inverse_variable: int | None = None,
    subtract_one: bool = False,
) -> str:
    terms: dict[RingExponent, int] = {}
    for exponent, coefficient in polynomial.items():
        ring_exponent = [0] * len(RING_VARIABLES)
        if inverse_variable is not None:
            ring_exponent[inverse_variable] = 1
        ring_exponent[2:] = exponent
        terms[tuple(ring_exponent)] = coefficient
    if subtract_one:
        zero = (0,) * len(RING_VARIABLES)
        terms[zero] = terms.get(zero, 0) - 1
        if not terms[zero]:
            del terms[zero]
    if not terms:
        return "0"
    serialized: list[str] = []
    for exponent in sorted(terms, reverse=True):
        coefficient = terms[exponent]
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        for name, power in zip(RING_VARIABLES, exponent):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        serialized.append("*".join(factors) if factors else "1")
    return "+".join(serialized).replace("+-", "-")


def exact_msolve_source(
    moments: dict[int, ParameterPolynomial],
    minor_name: str,
    lambda_value: int,
) -> str:
    prime = EXACT_LIFT_PRIME
    discriminant = exact_parameter_polynomial(
        discriminant_polynomial(prime), prime
    )
    minor = exact_parameter_polynomial(
        channel_minor_polynomial(minor_name, prime),
        prime,
        3 if minor_name == "03" else 1,
    )
    generators = [
        exact_fixed_ring_string(
            evaluate_lambda_exact(discriminant, lambda_value),
            inverse_variable=0,
            subtract_one=True,
        ),
        exact_fixed_ring_string(
            evaluate_lambda_exact(minor, lambda_value),
            inverse_variable=1,
            subtract_one=True,
        ),
    ]
    for polynomial in moments.values():
        generators.append(
            exact_fixed_ring_string(
                evaluate_lambda_exact(polynomial, lambda_value)
            )
        )
    return (
        ",".join(RING_VARIABLES)
        + "\n0\n"
        + ",\n".join(generators)
        + "\n"
    )


def ring_polynomial_string(polynomial: RingPolynomial, prime: int) -> str:
    if not polynomial:
        return "0"
    terms: list[str] = []
    for exponent in sorted(polynomial, reverse=True):
        coefficient = polynomial[exponent] % prime
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        for name, power in zip(RING_VARIABLES, exponent):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        terms.append("*".join(factors) if factors else "1")
    return "+".join(terms)


def generic_ring_polynomial_string(
    polynomial: GenericRingPolynomial, prime: int
) -> str:
    if not polynomial:
        return "0"
    terms: list[str] = []
    for exponent in sorted(polynomial, reverse=True):
        coefficient = polynomial[exponent] % prime
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        for name, power in zip(GENERIC_MSOLVE_VARIABLES, exponent):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        terms.append("*".join(factors) if factors else "1")
    return "+".join(terms)


def lift_fixed_polynomial(
    polynomial: FixedPolynomial,
    prime: int,
    inverse_variable: int | None = None,
    subtract_one: bool = False,
) -> RingPolynomial:
    answer: RingPolynomial = {}
    for exponent, coefficient in polynomial.items():
        ring_exponent = [0] * len(RING_VARIABLES)
        if inverse_variable is not None:
            ring_exponent[inverse_variable] = 1
        ring_exponent[2:] = exponent
        answer[tuple(ring_exponent)] = coefficient % prime
    if subtract_one:
        zero = (0,) * len(RING_VARIABLES)
        value = (answer.get(zero, 0) - 1) % prime
        if value:
            answer[zero] = value
        else:
            answer.pop(zero, None)
    return answer


def lift_parameter_polynomial(
    polynomial: ParameterPolynomial,
    prime: int,
    inverse_variable: int | None = None,
    subtract_one: bool = False,
) -> GenericRingPolynomial:
    answer: GenericRingPolynomial = {}
    for exponent, coefficient in polynomial.items():
        ring_exponent = [0] * len(GENERIC_MSOLVE_VARIABLES)
        if inverse_variable is not None:
            ring_exponent[inverse_variable] = 1
        # Parameter order is lambda,a1,a2,a3,b0,b1,b2.
        ring_exponent[2:8] = exponent[1:]
        ring_exponent[8] = exponent[0]
        answer[tuple(ring_exponent)] = coefficient % prime
    if subtract_one:
        zero = (0,) * len(GENERIC_MSOLVE_VARIABLES)
        value = (answer.get(zero, 0) - 1) % prime
        if value:
            answer[zero] = value
        else:
            answer.pop(zero, None)
    return answer


def discriminant_polynomial(prime: int) -> ParameterPolynomial:
    one = constant_parameter_polynomial(1, prime)
    lambda_polynomial = monomial(0, 1, prime)
    a1 = monomial(1, 1, prime)
    a2 = monomial(2, 1, prime)
    b1 = monomial(5, 1, prime)
    b2 = monomial(6, 1, prime)
    lambda_b1 = multiply_parameter_polynomials(
        lambda_polynomial, b1, prime
    )
    lambda_b2 = multiply_parameter_polynomials(
        lambda_polynomial, b2, prime
    )
    lambda_a1 = multiply_parameter_polynomials(
        lambda_polynomial, a1, prime
    )
    lambda_squared_b2 = multiply_parameter_polynomials(
        lambda_polynomial, lambda_b2, prime
    )
    first = sum_parameter_polynomials(
        prime, (9, one), (2, a1), (1, lambda_b2)
    )
    second = sum_parameter_polynomials(
        prime, (3, one), (2, lambda_b1), (3, b2)
    )
    third = sum_parameter_polynomials(
        prime,
        (3, a1),
        (-1, lambda_a1),
        (2, a2),
        (-3, lambda_polynomial),
        (-1, lambda_squared_b2),
    )
    return sum_parameter_polynomials(
        prime,
        (
            1,
            multiply_parameter_polynomials(first, first, prime),
        ),
        (
            1,
            multiply_parameter_polynomials(second, third, prime),
        ),
    )


def channel_minor_polynomial(name: str, prime: int) -> ParameterPolynomial:
    a1 = monomial(1, 1, prime)
    a2 = monomial(2, 1, prime)
    a3 = monomial(3, 1, prime)
    b0 = monomial(4, 1, prime)
    b1 = monomial(5, 1, prime)
    b2 = monomial(6, 1, prime)
    if name == "01":
        return sum_parameter_polynomials(
            prime,
            (1, b1),
            (-1, multiply_parameter_polynomials(a1, b0, prime)),
        )
    if name == "02":
        return sum_parameter_polynomials(
            prime,
            (1, b2),
            (-1, multiply_parameter_polynomials(a2, b0, prime)),
        )
    if name == "03":
        _, d_base = base_polynomials(prime)
        d3 = d_base[3]
        return sum_parameter_polynomials(
            prime,
            (1, d3),
            (-1, multiply_parameter_polynomials(a3, b0, prime)),
        )
    raise ValueError(name)


def msolve_source(
    moments: dict[int, ParameterPolynomial],
    prime: int,
    minor_name: str,
    lambda_value: int | None,
) -> str:
    if lambda_value is None:
        generators = [
            generic_ring_polynomial_string(
                lift_parameter_polynomial(
                    discriminant_polynomial(prime),
                    prime,
                    inverse_variable=0,
                    subtract_one=True,
                ),
                prime,
            ),
            generic_ring_polynomial_string(
                lift_parameter_polynomial(
                    channel_minor_polynomial(minor_name, prime),
                    prime,
                    inverse_variable=1,
                    subtract_one=True,
                ),
                prime,
            ),
        ]
        for polynomial in moments.values():
            generators.append(
                generic_ring_polynomial_string(
                    lift_parameter_polynomial(polynomial, prime), prime
                )
            )
        return (
            ",".join(GENERIC_MSOLVE_VARIABLES)
            + f"\n{prime}\n"
            + ",\n".join(generators)
            + "\n"
        )

    fixed_discriminant = evaluate_lambda(
        discriminant_polynomial(prime), lambda_value, prime
    )
    fixed_minor = evaluate_lambda(
        channel_minor_polynomial(minor_name, prime),
        lambda_value,
        prime,
    )
    generators = [
        ring_polynomial_string(
            lift_fixed_polynomial(
                fixed_discriminant,
                prime,
                inverse_variable=0,
                subtract_one=True,
            ),
            prime,
        ),
        ring_polynomial_string(
            lift_fixed_polynomial(
                fixed_minor,
                prime,
                inverse_variable=1,
                subtract_one=True,
            ),
            prime,
        ),
    ]
    for polynomial in moments.values():
        fixed = evaluate_lambda(polynomial, lambda_value, prime)
        generators.append(
            ring_polynomial_string(
                lift_fixed_polynomial(fixed, prime), prime
            )
        )
    return (
        ",".join(RING_VARIABLES)
        + f"\n{prime}\n"
        + ",\n".join(generators)
        + "\n"
    )


def combined_ring_polynomial_string(
    polynomial: ParameterPolynomial,
    prime: int,
    lambda_value: int | None,
    multiply_by_r: bool = False,
    subtract_one: bool = False,
) -> str:
    terms: dict[tuple[int, ...], int] = {}
    if lambda_value is None:
        variable_names = COMBINED_GENERIC_VARIABLES
        for exponent, coefficient in polynomial.items():
            ring_exponent = [0] * len(variable_names)
            ring_exponent[0] = int(multiply_by_r)
            ring_exponent[1:7] = exponent[1:]
            ring_exponent[7] = exponent[0]
            value = (
                terms.get(tuple(ring_exponent), 0) + coefficient
            ) % prime
            if value:
                terms[tuple(ring_exponent)] = value
            else:
                terms.pop(tuple(ring_exponent), None)
    else:
        variable_names = COMBINED_FIXED_VARIABLES
        fixed = evaluate_lambda(polynomial, lambda_value, prime)
        for exponent, coefficient in fixed.items():
            ring_exponent = [0] * len(variable_names)
            ring_exponent[0] = int(multiply_by_r)
            ring_exponent[1:] = exponent
            terms[tuple(ring_exponent)] = coefficient % prime
    if subtract_one:
        zero = (0,) * len(variable_names)
        value = (terms.get(zero, 0) - 1) % prime
        if value:
            terms[zero] = value
        else:
            terms.pop(zero, None)
    if not terms:
        return "0"
    serialized: list[str] = []
    for exponent in sorted(terms, reverse=True):
        coefficient = terms[exponent]
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        for name, power in zip(variable_names, exponent):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        serialized.append("*".join(factors) if factors else "1")
    return "+".join(serialized)


def combined_msolve_source(
    moments: dict[int, ParameterPolynomial],
    prime: int,
    minor_name: str,
    lambda_value: int | None,
) -> str:
    open_product = multiply_parameter_polynomials(
        discriminant_polynomial(prime),
        channel_minor_polynomial(minor_name, prime),
        prime,
    )
    generators = [
        combined_ring_polynomial_string(
            open_product,
            prime,
            lambda_value,
            multiply_by_r=True,
            subtract_one=True,
        )
    ]
    generators.extend(
        combined_ring_polynomial_string(
            polynomial, prime, lambda_value
        )
        for polynomial in moments.values()
    )
    variable_names = (
        COMBINED_GENERIC_VARIABLES
        if lambda_value is None
        else COMBINED_FIXED_VARIABLES
    )
    return (
        ",".join(variable_names)
        + f"\n{prime}\n"
        + ",\n".join(generators)
        + "\n"
    )


def multiply_fixed_polynomials(
    left: FixedPolynomial,
    right: FixedPolynomial,
    prime: int,
) -> FixedPolynomial:
    answer: FixedPolynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(len(left_exponent))
            )
            value = (
                answer.get(exponent, 0)
                + left_coefficient * right_coefficient
            ) % prime
            if value:
                answer[exponent] = value
            else:
                answer.pop(exponent, None)
    return answer


def parameter_polynomial_power(
    polynomial: ParameterPolynomial,
    exponent: int,
    prime: int,
) -> ParameterPolynomial:
    answer = constant_parameter_polynomial(1, prime)
    for _ in range(exponent):
        answer = multiply_parameter_polynomials(answer, polynomial, prime)
    return answer


def substitute_parameter_variable(
    polynomial: ParameterPolynomial,
    variable: int,
    replacement: ParameterPolynomial,
    prime: int,
) -> ParameterPolynomial:
    """Substitute one parameter by a polynomial in the other parameters."""

    maximum_power = max(
        (exponent[variable] for exponent in polynomial), default=0
    )
    replacement_powers = [
        parameter_polynomial_power(replacement, power, prime)
        for power in range(maximum_power + 1)
    ]
    answer: ParameterPolynomial = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[variable]
        residual_exponent = list(exponent)
        residual_exponent[variable] = 0
        residual = {tuple(residual_exponent): coefficient}
        summand = multiply_parameter_polynomials(
            residual, replacement_powers[power], prime
        )
        add_parameter_polynomial(answer, summand, 1, prime)
    return answer


def split_linear_parameter_polynomial(
    polynomial: ParameterPolynomial,
    variable: int,
    prime: int,
) -> tuple[ParameterPolynomial, ParameterPolynomial]:
    """Write a polynomial as variable*pivot + rest."""

    pivot: ParameterPolynomial = {}
    rest: ParameterPolynomial = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[variable]
        if power not in (0, 1):
            raise ValueError("requested polynomial is not linear")
        reduced_exponent = list(exponent)
        reduced_exponent[variable] = 0
        target = pivot if power else rest
        add_parameter_polynomial(
            target, {tuple(reduced_exponent): coefficient}, 1, prime
        )
    return pivot, rest


def linear_substitution_numerator(
    polynomial: ParameterPolynomial,
    variable: int,
    pivot: ParameterPolynomial,
    rest: ParameterPolynomial,
    prime: int,
) -> ParameterPolynomial:
    """Substitute variable=-rest/pivot and clear the minimal denominator."""

    groups: dict[int, ParameterPolynomial] = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[variable]
        reduced_exponent = list(exponent)
        reduced_exponent[variable] = 0
        add_parameter_polynomial(
            groups.setdefault(power, {}),
            {tuple(reduced_exponent): coefficient},
            1,
            prime,
        )
    clearing_power = max(groups, default=0)
    pivot_powers = [
        parameter_polynomial_power(pivot, power, prime)
        for power in range(clearing_power + 1)
    ]
    rest_powers = [
        parameter_polynomial_power(rest, power, prime)
        for power in range(clearing_power + 1)
    ]
    answer: ParameterPolynomial = {}
    for power, group in groups.items():
        multiplier = multiply_parameter_polynomials(
            rest_powers[power],
            pivot_powers[clearing_power - power],
            prime,
        )
        summand = multiply_parameter_polynomials(group, multiplier, prime)
        add_parameter_polynomial(
            answer, summand, -1 if power % 2 else 1, prime
        )
    return answer


def divide_parameter_polynomials(
    dividend: ParameterPolynomial,
    divisor: ParameterPolynomial,
    prime: int,
) -> ParameterPolynomial:
    """Return the exact sparse quotient for a known polynomial factor."""

    remainder = {
        exponent: coefficient % prime
        for exponent, coefficient in dividend.items()
        if coefficient % prime
    }
    quotient: ParameterPolynomial = {}
    divisor_lead = max(divisor)
    divisor_coefficient = divisor[divisor_lead] % prime
    while remainder:
        remainder_lead = max(remainder)
        if any(
            left < right
            for left, right in zip(remainder_lead, divisor_lead)
        ):
            raise ValueError("polynomial division has a remainder")
        exponent = tuple(
            left - right
            for left, right in zip(remainder_lead, divisor_lead)
        )
        coefficient = (
            remainder[remainder_lead]
            * pow(divisor_coefficient, -1, prime)
        ) % prime
        quotient[exponent] = (
            quotient.get(exponent, 0) + coefficient
        ) % prime
        for divisor_exponent, divisor_value in divisor.items():
            target = tuple(
                left + right
                for left, right in zip(exponent, divisor_exponent)
            )
            value = (
                remainder.get(target, 0)
                - coefficient * divisor_value
            ) % prime
            if value:
                remainder[target] = value
            else:
                remainder.pop(target, None)
    return quotient


def strip_parameter_factor(
    polynomial: ParameterPolynomial,
    factor: ParameterPolynomial,
    prime: int,
) -> tuple[ParameterPolynomial, int]:
    """Cancel every exact power of a known-nonzero modular factor."""

    answer = polynomial
    exponent = 0
    while answer:
        try:
            quotient = divide_parameter_polynomials(answer, factor, prime)
        except ValueError:
            break
        answer = quotient
        exponent += 1
    return answer, exponent


def coefficient_groups(
    polynomial: ParameterPolynomial,
    variable: int,
    prime: int,
) -> dict[int, ParameterPolynomial]:
    groups: dict[int, ParameterPolynomial] = {}
    for exponent, coefficient in polynomial.items():
        reduced_exponent = list(exponent)
        power = reduced_exponent[variable]
        reduced_exponent[variable] = 0
        add_parameter_polynomial(
            groups.setdefault(power, {}),
            {tuple(reduced_exponent): coefficient},
            1,
            prime,
        )
    return groups


def quadratic_remainder_numerator(
    polynomial: ParameterPolynomial,
    variable: int,
    quadratic: tuple[
        ParameterPolynomial,
        ParameterPolynomial,
        ParameterPolynomial,
    ],
    prime: int,
) -> tuple[ParameterPolynomial, ParameterPolynomial]:
    """Reduce modulo A*x^2+B*x+C, clearing a power of A."""

    a_polynomial, b_polynomial, c_polynomial = quadratic
    groups = coefficient_groups(polynomial, variable, prime)
    maximum = max(groups, default=0)
    zero: ParameterPolynomial = {}
    one = constant_parameter_polynomial(1, prime)
    if maximum == 0:
        return groups.get(0, {}), zero
    if maximum == 1:
        return groups.get(0, {}), groups.get(1, {})

    # R_k=A^(k-1)*x^k modulo A*x^2+B*x+C.
    remainders: list[
        tuple[ParameterPolynomial, ParameterPolynomial]
    ] = [(one, zero), (zero, one)]
    a_times_c = multiply_parameter_polynomials(
        a_polynomial, c_polynomial, prime
    )
    for _ in range(2, maximum + 1):
        previous_constant, previous_linear = remainders[-1]
        earlier_constant, earlier_linear = remainders[-2]
        constant = multiply_parameter_polynomials(
            b_polynomial, previous_constant, prime
        )
        add_parameter_polynomial(
            constant,
            multiply_parameter_polynomials(
                a_times_c, earlier_constant, prime
            ),
            1,
            prime,
        )
        constant = {
            exponent: (-coefficient) % prime
            for exponent, coefficient in constant.items()
        }
        linear = multiply_parameter_polynomials(
            b_polynomial, previous_linear, prime
        )
        add_parameter_polynomial(
            linear,
            multiply_parameter_polynomials(
                a_times_c, earlier_linear, prime
            ),
            1,
            prime,
        )
        linear = {
            exponent: (-coefficient) % prime
            for exponent, coefficient in linear.items()
        }
        remainders.append((constant, linear))

    a_powers = [
        parameter_polynomial_power(a_polynomial, power, prime)
        for power in range(maximum)
    ]
    answer_constant: ParameterPolynomial = {}
    answer_linear: ParameterPolynomial = {}
    for power, group in groups.items():
        if power <= 1:
            multiplier = a_powers[maximum - 1]
        else:
            multiplier = a_powers[maximum - power]
        group_multiplier = multiply_parameter_polynomials(
            group, multiplier, prime
        )
        reduced_constant, reduced_linear = remainders[power]
        add_parameter_polynomial(
            answer_constant,
            multiply_parameter_polynomials(
                group_multiplier, reduced_constant, prime
            ),
            1,
            prime,
        )
        add_parameter_polynomial(
            answer_linear,
            multiply_parameter_polynomials(
                group_multiplier, reduced_linear, prime
            ),
            1,
            prime,
        )
    return answer_constant, answer_linear


def affine_parameter_polynomial(
    constant: ParameterPolynomial,
    linear: ParameterPolynomial,
    variable: int,
    prime: int,
) -> ParameterPolynomial:
    answer = dict(constant)
    shifted: ParameterPolynomial = {}
    for exponent, coefficient in linear.items():
        shifted_exponent = list(exponent)
        shifted_exponent[variable] += 1
        shifted[tuple(shifted_exponent)] = coefficient
    add_parameter_polynomial(answer, shifted, 1, prime)
    return answer


def mu2_generic_reduced_msolve_source(
    moments: dict[int, ParameterPolynomial],
    prime: int,
    minor_name: str,
    branch: str,
    f4sat: bool,
) -> tuple[str, dict[str, object]]:
    """Reduce the first-mu2-pivot-open chart by mu3 and mu4."""

    if minor_name != "01":
        raise ValueError("the reduced generic modes use minor 01")
    first_pivot, first_rest = split_linear_parameter_polynomial(
        moments[2], 3, prime
    )
    eliminated = {
        order: linear_substitution_numerator(
            polynomial, 3, first_pivot, first_rest, prime
        )
        for order, polynomial in moments.items()
        if order > 2
    }
    mu3_groups = coefficient_groups(eliminated[3], 2, prime)
    if set(mu3_groups) != {0, 1, 2}:
        raise ValueError("mu3 is not a genuine quadratic in a2")
    c_polynomial = mu3_groups[0]
    b_polynomial = mu3_groups[1]
    a_polynomial = mu3_groups[2]
    discriminant = discriminant_polynomial(prime)
    minor = channel_minor_polynomial(minor_name, prime)

    if branch in ("quadratic", "quadratic-boundary"):
        reduced = {
            order: quadratic_remainder_numerator(
                polynomial,
                2,
                (a_polynomial, b_polynomial, c_polynomial),
                prime,
            )
            for order, polynomial in eliminated.items()
            if order > 3
        }
        u4, v4 = reduced[4]
        open_product = multiply_parameter_polynomials(
            first_pivot, a_polynomial, prime
        )
        open_product = multiply_parameter_polynomials(
            open_product, minor, prime
        )
        if branch == "quadratic":
            discriminant_constant, discriminant_linear = (
                quadratic_remainder_numerator(
                    discriminant,
                    2,
                    (a_polynomial, b_polynomial, c_polynomial),
                    prime,
                )
            )
            reduced_discriminant = multiply_parameter_polynomials(
                discriminant_constant, v4, prime
            )
            add_parameter_polynomial(
                reduced_discriminant,
                multiply_parameter_polynomials(
                    discriminant_linear, u4, prime
                ),
                -1,
                prime,
            )
            open_product = multiply_parameter_polynomials(
                open_product, reduced_discriminant, prime
            )
            open_product = multiply_parameter_polynomials(
                open_product, v4, prime
            )
            a_u4_squared = multiply_parameter_polynomials(
                a_polynomial,
                multiply_parameter_polynomials(u4, u4, prime),
                prime,
            )
            b_u4_v4 = multiply_parameter_polynomials(
                b_polynomial,
                multiply_parameter_polynomials(u4, v4, prime),
                prime,
            )
            c_v4_squared = multiply_parameter_polynomials(
                c_polynomial,
                multiply_parameter_polynomials(v4, v4, prime),
                prime,
            )
            quadratic_equation = dict(a_u4_squared)
            add_parameter_polynomial(
                quadratic_equation, b_u4_v4, -1, prime
            )
            add_parameter_polynomial(
                quadratic_equation, c_v4_squared, 1, prime
            )
            later_equations: dict[int, ParameterPolynomial] = {}
            for order, (constant, linear) in reduced.items():
                if order == 4:
                    continue
                equation = multiply_parameter_polynomials(
                    constant, v4, prime
                )
                add_parameter_polynomial(
                    equation,
                    multiply_parameter_polynomials(linear, u4, prime),
                    -1,
                    prime,
                )
                later_equations[order] = equation
            variables = MU2_GENERIC_REDUCED_VARIABLES
            active = (1, 4, 5, 6, 0)
            generators = [
                active_parameter_polynomial_string(
                    open_product,
                    prime,
                    variables,
                    active,
                    multiply_by_r=True,
                    subtract_one=True,
                ),
                active_parameter_polynomial_string(
                    quadratic_equation, prime, variables, active
                ),
            ]
            generators.extend(
                active_parameter_polynomial_string(
                    polynomial, prime, variables, active
                )
                for polynomial in later_equations.values()
            )
            metadata: dict[str, object] = {
                "branch": "first_pivot_open_mu3_quadratic_mu4_pivot_open",
                "mu3_leading_terms": len(a_polynomial),
                "mu4_pivot_terms": len(v4),
                "equation_terms": {
                    "mu3_after_mu4": len(quadratic_equation),
                    **{
                        str(order): len(polynomial)
                        for order, polynomial in later_equations.items()
                    },
                },
            }
        else:
            open_product = multiply_parameter_polynomials(
                open_product, discriminant, prime
            )
            variables = MU2_ELIMINATION_VARIABLES
            active = (1, 2, 4, 5, 6, 0)
            quadratic_equation = affine_parameter_polynomial(
                c_polynomial, b_polynomial, 2, prime
            )
            a2_squared = monomial(2, 1, prime)
            a2_squared = multiply_parameter_polynomials(
                a2_squared, a2_squared, prime
            )
            add_parameter_polynomial(
                quadratic_equation,
                multiply_parameter_polynomials(
                    a_polynomial, a2_squared, prime
                ),
                1,
                prime,
            )
            generators = [
                active_parameter_polynomial_string(
                    open_product,
                    prime,
                    variables,
                    active,
                    multiply_by_r=True,
                    subtract_one=True,
                ),
                active_parameter_polynomial_string(
                    quadratic_equation, prime, variables, active
                ),
                active_parameter_polynomial_string(
                    u4, prime, variables, active
                ),
                active_parameter_polynomial_string(
                    v4, prime, variables, active
                ),
            ]
            generators.extend(
                active_parameter_polynomial_string(
                    affine_parameter_polynomial(
                        constant, linear, 2, prime
                    ),
                    prime,
                    variables,
                    active,
                )
                for order, (constant, linear) in reduced.items()
                if order > 4
            )
            metadata = {
                "branch": (
                    "first_pivot_open_mu3_quadratic_mu4_pivot_zero"
                ),
                "mu3_leading_terms": len(a_polynomial),
                "mu4_constant_terms": len(u4),
                "mu4_pivot_terms": len(v4),
            }
    elif branch == "linear":
        linear_eliminated = {
            order: linear_substitution_numerator(
                polynomial, 2, b_polynomial, c_polynomial, prime
            )
            for order, polynomial in eliminated.items()
            if order > 3
        }
        reduced_discriminant = linear_substitution_numerator(
            discriminant, 2, b_polynomial, c_polynomial, prime
        )
        open_product = multiply_parameter_polynomials(
            first_pivot, b_polynomial, prime
        )
        open_product = multiply_parameter_polynomials(
            open_product, reduced_discriminant, prime
        )
        open_product = multiply_parameter_polynomials(
            open_product, minor, prime
        )
        variables = MU2_GENERIC_REDUCED_VARIABLES
        active = (1, 4, 5, 6, 0)
        generators = [
            active_parameter_polynomial_string(
                open_product,
                prime,
                variables,
                active,
                multiply_by_r=True,
                subtract_one=True,
            ),
            active_parameter_polynomial_string(
                a_polynomial, prime, variables, active
            ),
        ]
        generators.extend(
            active_parameter_polynomial_string(
                polynomial, prime, variables, active
            )
            for polynomial in linear_eliminated.values()
        )
        metadata = {
            "branch": "first_pivot_open_mu3_linear_pivot_open",
            "mu3_leading_terms": len(a_polynomial),
            "mu3_linear_terms": len(b_polynomial),
            "later_moment_terms": {
                str(order): len(polynomial)
                for order, polynomial in linear_eliminated.items()
            },
        }
    elif branch == "linear-boundary":
        open_product = multiply_parameter_polynomials(
            first_pivot, discriminant, prime
        )
        open_product = multiply_parameter_polynomials(
            open_product, minor, prime
        )
        variables = MU2_ELIMINATION_VARIABLES
        active = (1, 2, 4, 5, 6, 0)
        generators = [
            active_parameter_polynomial_string(
                open_product,
                prime,
                variables,
                active,
                multiply_by_r=True,
                subtract_one=True,
            ),
            active_parameter_polynomial_string(
                a_polynomial, prime, variables, active
            ),
            active_parameter_polynomial_string(
                b_polynomial, prime, variables, active
            ),
            active_parameter_polynomial_string(
                c_polynomial, prime, variables, active
            ),
        ]
        generators.extend(
            active_parameter_polynomial_string(
                polynomial, prime, variables, active
            )
            for order, polynomial in eliminated.items()
            if order > 3
        )
        metadata = {
            "branch": "first_pivot_open_mu3_degree_drop_boundary",
            "mu3_groups": {
                "quadratic": len(a_polynomial),
                "linear": len(b_polynomial),
                "constant": len(c_polynomial),
            },
        }
    else:
        raise ValueError(f"unknown reduced generic branch {branch}")

    if f4sat:
        open_generator = active_parameter_polynomial_string(
            open_product, prime, variables, active
        )
        generators = [*generators[1:], open_generator]
        source_variables = variables[1:]
    else:
        source_variables = variables
    source = (
        ",".join(source_variables)
        + f"\n{prime}\n"
        + ",\n".join(generators)
        + "\n"
    )
    metadata["f4sat"] = f4sat
    return source, metadata


def add_fixed_polynomial(
    target: FixedPolynomial,
    source: FixedPolynomial,
    scale: int,
    prime: int,
) -> None:
    for exponent, coefficient in source.items():
        value = (
            target.get(exponent, 0) + scale * coefficient
        ) % prime
        if value:
            target[exponent] = value
        else:
            target.pop(exponent, None)


def without_a3_exponent(exponent: Exponent) -> FixedExponent:
    return exponent[:3] + exponent[4:]


def eliminate_mu2_a3(
    moments: dict[int, ParameterPolynomial],
    prime: int,
) -> tuple[
    dict[int, FixedPolynomial],
    FixedPolynomial,
    FixedPolynomial,
]:
    """Substitute a3=-rest/pivot and clear only the necessary power."""

    mu2 = moments[2]
    pivot: FixedPolynomial = {}
    rest: FixedPolynomial = {}
    for exponent, coefficient in mu2.items():
        a3_power = exponent[3]
        if a3_power not in (0, 1):
            raise ValueError("mu2 is not linear in a3")
        target = pivot if a3_power else rest
        reduced_exponent = without_a3_exponent(exponent)
        value = (
            target.get(reduced_exponent, 0) + coefficient
        ) % prime
        if value:
            target[reduced_exponent] = value
        else:
            target.pop(reduced_exponent, None)

    maximum_power = max(
        exponent[3]
        for order, polynomial in moments.items()
        if order > 2
        for exponent in polynomial
    )
    one: FixedPolynomial = {(0,) * 6: 1}
    pivot_powers = [one]
    rest_powers = [one]
    for _ in range(maximum_power):
        pivot_powers.append(
            multiply_fixed_polynomials(pivot_powers[-1], pivot, prime)
        )
        rest_powers.append(
            multiply_fixed_polynomials(rest_powers[-1], rest, prime)
        )

    eliminated: dict[int, FixedPolynomial] = {}
    for order, polynomial in moments.items():
        if order == 2:
            continue
        groups: dict[int, FixedPolynomial] = {}
        for exponent, coefficient in polynomial.items():
            a3_power = exponent[3]
            reduced_exponent = without_a3_exponent(exponent)
            group = groups.setdefault(a3_power, {})
            value = (
                group.get(reduced_exponent, 0) + coefficient
            ) % prime
            if value:
                group[reduced_exponent] = value
            else:
                group.pop(reduced_exponent, None)
        clearing_power = max(groups)
        numerator: FixedPolynomial = {}
        for a3_power, group in groups.items():
            multiplier = multiply_fixed_polynomials(
                rest_powers[a3_power],
                pivot_powers[clearing_power - a3_power],
                prime,
            )
            summand = multiply_fixed_polynomials(
                group, multiplier, prime
            )
            add_fixed_polynomial(
                numerator,
                summand,
                -1 if a3_power % 2 else 1,
                prime,
            )
        eliminated[order] = numerator
    return eliminated, pivot, rest


def remove_a3(
    polynomial: ParameterPolynomial,
) -> FixedPolynomial:
    answer: FixedPolynomial = {}
    for exponent, coefficient in polynomial.items():
        if exponent[3]:
            raise ValueError("polynomial still depends on a3")
        reduced = without_a3_exponent(exponent)
        answer[reduced] = coefficient
    return answer


def mu2_eliminated_polynomial_string(
    polynomial: FixedPolynomial,
    prime: int,
    multiply_by_r: bool = False,
    subtract_one: bool = False,
) -> str:
    terms: dict[tuple[int, ...], int] = {}
    for exponent, coefficient in polynomial.items():
        ring_exponent = [0] * len(MU2_ELIMINATION_VARIABLES)
        ring_exponent[0] = int(multiply_by_r)
        # Reduced exponent order: lambda,a1,a2,b0,b1,b2.
        ring_exponent[1:6] = exponent[1:]
        ring_exponent[6] = exponent[0]
        terms[tuple(ring_exponent)] = coefficient % prime
    if subtract_one:
        zero = (0,) * len(MU2_ELIMINATION_VARIABLES)
        value = (terms.get(zero, 0) - 1) % prime
        if value:
            terms[zero] = value
        else:
            terms.pop(zero, None)
    serialized: list[str] = []
    for exponent in sorted(terms, reverse=True):
        coefficient = terms[exponent]
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        for name, power in zip(MU2_ELIMINATION_VARIABLES, exponent):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        serialized.append("*".join(factors) if factors else "1")
    return "+".join(serialized) or "0"


def active_parameter_polynomial_string(
    polynomial: ParameterPolynomial,
    prime: int,
    variable_names: tuple[str, ...],
    active_parameters: tuple[int, ...],
    multiply_by_r: bool = False,
    subtract_one: bool = False,
) -> str:
    """Serialize a polynomial after several variables were eliminated."""

    terms: dict[tuple[int, ...], int] = {}
    inactive = set(range(len(PARAMETER_NAMES))) - set(active_parameters)
    for exponent, coefficient in polynomial.items():
        if any(exponent[index] for index in inactive):
            raise ValueError("an eliminated parameter is still present")
        ring_exponent = [0] * len(variable_names)
        ring_exponent[0] = int(multiply_by_r)
        for target, source in enumerate(active_parameters, start=1):
            ring_exponent[target] = exponent[source]
        value = (
            terms.get(tuple(ring_exponent), 0) + coefficient
        ) % prime
        if value:
            terms[tuple(ring_exponent)] = value
        else:
            terms.pop(tuple(ring_exponent), None)
    if subtract_one:
        zero = (0,) * len(variable_names)
        value = (terms.get(zero, 0) - 1) % prime
        if value:
            terms[zero] = value
        else:
            terms.pop(zero, None)
    serialized: list[str] = []
    for exponent in sorted(terms, reverse=True):
        coefficient = terms[exponent]
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        for name, power in zip(variable_names, exponent):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        serialized.append("*".join(factors) if factors else "1")
    return "+".join(serialized) or "0"


def first_mu2_pivot_boundary_substitution(
    prime: int,
) -> ParameterPolynomial:
    """Return b1=-3*b0*(lambda+1)/4 on the first mu2 pivot boundary."""

    inverse_four = pow(4, -1, prime)
    b0 = monomial(4, 1, prime)
    lambda_plus_one = sum_parameter_polynomials(
        prime,
        (1, monomial(0, 1, prime)),
        (1, constant_parameter_polynomial(1, prime)),
    )
    return sum_parameter_polynomials(
        prime,
        (
            -3 * inverse_four,
            multiply_parameter_polynomials(b0, lambda_plus_one, prime),
        ),
    )


def mu2_pivot_boundary_reduced_msolve_source(
    moments: dict[int, ParameterPolynomial],
    prime: int,
    minor_name: str,
    branch: str,
    f4sat: bool,
) -> tuple[str, dict[str, object]]:
    """Split and reduce the exceptional a3-pivot divisor of mu2."""

    if minor_name != "01":
        raise ValueError("the reduced pivot-boundary modes use minor 01")
    b1_replacement = first_mu2_pivot_boundary_substitution(prime)
    boundary_moments = {
        order: substitute_parameter_variable(
            polynomial, 5, b1_replacement, prime
        )
        for order, polynomial in moments.items()
    }
    vanished_a3_pivot, _ = split_linear_parameter_polynomial(
        boundary_moments[2], 3, prime
    )
    if vanished_a3_pivot:
        raise ValueError(
            "the first mu2 pivot did not vanish after its parametrization"
        )
    a2_pivot, a2_rest = split_linear_parameter_polynomial(
        boundary_moments[2], 2, prime
    )
    discriminant = substitute_parameter_variable(
        discriminant_polynomial(prime), 5, b1_replacement, prime
    )
    minor = substitute_parameter_variable(
        channel_minor_polynomial(minor_name, prime),
        5,
        b1_replacement,
        prime,
    )

    if branch != "secondary":
        reduced_moments = {
            order: linear_substitution_numerator(
                polynomial, 2, a2_pivot, a2_rest, prime
            )
            for order, polynomial in boundary_moments.items()
            if order > 2
        }
        reduced_discriminant = linear_substitution_numerator(
            discriminant, 2, a2_pivot, a2_rest, prime
        )
        reduced_minor = linear_substitution_numerator(
            minor, 2, a2_pivot, a2_rest, prime
        )
        open_product = multiply_parameter_polynomials(
            a2_pivot, reduced_discriminant, prime
        )
        open_product = multiply_parameter_polynomials(
            open_product, reduced_minor, prime
        )
        if branch == "generic":
            variables = MU2_BOUNDARY_REDUCED_VARIABLES
            active = (1, 3, 4, 6, 0)
            generators = [
                active_parameter_polynomial_string(
                    open_product,
                    prime,
                    variables,
                    active,
                    multiply_by_r=True,
                    subtract_one=True,
                )
            ]
            generators.extend(
                active_parameter_polynomial_string(
                    polynomial, prime, variables, active
                )
                for polynomial in reduced_moments.values()
            )
            metadata: dict[str, object] = {
                "branch": "first_pivot_zero_second_pivot_open",
                "second_pivot_terms": len(a2_pivot),
                "later_moment_terms": {
                    str(order): len(polynomial)
                    for order, polynomial in reduced_moments.items()
                },
            }
        else:
            third_pivot, third_rest = split_linear_parameter_polynomial(
                reduced_moments[3], 3, prime
            )
            third_pivot_reduced = divide_parameter_polynomials(
                third_pivot, a2_pivot, prime
            )
            variables = MU2_BOUNDARY_REDUCED_VARIABLES
            active = (1, 3, 4, 6, 0)
            if branch == "tertiary":
                tertiary_moments: dict[int, ParameterPolynomial] = {}
                cancelled_factors: dict[str, dict[str, int]] = {}
                for order, polynomial in reduced_moments.items():
                    if order <= 3:
                        continue
                    reduced = linear_substitution_numerator(
                        polynomial,
                        3,
                        third_pivot,
                        third_rest,
                        prime,
                    )
                    reduced, second_exponent = strip_parameter_factor(
                        reduced, a2_pivot, prime
                    )
                    reduced, third_exponent = strip_parameter_factor(
                        reduced, third_pivot_reduced, prime
                    )
                    tertiary_moments[order] = reduced
                    cancelled_factors[str(order)] = {
                        "second_pivot": second_exponent,
                        "third_pivot": third_exponent,
                    }
                tertiary_open = multiply_parameter_polynomials(
                    open_product, third_pivot_reduced, prime
                )
                localized_open = tertiary_open
                variables = MU2_BOUNDARY_TERTIARY_VARIABLES
                active = (1, 4, 6, 0)
                generators = [
                    active_parameter_polynomial_string(
                        tertiary_open,
                        prime,
                        variables,
                        active,
                        multiply_by_r=True,
                        subtract_one=True,
                    )
                ]
                generators.extend(
                    active_parameter_polynomial_string(
                        polynomial, prime, variables, active
                    )
                    for polynomial in tertiary_moments.values()
                )
                metadata = {
                    "branch": (
                        "first_pivot_zero_second_and_third_pivots_open"
                    ),
                    "third_pivot_terms": len(third_pivot_reduced),
                    "cancelled_open_factor_powers": cancelled_factors,
                    "later_moment_terms": {
                        str(order): len(polynomial)
                        for order, polynomial in tertiary_moments.items()
                    },
                }
            elif branch == "tertiary-boundary":
                localized_open = open_product
                generators = [
                    active_parameter_polynomial_string(
                        open_product,
                        prime,
                        variables,
                        active,
                        multiply_by_r=True,
                        subtract_one=True,
                    ),
                    active_parameter_polynomial_string(
                        third_pivot_reduced, prime, variables, active
                    ),
                    active_parameter_polynomial_string(
                        third_rest, prime, variables, active
                    ),
                ]
                generators.extend(
                    active_parameter_polynomial_string(
                        polynomial, prime, variables, active
                    )
                    for order, polynomial in reduced_moments.items()
                    if order > 3
                )
                metadata = {
                    "branch": (
                        "first_pivot_zero_second_open_third_pivot_zero"
                    ),
                    "third_pivot_terms": len(third_pivot_reduced),
                    "third_rest_terms": len(third_rest),
                    "later_moment_terms": {
                        str(order): len(polynomial)
                        for order, polynomial in reduced_moments.items()
                        if order > 3
                    },
                }
            else:
                raise ValueError(f"unknown reduced boundary branch {branch}")
    else:
        b2_pivot, b2_rest = split_linear_parameter_polynomial(
            a2_pivot, 6, prime
        )
        if set(b2_pivot) != {zero_exponent()}:
            raise ValueError("the secondary b2 pivot is not constant")
        inverse_b2_pivot = pow(
            b2_pivot[zero_exponent()], -1, prime
        )
        b2_replacement = sum_parameter_polynomials(
            prime, (-inverse_b2_pivot, b2_rest)
        )
        reduced_moments = {
            order: substitute_parameter_variable(
                polynomial, 6, b2_replacement, prime
            )
            for order, polynomial in boundary_moments.items()
        }
        if substitute_parameter_variable(
            a2_pivot, 6, b2_replacement, prime
        ):
            raise ValueError(
                "the second mu2 pivot did not vanish after "
                "its parametrization"
            )
        reduced_discriminant = substitute_parameter_variable(
            discriminant, 6, b2_replacement, prime
        )
        reduced_minor = substitute_parameter_variable(
            minor, 6, b2_replacement, prime
        )
        open_product = multiply_parameter_polynomials(
            reduced_discriminant, reduced_minor, prime
        )
        localized_open = open_product
        variables = MU2_BOUNDARY_SECONDARY_VARIABLES
        active = (1, 2, 3, 4, 0)
        generators = [
            active_parameter_polynomial_string(
                open_product,
                prime,
                variables,
                active,
                multiply_by_r=True,
                subtract_one=True,
            )
        ]
        generators.extend(
            active_parameter_polynomial_string(
                polynomial, prime, variables, active
            )
            for polynomial in reduced_moments.values()
        )
        metadata = {
            "branch": "first_and_second_pivots_zero",
            "b2_pivot": b2_pivot[zero_exponent()],
            "moment_terms": {
                str(order): len(polynomial)
                for order, polynomial in reduced_moments.items()
            },
        }

    if branch in ("generic",):
        localized_open = open_product
    if f4sat:
        open_generator = active_parameter_polynomial_string(
            localized_open, prime, variables, active
        )
        generators = [*generators[1:], open_generator]
        source_variables = variables[1:]
    else:
        source_variables = variables
    source = (
        ",".join(source_variables)
        + f"\n{prime}\n"
        + ",\n".join(generators)
        + "\n"
    )
    metadata["f4sat"] = f4sat
    return source, metadata


def rational_polynomial(
    polynomial: ParameterPolynomial,
) -> RationalParameterPolynomial:
    return {
        exponent: Fraction(coefficient)
        for exponent, coefficient in polynomial.items()
    }


def add_rational_polynomial(
    target: RationalParameterPolynomial,
    source: RationalParameterPolynomial,
    scale: Fraction = Fraction(1),
) -> None:
    for exponent, coefficient in source.items():
        value = target.get(exponent, Fraction(0)) + scale * coefficient
        if value:
            target[exponent] = value
        else:
            target.pop(exponent, None)


def multiply_rational_polynomials(
    left: RationalParameterPolynomial,
    right: RationalParameterPolynomial,
) -> RationalParameterPolynomial:
    answer: RationalParameterPolynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(len(PARAMETER_NAMES))
            )
            value = (
                answer.get(exponent, Fraction(0))
                + left_coefficient * right_coefficient
            )
            if value:
                answer[exponent] = value
            else:
                answer.pop(exponent, None)
    return answer


def rational_polynomial_power(
    polynomial: RationalParameterPolynomial,
    exponent: int,
) -> RationalParameterPolynomial:
    answer: RationalParameterPolynomial = {
        zero_exponent(): Fraction(1)
    }
    for _ in range(exponent):
        answer = multiply_rational_polynomials(answer, polynomial)
    return answer


def substitute_rational_variable(
    polynomial: RationalParameterPolynomial,
    variable: int,
    replacement: RationalParameterPolynomial,
) -> RationalParameterPolynomial:
    maximum_power = max(
        (exponent[variable] for exponent in polynomial), default=0
    )
    replacement_powers = [
        rational_polynomial_power(replacement, power)
        for power in range(maximum_power + 1)
    ]
    answer: RationalParameterPolynomial = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[variable]
        residual_exponent = list(exponent)
        residual_exponent[variable] = 0
        summand = multiply_rational_polynomials(
            {tuple(residual_exponent): coefficient},
            replacement_powers[power],
        )
        add_rational_polynomial(answer, summand)
    return answer


def split_linear_rational_polynomial(
    polynomial: RationalParameterPolynomial,
    variable: int,
) -> tuple[RationalParameterPolynomial, RationalParameterPolynomial]:
    pivot: RationalParameterPolynomial = {}
    rest: RationalParameterPolynomial = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[variable]
        if power not in (0, 1):
            raise ValueError("requested rational polynomial is not linear")
        reduced_exponent = list(exponent)
        reduced_exponent[variable] = 0
        add_rational_polynomial(
            pivot if power else rest,
            {tuple(reduced_exponent): coefficient},
        )
    return pivot, rest


def rational_linear_substitution_numerator(
    polynomial: RationalParameterPolynomial,
    variable: int,
    pivot: RationalParameterPolynomial,
    rest: RationalParameterPolynomial,
) -> RationalParameterPolynomial:
    groups: dict[int, RationalParameterPolynomial] = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[variable]
        reduced_exponent = list(exponent)
        reduced_exponent[variable] = 0
        add_rational_polynomial(
            groups.setdefault(power, {}),
            {tuple(reduced_exponent): coefficient},
        )
    clearing_power = max(groups, default=0)
    pivot_powers = [
        rational_polynomial_power(pivot, power)
        for power in range(clearing_power + 1)
    ]
    rest_powers = [
        rational_polynomial_power(rest, power)
        for power in range(clearing_power + 1)
    ]
    answer: RationalParameterPolynomial = {}
    for power, group in groups.items():
        multiplier = multiply_rational_polynomials(
            rest_powers[power],
            pivot_powers[clearing_power - power],
        )
        summand = multiply_rational_polynomials(group, multiplier)
        add_rational_polynomial(
            answer, summand, Fraction(-1 if power % 2 else 1)
        )
    return answer


def divide_rational_polynomials(
    dividend: RationalParameterPolynomial,
    divisor: RationalParameterPolynomial,
) -> RationalParameterPolynomial:
    """Return the exact sparse quotient for a known rational factor."""

    remainder = dict(dividend)
    quotient: RationalParameterPolynomial = {}
    divisor_lead = max(divisor)
    divisor_coefficient = divisor[divisor_lead]
    while remainder:
        remainder_lead = max(remainder)
        if any(
            left < right
            for left, right in zip(remainder_lead, divisor_lead)
        ):
            raise ValueError("rational polynomial division has a remainder")
        exponent = tuple(
            left - right
            for left, right in zip(remainder_lead, divisor_lead)
        )
        coefficient = (
            remainder[remainder_lead] / divisor_coefficient
        )
        quotient[exponent] = (
            quotient.get(exponent, Fraction(0)) + coefficient
        )
        for divisor_exponent, divisor_value in divisor.items():
            target = tuple(
                left + right
                for left, right in zip(exponent, divisor_exponent)
            )
            value = (
                remainder.get(target, Fraction(0))
                - coefficient * divisor_value
            )
            if value:
                remainder[target] = value
            else:
                remainder.pop(target, None)
    return quotient


def strip_rational_factor(
    polynomial: RationalParameterPolynomial,
    factor: RationalParameterPolynomial,
) -> tuple[RationalParameterPolynomial, int]:
    """Cancel every exact power of a known-nonzero rational factor."""

    answer = polynomial
    exponent = 0
    while answer:
        try:
            quotient = divide_rational_polynomials(answer, factor)
        except ValueError:
            break
        answer = quotient
        exponent += 1
    return answer, exponent


def primitive_integer_coefficients(
    polynomial: RationalParameterPolynomial,
) -> dict[Exponent, int]:
    denominator = 1
    for coefficient in polynomial.values():
        denominator = lcm(denominator, coefficient.denominator)
    integer_polynomial = {
        exponent: int(coefficient * denominator)
        for exponent, coefficient in polynomial.items()
    }
    content = 0
    for coefficient in integer_polynomial.values():
        content = gcd(content, abs(coefficient))
    if content > 1:
        integer_polynomial = {
            exponent: coefficient // content
            for exponent, coefficient in integer_polynomial.items()
        }
    return integer_polynomial


def active_rational_polynomial_string(
    polynomial: RationalParameterPolynomial,
    variable_names: tuple[str, ...],
    active_parameters: tuple[int, ...],
    multiply_by_r: bool = False,
    subtract_one: bool = False,
) -> str:
    primitive = primitive_integer_coefficients(polynomial)
    terms: dict[tuple[int, ...], int] = {}
    inactive = set(range(len(PARAMETER_NAMES))) - set(active_parameters)
    for exponent, coefficient in primitive.items():
        if any(exponent[index] for index in inactive):
            raise ValueError("an eliminated parameter is still present")
        ring_exponent = [0] * len(variable_names)
        ring_exponent[0] = int(multiply_by_r)
        for target, source in enumerate(active_parameters, start=1):
            ring_exponent[target] = exponent[source]
        terms[tuple(ring_exponent)] = coefficient
    if subtract_one:
        zero = (0,) * len(variable_names)
        terms[zero] = terms.get(zero, 0) - 1
        if not terms[zero]:
            terms.pop(zero, None)
    serialized: list[str] = []
    for exponent in sorted(terms, reverse=True):
        coefficient = terms[exponent]
        factors: list[str] = []
        if coefficient != 1 or all(value == 0 for value in exponent):
            factors.append(str(coefficient))
        for name, power in zip(variable_names, exponent):
            factor = format_power(name, power)
            if factor:
                factors.append(factor)
        serialized.append("*".join(factors) if factors else "1")
    return "+".join(serialized).replace("+-", "-") or "0"


def exact_mu2_pivot_boundary_reduced_msolve_source(
    moments: dict[int, ParameterPolynomial],
    minor_name: str,
    branch: str,
    f4sat: bool,
) -> tuple[str, dict[str, object]]:
    """Characteristic-zero version of the two reduced boundary charts."""

    if minor_name != "01":
        raise ValueError("the exact reduced boundary modes use minor 01")
    b1_replacement: RationalParameterPolynomial = {
        (1, 0, 0, 0, 1, 0, 0): Fraction(-3, 4),
        (0, 0, 0, 0, 1, 0, 0): Fraction(-3, 4),
    }
    boundary_moments = {
        order: substitute_rational_variable(
            rational_polynomial(polynomial), 5, b1_replacement
        )
        for order, polynomial in moments.items()
    }
    vanished_a3_pivot, _ = split_linear_rational_polynomial(
        boundary_moments[2], 3
    )
    if vanished_a3_pivot:
        raise ValueError(
            "the exact first mu2 pivot did not vanish after "
            "its parametrization"
        )
    a2_pivot, a2_rest = split_linear_rational_polynomial(
        boundary_moments[2], 2
    )
    exact_discriminant = exact_parameter_polynomial(
        discriminant_polynomial(EXACT_LIFT_PRIME)
    )
    exact_minor = exact_parameter_polynomial(
        channel_minor_polynomial(minor_name, EXACT_LIFT_PRIME)
    )
    discriminant = substitute_rational_variable(
        rational_polynomial(exact_discriminant), 5, b1_replacement
    )
    minor = substitute_rational_variable(
        rational_polynomial(exact_minor), 5, b1_replacement
    )

    if branch != "secondary":
        reduced_moments = {
            order: rational_linear_substitution_numerator(
                polynomial, 2, a2_pivot, a2_rest
            )
            for order, polynomial in boundary_moments.items()
            if order > 2
        }
        reduced_discriminant = rational_linear_substitution_numerator(
            discriminant, 2, a2_pivot, a2_rest
        )
        reduced_minor = rational_linear_substitution_numerator(
            minor, 2, a2_pivot, a2_rest
        )
        open_product = multiply_rational_polynomials(
            a2_pivot, reduced_discriminant
        )
        open_product = multiply_rational_polynomials(
            open_product, reduced_minor
        )
        if branch == "generic":
            variables = MU2_BOUNDARY_REDUCED_VARIABLES
            active = (1, 3, 4, 6, 0)
            generators = [
                active_rational_polynomial_string(
                    open_product,
                    variables,
                    active,
                    multiply_by_r=True,
                    subtract_one=True,
                )
            ]
            generators.extend(
                active_rational_polynomial_string(
                    polynomial, variables, active
                )
                for polynomial in reduced_moments.values()
            )
            metadata: dict[str, object] = {
                "branch": "first_pivot_zero_second_pivot_open",
                "second_pivot_terms": len(a2_pivot),
                "later_moment_terms": {
                    str(order): len(polynomial)
                    for order, polynomial in reduced_moments.items()
                },
            }
        else:
            third_pivot, third_rest = (
                split_linear_rational_polynomial(
                    reduced_moments[3], 3
                )
            )
            third_pivot_reduced = divide_rational_polynomials(
                third_pivot, a2_pivot
            )
            variables = MU2_BOUNDARY_REDUCED_VARIABLES
            active = (1, 3, 4, 6, 0)
            if branch == "tertiary":
                tertiary_moments: dict[
                    int, RationalParameterPolynomial
                ] = {}
                cancelled_factors: dict[str, dict[str, int]] = {}
                for order, polynomial in reduced_moments.items():
                    if order <= 3:
                        continue
                    reduced = rational_linear_substitution_numerator(
                        polynomial, 3, third_pivot, third_rest
                    )
                    reduced, second_exponent = strip_rational_factor(
                        reduced, a2_pivot
                    )
                    reduced, third_exponent = strip_rational_factor(
                        reduced, third_pivot_reduced
                    )
                    tertiary_moments[order] = reduced
                    cancelled_factors[str(order)] = {
                        "second_pivot": second_exponent,
                        "third_pivot": third_exponent,
                    }
                tertiary_open = multiply_rational_polynomials(
                    open_product, third_pivot_reduced
                )
                localized_open = tertiary_open
                variables = MU2_BOUNDARY_TERTIARY_VARIABLES
                active = (1, 4, 6, 0)
                generators = [
                    active_rational_polynomial_string(
                        tertiary_open,
                        variables,
                        active,
                        multiply_by_r=True,
                        subtract_one=True,
                    )
                ]
                generators.extend(
                    active_rational_polynomial_string(
                        polynomial, variables, active
                    )
                    for polynomial in tertiary_moments.values()
                )
                metadata = {
                    "branch": (
                        "first_pivot_zero_second_and_third_pivots_open"
                    ),
                    "third_pivot_terms": len(third_pivot_reduced),
                    "cancelled_open_factor_powers": cancelled_factors,
                    "later_moment_terms": {
                        str(order): len(polynomial)
                        for order, polynomial in tertiary_moments.items()
                    },
                }
            elif branch == "tertiary-boundary":
                localized_open = open_product
                generators = [
                    active_rational_polynomial_string(
                        open_product,
                        variables,
                        active,
                        multiply_by_r=True,
                        subtract_one=True,
                    ),
                    active_rational_polynomial_string(
                        third_pivot_reduced, variables, active
                    ),
                    active_rational_polynomial_string(
                        third_rest, variables, active
                    ),
                ]
                generators.extend(
                    active_rational_polynomial_string(
                        polynomial, variables, active
                    )
                    for order, polynomial in reduced_moments.items()
                    if order > 3
                )
                metadata = {
                    "branch": (
                        "first_pivot_zero_second_open_third_pivot_zero"
                    ),
                    "third_pivot_terms": len(third_pivot_reduced),
                    "third_rest_terms": len(third_rest),
                    "later_moment_terms": {
                        str(order): len(polynomial)
                        for order, polynomial in reduced_moments.items()
                        if order > 3
                    },
                }
            else:
                raise ValueError(f"unknown reduced boundary branch {branch}")
    else:
        b2_pivot, b2_rest = split_linear_rational_polynomial(
            a2_pivot, 6
        )
        if set(b2_pivot) != {zero_exponent()}:
            raise ValueError("the exact secondary b2 pivot is not constant")
        b2_replacement = {
            exponent: -coefficient / b2_pivot[zero_exponent()]
            for exponent, coefficient in b2_rest.items()
        }
        expected_b2_replacement: RationalParameterPolynomial = {
            zero_exponent(): Fraction(-1),
            (0, 0, 0, 0, 1, 0, 0): Fraction(9, 16),
            (1, 0, 0, 0, 1, 0, 0): Fraction(1, 8),
            (2, 0, 0, 0, 1, 0, 0): Fraction(9, 16),
        }
        if b2_replacement != expected_b2_replacement:
            raise ValueError("unexpected exact secondary pivot formula")
        reduced_moments = {
            order: substitute_rational_variable(
                polynomial, 6, b2_replacement
            )
            for order, polynomial in boundary_moments.items()
        }
        if substitute_rational_variable(
            a2_pivot, 6, b2_replacement
        ):
            raise ValueError(
                "the exact second mu2 pivot did not vanish after "
                "its parametrization"
            )
        reduced_discriminant = substitute_rational_variable(
            discriminant, 6, b2_replacement
        )
        reduced_minor = substitute_rational_variable(
            minor, 6, b2_replacement
        )
        open_product = multiply_rational_polynomials(
            reduced_discriminant, reduced_minor
        )
        localized_open = open_product
        variables = MU2_BOUNDARY_SECONDARY_VARIABLES
        active = (1, 2, 3, 4, 0)
        generators = [
            active_rational_polynomial_string(
                open_product,
                variables,
                active,
                multiply_by_r=True,
                subtract_one=True,
            )
        ]
        generators.extend(
            active_rational_polynomial_string(
                polynomial, variables, active
            )
            for polynomial in reduced_moments.values()
        )
        metadata = {
            "branch": "first_and_second_pivots_zero",
            "b2_pivot": str(b2_pivot[zero_exponent()]),
            "moment_terms": {
                str(order): len(polynomial)
                for order, polynomial in reduced_moments.items()
            },
        }
    if branch in ("generic",):
        localized_open = open_product
    if f4sat:
        open_generator = active_rational_polynomial_string(
            localized_open, variables, active
        )
        generators = [*generators[1:], open_generator]
        source_variables = variables[1:]
    else:
        source_variables = variables
    source = (
        ",".join(source_variables)
        + "\n0\n"
        + ",\n".join(generators)
        + "\n"
    )
    metadata["f4sat"] = f4sat
    return source, metadata


def mu2_eliminated_msolve_source(
    moments: dict[int, ParameterPolynomial],
    prime: int,
    minor_name: str,
) -> tuple[str, dict[int, FixedPolynomial], FixedPolynomial]:
    eliminated, pivot, _ = eliminate_mu2_a3(moments, prime)
    open_product = multiply_parameter_polynomials(
        discriminant_polynomial(prime),
        channel_minor_polynomial(minor_name, prime),
        prime,
    )
    # The 01 and 02 channel minors and the discriminant do not involve a3.
    # For minor 03, keep the explicit a3 substitution for a later branch.
    if any(exponent[3] for exponent in open_product):
        raise ValueError(
            "--eliminate-mu2 currently supports minors 01 and 02"
        )
    reduced_open = remove_a3(open_product)
    localized_product = multiply_fixed_polynomials(
        reduced_open, pivot, prime
    )
    generators = [
        mu2_eliminated_polynomial_string(
            localized_product,
            prime,
            multiply_by_r=True,
            subtract_one=True,
        )
    ]
    generators.extend(
        mu2_eliminated_polynomial_string(polynomial, prime)
        for polynomial in eliminated.values()
    )
    source = (
        ",".join(MU2_ELIMINATION_VARIABLES)
        + f"\n{prime}\n"
        + ",\n".join(generators)
        + "\n"
    )
    return source, eliminated, pivot


def restore_a3_zero(
    polynomial: FixedPolynomial,
) -> ParameterPolynomial:
    return {
        exponent[:3] + (0,) + exponent[3:]: coefficient
        for exponent, coefficient in polynomial.items()
    }


def mu2_pivot_boundary_msolve_source(
    moments: dict[int, ParameterPolynomial],
    prime: int,
    minor_name: str,
) -> tuple[str, FixedPolynomial, FixedPolynomial]:
    _, pivot, rest = eliminate_mu2_a3(moments, prime)
    open_product = multiply_parameter_polynomials(
        discriminant_polynomial(prime),
        channel_minor_polynomial(minor_name, prime),
        prime,
    )
    generators = [
        combined_ring_polynomial_string(
            open_product,
            prime,
            None,
            multiply_by_r=True,
            subtract_one=True,
        ),
        combined_ring_polynomial_string(
            restore_a3_zero(pivot), prime, None
        ),
        combined_ring_polynomial_string(
            restore_a3_zero(rest), prime, None
        ),
    ]
    generators.extend(
        combined_ring_polynomial_string(
            polynomial, prime, None
        )
        for order, polynomial in moments.items()
        if order > 2
    )
    source = (
        ",".join(COMBINED_GENERIC_VARIABLES)
        + f"\n{prime}\n"
        + ",\n".join(generators)
        + "\n"
    )
    return source, pivot, rest


def normalized_orders(text: str) -> tuple[int, ...]:
    values = tuple(sorted({int(value) for value in text.split(",") if value}))
    if not values or values[0] < 1:
        raise ValueError("orders must be positive")
    return values


def channel_minor(name: str) -> str:
    d3 = "(-1-a1/3-lambda*b2/3)"
    choices = {
        "01": "b1-a1*b0",
        "02": "b2-a2*b0",
        "03": f"{d3}-a3*b0",
    }
    try:
        return choices[name]
    except KeyError as error:
        raise ValueError(f"unknown channel minor {name!r}") from error


def singular_source(
    moments: dict[int, ParameterPolynomial],
    prime: int,
    minor_name: str,
    backend: str,
    lambda_value: int | None,
) -> str:
    if lambda_value is None:
        ring_declaration = (
            f"ring r=({prime},lambda),({','.join(RING_VARIABLES)}),dp;"
        )
        lambda_declaration: list[str] = []
    else:
        ring_declaration = f"ring r={prime},({','.join(RING_VARIABLES)}),dp;"
        lambda_declaration = [f"number lambda={lambda_value % prime};"]
    definitions = [
        ring_declaration,
        *lambda_declaration,
        "poly discriminant=(9+2*a1+lambda*b2)^2"
        "+(3+2*lambda*b1+3*b2)"
        "*((3-lambda)*a1+2*a2-3*lambda-lambda^2*b2);",
        f"poly channelminor={channel_minor(minor_name)};",
    ]
    moment_names: list[str] = []
    for order, polynomial in moments.items():
        name = f"mu{order}"
        moment_names.append(name)
        definitions.append(
            f"poly {name}={singular_polynomial(polynomial, prime)};"
        )
    generators = ["h*discriminant-1", "k*channelminor-1", *moment_names]
    if backend == "slimgb":
        solve_command = "ideal G=slimgb(I);"
    else:
        solve_command = "ideal G=std(I);"
    definitions.extend(
        [
            f"ideal I={','.join(generators)};",
            "int start_timer=timer;",
            solve_command,
            'print("SIC33_HURWITZ_RESULT_BEGIN");',
            'print("elapsed_ticks="+string(timer-start_timer));',
            'print("basis_size="+string(size(G)));',
            'print("dimension="+string(dim(G)));',
            'print("contains_one="+string(reduce(1,G)==0));',
            'print("first_basis_element="+string(G[1]));',
            'print("SIC33_HURWITZ_RESULT_END");',
            "$",
        ]
    )
    return "\n".join(definitions) + "\n"


def limit_address_space(gigabytes: float):
    def apply_limit() -> None:
        byte_limit = int(gigabytes * 1024**3)
        resource.setrlimit(resource.RLIMIT_AS, (byte_limit, byte_limit))

    return apply_limit


@dataclass(frozen=True)
class SingularResult:
    returncode: int | None
    timed_out: bool
    elapsed_seconds: float
    stdout: str
    stderr: str


def run_singular(source: str, timeout: float, memory_gb: float) -> SingularResult:
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sing", prefix="sic33-hurwitz-", delete=False
    ) as handle:
        handle.write(source)
        input_path = Path(handle.name)
    started = time.monotonic()
    process = subprocess.Popen(
        ["Singular", "-q", str(input_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        preexec_fn=limit_address_space(memory_gb),
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout)
        return SingularResult(
            process.returncode,
            False,
            time.monotonic() - started,
            stdout,
            stderr,
        )
    except subprocess.TimeoutExpired:
        process.terminate()
        try:
            stdout, stderr = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
        return SingularResult(
            process.returncode,
            True,
            time.monotonic() - started,
            stdout,
            stderr,
        )
    finally:
        input_path.unlink(missing_ok=True)


def run_msolve(
    source: str,
    timeout: float,
    memory_gb: float,
    eliminate_lambda: bool,
    f4sat: bool = False,
) -> SingularResult:
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".ms", prefix="sic33-hurwitz-", delete=False
    ) as handle:
        handle.write(source)
        input_path = Path(handle.name)
    output_path = input_path.with_suffix(".out")
    started = time.monotonic()
    command = [
            "msolve",
            "-f",
            str(input_path),
            "-o",
            str(output_path),
            "-t",
            "4",
            "-l",
            "2",
            "-v",
            "1",
        ]
    if eliminate_lambda:
        command.extend(["-e", "8"])
    if f4sat:
        command.append("-S")
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        preexec_fn=limit_address_space(memory_gb),
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout)
        solver_output = (
            output_path.read_text() if output_path.exists() else ""
        )
        return SingularResult(
            process.returncode,
            False,
            time.monotonic() - started,
            stdout + solver_output,
            stderr,
        )
    except subprocess.TimeoutExpired:
        process.terminate()
        try:
            stdout, stderr = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
        solver_output = (
            output_path.read_text() if output_path.exists() else ""
        )
        return SingularResult(
            process.returncode,
            True,
            time.monotonic() - started,
            stdout + solver_output,
            stderr,
        )
    finally:
        input_path.unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)


def profile(polynomial: ParameterPolynomial) -> dict[str, int | str]:
    serialized = singular_polynomial(polynomial, 43)
    return {
        "terms": len(polynomial),
        "serialized_characters_at_prime_43": len(serialized),
        "sha256_at_prime_43": sha256(serialized.encode()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", type=int, default=43)
    parser.add_argument("--orders", default="2,3,4,5,6,7,8")
    parser.add_argument(
        "--minor", choices=("01", "02", "03"), default="01"
    )
    parser.add_argument(
        "--backend", choices=("std", "slimgb", "msolve"), default="msolve"
    )
    parser.add_argument("--lambda-value", type=int)
    parser.add_argument(
        "--lambda-sweep",
        action="store_true",
        help=(
            "solve every fixed lambda fibre in GF(p); currently supported "
            "only by the msolve backend"
        ),
    )
    parser.add_argument(
        "--characteristic-zero",
        action="store_true",
        help=(
            "recover the exact integer moments of 3F and solve one fixed "
            "lambda fibre over QQ"
        ),
    )
    parser.add_argument(
        "--combined-saturation",
        action="store_true",
        help=(
            "localize the discriminant and channel minor with one inverse "
            "of their product instead of two separate inverse variables"
        ),
    )
    parser.add_argument(
        "--eliminate-mu2",
        action="store_true",
        help=(
            "on the generic lambda chart, localize the a3 coefficient of "
            "mu2, substitute a3, and clear the minimal denominator powers"
        ),
    )
    parser.add_argument(
        "--mu2-pivot-boundary",
        action="store_true",
        help=(
            "solve the complementary boundary where the a3 coefficient "
            "of mu2 and the remaining part of mu2 both vanish"
        ),
    )
    parser.add_argument(
        "--mu2-pivot-boundary-reduced",
        choices=(
            "generic",
            "secondary",
            "tertiary",
            "tertiary-boundary",
        ),
        help=(
            "parametrize the first mu2 pivot boundary and either invert "
            "the resulting a2 pivot or impose its zero divisor"
        ),
    )
    parser.add_argument(
        "--mu2-generic-reduced",
        choices=(
            "quadratic",
            "quadratic-boundary",
            "linear",
            "linear-boundary",
        ),
        help=(
            "on the first mu2-pivot-open chart, reduce mu3 as a "
            "quadratic in a2 and select one of its four pivot branches"
        ),
    )
    parser.add_argument("--eliminate-lambda", action="store_true")
    parser.add_argument(
        "--f4sat",
        action="store_true",
        help=(
            "use msolve's native F4 saturation instead of an inverse "
            "variable; supported on the reduced mu2 branch modes"
        ),
    )
    parser.add_argument("--timeout", type=float, default=180.0)
    parser.add_argument("--memory-gb", type=float, default=5.0)
    parser.add_argument("--emit", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--no-solve", action="store_true")
    args = parser.parse_args()

    orders = normalized_orders(args.orders)
    if (
        not args.characteristic_zero
        and (args.prime <= 3 or max(orders) * 3 >= args.prime)
    ):
        raise SystemExit("use a prime > max(3, 3*max(orders))")
    if args.characteristic_zero and args.backend != "msolve":
        raise SystemExit("--characteristic-zero requires --backend msolve")
    if (
        args.characteristic_zero
        and args.lambda_value is None
        and not args.mu2_pivot_boundary_reduced
    ):
        raise SystemExit(
            "--characteristic-zero requires --lambda-value unless a "
            "reduced mu2 boundary branch is selected"
        )
    if args.characteristic_zero and args.lambda_sweep:
        raise SystemExit("--characteristic-zero does not support a sweep")
    if args.characteristic_zero and args.eliminate_lambda:
        raise SystemExit(
            "--characteristic-zero requires a fixed lambda fibre"
        )
    if args.characteristic_zero and args.combined_saturation:
        raise SystemExit(
            "--combined-saturation is not yet implemented for "
            "--characteristic-zero"
        )
    if args.eliminate_mu2 and args.characteristic_zero:
        raise SystemExit(
            "--eliminate-mu2 is currently a finite-field research mode"
        )
    if args.eliminate_mu2 and args.backend != "msolve":
        raise SystemExit("--eliminate-mu2 requires --backend msolve")
    if args.eliminate_mu2 and args.lambda_value is not None:
        raise SystemExit("--eliminate-mu2 currently keeps lambda generic")
    if args.eliminate_mu2 and args.lambda_sweep:
        raise SystemExit("--eliminate-mu2 is incompatible with a sweep")
    if args.eliminate_mu2 and 2 not in orders:
        raise SystemExit("--eliminate-mu2 requires moment order 2")
    if args.eliminate_mu2 and args.minor == "03":
        raise SystemExit(
            "--eliminate-mu2 currently supports minors 01 and 02"
        )
    if args.mu2_pivot_boundary and args.eliminate_mu2:
        raise SystemExit(
            "--mu2-pivot-boundary and --eliminate-mu2 are exclusive"
        )
    if args.mu2_pivot_boundary_reduced and (
        args.mu2_pivot_boundary or args.eliminate_mu2
    ):
        raise SystemExit(
            "the three mu2 elimination/boundary modes are exclusive"
        )
    if args.mu2_generic_reduced and (
        args.mu2_pivot_boundary_reduced
        or args.mu2_pivot_boundary
        or args.eliminate_mu2
    ):
        raise SystemExit("the mu2 reduction modes are exclusive")
    if args.mu2_pivot_boundary and args.characteristic_zero:
        raise SystemExit(
            "--mu2-pivot-boundary is currently a finite-field research mode"
        )
    if args.mu2_pivot_boundary and args.backend != "msolve":
        raise SystemExit(
            "--mu2-pivot-boundary requires --backend msolve"
        )
    if args.mu2_pivot_boundary and args.lambda_value is not None:
        raise SystemExit(
            "--mu2-pivot-boundary currently keeps lambda generic"
        )
    if args.mu2_pivot_boundary and args.lambda_sweep:
        raise SystemExit(
            "--mu2-pivot-boundary is incompatible with a sweep"
        )
    if args.mu2_pivot_boundary and 2 not in orders:
        raise SystemExit(
            "--mu2-pivot-boundary requires moment order 2"
        )
    if args.mu2_pivot_boundary_reduced and args.backend != "msolve":
        raise SystemExit(
            "--mu2-pivot-boundary-reduced requires --backend msolve"
        )
    if args.mu2_pivot_boundary_reduced and (
        args.lambda_value is not None or args.lambda_sweep
    ):
        raise SystemExit(
            "--mu2-pivot-boundary-reduced keeps lambda generic"
        )
    if args.mu2_pivot_boundary_reduced and 2 not in orders:
        raise SystemExit(
            "--mu2-pivot-boundary-reduced requires moment order 2"
        )
    if args.mu2_pivot_boundary_reduced and args.minor != "01":
        raise SystemExit(
            "--mu2-pivot-boundary-reduced currently uses minor 01"
        )
    if args.mu2_generic_reduced and args.characteristic_zero:
        raise SystemExit(
            "--mu2-generic-reduced is currently finite-field only"
        )
    if args.mu2_generic_reduced and args.backend != "msolve":
        raise SystemExit("--mu2-generic-reduced requires --backend msolve")
    if args.mu2_generic_reduced and (
        args.lambda_value is not None or args.lambda_sweep
    ):
        raise SystemExit("--mu2-generic-reduced keeps lambda generic")
    if args.mu2_generic_reduced and 2 not in orders:
        raise SystemExit("--mu2-generic-reduced requires moment order 2")
    if args.mu2_generic_reduced and args.minor != "01":
        raise SystemExit("--mu2-generic-reduced currently uses minor 01")
    if args.f4sat and not (
        args.mu2_generic_reduced or args.mu2_pivot_boundary_reduced
    ):
        raise SystemExit("--f4sat requires a reduced mu2 branch mode")
    if args.lambda_sweep and args.lambda_value is not None:
        raise SystemExit("--lambda-sweep and --lambda-value are exclusive")
    if args.lambda_sweep and args.backend != "msolve":
        raise SystemExit("--lambda-sweep currently requires --backend msolve")
    if args.lambda_sweep and args.eliminate_lambda:
        raise SystemExit("--lambda-sweep does not use --eliminate-lambda")

    if args.characteristic_zero:
        moment_polynomials: dict[int, ParameterPolynomial] = {}
        started = time.monotonic()
        recovered = exact_moment_polynomials(orders)
        for order in orders:
            moment_polynomials[order] = recovered[order]
            print(
                f"recovered exact mu_{order}(3F): "
                f"{len(recovered[order])} terms",
                flush=True,
            )
        print(
            "exact coefficient recovery completed in "
            f"{time.monotonic() - started:.3f}s",
            flush=True,
        )
        exact_boundary_metadata: dict[str, object] | None = None
        if args.mu2_pivot_boundary_reduced:
            source, exact_boundary_metadata = (
                exact_mu2_pivot_boundary_reduced_msolve_source(
                    moment_polynomials,
                    args.minor,
                    args.mu2_pivot_boundary_reduced,
                    args.f4sat,
                )
            )
        else:
            assert args.lambda_value is not None
            source = exact_msolve_source(
                moment_polynomials,
                args.minor,
                args.lambda_value,
            )
        if args.emit:
            args.emit.parent.mkdir(parents=True, exist_ok=True)
            args.emit.write_text(source)
        record: dict[str, object] = {
            "status": "exact characteristic-zero fixed-fibre computation",
            "chart": {
                "form": "(1+x)B(y)+x^2(lambda+x)D(y)",
                "scaled_form_used": "3F",
                "B": "1+a1*y+a2*y^2+a3*y^3",
                "D3_elimination": "-1-a1/3-lambda*b2/3",
                "channel_minor_open": args.minor,
                "quadratic_discriminant_open": True,
                "lambda_value": args.lambda_value,
                "mu2_pivot_boundary_reduced": (
                    args.mu2_pivot_boundary_reduced
                ),
            },
            "backend": args.backend,
            "characteristic": 0,
            "orders": list(orders),
            "exact_lift": {
                "prime": str(EXACT_LIFT_PRIME),
                "coefficient_bound": "(3m)!*52^m",
                "prime_exceeds_twice_bound": True,
            },
            "moment_profiles": {
                str(order): {
                    "terms": len(polynomial),
                    "sha256": sha256(
                        exact_polynomial_string(polynomial).encode()
                    ).hexdigest(),
                }
                for order, polynomial in moment_polynomials.items()
            },
            "msolve_input_sha256": sha256(source.encode()).hexdigest(),
            "reproduction_command": " ".join(sys.argv),
        }
        if exact_boundary_metadata is not None:
            record["mu2_elimination"] = exact_boundary_metadata
        if not args.no_solve:
            result = run_msolve(
                source,
                args.timeout,
                args.memory_gb,
                False,
                args.f4sat,
            )
            print(result.stdout, end="")
            if result.stderr:
                print(result.stderr, end="")
            record["solve"] = {
                "returncode": result.returncode,
                "timed_out": result.timed_out,
                "elapsed_seconds": result.elapsed_seconds,
                "memory_limit_gb": args.memory_gb,
                "unit_ideal": (
                    "Grobner basis has a single element" in result.stdout
                    and "No solution" in result.stdout
                    and (
                        "[1]:" in result.stdout
                        or "[-1]:" in result.stdout
                    )
                ),
                "stdout": result.stdout,
                "stderr": result.stderr,
            }
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(record, indent=2, sort_keys=True) + "\n"
            )
        return

    b_base, d_base = base_polynomials(args.prime)
    b_powers = powers(b_base, max(orders), args.prime)
    d_powers = powers(d_base, max(orders), args.prime)
    moment_polynomials: dict[int, ParameterPolynomial] = {}
    for order in orders:
        started = time.monotonic()
        polynomial = moment(
            order, b_powers, d_powers, args.prime
        )
        moment_polynomials[order] = polynomial
        print(
            f"generated mu_{order}: {len(polynomial)} terms "
            f"in {time.monotonic() - started:.3f}s",
            flush=True,
        )

    if args.lambda_sweep:
        fibres: list[dict[str, object]] = []
        for lambda_value in range(args.prime):
            if args.combined_saturation:
                source = combined_msolve_source(
                    moment_polynomials,
                    args.prime,
                    args.minor,
                    lambda_value,
                )
            else:
                source = msolve_source(
                    moment_polynomials,
                    args.prime,
                    args.minor,
                    lambda_value,
                )
            result = run_msolve(
                source,
                args.timeout,
                args.memory_gb,
                False,
            )
            no_solution = (
                "Grobner basis has a single element" in result.stdout
                and "No solution" in result.stdout
                and "[-1]:" in result.stdout
            )
            fibre = {
                "lambda_value": lambda_value,
                "returncode": result.returncode,
                "timed_out": result.timed_out,
                "elapsed_seconds": result.elapsed_seconds,
                "no_solution": no_solution,
                "input_sha256": sha256(source.encode()).hexdigest(),
                "stdout_sha256": sha256(result.stdout.encode()).hexdigest(),
                "stderr": result.stderr,
            }
            fibres.append(fibre)
            print(
                f"lambda={lambda_value}: "
                f"{'unit' if no_solution else 'unresolved'} "
                f"in {result.elapsed_seconds:.3f}s"
                + (" (timeout)" if result.timed_out else ""),
                flush=True,
            )

        record = {
            "status": (
                "bounded exact finite-field fibre sweep; "
                "not a characteristic-zero certificate"
            ),
            "chart": {
                "form": "(1+x)B(y)+x^2(lambda+x)D(y)",
                "B": "1+a1*y+a2*y^2+a3*y^3",
                "D3_elimination": "-1-a1/3-lambda*b2/3",
                "channel_minor_open": args.minor,
                "quadratic_discriminant_open": True,
                "combined_saturation": args.combined_saturation,
            },
            "backend": args.backend,
            "prime": args.prime,
            "orders": list(orders),
            "moment_profiles": {
                str(order): {
                    "terms": len(polynomial),
                    "sha256": sha256(
                        singular_polynomial(polynomial, args.prime).encode()
                    ).hexdigest(),
                }
                for order, polynomial in moment_polynomials.items()
            },
            "lambda_sweep": fibres,
            "all_fibres_no_solution": all(
                fibre["no_solution"] for fibre in fibres
            ),
            "reproduction_command": " ".join(sys.argv),
        }
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(
                json.dumps(record, indent=2, sort_keys=True) + "\n"
            )
        return

    if args.backend == "msolve":
        eliminated_profiles: dict[str, object] | None = None
        if args.eliminate_mu2:
            source, eliminated, mu2_pivot = mu2_eliminated_msolve_source(
                moment_polynomials,
                args.prime,
                args.minor,
            )
            eliminated_profiles = {
                "mu2_pivot_terms": len(mu2_pivot),
                "later_moment_terms": {
                    str(order): len(polynomial)
                    for order, polynomial in eliminated.items()
                },
            }
        elif args.mu2_pivot_boundary:
            source, mu2_pivot, mu2_rest = (
                mu2_pivot_boundary_msolve_source(
                    moment_polynomials,
                    args.prime,
                    args.minor,
                )
            )
            eliminated_profiles = {
                "boundary": "mu2_a3_pivot=mu2_rest=0",
                "mu2_pivot_terms": len(mu2_pivot),
                "mu2_rest_terms": len(mu2_rest),
            }
        elif args.mu2_pivot_boundary_reduced:
            source, eliminated_profiles = (
                mu2_pivot_boundary_reduced_msolve_source(
                    moment_polynomials,
                    args.prime,
                    args.minor,
                    args.mu2_pivot_boundary_reduced,
                    args.f4sat,
                )
            )
        elif args.mu2_generic_reduced:
            source, eliminated_profiles = (
                mu2_generic_reduced_msolve_source(
                    moment_polynomials,
                    args.prime,
                    args.minor,
                    args.mu2_generic_reduced,
                    args.f4sat,
                )
            )
        elif args.combined_saturation:
            source = combined_msolve_source(
                moment_polynomials,
                args.prime,
                args.minor,
                args.lambda_value,
            )
        else:
            source = msolve_source(
                moment_polynomials,
                args.prime,
                args.minor,
                args.lambda_value,
            )
    else:
        source = singular_source(
            moment_polynomials,
            args.prime,
            args.minor,
            args.backend,
            args.lambda_value,
        )
    if args.emit:
        args.emit.parent.mkdir(parents=True, exist_ok=True)
        args.emit.write_text(source)

    record: dict[str, object] = {
        "status": "exploratory",
        "chart": {
            "form": "(1+x)B(y)+x^2(lambda+x)D(y)",
            "B": "1+a1*y+a2*y^2+a3*y^3",
            "D3_elimination": "-1-a1/3-lambda*b2/3",
            "channel_minor_open": args.minor,
            "quadratic_discriminant_open": True,
            "lambda_is_coefficient_parameter": args.lambda_value is None,
            "lambda_value": args.lambda_value,
            "combined_saturation": args.combined_saturation,
            "mu2_eliminated": args.eliminate_mu2,
            "mu2_pivot_boundary": args.mu2_pivot_boundary,
            "mu2_pivot_boundary_reduced": (
                args.mu2_pivot_boundary_reduced
            ),
            "mu2_generic_reduced": args.mu2_generic_reduced,
        },
        "backend": args.backend,
        "prime": args.prime,
        "orders": list(orders),
        "moment_profiles": {
            str(order): {
                "terms": len(polynomial),
                "sha256": sha256(
                    singular_polynomial(polynomial, args.prime).encode()
                ).hexdigest(),
            }
            for order, polynomial in moment_polynomials.items()
        },
        "singular_input_sha256": sha256(source.encode()).hexdigest(),
    }
    if args.backend == "msolve" and eliminated_profiles is not None:
        record["mu2_elimination"] = eliminated_profiles

    if not args.no_solve:
        if args.backend == "msolve":
            result = run_msolve(
                source,
                args.timeout,
                args.memory_gb,
                args.eliminate_lambda,
                args.f4sat,
            )
        else:
            result = run_singular(source, args.timeout, args.memory_gb)
        print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, end="")
        record["solve"] = {
            "returncode": result.returncode,
            "timed_out": result.timed_out,
            "elapsed_seconds": result.elapsed_seconds,
            "memory_limit_gb": args.memory_gb,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
