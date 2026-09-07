#!/usr/bin/env python3
"""Exact discriminant geometry for Kihara's rank-at-least-14 family.

Kihara's printed quartic has rational-function coefficients in ``t``.  After
multiplication by ``t^4`` its binary-quartic discriminant factors into nine
low-degree, highly repeated factors and one squarefree even factor of degree
398.  The latter is the only high-degree conductor-frontier factor.  This
module derives that statement from the printed construction; no factor
coefficients are transcribed or fitted from specializations.

The degree-398 factor is represented as ``frontier(t) = f(t^2)``, where
``f`` is a primitive degree-199 integer polynomial.  Lightweight helpers
evaluate its homogenization and Hensel-lift simple roots for CRT searches.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gcd
from typing import Iterable, Sequence

import sympy as sp


FACTOR_SIGNATURE = (
    (1, 16),
    (2, 12),
    (2, 12),
    (2, 12),
    (2, 24),
    (2, 40),
    (4, 24),
    (8, 12),
    (8, 26),
    (398, 1),
)


@dataclass(frozen=True)
class DiscriminantGeometry:
    discriminant_degree: int
    factor_signature: tuple[tuple[int, int], ...]
    constant_factor: int
    frontier_coefficients_z: tuple[int, ...]

    @property
    def frontier_degree_z(self) -> int:
        return len(self.frontier_coefficients_z) - 1


def _kihara_symbolic_quartic() -> tuple[sp.Symbol, tuple[sp.Expr, ...]]:
    """Return ascending coefficients of ``t^4 r_t(x)`` over ``QQ[t]``."""

    x, t = sp.symbols("x t")
    field = sp.QQ.frac_field(t)
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
        4 * p**4
        + 8 * p**3 * q
        + 9 * p**2 * q**2
        - 2 * p * q**3
        + 2 * q**4,
    )
    product = sp.Poly(1, x, domain=field)
    for root in tuple(u + value for value in a_values) + tuple(
        -u + value for value in a_values
    ):
        product *= sp.Poly(x - root, x, domain=field)

    approximant = [sp.Integer(0)] * 7
    approximant[6] = sp.Integer(1)
    for index in range(5, -1, -1):
        polynomial = sp.Poly(
            sum(approximant[j] * x**j for j in range(7)), x, domain=field
        )
        degree = 6 + index
        approximant[index] = sp.cancel(
            (product.nth(degree) - (polynomial * polynomial).nth(degree)) / 2
        )
    polynomial = sp.Poly(
        sum(approximant[j] * x**j for j in range(7)), x, domain=field
    )
    remainder = polynomial * polynomial - product
    coefficients = tuple(sp.cancel(remainder.nth(index) * t**4) for index in range(5))
    if any(sp.denom(value) != 1 for value in coefficients):
        raise AssertionError("t^4 did not clear the symbolic quartic denominators")
    return t, coefficients


def _primitive_integer_coefficients(poly: sp.Poly) -> tuple[int, ...]:
    denominator, integer_poly = poly.clear_denoms(convert=True)
    del denominator
    values = [int(integer_poly.nth(index)) for index in range(integer_poly.degree() + 1)]
    content = 0
    for value in values:
        content = gcd(content, abs(value))
    if not content:
        raise ValueError("zero polynomial has no primitive normalization")
    values = [value // content for value in values]
    if values[-1] < 0:
        values = [-value for value in values]
    return tuple(values)


@lru_cache(maxsize=1)
def derive_discriminant_geometry() -> DiscriminantGeometry:
    """Derive and exactly factor the scaled quartic discriminant."""

    t, coefficients = _kihara_symbolic_quartic()
    e, d, c, b, a = coefficients
    invariant_i = sp.Poly(sp.expand(12 * a * e - 3 * b * d + c**2), t, domain=sp.QQ)
    invariant_j = sp.Poly(
        sp.expand(
            72 * a * c * e
            + 9 * b * c * d
            - 27 * a * d**2
            - 27 * b**2 * e
            - 2 * c**3
        ),
        t,
        domain=sp.QQ,
    )
    discriminant = 4 * invariant_i**3 - invariant_j**2
    constant, raw_factors = sp.factor_list(discriminant.as_expr(), t)
    factors = tuple((sp.Poly(factor, t, domain=sp.QQ), exponent) for factor, exponent in raw_factors)
    signature = tuple(
        (int(factor.degree()), int(exponent)) for factor, exponent in factors
    )
    if signature != FACTOR_SIGNATURE:
        raise AssertionError(f"unexpected discriminant factor signature {signature}")
    frontier = next(factor for factor, exponent in factors if factor.degree() == 398 and exponent == 1)
    if any(frontier.nth(index) for index in range(1, 399, 2)):
        raise AssertionError("the degree-398 frontier factor is not even")
    z = sp.symbols("z")
    frontier_z = sp.Poly(
        sum(frontier.nth(2 * index) * z**index for index in range(200)),
        z,
        domain=sp.QQ,
    )
    if frontier_z.degree() != 199 or not frontier_z.is_irreducible:
        raise AssertionError("the frontier factor did not remain irreducible in t^2")
    return DiscriminantGeometry(
        discriminant_degree=discriminant.degree(),
        factor_signature=signature,
        constant_factor=int(constant),
        frontier_coefficients_z=_primitive_integer_coefficients(frontier_z),
    )


def polynomial_value(coefficients: Sequence[int], value: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = answer * value + int(coefficient)
    return answer


def polynomial_derivative_value(coefficients: Sequence[int], value: int) -> int:
    answer = 0
    for degree in range(len(coefficients) - 1, 0, -1):
        answer = answer * value + degree * int(coefficients[degree])
    return answer


def frontier_value_t(coefficients_z: Sequence[int], parameter_t: int) -> int:
    return polynomial_value(coefficients_z, int(parameter_t) ** 2)


def homogeneous_frontier_value(
    coefficients_z: Sequence[int], numerator: int, denominator: int
) -> int:
    """Evaluate ``denominator^398 f((numerator/denominator)^2)`` exactly."""

    if denominator == 0:
        raise ValueError("a rational parameter needs nonzero denominator")
    degree = len(coefficients_z) - 1
    numerator_squared = int(numerator) ** 2
    denominator_squared = int(denominator) ** 2
    answer = 0
    numerator_power = 1
    # Direct homogeneous evaluation avoids constructing a Fraction with
    # thousand-digit numerator and denominator only to cancel it again.
    for index, coefficient in enumerate(coefficients_z):
        answer += int(coefficient) * numerator_power * denominator_squared ** (degree - index)
        numerator_power *= numerator_squared
    return answer


def roots_mod_prime_t(coefficients_z: Sequence[int], prime: int) -> tuple[int, ...]:
    if prime < 2:
        raise ValueError("prime must be at least two")
    return tuple(
        residue
        for residue in range(prime)
        if frontier_value_t(coefficients_z, residue) % prime == 0
    )


def hensel_lift_simple_t_root(
    coefficients_z: Sequence[int], residue: int, prime: int, exponent: int
) -> int:
    """Lift a simple root of ``f(t^2)`` from ``p`` to ``p^exponent``."""

    if prime < 2 or exponent < 1:
        raise ValueError("invalid Hensel modulus")
    residue %= prime
    value = frontier_value_t(coefficients_z, residue)
    derivative = (
        2
        * residue
        * polynomial_derivative_value(coefficients_z, residue * residue)
    )
    if value % prime or derivative % prime == 0:
        raise ValueError("the supplied residue is not a simple root")
    modulus = prime
    lifted = residue
    for _ in range(1, exponent):
        quotient = frontier_value_t(coefficients_z, lifted) // modulus
        correction = -quotient * pow(derivative, -1, prime) % prime
        lifted += correction * modulus
        modulus *= prime
        derivative = (
            2
            * lifted
            * polynomial_derivative_value(coefficients_z, lifted * lifted)
        )
    if frontier_value_t(coefficients_z, lifted) % modulus:
        raise AssertionError("Hensel lifting failed")
    return lifted % modulus


def distinct_prime_divisors_through(value: int, primes: Iterable[int]) -> tuple[int, ...]:
    value = abs(int(value))
    return tuple(prime for prime in primes if value % prime == 0)
