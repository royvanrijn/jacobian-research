#!/usr/bin/env python3
"""Exact independence certificates from reductions modulo good primes.

Let ``P_1,...,P_r`` be rational points on an elliptic curve with no rational
2-torsion.  A rational integral relation reduces at every good prime.  If the
images of the points in the product of finite quotients

``prod_p E(F_p) / 2 E(F_p)``

are linearly independent over ``F_2``, every coefficient in such a relation
is even.  Dividing the relation by two and using ``E(Q)[2]=0`` repeats the
argument indefinitely, so every coefficient is zero.

The code below constructs the finite quotients by exhaustive group
enumeration.  It is intended for small certificate primes, not for point
counting at cryptographic sizes.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from typing import Iterable, Sequence


Q = Fraction
FinitePoint = tuple[int, int] | None
RationalPoint = tuple[Fraction, Fraction]


@dataclass(frozen=True)
class Mod2ReductionSignature:
    prime: int
    group_order: int
    doubled_subgroup_order: int
    quotient_dimension: int
    rows: tuple[tuple[int, ...], ...]


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def _reduce_rational(value: Fraction, prime: int) -> int:
    value = Q(value)
    if value.denominator % prime == 0:
        raise ValueError(f"denominator is not invertible modulo {prime}")
    return value.numerator * pow(value.denominator, -1, prime) % prime


def finite_add(
    left: FinitePoint, right: FinitePoint, coefficient_a: int, prime: int
) -> FinitePoint:
    """Add points on ``y^2=x^3+A*x+B`` over ``F_prime``."""

    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if (y1 + y2) % prime == 0:
            return None
        slope = (3 * x1 * x1 + coefficient_a) * pow(2 * y1, -1, prime)
    else:
        slope = (y2 - y1) * pow((x2 - x1) % prime, -1, prime)
    slope %= prime
    x3 = (slope * slope - x1 - x2) % prime
    y3 = (-y1 + slope * (x1 - x3)) % prime
    return x3, y3


def finite_negate(point: FinitePoint, prime: int) -> FinitePoint:
    if point is None:
        return None
    return point[0], (-point[1]) % prime


def finite_subtract(
    left: FinitePoint, right: FinitePoint, coefficient_a: int, prime: int
) -> FinitePoint:
    return finite_add(left, finite_negate(right, prime), coefficient_a, prime)


def finite_multiply(
    point: FinitePoint, scalar: int, coefficient_a: int, prime: int
) -> FinitePoint:
    if scalar < 0:
        return finite_multiply(finite_negate(point, prime), -scalar, coefficient_a, prime)
    answer: FinitePoint = None
    addend = point
    while scalar:
        if scalar & 1:
            answer = finite_add(answer, addend, coefficient_a, prime)
        addend = finite_add(addend, addend, coefficient_a, prime)
        scalar >>= 1
    return answer


def finite_curve_points(
    coefficient_a: int, coefficient_b: int, prime: int
) -> tuple[FinitePoint, ...]:
    """Enumerate every point on a nonsingular short curve over ``F_prime``."""

    if prime <= 2 or not _is_prime(prime):
        raise ValueError("certificate primes must be odd")
    square_roots: dict[int, list[int]] = {}
    for ordinate in range(prime):
        square_roots.setdefault(ordinate * ordinate % prime, []).append(ordinate)
    points: list[FinitePoint] = [None]
    for abscissa in range(prime):
        rhs = (abscissa**3 + coefficient_a * abscissa + coefficient_b) % prime
        points.extend((abscissa, ordinate) for ordinate in square_roots.get(rhs, ()))
    return tuple(points)


def mod2_reduction_signature(
    coefficients: Sequence[Fraction],
    points: Sequence[RationalPoint],
    prime: int,
) -> Mod2ReductionSignature:
    """Return exact coordinate rows for the images in ``E(F_p)/2E(F_p)``."""

    from research_runtime.finite_reduction import default_reduction_cache
    return default_reduction_cache().signature(coefficients, points, prime)


def gf2_rank(rows: Iterable[Sequence[int]], column_count: int) -> int:
    """Return the row rank of a binary matrix by exact bit elimination."""

    pivots: dict[int, int] = {}
    for row in rows:
        row = tuple(int(value) for value in row)
        if len(row) != column_count or any(value not in (0, 1) for value in row):
            raise ValueError("rows must be binary and have the declared width")
        packed = sum(value << index for index, value in enumerate(row))
        while packed:
            pivot = packed.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = packed
                break
            packed ^= pivots[pivot]
    return len(pivots)


def combined_mod2_rank(
    signatures: Iterable[Mod2ReductionSignature], column_count: int
) -> int:
    return gf2_rank(
        (row for signature in signatures for row in signature.rows), column_count
    )


def _primes_up_to(bound: int) -> tuple[int, ...]:
    primes: list[int] = []
    for candidate in range(2, bound + 1):
        if all(candidate % prime for prime in primes if prime * prime <= candidate):
            primes.append(candidate)
    return tuple(primes)


def find_mod2_reduction_certificate(
    coefficients: Sequence[Fraction],
    points: Sequence[RationalPoint],
    *,
    prime_bound: int = 1000,
) -> tuple[Mod2ReductionSignature, ...]:
    """Greedily retain small good primes that increase the exact binary rank.

    The returned sequence is deterministic.  An empty or rank-deficient
    sequence is a valid negative bounded-search result, not evidence of point
    dependence: the supplied subgroup may simply be nonsaturated at 2.
    """

    if prime_bound < 3:
        raise ValueError("prime_bound must be at least 3")
    selected: list[Mod2ReductionSignature] = []
    current_rank = 0
    for prime in _primes_up_to(prime_bound):
        if prime == 2:
            continue
        try:
            signature = mod2_reduction_signature(coefficients, points, prime)
        except ValueError:
            continue
        candidate_rank = combined_mod2_rank((*selected, signature), len(points))
        if candidate_rank > current_rank:
            selected.append(signature)
            current_rank = candidate_rank
        if current_rank == len(points):
            break
    return tuple(selected)


def short_curve_has_no_rational_2_torsion_modular_certificate(
    coefficients: Sequence[Fraction], prime: int
) -> bool:
    """Certify irreducibility of the 2-division cubic by reduction mod ``prime``.

    A cubic over a field is irreducible exactly when it has no root.  If the
    primitive rational 2-division polynomial remains cubic and irreducible
    modulo one prime, it is irreducible over Q by Gauss's lemma.
    """

    if len(coefficients) != 5 or any(Q(value) for value in coefficients[:3]):
        raise ValueError("the certificate currently requires a short Weierstrass model")
    if prime <= 2 or not _is_prime(prime):
        raise ValueError("the certificate modulus must be an odd prime")
    coefficient_a = _reduce_rational(Q(coefficients[3]), prime)
    coefficient_b = _reduce_rational(Q(coefficients[4]), prime)
    return all(
        (value**3 + coefficient_a * value + coefficient_b) % prime
        for value in range(prime)
    )


def find_two_torsion_certificate_prime(
    coefficients: Sequence[Fraction], *, prime_bound: int = 200
) -> int:
    """Return the least usable prime proving ``E(Q)[2]=0`` by reduction."""

    for prime in _primes_up_to(prime_bound):
        if prime == 2:
            continue
        try:
            certified = short_curve_has_no_rational_2_torsion_modular_certificate(
                coefficients, prime
            )
        except ValueError:
            continue
        if certified:
            return prime
    raise ValueError(
        f"no rational-2-torsion certificate prime was found through {prime_bound}"
    )
