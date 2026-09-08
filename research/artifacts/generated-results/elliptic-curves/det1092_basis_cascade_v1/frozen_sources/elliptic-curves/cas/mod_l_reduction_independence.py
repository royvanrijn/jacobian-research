#!/usr/bin/env python3
"""Exact independence certificates modulo an arbitrary small prime ``ell``.

Let ``P_1,...,P_r`` be rational points on a short elliptic curve and let
``ell`` be prime.  If their images in a product of finite quotients

``prod_p E(F_p) / ell E(F_p)``

have column rank ``r`` over ``F_ell``, every integral relation among the
points has all coefficients divisible by ``ell``.  If ``E(Q)[ell]=0``,
infinite descent then proves that the points are independent.

This is the prime-modulus generalization of :mod:`mod2_reduction_independence`.
It deliberately uses exhaustive finite-group enumeration and is intended for
small certificate primes, not cryptographic-size point counting.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from typing import Iterable, Sequence

from mod2_reduction_independence import (
    FinitePoint,
    RationalPoint,
    finite_add,
    finite_curve_points,
    finite_multiply,
    finite_subtract,
)


Q = Fraction


@dataclass(frozen=True)
class ModLReductionSignature:
    """Coordinates of rational points in one finite ``ell``-quotient."""

    modulus: int
    prime: int
    group_order: int
    multiple_subgroup_order: int
    quotient_dimension: int
    rows: tuple[tuple[int, ...], ...]


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def _primes_up_to(bound: int) -> tuple[int, ...]:
    return tuple(value for value in range(2, bound + 1) if _is_prime(value))


def _reduce_rational(value: Fraction, prime: int) -> int:
    value = Q(value)
    if value.denominator % prime == 0:
        raise ValueError(f"denominator is not invertible modulo {prime}")
    return value.numerator * pow(value.denominator, -1, prime) % prime


def _short_curve_data(
    coefficients: Sequence[Fraction],
    points: Sequence[RationalPoint],
    prime: int,
) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    if len(coefficients) != 5 or any(Q(value) for value in coefficients[:3]):
        raise ValueError("the certificate currently requires a short Weierstrass model")
    coefficient_a_q = Q(coefficients[3])
    coefficient_b_q = Q(coefficients[4])
    for point in points:
        x_q, y_q = (Q(value) for value in point)
        if y_q**2 != x_q**3 + coefficient_a_q * x_q + coefficient_b_q:
            raise ValueError("a supplied rational point is not on the short curve")
    coefficient_a = _reduce_rational(coefficient_a_q, prime)
    coefficient_b = _reduce_rational(coefficient_b_q, prime)
    discriminant = -16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
    if discriminant % prime == 0:
        raise ValueError(f"the curve has bad reduction at {prime}")
    reduced_points = tuple(
        (_reduce_rational(Q(point[0]), prime), _reduce_rational(Q(point[1]), prime))
        for point in points
    )
    return coefficient_a, coefficient_b, reduced_points


def mod_l_reduction_signature(
    coefficients: Sequence[Fraction],
    points: Sequence[RationalPoint],
    prime: int,
    modulus: int,
) -> ModLReductionSignature:
    """Return exact coordinates in ``E(F_prime)/modulus*E(F_prime)``."""

    if prime <= 2 or not _is_prime(prime):
        raise ValueError("certificate primes must be odd")
    if not _is_prime(modulus):
        raise ValueError("the quotient modulus must be prime")
    coefficient_a, coefficient_b, reduced_points = _short_curve_data(
        coefficients, points, prime
    )
    finite_points = finite_curve_points(coefficient_a, coefficient_b, prime)
    multiples = {
        finite_multiply(point, modulus, coefficient_a, prime)
        for point in finite_points
    }

    # Build a basis and a coordinate-labelled list of representatives for the
    # elementary abelian group E(F_p)/ell E(F_p).
    basis: list[FinitePoint] = []
    span: list[FinitePoint] = [None]
    coordinates: list[tuple[int, ...]] = [()]
    for point in finite_points:
        if any(
            finite_subtract(point, representative, coefficient_a, prime)
            in multiples
            for representative in span
        ):
            continue
        basis.append(point)
        old_span = tuple(span)
        old_coordinates = tuple(coordinates)
        span = []
        coordinates = []
        for scalar in range(modulus):
            summand = finite_multiply(point, scalar, coefficient_a, prime)
            for representative, coordinate in zip(old_span, old_coordinates):
                span.append(
                    finite_add(representative, summand, coefficient_a, prime)
                )
                coordinates.append((*coordinate, scalar))

    if len(span) * len(multiples) != len(finite_points):
        raise AssertionError("the quotient representatives do not cover the group")
    if len(span) != modulus ** len(basis):
        raise AssertionError("the quotient does not have the expected vector-space size")

    rows = [[0] * len(points) for _ in basis]
    for point_index, reduced in enumerate(reduced_points):
        coordinate = next(
            (
                coordinate
                for representative, coordinate in zip(span, coordinates)
                if finite_subtract(reduced, representative, coefficient_a, prime)
                in multiples
            ),
            None,
        )
        if coordinate is None:
            raise AssertionError("a reduced point missed every quotient coset")
        for basis_index, value in enumerate(coordinate):
            rows[basis_index][point_index] = value

    return ModLReductionSignature(
        modulus=modulus,
        prime=prime,
        group_order=len(finite_points),
        multiple_subgroup_order=len(multiples),
        quotient_dimension=len(basis),
        rows=tuple(tuple(row) for row in rows),
    )


def gf_l_rank(
    rows: Iterable[Sequence[int]], column_count: int, modulus: int
) -> int:
    """Return the row rank over the prime field ``F_modulus``."""

    if not _is_prime(modulus):
        raise ValueError("the matrix modulus must be prime")
    matrix: list[list[int]] = []
    for row in rows:
        if len(row) != column_count:
            raise ValueError("a row has the wrong declared width")
        reduced = [int(value) % modulus for value in row]
        if any(reduced):
            matrix.append(reduced)
    rank = 0
    for column in range(column_count):
        pivot = next(
            (
                index
                for index in range(rank, len(matrix))
                if matrix[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, modulus)
        matrix[rank] = [(value * inverse) % modulus for value in matrix[rank]]
        for index, row in enumerate(matrix):
            if index == rank or row[column] == 0:
                continue
            multiple = row[column]
            matrix[index] = [
                (left - multiple * right) % modulus
                for left, right in zip(row, matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def combined_mod_l_rank(
    signatures: Iterable[ModLReductionSignature],
    column_count: int,
    modulus: int,
) -> int:
    signatures = tuple(signatures)
    if any(signature.modulus != modulus for signature in signatures):
        raise ValueError("all signatures must use the declared modulus")
    return gf_l_rank(
        (row for signature in signatures for row in signature.rows),
        column_count,
        modulus,
    )


def find_mod_l_reduction_certificate(
    coefficients: Sequence[Fraction],
    points: Sequence[RationalPoint],
    *,
    modulus: int,
    prime_bound: int = 1000,
) -> tuple[ModLReductionSignature, ...]:
    """Greedily retain good primes that increase the exact column rank."""

    if not _is_prime(modulus):
        raise ValueError("the quotient modulus must be prime")
    if prime_bound < 3:
        raise ValueError("prime_bound must be at least 3")
    selected: list[ModLReductionSignature] = []
    current_rank = 0
    for prime in _primes_up_to(prime_bound):
        if prime in (2, modulus):
            continue
        try:
            signature = mod_l_reduction_signature(
                coefficients, points, prime, modulus
            )
        except ValueError:
            continue
        candidate_rank = combined_mod_l_rank(
            (*selected, signature), len(points), modulus
        )
        if candidate_rank > current_rank:
            selected.append(signature)
            current_rank = candidate_rank
        if current_rank == len(points):
            break
    return tuple(selected)


def no_rational_l_torsion_reduction_certificate(
    coefficients: Sequence[Fraction], prime: int, modulus: int
) -> bool:
    """Prove ``E(Q)[ell]=0`` when a good finite group has order prime to ell."""

    if prime <= 2 or not _is_prime(prime) or prime == modulus:
        raise ValueError("use an odd certificate prime different from the modulus")
    if not _is_prime(modulus):
        raise ValueError("the torsion modulus must be prime")
    coefficient_a, coefficient_b, _ = _short_curve_data(coefficients, (), prime)
    group_order = len(finite_curve_points(coefficient_a, coefficient_b, prime))
    return group_order % modulus != 0


def find_no_rational_l_torsion_prime(
    coefficients: Sequence[Fraction], *, modulus: int, prime_bound: int = 200
) -> int:
    """Return the least good prime whose group order is not divisible by ell."""

    for prime in _primes_up_to(prime_bound):
        if prime in (2, modulus):
            continue
        try:
            certified = no_rational_l_torsion_reduction_certificate(
                coefficients, prime, modulus
            )
        except ValueError:
            continue
        if certified:
            return prime
    raise ValueError(
        f"no rational-{modulus}-torsion certificate prime was found through "
        f"{prime_bound}"
    )
