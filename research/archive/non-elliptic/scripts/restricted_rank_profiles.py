#!/usr/bin/env python3
"""Rank and power-rank profiles used by restricted-conjecture searches.

The search layer deliberately works over a fixed good prime.  A modular rank
is a lower bound for the characteristic-zero generic rank, and a nonzero
specialized matrix power proves that the corresponding polynomial matrix
power is nonzero.  Vanishing at one specialization proves neither generic
rank nor polynomial nilpotence, so callers must label these values as search
diagnostics until a shortlisted map receives an exact QQ certificate.

The functions in this module contain no BCW-specific logic.  They are shared
by circuit searches for cubic Jacobian rank/index and quartic cotangent
Hessian rank/index.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import random
from typing import Iterable, Sequence

import sympy as sp


DEFAULT_PRIME = 1_000_003
COMPLEX_PRIME = 1_000_033  # 1 modulo 4
SQRT_MINUS_ONE = 350_504


@dataclass(frozen=True)
class PowerRankProfile:
    """Ranks of ``M, M^2, ...`` at one deterministic specialization."""

    ranks: tuple[int, ...]
    terminated: bool

    @property
    def rank(self) -> int:
        return self.ranks[0] if self.ranks else 0

    @property
    def sampled_index(self) -> int | None:
        """Return the specialized nilpotency index when a zero power occurs."""

        return len(self.ranks) if self.terminated else None

    @property
    def area(self) -> int:
        """A monotone scalar tie-breaker; the full tuple remains authoritative."""

        return sum(self.ranks)


def residue(value: object, prime: int) -> int:
    """Map an integer/rational SymPy or stdlib value to ``GF(prime)``."""

    rational = sp.Rational(value)
    numerator, denominator = map(int, rational.as_numer_denom())
    return numerator % prime * pow(denominator % prime, -1, prime) % prime


def deterministic_point(dimension: int, prime: int, seed: int) -> tuple[int, ...]:
    generator = random.Random(seed)
    return tuple(generator.randrange(1, prime) for _ in range(dimension))


def modular_rank(matrix: Sequence[Sequence[int]], prime: int) -> int:
    work = [[value % prime for value in row] for row in matrix]
    if not work:
        return 0
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [
            value * inverse % prime for value in work[pivot_row]
        ]
        for row in range(pivot_row + 1, row_count):
            if work[row][column]:
                scale = work[row][column]
                work[row] = [
                    (left - scale * right) % prime
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def modular_product(
    left: Sequence[Sequence[int]],
    right: Sequence[Sequence[int]],
    prime: int,
) -> list[list[int]]:
    if not left:
        return []
    row_count = len(left)
    inner = len(right)
    column_count = len(right[0]) if right else 0
    answer = [[0] * column_count for _ in range(row_count)]
    nonzero_right = [
        [(column, value) for column, value in enumerate(row) if value % prime]
        for row in right
    ]
    for row_index, row in enumerate(left):
        for inner_index, value in enumerate(row):
            value %= prime
            if value:
                for column, right_value in nonzero_right[inner_index]:
                    answer[row_index][column] = (
                        answer[row_index][column] + value * right_value
                    ) % prime
    assert inner == len(left[0])
    return answer


def power_rank_profile(
    matrix: Sequence[Sequence[int]],
    prime: int,
    max_power: int | None = None,
) -> PowerRankProfile:
    """Compute ranks through the first zero power or ``max_power``."""

    if not matrix:
        return PowerRankProfile((0,), True)
    if max_power is None:
        max_power = len(matrix)
    base = [[value % prime for value in row] for row in matrix]
    power = base
    ranks: list[int] = []
    for _ in range(max_power):
        rank = modular_rank(power, prime)
        ranks.append(rank)
        if rank == 0:
            return PowerRankProfile(tuple(ranks), True)
        power = modular_product(power, base, prime)
    return PowerRankProfile(tuple(ranks), False)


def evaluated_jacobian(
    expressions: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    point: Sequence[int],
    prime: int,
    *,
    subtract_identity: bool = False,
) -> list[list[int]]:
    if len(expressions) != len(variables) or len(point) != len(variables):
        raise ValueError("a square map and a point of matching dimension are required")
    substitution = dict(zip(variables, point))
    rows: list[list[int]] = []
    for row_index, expression in enumerate(expressions):
        row: list[int] = []
        for column_index, variable in enumerate(variables):
            value = sp.Poly(
                sp.diff(expression, variable), *variables, domain=sp.QQ
            ).eval(substitution)
            entry = residue(value, prime)
            if subtract_identity and row_index == column_index:
                entry = (entry - 1) % prime
            row.append(entry)
        rows.append(row)
    return rows


def polynomial_jacobian_profile(
    components: Sequence[sp.Poly],
    variables: Sequence[sp.Symbol],
    *,
    prime: int = DEFAULT_PRIME,
    seed: int = 20_260_723,
    max_power: int | None = None,
) -> PowerRankProfile:
    point = deterministic_point(len(variables), prime, seed)
    expressions = [poly.as_expr() for poly in components]
    matrix = evaluated_jacobian(expressions, variables, point, prime)
    return power_rank_profile(matrix, prime, max_power)


def correction_profile(
    expressions: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    prime: int = DEFAULT_PRIME,
    seed: int = 20_260_723,
    max_power: int = 8,
) -> PowerRankProfile:
    """Profile ``J(F-id)`` for a partial, not-yet-homogenized circuit."""

    point = deterministic_point(len(variables), prime, seed)
    matrix = evaluated_jacobian(
        expressions,
        variables,
        point,
        prime,
        subtract_identity=True,
    )
    return power_rank_profile(matrix, prime, max_power)


def _fraction_polynomial(poly: sp.Poly) -> dict[tuple[int, ...], Fraction]:
    return {
        tuple(exponents): Fraction(int(coefficient.p), int(coefficient.q))
        for exponents, coefficient in poly.terms()
        if coefficient
    }


def cotangent_hessian_rank(
    components: Sequence[sp.Poly],
    *,
    seed: int = 20_260_723,
) -> tuple[int, int, int]:
    """Return sampled ranks ``(JH, A, Hess(y.H))`` over a good prime."""

    # Import lazily: the audit module carries the sparse exact block
    # construction, while keeping the lightweight matrix helpers here usable
    # without the Hessian search.
    from audit_fixed_rank_hessian_witness import (  # pylint: disable=import-outside-toplevel
        cotangent_hessian,
        specialization_profile,
    )

    sparse = [_fraction_polynomial(poly) for poly in components]
    jacobian, upper_left, hessian = cotangent_hessian(sparse)
    return specialization_profile(jacobian, upper_left, hessian, seed)


def transformed_hessian_power_profile(
    components: Sequence[sp.Poly],
    *,
    seed: int = 20_260_723,
    max_power: int | None = None,
) -> PowerRankProfile:
    """Sample powers after the orthogonal change to the HN quartic.

    The Witt-block Hessian is transformed by congruence to the Hessian of the
    quartic potential.  Congruence preserves rank but not powers, so this
    routine uses the explicit transformed matrix rather than powers of the
    Witt block.
    """

    from audit_fixed_rank_hessian_witness import (  # pylint: disable=import-outside-toplevel
        cotangent_blocks,
        evaluate_matrix,
        transformed_hessian,
    )

    sparse = [_fraction_polynomial(poly) for poly in components]
    jacobian, upper_left = cotangent_blocks(sparse)
    dimension = len(components)
    point = deterministic_point(2 * dimension, COMPLEX_PRIME, seed)
    jacobian_value = evaluate_matrix(jacobian, list(point), COMPLEX_PRIME)
    upper_left_value = evaluate_matrix(upper_left, list(point), COMPLEX_PRIME)
    matrix = transformed_hessian(
        jacobian_value, upper_left_value, COMPLEX_PRIME
    )
    return power_rank_profile(
        matrix,
        COMPLEX_PRIME,
        max_power=max_power or 2 * dimension,
    )


def profile_key(profile: PowerRankProfile) -> tuple[object, ...]:
    """Lexicographic minimization key retaining the entire power-rank tuple."""

    sentinel = len(profile.ranks) + 1
    sampled_index = profile.sampled_index or sentinel
    return (
        profile.rank,
        sampled_index,
        profile.area,
        profile.ranks,
        not profile.terminated,
    )


def fraction_vector(values: Iterable[object]) -> tuple[Fraction, ...]:
    """Canonical helper used by collision-aware search states."""

    return tuple(Fraction(sp.Rational(value)) for value in values)
