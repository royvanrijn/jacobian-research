"""Small exact integer-lattice helpers with no CAS dependency."""

from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm
from typing import Iterable, Sequence


def content(values: Iterable[int]) -> int:
    """Return the nonnegative gcd of an integer vector."""

    answer = 0
    for value in values:
        answer = gcd(answer, abs(int(value)))
    return answer


def primitive_coordinates(values: Sequence[int]) -> tuple[int, ...]:
    """Return the primitive vector on the same oriented rational ray."""

    coordinates = tuple(int(value) for value in values)
    divisor = content(coordinates)
    if divisor == 0:
        raise ValueError("the zero vector has no primitive representative")
    return tuple(value // divisor for value in coordinates)


def canonical_unoriented(values: Sequence[int]) -> tuple[int, ...]:
    """Return the primitive canonical representative of ``{v,-v}``."""

    coordinates = primitive_coordinates(values)
    first = next(value for value in coordinates if value)
    if first < 0:
        coordinates = tuple(-value for value in coordinates)
    return coordinates


def rational_rank(rows: Sequence[Sequence[int]]) -> int:
    """Compute row rank over Q by exact fraction-free-sized elimination."""

    matrix = [[Fraction(int(value)) for value in row] for row in rows]
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix rows have inconsistent widths")
    rank = 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                value - multiplier * pivot_entry
                for value, pivot_entry in zip(matrix[index], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def rational_nullspace(
    rows: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    """Return primitive integer rows spanning the right nullspace over Q."""

    if not rows:
        raise ValueError("nullspace input must be nonempty")
    matrix = [[Fraction(int(value)) for value in row] for row in rows]
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix rows have inconsistent widths")
    rank = 0
    pivots: list[int] = []
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                value - multiplier * pivot_entry
                for value, pivot_entry in zip(matrix[index], matrix[rank])
            ]
        pivots.append(column)
        rank += 1
        if rank == len(matrix):
            break
    pivot_set = set(pivots)
    free_columns = [column for column in range(width) if column not in pivot_set]
    answer = []
    for free_column in free_columns:
        vector = [Fraction(0) for _ in range(width)]
        vector[free_column] = Fraction(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -matrix[row][free_column]
        denominator = 1
        for value in vector:
            denominator = lcm(denominator, value.denominator)
        integral = [int(value * denominator) for value in vector]
        divisor = content(integral)
        integral = [value // divisor for value in integral]
        first = next(value for value in integral if value)
        if first < 0:
            integral = [-value for value in integral]
        answer.append(tuple(integral))
    return tuple(answer)


def add(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    if len(left) != len(right):
        raise ValueError("vector widths differ")
    return tuple(int(a) + int(b) for a, b in zip(left, right))


def subtract(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    if len(left) != len(right):
        raise ValueError("vector widths differ")
    return tuple(int(a) - int(b) for a, b in zip(left, right))
