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


def canonical_rational_unoriented(
    values: Sequence[int | Fraction],
) -> tuple[int, ...]:
    """Return the primitive integer representative of a rational line.

    This is the finite-index-aware counterpart of
    :func:`canonical_unoriented`: denominators from a nonsaturated reference
    subgroup are cleared exactly before content and orientation are removed.
    The cleared denominator should be recorded separately when its saturation
    class is part of the arithmetic metadata.
    """

    rational = tuple(Fraction(value) for value in values)
    if not any(rational):
        raise ValueError("the zero vector has no primitive representative")
    denominator = 1
    for value in rational:
        denominator = lcm(denominator, value.denominator)
    return canonical_unoriented(
        tuple(int(value * denominator) for value in rational)
    )


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


def modular_rank(rows: Sequence[Sequence[int]], prime: int) -> int:
    """Return exact row rank over ``F_prime`` for a prime modulus."""

    prime = int(prime)
    if prime < 2:
        raise ValueError("modular rank needs a prime modulus")
    matrix = [[int(value) % prime for value in row] for row in rows]
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
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(value * inverse) % prime for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                (value - multiplier * pivot_value) % prime
                for value, pivot_value in zip(matrix[index], matrix[rank])
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


def row_basis_coordinates(
    vectors: Sequence[Sequence[int]],
    basis_rows: Sequence[Sequence[int]],
    *,
    require_integral: bool = True,
) -> tuple[tuple[int | Fraction, ...], ...]:
    """Express vectors in an exact rational row basis.

    A full-rank pivot minor is inverted once.  Every answer is replayed in all
    ambient coordinates; when ``require_integral`` is true, the routine also
    certifies membership in the integral lattice generated by the rows.
    """

    if not basis_rows:
        raise ValueError("coordinate basis must be nonempty")
    width = len(basis_rows[0])
    rank = len(basis_rows)
    if any(len(row) != width for row in basis_rows):
        raise ValueError("basis row widths differ")
    if rational_rank(basis_rows) != rank:
        raise ValueError("coordinate basis is not rationally independent")
    matrix = [[Fraction(int(value)) for value in row] for row in basis_rows]
    echelon = [row[:] for row in matrix]
    pivots = []
    current = 0
    for column in range(width):
        pivot = next(
            (index for index in range(current, rank) if echelon[index][column]),
            None,
        )
        if pivot is None:
            continue
        echelon[current], echelon[pivot] = echelon[pivot], echelon[current]
        value = echelon[current][column]
        echelon[current] = [entry / value for entry in echelon[current]]
        for index in range(rank):
            if index == current or not echelon[index][column]:
                continue
            value = echelon[index][column]
            echelon[index] = [
                entry - value * pivot_entry
                for entry, pivot_entry in zip(echelon[index], echelon[current])
            ]
        pivots.append(column)
        current += 1
        if current == rank:
            break
    if len(pivots) != rank:
        raise ArithmeticError("failed to locate a full-rank basis minor")

    # Solve transpose(B[:,pivots]) * c^T = v[pivots]^T by one augmented RREF.
    system = [
        [matrix[column][pivot] for column in range(rank)]
        + [Fraction(1 if row == column else 0) for column in range(rank)]
        for row, pivot in enumerate(pivots)
    ]
    for column in range(rank):
        pivot = next(index for index in range(column, rank) if system[index][column])
        system[column], system[pivot] = system[pivot], system[column]
        value = system[column][column]
        system[column] = [entry / value for entry in system[column]]
        for index in range(rank):
            if index == column or not system[index][column]:
                continue
            value = system[index][column]
            system[index] = [
                entry - value * pivot_entry
                for entry, pivot_entry in zip(system[index], system[column])
            ]
    inverse = [row[rank:] for row in system]
    answers = []
    for vector in vectors:
        if len(vector) != width:
            raise ValueError("vector width differs from coordinate basis")
        restricted = [Fraction(int(vector[pivot])) for pivot in pivots]
        coordinates = tuple(
            sum(inverse[row][column] * restricted[column] for column in range(rank))
            for row in range(rank)
        )
        if any(
            sum(coordinates[row] * matrix[row][column] for row in range(rank))
            != int(vector[column])
            for column in range(width)
        ):
            raise ValueError("vector is outside the rational row span")
        if require_integral:
            if any(value.denominator != 1 for value in coordinates):
                raise ValueError("vector is not in the integral row lattice")
            answers.append(tuple(int(value) for value in coordinates))
        else:
            answers.append(coordinates)
    return tuple(answers)


def add(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    if len(left) != len(right):
        raise ValueError("vector widths differ")
    return tuple(int(a) + int(b) for a, b in zip(left, right))


def subtract(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    if len(left) != len(right):
        raise ValueError("vector widths differ")
    return tuple(int(a) - int(b) for a, b in zip(left, right))
