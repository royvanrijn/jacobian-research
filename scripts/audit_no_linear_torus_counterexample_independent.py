#!/usr/bin/env python3
"""Dependency-free replay of the linear-torus-free Keller certificate.

This script intentionally does not import SymPy or reuse the construction
code.  It implements the needed sparse polynomial arithmetic over Q directly
with Python's standard library.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm


Exponent = tuple[int, int, int]


class Poly:
    """A sparse polynomial in x,y,z with rational coefficients."""

    def __init__(self, terms: dict[Exponent, Fraction | int] | None = None):
        self.terms = {
            exponent: Fraction(coefficient)
            for exponent, coefficient in (terms or {}).items()
            if coefficient
        }

    @staticmethod
    def constant(value: Fraction | int) -> "Poly":
        return Poly({(0, 0, 0): Fraction(value)}) if value else Poly()

    @staticmethod
    def variable(index: int) -> "Poly":
        exponent = [0, 0, 0]
        exponent[index] = 1
        return Poly({tuple(exponent): Fraction(1)})

    def __add__(self, other: "Poly | Fraction | int") -> "Poly":
        other = as_poly(other)
        result = dict(self.terms)
        for exponent, coefficient in other.terms.items():
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
            if not result[exponent]:
                del result[exponent]
        return Poly(result)

    __radd__ = __add__

    def __neg__(self) -> "Poly":
        return Poly(
            {exponent: -coefficient for exponent, coefficient in self.terms.items()}
        )

    def __sub__(self, other: "Poly | Fraction | int") -> "Poly":
        return self + (-as_poly(other))

    def __rsub__(self, other: "Poly | Fraction | int") -> "Poly":
        return as_poly(other) - self

    def __mul__(self, other: "Poly | Fraction | int") -> "Poly":
        other = as_poly(other)
        result: dict[Exponent, Fraction] = {}
        for left_exponent, left_coefficient in self.terms.items():
            for right_exponent, right_coefficient in other.terms.items():
                exponent = tuple(
                    left + right
                    for left, right in zip(left_exponent, right_exponent)
                )
                result[exponent] = (
                    result.get(exponent, Fraction(0))
                    + left_coefficient * right_coefficient
                )
        return Poly(result)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Poly":
        assert exponent >= 0
        result = Poly.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def derivative(self, index: int) -> "Poly":
        result = {}
        for exponent, coefficient in self.terms.items():
            if exponent[index]:
                new_exponent = list(exponent)
                new_exponent[index] -= 1
                result[tuple(new_exponent)] = coefficient * exponent[index]
        return Poly(result)

    def evaluate(self, point: tuple[Fraction, Fraction, Fraction]) -> Fraction:
        return sum(
            coefficient
            * point[0] ** exponent[0]
            * point[1] ** exponent[1]
            * point[2] ** exponent[2]
            for exponent, coefficient in self.terms.items()
        )


def as_poly(value: Poly | Fraction | int) -> Poly:
    return value if isinstance(value, Poly) else Poly.constant(value)


def determinant_3_by_3(matrix: list[list[Poly]]) -> Poly:
    return (
        matrix[0][0] * matrix[1][1] * matrix[2][2]
        + matrix[0][1] * matrix[1][2] * matrix[2][0]
        + matrix[0][2] * matrix[1][0] * matrix[2][1]
        - matrix[0][2] * matrix[1][1] * matrix[2][0]
        - matrix[0][1] * matrix[1][0] * matrix[2][2]
        - matrix[0][0] * matrix[1][2] * matrix[2][1]
    )


def primitive_integer_row(row: list[Fraction]) -> list[int]:
    denominator_lcm = lcm(*(entry.denominator for entry in row))
    integers = [entry.numerator * (denominator_lcm // entry.denominator) for entry in row]
    common_factor = gcd(*(abs(entry) for entry in integers))
    integers = [entry // common_factor for entry in integers]
    first_nonzero = next(entry for entry in integers if entry)
    return [-entry for entry in integers] if first_nonzero < 0 else integers


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant with row pivoting."""

    work = [row[:] for row in matrix]
    size = len(work)
    assert all(len(row) == size for row in work)
    previous_pivot = 1
    sign = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            row
            for row in range(pivot_index, size)
            if work[row][pivot_index]
        )
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                quotient, remainder = divmod(numerator, previous_pivot)
                assert remainder == 0
                work[row][column] = quotient
            work[row][pivot_index] = 0
        previous_pivot = pivot
    return sign * work[-1][-1]


def rational_rank(matrix: list[list[Fraction]]) -> int:
    """Exact row rank over Q."""

    work = [row[:] for row in matrix if any(row)]
    if not work:
        return 0
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [
            entry / pivot_value
            for entry in work[pivot_row]
        ]
        for row in range(pivot_row + 1, row_count):
            factor = work[row][column]
            if factor:
                work[row] = [
                    entry - factor * pivot_entry
                    for entry, pivot_entry in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


x, y, z = (Poly.variable(index) for index in range(3))
variables = (x, y, z)
t = 1 + x * y
q = t**2 * z - y**2 * (1 + 3 * t)
F = (
    -Fraction(1, 2) * t * q,
    y
    - 3 * x * q
    - t * q
    + 2 * t**2 * x**2 * q**4,
    x * (5 - 3 * t)
    + x**3 * z
    - (x * q) ** 4,
)
jacobian = [
    [component.derivative(index) for index in range(3)]
    for component in F
]

# Independent determinant-one and collision checks.
assert determinant_3_by_3(jacobian).terms == {(0, 0, 0): Fraction(1)}
points = (
    (Fraction(0), Fraction(1), Fraction(5)),
    (-Fraction(1), Fraction(2), -Fraction(9)),
    (Fraction(1, 3), -Fraction(4), -Fraction(27)),
    (Fraction(2, 3), -Fraction(1), Fraction(45)),
)
for point in points:
    assert tuple(component.evaluate(point) for component in F) == (
        -Fraction(1, 2),
        Fraction(0),
        Fraction(0),
    )

# Construct the 18 columns of B F - JF A x without symbolic matrix variables.
columns: list[tuple[Poly, Poly, Poly]] = []
for source_row in range(3):
    for source_column in range(3):
        columns.append(
            tuple(
                -jacobian[component][source_row] * variables[source_column]
                for component in range(3)
            )
        )
for target_row in range(3):
    for target_column in range(3):
        columns.append(
            tuple(
                F[target_column] if component == target_row else Poly()
                for component in range(3)
            )
        )
assert len(columns) == 18

coefficient_rows: dict[tuple[int, Exponent], list[Fraction]] = {}
for column, residual in enumerate(columns):
    for component, polynomial in enumerate(residual, start=1):
        for exponent, coefficient in polynomial.terms.items():
            label = (component, exponent)
            coefficient_rows.setdefault(label, [Fraction(0)] * 18)[column] = coefficient
assert len(coefficient_rows) == 734

# Affine-linear vector fields add three source and three target translation
# columns.  Their exact coefficient matrix also has full column rank.
affine_columns = columns[:9]
affine_columns.extend(
    tuple(-jacobian[component][source_row] for component in range(3))
    for source_row in range(3)
)
affine_columns.extend(columns[9:])
affine_columns.extend(
    tuple(
        Poly.constant(1) if component == target_row else Poly()
        for component in range(3)
    )
    for target_row in range(3)
)
assert len(affine_columns) == 24
affine_coefficient_rows: dict[tuple[int, Exponent], list[Fraction]] = {}
for column, residual in enumerate(affine_columns):
    for component, polynomial in enumerate(residual, start=1):
        for exponent, coefficient in polynomial.terms.items():
            label = (component, exponent)
            affine_coefficient_rows.setdefault(
                label, [Fraction(0)] * 24
            )[column] = coefficient
assert len(affine_coefficient_rows) == 785
assert rational_rank(list(affine_coefficient_rows.values())) == 24

certificate_labels = (
    (1, (12, 10, 4)),
    (1, (12, 8, 4)),
    (1, (4, 3, 0)),
    (1, (4, 2, 1)),
    (1, (3, 4, 0)),
    (1, (3, 3, 1)),
    (1, (3, 2, 2)),
    (1, (2, 4, 1)),
    (1, (2, 4, 0)),
    (1, (2, 3, 2)),
    (1, (2, 2, 1)),
    (2, (12, 10, 4)),
    (2, (12, 8, 4)),
    (2, (3, 3, 1)),
    (2, (3, 2, 1)),
    (3, (12, 10, 4)),
    (3, (12, 8, 4)),
    (3, (3, 3, 1)),
)
certificate_matrix = [
    primitive_integer_row(coefficient_rows[label])
    for label in certificate_labels
]
expected_matrix = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [3, 0, 0, 0, 3, 0, 0, 0, 1, -1, -2, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [2, 0, 0, 0, 4, 0, 0, 0, 0, -1, -2, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [2, 0, 0, 0, 2, 0, 0, 0, 1, -1, -2, 0, 0, 0, 0, 0, 0, 0],
    [12, 0, 0, 0, 10, 0, 0, 0, 4, 0, 0, 0, 0, -1, 0, 0, 0, 0],
    [0, 0, 0, 180, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [6, 0, 0, 0, 6, 0, 0, 0, 2, 0, 0, 0, -1, -2, 0, 0, 0, 0],
    [3, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [12, 0, 0, 0, 8, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, -1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0],
]
assert certificate_matrix == expected_matrix
certificate_determinant = bareiss_determinant(certificate_matrix)
assert certificate_determinant == -5

print("PASS independent sparse-Q replay: det(JF) = 1")
print("PASS independent sparse-Q replay: four-point rational collision")
print("PASS independent sparse-Q replay: 734-by-18 coefficient system rebuilt")
print("PASS independent sparse-Q replay: 785-by-24 affine system has full rank")
print(
    "PASS independent Bareiss determinant:",
    certificate_determinant,
)
