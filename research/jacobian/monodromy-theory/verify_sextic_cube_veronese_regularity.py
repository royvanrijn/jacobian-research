"""Exact rational ranks for the algebraic cube edges in Veronese descent."""

from fractions import Fraction
from math import factorial


# (weight of t in the original frame, determinant order, edge exponents,
#  degree of the algebraic coefficient, expected number of nonempty matrices)
EDGE_CASES = (
    (1, 15, 1, 0, 3, 228),
    (1, 15, 1, 1, 3, 192),
    (1, 15, 2, 0, 3, 334),
    (2, 29, 1, 1, 2, 176),
    (2, 29, 3, 0, 2, 384),
)


def binomial(exponent, index):
    if index < 0:
        return Fraction(0)
    numerator = Fraction(1)
    for step in range(index):
        numerator *= exponent - step
    return numerator / factorial(index)


def rank(matrix, columns):
    """Row reduction over the rational numbers."""
    pivot_row = 0
    for col in range(columns):
        pivot = next((row for row in range(pivot_row, len(matrix))
                      if matrix[row][col]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        leading = matrix[pivot_row]
        for row in matrix[pivot_row + 1:]:
            if row[col]:
                multiplier = row[col] / leading[col]
                for j in range(col, columns):
                    row[j] -= multiplier * leading[j]
        pivot_row += 1
    return pivot_row


def verify_case(d, determinant_order, a, b, degree, expected_count):
    epsilon_weight = 3 - b
    count = 0
    for loss in range(d, 6 * d, d):
        if (loss // d) % 3 == 0:
            continue
        maximum_weight = epsilon_weight * loss + a * (12 - loss // d)
        for weight in range(maximum_weight + 1):
            for residue in range(a * degree):
                columns = []
                for i in range(loss, determinant_order):
                    if i % (a * degree) != residue:
                        continue
                    exponent, remainder = divmod(weight - epsilon_weight * i, a)
                    if remainder or (exponent >= 0 and exponent % 3 == 0):
                        continue
                    columns.append(i)
                if not columns:
                    continue
                rows = [
                    n for n in range(determinant_order)
                    if epsilon_weight * n > weight and n % a == residue % a
                ]
                matrix = [
                    [binomial(Fraction(weight - epsilon_weight * i, 3 * a),
                              (n - i) // a)
                     for i in columns]
                    for n in rows
                ]
                assert rank(matrix, len(columns)) == len(columns), (
                    d, determinant_order, a, b, degree, loss, weight, residue,
                )
                count += 1
    assert count == expected_count, (d, determinant_order, a, b, degree, count)
    print(f"{(d, determinant_order, a, b, degree)}: {count} full-column-rank matrices")


if __name__ == "__main__":
    for case in EDGE_CASES:
        verify_case(*case)
