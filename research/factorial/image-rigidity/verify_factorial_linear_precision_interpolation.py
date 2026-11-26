"""Replay the linear denominator bound from frozen partition weights."""

from itertools import product
from math import comb, factorial

from sympy import Matrix

from verify_factorial_moment_prime_power_rationality import (
    finite_transfer,
    native_series,
    valuation,
)


def frozen_transfer(columns, weights, prime):
    dimension = len(columns[0])
    maximum = max(map(max, columns))
    period = prime ** (1 + valuation(factorial(maximum), prime))
    while period <= maximum:
        period *= prime

    states = []
    for coordinates in product(range(-1, period), repeat=dimension):
        residues = tuple(max(0, digit) for digit in coordinates)
        positive_mask = sum(1 << i for i, digit in enumerate(coordinates) if digit >= 0)
        states.append((residues, positive_mask))
    index = {state: i for i, state in enumerate(states)}
    matrix = [[0] * len(states) for _ in states]
    for source, (residues, mask) in enumerate(states):
        for column, weight in zip(columns, weights):
            active = mask | sum(1 << i for i, exponent in enumerate(column) if exponent)
            first = (active & -active).bit_length() - 1
            if not column[first]:
                continue
            frozen = tuple(r + n for r, n in zip(residues, column))
            target = index[tuple(t % period for t in frozen), active]
            factor = comb(frozen[first] - 1, column[first] - 1)
            for i, exponent in enumerate(column):
                if i != first:
                    factor *= comb(frozen[i], exponent)
            matrix[target][source] += weight * factor
    return matrix


def multiply(left, right, modulus):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                result[i + j] += a * b
    return [entry % modulus for entry in result]


def verify():
    cases = (
        (2, (2, 3), ((2,), (1,), (3,)), (1, 2, 1)),
        (2, (2, 3), ((2, 0), (1, 1), (0, 1)), (1, 2, 1)),
        (3, (2, 3), ((3,), (1,), (2,)), (2, 3, 1)),
        (5, (2,), ((1, 0), (0, 1), (1, 1)), (1, 3, 2)),
    )
    checks = 0
    native_checks = 0
    for prime, precisions, columns, source_weights in cases:
        for weights in (source_weights, tuple(weight + prime for weight in source_weights)):
            matrix = frozen_transfer(columns, weights, prime)
            characteristic = [int(coefficient) for coefficient in Matrix(matrix).charpoly().all_coeffs()]
            highest_degree = max(map(sum, columns))
            for precision in precisions:
                modulus = prime**precision
                exponent = 1 + highest_degree * (precision - 1)
                numerator_limit = len(matrix) * exponent - 1
                denominator = [1]
                for _ in range(exponent):
                    denominator = multiply(denominator, characteristic, modulus)
                moments = finite_transfer(columns, weights, prime, precision, numerator_limit + 12)
                assert moments[:9] == native_series(columns, weights, 8, modulus)
                native_checks += 9
                for degree in range(numerator_limit + 1, len(moments)):
                    assert sum(a * moments[degree - j]
                               for j, a in enumerate(denominator[:degree + 1])) % modulus == 0
                    checks += 1
    print(f"{native_checks} literal moments and {checks} denominator tail coefficients verified.")


if __name__ == "__main__":
    verify()
