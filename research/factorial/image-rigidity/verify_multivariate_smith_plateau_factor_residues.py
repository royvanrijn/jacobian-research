#!/usr/bin/env python3
"""Replay a monic residue from a stable Smith interval of ordinary parents."""

from collections import Counter
from itertools import product

from verify_multivariate_height_capped_native_unit_interpolation import coefficient_data
from verify_multivariate_joint_native_shape_reconstruction import residue_gcd
from verify_multivariate_shared_residual_block_lifting import divide_mod


def sylvester(first, second):
    first_degree, second_degree = len(first) - 1, len(second) - 1
    return [
        [first[i - k] if 0 <= i - k <= first_degree else 0 for k in range(second_degree)]
        + [second[i - k] if 0 <= i - k <= second_degree else 0 for k in range(first_degree)]
        for i in range(first_degree + second_degree)
    ]


def diagonal_image(matrix, prime, depth):
    size, modulus = len(matrix), prime**depth
    matrix = [[entry % modulus for entry in row] for row in matrix]
    row_change = [[int(i == j) for j in range(size)] for i in range(size)]
    depths = []

    def valuation(entry):
        current = 0
        while current < depth and entry % prime == 0:
            entry //= prime
            current += 1
        return current

    for k in range(size):
        current, i, j = min(
            (valuation(matrix[i][j]), i, j) for i in range(k, size) for j in range(k, size)
        )
        if current == depth:
            depths.extend([depth] * (size - k))
            break
        depths.append(current)
        matrix[k], matrix[i] = matrix[i], matrix[k]
        row_change[k], row_change[i] = row_change[i], row_change[k]
        for row in matrix:
            row[k], row[j] = row[j], row[k]
        divisor = prime**current
        inverse = pow(matrix[k][k] // divisor, -1, modulus // divisor)
        for i in range(k + 1, size):
            ratio = (matrix[i][k] // divisor) * inverse % (modulus // divisor)
            for array in (matrix, row_change):
                array[i] = [
                    (value - ratio * pivot) % modulus
                    for value, pivot in zip(array[i], array[k], strict=True)
                ]
        for j in range(k + 1, size):
            ratio = (matrix[k][j] // divisor) * inverse % (modulus // divisor)
            for row in matrix:
                row[j] = (row[j] - ratio * row[k]) % modulus

    assert all(value == 0 for i, row in enumerate(matrix)
               for j, value in enumerate(row) if i != j)

    def contains(vector):
        return all(
            sum(a * b for a, b in zip(row, vector, strict=True)) % prime**exponent == 0
            for row, exponent in zip(row_change, depths, strict=True)
        )

    return depths, contains


def main():
    rows = Counter({(2, 3): 4, (3, 2): 4})
    prime, input_modulus = 5, 125
    parents = []
    for order in (6, 7):
        data = coefficient_data(rows, order, prime, input_modulus)
        minimum = min(valuation for valuation, _ in data)
        scaled = [prime**(valuation - minimum) * unit % input_modulus
                  for valuation, unit in data]
        inverse = pow(scaled[-1], -1, input_modulus)
        parents.append([coefficient * inverse % input_modulus for coefficient in scaled])
    assert parents == [[1, 16, 90, 0, 90, 16, 1], [1, 0, 101, 85, 85, 101, 0, 1]]
    assert residue_gcd(*parents, prime) == [1, 0, 0, 0, 0, 1]

    matrix = sylvester(*parents)
    depths, contains = diagonal_image(matrix, prime, 3)
    assert [sum(value >= threshold for value in depths) for threshold in (1, 2, 3)] == [5, 2, 2]

    quadratics = [[1 + 5 * a, 2 + 5 * b, 1] for a, b in product(range(5), repeat=2)]
    assert all(not divide_mod(parent, quadratic, 25)[1]
               for parent in parents for quadratic in quadratics)
    saturated = [
        quadratic for quadratic in quadratics
        if contains([prime * coefficient for coefficient in quadratic] + [0] * (len(matrix) - 3))
    ]
    assert saturated == [[1, 12, 1]]
    lifts = [[1 + 25 * a, 12 + 25 * b, 1] for a, b in product(range(5), repeat=2)]
    assert all(any(divide_mod(parent, lift, input_modulus)[1] for parent in parents)
               for lift in lifts)
    print("Verified a unique saturated quadratic from 25 repeated-residue candidates.")


if __name__ == "__main__":
    main()
