#!/usr/bin/env python3
"""Dependency-free audit of the all-order rooted-tree recurrence certificate."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


def fraction(value: str) -> Fraction:
    return Fraction(value)


def matrix_multiply(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [
            sum(left[row][index] * right[index][column] for index in range(3))
            for column in range(3)
        ]
        for row in range(3)
    ]


def matrix_add(
    *matrices: list[list[Fraction]],
) -> list[list[Fraction]]:
    return [
        [sum(matrix[row][column] for matrix in matrices) for column in range(3)]
        for row in range(3)
    ]


def matrix_scale(
    scalar: Fraction, matrix: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [[scalar * entry for entry in row] for row in matrix]


def matrix_vector(
    matrix: list[list[Fraction]], vector: list[Fraction]
) -> list[Fraction]:
    return [
        sum(matrix[row][column] * vector[column] for column in range(3))
        for row in range(3)
    ]


def identity() -> list[list[Fraction]]:
    return [
        [Fraction(int(row == column)) for column in range(3)] for row in range(3)
    ]


def is_zero(matrix: list[list[Fraction]]) -> bool:
    return all(entry == 0 for row in matrix for entry in row)


def third_sequence(
    matrix: list[list[Fraction]], seed: list[Fraction], terms: int
) -> list[Fraction]:
    values = []
    vector = seed
    for _ in range(terms):
        values.append(vector[2])
        vector = matrix_vector(matrix, vector)
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "certificate",
        nargs="?",
        type=Path,
        default=Path(
            "artifacts/generated-results/lr_rooted_tree_normal_classes.json"
        ),
    )
    arguments = parser.parse_args()

    data = json.loads(arguments.certificate.read_text(encoding="utf-8"))
    family = data["all_order_family"]
    matrix = [
        [fraction(entry) for entry in row]
        for row in family["transfer_matrix_at_gamma_0_u_1_over_6"]
    ]
    recurrence = [
        fraction(value)
        for value in family["cayley_hamilton_recurrence_coefficients"]
    ]
    assert all(value > 0 for value in recurrence)

    square = matrix_multiply(matrix, matrix)
    cube = matrix_multiply(square, matrix)
    cayley_hamilton = matrix_add(
        cube,
        matrix_scale(-recurrence[0], square),
        matrix_scale(-recurrence[1], matrix),
        matrix_scale(-recurrence[2], identity()),
    )
    assert is_zero(cayley_hamilton)

    even_seed = [
        fraction(value)
        for value in family["even_seed_vector_at_gamma_0_u_1_over_6"]
    ]
    odd_seed = [
        fraction(value)
        for value in family["odd_seed_vector_at_gamma_0_u_1_over_6"]
    ]
    even_values = third_sequence(matrix, even_seed, 4)
    odd_values = third_sequence(matrix, odd_seed, 3)
    assert [str(value) for value in even_values] == family[
        "even_initial_third_values_k_0_through_3"
    ]
    assert [str(value) for value in odd_values] == family[
        "odd_initial_third_values_k_0_through_2"
    ]
    assert even_values[0] > 0
    assert all(value < 0 for value in even_values[1:])
    assert all(value < 0 for value in odd_values)

    # Cayley-Hamilton gives s_(k+3)=a*s_(k+2)+b*s_(k+1)+c*s_k.
    # Positive a,b,c preserve strict negativity.  The checked seeds therefore
    # prove every odd term and every even term after the exceptional k=0 term.
    print("PASS: exact 3x3 Cayley-Hamilton identity")
    print("PASS: all recurrence coefficients are positive")
    print("PASS: seed signs prove every separator value is nonzero")
    print("PASS: all displayed rooted-tree normal classes are nonzero for n>=2")


if __name__ == "__main__":
    main()
