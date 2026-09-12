#!/usr/bin/env python3
"""Dependency-free replay of the rank-35 identity-slice witness."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import random


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_rank_reduced_bcw_22_counterexample.json"
)
PRIME = 1_000_003
SEEDS = (20_260_728, 20_260_729, 20_260_730)

Exponent = tuple[int, ...]
Polynomial = dict[Exponent, Fraction]


def add_term(
    polynomial: Polynomial,
    exponent: Exponent,
    coefficient: Fraction,
) -> None:
    value = polynomial.get(exponent, Fraction(0)) + coefficient
    if value:
        polynomial[exponent] = value
    else:
        polynomial.pop(exponent, None)


def derivative(polynomial: Polynomial, variable: int) -> Polynomial:
    result: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[variable]
        if power:
            reduced = list(exponent)
            reduced[variable] -= 1
            add_term(result, tuple(reduced), coefficient * power)
    return result


def decode_slice(stored: dict[str, object]) -> list[Polynomial]:
    result: list[Polynomial] = []
    for component in stored["H"][:-1]:
        polynomial: Polynomial = {}
        for term in component:
            exponent = [0] * 21
            for variable, power in term["monomial"]:
                if variable < 21:
                    exponent[variable] = power
                else:
                    assert variable == 21
            add_term(
                polynomial,
                tuple(exponent),
                Fraction(term["coefficient"]),
            )
        result.append(polynomial)
    return result


def evaluate_rational(polynomial: Polynomial, point: list[Fraction]) -> Fraction:
    result = Fraction(0)
    for exponent, coefficient in polynomial.items():
        term = coefficient
        for value, power in zip(point, exponent):
            term *= value**power
        result += term
    return result


def residue(value: Fraction) -> int:
    return (
        value.numerator
        * pow(value.denominator % PRIME, -1, PRIME)
        % PRIME
    )


def evaluate_modular(polynomial: Polynomial, point: list[int]) -> int:
    result = 0
    for exponent, coefficient in polynomial.items():
        term = residue(coefficient)
        for value, power in zip(point, exponent):
            term = term * pow(value, power, PRIME) % PRIME
        result = (result + term) % PRIME
    return result


def rank(matrix: list[list[int]]) -> int:
    work = [[value % PRIME for value in row] for row in matrix]
    pivot_row = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next(
            (
                row
                for row in range(pivot_row, len(work))
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, PRIME)
        work[pivot_row] = [
            value * inverse % PRIME for value in work[pivot_row]
        ]
        for row in range(pivot_row + 1, len(work)):
            if work[row][column]:
                factor = work[row][column]
                work[row] = [
                    (left - factor * right) % PRIME
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
    return pivot_row


def main() -> None:
    stored = json.loads(SOURCE.read_text())
    assert stored["dimension"] == 22 and stored["H"][-1] == []
    assert stored["statistics"]["nilpotency_index_JH"] == 18
    components = decode_slice(stored)
    assert len(components) == 21

    points = [
        [Fraction(value) for value in point[:-1]]
        for point in stored["collision_points"]
    ]
    assert [Fraction(point[-1]) for point in stored["collision_points"]] == [
        1,
        1,
        1,
    ]
    images = [
        [
            coordinate + evaluate_rational(component, point)
            for coordinate, component in zip(point, components)
        ]
        for point in points
    ]
    assert images[0] == images[1] == images[2]
    assert [point[0] for point in points] == [0, 1, -1]

    jacobian = [
        [derivative(component, variable) for variable in range(21)]
        for component in components
    ]
    second = [
        [
            [
                derivative(jacobian[output][first], second_variable)
                for output in range(21)
            ]
            for second_variable in range(21)
        ]
        for first in range(21)
    ]

    profiles = []
    for seed in SEEDS:
        generator = random.Random(seed)
        x = [generator.randrange(1, PRIME) for _ in range(21)]
        y = [generator.randrange(1, PRIME) for _ in range(21)]
        jacobian_value = [
            [evaluate_modular(entry, x) for entry in row]
            for row in jacobian
        ]
        upper_left = [
            [
                sum(
                    y[output] * evaluate_modular(second[first][column][output], x)
                    for output in range(21)
                )
                % PRIME
                for column in range(21)
            ]
            for first in range(21)
        ]
        transpose = [list(column) for column in zip(*jacobian_value)]
        hessian = [
            upper_left[row] + transpose[row] for row in range(21)
        ]
        hessian.extend(
            jacobian_value[row] + [0] * 21 for row in range(21)
        )
        profile = (
            rank(jacobian_value),
            rank(upper_left),
            rank(hessian),
        )
        assert profile == (17, 14, 35)
        profiles.append(profile)

    print("PASS independent rank-35 slice: exact rational three-point collision")
    print("PASS independent rank-35 slice: coordinate 0 separates as 0,1,-1")
    print("PASS independent rank-35 slice: three modular block profiles are (17,14,35)")
    print("PASS independent rank-35 slice: no SymPy or Singular dependency")


if __name__ == "__main__":
    main()
