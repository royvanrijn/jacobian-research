#!/usr/bin/env python3
"""Dependency-free replay of the rank-34 double identity slice."""

from __future__ import annotations

from fractions import Fraction
from math import comb
import json
from pathlib import Path
import random


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_rank_35_identity_slice_counterexample.json"
)
PRIME = 1_000_003
SEEDS = (20_260_728, 20_260_729, 20_260_730)
KEPT = tuple(index for index in range(21) if index != 9)

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


def decode(stored: dict[str, object]) -> list[Polynomial]:
    result = []
    for component in stored["K"]:
        polynomial: Polynomial = {}
        for term in component:
            exponent = [0] * 21
            for variable, power in term["monomial"]:
                exponent[variable] = power
            add_term(
                polynomial,
                tuple(exponent),
                Fraction(term["coefficient"]),
            )
        result.append(polynomial)
    return result


def eliminate(polynomial: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[9]
        reduced = [exponent[index] for index in KEPT]
        for x6_power in range(power + 1):
            expanded = reduced[:]
            expanded[1] += power - x6_power
            expanded[6] += x6_power
            add_term(
                result,
                tuple(expanded),
                coefficient
                * comb(power, x6_power)
                * 3 ** (power - x6_power),
            )
    return result


def derivative(polynomial: Polynomial, variable: int) -> Polynomial:
    result: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        power = exponent[variable]
        if power:
            reduced = list(exponent)
            reduced[variable] -= 1
            add_term(result, tuple(reduced), coefficient * power)
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
    k = decode(stored)
    relation: Polynomial = {}
    for scalar, component in ((1, k[9]), (-3, k[1]), (-1, k[6])):
        for exponent, coefficient in component.items():
            add_term(relation, exponent, scalar * coefficient)
    assert relation == {}

    points_21 = [
        [Fraction(value) for value in point]
        for point in stored["collision_points"]
    ]
    assert [
        -3 * point[1] - point[6] + point[9] for point in points_21
    ] == [0, 0, 0]
    points = [[point[index] for index in KEPT] for point in points_21]
    ell = [eliminate(k[index]) for index in KEPT]
    images = [
        [
            coordinate + evaluate_rational(component, point)
            for coordinate, component in zip(point, ell)
        ]
        for point in points
    ]
    assert len(set(map(tuple, points))) == 3
    assert images[0] == images[1] == images[2]
    assert [point[0] for point in points] == [0, 1, -1]

    jacobian = [
        [derivative(component, variable) for variable in range(20)]
        for component in ell
    ]
    second = [
        [
            [
                derivative(jacobian[output][first], second_variable)
                for output in range(20)
            ]
            for second_variable in range(20)
        ]
        for first in range(20)
    ]

    profiles = []
    for seed in SEEDS:
        generator = random.Random(seed)
        x = [generator.randrange(1, PRIME) for _ in range(20)]
        y = [generator.randrange(1, PRIME) for _ in range(20)]
        jacobian_value = [
            [evaluate_modular(entry, x) for entry in row]
            for row in jacobian
        ]
        upper_left = [
            [
                sum(
                    y[output]
                    * evaluate_modular(second[first][column][output], x)
                    for output in range(20)
                )
                % PRIME
                for column in range(20)
            ]
            for first in range(20)
        ]
        transpose = [list(column) for column in zip(*jacobian_value)]
        hessian = [
            upper_left[row] + transpose[row] for row in range(20)
        ]
        hessian.extend(
            jacobian_value[row] + [0] * 20 for row in range(20)
        )
        profile = (
            rank(jacobian_value),
            rank(upper_left),
            rank(hessian),
        )
        assert profile == (17, 12, 34)
        profiles.append(profile)

    # Build the coefficient equations for a constant Hessian-kernel vector.
    zero_20 = (0,) * 20
    jacobian_lifted = [
        [
            {exponent + zero_20: coefficient for exponent, coefficient in entry.items()}
            for entry in row
        ]
        for row in jacobian
    ]
    upper_polynomials: list[list[Polynomial]] = []
    for first in range(20):
        row = []
        for column in range(20):
            entry: Polynomial = {}
            for output in range(20):
                for exponent, coefficient in second[first][column][output].items():
                    lifted = list(exponent) + [0] * 20
                    lifted[20 + output] = 1
                    add_term(entry, tuple(lifted), coefficient)
            row.append(entry)
        upper_polynomials.append(row)
    zero: Polynomial = {}
    hessian_polynomials = [
        upper_polynomials[row]
        + [jacobian_lifted[column][row] for column in range(20)]
        for row in range(20)
    ]
    hessian_polynomials.extend(
        jacobian_lifted[row] + [zero] * 20 for row in range(20)
    )
    monomials = sorted(
        {
            exponent
            for row in hessian_polynomials
            for entry in row
            for exponent in entry
        }
    )
    coefficient_rows = [
        [
            residue(entry.get(monomial, Fraction(0)))
            for entry in row
        ]
        for row in hessian_polynomials
        for monomial in monomials
        if any(entry.get(monomial, Fraction(0)) for entry in row)
    ]
    assert rank(coefficient_rows) == 40

    print("PASS independent rank-34 slice: exact relation K_9-3*K_1-K_6=0")
    print("PASS independent rank-34 slice: exact 20D rational collision")
    print("PASS independent rank-34 slice: three modular profiles are (17,12,34)")
    print("PASS independent rank-34 slice: constant Hessian kernel is zero")
    print("PASS independent rank-34 slice: no SymPy or Singular dependency")


if __name__ == "__main__":
    main()
