#!/usr/bin/env python3
"""Audit nonlinear semiconjugacies carrying a collision observable.

On the 20-variable identity slice W, both the multiplier coordinate X_0 and
the quadratic observable

    Q = X_18 - X_6*X_8

separate the three stored collision points before W and are constant after
W.  For each observable q, this script evaluates the gradients of

    q, q o W, ..., q o W^12

at one integral point modulo the good prime 1000003.  Their row rank is 13.
A nonzero modular Jacobian minor proves algebraic independence over Q, so any
rational semiconjugate quotient through which either observable factors has
dimension at least 13.

Ranks of longer prefixes at the chosen modular point are printed only as
exploratory diagnostics; they are not upper bounds on generic rank.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "essential_bcw_21_counterexample.json"
)
PRIME = 1_000_003
DIMENSION = 20

StoredTerm = tuple[int, tuple[int, ...]]
Matrix = list[list[int]]


def rational_mod(value: str) -> int:
    number = Fraction(value)
    return number.numerator * pow(number.denominator, -1, PRIME) % PRIME


def decode_slice(source: dict[str, object]) -> list[list[StoredTerm]]:
    components = []
    for component in source["H"][:DIMENSION]:
        terms = []
        for term in component:
            exponent = [0] * 21
            for variable, power in term["monomial"]:
                exponent[variable] = power
            terms.append((rational_mod(term["coefficient"]), tuple(exponent)))
        components.append(terms)
    return components


def evaluate_monomial(exponent: tuple[int, ...], point: list[int]) -> int:
    value = 1
    for variable, coordinate in enumerate(point):
        value = value * pow(coordinate, exponent[variable], PRIME) % PRIME
    # X_20=1 on the identity slice.
    return value


def map_and_jacobian(
    components: list[list[StoredTerm]], point: list[int]
) -> tuple[list[int], Matrix]:
    image = list(point)
    jacobian = [
        [int(row == column) for column in range(DIMENSION)]
        for row in range(DIMENSION)
    ]
    for output, component in enumerate(components):
        for coefficient, exponent in component:
            image[output] = (
                image[output] + coefficient * evaluate_monomial(exponent, point)
            ) % PRIME
            for variable in range(DIMENSION):
                power = exponent[variable]
                if not power:
                    continue
                derivative = list(exponent)
                derivative[variable] -= 1
                jacobian[output][variable] = (
                    jacobian[output][variable]
                    + coefficient
                    * power
                    * evaluate_monomial(tuple(derivative), point)
                ) % PRIME
    return image, jacobian


def matrix_product(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(
                left[row][middle] * right[middle][column]
                for middle in range(DIMENSION)
            )
            % PRIME
            for column in range(DIMENSION)
        ]
        for row in range(DIMENSION)
    ]


def row_times_matrix(row: list[int], matrix: Matrix) -> list[int]:
    return [
        sum(row[middle] * matrix[middle][column] for middle in range(DIMENSION))
        % PRIME
        for column in range(DIMENSION)
    ]


def row_rank(rows: list[list[int]]) -> int:
    matrix = [row[:] for row in rows]
    rank = 0
    for column in range(DIMENSION):
        pivot = next(
            (
                row
                for row in range(rank, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, PRIME)
        matrix[rank] = [value * inverse % PRIME for value in matrix[rank]]
        for row in range(rank + 1, len(matrix)):
            if matrix[row][column]:
                factor = matrix[row][column]
                matrix[row] = [
                    (left - factor * right) % PRIME
                    for left, right in zip(matrix[row], matrix[rank])
                ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def coordinate_zero_gradient(_point: list[int]) -> list[int]:
    return [int(index == 0) for index in range(DIMENSION)]


def quadratic_gradient(point: list[int]) -> list[int]:
    gradient = [0] * DIMENSION
    gradient[18] = 1
    gradient[6] = -point[8] % PRIME
    gradient[8] = -point[6] % PRIME
    return gradient


def observable_rows(
    components: list[list[StoredTerm]],
    gradient,
    length: int,
    seed: int | None = None,
) -> list[list[int]]:
    if seed is None:
        point = [
            (index * index + 3 * index + 5) % PRIME
            for index in range(DIMENSION)
        ]
    else:
        point = [
            int.from_bytes(
                hashlib.sha256(f"{seed}:{index}".encode()).digest()[:8],
                "big",
            )
            % PRIME
            for index in range(DIMENSION)
        ]
    derivative = [
        [int(row == column) for column in range(DIMENSION)]
        for row in range(DIMENSION)
    ]
    rows = []
    for _order in range(length):
        rows.append(row_times_matrix(gradient(point), derivative))
        point, next_jacobian = map_and_jacobian(components, point)
        derivative = matrix_product(next_jacobian, derivative)
    return rows


def observable_ranks(
    components: list[list[StoredTerm]],
    gradient,
    length: int,
    seed: int | None = None,
) -> list[int]:
    rows = observable_rows(components, gradient, length, seed)
    ranks = []
    for prefix in range(1, len(rows) + 1):
        ranks.append(row_rank(rows[:prefix]))
    return ranks


def main() -> None:
    source = json.loads(SOURCE.read_text())
    assert source["dimension"] == 21
    assert source["H"][20] == []
    assert all(point[20] == "1" for point in source["collision_points"])
    components = decode_slice(source)

    collision = [[Fraction(value) for value in point] for point in source["collision_points"]]
    q_values = [point[18] - point[6] * point[8] for point in collision]
    assert q_values == [Fraction(0), Fraction(-39, 16), Fraction(39, 16)]
    x0_values = [point[0] for point in collision]
    assert x0_values == [Fraction(0), Fraction(1), Fraction(-1)]
    target = [Fraction(value) for value in source["common_image"]]
    assert target[0] == 0
    assert target[18] - target[6] * target[8] == 0

    ranks_x0 = observable_ranks(components, coordinate_zero_gradient, 25)
    ranks_q = observable_ranks(components, quadratic_gradient, 25)
    assert ranks_x0[:13] == list(range(1, 14))
    assert ranks_q[:13] == list(range(1, 14))
    seeds = (11, 23, 47, 89, 131, 197)
    plateau_x0 = [
        observable_ranks(components, coordinate_zero_gradient, 25, seed)[-1]
        for seed in seeds
    ]
    plateau_q = [
        observable_ranks(components, quadratic_gradient, 25, seed)[-1]
        for seed in seeds
    ]
    assert plateau_x0 == plateau_q == [13] * len(seeds)
    stacked_x0 = [
        row
        for seed in seeds[:3]
        for row in observable_rows(
            components, coordinate_zero_gradient, 25, seed
        )
    ]
    stacked_q = [
        row
        for seed in seeds[:3]
        for row in observable_rows(
            components, quadratic_gradient, 25, seed
        )
    ]
    assert row_rank(stacked_x0) == row_rank(stacked_q) == DIMENSION

    print("PASS observable quotient audit: X_0 values on collision are 0,1,-1")
    print("PASS observable quotient audit: Q values on collision are 0,-39/16,39/16")
    print("PASS observable quotient audit: both observables vanish on the common image")
    print(f"PASS observable quotient audit: first 13 iterates of X_0 have Jacobian rank 13 mod {PRIME}")
    print(f"PASS observable quotient audit: first 13 iterates of Q have Jacobian rank 13 mod {PRIME}")
    print("THEOREM: any rational semiconjugate quotient carrying X_0 or Q has dimension at least 13")
    print("EXPERIMENT ONLY: ranks through 25 iterates at the base modular point")
    print("  X_0:", ranks_x0)
    print("  Q:", ranks_q)
    print(f"EXPERIMENT ONLY: final ranks at independent seeds {seeds}")
    print("  X_0:", plateau_x0)
    print("  Q:", plateau_q)
    print("PASS observable quotient audit: stacked three-point codistributions have rank 20")
    print("THEOREM: the iterate families have no common constant translation direction")


if __name__ == "__main__":
    main()
