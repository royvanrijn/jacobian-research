#!/usr/bin/env python3
"""Search the first migrating binary GVC jets over faithful finite fields.

This is an experiment, not a theorem checker.  It studies the repeated-root
quartic leading orbit

    Lambda_4 = d_x^4,        P_5 in y^2 Sym^3(x,y),

and the analogous (3,1) and (2,2) monomial orbits.  The first search uses

    Lambda = Lambda_4 + Lambda_5,   P = P_5 + P_4.

The second, conditioned search keeps the (4) orbit and adds the complete
Weierstrass-normalized defect-two data Lambda_6 and P_3.  It solves the
moment-one and moment-two affine equations before sampling the later moments.

Every prime is larger than the maximum input degree in the moments it checks,
so reduction of the characteristic-zero differentiation formulas is faithful.
Passing a bounded modular search is not an all-order GVC statement.
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "binary_repeated_quartic_gvc_jet_search.json"
)

Exponent = tuple[int, int]
Polynomial = dict[Exponent, int]
Operator = dict[Exponent, int]


@dataclass(frozen=True)
class Orbit:
    name: str
    leading_operator: Exponent
    leading_polynomial_branches: tuple[tuple[Exponent, ...], ...]
    quintic_jet_exponents: tuple[Exponent, ...]


ORBITS = (
    Orbit(
        "4",
        (4, 0),
        (tuple((i, 5 - i) for i in range(4)),),
        tuple((i, 5 - i) for i in range(4)),
    ),
    Orbit(
        "3+1",
        (3, 1),
        (((5, 0),), tuple((i, 5 - i) for i in range(3))),
        tuple((i, 5 - i) for i in (0, 1, 2, 5)),
    ),
    Orbit(
        "2+2",
        (2, 2),
        (((4, 1), (5, 0)), ((0, 5), (1, 4))),
        ((0, 5), (1, 4), (4, 1), (5, 0)),
    ),
)


class ModDifferentialRing:
    def __init__(self, prime: int, max_degree: int, max_order: int) -> None:
        self.prime = prime
        self.falling = {
            (degree, order): self._falling_factorial(degree, order)
            for degree in range(max_degree + 1)
            for order in range(max_order + 1)
        }

    def _falling_factorial(self, degree: int, order: int) -> int:
        result = 1
        for offset in range(order):
            result = result * (degree - offset) % self.prime
        return result

    def multiply(self, left: Polynomial, right: Polynomial) -> Polynomial:
        result: Polynomial = {}
        for (left_x, left_y), left_coefficient in left.items():
            for (right_x, right_y), right_coefficient in right.items():
                exponent = (left_x + right_x, left_y + right_y)
                result[exponent] = (
                    result.get(exponent, 0)
                    + left_coefficient * right_coefficient
                ) % self.prime
        return {
            exponent: coefficient
            for exponent, coefficient in result.items()
            if coefficient
        }

    def apply(self, polynomial: Polynomial, operator: Operator) -> Polynomial:
        result: Polynomial = {}
        for (x_degree, y_degree), coefficient in polynomial.items():
            for (x_order, y_order), operator_coefficient in operator.items():
                if x_degree < x_order or y_degree < y_order:
                    continue
                exponent = (x_degree - x_order, y_degree - y_order)
                contribution = (
                    coefficient
                    * operator_coefficient
                    * self.falling[x_degree, x_order]
                    * self.falling[y_degree, y_order]
                )
                result[exponent] = (
                    result.get(exponent, 0) + contribution
                ) % self.prime
        return {
            exponent: coefficient
            for exponent, coefficient in result.items()
            if coefficient
        }

    def moment(
        self,
        polynomial: Polynomial,
        operator: Operator,
        order: int,
        multiplier: Polynomial | None = None,
    ) -> Polynomial:
        result: Polynomial = {(0, 0): 1}
        for _ in range(order):
            result = self.multiply(result, polynomial)
        if multiplier is not None:
            result = self.multiply(multiplier, result)
        for _ in range(order):
            result = self.apply(result, operator)
        return result


def random_projective(
    rng: random.Random,
    prime: int,
    exponents: tuple[Exponent, ...],
) -> Polynomial:
    coefficients = [rng.randrange(prime) for _ in exponents]
    while not any(coefficients):
        coefficients = [rng.randrange(prime) for _ in exponents]
    first = next(coefficient for coefficient in coefficients if coefficient)
    inverse = pow(first, -1, prime)
    return {
        exponent: coefficient * inverse % prime
        for exponent, coefficient in zip(exponents, coefficients)
        if coefficient
    }


def random_affine_nonzero(
    rng: random.Random,
    prime: int,
    exponents: tuple[Exponent, ...],
) -> Polynomial:
    coefficients = [rng.randrange(prime) for _ in exponents]
    while not any(coefficients):
        coefficients = [rng.randrange(prime) for _ in exponents]
    return {
        exponent: coefficient
        for exponent, coefficient in zip(exponents, coefficients)
        if coefficient
    }


def serialize_polynomial(polynomial: Polynomial) -> dict[str, int]:
    return {
        f"{x_degree},{y_degree}": coefficient
        for (x_degree, y_degree), coefficient in sorted(polynomial.items())
    }


def solve_affine_system(
    rng: random.Random,
    prime: int,
    constant: list[int],
    columns: list[list[int]],
) -> list[int] | None:
    """Choose a random point of constant + span(columns) = 0."""
    matrix = [
        [
            *(columns[column][row] % prime for column in range(len(columns))),
            -constant[row] % prime,
        ]
        for row in range(len(constant))
    ]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(len(columns)):
        candidate = next(
            (
                row
                for row in range(pivot_row, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if candidate is None:
            continue
        matrix[pivot_row], matrix[candidate] = (
            matrix[candidate],
            matrix[pivot_row],
        )
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [
            value * inverse % prime for value in matrix[pivot_row]
        ]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            scalar = matrix[row][column]
            matrix[row] = [
                (left - scalar * right) % prime
                for left, right in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    if any(not any(row[:-1]) and row[-1] for row in matrix):
        return None

    free_columns = [
        column
        for column in range(len(columns))
        if column not in pivot_columns
    ]
    solution = [0] * len(columns)
    for column in free_columns:
        solution[column] = rng.randrange(prime)
    for row, column in enumerate(pivot_columns):
        solution[column] = (
            matrix[row][-1]
            - sum(
                matrix[row][free] * solution[free]
                for free in free_columns
            )
        ) % prime
    return solution


def depth_one_search(
    orbit: Orbit,
    *,
    prime: int,
    samples: int,
    seed: int,
) -> dict[str, object]:
    max_moment = 4
    assert prime > 5 * max_moment
    ring = ModDifferentialRing(prime, 5 * max_moment, 5)
    rng = random.Random(seed)
    moment_counts = {str(moment): 0 for moment in range(2, max_moment + 1)}
    survivors: list[dict[str, object]] = []
    leading_factorial = (
        ring.falling[
            orbit.leading_operator[0], orbit.leading_operator[0]
        ]
        * ring.falling[
            orbit.leading_operator[1], orbit.leading_operator[1]
        ]
    ) % prime

    for sample in range(samples):
        branch = orbit.leading_polynomial_branches[
            sample % len(orbit.leading_polynomial_branches)
        ]
        polynomial_5 = random_projective(rng, prime, branch)
        operator_5 = random_affine_nonzero(
            rng, prime, orbit.quintic_jet_exponents
        )
        polynomial_4 = {
            (x_degree, 4 - x_degree): rng.randrange(prime)
            for x_degree in range(5)
            if (x_degree, 4 - x_degree) != orbit.leading_operator
        }
        apolar_value = sum(
            coefficient
            * operator_5.get(exponent, 0)
            * ring.falling[exponent[0], exponent[0]]
            * ring.falling[exponent[1], exponent[1]]
            for exponent, coefficient in polynomial_5.items()
        ) % prime
        polynomial_4[orbit.leading_operator] = (
            -apolar_value * pow(leading_factorial, -1, prime)
        ) % prime

        polynomial = polynomial_5 | {
            exponent: coefficient
            for exponent, coefficient in polynomial_4.items()
            if coefficient
        }
        operator = {orbit.leading_operator: 1} | operator_5
        survived = True
        for moment_order in range(2, max_moment + 1):
            if ring.moment(polynomial, operator, moment_order):
                survived = False
                break
            moment_counts[str(moment_order)] += 1
        if survived:
            mixed = {}
            for name, multiplier in {
                "x": {(1, 0): 1},
                "y": {(0, 1): 1},
            }.items():
                mixed[name] = [
                    bool(
                        ring.moment(
                            polynomial,
                            operator,
                            moment_order,
                            multiplier,
                        )
                    )
                    for moment_order in range(1, max_moment + 1)
                ]
            survivors.append(
                {
                    "P5": serialize_polynomial(polynomial_5),
                    "P4": serialize_polynomial(
                        {
                            exponent: coefficient
                            for exponent, coefficient in polynomial_4.items()
                            if coefficient
                        }
                    ),
                    "Lambda5": serialize_polynomial(operator_5),
                    "boundary_coefficients": {
                        "P4_dual_to_Lambda4": polynomial_4.get(
                            orbit.leading_operator, 0
                        ),
                        "Lambda5_dual_to_lowest_P5_tip": operator_5.get(
                            branch[0], 0
                        ),
                    },
                    "mixed_nonzero": mixed,
                }
            )

    return {
        "orbit": orbit.name,
        "prime": prime,
        "samples": samples,
        "seed": seed,
        "maximum_moment": max_moment,
        "faithful_degree_bound": 5 * max_moment,
        "moment_survivor_counts": moment_counts,
        "survivor_count": len(survivors),
        "survivors": survivors,
    }


def _early_vector(
    ring: ModDifferentialRing,
    polynomial: Polynomial,
    operator: Operator,
) -> list[int]:
    first = ring.moment(polynomial, operator, 1)
    second = ring.moment(polynomial, operator, 2)
    return [
        first.get((0, 0), 0),
        second.get((0, 1), 0),
        second.get((1, 0), 0),
    ]


def depth_one_conditioned_orbit_four_search(
    *,
    prime: int,
    samples: int,
    seed: int,
) -> dict[str, object]:
    """Push the minimal (4)-orbit ansatz through a faithful fifth moment."""
    max_moment = 5
    assert prime > 5 * max_moment
    ring = ModDifferentialRing(prime, 5 * max_moment + 1, 5)
    rng = random.Random(seed)
    polynomial_5_exponents = tuple((i, 5 - i) for i in range(4))
    polynomial_4_exponents = tuple((i, 4 - i) for i in range(5))
    leading_operator = {(4, 0): 1}
    counts = {
        "solvable_through_defect_one": 0,
        "moment_two": 0,
        "moment_three": 0,
        "moment_four": 0,
        "moment_five": 0,
    }
    survivors: list[dict[str, object]] = []

    for _sample in range(samples):
        polynomial_5 = random_projective(
            rng, prime, polynomial_5_exponents
        )
        operator_5 = random_affine_nonzero(
            rng, prime, polynomial_5_exponents
        )
        operator = leading_operator | operator_5
        constant = _early_vector(ring, polynomial_5, operator)
        columns = []
        for exponent in polynomial_4_exponents:
            value = _early_vector(
                ring, polynomial_5 | {exponent: 1}, operator
            )
            columns.append(
                [
                    (entry - base) % prime
                    for entry, base in zip(value, constant)
                ]
            )
        polynomial_4_solution = solve_affine_system(
            rng, prime, constant, columns
        )
        if polynomial_4_solution is None:
            continue
        counts["solvable_through_defect_one"] += 1
        polynomial_4 = {
            exponent: coefficient
            for exponent, coefficient in zip(
                polynomial_4_exponents, polynomial_4_solution
            )
            if coefficient
        }
        polynomial = polynomial_5 | polynomial_4
        survived = True
        for moment_order in range(2, max_moment + 1):
            if ring.moment(polynomial, operator, moment_order):
                survived = False
                break
            moment_name = ("two", "three", "four", "five")[
                moment_order - 2
            ]
            counts[f"moment_{moment_name}"] += 1
        if not survived:
            continue
        mixed = {}
        for name, multiplier in {
            "x": {(1, 0): 1},
            "y": {(0, 1): 1},
        }.items():
            mixed[name] = [
                bool(
                    ring.moment(
                        polynomial,
                        operator,
                        moment_order,
                        multiplier,
                    )
                )
                for moment_order in range(1, max_moment + 1)
            ]
        survivors.append(
            {
                "P5": serialize_polynomial(polynomial_5),
                "P4": serialize_polynomial(polynomial_4),
                "Lambda5": serialize_polynomial(operator_5),
                "support_extrema": {
                    "P5_max_x": max(exponent[0] for exponent in polynomial_5),
                    "P4_max_x": max(
                        (exponent[0] for exponent in polynomial_4),
                        default=-1,
                    ),
                    "Lambda5_min_x": min(
                        exponent[0] for exponent in operator_5
                    ),
                },
                "mixed_nonzero": mixed,
            }
        )

    return {
        "orbit": "4",
        "prime": prime,
        "samples": samples,
        "seed": seed,
        "maximum_moment": max_moment,
        "faithful_degree_bound": 5 * max_moment,
        "conditioning": (
            "solve moment one and both degree-one coefficients of moment "
            "two in P4, then test the scalar second moment and moments "
            "three through five"
        ),
        "counts": counts,
        "survivor_count": len(survivors),
        "survivors": survivors,
    }


def depth_two_conditioned_search(
    *,
    prime: int,
    samples: int,
    seed: int,
) -> dict[str, object]:
    max_moment = 5
    assert prime > 5 * max_moment
    ring = ModDifferentialRing(prime, 5 * max_moment, 6)
    rng = random.Random(seed)
    polynomial_5_exponents = tuple((i, 5 - i) for i in range(4))
    polynomial_4_exponents = tuple((i, 4 - i) for i in range(5))
    polynomial_3_exponents = tuple((i, 3 - i) for i in range(4))
    operator_6_exponents = tuple((i, 6 - i) for i in range(4))
    leading_operator = {(4, 0): 1}

    counts = {
        "solvable_through_defect_one": 0,
        "solvable_through_moment_two": 0,
        "moment_three": 0,
        "moment_four": 0,
        "moment_five": 0,
    }
    delayed: list[dict[str, object]] = []

    for _sample in range(samples):
        polynomial_5 = random_projective(
            rng, prime, polynomial_5_exponents
        )
        operator_5 = random_affine_nonzero(
            rng, prime, polynomial_5_exponents
        )
        operator_through_5 = leading_operator | operator_5

        constant = _early_vector(
            ring, polynomial_5, operator_through_5
        )
        columns = []
        for exponent in polynomial_4_exponents:
            value = _early_vector(
                ring,
                polynomial_5 | {exponent: 1},
                operator_through_5,
            )
            columns.append(
                [
                    (entry - base) % prime
                    for entry, base in zip(value, constant)
                ]
            )
        polynomial_4_solution = solve_affine_system(
            rng, prime, constant, columns
        )
        if polynomial_4_solution is None:
            continue
        counts["solvable_through_defect_one"] += 1
        polynomial_4 = {
            exponent: coefficient
            for exponent, coefficient in zip(
                polynomial_4_exponents, polynomial_4_solution
            )
            if coefficient
        }
        base_polynomial = polynomial_5 | polynomial_4

        second_constant = ring.moment(
            base_polynomial, operator_through_5, 2
        ).get((0, 0), 0)
        second_columns: list[list[int]] = []
        for exponent in polynomial_3_exponents:
            value = ring.moment(
                base_polynomial | {exponent: 1},
                operator_through_5,
                2,
            ).get((0, 0), 0)
            second_columns.append([(value - second_constant) % prime])
        for exponent in operator_6_exponents:
            value = ring.moment(
                base_polynomial,
                operator_through_5 | {exponent: 1},
                2,
            ).get((0, 0), 0)
            second_columns.append([(value - second_constant) % prime])
        depth_two_solution = solve_affine_system(
            rng, prime, [second_constant], second_columns
        )
        if depth_two_solution is None:
            continue
        counts["solvable_through_moment_two"] += 1
        polynomial_3 = {
            exponent: coefficient
            for exponent, coefficient in zip(
                polynomial_3_exponents, depth_two_solution[:4]
            )
            if coefficient
        }
        operator_6 = {
            exponent: coefficient
            for exponent, coefficient in zip(
                operator_6_exponents, depth_two_solution[4:]
            )
            if coefficient
        }
        polynomial = base_polynomial | polynomial_3
        operator = operator_through_5 | operator_6

        third = ring.moment(polynomial, operator, 3)
        if third:
            continue
        counts["moment_three"] += 1
        fourth = ring.moment(polynomial, operator, 4)
        if fourth:
            delayed.append(
                {
                    "P5": serialize_polynomial(polynomial_5),
                    "P4": serialize_polynomial(polynomial_4),
                    "P3": serialize_polynomial(polynomial_3),
                    "Lambda5": serialize_polynomial(operator_5),
                    "Lambda6": serialize_polynomial(operator_6),
                    "first_failure": 4,
                    "fourth_moment": serialize_polynomial(fourth),
                }
            )
            continue
        counts["moment_four"] += 1
        fifth = ring.moment(polynomial, operator, 5)
        if fifth:
            delayed.append(
                {
                    "P5": serialize_polynomial(polynomial_5),
                    "P4": serialize_polynomial(polynomial_4),
                    "P3": serialize_polynomial(polynomial_3),
                    "Lambda5": serialize_polynomial(operator_5),
                    "Lambda6": serialize_polynomial(operator_6),
                    "first_failure": 5,
                    "fifth_moment": serialize_polynomial(fifth),
                }
            )
            continue
        counts["moment_five"] += 1

    return {
        "orbit": "4",
        "prime": prime,
        "samples": samples,
        "seed": seed,
        "maximum_moment": max_moment,
        "faithful_degree_bound": 5 * max_moment,
        "conditioning": (
            "solve moment one and both degree-one coefficients of moment "
            "two in P4, then solve the scalar part of moment two in P3 "
            "and Lambda6"
        ),
        "counts": counts,
        "delayed_failures": delayed,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="use reduced sample counts and do not overwrite the pinned result",
    )
    args = parser.parse_args()

    if args.quick:
        depth_one_counts = {"4": 20_000, "3+1": 20_000, "2+2": 20_000}
        conditioned_depth_one_count = 20_000
        depth_two_count = 20_000
    else:
        depth_one_counts = {"4": 100_000, "3+1": 100_000, "2+2": 100_000}
        conditioned_depth_one_count = 500_000
        depth_two_count = 500_000

    seeds = {"4": 20260733, "3+1": 20260732, "2+2": 20260733}
    depth_one = [
        depth_one_search(
            orbit,
            prime=23,
            samples=depth_one_counts[orbit.name],
            seed=seeds[orbit.name],
        )
        for orbit in ORBITS
    ]
    conditioned_depth_one = depth_one_conditioned_orbit_four_search(
        prime=29,
        samples=conditioned_depth_one_count,
        seed=20260738,
    )
    depth_two = depth_two_conditioned_search(
        prime=29,
        samples=depth_two_count,
        seed=20260737,
    )

    artifact = {
        "format": "binary-repeated-quartic-gvc-jet-search-v1",
        "status": "experiment",
        "field": "finite fields of faithful characteristic for each window",
        "scope": (
            "bounded modular samples on repeated-root quartic-leading "
            "degree-five binary GVC jets; not an exhaustive search, "
            "all-order proof, or counterexample"
        ),
        "minimal_migrating_architecture": {
            "operator": "Lambda4 + Lambda5",
            "polynomial": "P5 + P4",
            "defect_formula": (
                "number of Lambda5 selections plus number of P4 selections"
            ),
        },
        "depth_one_searches": depth_one,
        "conditioned_depth_one_fifth_moment_search": conditioned_depth_one,
        "depth_two_conditioned_search": depth_two,
        "conclusion": (
            "No sampled point survives the declared final moment outside "
            "the recorded boundary strata. The search does not exclude "
            "points outside the samples or failure at a later moment."
        ),
    }

    if args.quick:
        print(json.dumps(artifact, indent=2))
        print("PASS quick repeated-quartic GVC jet experiment")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS repeated-quartic GVC jet experiment")
    print(f"PASS wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
