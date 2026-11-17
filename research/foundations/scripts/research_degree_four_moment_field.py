#!/usr/bin/env python3
"""Bounded exact searches for the degree-four moment field.

This is an exploratory checker, not a proof that the full moment field has
generic degree two.  It searches, modulo a good prime, for homogeneous
relations of the form

    Q(mu_1, mu_2, ...) + c_234^2 P(mu_1, mu_2, ...) = 0,

where c_234 is the first apolar-adjoint-odd primitive invariant.  A full
column rank proves that no characteristic-zero relation with that support
exists.  A rank defect is only a candidate until its coefficients have
been reconstructed over QQ and the resulting polynomial identity verified.
"""

from __future__ import annotations

import argparse
import json
from math import factorial
from pathlib import Path

import sympy as sp

import verify_two_pair_counterexample_missing_invariant as base


DEFAULT_PRIME = 1_000_003
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_four_moment_field_bounded_relations.json"
)


def weighted_exponents(total: int) -> list[tuple[int, ...]]:
    """Exponent vectors for variables of weights 1,...,total."""
    result: list[tuple[int, ...]] = []

    def recurse(weight: int, remainder: int, reversed_prefix: list[int]) -> None:
        if weight == 0:
            if remainder == 0:
                result.append(tuple(reversed(reversed_prefix)))
            return
        for exponent in range(remainder // weight + 1):
            recurse(
                weight - 1,
                remainder - exponent * weight,
                reversed_prefix + [exponent],
            )

    recurse(total, total, [])
    return result


def reduce_rational(value: sp.Rational, prime: int) -> int:
    numerator, denominator = sp.fraction(value)
    return int(numerator) % prime * pow(int(denominator) % prime, -1, prime) % prime


def projector_matrices_mod(prime: int) -> dict[int, list[list[int]]]:
    return {
        component: [
            [reduce_rational(entry, prime) for entry in row]
            for row in projector.tolist()
        ]
        for component, projector in base.casimir_projectors().items()
    }


def low_degree_invariants_mod(
    point: list[list[int]],
    projectors: dict[int, list[list[int]]],
    prime: int,
) -> tuple[list[int], int]:
    factorial_diagonal = [
        factorial(index) * factorial(4 - index) % prime
        for index in range(5)
    ]
    # Column-major vectorization of A=C^T D.
    vector = [
        point[column][row] * factorial_diagonal[column] % prime
        for column in range(5)
        for row in range(5)
    ]
    components: dict[int, list[list[int]]] = {}
    for component in range(5):
        projected = [
            sum(left * right for left, right in zip(row, vector)) % prime
            for row in projectors[component]
        ]
        components[component] = [
            [projected[row + 5 * column] for column in range(5)]
            for row in range(5)
        ]

    quadratics = []
    for component in range(5):
        matrix = components[component]
        quadratics.append(
            sum(
                matrix[row][column] * matrix[column][row]
                for row in range(5)
                for column in range(5)
            )
            % prime
        )

    product_23 = [
        [
            sum(
                components[2][row][middle] * components[3][middle][column]
                for middle in range(5)
            )
            % prime
            for column in range(5)
        ]
        for row in range(5)
    ]
    odd_cubic = (
        sum(
            product_23[row][middle] * components[4][middle][row]
            for row in range(5)
            for middle in range(5)
        )
        % prime
    )
    return quadratics, odd_cubic


def moments_mod(
    point: list[list[int]],
    cutoff: int,
    prime: int,
) -> list[int]:
    factorials = [1]
    for value in range(1, 4 * cutoff + 1):
        factorials.append(factorials[-1] * value % prime)
    power = [[1]]
    values = []
    for order in range(1, cutoff + 1):
        size = 4 * order + 1
        current = [[0] * size for _ in range(size)]
        for left, previous_row in enumerate(power):
            for right, value in enumerate(previous_row):
                if not value:
                    continue
                for row in range(5):
                    for column in range(5):
                        current[left + row][right + column] = (
                            current[left + row][right + column]
                            + value * point[row][column]
                        ) % prime
        power = current
        values.append(
            sum(
                factorials[index]
                * factorials[4 * order - index]
                * power[index][index]
                for index in range(size)
            )
            % prime
        )
    return values


def monomial_value(
    values: list[int],
    exponents: tuple[int, ...],
    prime: int,
) -> int:
    result = 1
    for value, exponent in zip(values, exponents, strict=True):
        if exponent:
            result = result * pow(value, exponent, prime) % prime
    return result


def deterministic_point(sample: int, prime: int) -> list[list[int]]:
    # SplitMix64 gives independently varying, reproducible coordinates.
    state = sample & ((1 << 64) - 1)
    result = []
    for _row in range(5):
        row_values = []
        for _column in range(5):
            state = (state + 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
            value = state
            value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9
            value &= (1 << 64) - 1
            value = (value ^ (value >> 27)) * 0x94D049BB133111EB
            value &= (1 << 64) - 1
            value ^= value >> 31
            row_values.append(value % prime)
        result.append(row_values)
    return result


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    reduced = [row[:] for row in matrix]
    row_count = len(reduced)
    column_count = len(reduced[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(rank, row_count)
                if reduced[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
        inverse = pow(reduced[rank][column], -1, prime)
        for row in range(rank + 1, row_count):
            if not reduced[row][column]:
                continue
            scale = reduced[row][column] * inverse % prime
            reduced[row] = [
                (left - scale * right) % prime
                for left, right in zip(reduced[row], reduced[rank], strict=True)
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def target_relation_rank(
    sample_data: list[tuple[list[int], list[int], int]],
    exponent_cache: dict[int, list[tuple[int, ...]]],
    weight: int,
    target_degree: int,
    target_index: int | None,
    extra_samples: int,
    prime: int,
) -> tuple[int, int]:
    pure_exponents = exponent_cache[weight]
    target_exponents = exponent_cache[weight - target_degree]
    matrix = []
    for moments, quadratics, odd_cubic in sample_data[
        : len(pure_exponents) + len(target_exponents) + extra_samples
    ]:
        pure_values = [
            monomial_value(moments[:weight], exponents, prime)
            for exponents in pure_exponents
        ]
        if target_index is None:
            target = odd_cubic * odd_cubic % prime
        else:
            target = quadratics[target_index]
        target_values = [
            target
            * monomial_value(
                moments[: weight - target_degree],
                exponents,
                prime,
            )
            % prime
            for exponents in target_exponents
        ]
        matrix.append(pure_values + target_values)
    return rank_mod(matrix, prime), len(matrix[0])


def search(
    max_weight: int,
    extra_samples: int,
    prime: int,
    targets: str,
) -> dict[str, list[dict[str, int]]]:
    projectors = projector_matrices_mod(prime)
    exponent_cache = {
        weight: weighted_exponents(weight) for weight in range(max_weight + 1)
    }
    largest_column_count = len(exponent_cache[max_weight])
    if max_weight >= 2:
        largest_column_count += len(exponent_cache[max_weight - 2])
    sample_count = largest_column_count + extra_samples

    sample_data = []
    for sample in range(sample_count):
        point = deterministic_point(sample + 1, prime)
        sample_data.append(
            (
                moments_mod(point, max_weight, prime),
                *low_degree_invariants_mod(point, projectors, prime),
            )
        )

    results: dict[str, list[dict[str, int]]] = {}
    if targets in ("all", "quadratics"):
        for component in range(1, 5):
            label = f"q_{2 * component}"
            results[label] = []
            for weight in range(2, max_weight + 1):
                rank, columns = target_relation_rank(
                    sample_data,
                    exponent_cache,
                    weight,
                    2,
                    component,
                    extra_samples,
                    prime,
                )
                record = {
                    "weight": weight,
                    "columns": columns,
                    "rank": rank,
                    "nullity": columns - rank,
                }
                results[label].append(record)
                print(
                    f"target={label} weight={weight} "
                    f"columns={columns} rank={rank} "
                    f"nullity={columns - rank}"
                )
                if rank < columns:
                    break

    if targets in ("all", "odd-square"):
        label = "c_234^2"
        results[label] = []
        for weight in range(6, max_weight + 1):
            rank, columns = target_relation_rank(
                sample_data,
                exponent_cache,
                weight,
                6,
                None,
                extra_samples,
                prime,
            )
            record = {
                "weight": weight,
                "columns": columns,
                "rank": rank,
                "nullity": columns - rank,
            }
            results[label].append(record)
            print(
                f"target={label} weight={weight} "
                f"columns={columns} rank={rank} "
                f"nullity={columns - rank}"
            )
            if rank < columns:
                break
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-weight", type=int, default=12)
    parser.add_argument("--extra-samples", type=int, default=3)
    parser.add_argument("--prime", type=int, default=DEFAULT_PRIME)
    parser.add_argument(
        "--targets",
        choices=("all", "quadratics", "odd-square"),
        default="all",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = search(
        args.max_weight,
        args.extra_samples,
        args.prime,
        args.targets,
    )
    payload = {
        "format": "degree-four-moment-field-bounded-relations-v1",
        "prime": args.prime,
        "max_weight": args.max_weight,
        "extra_samples": args.extra_samples,
        "targets": args.targets,
        "relation_ansatz": (
            "Q(mu)+target*P(mu)=0, homogeneous in invariant weight"
        ),
        "results": results,
        "status": (
            "full column rank excludes the displayed characteristic-zero "
            "relation support; a rank defect would require rational "
            "reconstruction and identity verification"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
