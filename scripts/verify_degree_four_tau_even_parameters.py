#!/usr/bin/env python3
"""Exact low-degree parameters for the apolar-even quartic quotient.

The script constructs a concrete set of 22 algebraically independent
tau-even SL2 invariants on

    End(Sym^4) = Sym^0 + Sym^2 + Sym^4 + Sym^6 + Sym^8.

It uses the scalar invariant, four primitive quadratic Casimir contractions,
nine tau-even primitive cubic contractions, and greedily selected
tau-even quartic trace words.  A nonzero 22-by-22 Jacobian minor
modulo 1000003 is an exact characteristic-zero independence certificate.

This produces rational parameters for a finite cover of the tau-fixed
quotient.  It does not prove that these parameters, or the moments, generate
the full tau-fixed function field.
"""

from __future__ import annotations

import itertools
import json
from math import factorial
from pathlib import Path

import sympy as sp

import verify_two_pair_counterexample_missing_invariant as base


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_four_tau_even_parameters.json"
)
PRIME = 1_000_003
EVEN_CUBIC_TRIPLES = tuple(
    triple for triple in base.CUBIC_TRIPLES if triple != (2, 3, 4)
)


def reduce_rational(value: sp.Rational) -> int:
    numerator, denominator = sp.fraction(value)
    return (
        int(numerator)
        * pow(int(denominator) % PRIME, -1, PRIME)
        % PRIME
    )


def matrix_multiply(
    left: list[list[int]],
    right: list[list[int]],
) -> list[list[int]]:
    return [
        [
            sum(left[row][middle] * right[middle][column] for middle in range(5))
            % PRIME
            for column in range(5)
        ]
        for row in range(5)
    ]


def trace(matrix: list[list[int]]) -> int:
    return sum(matrix[index][index] for index in range(5)) % PRIME


def trace_word(matrices: list[list[list[int]]]) -> int:
    product = matrices[0]
    for matrix in matrices[1:]:
        product = matrix_multiply(product, matrix)
    return trace(product)


def rank_and_pivots(
    rows: list[list[int]],
) -> tuple[int, list[int], int | None]:
    reduced = [[value % PRIME for value in row] for row in rows]
    rank = 0
    pivots = []
    determinant = 1
    for column in range(len(reduced[0])):
        pivot = next(
            (
                row
                for row in range(rank, len(reduced))
                if reduced[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        if pivot != rank:
            reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
            determinant = -determinant
        pivot_value = reduced[rank][column]
        determinant = determinant * pivot_value % PRIME
        inverse = pow(pivot_value, -1, PRIME)
        for row in range(rank + 1, len(reduced)):
            if not reduced[row][column]:
                continue
            scale = reduced[row][column] * inverse % PRIME
            reduced[row] = [
                (left - scale * right) % PRIME
                for left, right in zip(reduced[row], reduced[rank], strict=True)
            ]
        pivots.append(column)
        rank += 1
        if rank == len(reduced):
            break
    if rank != len(rows):
        return rank, pivots, None
    return rank, pivots, determinant % PRIME


def generic_point() -> list[list[int]]:
    return [
        [
            (
                19
                + 37 * row
                + 61 * column
                + 11 * row * column
                + 7 * row**2
                + 13 * column**2
            )
            % PRIME
            for column in range(5)
        ]
        for row in range(5)
    ]


def moment_jacobian_rows(point: list[list[int]], cutoff: int) -> list[list[int]]:
    factorials = [1]
    for value in range(1, 4 * cutoff + 1):
        factorials.append(factorials[-1] * value % PRIME)
    powers = [[[1]]]
    for order in range(1, cutoff):
        previous = powers[-1]
        size = 4 * order + 1
        current = [[0] * size for _ in range(size)]
        for left, previous_row in enumerate(previous):
            for right, value in enumerate(previous_row):
                if not value:
                    continue
                for row in range(5):
                    for column in range(5):
                        current[left + row][right + column] = (
                            current[left + row][right + column]
                            + value * point[row][column]
                        ) % PRIME
        powers.append(current)

    result = []
    for order in range(1, cutoff + 1):
        previous = powers[order - 1]
        row_values = []
        for dual_index in range(5):
            for coordinate_index in range(5):
                value = 0
                for total in range(4 * order + 1):
                    left = total - dual_index
                    right = total - coordinate_index
                    if (
                        0 <= left < len(previous)
                        and 0 <= right < len(previous)
                    ):
                        value += (
                            factorials[total]
                            * factorials[4 * order - total]
                            * previous[left][right]
                        )
                row_values.append(order * value % PRIME)
        result.append(row_values)
    return result


def component_data() -> tuple[
    dict[int, list[list[int]]],
    list[dict[int, list[list[int]]]],
]:
    projectors = {
        component: [
            [reduce_rational(entry) for entry in row]
            for row in projector.tolist()
        ]
        for component, projector in base.casimir_projectors().items()
    }
    diagonal = [
        factorial(index) * factorial(4 - index) % PRIME
        for index in range(5)
    ]

    def operator_vector(point: list[list[int]]) -> list[int]:
        return [
            point[column][row] * diagonal[column] % PRIME
            for column in range(5)
            for row in range(5)
        ]

    def project(vector: list[int]) -> dict[int, list[list[int]]]:
        result = {}
        for component, projector in projectors.items():
            projected = [
                sum(a * b for a, b in zip(row, vector, strict=True)) % PRIME
                for row in projector
            ]
            result[component] = [
                [projected[row + 5 * column] for column in range(5)]
                for row in range(5)
            ]
        return result

    point_components = project(operator_vector(generic_point()))
    direction_components = []
    for coefficient in range(25):
        direction = [[0] * 5 for _ in range(5)]
        direction[coefficient // 5][coefficient % 5] = 1
        direction_components.append(project(operator_vector(direction)))
    return point_components, direction_components


def trace_word_derivative(
    word: tuple[int, ...],
    point: dict[int, list[list[int]]],
    direction: dict[int, list[list[int]]],
) -> int:
    result = 0
    for differentiated in range(len(word)):
        matrices = [
            (
                direction[component]
                if index == differentiated
                else point[component]
            )
            for index, component in enumerate(word)
        ]
        result += trace_word(matrices)
    return result % PRIME


def symmetrized_trace_derivative(
    words: list[tuple[int, ...]],
    point: dict[int, list[list[int]]],
    direction: dict[int, list[list[int]]],
) -> int:
    return sum(
        trace_word_derivative(word, point, direction) for word in words
    ) % PRIME


def canonical_cyclic_word(word: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [
        word[index:] + word[:index] for index in range(len(word))
    ]
    return min(rotations)


def main() -> None:
    point, directions = component_data()
    rows: list[list[int]] = []
    labels: list[str] = []
    degrees: list[int] = []

    # The scalar component is a nonzero multiple of mu_1.
    rows.append(
        [
            trace(direction[0])
            for direction in directions
        ]
    )
    labels.append("tr(A_0)")
    degrees.append(1)

    for component in range(1, 5):
        rows.append(
            [
                trace_word_derivative(
                    (component, component),
                    point,
                    direction,
                )
                for direction in directions
            ]
        )
        labels.append(f"tr(A_{2 * component}^2)")
        degrees.append(2)

    for triple in EVEN_CUBIC_TRIPLES:
        permutations = sorted(set(itertools.permutations(triple)))
        rows.append(
            [
                symmetrized_trace_derivative(
                    permutations,
                    point,
                    direction,
                )
                for direction in directions
            ]
        )
        labels.append("symtr(" + ",".join(map(str, triple)) + ")")
        degrees.append(3)

    rank, _, _ = rank_and_pivots(rows)
    assert rank == 14

    cyclic_words = sorted(
        {
            canonical_cyclic_word(word)
            for word in itertools.product(range(1, 5), repeat=4)
        }
    )
    selected_quartics = []
    for word in cyclic_words:
        if sum(word) % 2:
            continue
        row = [
            trace_word_derivative(word, point, direction)
            for direction in directions
        ]
        if not any(row):
            continue
        candidate_rank, _, _ = rank_and_pivots(rows + [row])
        if candidate_rank > rank:
            rows.append(row)
            labels.append(
                "tr("
                + "".join(map(str, word))
                + ")"
            )
            degrees.append(4)
            selected_quartics.append(
                {
                    "word": list(word),
                    "component_index_sum": sum(word),
                }
            )
            rank = candidate_rank
        if rank == 22:
            break

    assert rank == 22
    assert len(selected_quartics) == 8
    rank, pivots, determinant = rank_and_pivots(rows)
    assert rank == 22
    assert determinant not in (None, 0)

    # The first odd cubic has nonzero differential but cannot increase the
    # quotient rank beyond 22.
    odd_word = (2, 3, 4)
    odd_row = [
        trace_word_derivative(odd_word, point, direction)
        for direction in directions
    ]
    odd_augmented_rank, _, _ = rank_and_pivots(rows + [odd_row])
    assert odd_augmented_rank == 22

    moment_rows = moment_jacobian_rows(generic_point(), 22)
    moment_rank, moment_pivots, moment_determinant = rank_and_pivots(
        moment_rows
    )
    assert moment_rank == 22
    assert moment_determinant not in (None, 0)
    combined_cotangent_rank, _, _ = rank_and_pivots(rows + moment_rows)
    assert combined_cotangent_rank == 22

    payload = {
        "format": "degree-four-tau-even-parameters-v1",
        "prime": PRIME,
        "point": generic_point(),
        "parameter_count": len(rows),
        "degrees": degrees,
        "labels": labels,
        "selected_tau_even_quartic_trace_words": selected_quartics,
        "jacobian_rank": rank,
        "pivot_columns_zero_based": pivots,
        "minor_determinant_mod_prime": determinant,
        "first_tau_odd_generator": "tr(A_4 A_6 A_8)",
        "rank_after_appending_odd_generator": odd_augmented_rank,
        "first_22_moment_jacobian": {
            "rank": moment_rank,
            "pivot_columns_zero_based": moment_pivots,
            "minor_determinant_mod_prime": moment_determinant,
        },
        "combined_even_parameter_and_moment_cotangent_rank": (
            combined_cotangent_rank
        ),
        "status": (
            "exact algebraic-independence certificate; no assertion that "
            "these parameters or the moments generate the tau-fixed field"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(
        "PASS tau-even parameters: degrees",
        degrees,
        "have Jacobian rank",
        rank,
        "modulo",
        PRIME,
    )
    print("PASS tau-even parameters: selected quartic words", selected_quartics)
    print("PASS tau-even parameters: minor determinant", determinant)
    print(
        "PASS tau-even parameters: first 22 moments have the same "
        "22-dimensional generic cotangent space"
    )


if __name__ == "__main__":
    main()
