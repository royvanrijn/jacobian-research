#!/usr/bin/env python3
"""Dependency-free audit of the frozen 22D index-reduced cubic witness.

This replay intentionally starts from the final sparse artifact rather than
importing the SymPy generator.  Exact ``Fraction`` polynomial arithmetic
checks the collision and multiplies the full polynomial Jacobian through
powers 17 and 18.  Thus nilpotence, the Keller property, and index 18 do not
depend on a finite-field sample or on the BCW determinant bridge.
"""

from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path

from audit_shared_bcw_33_independent import (
    Poly,
    add,
    dense,
    derivative,
    evaluate,
    multiply,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "index_reduced_bcw_22_counterexample.json"
)


def decode_h(stored: dict[str, object]) -> list[Poly]:
    dimension = int(stored["dimension"])
    return [
        {
            dense(term["monomial"], dimension): Q(term["coefficient"])
            for term in component
        }
        for component in stored["H"]
    ]


def matrix_product(left: list[list[Poly]], right: list[list[Poly]]) -> list[list[Poly]]:
    dimension = len(left)
    answer: list[list[Poly]] = [
        [{} for _ in range(dimension)] for _ in range(dimension)
    ]
    nonzero_right = [
        [(column, poly) for column, poly in enumerate(row) if poly]
        for row in right
    ]
    for row_index, row in enumerate(left):
        for inner, left_poly in enumerate(row):
            if not left_poly:
                continue
            for column, right_poly in nonzero_right[inner]:
                answer[row_index][column] = add(
                    answer[row_index][column],
                    multiply(left_poly, right_poly),
                )
    return answer


def rational_rank(matrix: list[list[Q]]) -> int:
    work = [row.copy() for row in matrix]
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(pivot_row + 1, len(work)):
            if work[row][column]:
                scale = work[row][column]
                work[row] = [
                    left - scale * right
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def main() -> None:
    stored = json.loads(ARTIFACT.read_text())
    assert stored["format"] == "index-reduced-bcw-sparse-cubic-homogeneous-map-v1"
    dimension = int(stored["dimension"])
    assert dimension == 22
    h = decode_h(stored)
    assert len(h) == dimension
    assert all(sum(exponents) == 3 for poly in h for exponents in poly)

    points = [[Q(value) for value in point] for point in stored["collision_points"]]
    common_image = [Q(value) for value in stored["common_image"]]
    assert len({tuple(point) for point in points}) == 3
    for point in points:
        image = [
            coordinate + evaluate(poly, point)
            for coordinate, poly in zip(point, h)
        ]
        assert image == common_image

    jacobian = [
        [derivative(h[row], column) for column in range(dimension)]
        for row in range(dimension)
    ]
    point = [Q(index * index + 3 * index + 5) for index in range(dimension)]
    evaluated = [
        [evaluate(entry, point) for entry in row] for row in jacobian
    ]
    assert rational_rank(evaluated) == 18

    power = jacobian
    nonzero_counts: list[int] = []
    term_counts: list[int] = []
    for _ in range(1, 19):
        nonzero_counts.append(sum(bool(entry) for row in power for entry in row))
        term_counts.append(sum(len(entry) for row in power for entry in row))
        if len(nonzero_counts) < 18:
            power = matrix_product(power, jacobian)
    assert nonzero_counts[-2:] == [6, 0]
    assert term_counts[-2:] == [8, 0]

    statistics = stored["statistics"]
    assert statistics["generic_rank_JH_over_QQ_x"] == 18
    assert statistics["nilpotency_index_JH"] == 18
    certificate = statistics["exact_certificate"]
    assert certificate["independent_generic_kernel_columns"] == 4
    assert certificate["specialized_rank_lower_bound"] == 18
    assert certificate["nonzero_entries_JH_power_17"] == 6
    assert certificate["nonzero_entries_JH_power_18"] == 0

    print("PASS (stdlib) index-reduced BCW 22: parsed a cubic-homogeneous 22D correction")
    print("PASS (stdlib) index-reduced BCW 22: verified three exact collision points")
    print("PASS (stdlib) index-reduced BCW 22: exact specialized Jacobian rank is 18")
    print("PASS (stdlib) index-reduced BCW 22: (JH)^17 has 6 nonzero entries and (JH)^18=0")
    print("PASS (stdlib) index-reduced BCW 22: direct nilpotence certifies the Keller property")


if __name__ == "__main__":
    main()
