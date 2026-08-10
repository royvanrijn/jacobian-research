#!/usr/bin/env python3
"""Fast regressions for the generic sparse finite-field circuit backend."""

from __future__ import annotations

from dataclasses import dataclass

from sparse_circuit_modp import (
    apply_scaled_correction,
    apply_row_functionals,
    evaluate,
    evaluate_truncated_series,
    evaluate_univariate,
    interpolate_consecutive_values,
    interpolate_corrections,
    left_cokernel_basis,
    polynomial_gcd_modp,
    solve_linearization,
)


@dataclass
class TinyDAG:
    nodes: list[tuple[object, ...]]


def rational_embedding(key: object) -> int:
    return int(key[0][0]) % 31


def main() -> None:
    zero = ((0, 1), (0, 1))
    minus_one = ((-1, 1), (0, 1))
    dag = TinyDAG(
        [
            ("const", zero),
            ("var", "x"),
            ("var", "y"),
            ("mul", 1, 1),
            ("add", 3, 2),       # x^2+y
            ("scale", minus_one, 1),
            ("add", 2, 5),       # y-x
        ]
    )
    point = {"x": 1}
    evaluated = evaluate(
        dag, point, 31, rational_embedding, with_jacobian=True
    )
    assert [evaluated.values[root] for root in (4, 6)] == [1, 30]

    full = solve_linearization(evaluated, [4, 6], 31)
    assert full.rank == 2 and not full.inconsistent_rows
    corrected = apply_scaled_correction(point, full.correction, 1, 31)
    for root in (4, 6):
        gradient = evaluated.gradients[root]
        linear_value = evaluated.values[root] + sum(
            value * full.correction.get(evaluated.variable_names[index], 0)
            for index, value in gradient.items()
        )
        assert linear_value % 31 == 0

    restricted = solve_linearization(
        evaluated, [4, 6], 31, allowed_variables={"x"}
    )
    assert restricted.inconsistent_rows
    restricted_cokernel = left_cokernel_basis(
        evaluated, [4, 6], 31, allowed_variables={"x"}
    )
    assert restricted_cokernel.rank == 1
    assert len(restricted_cokernel.functionals) == 1
    projected = apply_row_functionals(
        restricted_cokernel.functionals,
        [evaluated.values[4], evaluated.values[6]],
        31,
    )
    assert projected != [0]

    affine = solve_linearization(
        evaluated, [4], 31, free_values={"y": 4}
    )
    assert affine.correction["y"] == 4
    assert (
        evaluated.values[4]
        + 2 * affine.correction["x"]
        + affine.correction["y"]
    ) % 31 == 0
    forced = solve_linearization(
        evaluated, [4], 31, right_hand_side=[7], free_values={"y": 3}
    )
    assert (2 * forced.correction["x"] + forced.correction["y"]) % 31 == 7
    nonmonotone_pivots = solve_linearization(
        evaluated,
        [2, 4],
        31,
        right_hand_side=[3, 7],
    )
    assert nonmonotone_pivots.pivot_variables == ("y", "x")
    assert nonmonotone_pivots.correction == {"x": 2, "y": 3}
    gauged = solve_linearization(
        evaluated,
        [2, 4],
        31,
        right_hand_side=[3, 7],
        prescribed_values={"x": 2},
    )
    assert not gauged.inconsistent_rows
    assert gauged.correction == {"x": 2, "y": 3}
    assert interpolate_corrections(
        {"x": 2, "z": 4}, {"x": 5, "y": 1}, 2, 31
    ) == {"x": 8, "y": 2, "z": 27}

    polynomial = [7, 3, 0, 5]
    samples = [
        evaluate_univariate(polynomial, value, 31) for value in range(31)
    ]
    assert interpolate_consecutive_values(samples, 31, 7) == polynomial
    assert polynomial_gcd_modp([30, 0, 1], [30, 1], 31) == [30, 1]
    assert polynomial_gcd_modp([1, 1], [1, 0, 1], 31) == [1]
    assert corrected != point
    series = evaluate_truncated_series(
        dag,
        {"x": [1, 2], "y": [0, 3, 4]},
        31,
        rational_embedding,
        2,
    )
    assert series[4] == (1, 7, 8)  # x^2+y
    assert series[6] == (30, 1, 4)  # y-x
    print("SPARSE_CIRCUIT_MODP_PASS")


if __name__ == "__main__":
    main()
