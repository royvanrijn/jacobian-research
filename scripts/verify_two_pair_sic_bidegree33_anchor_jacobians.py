#!/usr/bin/env python3
"""Exact Jacobian certificates on the five full non-null quadratic charts.

After normalizing the Sym^2 component to 2*X*T and quotienting its residual
torus on a nonzero-weight coefficient chart, moments 2 through 12 form an
eleven-by-eleven system.  This checker evaluates its Jacobian exactly at
one rational point on each of the five Weyl-orbit representative charts.
"""

from __future__ import annotations

from math import factorial
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_anchor_jacobians.json"
)
DEGREE = 3
MONOMIALS = tuple(
    (i, j) for i in range(DEGREE + 1) for j in range(DEGREE + 1)
)
PARAMETERS = (
    "s0", "s1", "s2", "s3", "s4", "s5", "s6",
    "t0", "t1", "t2", "t3", "t4",
)
SEXTIC_MAP = (
    (0, 0, 3, -1), (0, 1, 4, -3), (0, 2, 5, -3), (0, 3, 6, -1),
    (1, 0, 2, 3), (1, 1, 3, 9), (1, 2, 4, 9), (1, 3, 5, 3),
    (2, 0, 1, -3), (2, 1, 2, -9), (2, 2, 3, -9), (2, 3, 4, -3),
    (3, 0, 0, 1), (3, 1, 1, 3), (3, 2, 2, 3), (3, 3, 3, 1),
)
QUARTIC_MAP_LOCAL = (
    (0, 0, 2, 1), (0, 1, 3, 2), (0, 2, 4, 1),
    (1, 0, 1, -2), (1, 1, 2, -3), (1, 3, 4, 1),
    (2, 0, 0, 1), (2, 2, 2, -3), (2, 3, 3, -2),
    (3, 1, 0, 1), (3, 2, 1, 2), (3, 3, 2, 1),
)
QUARTIC_MAP = tuple(
    (i, j, parameter + 7, coefficient)
    for i, j, parameter, coefficient in QUARTIC_MAP_LOCAL
)
COEFFICIENT_MAP = SEXTIC_MAP + QUARTIC_MAP
NORMALIZED_QUADRATIC = (
    (0, 0, -1), (1, 1, -1), (2, 2, 1), (3, 3, 1),
)
MOMENT_ORDERS = tuple(range(2, 13))
CHART_POINTS = {
    0: (1, 0, 2, 0, 0, 2, 3, -2, -1, 3, -2, -2),
    1: (2, 1, -2, -1, 1, -2, -3, 2, 1, 0, 1, 1),
    2: (-3, -1, 1, -2, 1, -2, -3, 0, -2, -1, 1, 3),
    7: (2, 3, -1, -3, 0, 0, -1, 1, -1, -3, 0, -2),
    8: (3, -2, 3, 1, -1, 1, -2, 1, 1, 1, -3, 0),
}


def multiply_dense_bivariate(
    left: list[list[int]], right: tuple[tuple[int, ...], ...]
) -> list[list[int]]:
    answer = [
        [0] * (len(left[0]) + len(right[0]) - 1)
        for _ in range(len(left) + len(right) - 1)
    ]
    for left_i, row in enumerate(left):
        for left_j, left_coefficient in enumerate(row):
            if left_coefficient == 0:
                continue
            for right_i, right_row in enumerate(right):
                for right_j, right_coefficient in enumerate(right_row):
                    answer[left_i + right_i][left_j + right_j] += (
                        left_coefficient * right_coefficient
                    )
    return answer


def coefficient_grid(point: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    grid = [[0] * (DEGREE + 1) for _ in range(DEGREE + 1)]
    for i, j, parameter, coefficient in COEFFICIENT_MAP:
        grid[i][j] += coefficient * point[parameter]
    for i, j, coefficient in NORMALIZED_QUADRATIC:
        grid[i][j] += coefficient
    return tuple(tuple(row) for row in grid)


def chart_jacobian(
    fixed_parameter: int,
    point: tuple[int, ...],
) -> sp.Matrix:
    assert point[fixed_parameter] == 1
    variables = tuple(
        parameter
        for parameter in range(len(PARAMETERS))
        if parameter != fixed_parameter
    )
    grid = coefficient_grid(point)
    power: list[list[int]] = [[1]]
    rows: list[list[int]] = []

    for order in MOMENT_ORDERS:
        power = multiply_dense_bivariate(power, grid)
        coefficient_gradient: list[int] = []
        for coefficient_i, coefficient_j in MONOMIALS:
            value = 0
            for diagonal in range(DEGREE * order + 1):
                power_i = diagonal - coefficient_i
                power_j = diagonal - coefficient_j
                if (
                    0 <= power_i < len(power)
                    and 0 <= power_j < len(power[0])
                ):
                    value += (
                        factorial(DEGREE * order - diagonal)
                        * factorial(diagonal)
                        * power[power_i][power_j]
                    )
            coefficient_gradient.append(order * value)

        rows.append([
            sum(
                coefficient
                * coefficient_gradient[MONOMIALS.index((i, j))]
                for i, j, parameter_index, coefficient in COEFFICIENT_MAP
                if parameter_index == parameter
            )
            for parameter in variables
        ])

    return sp.Matrix(rows)


def main() -> None:
    certificates: list[dict[str, object]] = []
    for fixed_parameter, point in CHART_POINTS.items():
        jacobian = chart_jacobian(fixed_parameter, point)
        determinant = int(jacobian.det(method="domain-ge"))
        assert determinant != 0
        assert jacobian.rank() == 11
        certificates.append({
            "fixed_parameter": PARAMETERS[fixed_parameter],
            "point": list(point),
            "variable_order": [
                parameter
                for index, parameter in enumerate(PARAMETERS)
                if index != fixed_parameter
            ],
            "moment_orders": list(MOMENT_ORDERS),
            "determinant": str(determinant),
            "rank": 11,
        })
        print(
            f"PASS full anchor chart {PARAMETERS[fixed_parameter]}=1: "
            "moments 2..12 have exact Jacobian rank 11"
        )

    artifact = {
        "format": "two-pair-sic-bidegree33-anchor-jacobians-v1",
        "field": "characteristic zero",
        "normalized_quadratic": "2*X*T",
        "higher_components": "full Sym^6+Sym^4",
        "representative_charts": [
            "s0=1", "s1=1", "s2=1", "t0=1", "t1=1"
        ],
        "weyl_reflections": [
            "s6=1", "s5=1", "s4=1", "t4=1", "t3=1"
        ],
        "certificates": certificates,
        "conclusion": (
            "on every nonzero-weight chart orbit, moments 2 through 12 "
            "are algebraically independent; this does not determine their "
            "common zero fiber"
        ),
        "reproduce": (
            ".venv/bin/python "
            "scripts/verify_two_pair_sic_bidegree33_anchor_jacobians.py"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
