#!/usr/bin/env python3
"""Exact Jacobian certificates on the five full non-null quadratic charts.

After normalizing the Sym^2 component to 2*X*T and quotienting its residual
torus on a nonzero-weight coefficient chart, moments 2 through 12 form an
eleven-by-eleven system.  This checker evaluates its Jacobian exactly at
one rational point on each of the five Weyl-orbit representative charts.
It also proves that moment 2 eliminates the opposite-weight coordinate
with a constant nonzero pivot on every chart and records the next exact
two-open-plus-boundary split of moment 3 on the s0 chart.  Finally it
closes a natural two-parameter plane in the common boundary by an exact
moment-3/moment-4 unit certificate.
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
OPPOSITE_PIVOTS = {
    0: 6,
    1: 5,
    2: 4,
    7: 11,
    8: 10,
}
SECOND_MOMENT_PIVOT_DERIVATIVES = {
    0: -72,
    1: 432,
    2: -1080,
    7: 336,
    8: -1344,
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


def symbolic_moment(order: int) -> tuple[tuple[sp.Symbol, ...], sp.Expr]:
    """Construct one exact normalized full-anchor moment."""

    parameters = sp.symbols(" ".join(PARAMETERS))
    x, y = sp.symbols("x y")
    grid: list[list[sp.Expr]] = [
        [sp.Integer(0)] * (DEGREE + 1)
        for _ in range(DEGREE + 1)
    ]
    for i, j, parameter, coefficient in COEFFICIENT_MAP:
        grid[i][j] += coefficient * parameters[parameter]
    for i, j, coefficient in NORMALIZED_QUADRATIC:
        grid[i][j] += coefficient
    polynomial = sum(
        grid[i][j] * x**i * y**j
        for i, j in MONOMIALS
    )
    power = sp.Poly(sp.expand(polynomial**order), x, y)
    moment = sp.expand(sum(
        factorial(DEGREE * order - diagonal)
        * factorial(diagonal)
        * power.coeff_monomial(x**diagonal * y**diagonal)
        for diagonal in range(DEGREE * order + 1)
    ))
    return parameters, moment


def check_exact_pivot_reduction() -> dict[str, object]:
    """Eliminate one variable exactly with mu_2 on every anchor chart."""

    parameters, second_moment = symbolic_moment(2)
    expected = -24 * (
        3 * parameters[0] * parameters[6]
        - 18 * parameters[1] * parameters[5]
        + 45 * parameters[2] * parameters[4]
        - 30 * parameters[3] ** 2
        - 14 * parameters[7] * parameters[11]
        + 56 * parameters[8] * parameters[10]
        - 42 * parameters[9] ** 2
        - 70
    )
    assert sp.expand(second_moment - expected) == 0

    pivots: list[dict[str, object]] = []
    for fixed_parameter, pivot_parameter in OPPOSITE_PIVOTS.items():
        chart_moment = sp.expand(
            second_moment.subs(parameters[fixed_parameter], 1)
        )
        derivative = sp.diff(chart_moment, parameters[pivot_parameter])
        assert derivative == SECOND_MOMENT_PIVOT_DERIVATIVES[fixed_parameter]
        assert sp.diff(
            chart_moment,
            parameters[pivot_parameter],
            2,
        ) == 0
        pivots.append({
            "chart": f"{PARAMETERS[fixed_parameter]}=1",
            "eliminated_variable": PARAMETERS[pivot_parameter],
            "constant_derivative": int(derivative),
        })

    _, third_moment = symbolic_moment(3)
    s0_chart_second = sp.expand(second_moment.subs(parameters[0], 1))
    s6_coefficient = sp.diff(s0_chart_second, parameters[6])
    s6_value = sp.expand(
        -(
            s0_chart_second
            - s6_coefficient * parameters[6]
        )
        / s6_coefficient
    )
    reduced_third = sp.expand(
        third_moment.subs(parameters[0], 1).subs(parameters[6], s6_value)
    )
    s5_coefficient = sp.factor(sp.diff(reduced_third, parameters[5]))
    t4_coefficient = sp.factor(sp.diff(reduced_third, parameters[11]))
    expected_s5_coefficient = -103680 * (
        6 * parameters[1] ** 2 * parameters[8]
        - 3 * parameters[1] * parameters[2] * parameters[7]
        - 3 * parameters[1] * parameters[9]
        - 3 * parameters[2] * parameters[8]
        + 2 * parameters[3] * parameters[7]
        - 3 * parameters[7]
        + parameters[10]
    )
    expected_t4_coefficient = -17280 * (
        12 * parameters[1] * parameters[3]
        + 28 * parameters[1] * parameters[7] * parameters[8]
        - 18 * parameters[1]
        - 9 * parameters[2] ** 2
        - 14 * parameters[2] * parameters[7] ** 2
        - 3 * parameters[4]
        - 2 * parameters[7] * parameters[9]
        - 12 * parameters[8] ** 2
    )
    assert sp.expand(s5_coefficient - expected_s5_coefficient) == 0
    assert sp.expand(t4_coefficient - expected_t4_coefficient) == 0
    assert sp.diff(reduced_third, parameters[5], 2) == 0
    assert sp.diff(reduced_third, parameters[11], 2) == 0
    assert sp.diff(
        reduced_third,
        parameters[5],
        parameters[11],
    ) == 0

    return {
        "normalized_second_moment": str(expected),
        "chart_pivots": pivots,
        "conclusion": (
            "mu_2 eliminates one opposite-weight variable with a "
            "nonzero constant pivot on every full anchor chart"
        ),
        "s0_chart_next_split": {
            "after_eliminating": "s6 with mu_2",
            "mu_3_is_affine_in": ["s5", "t4"],
            "s5_coefficient": str(expected_s5_coefficient),
            "t4_coefficient": str(expected_t4_coefficient),
            "next_branches": (
                "s5 coefficient nonzero; then t4 coefficient nonzero; "
                "then both coefficients zero"
            ),
        },
    }


def check_s0_sparse_boundary_slice() -> dict[str, object]:
    """Close a natural two-parameter slice of the A=B=0 boundary."""

    a, b, x, y = sp.symbols("a b x y")
    values: dict[str, sp.Expr] = {
        parameter: sp.Integer(0) for parameter in PARAMETERS
    }
    values.update({
        "s0": sp.Integer(1),
        "s6": (14 * a * b + 70) / 3,
        "t0": a,
        "t3": 3 * a,
        "t4": b,
    })
    grid: list[list[sp.Expr]] = [
        [sp.Integer(0)] * (DEGREE + 1)
        for _ in range(DEGREE + 1)
    ]
    for i, j, parameter, coefficient in COEFFICIENT_MAP:
        grid[i][j] += coefficient * values[PARAMETERS[parameter]]
    for i, j, coefficient in NORMALIZED_QUADRATIC:
        grid[i][j] += coefficient
    polynomial = sum(
        grid[i][j] * x**i * y**j
        for i, j in MONOMIALS
    )

    moments: dict[int, sp.Expr] = {}
    for order in (2, 3, 4):
        power = sp.Poly(sp.expand(polynomial**order), x, y)
        moments[order] = sp.factor(sum(
            factorial(DEGREE * order - diagonal)
            * factorial(diagonal)
            * power.coeff_monomial(x**diagonal * y**diagonal)
            for diagonal in range(DEGREE * order + 1)
        ))
    assert moments[2] == 0
    assert moments[3] == 1866240 * a**3
    quartic_factor = 11249 - 8776 * a * b - 901 * a**2 * b**2
    assert sp.expand(moments[4] - 138240 * quartic_factor) == 0, (
        moments[4],
        quartic_factor,
    )

    constant = sp.Integer(11249)
    u = a * b
    inverse_mod_a3 = (
        1 / constant
        + sp.Rational(8776, constant**2) * u
        + sp.Rational(
            8776**2 + 901 * constant,
            constant**3,
        )
        * u**2
    )
    a3_coefficient = sp.cancel(
        (1 - inverse_mod_a3 * quartic_factor) / a**3
    )
    assert sp.expand(
        inverse_mod_a3 * moments[4] / 138240
        + a3_coefficient * moments[3] / 1866240
        - 1
    ) == 0

    return {
        "slice": {
            "s0": "1",
            "s1,s2,s3,s4,s5,t1,t2": "0",
            "t0": "a",
            "t3": "3*a",
            "t4": "b",
            "s6": "(14*a*b+70)/3",
        },
        "boundary_conditions": "A=B=0 and mu_2=0 identically",
        "mu_3": str(moments[3]),
        "mu_4": str(moments[4]),
        "inverse_of_mu4_factor_mod_a3": str(inverse_mod_a3),
        "a3_certificate_coefficient": str(a3_coefficient),
        "conclusion": (
            "mu_3 and mu_4 generate the unit ideal on this "
            "two-parameter sparse boundary slice"
        ),
    }


def check_s0_stratum_jacobians() -> dict[str, object]:
    """Prove maximal differential rank on all three mu_3 pivot strata."""

    points: dict[str, tuple[sp.Expr, ...]] = {
        "A_nonzero": tuple(map(sp.Integer, CHART_POINTS[0])),
        "A_zero_B_nonzero": tuple(map(
            sp.Integer,
            (1, 0, 0, -3, 3, 2, 3, -2, 1, 2, -18, 2),
        )),
        "A_zero_B_zero": (
            sp.Integer(1),
            sp.Integer(0),
            sp.Integer(2),
            sp.Integer(0),
            sp.Rational(-148, 3),
            sp.Integer(2),
            sp.Integer(3),
            sp.Integer(-2),
            sp.Integer(-1),
            sp.Integer(3),
            sp.Integer(-12),
            sp.Integer(-2),
        ),
    }
    expected: dict[str, tuple[sp.Expr, sp.Expr, int]] = {
        "A_nonzero": (sp.Integer(10), sp.Integer(-148), 11),
        "A_zero_B_nonzero": (sp.Integer(0), sp.Integer(-13), 10),
        "A_zero_B_zero": (sp.Integer(0), sp.Integer(0), 9),
    }
    certificates: list[dict[str, object]] = []
    for stratum, point in points.items():
        (
            _,
            s1,
            s2,
            s3,
            s4,
            _,
            _,
            t0,
            t1,
            t2,
            t3,
            _,
        ) = point
        coefficient_a = (
            6 * s1**2 * t1
            - 3 * s1 * s2 * t0
            - 3 * s1 * t2
            - 3 * s2 * t1
            + 2 * s3 * t0
            - 3 * t0
            + t3
        )
        coefficient_b = (
            12 * s1 * s3
            + 28 * s1 * t0 * t1
            - 18 * s1
            - 9 * s2**2
            - 14 * s2 * t0**2
            - 3 * s4
            - 2 * t0 * t2
            - 12 * t1**2
        )
        expected_a, expected_b, restricted_rank = expected[stratum]
        assert coefficient_a == expected_a
        assert coefficient_b == expected_b
        jacobian = chart_jacobian(0, point)
        determinant = sp.factor(jacobian.det(method="domain-ge"))
        assert determinant != 0
        certificates.append({
            "stratum": stratum,
            "point": [str(coordinate) for coordinate in point],
            "A": str(coefficient_a),
            "B": str(coefficient_b),
            "full_chart_determinant": str(determinant),
            "full_chart_rank": 11,
            "restricted_stratum_rank": restricted_rank,
        })
    return {
        "certificates": certificates,
        "conclusion": (
            "moments 2 through 12 have maximal differential rank on "
            "each of the A-nonzero, A-zero B-nonzero, and A=B=0 strata"
        ),
    }


def main() -> None:
    pivot_reduction = check_exact_pivot_reduction()
    sparse_boundary_slice = check_s0_sparse_boundary_slice()
    stratum_jacobians = check_s0_stratum_jacobians()
    print(
        "PASS full anchor charts: mu_2 has a constant nonzero "
        "opposite-weight pivot on all five chart orbits"
    )
    print(
        "PASS s0 anchor chart: after eliminating s6, mu_3 gives two "
        "explicit principal-open pivots and their common boundary"
    )
    print(
        "PASS s0 sparse common-boundary slice: mu_3 and mu_4 "
        "generate the unit ideal"
    )
    print(
        "PASS s0 pivot strata: moments 2..12 have maximal restricted "
        "differential ranks 11,10,9"
    )
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
        "format": "two-pair-sic-bidegree33-anchor-jacobians-v3",
        "field": "characteristic zero",
        "normalized_quadratic": "2*X*T",
        "higher_components": "full Sym^6+Sym^4",
        "representative_charts": [
            "s0=1", "s1=1", "s2=1", "t0=1", "t1=1"
        ],
        "weyl_reflections": [
            "s6=1", "s5=1", "s4=1", "t4=1", "t3=1"
        ],
        "exact_pivot_reduction": pivot_reduction,
        "s0_sparse_boundary_slice": sparse_boundary_slice,
        "s0_stratum_jacobians": stratum_jacobians,
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
