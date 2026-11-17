#!/usr/bin/env python3
"""Verify the general-position squarefree quartic synchronization closure.

For a projective frame of four lines this checker verifies:

* absence of quartic polar syzygies with zero or one tangent flag;
* the four-chart rank classification for two tangent and one transverse
  polar on three coordinate lines;
* the eight-chart rank classification for three tangent polars on three
  coordinate lines.

Together with the written polar-injectivity arguments, these calculations
close all sixteen flag patterns in the no-three-concurrent arrangement.
"""

from __future__ import annotations

import itertools

import sympy as sp


x, y, z, normal = sp.symbols("x y z normal")
variables = (x, y, z)
quintic_exponents = [
    (first, second, 5 - first - second)
    for first in range(6)
    for second in range(6 - first)
]
quintic_monomials = [
    x**first * y**second * z**third
    for first, second, third in quintic_exponents
]
coefficients = sp.symbols("c0:21")
h_five = sum(
    coefficient * monomial
    for coefficient, monomial in zip(coefficients, quintic_monomials)
)


def constraint_block(
    line: sp.Expr,
    direction: tuple[sp.Expr | int, ...],
    order: int,
) -> sp.Matrix:
    derivative = sp.expand(
        sum(
            sp.sympify(entry) * sp.diff(h_five, variable)
            for entry, variable in zip(direction, variables)
        )
    )
    if line == x:
        substitution, chart = {x: normal}, (normal, y, z)
    elif line == y:
        substitution, chart = {y: normal}, (x, normal, z)
    elif line == z:
        substitution, chart = {z: normal}, (x, y, normal)
    else:
        raise ValueError(f"unsupported coordinate line {line}")
    polynomial = sp.Poly(sp.expand(derivative.subs(substitution)), *chart)
    normal_index = chart.index(normal)
    return sp.Matrix(
        [
            [sp.expand(term_coefficient).coeff(coefficient) for coefficient in coefficients]
            for monomial, term_coefficient in polynomial.terms()
            if monomial[normal_index] < order
        ]
    )


def constraint_matrix(
    directions: tuple[tuple[sp.Expr | int, ...], ...],
    orders: tuple[int, int, int],
) -> sp.Matrix:
    return sp.Matrix.vstack(
        *(
            constraint_block(line, direction, order)
            for line, direction, order in zip(variables, directions, orders)
        )
    )


def pivot_minor(
    matrix: sp.Matrix,
    sample: dict[sp.Symbol, int],
    size: int,
) -> sp.Expr:
    numeric = matrix.subs(sample)
    assert numeric.rank() >= size
    _, pivot_rows = numeric.T.rref()
    _, pivot_columns = numeric.rref()
    rows = list(pivot_rows[:size])
    columns = list(pivot_columns[:size])
    return sp.factor(matrix.extract(rows, columns).det(method="domain-ge"))


def assert_associate(actual: sp.Expr, expected: sp.Expr) -> None:
    quotient = sp.cancel(actual / expected)
    assert sp.denom(quotient) == 1
    assert quotient != 0
    assert not sp.sympify(quotient).free_symbols


def nullspace_forms(matrix: sp.Matrix) -> list[sp.Expr]:
    return [
        sp.factor(
            sum(entry * monomial for entry, monomial in zip(vector, quintic_monomials))
        )
        for vector in matrix.nullspace()
    ]


def quartic_polar_syzygy_dimension(
    lines: tuple[sp.Expr, ...], pattern: str
) -> int:
    quartic_monomials = [
        x**first * y**second * z ** (4 - first - second)
        for first in range(5)
        for second in range(5 - first)
    ]
    columns: list[sp.Poly] = []
    for line, flag in zip(lines, pattern):
        exponent = 2 if flag == "T" else 3
        residual_degree = 4 - exponent
        residuals = [
            x**first * y**second * z ** (residual_degree - first - second)
            for first in range(residual_degree + 1)
            for second in range(residual_degree + 1 - first)
        ]
        columns.extend(
            sp.Poly(sp.expand(line**exponent * residual), x, y, z)
            for residual in residuals
        )
    matrix = sp.Matrix(
        [
            [column.coeff_monomial(monomial) for column in columns]
            for monomial in quartic_monomials
        ]
    )
    return len(matrix.nullspace())


def verify_zero_and_one_tangent() -> None:
    lines = (x, y, z, x + y + z)
    assert quartic_polar_syzygy_dimension(lines, "RRRR") == 0
    assert quartic_polar_syzygy_dimension(lines, "TRRR") == 0


def verify_two_tangent_coordinate_lemma() -> None:
    """Classify D_v h in (x^2),(y^2),(z^3) in four direction charts."""
    a, b, c, d = sp.symbols("a b c d")

    # Both tangent directions in their finite charts.
    main = constraint_matrix(
        ((0, a, 1), (b, 0, 1), (c, d, 1)), (2, 2, 3)
    )
    direction_determinant = a * b - a * c - b * d
    assert_associate(
        pivot_minor(main, {a: 2, b: 3, c: 5, d: 7}, 21),
        a**7 * b**3 * d**3 * direction_determinant,
    )
    main_a_zero = main.subs(a, 0)
    assert_associate(
        pivot_minor(main_a_zero, {b: 2, c: 3, d: 5}, 21), b**7 * d**5
    )
    main_b_zero = main.subs(b, 0)
    assert_associate(
        pivot_minor(main_b_zero, {a: 2, c: 3, d: 5}, 21),
        a**7 * c**4 * d,
    )
    main_b_d_zero = main.subs({b: 0, d: 0})
    assert_associate(
        pivot_minor(main_b_d_zero, {a: 2, c: 3}, 21), a**7 * c**5
    )
    main_d_zero = main.subs(d, 0)
    assert_associate(
        pivot_minor(main_d_zero, {a: 2, b: 3, c: 5}, 21),
        a**7 * b**3 * c**3 * (a * b - a * c),
    )
    binary_component = main.subs({c: 0, d: 0})
    assert nullspace_forms(binary_component) == [x**2 * y**3, x**3 * y**2]
    assert_associate(
        pivot_minor(binary_component, {a: 2, b: 3}, 19), a**8 * b**4
    )

    # Boundary tangent direction on x=0.
    x_boundary = constraint_matrix(
        ((0, 1, 0), (b, 0, 1), (c, d, 1)), (2, 2, 3)
    )
    assert_associate(
        pivot_minor(x_boundary, {b: 2, c: 3, d: 5}, 21),
        b**3 * d**3 * (b - c),
    )
    assert_associate(
        pivot_minor(x_boundary.subs(b, 0), {c: 3, d: 5}, 21), c**4 * d
    )
    assert_associate(
        pivot_minor(x_boundary.subs({b: 0, d: 0}), {c: 3}, 21), c**5
    )
    assert_associate(
        pivot_minor(x_boundary.subs(d, 0), {b: 2, c: 3}, 21),
        b**3 * c**3 * (b - c),
    )
    x_boundary_binary = x_boundary.subs({c: 0, d: 0})
    assert nullspace_forms(x_boundary_binary) == [x**2 * y**3, x**3 * y**2]
    assert_associate(pivot_minor(x_boundary_binary, {b: 3}, 19), b**4)

    # Boundary tangent direction on y=0.
    y_boundary = constraint_matrix(
        ((0, a, 1), (1, 0, 0), (c, d, 1)), (2, 2, 3)
    )
    assert_associate(
        pivot_minor(y_boundary, {a: 2, c: 3, d: 5}, 21),
        a**3 * d**3 * (a - d),
    )
    assert_associate(
        pivot_minor(y_boundary.subs(a, 0), {c: 3, d: 5}, 21), d**5
    )
    assert_associate(
        pivot_minor(y_boundary.subs(d, 0), {a: 2, c: 3}, 21), a**4 * c**3
    )
    y_boundary_binary = y_boundary.subs({c: 0, d: 0})
    assert nullspace_forms(y_boundary_binary) == [x**2 * y**3, x**3 * y**2]
    assert_associate(pivot_minor(y_boundary_binary, {a: 3}, 19), a**4)

    # Both tangent directions at infinity.  The whole solution is z^5.
    both_boundary = constraint_matrix(
        ((0, 1, 0), (1, 0, 0), (c, d, 1)), (2, 2, 3)
    )
    assert nullspace_forms(both_boundary) == [z**5]
    assert_associate(pivot_minor(both_boundary, {c: 2, d: 3}, 20), d**3)
    assert_associate(
        pivot_minor(both_boundary.subs(d, 0), {c: 2}, 20), c**3
    )
    both_zero = both_boundary.subs({c: 0, d: 0})
    assert nullspace_forms(both_zero) == [z**5, x**2 * y**3, x**3 * y**2]
    assert pivot_minor(both_zero, {}, 18) != 0

    # The only quartic syzygy among x^2*S2, y^2*S2, z^3*S1 is the repeated
    # tangent-pair relation (-x^2*y^2, x^2*y^2, 0).
    assert quartic_polar_syzygy_dimension((x, y, z), "TTR") == 1


def verify_three_tangent_coordinate_lemma() -> None:
    """Classify D_v h in (x^2),(y^2),(z^2) on all eight P1 charts."""
    a, b, c = sp.symbols("a b c")
    direction_charts = (
        ((0, a, 1), (0, 1, 0)),
        ((b, 0, 1), (1, 0, 0)),
        ((c, 1, 0), (1, 0, 0)),
    )

    main = constraint_matrix(
        tuple(chart[0] for chart in direction_charts), (2, 2, 2)
    )
    direction_determinant = a * c + b
    assert_associate(
        pivot_minor(main, {a: 2, b: 3, c: 5}, 21),
        a**7 * b**3 * direction_determinant,
    )
    assert_associate(
        pivot_minor(main.subs(a, 0), {b: 3, c: 5}, 21), b**7
    )
    assert_associate(
        pivot_minor(main.subs(b, 0), {a: 2, c: 3}, 21), a**7 * c**4
    )
    dependent_main = main.subs(b, -a * c)
    assert nullspace_forms(dependent_main) == [(a * c * z - c * y + x) ** 5]
    assert pivot_minor(dependent_main, {a: 0, c: 0}, 20) != 0

    # The other seven projective charts.
    matrices: dict[tuple[int, int, int], sp.Matrix] = {}
    for bits in itertools.product((0, 1), repeat=3):
        if bits == (0, 0, 0):
            continue
        directions = tuple(
            direction_charts[index][bit] for index, bit in enumerate(bits)
        )
        matrices[bits] = constraint_matrix(directions, (2, 2, 2))

    chart_001 = matrices[(0, 0, 1)]
    assert_associate(
        pivot_minor(chart_001, {a: 2, b: 3}, 21), a**8 * b**3
    )
    assert_associate(pivot_minor(chart_001.subs(b, 0), {a: 2}, 21), a**7)
    dependent_001 = chart_001.subs(a, 0)
    assert nullspace_forms(dependent_001) == [y**5]
    assert_associate(pivot_minor(dependent_001, {b: 2}, 20), b**7)
    assert pivot_minor(dependent_001.subs(b, 0), {}, 20) != 0

    chart_010 = matrices[(0, 1, 0)]
    assert_associate(pivot_minor(chart_010, {a: 2, c: 3}, 21), a**3)
    assert pivot_minor(chart_010.subs(a, 0), {c: 2}, 21) != 0

    chart_011 = matrices[(0, 1, 1)]
    assert nullspace_forms(chart_011) == [(-a * z + y) ** 5]
    assert_associate(pivot_minor(chart_011, {a: 2}, 20), a**3)
    assert pivot_minor(chart_011.subs(a, 0), {}, 20) != 0

    chart_100 = matrices[(1, 0, 0)]
    assert_associate(pivot_minor(chart_100, {b: 3, c: 2}, 21), b**3 * c)
    assert_associate(pivot_minor(chart_100.subs(b, 0), {c: 2}, 21), c**4)
    dependent_100 = chart_100.subs(c, 0)
    assert nullspace_forms(dependent_100) == [(-b * z + x) ** 5]
    assert_associate(pivot_minor(dependent_100, {b: 2}, 20), b**3)
    assert pivot_minor(dependent_100.subs(b, 0), {}, 20) != 0

    chart_101 = matrices[(1, 0, 1)]
    assert_associate(pivot_minor(chart_101, {b: 2}, 21), b**3)
    assert pivot_minor(chart_101.subs(b, 0), {}, 21) != 0

    chart_110 = matrices[(1, 1, 0)]
    assert nullspace_forms(chart_110) == [z**5]
    assert pivot_minor(chart_110, {c: 2}, 20) != 0

    chart_111 = matrices[(1, 1, 1)]
    assert nullspace_forms(chart_111) == [z**5]
    assert pivot_minor(chart_111, {}, 20) != 0


def main() -> None:
    verify_zero_and_one_tangent()
    print("PASS: no general-position polar syzygy with zero or one tangent")
    verify_two_tangent_coordinate_lemma()
    print("PASS: two-tangent coordinate-line atlas has only cones")
    verify_three_tangent_coordinate_lemma()
    print("PASS: three-tangent coordinate-line atlas has only fifth powers")
    print("PASS: squarefree quartic general-position closure")


if __name__ == "__main__":
    main()
