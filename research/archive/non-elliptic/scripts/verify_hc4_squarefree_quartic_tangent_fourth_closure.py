#!/usr/bin/env python3
"""Verify the last eight clean squarefree quartic-denominator rows.

The four representative flag patterns are RRRT, TRRT, TTRT, and TTTT on

    (x, y, x+y, z).

For a polar relation write q_i=D_{w_i}h.  Exact Singular saturation proves
that the mixed-partial equations have no component on which both the q_i
and the w_i have rank at least two.  A noncone quintic has injective polar
map, so the two ranks must agree.  The remaining rank-one relations are
repeated tangent pairs; the final SymPy checks replay the elementary
blocking-polar calculation.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp


x, y, z = sp.symbols("x y z")
variables = (x, y, z)
lines = (x, y, x + y, z)
quartic_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]


def residual_monomials(degree: int) -> list[sp.Expr]:
    return [
        x**i * y**j * z ** (degree - i - j)
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
    ]


def polar_syzygies(pattern: str) -> list[tuple[sp.Expr, ...]]:
    """Compute sum(q_i)=0 with q_i in (L_i^2) or (L_i^3)."""
    columns: list[sp.Poly] = []
    blocks: list[tuple[sp.Expr, list[sp.Expr]]] = []
    for line, flag in zip(lines, pattern):
        exponent = 2 if flag == "T" else 3
        residuals = residual_monomials(4 - exponent)
        blocks.append((line**exponent, residuals))
        columns.extend(
            sp.Poly(sp.expand(line**exponent * residual), *variables)
            for residual in residuals
        )
    matrix = sp.Matrix(
        [
            [column.coeff_monomial(monomial) for column in columns]
            for monomial in quartic_monomials
        ]
    )
    result: list[tuple[sp.Expr, ...]] = []
    for vector in matrix.nullspace():
        offset = 0
        components: list[sp.Expr] = []
        for line_power, residuals in blocks:
            components.append(
                sp.factor(
                    line_power
                    * sum(
                        vector[offset + index] * residual
                        for index, residual in enumerate(residuals)
                    )
                )
            )
            offset += len(residuals)
        result.append(tuple(components))
    return result


def coefficient_vector(form: sp.Expr) -> tuple[sp.Expr, ...]:
    polynomial = sp.Poly(sp.expand(form), *variables)
    return tuple(
        polynomial.coeff_monomial(monomial) for monomial in quartic_monomials
    )


def directional_derivative(
    direction: tuple[sp.Expr, ...], form: sp.Expr
) -> sp.Expr:
    return sp.expand(
        sum(entry * sp.diff(form, variable) for entry, variable in zip(direction, variables))
    )


def polynomial_coefficients(form: sp.Expr) -> list[sp.Expr]:
    return [coefficient for _, coefficient in sp.Poly(sp.expand(form), *variables).terms()]


def incidence_data(pattern: str):
    basis = polar_syzygies(pattern)
    parameters = sp.symbols(f"s0:{len(basis)}")
    syzygy = [
        sp.expand(
            sum(parameter * basis_vector[index] for parameter, basis_vector in zip(parameters, basis))
        )
        for index in range(4)
    ]

    direction_variables = sp.symbols("a0:9")
    w1 = tuple(direction_variables[0:3])
    w2 = tuple(direction_variables[3:6])
    w3 = tuple(direction_variables[6:9])
    w4 = tuple(-w1[index] - w2[index] - w3[index] for index in range(3))
    directions = (w1, w2, w3, w4)

    equations: list[sp.Expr] = []
    for first in range(4):
        for second in range(first + 1, 4):
            equations.extend(
                polynomial_coefficients(
                    directional_derivative(directions[first], syzygy[second])
                    - directional_derivative(directions[second], syzygy[first])
                )
            )
    for index, flag in enumerate(pattern):
        if flag == "T":
            equations.append(
                sp.expand(lines[index].subs(dict(zip(variables, directions[index]))))
            )
    equations = [equation for equation in equations if equation != 0]
    return parameters, direction_variables, syzygy, directions, equations


def two_by_two_minors(columns: list[tuple[sp.Expr, ...]]) -> list[sp.Expr]:
    minors: list[sp.Expr] = []
    seen: set[str] = set()
    row_count = len(columns[0])
    for first_row in range(row_count):
        for second_row in range(first_row + 1, row_count):
            for first_column in range(len(columns)):
                for second_column in range(first_column + 1, len(columns)):
                    minor = sp.factor(
                        columns[first_column][first_row]
                        * columns[second_column][second_row]
                        - columns[first_column][second_row]
                        * columns[second_column][first_row]
                    )
                    if minor == 0:
                        continue
                    key = str(minor)
                    negative_key = str(-minor)
                    if key in seen or negative_key in seen:
                        continue
                    seen.add(key)
                    minors.append(minor)
    return minors


def singular_polynomial(
    expression: sp.Expr, ring_variables: tuple[sp.Symbol, ...]
) -> str:
    _, polynomial = sp.Poly(
        expression, *ring_variables, domain=sp.QQ
    ).clear_denoms(convert=True)
    return str(sp.expand(polynomial.as_expr())).replace("**", "^")


def verify_rank_at_least_two_empty(pattern: str) -> None:
    parameters, direction_variables, syzygy, directions, equations = incidence_data(pattern)
    ring_variables = parameters + direction_variables
    q_columns = [coefficient_vector(component) for component in syzygy]
    q_minors = two_by_two_minors(q_columns)
    w_minors = two_by_two_minors([tuple(direction) for direction in directions])
    assert q_minors and w_minors

    def ideal(expressions: list[sp.Expr]) -> str:
        return ",".join(
            singular_polynomial(expression, ring_variables) for expression in expressions
        )

    source = "\n".join(
        [
            'LIB "elim.lib";',
            f"ring r=0,({','.join(map(str, ring_variables))}),dp;",
            "option(redSB);",
            f"ideal I={ideal(equations)};",
            f"ideal Q={ideal(q_minors)};",
            f"ideal W={ideal(w_minors)};",
            "ideal IQ=sat(I,Q);",
            "ideal IQW=sat(IQ,W);",
            'print("RANK_TWO_SATURATION");',
            "IQW;",
            "exit;",
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"],
        input=source,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "RANK_TWO_SATURATION\nIQW[1]=1" in completed.stdout


def coefficient_ideal(form: sp.Expr, generators: tuple[sp.Symbol, ...]) -> sp.GroebnerBasis:
    coefficients = [
        coefficient for _, coefficient in sp.Poly(sp.expand(form), x, y, z).terms()
    ]
    return sp.groebner(coefficients, *generators, order="grevlex")


def verify_rank_one_blockers() -> None:
    """Replay the two normal forms for every repeated tangent pair."""
    a, b, c = sp.symbols("a b c")

    # Two lines in the concurrent pencil repeat the common direction d/dz.
    concurrent_core = x**2 * y**2 * z
    assert sp.diff(concurrent_core, z) == x**2 * y**2
    z_coefficient = sp.diff(
        directional_derivative((a, b, 0), concurrent_core), z
    )
    assert coefficient_ideal(z_coefficient, (a, b)).contains(a)
    assert coefficient_ideal(z_coefficient, (a, b)).contains(b)

    # If the third pencil line is transverse, divisibility by (x+y)^3 of
    # the same z coefficient also forces a=b=0.  Substitute y=n-x.
    n = sp.symbols("n")
    normal_form = sp.expand(z_coefficient.subs(y, n - x))
    low_normal_coefficients = [
        sp.Poly(normal_form, x, n).coeff_monomial(x ** (3 - degree) * n**degree)
        for degree in range(3)
    ]
    transverse_basis = sp.groebner(low_normal_coefficients, a, b, order="grevlex")
    assert transverse_basis.contains(a)
    assert transverse_basis.contains(b)

    # A pencil line and z repeat their common tangent direction.  For the
    # normalized pair (x,z), integration gives x^2*y*z^2 plus F(x,z).
    cross_core = x**2 * y * z**2
    assert sp.diff(cross_core, y) == x**2 * z**2
    y_coefficient = sp.diff(
        directional_derivative((a, b, c), cross_core), y
    )
    blocker_basis = coefficient_ideal(y_coefficient, (a, c))
    assert blocker_basis.contains(a)
    assert blocker_basis.contains(c)
    assert directional_derivative((0, b, 0), cross_core).subs(y, 0) == b * x**2 * z**2


def main() -> None:
    assert shutil.which("Singular") is not None, "Singular is required"
    expected_dimensions = {"RRRT": 1, "TRRT": 3, "TTRT": 6, "TTTT": 9}
    for pattern, expected_dimension in expected_dimensions.items():
        assert len(polar_syzygies(pattern)) == expected_dimension
        verify_rank_at_least_two_empty(pattern)
        print(f"PASS: {pattern} has no equal-rank relation of rank at least two")
    verify_rank_one_blockers()
    print("PASS: every rank-one relation is blocked by an unused flag")
    print("PASS: squarefree quartic tangent-fourth closure")


if __name__ == "__main__":
    main()
