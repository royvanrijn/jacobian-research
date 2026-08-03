#!/usr/bin/env python3
"""Replay the first Ferrers Hall-face radicals over characteristic zero.

These are regressions for Section 5 of
``extended-geometry/BINARY_GVC_ENVELOPE_CLOSURE.md``.  They are not used in
the unrestricted proof: shifted-ray endpoint rigidity supplies the general
horizontal-separation theorem directly.

The default octic gap-three check takes a few seconds.  ``--gap-four`` adds
the first degree-nine staircase; its eight exact msolve saturations can take
several minutes.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from collections.abc import Iterable

import sympy as sp

from verify_binary_quartic_quadruple_root_gvc import (
    moment,
    singular_expression,
)


def primitive(
    expression: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> sp.Expr:
    """Return the primitive rational polynomial underlying ``expression``."""

    return sp.Poly(expression, *variables, domain=sp.QQ).primitive()[1].as_expr()


def moment_equations(
    polynomial: dict[tuple[int, int], sp.Expr],
    operator: dict[tuple[int, int], sp.Expr],
    variables: tuple[sp.Symbol, ...],
    bound: int = 6,
) -> tuple[sp.Expr, ...]:
    """Return distinct primitive coefficients through the requested moment."""

    equations: list[sp.Expr] = []
    for order in range(1, bound + 1):
        for coefficient in moment(polynomial, operator, order).values():
            value = primitive(coefficient, variables)
            if value and value not in equations:
                equations.append(value)
    return tuple(equations)


def in_monomial_ideal(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    generators: tuple[sp.Expr, ...],
) -> bool:
    """Test membership in a monomial ideal term by term."""

    generator_exponents = tuple(
        sp.Poly(generator, *variables).monoms()[0] for generator in generators
    )
    for exponent, _coefficient in sp.Poly(expression, *variables).terms():
        if not any(
            all(left >= right for left, right in zip(exponent, candidate))
            for candidate in generator_exponents
        ):
            return False
    return True


def msolve_empty(
    variables: tuple[sp.Symbol, ...], equations: Iterable[sp.Expr]
) -> None:
    """Certify that one rational affine saturation has no complex point."""

    executable = shutil.which("msolve")
    if executable is None:
        raise RuntimeError("msolve is required for the Ferrers regression")

    source = (
        ",".join(map(str, variables))
        + "\n0\n"
        + ",\n".join(singular_expression(value) for value in equations)
        + "\n"
    )
    input_name = output_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".msolve", delete=False
        ) as handle:
            handle.write(source)
            input_name = handle.name
        output_name = input_name + ".out"
        subprocess.run(
            [executable, "-f", input_name, "-o", output_name, "-t", "8"],
            check=True,
            capture_output=True,
            text=True,
        )
        with open(output_name, encoding="utf-8") as handle:
            assert handle.read().strip() == "[-1]:"
    finally:
        for name in (input_name, output_name):
            if name and os.path.exists(name):
                os.unlink(name)


def verify_octic_gap_three() -> None:
    """Verify ``(A,S,T,BP,BQ,CQ)`` for the first octic staircase."""

    A, B, C, P, Q, T, S = sp.symbols("A B C P Q T S")
    variables = (A, B, C, P, Q, T, S)
    operator = {
        (3, 1): 1,
        (0, 7): A,
        (1, 5): B,
        (2, 3): C,
    }
    polynomial = {
        (0, 8): 1,
        (1, 6): P,
        (2, 4): Q,
        (3, 2): T,
        (4, 0): S,
    }
    equations = moment_equations(polynomial, operator, variables)
    assert len(equations) == 15

    linear = sp.solve(equations[0], T)
    assert linear == [-3360 * A - 60 * B * P - 4 * C * Q]

    reduced_variables = (A, B, C, P, Q, S)
    reduced_equations: list[sp.Expr] = []
    for equation in equations[1:]:
        value = primitive(
            sp.expand(equation.subs(T, linear[0])), reduced_variables
        )
        if value and value not in reduced_equations:
            reduced_equations.append(value)
    assert len(reduced_equations) == 14

    expected = (A, S, B * P, B * Q, C * Q)
    assert all(
        in_monomial_ideal(value, reduced_variables, expected)
        for value in reduced_equations
    )

    z = sp.Symbol("z")
    for generator in expected:
        msolve_empty(
            (*reduced_variables, z),
            (*reduced_equations, z * generator - 1),
        )

    # The eliminated equation puts T in the same radical.  The squarefree
    # ideal has the three displayed minimal vertex covers.
    covers = (
        (A, S, P, Q),
        (A, S, B, Q),
        (A, S, B, C),
    )
    assert all(
        all(
            in_monomial_ideal(generator, reduced_variables, cover)
            for generator in expected
        )
        for cover in covers
    )
    print("PASS octic gap three: (A,S,T,BP,BQ,CQ)")


def gap_four_data():
    """Construct the degree-nine/order-four slope-two face."""

    degree = 9
    order = multiplicity = gap = 4
    operator_base = (4, 0)
    polynomial_base = (0, 9)

    def normal_operator(exponent: tuple[int, int]) -> bool:
        x_order, y_order = exponent
        total = x_order + y_order
        if total == order:
            return x_order >= multiplicity
        if total < order + 1:
            return False
        excess = total - order
        return x_order < multiplicity or x_order > multiplicity + excess

    def weight(exponent: tuple[int, int]) -> int:
        return 2 * exponent[0] + exponent[1]

    operator_support = [operator_base]
    operator_support.extend(
        (x_order, y_order)
        for x_order in range(multiplicity)
        for y_order in range(weight(operator_base) + 1)
        if normal_operator((x_order, y_order))
        and weight((x_order, y_order)) == weight(operator_base)
    )
    polynomial_support = [polynomial_base]
    polynomial_support.extend(
        (x_degree, y_degree)
        for x_degree in range(degree + 1)
        for y_degree in range(degree + 1 - x_degree)
        if (x_degree, y_degree) != polynomial_base
        and weight((x_degree, y_degree)) == weight(polynomial_base)
    )

    operator_variables = sp.symbols(f"a1:{len(operator_support)}")
    polynomial_variables = sp.symbols(f"b1:{len(polynomial_support)}")
    variables = (*operator_variables, *polynomial_variables)
    operator = {
        operator_base: 1,
        **dict(zip(operator_support[1:], operator_variables, strict=True)),
    }
    polynomial = {
        polynomial_base: 1,
        **dict(
            zip(polynomial_support[1:], polynomial_variables, strict=True)
        ),
    }
    equations = moment_equations(polynomial, operator, variables)

    operator_by_deficit = {
        multiplicity - exponent[0]: variable
        for exponent, variable in zip(
            operator_support[1:], operator_variables, strict=True
        )
    }
    polynomial_by_advance = {
        exponent[0]: variable
        for exponent, variable in zip(
            polynomial_support[1:], polynomial_variables, strict=True
        )
    }
    expected: list[sp.Expr] = []
    expected.extend(
        variable
        for index, variable in operator_by_deficit.items()
        if index >= gap
    )
    expected.extend(
        variable
        for index, variable in polynomial_by_advance.items()
        if index >= gap
    )
    expected.extend(
        left * right
        for left_index, left in operator_by_deficit.items()
        for right_index, right in polynomial_by_advance.items()
        if left_index < gap
        and right_index < gap
        and left_index + right_index >= gap
    )
    return variables, equations, tuple(expected)


def verify_degree_nine_gap_four() -> None:
    """Verify the first eight-generator Ferrers staircase."""

    variables, equations, expected = gap_four_data()
    assert len(variables) == 8
    assert len(equations) == 15
    assert len(expected) == 8
    assert all(
        in_monomial_ideal(value, variables, expected) for value in equations
    )

    z = sp.Symbol("z")
    for generator in expected:
        msolve_empty((*variables, z), (*equations, z * generator - 1))
        print(f"PASS gap-four saturation: {generator}")
    print("PASS degree-nine gap four: Ferrers staircase")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gap-four",
        action="store_true",
        help="also run the several-minute degree-nine gap-four saturation",
    )
    args = parser.parse_args()
    verify_octic_gap_three()
    if args.gap_four:
        verify_degree_nine_gap_four()


if __name__ == "__main__":
    main()
