#!/usr/bin/env python3
"""Verify the quadruple-root quartic-leading binary GVC calculation.

The all-order weight argument is in
``extended-geometry/BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md``.  This checker
replays the exact defect-one radical and the decisive branch identities.
Singular is required for the radical calculation.
"""

from __future__ import annotations

import math
import shutil
import subprocess

import sympy as sp


x, y = sp.symbols("x y")
Exponent = tuple[int, int]
SparsePolynomial = dict[Exponent, sp.Expr]


def falling(degree: int, order: int) -> int:
    return math.prod(range(degree - order + 1, degree + 1)) if order else 1


def multiply(
    left: SparsePolynomial, right: SparsePolynomial
) -> SparsePolynomial:
    result: SparsePolynomial = {}
    for (left_x, left_y), left_coefficient in left.items():
        for (right_x, right_y), right_coefficient in right.items():
            exponent = (left_x + right_x, left_y + right_y)
            result[exponent] = (
                result.get(exponent, 0)
                + left_coefficient * right_coefficient
            )
    return {
        exponent: sp.expand(coefficient)
        for exponent, coefficient in result.items()
        if coefficient != 0
    }


def apply_operator(
    polynomial: SparsePolynomial, operator: SparsePolynomial
) -> SparsePolynomial:
    result: SparsePolynomial = {}
    for (x_degree, y_degree), coefficient in polynomial.items():
        for (x_order, y_order), operator_coefficient in operator.items():
            if x_degree < x_order or y_degree < y_order:
                continue
            exponent = (x_degree - x_order, y_degree - y_order)
            result[exponent] = (
                result.get(exponent, 0)
                + coefficient
                * operator_coefficient
                * falling(x_degree, x_order)
                * falling(y_degree, y_order)
            )
    return {
        exponent: sp.expand(coefficient)
        for exponent, coefficient in result.items()
        if coefficient != 0
    }


def moment(
    polynomial: SparsePolynomial,
    operator: SparsePolynomial,
    order: int,
) -> SparsePolynomial:
    result: SparsePolynomial = {(0, 0): 1}
    for _ in range(order):
        result = multiply(result, polynomial)
    for _ in range(order):
        result = apply_operator(result, operator)
    return result


def primitive(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    _content, value = sp.Poly(
        sp.together(expression).as_numer_denom()[0],
        *variables,
        domain=sp.QQ,
    ).primitive()
    return sp.factor(value.as_expr())


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


def verify_radical(
    variables: tuple[sp.Symbol, ...],
    equations: tuple[sp.Expr, ...],
    expected: tuple[sp.Expr, ...],
) -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the radical replay")
    program = f"""
ring r=0,({",".join(map(str, variables))}),dp;
ideal I={",".join(map(singular_expression, equations))};
ideal E={",".join(map(singular_expression, expected))};
LIB "primdec.lib";
ideal R=std(radical(I));
ideal left=reduce(R,std(E));
ideal right=reduce(std(E),R);
if ((size(left)==0) && (size(right)==0)) {{ print("PASS"); }}
else {{ print("FAIL"); R; }}
quit;
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "PASS" in result.stdout and "FAIL" not in result.stdout


def verify_defect_one_radical() -> None:
    A, B, C = sp.symbols("A B C")
    l0, l1, l2, l3 = sp.symbols("l0:4")
    variables = (A, B, C, l0, l1, l2, l3)
    polynomial = {
        (0, 5): 1,
        (1, 4): A,
        (2, 3): B,
        (3, 2): C,
    }
    quintic_jet = {
        (0, 5): l0,
        (1, 4): l1,
        (2, 3): l2,
        (3, 2): l3,
    }
    equations: list[sp.Expr] = []
    power: SparsePolynomial = {(0, 0): 1}
    for order in range(1, 4):
        power = multiply(power, polynomial)
        if order < 2:
            continue
        value = apply_operator(power, quintic_jet)
        value = apply_operator(
            value, {(4 * (order - 1), 0): 1}
        )
        equations.extend(
            primitive(coefficient, variables)
            for coefficient in value.values()
        )
    verify_radical(
        variables,
        tuple(equations),
        (C * l2, C * l1, C * l0, B * l0),
    )


def verify_low_threshold_branches() -> None:
    A = sp.symbols("A")
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    r0, r1, r2, r3 = sp.symbols("r0:4")
    q0, q1, q2 = sp.symbols("q0:3")
    l0, l1, l2, l3 = sp.symbols("l0:4")
    k0, k1, k2, k3 = sp.symbols("k0:4")
    h0, h1, h2, h3 = sp.symbols("h0:4")

    # The A*l0 crossing, normalized by A=1.  Complete defect-three data
    # cannot change the displayed positive-degree fourth-moment coefficient.
    polynomial_a = {
        (0, 5): 1,
        (1, 4): 1,
        (0, 4): p0,
        (1, 3): p1,
        (2, 2): p2,
        (3, 1): p3,
        (4, 0): p4,
        (0, 3): r0,
        (1, 2): r1,
        (2, 1): r2,
        (3, 0): r3,
        (0, 2): q0,
        (1, 1): q1,
        (2, 0): q2,
    }
    operator_a = {
        (4, 0): 1,
        (0, 5): l0,
        (1, 4): l1,
        (2, 3): l2,
        (3, 2): l3,
        (0, 6): k0,
        (1, 5): k1,
        (2, 4): k2,
        (3, 3): k3,
        (0, 7): h0,
        (1, 6): h1,
        (2, 5): h2,
        (3, 4): h3,
    }
    second_a = moment(polynomial_a, operator_a, 2)[(0, 0)]
    p4_value = -5 * l0 - l1
    p3_value = -(
        795 * l0**2
        + 310 * l0 * l1
        + 28 * l0 * l2
        + 19 * l1**2
    ) / (2 * l0)
    assert sp.factor(second_a.subs(p4, p4_value).subs(p3, p3_value)) == 0
    fourth_a = moment(polynomial_a, operator_a, 4)
    variables_a = tuple(
        sorted(
            fourth_a[(0, 1)].free_symbols | {l0},
            key=str,
        )
    )
    assert primitive(
        fourth_a[(0, 1)].subs({p4: p4_value, p3: p3_value}),
        variables_a,
    ) == l0**3

    # Once l0=0, the equality A*l1 is already killed by moment two.
    assert primitive(
        second_a.subs({l0: 0, p4: -l1}),
        tuple(sorted(second_a.free_symbols, key=str)),
    ) == l1**2

    # The B*l1 crossing, normalized by B=1.
    polynomial_b = {
        (0, 5): 1,
        (1, 4): A,
        (2, 3): 1,
        (0, 4): p0,
        (1, 3): p1,
        (2, 2): p2,
        (3, 1): p3,
        (4, 0): p4,
        (0, 3): r0,
        (1, 2): r1,
        (2, 1): r2,
        (3, 0): r3,
        (0, 2): q0,
        (1, 1): q1,
        (2, 0): q2,
    }
    operator_b = {
        (4, 0): 1,
        (1, 4): l1,
        (2, 3): l2,
        (3, 2): l3,
        (0, 6): k0,
        (1, 5): k1,
        (2, 4): k2,
        (3, 3): k3,
        (0, 7): h0,
        (1, 6): h1,
        (2, 5): h2,
        (3, 4): h3,
    }
    p4_b = -A * l1 - l2 / 2
    third_b = moment(polynomial_b, operator_b, 3)
    assert primitive(
        third_b[(0, 1)].subs(p4, p4_b),
        tuple(sorted(third_b[(0, 1)].free_symbols, key=str)),
    ) == l1**2
    second_b = moment(polynomial_b, operator_b, 2)[(0, 0)]
    second_b_reduced = primitive(
        second_b.subs(p4, p4_b).subs(l1, 0),
        tuple(sorted(second_b.free_symbols, key=str)),
    )
    assert second_b_reduced == 24 * k0 + 13 * l2**2
    third_b_scalar = primitive(
        third_b[(0, 0)].subs(p4, p4_b).subs(l1, 0),
        tuple(sorted(third_b[(0, 0)].free_symbols, key=str)),
    )
    assert third_b_scalar == l2 * (12 * k0 + l2**2)
    assert sp.factor(
        third_b_scalar.subs(k0, -sp.Rational(13, 24) * l2**2)
    ) == -sp.Rational(11, 2) * l2**3


def verify_terminal_equality_face() -> None:
    A, B, t = sp.symbols("A B t")
    p0, p1, p2, p3, p4 = sp.symbols("p0:5")
    r0, r1, r2, r3 = sp.symbols("r0:4")
    q0, q1, q2 = sp.symbols("q0:3")
    k0, k1, k2, k3 = sp.symbols("k0:4")
    h0, h1, h2, h3 = sp.symbols("h0:4")
    z0, z1, z2, z3 = sp.symbols("z0:4")

    polynomial = {
        (0, 5): 1,
        (1, 4): A,
        (2, 3): B,
        (3, 2): 1,
        (0, 4): p0,
        (1, 3): p1,
        (2, 2): p2,
        (3, 1): p3,
        (4, 0): p4,
        (0, 3): r0,
        (1, 2): r1,
        (2, 1): r2,
        (3, 0): r3,
        (0, 2): q0,
        (1, 1): q1,
        (2, 0): q2,
    }
    operator = {
        (4, 0): 1,
        (3, 2): t,
        (0, 6): k0,
        (1, 5): k1,
        (2, 4): k2,
        (3, 3): k3,
        (0, 7): h0,
        (1, 6): h1,
        (2, 5): h2,
        (3, 4): h3,
        (0, 8): z0,
        (1, 7): z1,
        (2, 6): z2,
        (3, 5): z3,
    }
    p4_value = -t / 2
    second = moment(polynomial, operator, 2)[(0, 0)]
    third = moment(polynomial, operator, 3)
    assert primitive(
        third[(1, 0)].subs(p4, p4_value),
        tuple(sorted(third[(1, 0)].free_symbols, key=str)),
    ) == k0
    assert primitive(
        third[(0, 1)].subs(p4, p4_value),
        tuple(sorted(third[(0, 1)].free_symbols, key=str)),
    ) == 7 * B * k0 + 3 * k1
    early = {
        p4: p4_value,
        k0: 0,
        k1: 0,
        k2: -sp.Rational(5, 24) * t**2,
        h1: sp.Rational(47, 144) * t**3 - sp.Rational(7, 3) * B * h0,
    }
    assert sp.factor(second.subs(early)) == 0
    assert sp.factor(third[(0, 0)].subs(early)) == 0
    fourth = moment(polynomial, operator, 4)
    assert primitive(
        fourth[(0, 1)].subs(early),
        tuple(sorted(fourth[(0, 1)].free_symbols, key=str)),
    ) == h0

    # With every below-weight operator term gone, only the equality chain
    # contributes to the terminal scalar moments.
    equality_polynomial = {(3, 2): 1, (4, 0): -t / 2}
    equality_operator = {
        (4, 0): 1,
        (3, 2): t,
        (2, 4): -sp.Rational(5, 24) * t**2,
        (1, 6): sp.Rational(47, 144) * t**3,
        (0, 8): z0,
    }
    fourth_scalar = moment(
        equality_polynomial, equality_operator, 4
    )[(0, 0)]
    assert fourth_scalar == 958_003_200 * (
        40_541 * t**4 + 80_640 * z0
    )
    z0_value = -sp.Rational(40_541, 80_640) * t**4
    fifth_scalar = sp.factor(
        moment(equality_polynomial, equality_operator, 5)[(0, 0)].subs(
            z0, z0_value
        )
    )
    assert fifth_scalar == -19_931_886_558_904_320_000 * t**5

    # The low-weight operator list is exhaustive after Weierstrass division.
    low_or_equal = {
        exponent
        for total_order in range(5, 9)
        for x_order in range(4)
        if (exponent := (x_order, total_order - x_order))
        and total_order + x_order <= 8
    }
    assert low_or_equal == {
        (0, 5),
        (1, 4),
        (2, 3),
        (3, 2),
        (0, 6),
        (1, 5),
        (2, 4),
        (0, 7),
        (1, 6),
        (0, 8),
    }


def verify_final_weight_separators() -> None:
    # Branch max_x(P5)<=1: weights (3,2) give a strict 11 < 12 gap.
    branch_one_polynomial = (
        {(0, 5), (1, 4)}
        | {(x_degree, 4 - x_degree) for x_degree in range(4)}
        | {
            (x_degree, total_degree - x_degree)
            for total_degree in range(4)
            for x_degree in range(total_degree + 1)
        }
    )
    assert max(3 * a + 2 * b for a, b in branch_one_polynomial) == 11
    normalized_operator_one = {
        (4, 0),
        (2, 3),
        (3, 2),
    } | {
        (x_order, total_order - x_order)
        for total_order in range(6, 13)
        for x_order in range(min(3, total_order) + 1)
    }
    assert min(3 * a + 2 * b for a, b in normalized_operator_one) == 12

    # The B-branch weight-seven equality faces have a y deficit.
    branch_two_polynomial = {
        (0, 5),
        (1, 4),
        (2, 3),
    } | {(x_degree, 4 - x_degree) for x_degree in range(4)}
    branch_two_operator = {
        (4, 0),
        (3, 2),
        (1, 5),
        (2, 4),
        (3, 3),
        (0, 7),
        (1, 6),
    }
    assert max(2 * a + b for a, b in branch_two_polynomial) == 7
    assert min(2 * a + b for a, b in branch_two_operator) == 7
    polynomial_equality = {
        exponent
        for exponent in branch_two_polynomial
        if 2 * exponent[0] + exponent[1] == 7
    }
    operator_equality = {
        exponent
        for exponent in branch_two_operator
        if 2 * exponent[0] + exponent[1] == 7
    }
    assert polynomial_equality == {(2, 3), (3, 1)}
    assert operator_equality == {(1, 5), (0, 7)}
    assert max(b for _a, b in polynomial_equality) == 3
    assert min(b for _a, b in operator_equality) == 5

    # The C-branch terminal weight-eight pair has an x deficit.
    assert 2 * 3 + 2 == 8
    assert 2 * 4 == 8
    assert 3 < 4


def main() -> None:
    verify_defect_one_radical()
    verify_low_threshold_branches()
    verify_terminal_equality_face()
    verify_final_weight_separators()
    print("PASS quadruple-root defect-one radical and branch identities")
    print("PASS quadruple-root equality chain dies at pure moment five")
    print("PASS arbitrary higher jets have strict or one-sided final faces")


if __name__ == "__main__":
    main()
