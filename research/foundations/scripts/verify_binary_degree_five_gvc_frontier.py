#!/usr/bin/env python3
"""Verify the exact faces on the binary degree-five GVC frontier.

The all-order argument is in
``extended-geometry/BINARY_DEGREE_FIVE_GVC_FRONTIER.md``.  This checker
replays the decisive characteristic-zero moment identities for the eight
order-three normal forms and the order-four squarefree row.  With
``--singular`` it also verifies the two small exact radical calculations.

The uniform squarefree-quartic top-form saturation is intentionally a
separate, expensive command (``--singular-top``); it takes several minutes.
Finite-field discovery is reproduced independently by
``search_binary_degree_five_gvc_faces_mod_p.py``.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess

import sympy as sp


x, y = sp.symbols("x y")


def apply_operator(
    polynomial: sp.Expr,
    terms: dict[tuple[int, int], sp.Expr],
) -> sp.Expr:
    return sp.expand(
        sum(
            coefficient * sp.diff(polynomial, x, i, y, j)
            for (i, j), coefficient in terms.items()
        )
    )


def moment(
    polynomial: sp.Expr,
    terms: dict[tuple[int, int], sp.Expr],
    order: int,
) -> sp.Expr:
    result = sp.expand(polynomial**order)
    for _ in range(order):
        result = apply_operator(result, terms)
    return sp.expand(result)


def primitive_scalar(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    _content, primitive = sp.Poly(
        sp.expand(expression), *variables, domain=sp.QQ
    ).primitive()
    return sp.expand(primitive.as_expr())


def singular_polynomial(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


def verify_radical_with_singular(
    variables: tuple[sp.Symbol, ...],
    equations: tuple[sp.Expr, ...],
    expected: tuple[sp.Expr, ...],
) -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for --singular")
    program = f"""
ring r=0,({",".join(map(str, variables))}),dp;
ideal I={",".join(map(singular_polynomial, equations))};
ideal E={",".join(map(singular_polynomial, expected))};
LIB "primdec.lib";
ideal R=std(radical(I));
ideal left=reduce(R,std(E));
ideal right=reduce(std(E),R);
if ((size(left)==0) && (size(right)==0)) {{ print("PASS"); }}
else {{ print("FAIL"); }}
quit;
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "PASS" in result.stdout and "FAIL" not in result.stdout, result.stdout


def assert_mixed_tail(
    operator: dict[tuple[int, int], sp.Expr],
    polynomial: sp.Expr,
    multiplier: sp.Expr,
    *,
    tail_start: int,
    through: int = 6,
) -> tuple[sp.Expr, ...]:
    values = tuple(
        moment(polynomial, operator, order)
        if multiplier == 1
        else mixed_moment(polynomial, operator, multiplier, order)
        for order in range(1, through + 1)
    )
    assert any(value != 0 for value in values[: tail_start - 1])
    assert all(value == 0 for value in values[tail_start - 1 :])
    return values


def mixed_moment(
    polynomial: sp.Expr,
    terms: dict[tuple[int, int], sp.Expr],
    multiplier: sp.Expr,
    order: int,
) -> sp.Expr:
    result = sp.expand(multiplier * polynomial**order)
    for _ in range(order):
        result = apply_operator(result, terms)
    return sp.expand(result)


def triple_root_x2y3() -> None:
    a2, a3 = sp.symbols("a2 a3")
    b3, b4, b5, b6 = sp.symbols("b3 b4 b5 b6")
    c4, c5, c6, c7, c8, c9 = sp.symbols("c4:10")
    R, S, T = sp.symbols("R S T")
    polynomial = x**2 * y**3 + R * x**4 + S * x**3 * y + T * x**3
    operator = {
        (3, 0): 1,
        (2, 2): a2,
        (2, 3): a3,
        (1, 3): b3,
        (1, 4): b4,
        (1, 5): b5,
        (1, 6): b6,
        (0, 4): c4,
        (0, 5): c5,
        (0, 6): c6,
        (0, 7): c7,
        (0, 8): c8,
        (0, 9): c9,
    }
    first = {R: -b3 / 2, S: -2 * a2, T: -2 * a3}
    assert apply_operator(polynomial, operator).subs(first) == 0
    second = {
        c4: 0,
        b3: 0,
        b4: 0,
        c5: 0,
        c6: 0,
        b5: -a2 * a3 / 3,
        b6: -sp.Rational(5, 12) * a3**2,
    }
    normalized_polynomial = sp.expand(polynomial.subs(first).subs(second))
    normalized_operator = {
        exponent: sp.expand(sp.sympify(coefficient).subs(second))
        for exponent, coefficient in operator.items()
    }
    assert moment(normalized_polynomial, normalized_operator, 2) == 0
    third = sp.Poly(
        moment(normalized_polynomial, normalized_operator, 3), x, y
    )
    assert third.coeff_monomial(y**3) == -6_531_840 * a2**3
    third_branch = {
        a2: 0,
        c7: 0,
        c8: 0,
        c9: sp.Rational(1169, 2160) * a3**3,
    }
    terminal_polynomial = sp.expand(normalized_polynomial.subs(third_branch))
    terminal_operator = {
        exponent: sp.expand(sp.sympify(coefficient).subs(third_branch))
        for exponent, coefficient in normalized_operator.items()
    }
    assert moment(terminal_polynomial, terminal_operator, 3) == 0
    assert (
        moment(terminal_polynomial, terminal_operator, 4)
        == 67_315_784_417_280 * a3**4
    )

    base = {(3, 0): 1}
    mixed = assert_mixed_tail(
        base, x**2 * y**3, x, tail_start=2, through=5
    )
    assert mixed[0] == 6 * y**3
    # The distinct quadratic top form shears to this face plus -y^5/4,
    # which is four weight units lower for weights (3,1).
    assert sp.expand(
        (x + y / 2) * (x - y / 2) * y**3
        - (x**2 * y**3 - y**5 / 4)
    ) == 0


def triple_root_xy4(
    *,
    run_singular: bool,
) -> None:
    A, B, C, u, v = sp.symbols("A B C u v")
    polynomial = x * y**4 + u * x**2 * y**2 + v * x**3
    operator = {(3, 0): 1, (2, 2): A, (1, 4): B, (0, 6): C}
    branch = {v: -sp.Rational(2, 3) * A * u - 4 * B}
    scalar_moments = tuple(
        primitive_scalar(
            moment(polynomial.subs(branch), operator, order),
            (A, B, C, u),
        )
        for order in range(2, 6)
    )
    assert scalar_moments[0] == (
        -6 * A**3 * u
        - 36 * A**2 * B
        + 2 * A**2 * u**2
        + 123 * A * B * u
        + 1260 * A * C
        + 648 * B**2
        + 9 * B * u**2
        + 135 * C * u
    )
    if run_singular:
        verify_radical_with_singular(
            (A, B, C, u),
            scalar_moments,
            (B, C, A * u),
        )

    a_axis_operator = {(3, 0): 1, (2, 2): A}
    u_axis_polynomial = x * y**4 + u * x**2 * y**2
    assert_mixed_tail(
        a_axis_operator, x * y**4, x**2, tail_start=3, through=5
    )
    assert_mixed_tail(
        {(3, 0): 1}, u_axis_polynomial, x**2, tail_start=3, through=6
    )

    # Low transverse jets preceding the equality face.
    b, c4, d = sp.symbols("b c4 d")
    R, S = sp.symbols("R S")
    full_polynomial = (
        x * y**4
        + R * x**4
        + S * x**3 * y
        + u * x**2 * y**2
        + v * x**3
    )
    full_operator = {
        (3, 0): 1,
        (2, 2): A,
        (1, 3): b,
        (1, 4): B,
        (0, 4): c4,
        (0, 5): d,
        (0, 6): C,
    }
    low_branch = {
        R: -c4,
        S: -4 * b,
        v: -sp.Rational(2, 3) * A * u - 4 * B,
    }
    assert apply_operator(full_polynomial, full_operator).subs(low_branch) == 0
    low_second = sp.Poly(
        moment(full_polynomial.subs(low_branch), full_operator, 2),
        x,
        y,
    )
    assert low_second.coeff_monomial(x**2) == 54_720 * c4**2
    assert (
        low_second.coeff_monomial(y**2).subs(c4, 0)
        == 28_800 * b**2
    )
    reduced_operator = {
        exponent: sp.sympify(coefficient).subs({c4: 0, b: 0})
        for exponent, coefficient in full_operator.items()
    }
    reduced_polynomial = full_polynomial.subs(low_branch).subs(
        {c4: 0, b: 0}
    )
    low_third = sp.Poly(
        moment(reduced_polynomial, reduced_operator, 3), x, y
    )
    assert low_third.coeff_monomial(y**2) == 4_311_014_400 * d**2


def triple_root_y5() -> None:
    T, b, Z, c5, H = sp.symbols("T b Z c5 H")
    A0, B0, C0, c4 = sp.symbols("A0 B0 C0 c4")
    R, S, v = sp.symbols("R S v")
    full_polynomial = (
        y**5
        + R * x**4
        + S * x**3 * y
        + T * x**2 * y**2
        + Z * x * y**3
        + v * x**3
        + H * x**2 * y
    )
    full_operator = {
        (3, 0): 1,
        (2, 2): A0,
        (1, 3): b,
        (1, 4): B0,
        (0, 4): c4,
        (0, 5): c5,
        (0, 6): C0,
    }
    first_branch = {
        R: 0,
        S: -20 * c4,
        v: -sp.Rational(2, 3) * A0 * T - Z * b - 20 * c5,
    }
    assert apply_operator(full_polynomial, full_operator).subs(first_branch) == 0
    second = sp.Poly(
        moment(full_polynomial.subs(first_branch), full_operator, 2),
        x,
        y,
    )
    assert second.coeff_monomial(y**2) == 1_929_600 * c4**2
    assert (
        sp.factor(second.coeff_monomial(y).subs(c4, 0))
        == 576 * T * b * (2 * T + 35 * b)
    )
    reduced_operator = {
        exponent: sp.sympify(coefficient).subs(c4, 0)
        for exponent, coefficient in full_operator.items()
    }
    reduced_polynomial = full_polynomial.subs(first_branch).subs(c4, 0)
    third = sp.Poly(
        moment(reduced_polynomial, reduced_operator, 3), x, y
    )
    assert (
        sp.factor(third.coeff_monomial(x))
        == 311_040 * T**2 * b**2 * (5 * T + 84 * b)
    )
    # The two nonzero linear ratios are incompatible.
    ratio_resultant = sp.resultant(
        2 * T + 35 * b,
        5 * T + 84 * b,
        T,
    )
    assert ratio_resultant == -7 * b

    A = sp.symbols("A")
    quartic_face_polynomial = x**2 * y**2 - sp.Rational(2, 3) * A * x**3
    quartic_face_operator = {
        (3, 0): 1,
        (2, 2): A,
        (1, 4): -sp.Rational(2, 9) * A**2,
        (0, 6): sp.Rational(136, 405) * A**3,
    }
    assert (
        moment(quartic_face_polynomial, quartic_face_operator, 4)
        == 3_361_505_280 * A**4
    )

    b_face_polynomial = (
        y**5 + Z * x * y**3 - sp.Rational(17, 20) * Z**2 * x**2 * y
    )
    b_face_operator = {
        (1, 3): b,
        (0, 5): -Z * b / 20,
    }
    assert (
        moment(b_face_polynomial, b_face_operator, 3)
        == -37_666_944 * Z**3 * b**3
    )
    assert second.as_expr().subs({c4: 0, T: 0, b: 0}) == 3_859_200 * c5**2

    assert_mixed_tail({(1, 3): b}, y**5, x, tail_start=2, through=5)


def double_and_squarefree_rows(
    *,
    run_singular: bool,
) -> None:
    U, A, B = sp.symbols("U A B")
    operator = {(2, 1): 1, (4, 0): U}
    polynomial = x**5 + A * x**3 * y + B * x * y**2
    assert sp.expand(
        apply_operator(polynomial, operator) - 6 * x * (A + 20 * U)
    ) == 0
    second = moment(polynomial.subs(A, -20 * U), operator, 2)
    assert sp.expand(
        second
        - 144 * (
            10 * (B + 340 * U**2) * x**2 - 3200 * U**3 * y
        )
    ) == 0

    D, E, T, C = sp.symbols("D E T C")
    assert (
        moment(x**5 + D * y**3, {(2, 1): 1}, 2)
        == 1440 * D * x * y
    )
    assert (
        moment(x**5 + E * y**2, {(2, 1): 1}, 2)
        == 480 * E * x
    )
    last_operator = {(2, 1): 1, (5, 0): T}
    last_polynomial = x**5 + C * x**2 * y
    assert apply_operator(last_polynomial, last_operator) == 2 * (C + 60 * T)
    assert (
        moment(last_polynomial.subs(C, -60 * T), last_operator, 2)
        == 2_592_000 * T**2
    )
    # The extra squarefree term d_x d_y^2 is strictly above all four faces.
    weights = ((3, 5), (1, 2), (2, 5), (1, 3))
    for wx, wy in weights:
        assert wx + 2 * wy > 2 * wx + wy

    Hx, Jx = sp.symbols("Hx Jx")
    xy4_operator = {(2, 1): 1, (1, 3): Hx, (0, 5): Jx}
    xy4_polynomial = (
        x * y**4 - 6 * Hx * x**2 * y**2 + B * x**3
    )
    xy4_second = sp.Poly(moment(xy4_polynomial, xy4_operator, 2), x, y)
    assert sp.expand(
        xy4_second.coeff_monomial(y**2) / 576
        - (B - 2 * Hx**2 + 140 * Jx)
    ) == 0
    assert sp.expand(
        xy4_second.coeff_monomial(x) / 576
        + Hx * (B + 18 * Hx**2 - 100 * Jx)
    ) == 0
    xy4_third = moment(xy4_polynomial, xy4_operator, 3)
    assert xy4_third.subs(
        {B: -sp.Rational(29, 3) * Hx**2, Jx: Hx**2 / 12}
    ).coeff(y, 3) == 71_884_800 * Hx**3
    assert xy4_third.subs(
        {Hx: 0, B: -140 * Jx}
    ).coeff(x, 1).coeff(y, 1) == 9_638_092_800 * Jx**2

    M = sp.symbols("M")
    assert (
        moment(x * y**4, {(2, 1): 1, (0, 6): M}, 2)
        == 161_280 * M * y
    )
    K, L = sp.symbols("K L")
    terminal_xy4_operator = {
        (2, 1): 1,
        (1, 4): K,
        (0, 7): L,
    }
    terminal_xy4_polynomial = x * y**4 - 12 * K * x**2 * y
    terminal_xy4_branch = {L: -sp.Rational(23, 70) * K**2}
    assert (
        moment(
            terminal_xy4_polynomial,
            {
                exponent: sp.sympify(coefficient).subs(terminal_xy4_branch)
                for exponent, coefficient in terminal_xy4_operator.items()
            },
            3,
        )
        == -3_318_921_216 * K**3
    )

    Ay, Hy, Jy, By = sp.symbols("Ay Hy Jy By")
    y5_operator = {(2, 1): 1, (1, 3): Hy, (0, 5): Jy}
    y5_polynomial = y**5 + Ay * x * y**3 + By * x**2 * y
    y5_branch = {By: -3 * Ay * Hy - 60 * Jy}
    y5_scalars = tuple(
        primitive_scalar(
            moment(y5_polynomial.subs(y5_branch), y5_operator, order),
            (Ay, Hy, Jy),
        )
        for order in (2, 3, 4)
    )
    groebner = sp.groebner(
        y5_scalars,
        Ay,
        Jy,
        Hy,
        order="grevlex",
    )
    assert any(
        sp.expand(polynomial.as_expr() - Jy**6) == 0
        for polynomial in groebner.polys
    )
    if run_singular:
        verify_radical_with_singular(
            (Ay, By, Hy, Jy),
            (
                apply_operator(y5_polynomial, y5_operator),
                moment(y5_polynomial, y5_operator, 2),
                moment(y5_polynomial, y5_operator, 3),
                moment(y5_polynomial, y5_operator, 4),
            ),
            (Jy, By, Ay * Hy),
        )

    a_line = y**5 + Ay * x * y**3
    h_line_operator = {(2, 1): 1, (1, 3): Hy}
    a_values = assert_mixed_tail(
        {(2, 1): 1}, a_line, x**2, tail_start=3, through=6
    )
    h_values = assert_mixed_tail(
        h_line_operator, y**5, x**2, tail_start=3, through=6
    )
    assert sp.expand(
        a_values[0] - 2 * y**2 * (9 * Ay * x + 5 * y**2)
    ) == 0
    assert a_values[1] == 720 * Ay**2 * y**4
    assert sp.expand(
        h_values[0] - 10 * y**2 * (12 * Hy * x + y**2)
    ) == 0
    assert h_values[1] == 302_400 * Hy**2 * y**4


def squarefree_order_four_faces() -> None:
    lam = sp.symbols("lam")
    leading_operator = {
        (3, 1): 1,
        (2, 2): -(lam + 1),
        (1, 3): lam,
    }
    X, Y = sp.symbols("X Y")
    leading_symbol = X * Y * (X - Y) * (X - lam * Y)
    assert sp.expand(
        leading_symbol
        - (
            X**3 * Y
            - (lam + 1) * X**2 * Y**2
            + lam * X * Y**3
        )
    ) == 0

    b0, b1, b2, b3, b4, a5 = sp.symbols("b0 b1 b2 b3 b4 a5")
    quartic_correction = (
        b0 * x**4
        + b1 * x**3 * y
        + b2 * x**2 * y**2
        + b3 * x * y**3
        + b4 * y**4
    )
    fifth_operator = {(5, 0): a5}

    def correction_layer(order: int) -> sp.Poly:
        polynomial_term = order * x ** (5 * (order - 1)) * quartic_correction
        operator_term = order * apply_operator(x ** (5 * order), fifth_operator)
        for _ in range(order):
            polynomial_term = apply_operator(polynomial_term, leading_operator)
        for _ in range(order - 1):
            operator_term = apply_operator(operator_term, leading_operator)
        branch = {b0: 0, a5: -b1 / 20}
        return sp.Poly(
            sp.expand((polynomial_term + operator_term).subs(branch)),
            x,
            y,
        )

    fourth_layer = correction_layer(4)
    assert (
        fourth_layer.coeff_monomial(x**3)
        == 20_922_789_888_000 * b4
    )
    third_layer = correction_layer(3)
    assert sp.expand(
        third_layer.coeff_monomial(x**2)
        - 32_659_200 * (11 * b3 - 12 * (lam + 1) * b4)
    ) == 0
    assert third_layer.coeff_monomial(x * y) == 261_273_600 * b4
    second_layer = correction_layer(2)
    assert sp.expand(
        second_layer.coeff_monomial(y)
        - 2880 * (3 * b3 - 4 * (lam + 1) * b4)
    ) == 0
    assert sp.expand(
        second_layer.coeff_monomial(x)
        - 2880
        * (
            7 * b2
            - 6 * (lam + 1) * b3
            + 2 * (lam**2 + 4 * lam + 1) * b4
        )
    ) == 0
    first_layer = sp.expand(
        apply_operator(quartic_correction, leading_operator)
        + apply_operator(x**5, fifth_operator)
    )
    assert sp.expand(
        first_layer.subs({b0: 0, b2: 0, b3: 0, b4: 0})
        - 6 * (20 * a5 + b1)
    ) == 0

    d = sp.symbols("d")
    transverse = moment(x**5 + d * y**3, {(3, 1): 1}, 3)
    assert transverse == 65_318_400 * d * x

    t, s = sp.symbols("t s")
    a = -t / 20
    equality_operator = {(3, 1): 1, (5, 0): a}
    equality_polynomial = x**5 + t * x**3 * y + s * x * y**2
    assert apply_operator(equality_polynomial, equality_operator) == 0
    second = moment(equality_polynomial, equality_operator, 2)
    third = moment(equality_polynomial, equality_operator, 3)
    assert second == 144 * (17 * t**2 + 20 * s)
    assert (
        third.subs(s, -sp.Rational(17, 20) * t**2)
        == -37_666_944 * t**3
    )

    base_values = assert_mixed_tail(
        {(3, 1): 1}, x**5, y, tail_start=2, through=5
    )
    assert base_values[0] == 60 * x**2
    assert leading_operator[(3, 1)] == 1


def squarefree_top_saturation_with_singular() -> None:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for --singular-top")
    lam = sp.symbols("lam")
    coefficients = sp.symbols("a0:6")
    polynomial = sum(
        coefficients[index] * x ** (5 - index) * y**index
        for index in range(6)
    )
    operator = {
        (3, 1): 1,
        (2, 2): -(lam + 1),
        (1, 3): lam,
    }
    first_branch = {
        coefficients[1]: (
            (lam + 1) * coefficients[2] - lam * coefficients[3]
        )
        / 2,
        coefficients[4]: (
            (lam + 1) * coefficients[3] - coefficients[2]
        )
        / (2 * lam),
    }
    equations = []
    for order in (2, 3):
        contracted = sp.Poly(
            moment(polynomial.subs(first_branch), operator, order),
            x,
            y,
        )
        for coefficient in contracted.coeffs():
            numerator = sp.together(coefficient).as_numer_denom()[0]
            equations.append(sp.expand(numerator))
    a0, _a1, a2, a3, _a4, a5 = coefficients
    ring_variables = ("z", "lam", "a0", "a2", "a3", "a5")
    generators = ",".join(map(singular_polynomial, equations))
    program = f"""
LIB "elim.lib";
ring r=0,({",".join(ring_variables)}),dp;
option(redSB);
ideal I={generators},z*lam*(lam-1)-1;
ideal M=a0,a2,a3,a5;
list sat_data=sat(I,M);
ideal S=std(sat_data[1]);
ideal J1=a2,a3,a5;
ideal J2=a0,a2,a3;
ideal J3=a0-a5,a2-10*a5,a3-10*a5;
ideal J4=a0-lam^5*a5,a2-10*lam^3*a5,a3-10*lam^2*a5;
ideal J=std(intersect(intersect(J1,J2),intersect(J3,J4))
            +ideal(z*lam*(lam-1)-1));
ideal left=reduce(S,J);
ideal right=reduce(J,S);
if ((size(left)==0) && (size(right)==0)) {{ print("PASS"); }}
else {{ print("FAIL"); }}
quit;
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "PASS" in result.stdout and "FAIL" not in result.stdout, result.stdout


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--singular",
        action="store_true",
        help="verify the two small face radicals with Singular",
    )
    parser.add_argument(
        "--singular-top",
        action="store_true",
        help="run the several-minute uniform quartic top-form saturation",
    )
    args = parser.parse_args()

    triple_root_x2y3()
    triple_root_xy4(run_singular=args.singular)
    triple_root_y5()
    double_and_squarefree_rows(run_singular=args.singular)
    squarefree_order_four_faces()
    if args.singular_top:
        squarefree_top_saturation_with_singular()

    print("PASS r=3: all eight cubic-leading quintic normal forms close")
    print("PASS r=4: squarefree cross-ratio row has no exceptional parameter")
    print("PASS terminal components are one-sided with mixed cutoffs")
    if args.singular:
        print("PASS Singular: exact residual face radicals")
    if args.singular_top:
        print("PASS Singular: uniform squarefree-quartic top saturation")


if __name__ == "__main__":
    main()
