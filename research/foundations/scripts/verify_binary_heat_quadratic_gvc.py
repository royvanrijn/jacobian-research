#!/usr/bin/env python3
"""Exact symbolic checks for low-order binary GVC theorems.

The all-order cutoff is the derivative-count proof in
extended-geometry/SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md.  This script
derives the quadratic and cubic first-moment normal forms, verifies their
second moments (including every operator term that can act), and checks the
heat square identity on the generic degree-six heat polynomial.
"""

from __future__ import annotations

import sympy as sp


def main() -> None:
    x, y = sp.symbols("x y")
    A, B, C = sp.symbols("A B C")
    a, b, c, d, e, f = sp.symbols("a b c d e f")

    def operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x)
            + A * sp.diff(polynomial, x, 2)
            + B * sp.diff(polynomial, x, y)
            + C * sp.diff(polynomial, y, 2)
        )

    general = a * x**2 + b * x * y + c * y**2 + d * x + e * y + f
    first = sp.Poly(operator(general), x, y)
    assert first.coeff_monomial(x) == 2 * a
    assert first.coeff_monomial(y) == b
    assert first.coeff_monomial(1) == d + 2 * A * a + B * b + 2 * C * c

    normal = c * (y**2 - 2 * C * x) + e * y + f
    assert operator(normal) == 0
    second = sp.factor(operator(operator(sp.expand(normal**2))))
    assert second == 16 * C**2 * c**2

    D, E, F, G = sp.symbols("D E F G")

    def cubic_coupled_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            operator(polynomial)
            + D * sp.diff(polynomial, x, 3)
            + E * sp.diff(polynomial, x, 2, y)
            + F * sp.diff(polynomial, x, y, 2)
            + G * sp.diff(polynomial, y, 3)
        )

    assert cubic_coupled_operator(normal) == 0
    coupled_second = sp.factor(
        cubic_coupled_operator(
            cubic_coupled_operator(sp.expand(normal**2))
        )
    )
    assert coupled_second == 16 * C**2 * c**2

    c0, c2, c5, c9 = sp.symbols("c0 c2 c5 c9")
    coupled_cubic_normal = (
        c0
        + c2 * y
        + c5 * y**2
        + c9 * y**3
        + x * (6 * B * C * c9 - 2 * C * c5 - 6 * C * c9 * y - 6 * G * c9)
    )
    assert cubic_coupled_operator(coupled_cubic_normal) == 0
    coupled_cubic_second = sp.Poly(
        cubic_coupled_operator(
            cubic_coupled_operator(sp.expand(coupled_cubic_normal**2))
        ),
        x,
        y,
    )
    assert (
        coupled_cubic_second.coeff_monomial(y**2)
        == 144 * C**2 * c9**2
    )
    assert sp.factor(
        coupled_cubic_second.as_expr().subs(c9, 0)
    ) == 16 * C**2 * c5**2
    assert sp.factor(
        coupled_cubic_second.as_expr().subs(C, 0)
    ) == 648 * G**2 * c9**2

    # For cubic P, operator terms of order at least six cannot occur in the
    # second moment. Add the complete order-four and order-five pieces.
    H, I, J, K, L = sp.symbols("H I J K L")
    U, V, W, X, Y, Z = sp.symbols("U V W X Y Z")

    def full_relevant_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            cubic_coupled_operator(polynomial)
            + H * sp.diff(polynomial, x, 4)
            + I * sp.diff(polynomial, x, 3, y)
            + J * sp.diff(polynomial, x, 2, y, 2)
            + K * sp.diff(polynomial, x, y, 3)
            + L * sp.diff(polynomial, y, 4)
            + U * sp.diff(polynomial, x, 5)
            + V * sp.diff(polynomial, x, 4, y)
            + W * sp.diff(polynomial, x, 3, y, 2)
            + X * sp.diff(polynomial, x, 2, y, 3)
            + Y * sp.diff(polynomial, x, y, 4)
            + Z * sp.diff(polynomial, y, 5)
        )

    assert full_relevant_operator(coupled_cubic_normal) == 0
    full_cubic_second = sp.Poly(
        full_relevant_operator(
            full_relevant_operator(sp.expand(coupled_cubic_normal**2))
        ),
        x,
        y,
    )
    assert (
        full_cubic_second.coeff_monomial(y**2)
        == 144 * C**2 * c9**2
    )
    assert sp.factor(
        full_cubic_second.as_expr().subs(c9, 0)
    ) == 16 * C**2 * c5**2
    assert sp.factor(
        full_cubic_second.as_expr().subs(C, 0)
    ) == 648 * G**2 * c9**2

    # Quartic P only sees the operator 4-jet in the first equation and the
    # 7-jet in the second. These are the three successive transverse
    # branches predicted by formal drift straightening.
    jet_coefficients: dict[tuple[int, int], sp.Expr] = {
        (1, 0): sp.Integer(1)
    }
    for total_order in range(2, 8):
        for y_order in range(total_order + 1):
            x_order = total_order - y_order
            jet_coefficients[(x_order, y_order)] = sp.symbols(
                f"j{x_order}{y_order}"
            )

    def jet_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in jet_coefficients.items()
            )
        )

    quartic_coefficients = sp.symbols("q0:5")
    quartic_boundary = sum(
        coefficient * y**degree
        for degree, coefficient in enumerate(quartic_coefficients)
    )
    j20 = jet_coefficients[(2, 0)]
    j11 = jet_coefficients[(1, 1)]
    j02 = jet_coefficients[(0, 2)]
    j12 = jet_coefficients[(1, 2)]
    j03 = jet_coefficients[(0, 3)]
    j04 = jet_coefficients[(0, 4)]
    effective_order_two = j02
    effective_order_three = j03 - j11 * j02
    effective_order_four = (
        j04
        - j11 * j03
        + (j11**2 - j12 + j20 * j02) * j02
    )
    quartic_normal = sp.expand(
        quartic_boundary
        - x
        * (
            effective_order_two * sp.diff(quartic_boundary, y, 2)
            + effective_order_three * sp.diff(quartic_boundary, y, 3)
            + effective_order_four * sp.diff(quartic_boundary, y, 4)
        )
        + 12 * x**2 * j02**2 * quartic_coefficients[4]
    )
    assert jet_operator(quartic_normal) == 0
    quartic_second = sp.Poly(
        jet_operator(jet_operator(sp.expand(quartic_normal**2))),
        x,
        y,
    )
    top_quartic = quartic_coefficients[4]
    assert (
        quartic_second.coeff_monomial(x**2)
        == 2304 * j02**4 * top_quartic**2
    )
    assert (
        sp.Poly(quartic_second.as_expr().subs(j02, 0), x, y)
        .coeff_monomial(y**2)
        == 15552 * j03**2 * top_quartic**2
    )
    assert (
        sp.factor(
            quartic_second.as_expr().subs({j02: 0, j03: 0})
        )
        == 39168 * j04**2 * top_quartic**2
    )

    # Lowest positive order two, cubic P: double-line quadratic symbol.
    s0, s1, s2, s4, s5, s8, s9 = sp.symbols(
        "s0 s1 s2 s4 s5 s8 s9"
    )
    sa, sb, sc, sd = sp.symbols("sa sb sc sd")
    square_four = sp.symbols("sf0:5")
    square_five = sp.symbols("sg0:6")
    square_normal = (
        -sc * s8 * x**2
        - 3 * sd * s9 * x**2
        + s0
        + s1 * x
        + s2 * y
        + s4 * x * y
        + s5 * y**2
        + s8 * x * y**2
        + s9 * y**3
    )

    def square_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2)
            + sa * sp.diff(polynomial, x, 3)
            + sb * sp.diff(polynomial, x, 2, y)
            + sc * sp.diff(polynomial, x, y, 2)
            + sd * sp.diff(polynomial, y, 3)
            + sum(
                square_four[y_order]
                * sp.diff(
                    polynomial,
                    x,
                    4 - y_order,
                    y,
                    y_order,
                )
                for y_order in range(5)
            )
            + sum(
                square_five[y_order]
                * sp.diff(
                    polynomial,
                    x,
                    5 - y_order,
                    y,
                    y_order,
                )
                for y_order in range(6)
            )
        )

    assert square_operator(square_normal) == 0
    square_second = sp.Poly(
        square_operator(
            square_operator(sp.expand(square_normal**2))
        ),
        x,
        y,
    )
    assert square_second.coeff_monomial(y) == 96 * sd * s8**2
    assert (
        sp.factor(square_second.as_expr().subs(sd, 0))
        == 24 * s8**2 * (sc**2 + 4 * square_four[4])
    )
    square_third_value = sp.expand(square_normal**3)
    for _ in range(3):
        square_third_value = square_operator(square_third_value)
    assert (
        sp.factor(
            square_third_value.subs(
                {
                    sd: 0,
                    square_four[4]: -sc**2 / 4,
                }
            )
        )
        == -4608 * sc**3 * s8**3
    )

    # The distinct-line orbit, normalized to dx*dy.
    t0, t1, t2, t3, t5, t6, t9 = sp.symbols(
        "t0 t1 t2 t3 t5 t6 t9"
    )
    ta, tb, tc, td = sp.symbols("ta tb tc td")
    split_four = sp.symbols("tf0:5")
    split_normal = (
        -6 * (ta * t6 + td * t9) * x * y
        + t0
        + t1 * x
        + t2 * y
        + t3 * x**2
        + t5 * y**2
        + t6 * x**3
        + t9 * y**3
    )

    def split_quadratic_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, y)
            + ta * sp.diff(polynomial, x, 3)
            + tb * sp.diff(polynomial, x, 2, y)
            + tc * sp.diff(polynomial, x, y, 2)
            + td * sp.diff(polynomial, y, 3)
            + sum(
                split_four[y_order]
                * sp.diff(
                    polynomial,
                    x,
                    4 - y_order,
                    y,
                    y_order,
                )
                for y_order in range(5)
            )
        )

    assert split_quadratic_operator(split_normal) == 0
    split_second = sp.Poly(
        split_quadratic_operator(
            split_quadratic_operator(sp.expand(split_normal**2))
        ),
        x,
        y,
    )
    assert split_second.coeff_monomial(x * y) == 72 * t6 * t9
    assert (
        sp.factor(
            split_second.as_expr().subs({t9: 0, t5: 0})
        )
        == 288 * ta**2 * t6**2
    )
    assert (
        sp.factor(
            split_second.as_expr().subs({t6: 0, t3: 0})
        )
        == 288 * td**2 * t9**2
    )

    # Lowest order two, quartic P. First the distinct-line orbit with
    # leading quartic x^4.
    split_quartic_jet = {(1, 1): sp.Integer(1)}
    for total_order in range(3, 7):
        for y_order in range(total_order + 1):
            x_order = total_order - y_order
            split_quartic_jet[(x_order, y_order)] = sp.symbols(
                f"k{x_order}{y_order}"
            )

    def split_quartic_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in split_quartic_jet.items()
            )
        )

    kp00, kp10, kp11, kp20, kp22, kp30, kp33 = sp.symbols(
        "kp00 kp10 kp11 kp20 kp22 kp30 kp33"
    )
    k30 = split_quartic_jet[(3, 0)]
    k21 = split_quartic_jet[(2, 1)]
    k03 = split_quartic_jet[(0, 3)]
    k40 = split_quartic_jet[(4, 0)]
    split_quartic_normal = (
        x**4
        - 12 * k30 * x**2 * y
        + (
            -6 * k03 * kp33
            + 24 * k21 * k30
            - 6 * k30 * kp30
            - 24 * k40
        )
        * x
        * y
        + kp00
        + kp10 * x
        + kp11 * y
        + kp20 * x**2
        + kp22 * y**2
        + kp30 * x**3
        + kp33 * y**3
    )
    assert split_quartic_operator(split_quartic_normal) == 0
    split_quartic_second = sp.Poly(
        split_quartic_operator(
            split_quartic_operator(
                sp.expand(split_quartic_normal**2)
            )
        ),
        x,
        y,
    )
    assert split_quartic_second.coeff_monomial(x**2 * y) == 144 * kp33
    assert (
        split_quartic_second.coeff_monomial(x**2).subs(kp33, 0)
        == 48 * (132 * k30**2 + kp22)
    )
    split_quartic_branch = {
        kp33: 0,
        kp22: -132 * k30**2,
    }
    assert (
        sp.factor(
            split_quartic_second.coeff_monomial(y).subs(
                split_quartic_branch
            )
        )
        == 9216 * k30**3
    )
    assert (
        sp.factor(
            split_quartic_second.as_expr().subs(
                {
                    **split_quartic_branch,
                    k30: 0,
                    kp22: 0,
                }
            )
        )
        == 31104 * k40**2
    )

    # Double-line orbit, leading quartic x*y^3.
    double_quartic_jet = {(2, 0): sp.Integer(1)}
    for total_order in range(3, 9):
        for y_order in range(total_order + 1):
            x_order = total_order - y_order
            double_quartic_jet[(x_order, y_order)] = sp.symbols(
                f"r{x_order}{y_order}"
            )

    def double_quartic_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in double_quartic_jet.items()
            )
        )

    def double_quartic_operator_up_to_six(
        polynomial: sp.Expr,
    ) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in double_quartic_jet.items()
                if x_order + y_order <= 6
            )
        )

    def specialized_double_quartic_operator(
        polynomial: sp.Expr,
        substitutions: dict[sp.Expr, sp.Expr],
    ) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient.subs(substitutions)
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in double_quartic_jet.items()
                if coefficient.subs(substitutions) != 0
            )
        )

    rp00, rp10, rp11, rp21, rp22, rp32, rp33 = sp.symbols(
        "rp00 rp10 rp11 rp21 rp22 rp32 rp33"
    )
    r30 = double_quartic_jet[(3, 0)]
    r21 = double_quartic_jet[(2, 1)]
    r12 = double_quartic_jet[(1, 2)]
    r03 = double_quartic_jet[(0, 3)]
    r13 = double_quartic_jet[(1, 3)]
    r04 = double_quartic_jet[(0, 4)]
    r05 = double_quartic_jet[(0, 5)]
    r06 = double_quartic_jet[(0, 6)]
    double_xy3_normal = (
        x * y**3
        + (
            3 * r03 * r30
            - 3 * r03 * rp33
            + 3 * r12 * r21
            - r12 * rp32
            - 3 * r13
        )
        * x**2
        - r03 * x**3
        - 3 * r12 * x**2 * y
        + rp00
        + rp10 * x
        + rp11 * y
        + rp21 * x * y
        + rp22 * y**2
        + rp32 * x * y**2
        + rp33 * y**3
    )
    assert double_quartic_operator(double_xy3_normal) == 0
    double_xy3_second = sp.Poly(
        double_quartic_operator_up_to_six(
            double_quartic_operator_up_to_six(
                sp.expand(double_xy3_normal**2)
            )
        ),
        x,
        y,
    )
    assert double_xy3_second.coeff_monomial(x**2) == 792 * r03**2
    assert sp.expand(
        double_xy3_second.coeff_monomial(y**2).subs(r03, 0)
        - 72 * (20 * r04 + r12**2)
    ) == 0
    xy3_cancellation = {
        r03: 0,
        r04: -r12**2 / 20,
    }
    assert (
        sp.factor(
            double_xy3_second.coeff_monomial(x).subs(
                xy3_cancellation
            )
        )
        == -sp.Rational(288, 5) * r12**3
    )
    assert (
        sp.factor(
            double_xy3_second.coeff_monomial(y).subs(
                {r03: 0, r12: 0, r04: 0}
            )
        )
        == 2880 * r05
    )
    assert sp.expand(
        double_xy3_second.as_expr().subs(
            {
                r03: 0,
                r12: 0,
                r04: 0,
                r05: 0,
            }
        )
        - 24 * (120 * r06 + 51 * r13**2)
    ) == 0
    xy3_third_branch = {
        r03: 0,
        r12: 0,
        r04: 0,
        r05: 0,
        r06: -sp.Rational(17, 40) * r13**2,
    }
    double_xy3_third = sp.expand(
        double_xy3_normal.subs(xy3_third_branch) ** 3
    )
    for _ in range(3):
        double_xy3_third = specialized_double_quartic_operator(
            double_xy3_third,
            xy3_third_branch,
        )
    assert sp.factor(double_xy3_third) == -3604176 * r13**3

    # Double-line orbit, leading quartic y^4. The first three moments
    # leave one finite resultant branch.
    zp00, zp10, zp11, zp21, zp22, zp32, zp33 = sp.symbols(
        "zp00 zp10 zp11 zp21 zp22 zp32 zp33"
    )
    double_y4_normal = (
        y**4
        + (
            12 * r03 * r21
            - 3 * r03 * zp33
            - 12 * r04
            - r12 * zp32
        )
        * x**2
        - 12 * r03 * x**2 * y
        + zp00
        + zp10 * x
        + zp11 * y
        + zp21 * x * y
        + zp22 * y**2
        + zp32 * x * y**2
        + zp33 * y**3
    )
    assert double_quartic_operator(double_y4_normal) == 0
    double_y4_second = sp.Poly(
        double_quartic_operator_up_to_six(
            double_quartic_operator_up_to_six(
                sp.expand(double_y4_normal**2)
            )
        ),
        x,
        y,
    )
    assert double_y4_second.coeff_monomial(y**2) == 17856 * r03**2
    residual_second = sp.factor(
        double_y4_second.as_expr().subs(r03, 0) / 24
    )
    residual_s = (
        1728 * r04**2
        - 48 * r12**2 * r04
        + 112 * r12 * r04 * zp32
        + 4 * r04 * zp32**2
        - 4 * r12**3 * zp32
        + r12**2 * zp32**2
    )
    assert sp.expand(residual_second - residual_s) == 0
    y4_third_branch = {r03: 0}
    double_y4_third = sp.expand(
        double_y4_normal.subs(y4_third_branch) ** 3
    )
    for _ in range(3):
        double_y4_third = specialized_double_quartic_operator(
            double_y4_third,
            y4_third_branch,
        )
    residual_t = (
        815616 * r04**3
        - 13824 * r12**2 * r04**2
        + 52416 * r12 * r04**2 * zp32
        + 1152 * r04**2 * zp32**2
        - 1584 * r12**3 * r04 * zp32
        + 936 * r12**2 * r04 * zp32**2
        + 36 * r12 * r04 * zp32**3
        - 36 * r12**4 * zp32**2
        + r12**3 * zp32**3
    )
    assert sp.expand(double_y4_third - 576 * residual_t) == 0
    residual_resultant = sp.factor(
        sp.resultant(residual_s, residual_t, r04)
    )
    residual_sextic = (
        -2583360 * r12**6
        + 1828368 * r12**5 * zp32
        + 1514304 * r12**4 * zp32**2
        + 502328 * r12**3 * zp32**3
        + 80916 * r12**2 * zp32**4
        + 6117 * r12 * zp32**5
        + 92 * zp32**6
    )
    assert sp.expand(
        residual_resultant
        - 1769472 * r12**3 * zp32**3 * residual_sextic
    ) == 0

    # Minimal three-term realization of the residual: moment four gives
    # a coprime degree-eight resultant.
    minimal_residual_polynomial = (
        y**4
        + zp32 * x * y**2
        - (12 * r04 + r12 * zp32) * x**2
    )

    def minimal_residual_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2)
            + r12 * sp.diff(polynomial, x, y, 2)
            + r04 * sp.diff(polynomial, y, 4)
        )

    minimal_fourth = sp.expand(minimal_residual_polynomial**4)
    for _ in range(4):
        minimal_fourth = minimal_residual_operator(minimal_fourth)
    residual_m = sp.factor(minimal_fourth / 17280)
    residual_sm_resultant = sp.factor(
        sp.resultant(residual_s, residual_m, r04)
    )
    residual_octavic = (
        -2557094400 * r12**8
        + 639596160 * r12**7 * zp32
        + 532247424 * r12**6 * zp32**2
        + 236419896 * r12**5 * zp32**3
        + 52199245 * r12**4 * zp32**4
        + 7972150 * r12**3 * zp32**5
        + 720528 * r12**2 * zp32**6
        + 39490 * r12 * zp32**7
        + 550 * zp32**8
    )
    sm_scalar = -28991029248
    assert sp.expand(
        residual_sm_resultant
        - sm_scalar * r12**4 * zp32**4 * residual_octavic
    ) == 0
    residual_t_variable = sp.symbols("residual_t")
    residual_sextic_univariate = residual_sextic.subs(
        {r12: 1, zp32: residual_t_variable}
    )
    residual_octavic_univariate = residual_octavic.subs(
        {r12: 1, zp32: residual_t_variable}
    )
    assert sp.gcd(
        residual_sextic_univariate,
        residual_octavic_univariate,
    ) == 1
    assert sp.resultant(
        residual_sextic_univariate,
        residual_octavic_univariate,
        residual_t_variable,
    ) == (
        -22002331580862445954532620608845574194895939575073794373253473026129579212800
    )

    # Lowest order three, quartic P, triple-root leading symbol.
    tr_u, tr_j, tr_k, tr_e, tr_d = sp.symbols(
        "tr_u tr_j tr_k tr_e tr_d"
    )
    triple_a_polynomial = (
        x**2 * y**2 - sp.Rational(2, 3) * tr_u * x**3
    )

    def triple_a_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 3)
            + tr_u * sp.diff(polynomial, x, 2, y, 2)
            - sp.Rational(2, 9)
            * tr_u**2
            * sp.diff(polynomial, x, y, 4)
            + tr_j * sp.diff(polynomial, y, 5)
            + tr_k * sp.diff(polynomial, y, 6)
        )

    assert triple_a_operator(triple_a_polynomial) == 0
    triple_a_third = sp.expand(triple_a_polynomial**3)
    for _ in range(3):
        triple_a_third = triple_a_operator(triple_a_third)
    assert sp.expand(
        triple_a_third
        - 3840
        * (
            405 * tr_j * y
            + 405 * tr_k
            - 136 * tr_u**3
        )
    ) == 0
    triple_a_fourth = sp.expand(triple_a_polynomial**4)
    triple_a_branch = {
        tr_j: 0,
        tr_k: sp.Rational(136, 405) * tr_u**3,
    }
    for _ in range(4):
        triple_a_fourth = triple_a_operator(
            triple_a_fourth
        ).subs(triple_a_branch)
    assert sp.factor(triple_a_fourth) == 3361505280 * tr_u**4

    triple_b_polynomial = x * y**3 - tr_d * x**3

    def triple_b_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 3)
            + tr_d * sp.diff(polynomial, x, y, 3)
            + tr_e * sp.diff(polynomial, y, 4)
        )

    assert triple_b_operator(triple_b_polynomial) == 0
    triple_b_second = sp.expand(triple_b_polynomial**2)
    for _ in range(2):
        triple_b_second = triple_b_operator(triple_b_second)
    assert sp.factor(triple_b_second) == 1584 * tr_d**2
    triple_b_third = sp.expand(triple_b_polynomial**3)
    for _ in range(3):
        triple_b_third = triple_b_operator(triple_b_third)
    assert sp.expand(
        triple_b_third
        - 12960 * (101 * tr_d**3 + 504 * tr_e**2 * y)
    ) == 0

    triple_c_polynomial = y**4 - 4 * tr_e * x**3

    def triple_c_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 3)
            + tr_e * sp.diff(polynomial, y, 4)
        )

    assert triple_c_operator(triple_c_polynomial) == 0
    triple_c_second = sp.expand(triple_c_polynomial**2)
    for _ in range(2):
        triple_c_second = triple_c_operator(triple_c_second)
    assert sp.factor(triple_c_second) == 49536 * tr_e**2

    # Lowest order three, quartic P, double-root leading symbol:
    # all branches except the pure y^4 endpoint.
    dr_t = sp.symbols("dr_t")
    double_root_x4_polynomial = (
        x**4 - 12 * dr_t * x**2 * y - 132 * dr_t**2 * y**2
    )

    def double_root_x4_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2, y)
            + dr_t * sp.diff(polynomial, x, 4)
        )

    assert double_root_x4_operator(double_root_x4_polynomial) == 0
    double_root_x4_third = sp.expand(double_root_x4_polynomial**3)
    for _ in range(3):
        double_root_x4_third = double_root_x4_operator(
            double_root_x4_third
        )
    assert (
        sp.factor(double_root_x4_third)
        == 129392640 * dr_t**3
    )

    dr_e, dr_h, dr_b = sp.symbols("dr_e dr_h dr_b")
    dr_j = -(
        248 * dr_e**2 * dr_b**2
        + 56 * dr_e * dr_h * dr_b
        + dr_h**2
    ) / 20
    double_root_xy3_polynomial = (
        x * y**3
        + dr_b * y**4
        + (-12 * dr_e * dr_b - 3 * dr_h) * x**2 * y
        - 10 * dr_e * x**3
    )

    def double_root_xy3_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2, y)
            + dr_e * sp.diff(polynomial, y, 4)
            + dr_h * sp.diff(polynomial, x, y, 3)
            + dr_j * sp.diff(polynomial, y, 5)
        )

    assert double_root_xy3_operator(double_root_xy3_polynomial) == 0
    double_root_xy3_second = sp.expand(double_root_xy3_polynomial**2)
    for _ in range(2):
        double_root_xy3_second = double_root_xy3_operator(
            double_root_xy3_second
        )
    assert double_root_xy3_second == 0
    double_root_xy3_third = sp.expand(double_root_xy3_polynomial**3)
    for _ in range(3):
        double_root_xy3_third = double_root_xy3_operator(
            double_root_xy3_third
        )
    assert (
        sp.Poly(double_root_xy3_third, x, y).coeff_monomial(x)
        == 7827840 * dr_e**2
    )
    assert (
        sp.factor(double_root_xy3_third.subs(dr_e, 0))
        == -528768 * dr_h**3
    )

    dr_r, dr_z, dr_u, dr_p33 = sp.symbols(
        "dr_r dr_z dr_u dr_p33"
    )
    double_root_y4_polynomial = (
        y**4
        - 12 * dr_e * x**2 * y
        + dr_r * x**3
        + dr_z * x * y**2
        + dr_p33 * y**3
    )

    def double_root_y4_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2, y)
            + dr_e * sp.diff(polynomial, y, 4)
            + dr_h * sp.diff(polynomial, x, y, 3)
            + dr_u * sp.diff(polynomial, x, 2, y, 2)
        )

    assert double_root_y4_operator(double_root_y4_polynomial) == 0
    double_root_y4_second = sp.expand(double_root_y4_polynomial**2)
    for _ in range(2):
        double_root_y4_second = double_root_y4_operator(
            double_root_y4_second
        )
    assert sp.expand(
        double_root_y4_second
        - 96
        * (
            372 * dr_e**2
            + 6 * dr_h * dr_r
            + dr_r * dr_z
        )
    ) == 0
    double_root_y4_third = sp.expand(double_root_y4_polynomial**3)
    for _ in range(3):
        double_root_y4_third = double_root_y4_operator(
            double_root_y4_third
        )
    assert (
        sp.Poly(double_root_y4_third, x, y).coeff_monomial(y)
        == 51840 * dr_r**2
    )

    # Lowest order three, quartic P, squarefree leading symbol.  The
    # leading nullcone consists of three fourth powers in one stabilizer
    # orbit, and its x^4 tip has the double-root weighted correction face.
    sf_a, sf_b, sf_c = sp.symbols("sf_a sf_b sf_c")
    squarefree_leading_polynomial = (
        sf_a * x**4
        - sp.Rational(2, 3) * sf_b * x**3 * y
        + sf_b * x**2 * y**2
        - sp.Rational(2, 3) * sf_b * x * y**3
        + sf_c * y**4
    )

    def squarefree_leading_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2, y)
            + sp.diff(polynomial, x, y, 2)
        )

    squarefree_leading_second = sp.Poly(
        squarefree_leading_operator(
            squarefree_leading_operator(
                sp.expand(squarefree_leading_polynomial**2)
            )
        ),
        x,
        y,
    )
    squarefree_leading_ideal = squarefree_leading_second.coeffs()
    squarefree_groebner = sp.groebner(
        squarefree_leading_ideal,
        sf_a,
        sf_b,
        sf_c,
        order="lex",
    )
    assert [polynomial.as_expr() for polynomial in squarefree_groebner.polys] == [
        sf_a * sf_b - sf_b * sf_c,
        6 * sf_a * sf_c - sf_b * sf_c,
        sf_b**2 - 6 * sf_b * sf_c,
    ]

    sf_operator_coefficients: dict[tuple[int, int], sp.Expr] = {}
    for total_order in (4, 5):
        for y_order in range(total_order + 1):
            x_order = total_order - y_order
            sf_operator_coefficients[(x_order, y_order)] = sp.symbols(
                f"sf_l{x_order}{y_order}"
            )
    sf_polynomial_coefficients: dict[tuple[int, int], sp.Expr] = {}
    for total_degree in range(4):
        for y_degree in range(total_degree + 1):
            x_degree = total_degree - y_degree
            sf_polynomial_coefficients[(x_degree, y_degree)] = sp.symbols(
                f"sf_p{x_degree}{y_degree}"
            )

    def squarefree_jet_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            squarefree_leading_operator(polynomial)
            + sum(
                coefficient
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in sf_operator_coefficients.items()
            )
        )

    sf_t = sf_operator_coefficients[(4, 0)]
    sf_p21 = sf_polynomial_coefficients[(2, 1)]
    sf_p12 = sf_polynomial_coefficients[(1, 2)]
    sf_p03 = sf_polynomial_coefficients[(0, 3)]
    sf_p02 = sf_polynomial_coefficients[(0, 2)]
    sf_p30 = sf_polynomial_coefficients[(3, 0)]
    squarefree_x4_polynomial = x**4 + sum(
        coefficient * x**x_degree * y**y_degree
        for (x_degree, y_degree), coefficient
        in sf_polynomial_coefficients.items()
    )
    squarefree_first = squarefree_jet_operator(
        squarefree_x4_polynomial
    )
    assert squarefree_first == 24 * sf_t + 2 * sf_p21 + 2 * sf_p12
    squarefree_first_branch = {sf_p21: -12 * sf_t - sf_p12}
    squarefree_second = sp.Poly(
        squarefree_jet_operator(
            squarefree_jet_operator(
                sp.expand(
                    squarefree_x4_polynomial.subs(
                        squarefree_first_branch
                    )
                    ** 2
                )
            )
        ),
        x,
        y,
    )
    assert (
        squarefree_second.coeff_monomial(x)
        == 96 * (5 * sf_p12 + 6 * sf_p03)
    )
    assert squarefree_second.coeff_monomial(y) == 288 * sf_p03
    squarefree_second_constant = squarefree_second.coeff_monomial(1)
    sf_l31 = sf_operator_coefficients[(3, 1)]
    sf_l22 = sf_operator_coefficients[(2, 2)]
    assert sp.expand(
        squarefree_second_constant
        - 48
        * (
            12 * sf_l22 * sf_p03
            + 12 * sf_l31 * sf_p03
            + 20 * sf_l31 * sf_p12
            + 264 * sf_t**2
            - 24 * sf_t * sf_p03
            - 52 * sf_t * sf_p12
            + 2 * sf_p02
            - 2 * sf_p03 * sf_p12
            + 3 * sf_p03 * sf_p30
            - sf_p12**2
            + 2 * sf_p12 * sf_p30
        )
    ) == 0

    squarefree_face_polynomial = (
        x**4 - 12 * sf_t * x**2 * y - 132 * sf_t**2 * y**2
    )

    def squarefree_face_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2, y)
            + sp.diff(polynomial, x, y, 2)
            + sf_t * sp.diff(polynomial, x, 4)
        )

    assert squarefree_face_operator(squarefree_face_polynomial) == 0
    squarefree_face_second = sp.expand(squarefree_face_polynomial**2)
    for _ in range(2):
        squarefree_face_second = squarefree_face_operator(
            squarefree_face_second
        )
    assert squarefree_face_second == 0
    squarefree_face_third = sp.expand(squarefree_face_polynomial**3)
    for _ in range(3):
        squarefree_face_third = squarefree_face_operator(
            squarefree_face_third
        )
    assert sp.factor(squarefree_face_third) == 129392640 * sf_t**3

    # First unresolved r=2, deg(P)=5 row: exact top-form reduction.
    qf_coefficients = sp.symbols("qf0:6")
    quadratic_quintic = sum(
        coefficient * x ** (5 - y_degree) * y**y_degree
        for y_degree, coefficient in enumerate(qf_coefficients)
    )
    distinct_quintic_first = sp.Poly(
        sp.diff(quadratic_quintic, x, y),
        x,
        y,
    )
    distinct_quintic_branch = {
        qf_coefficients[1]: 0,
        qf_coefficients[2]: 0,
        qf_coefficients[3]: 0,
        qf_coefficients[4]: 0,
    }
    assert distinct_quintic_first.as_expr().subs(
        distinct_quintic_branch
    ) == 0
    assert all(
        coefficient == 0
        for coefficient in sp.solve(
            distinct_quintic_first.coeffs(),
            qf_coefficients[1:5],
            dict=True,
        )[0].values()
    )
    distinct_quintic_top = (
        qf_coefficients[0] * x**5 + qf_coefficients[5] * y**5
    )
    distinct_quintic_second = sp.diff(
        distinct_quintic_top**2,
        x,
        2,
        y,
        2,
    )
    assert (
        sp.factor(distinct_quintic_second)
        == 800
        * qf_coefficients[0]
        * qf_coefficients[5]
        * x**3
        * y**3
    )
    double_quintic_first = sp.Poly(
        sp.diff(quadratic_quintic, x, 2),
        x,
        y,
    )
    assert sp.solve(
        double_quintic_first.coeffs(),
        qf_coefficients[:4],
        dict=True,
    ) == [
        {
            qf_coefficients[0]: 0,
            qf_coefficients[1]: 0,
            qf_coefficients[2]: 0,
            qf_coefficients[3]: 0,
        }
    ]

    # The distinct-root x^5 tip closes after the full second moment.
    dq_operator_coefficients: dict[tuple[int, int], sp.Expr] = {
        (1, 1): sp.Integer(1)
    }
    for total_order in range(3, 9):
        for y_order in range(total_order + 1):
            x_order = total_order - y_order
            dq_operator_coefficients[(x_order, y_order)] = sp.symbols(
                f"dq_l{x_order}{y_order}"
            )
    dq_polynomial_coefficients: dict[tuple[int, int], sp.Expr] = {}
    for total_degree in range(5):
        for y_degree in range(total_degree + 1):
            x_degree = total_degree - y_degree
            dq_polynomial_coefficients[(x_degree, y_degree)] = sp.symbols(
                f"dq_p{x_degree}{y_degree}"
            )

    def distinct_quintic_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in dq_operator_coefficients.items()
            )
        )

    distinct_quintic_polynomial = x**5 + sum(
        coefficient * x**x_degree * y**y_degree
        for (x_degree, y_degree), coefficient
        in dq_polynomial_coefficients.items()
    )
    distinct_quintic_first_full = sp.Poly(
        distinct_quintic_operator(distinct_quintic_polynomial),
        x,
        y,
    )
    distinct_quintic_first_solution = sp.solve(
        distinct_quintic_first_full.coeffs(),
        list(dq_polynomial_coefficients.values()),
        dict=True,
        simplify=False,
    )[0]
    distinct_quintic_normal = sp.expand(
        distinct_quintic_polynomial.subs(
            distinct_quintic_first_solution
        )
    )
    distinct_quintic_second_full = sp.Poly(
        distinct_quintic_operator(
            distinct_quintic_operator(
                sp.expand(distinct_quintic_normal**2)
            )
        ),
        x,
        y,
    )
    dq_u = dq_polynomial_coefficients[(0, 4)]
    dq_v = dq_polynomial_coefficients[(0, 3)]
    dq_w = dq_polynomial_coefficients[(0, 2)]
    dq_a = dq_operator_coefficients[(3, 0)]
    dq_b = dq_operator_coefficients[(4, 0)]
    dq_c = dq_operator_coefficients[(5, 0)]
    assert (
        distinct_quintic_second_full.coeff_monomial(x**3 * y**2)
        == 480 * dq_u
    )
    assert (
        distinct_quintic_second_full.coeff_monomial(x**4).subs(dq_u, 0)
        == 48000 * dq_a**2
    )
    assert (
        distinct_quintic_second_full
        .coeff_monomial(x**3 * y)
        .subs({dq_u: 0, dq_a: 0})
        == 240 * dq_v
    )
    assert (
        distinct_quintic_second_full
        .coeff_monomial(x**3)
        .subs({dq_u: 0, dq_a: 0, dq_v: 0})
        == 80 * dq_w
    )
    assert (
        distinct_quintic_second_full
        .coeff_monomial(x**2)
        .subs({dq_u: 0, dq_a: 0, dq_v: 0, dq_w: 0})
        == 1296000 * dq_b**2
    )
    distinct_quintic_closure = {
        dq_u: 0,
        dq_a: 0,
        dq_v: 0,
        dq_w: 0,
        dq_b: 0,
    }
    assert (
        distinct_quintic_second_full
        .coeff_monomial(1)
        .subs(distinct_quintic_closure)
        == 3340800 * dq_c**2
    )
    distinct_quintic_final_polynomial = sp.expand(
        distinct_quintic_normal.subs(
            {**distinct_quintic_closure, dq_c: 0}
        )
    )
    assert sp.diff(distinct_quintic_final_polynomial, y, 2) == 0
    assert sp.diff(distinct_quintic_final_polynomial, x, y) == 0

    # The double-line x*y^4 tip closes through two nested weighted faces.
    xy4_operator_coefficients: dict[tuple[int, int], sp.Expr] = {
        (2, 0): sp.Integer(1)
    }
    for total_order in range(3, 9):
        for y_order in range(total_order + 1):
            x_order = total_order - y_order
            xy4_operator_coefficients[(x_order, y_order)] = sp.symbols(
                f"xy4_l{x_order}{y_order}"
            )
    xy4_polynomial_coefficients: dict[tuple[int, int], sp.Expr] = {}
    for total_degree in range(5):
        for y_degree in range(total_degree + 1):
            x_degree = total_degree - y_degree
            xy4_polynomial_coefficients[(x_degree, y_degree)] = sp.symbols(
                f"xy4_p{x_degree}{y_degree}"
            )

    def xy4_quintic_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in xy4_operator_coefficients.items()
            )
        )

    xy4_quintic_polynomial = x * y**4 + sum(
        coefficient * x**x_degree * y**y_degree
        for (x_degree, y_degree), coefficient
        in xy4_polynomial_coefficients.items()
    )
    xy4_quintic_first = sp.Poly(
        xy4_quintic_operator(xy4_quintic_polynomial),
        x,
        y,
    )
    xy4_quintic_first_solution = sp.solve(
        xy4_quintic_first.coeffs(),
        list(xy4_polynomial_coefficients.values()),
        dict=True,
        simplify=False,
    )[0]
    xy4_quintic_normal = sp.expand(
        xy4_quintic_polynomial.subs(xy4_quintic_first_solution)
    )
    xy4_quintic_second = sp.Poly(
        xy4_quintic_operator(
            xy4_quintic_operator(
                sp.expand(xy4_quintic_normal**2)
            )
        ),
        x,
        y,
    )
    xy4_l03 = xy4_operator_coefficients[(0, 3)]
    xy4_l12 = xy4_operator_coefficients[(1, 2)]
    xy4_l04 = xy4_operator_coefficients[(0, 4)]
    xy4_l05 = xy4_operator_coefficients[(0, 5)]
    xy4_l06 = xy4_operator_coefficients[(0, 6)]
    xy4_l07 = xy4_operator_coefficients[(0, 7)]
    xy4_l08 = xy4_operator_coefficients[(0, 8)]
    xy4_h = xy4_operator_coefficients[(1, 3)]
    xy4_j = xy4_operator_coefficients[(1, 4)]
    xy4_l21 = xy4_operator_coefficients[(2, 1)]
    assert xy4_quintic_second.coeff_monomial(y**5) == 1152 * xy4_l03
    assert (
        xy4_quintic_second
        .coeff_monomial(y**4)
        .subs(xy4_l03, 0)
        == 96 * (68 * xy4_l04 + xy4_l12**2)
    )
    assert sp.factor(
        xy4_quintic_second
        .coeff_monomial(x * y**2)
        .subs(
            {
                xy4_l03: 0,
                xy4_l04: -xy4_l12**2 / 68,
            }
        )
    ) == -sp.Rational(115200, 17) * xy4_l12**3
    xy4_initial_branch = {
        xy4_l03: 0,
        xy4_l12: 0,
        xy4_l04: 0,
    }
    assert (
        xy4_quintic_second.coeff_monomial(y**3).subs(
            xy4_initial_branch
        )
        == 26880 * xy4_l05
    )
    assert (
        xy4_quintic_second
        .coeff_monomial(y**2)
        .subs({**xy4_initial_branch, xy4_l05: 0})
        == 576 * (140 * xy4_l06 + 46 * xy4_h**2)
    )
    xy4_l06_branch = -sp.Rational(23, 70) * xy4_h**2
    xy4_l07_branch = (
        sp.Rational(16, 35) * xy4_h**2 * xy4_l21
        - sp.Rational(11, 14) * xy4_h * xy4_j
    )
    assert sp.expand(
        xy4_quintic_second
        .coeff_monomial(y)
        .subs(
            {
                **xy4_initial_branch,
                xy4_l05: 0,
                xy4_l06: xy4_l06_branch,
                xy4_l07: xy4_l07_branch,
            }
        )
    ) == 0
    xy4_second_after_h = {
        **xy4_initial_branch,
        xy4_l05: 0,
        xy4_h: 0,
        xy4_l06: 0,
        xy4_l07: 0,
    }
    assert (
        xy4_quintic_second.coeff_monomial(1).subs(
            xy4_second_after_h
        )
        == 24 * (6720 * xy4_l08 + 3216 * xy4_j**2)
    )

    xy4_first_face_polynomial = x * y**4 - 12 * xy4_h * x**2 * y

    def xy4_first_face_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2)
            + xy4_h * sp.diff(polynomial, x, y, 3)
            - sp.Rational(23, 70)
            * xy4_h**2
            * sp.diff(polynomial, y, 6)
        )

    xy4_first_face_third = sp.expand(xy4_first_face_polynomial**3)
    for _ in range(3):
        xy4_first_face_third = xy4_first_face_operator(
            xy4_first_face_third
        )
    assert (
        sp.Poly(xy4_first_face_third, x, y).coeff_monomial(y**3)
        == -553153536 * xy4_h**3
    )

    xy4_second_face_polynomial = x * y**4 - 12 * xy4_j * x**2

    def xy4_second_face_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2)
            + xy4_j * sp.diff(polynomial, x, y, 4)
            - sp.Rational(67, 140)
            * xy4_j**2
            * sp.diff(polynomial, y, 8)
        )

    xy4_second_face_third = sp.expand(xy4_second_face_polynomial**3)
    for _ in range(3):
        xy4_second_face_third = xy4_second_face_operator(
            xy4_second_face_third
        )
    assert (
        sp.factor(xy4_second_face_third)
        == -5430509568 * xy4_j**3
    )

    # The final double-line y^5 tip has six second-moment face ratios.
    y5_operator_coefficients: dict[tuple[int, int], sp.Expr] = {
        (2, 0): sp.Integer(1)
    }
    for total_order in range(3, 9):
        for y_order in range(total_order + 1):
            x_order = total_order - y_order
            y5_operator_coefficients[(x_order, y_order)] = sp.symbols(
                f"y5_l{x_order}{y_order}"
            )
    y5_polynomial_coefficients: dict[tuple[int, int], sp.Expr] = {}
    for total_degree in range(5):
        for y_degree in range(total_degree + 1):
            x_degree = total_degree - y_degree
            y5_polynomial_coefficients[(x_degree, y_degree)] = sp.symbols(
                f"y5_p{x_degree}{y_degree}"
            )

    def y5_quintic_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                coefficient
                * sp.diff(polynomial, x, x_order, y, y_order)
                for (x_order, y_order), coefficient
                in y5_operator_coefficients.items()
            )
        )

    y5_quintic_polynomial = y**5 + sum(
        coefficient * x**x_degree * y**y_degree
        for (x_degree, y_degree), coefficient
        in y5_polynomial_coefficients.items()
    )
    y5_quintic_first = sp.Poly(
        y5_quintic_operator(y5_quintic_polynomial),
        x,
        y,
    )
    y5_quintic_first_solution = sp.solve(
        y5_quintic_first.coeffs(),
        list(y5_polynomial_coefficients.values()),
        dict=True,
        simplify=False,
    )[0]
    y5_quintic_normal = sp.expand(
        y5_quintic_polynomial.subs(y5_quintic_first_solution)
    )
    y5_quintic_second = sp.Poly(
        y5_quintic_operator(
            y5_quintic_operator(
                sp.expand(y5_quintic_normal**2)
            )
        ),
        x,
        y,
    )
    y5_l03 = y5_operator_coefficients[(0, 3)]
    assert (
        y5_quintic_second.coeff_monomial(y**4)
        == 122400 * y5_l03**2
    )

    y5_a, y5_b, y5_z = sp.symbols("y5_a y5_b y5_z")
    y5_face_polynomial = (
        y**5
        + y5_z * x * y**3
        - 3 * (20 * y5_a + y5_b * y5_z) * x**2 * y
    )

    def y5_face_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2)
            + y5_b * sp.diff(polynomial, x, y, 2)
            + y5_a * sp.diff(polynomial, y, 4)
        )

    assert y5_face_operator(y5_face_polynomial) == 0
    y5_face_second = sp.expand(y5_face_polynomial**2)
    for _ in range(2):
        y5_face_second = y5_face_operator(y5_face_second)
    y5_residual_s = (
        24000 * y5_a**2
        - 1200 * y5_a * y5_b**2
        + 880 * y5_a * y5_b * y5_z
        + 20 * y5_a * y5_z**2
        - 60 * y5_b**3 * y5_z
        + y5_b**2 * y5_z**2
    )
    assert sp.expand(
        y5_face_second
        - 72
        * (
            16
            * y5_a
            * y5_z
            * (40 * y5_a + y5_b * y5_z)
            * x
            + y5_residual_s * y**2
        )
    ) == 0
    assert sp.expand(
        y5_residual_s.subs(y5_a, 0)
        - y5_b**2 * y5_z * (y5_z - 60 * y5_b)
    ) == 0
    assert sp.expand(
        y5_residual_s.subs(y5_z, 0)
        - 1200 * y5_a * (20 * y5_a - y5_b**2)
    ) == 0
    assert sp.expand(
        y5_residual_s.subs(
            y5_a,
            -y5_b * y5_z / 40,
        )
        + y5_b
        * y5_z
        * (y5_z**2 + 12 * y5_b * y5_z + 60 * y5_b**2)
        / 2
    ) == 0
    y5_face_third = sp.expand(y5_face_polynomial**3)
    for _ in range(3):
        y5_face_third = y5_face_operator(y5_face_third)
    assert sp.factor(
        y5_face_third.subs({y5_a: 0, y5_b: 1, y5_z: 60})
    ) == -1119744000 * y * (42 * x + 17 * y**2)
    assert sp.factor(
        y5_face_third.subs(
            {y5_a: sp.Rational(1, 20), y5_b: 1, y5_z: 0}
        )
    ) == -373248 * y * (4 * x - 27 * y**2)
    y5_t = sp.symbols("y5_t")
    y5_quadratic_branch = y5_face_third.subs(
        {
            y5_a: -y5_t / 40,
            y5_b: 1,
            y5_z: y5_t,
        }
    )
    assert sp.expand(
        sp.rem(
            y5_quadratic_branch,
            y5_t**2 + 12 * y5_t + 60,
            y5_t,
        )
        - (
            124416
            * y
            * (
                (51 * y5_t + 1080) * x
                - (259 * y5_t + 2220) * y**2
            )
        )
    ) == 0

    # The next r=3, deg(P)=5 row has eight top-form normal forms.
    c5_coefficients = sp.symbols("c5_0:6")
    cubic_quintic = sum(
        coefficient * x ** (5 - y_degree) * y**y_degree
        for y_degree, coefficient in enumerate(c5_coefficients)
    )
    triple_quintic_first = sp.Poly(
        sp.diff(cubic_quintic, x, 3),
        x,
        y,
    )
    assert sp.solve(
        triple_quintic_first.coeffs(),
        c5_coefficients[:3],
        dict=True,
    ) == [
        {
            c5_coefficients[0]: 0,
            c5_coefficients[1]: 0,
            c5_coefficients[2]: 0,
        }
    ]

    c5_a, c5_b, c5_c = sp.symbols("c5_a c5_b c5_c")
    double_cubic_quintic = (
        c5_a * x**5 + c5_b * x * y**4 + c5_c * y**5
    )
    double_cubic_quintic_second = sp.expand(
        double_cubic_quintic**2
    )
    for _ in range(2):
        double_cubic_quintic_second = sp.diff(
            double_cubic_quintic_second,
            x,
            2,
            y,
        )
    assert sp.factor(double_cubic_quintic_second) == (
        960
        * c5_a
        * x
        * y**2
        * (9 * c5_b * x + 5 * c5_c * y)
    )

    c5_d = sp.symbols("c5_d")
    squarefree_cubic_quintic = (
        c5_a * x**5
        - c5_d * x**4 * y
        + 2 * c5_d * x**3 * y**2
        - 2 * c5_d * x**2 * y**3
        + c5_d * x * y**4
        + c5_c * y**5
    )

    def squarefree_cubic_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x, 2, y)
            + sp.diff(polynomial, x, y, 2)
        )

    assert squarefree_cubic_operator(squarefree_cubic_quintic) == 0
    squarefree_cubic_quintic_second = sp.Poly(
        squarefree_cubic_operator(
            squarefree_cubic_operator(
                sp.expand(squarefree_cubic_quintic**2)
            )
        ),
        x,
        y,
    )
    squarefree_cubic_quintic_groebner = sp.groebner(
        squarefree_cubic_quintic_second.coeffs(),
        c5_a,
        c5_d,
        c5_c,
        order="lex",
    )
    assert [
        polynomial.as_expr()
        for polynomial in squarefree_cubic_quintic_groebner.polys
    ] == [
        c5_a * c5_d + c5_c * c5_d,
        5 * c5_a * c5_c - c5_c * c5_d,
        c5_d**2 + 5 * c5_c * c5_d,
    ]

    cubic_normal = (
        c0
        + c2 * y
        + c5 * (y**2 - 2 * C * x)
        + c9 * (y**3 - 6 * C * x * y + 6 * B * C * x)
    )
    assert operator(cubic_normal) == 0
    cubic_second = sp.Poly(
        operator(operator(sp.expand(cubic_normal**2))),
        x,
        y,
    )
    assert cubic_second.coeff_monomial(y**2) == 144 * C**2 * c9**2
    assert sp.factor(
        cubic_second.as_expr().subs(c9, 0)
    ) == 16 * C**2 * c5**2

    # A general heat-harmonic polynomial is the finite heat evolution of
    # arbitrary boundary data p(y).
    boundary_coefficients = sp.symbols("h0:7")
    boundary = sum(
        coefficient * y**degree
        for degree, coefficient in enumerate(boundary_coefficients)
    )
    heat_polynomial = 0
    derivative = boundary
    heat_order = 0
    while derivative != 0:
        heat_polynomial += (
            (-C * x) ** heat_order
            / sp.factorial(heat_order)
            * derivative
        )
        derivative = sp.diff(derivative, y, 2)
        heat_order += 1
    heat_polynomial = sp.expand(heat_polynomial)

    def heat_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x) + C * sp.diff(polynomial, y, 2)
        )

    assert heat_operator(heat_polynomial) == 0
    heat_square = sp.factor(
        heat_operator(heat_operator(sp.expand(heat_polynomial**2)))
        - 4 * C**2 * sp.diff(heat_polynomial, y, 2) ** 2
    )
    assert heat_square == 0

    E = sp.symbols("E")

    def transverse_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            C * sp.diff(polynomial, y, 2)
            + E * sp.diff(polynomial, y, 3)
        )

    separated_polynomial = 0
    transverse_derivative = boundary
    separated_order = 0
    while transverse_derivative != 0:
        separated_polynomial += (
            (-x) ** separated_order
            / sp.factorial(separated_order)
            * transverse_derivative
        )
        transverse_derivative = transverse_operator(transverse_derivative)
        separated_order += 1
    separated_polynomial = sp.expand(separated_polynomial)

    def separated_operator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            sp.diff(polynomial, x) + transverse_operator(polynomial)
        )

    assert separated_operator(separated_polynomial) == 0
    separated_second = sp.Poly(
        separated_operator(
            separated_operator(sp.expand(separated_polynomial**2))
        ).subs(x, 0),
        y,
    )
    top_boundary = boundary_coefficients[6]
    assert (
        separated_second.coeff_monomial(y**8)
        == 3600 * C**2 * top_boundary**2
    )

    print("PASS heat class: first moment gives a=b=0 and d=-2*C*c")
    print("PASS heat class: Lambda^2(P^2)=16*C^2*c^2")
    print("PASS coupled cubic: quadratic-P second moment is unchanged")
    print("PASS cubic P: exact all-operator-order three-branch obstruction")
    print("PASS quartic P: complete 7-jet three-branch obstruction")
    print("PASS lowest order two: cubic double-line third-moment closure")
    print("PASS lowest order two: cubic distinct-line second-moment closure")
    print("PASS lowest order two: quartic distinct-line second-moment closure")
    print("PASS lowest order two: quartic double-line finite residual sextic")
    print("PASS lowest order two: full quartic residual dies at moment four")
    print("PASS lowest order three: triple-root quartic cell closes")
    print("PASS lowest order three: double-root quartic cell closes")
    print("PASS lowest order three: squarefree quartic cell closes")
    print("PASS first degree-five row: three quadratic-leading top forms")
    print("PASS quadratic-leading quintic: distinct-root branch closes")
    print("PASS quadratic-leading quintic: x*y^4 branch closes")
    print("PASS quadratic-leading quintic: y^5 branch closes")
    print("PASS cubic-leading quintic: eight top-form normal forms")
    print("PASS drift--diffusion: cubic top coefficients vanish at moment two")
    print("PASS heat class: degree-six heat square identity is exact")
    print("PASS separated drift: degree-six D^2+D^3 leading obstruction is exact")
    print("SCOPE: the all-order mixed cutoff is the written derivative count")


if __name__ == "__main__":
    main()
