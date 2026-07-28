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
    print("PASS drift--diffusion: cubic top coefficients vanish at moment two")
    print("PASS heat class: degree-six heat square identity is exact")
    print("PASS separated drift: degree-six D^2+D^3 leading obstruction is exact")
    print("SCOPE: the all-order mixed cutoff is the written derivative count")


if __name__ == "__main__":
    main()
