#!/usr/bin/env python3
"""Exact reductions for extending the Q(sqrt(-31)) fixed-quintic Hasse row.

The script verifies two routes through the certified point:

1. fixing its normalized quadratic factor parameters produces a genus-two
   cube curve; and
2. varying the two field generators affinely reduces the trace normalization
   to one explicit quadric and the cube condition to one explicit quartic.

These identities are a frontier reduction, not an infinitude theorem.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt

import sympy as sp


def main() -> None:
    A, V, W, R, K, P = sp.symbols("A V W R K P")
    M = (
        3 * A**8
        - 24 * A**6 * V
        - 50 * A**6
        + 70 * A**4 * V**2
        + 270 * A**4 * V
        + 56 * A**4 * W
        + 275 * A**4
        - 76 * A**2 * V**3
        - 510 * A**2 * V**2
        - 288 * A**2 * V * W
        - 750 * A**2 * V
        - 360 * A**2 * W
        - 500 * A**2
        + 27 * V**4
        + 270 * V**3
        + 216 * V**2 * W
        + 675 * V**2
        + 1080 * V * W
        + 432 * W**2
    )

    specialized = sp.factor(
        M.subs(V, (A**2 + 31 * R**2) / 4).subs({A: -8, R: 2})
    )
    assert sp.factor(
        specialized - 16 * (27 * W**2 - 8254 * W + 617811)
    ) == 0
    assert 31 * specialized.subs(W, 125) == 1984**2

    parametrized_w = sp.factor(
        (125 * K**2 - 3968 * K - 2419984) / (K**2 - 13392)
    )
    parametrized_y = 1984 + K * (parametrized_w - 125)
    assert sp.factor(parametrized_y**2 - 31 * specialized.subs(W, parametrized_w)) == 0

    cube_curve = 1984 * (27 * P**6 - 8254 * P**3 + 617811)
    conic_discriminant = sp.factor(
        3968**2
        - 4 * (P**3 - 125) * (2419984 - 13392 * P**3)
    )
    assert sp.factor(conic_discriminant - cube_curve) == 0
    assert cube_curve.subs(P, 5) == 3968**2

    # Affine variation in
    # Q(sqrt(-31)) x Q(theta), theta^3+8 theta^2+12 theta+8=0.
    u, v, w, s = sp.symbols("u v w s")
    x = sp.symbols("x")
    cubic = x**3 + 8 * x**2 + 12 * x + 8

    def cubic_trace(poly: sp.Expr) -> sp.Expr:
        remainder = sp.rem(sp.Poly(sp.expand(poly), x), sp.Poly(cubic, x)).as_expr()
        return sp.expand(
            3 * remainder.coeff(x, 0)
            - 8 * remainder.coeff(x, 1)
            + 40 * remainder.coeff(x, 2)
        )

    quadratic_p2 = 2 * (v**2 - 31 * u**2)
    quadratic_p4 = 2 * (v**4 - 186 * v**2 * u**2 + 961 * u**4)
    eta3 = w * x + s
    trace = sp.expand(2 * v + cubic_trace(eta3))
    p2 = sp.expand(quadratic_p2 + cubic_trace(eta3**2))
    p4 = sp.expand(quadratic_p4 + cubic_trace(eta3**4))

    trace_zero_s = (8 * w - 2 * v) / 3
    assert sp.factor(trace.subs(s, trace_zero_s)) == 0
    assert sp.factor(
        p2.subs(s, trace_zero_s)
        - 10
        - 2 * (5 * v**2 + 28 * w**2 - 93 * u**2 - 15) / 3
    ) == 0

    cube_numerator = (
        25947 * u**4
        - 5022 * u**2 * v**2
        + 35 * v**4
        + 672 * v**2 * w**2
        + 1504 * v * w**3
        + 2352 * w**4
        - 675
    )
    assert sp.factor(
        (50 - p4.subs(s, trace_zero_s)) / 16
        + cube_numerator / 216
    ) == 0

    known = {u: 1, v: 4, w: 1}
    assert (5 * v**2 + 28 * w**2 - 93 * u**2).subs(known) == 15
    assert (-cube_numerator / 216).subs(known) == 125

    rational_points: list[tuple[Fraction, Fraction]] = []
    height_bound = 600
    for denominator in range(1, height_bound + 1):
        for numerator in range(-height_bound, height_bound + 1):
            if gcd(abs(numerator), denominator) != 1:
                continue
            pi_value = Fraction(numerator, denominator)
            square_value = 1984 * (
                27 * pi_value**6
                - 8254 * pi_value**3
                + 617811
            )
            if square_value < 0:
                continue
            root_numerator = isqrt(square_value.numerator)
            root_denominator = isqrt(square_value.denominator)
            if (
                root_numerator**2 == square_value.numerator
                and root_denominator**2 == square_value.denominator
            ):
                rational_points.append(
                    (
                        pi_value,
                        Fraction(root_numerator, root_denominator),
                    )
                )
    assert rational_points == [(Fraction(5), Fraction(3968))]

    print("PASS: the fixed (A,R)=(-8,2) conic parametrization is exact")
    print("PASS: its remaining cube condition is the displayed genus-two curve")
    print("PASS: its only rational x-coordinate of height <=600 is Pi=5")
    print("PASS: affine generator variation reduces to an explicit rational quadric")
    print("PASS: its fourth-moment cube condition is the displayed quartic")


if __name__ == "__main__":
    main()
