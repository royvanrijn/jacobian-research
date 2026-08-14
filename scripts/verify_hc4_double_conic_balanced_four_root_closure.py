#!/usr/bin/env python3
"""Verify the exact closure of the double-conic (3,3,2,2) row.

The normal-layer checker leaves this row generically empty over QQ(cr).
Here a small coefficient calculation removes every exceptional fiber.  The
only assumption on the cross-ratio is cr != 0,1, as required for four
distinct support points.
"""

from __future__ import annotations

import sympy as sp

from verify_hc4_double_conic_normal_layers import (
    canonical_lift,
    cross_ratio,
    g_coefficients,
    general_cubic,
    q,
    s,
    t,
    x,
    y,
    z,
)


def coefficient(polynomial: sp.Poly, exponents: tuple[int, int, int]) -> sp.Expr:
    """Return the coefficient of x^a*y^b*z^c."""

    a, b, c = exponents
    return polynomial.coeff_monomial(x**a * y**b * z**c)


def main() -> None:
    cr = cross_ratio
    g = g_coefficients
    binary_form = s**3 * t**3 * (s - t) ** 2 * (s - cr * t) ** 2
    h5 = canonical_lift(binary_form) + q * general_cubic
    determinant = sp.Poly(
        sp.expand(sp.hessian(h5, (x, y, z)).det()), x, y, z
    )

    line_x = coefficient(determinant, (5, 0, 4))
    line_y = coefficient(determinant, (4, 1, 4))
    line_z = coefficient(determinant, (4, 0, 5))
    remainder = sp.Poly(
        sp.expand(
            determinant.as_expr()
            - q**4 * (line_x * x + line_y * y + line_z * z)
        ),
        x,
        y,
        z,
    )

    # The x-endpoint normal layers successively force
    # g9=0, g8=-1, and g7=2*cr+2.  The residual x-line coefficient then
    # vanishes identically.
    assert coefficient(remainder, (9, 0, 0)) == 32 * g[9] ** 3
    x_first = {g[9]: 0}
    assert sp.expand(
        coefficient(remainder, (7, 1, 1)).subs(x_first)
        - 12 * (g[8] + 1) ** 3
    ) == 0
    x_second = x_first | {g[8]: -1}
    assert sp.expand(
        coefficient(remainder, (6, 1, 2)).subs(x_second)
        + 144 * (g[7] - 2 * cr - 2) ** 2
    ) == 0
    x_endpoint = x_second | {g[7]: 2 * cr + 2}
    assert sp.expand(line_x.subs(x_endpoint)) == 0

    # At the z-endpoint, cr != 0 gives the symmetric three-step chain.
    assert coefficient(remainder, (0, 0, 9)) == 32 * g[0] ** 3
    z_first = {g[0]: 0}
    assert sp.expand(
        coefficient(remainder, (1, 1, 7)).subs(z_first)
        - 12 * (g[1] + cr**2) ** 3
    ) == 0
    z_second = z_first | {g[1]: -(cr**2)}
    assert sp.expand(
        coefficient(remainder, (2, 1, 6)).subs(z_second)
        + 144 * cr**2 * (g[4] - 2 * cr**2 - 2 * cr) ** 2
    ) == 0
    z_endpoint = z_second | {g[4]: 2 * cr**2 + 2 * cr}
    assert sp.expand(line_z.subs(z_endpoint)) == 0

    endpoint_substitution = x_endpoint | z_endpoint
    assert sp.expand(line_x.subs(endpoint_substitution)) == 0
    assert sp.expand(line_z.subs(endpoint_substitution)) == 0

    u, v = sp.symbols("u v")
    middle_coordinates = {
        g[5]: u - cr**2 - 4 * cr - 1,
        g[6]: v + 2 * cr + 2,
    }

    def middle(expression: sp.Expr) -> sp.Expr:
        return sp.expand(
            expression.subs(endpoint_substitution).subs(middle_coordinates)
        )

    assert sp.expand(middle(line_y) - 16 * u**3) == 0
    assert sp.expand(
        middle(coefficient(remainder, (4, 5, 0)))
        + 16 * (3 * u - v**2)
    ) == 0
    assert sp.expand(
        middle(coefficient(remainder, (4, 4, 1)))
        + 16 * v * (2 * u - v**2)
    ) == 0

    # Two redundant neighboring layers are useful regression checks.
    assert sp.expand(
        middle(coefficient(remainder, (4, 3, 2)))
        + 32 * u * (u - v**2)
    ) == 0
    assert sp.expand(
        middle(coefficient(remainder, (4, 2, 3))) - 32 * u**2 * v
    ) == 0

    print("PASS: endpoint normal layers force both endpoint line coefficients to zero")
    print("PASS: middle layers are B=16*u^3, 3*u-v^2, and v*(2*u-v^2)")
    print("THEOREM: the complete four-root partition (3,3,2,2) is empty")


if __name__ == "__main__":
    main()
