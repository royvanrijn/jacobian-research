#!/usr/bin/env python3
"""Exact verification of the prime-to-characteristic realization theorem.

The theorem uses two explicit maps:

* an integral cubic map with determinant 1, generic separable degree 3,
  and a complete three-point split fiber over every field;
* a ten-monomial characteristic-3 quartic with determinant 1, generic
  separable degree 4, and a visible rational collision.

The irreducibility arguments are structural and are documented in the
companion Markdown note. This script checks every polynomial identity used
by those arguments.
"""

from __future__ import annotations

import sympy as sp


x, y, z, T = sp.symbols("x y z T")


def numerator(expr: sp.Expr) -> sp.Expr:
    """Return the expanded numerator of a rational expression."""
    return sp.expand(sp.cancel(expr).as_numer_denom()[0])


def poly_mod(expr: sp.Expr, variables: tuple[sp.Symbol, ...], p: int) -> sp.Poly:
    return sp.Poly(sp.expand(expr), *variables, modulus=p)


def verify_integral_cubic() -> None:
    beta = x * z - 2 * y
    gamma = -x * z + 3 * y

    a = x
    b = 1 + x * beta
    c = 1 + x * gamma
    d = -(beta + gamma) - x * beta * gamma
    e = -z - 2 * beta**2 - 4 * beta * gamma - 2 * x * beta**2 * gamma

    P = sp.expand(a * c)
    Q = sp.expand(a * e + b * d)
    R = sp.expand(b * e)
    F = (P, Q, R)

    # The map is Keller over Z, hence after base change to every field.
    jacobian = sp.Matrix(F).jacobian((x, y, z)).det()
    assert sp.expand(jacobian - 1) == 0

    # Normalized linear-times-quadratic factorization.
    assert sp.expand(a * d + b * c - 1) == 0
    assert sp.expand(a**2 * e - a * b * d + b**2 * c - 1) == 0
    cubic = P * T**3 + T**2 + Q * T + R
    assert sp.expand(cubic - (a * T + b) * (c * T**2 + d * T + e)) == 0

    # The full coefficient-resultant map has Jacobian -Res^2.
    aa, bb, cc, dd, ee = sp.symbols("a b c d e")
    resultant = aa**2 * ee - aa * bb * dd + bb**2 * cc
    theta = (
        aa * cc,
        aa * dd + bb * cc,
        aa * ee + bb * dd,
        bb * ee,
        resultant,
    )
    theta_jacobian = sp.Matrix(theta).jacobian((aa, bb, cc, dd, ee)).det()
    assert sp.expand(theta_jacobian + resultant**2) == 0

    # The displayed source chart is an integral isomorphism A^3_Z -> X.
    inverse_y = (
        2 * aa * dd**2
        + aa * cc * ee
        + 6 * aa * bb * dd**2
        + 3 * aa * bb * cc * ee
        - 4 * aa * ee
        - bb * dd
    )
    inverse_z = (
        4 * dd**2
        + 2 * cc * ee
        + 12 * bb * dd**2
        + 6 * bb * cc * ee
        - 9 * ee
    )
    forward_substitution = {aa: a, bb: b, cc: c, dd: d, ee: e}
    assert sp.expand(inverse_y.subs(forward_substitution) - y) == 0
    assert sp.expand(inverse_z.subs(forward_substitution) - z) == 0

    relation_1 = aa * dd + bb * cc - 1
    relation_2 = aa**2 * ee - aa * bb * dd + bb**2 * cc - 1
    groebner = sp.groebner(
        (relation_1, relation_2), ee, dd, cc, bb, aa, order="lex", domain=sp.ZZ
    )
    reverse_beta = aa * inverse_z - 2 * inverse_y
    reverse_gamma = -aa * inverse_z + 3 * inverse_y
    reverse = (
        1 + aa * reverse_beta,
        1 + aa * reverse_gamma,
        -(reverse_beta + reverse_gamma) - aa * reverse_beta * reverse_gamma,
        -inverse_z
        - 2 * reverse_beta**2
        - 4 * reverse_beta * reverse_gamma
        - 2 * aa * reverse_beta**2 * reverse_gamma,
    )
    for expected, reconstructed in zip((bb, cc, dd, ee), reverse):
        remainder = groebner.reduce(sp.expand(expected - reconstructed))[1]
        assert remainder == 0

    # Rational reconstruction from the marked root tau=-b/a.
    tau = -b / a
    delta = 3 * P * tau**2 + 2 * tau + Q
    assert numerator(a * delta - 1) == 0
    assert numerator(P * delta**2 - 2 * delta - tau - y) == 0
    assert numerator(delta * (2 * P * delta**2 - 5 * delta - 3 * tau) - z) == 0
    assert numerator(cubic.subs(T, tau)) == 0
    assert numerator(sp.diff(cubic, T).subs(T, tau) - delta) == 0

    # Complete split fiber over (0,1,0), valid over every field.
    target = (0, 1, 0)
    points = ((-1, 3, -8), (0, -1, 16), (1, -2, -5))
    factor_tuples = (
        (-1, -1, 0, -1, 0),
        (0, 1, 1, 1, 0),
        (1, 0, 0, 1, 1),
    )
    for point, factor_tuple in zip(points, factor_tuples):
        substitution = dict(zip((x, y, z), point))
        assert tuple(int(value.subs(substitution)) for value in F) == target
        assert tuple(
            int(value.subs(substitution)) for value in (a, b, c, d, e)
        ) == factor_tuple

    # The binary target cubic is TS(T+S), so the three marked linear factors
    # are T, S, and T+S in every characteristic.
    U, V = sp.symbols("U V")
    binary_cubic = U**2 * V + U * V**2
    assert sp.expand(binary_cubic - U * V * (U + V)) == 0


def verify_characteristic_three_quartic() -> None:
    p = 3

    A = -y - x**2 * y**2
    B = (
        x * y
        - z
        + x * z**2
        + x**5 * y**3
        - x**6 * y**3 * z
        + x**7 * y**3 * z**2
    )
    C = x + x**2 * z
    F = (A, B, C)

    jacobian = sp.Matrix(F).jacobian((x, y, z)).det()
    assert poly_mod(jacobian - 1, (x, y, z), p).is_zero
    assert sum(len(poly_mod(value, (x, y, z), p).terms()) for value in F) == 10

    # Visible prime-field collision.
    for point in ((1, 0, 0), (-1, 0, -1)):
        substitution = dict(zip((x, y, z), point))
        values = tuple(int(value.subs(substitution)) % p for value in F)
        assert values == (0, 0, 1)

    # Weight-redistributed suspension certificate.
    u = 1 + x**2 * y
    gamma = 1 + x * z
    W = u * gamma
    H = 2 * W**2 + W**4
    H_prime = W + W**3  # dH/dW in characteristic 3
    S = H_prime + gamma
    core_T = W * S - H

    assert poly_mod(A * C**2 - core_T, (x, y, z), p).is_zero
    assert poly_mod(B * C - S, (x, y, z), p).is_zero

    # Three determinant factors in the suspension square.
    source_jacobian = sp.Matrix((W, gamma, C)).jacobian((x, y, z)).det()
    assert poly_mod(source_jacobian - x**3 * gamma**2, (x, y, z), p).is_zero

    w, g, c0 = sp.symbols("w g c0")
    h = 2 * w**2 + w**4
    h_prime = w + w**3
    s = h_prime + g
    t0 = w * s - h
    core_jacobian = sp.Matrix((t0, s, c0)).jacobian((w, g, c0)).det()
    assert poly_mod(core_jacobian - g, (w, g, c0), p).is_zero

    tt, ss, cc0 = sp.symbols("tt ss cc0")
    target_jacobian = sp.Matrix((tt / cc0**2, ss / cc0, cc0)).jacobian(
        (tt, ss, cc0)
    ).det()
    assert sp.cancel(target_jacobian - cc0**-3) == 0

    # Inverse quartic and its marked-root derivative.
    inverse_quartic = W**4 + 2 * W**2 - B * C * W + A * C**2
    assert poly_mod(inverse_quartic, (x, y, z), p).is_zero
    inverse_derivative = 4 * W**3 + 4 * W - B * C
    assert poly_mod(inverse_derivative + gamma, (x, y, z), p).is_zero

    # Denominator-cleared reconstruction identities.
    reconstructed_gamma = B * C - (W + W**3)
    assert poly_mod(reconstructed_gamma - gamma, (x, y, z), p).is_zero
    assert poly_mod(x * reconstructed_gamma - C, (x, y, z), p).is_zero
    assert poly_mod(u * reconstructed_gamma - W, (x, y, z), p).is_zero
    assert poly_mod(x**2 * y - (u - 1), (x, y, z), p).is_zero
    assert poly_mod(x * z - (gamma - 1), (x, y, z), p).is_zero

    # Structural irreducibility checks: over k(B,C), the inverse equation is
    # degree one in A, and it has exact W-degree four.
    AA, BB, CC, WW = sp.symbols("A B C W")
    generic_quartic = WW**4 + 2 * WW**2 - BB * CC * WW + AA * CC**2
    generic_poly = sp.Poly(
        generic_quartic, AA, WW, domain=sp.GF(p).frac_field(BB, CC)
    )
    assert generic_poly.degree(AA) == 1
    assert (
        sp.Poly(
            generic_quartic, WW, domain=sp.GF(p).frac_field(AA, BB, CC)
        ).degree()
        == 4
    )


def main() -> None:
    verify_integral_cubic()
    print("integral cubic: all exact identities passed")
    verify_characteristic_three_quartic()
    print("characteristic-3 quartic: all exact identities passed")
    print("prime-to-characteristic realization verification passed")


if __name__ == "__main__":
    main()
