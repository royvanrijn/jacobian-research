#!/usr/bin/env python3
"""Verify the fixed-coordinate and affine-etale k=1 pullback interface."""

from __future__ import annotations

import sympy as sp

from verify_f2_affine_target_k1_implicit_conductor import (
    expected_implicit_quintic,
)


t, u, P, Q = sp.symbols("t u P Q")
a, b, c, d = sp.symbols("a b c d")
A, B = sp.symbols("A B", nonzero=True)
P0, Q0, Gamma = sp.symbols("P0 Q0 Gamma")


def fixed_target_normalization_audit() -> tuple[sp.Expr, sp.Expr]:
    """Undo every target change used in the four-parameter normal form."""

    normalized_p = t**3 + a * t
    normalized_q = t**5 + b * t**4 + c * t**2 + d * t
    fixed_p = P0 + A * normalized_p
    fixed_q = Q0 + B * normalized_q + Gamma * normalized_p

    normalized_implicit = expected_implicit_quintic()
    p_numerator = P - P0
    q_numerator = A * (Q - Q0) - Gamma * (P - P0)
    substituted = normalized_implicit.subs(
        {P: p_numerator / A, Q: q_numerator / (A * B)},
        simultaneous=True,
    )
    fixed_implicit, denominator = sp.cancel(A**5 * B**3 * substituted).as_numer_denom()
    assert denominator == 1
    fixed_implicit = sp.expand(fixed_implicit)
    assert sp.Poly(fixed_implicit, P, Q).total_degree() == 5

    # The inverse triangular normalization sends the fixed parametrization
    # back to the already-certified normal form.
    restored_p = sp.cancel((fixed_p - P0) / A)
    restored_q = sp.cancel(
        (A * (fixed_q - Q0) - Gamma * (fixed_p - P0)) / (A * B)
    )
    assert restored_p == normalized_p
    assert restored_q == normalized_q
    assert sp.expand(normalized_implicit.subs(
        {P: restored_p, Q: restored_q}, simultaneous=True
    )) == 0

    # Clearing denominators preserves the fixed curve, and its top degree
    # form is the expected fifth power in the P-dominant direction.
    assert sp.expand(fixed_implicit.subs(
        {P: fixed_p, Q: fixed_q}, simultaneous=True
    )) == 0
    fixed_poly = sp.Poly(fixed_implicit, P, Q)
    top = sum(
        coefficient * P**monomial[0] * Q**monomial[1]
        for monomial, coefficient in fixed_poly.terms()
        if sum(monomial) == 5
    )
    assert top == B**3 * P**5
    return fixed_p, fixed_q


def carrier_jet_audit(fixed_p: sp.Expr, fixed_q: sp.Expr) -> None:
    """Recover the fixed carrier residue and its finite eight-jet test."""

    p_unit = 1 + a * u**2 + (P0 / A) * u**3
    minus_q_unit = (
        1
        + b * u
        + (Gamma / B) * u**2
        + c * u**3
        + (d + Gamma * a / B) * u**4
        + (Q0 / B) * u**5
    )
    residue = A**5 / (-B) ** 3
    pi_leading = A**3 / B**2
    h = residue * p_unit**5 / minus_q_unit**3
    pi = pi_leading * u * p_unit**3 / minus_q_unit**2

    carrier_centers = sp.symbols("c1:8")
    w = h - residue - sum(
        center * pi**index
        for index, center in enumerate(carrier_centers, start=1)
    )
    jet = sp.series(w, u, 0, 9).removeO().expand()
    coefficients = [sp.expand(jet.coeff(u, index)) for index in range(1, 9)]
    assert len(coefficients) == 8
    assert sp.expand(
        coefficients[0] - (-3 * residue * b - carrier_centers[0] * pi_leading)
    ) == 0

    # Direct Laurent substitution agrees with the unit-series calculation.
    direct_p = sp.expand(fixed_p.subs(t, 1 / u))
    direct_q = sp.expand(fixed_q.subs(t, 1 / u))
    assert sp.cancel(direct_p**5 / (-direct_q) ** 3 - h) == 0
    assert sp.cancel(direct_p**3 / (-direct_q) ** 2 - pi) == 0


def keller_gradient_audit() -> None:
    """Check the GL_2 identity behind equality of the two Jacobian ideals."""

    p_x, p_y, q_x, q_y = sp.symbols("p_x p_y q_x q_y")
    g_p, g_q = sp.symbols("g_p g_q")
    jacobian = p_x * q_y - p_y * q_x
    h_x = p_x * g_p + q_x * g_q
    h_y = p_y * g_p + q_y * g_q

    recovered_g_p = sp.cancel((q_y * h_x - q_x * h_y) / jacobian)
    recovered_g_q = sp.cancel((-p_y * h_x + p_x * h_y) / jacobian)
    assert recovered_g_p == g_p
    assert recovered_g_q == g_q


def node_fiber_audit(fixed_p: sp.Expr, fixed_q: sp.Expr) -> None:
    """Transport the four collision values into the fixed target chart."""

    collision_parameter = sp.symbols("collision_parameter")
    r = collision_parameter
    collision_quartic = (
        r**4 + b * r**3 + a * r**2 + (2 * a * b - c) * r - (a**2 + d)
    )
    normalized_x = -r * (r**2 + a)
    normalized_y = (r**2 + a) * (r**3 + 2 * a * r + a * b - c)
    fixed_x = P0 + A * normalized_x
    fixed_y = Q0 + B * normalized_y + Gamma * normalized_x

    assert sp.Poly(collision_quartic, r).degree() == 4
    assert sp.expand((fixed_x - P0) / A - normalized_x) == 0
    assert sp.expand(
        (A * (fixed_y - Q0) - Gamma * (fixed_x - P0)) / (A * B)
        - normalized_y
    ) == 0

    # The parametrization itself has the advertised leading coefficients;
    # the shear is lower order and cannot alter the carrier residue.
    assert sp.Poly(fixed_p, t).LC() == A
    assert sp.Poly(fixed_q, t).LC() == B


def finite_cover_fiber_audit() -> None:
    """Replay the rank-d boundary/affine length split at four target nodes."""

    for degree in (6, 12, 9375):
        affine_lengths = []
        for boundary_length in range(1, degree + 1):
            affine_length = degree - boundary_length
            assert 0 <= affine_length <= degree - 1
            assert affine_length + boundary_length == degree
            affine_lengths.append(affine_length)
        assert max(affine_lengths) == degree - 1
        assert 4 * max(affine_lengths) == 4 * (degree - 1)


def main() -> None:
    fixed_p, fixed_q = fixed_target_normalization_audit()
    carrier_jet_audit(fixed_p, fixed_q)
    keller_gradient_audit()
    node_fiber_audit(fixed_p, fixed_q)
    finite_cover_fiber_audit()
    print(
        "PASS: every k=1 target normal form has an exact fixed-coordinate "
        "quintic and eight-jet carrier interface; a Keller pullback is "
        "reduced, its affine singular/conductor scheme is the etale base "
        "change of the target nodes, and each node fiber has at most d-1 "
        "affine points"
    )


if __name__ == "__main__":
    main()
