#!/usr/bin/env python3
"""Exact checks for finite etale algebras as Keller fibers.

The checker has three independent layers:

1. direct Jacobian and inverse-polynomial checks for translated seeds;
2. quotient-ring reconstruction, using an explicit Bezout inverse of E';
3. an exact target-scaling audit of the displayed Berend--Bilu map.
"""
from __future__ import annotations

import sympy as sp

x, y, z, S, ell = sp.symbols("x y z S ell")


def quadratic_gauge_map(G: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the normalized determinant-minus-two gauge for a rooted G."""
    poly = sp.Poly(sp.expand(G), S, domain=sp.QQ)
    degree = poly.degree()
    coeff = {k: poly.coeff_monomial(S**k) for k in range(1, degree + 1)}
    g1 = coeff[1]
    g3 = coeff[3]
    assert g1 != 0 and g3 != 0 and coeff[degree] != 0

    t = 1 + x * y
    q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
    pi = t * q
    b = y + 3 * (g3 / g1) * x * q + 2 * (coeff.get(2, 0) / g1) * t * q
    c = x * (5 - 3 * t) - (g3 / g1) * x**3 * z
    for k in range(4, degree + 1):
        b += k * (coeff[k] / g1) * t**2 * x ** (k - 2) * q**k
        c -= (k - 2) * (coeff[k] / g1) * (x * q) ** k
    return tuple(sp.cancel(component) for component in (pi, b, c))


def quotient_reduce(expression: sp.Expr, modulus: sp.Expr) -> sp.Expr:
    """Reduce a rational expression in QQ(S) modulo a squarefree modulus."""
    expression = sp.cancel(expression)
    numerator, denominator = expression.as_numer_denom()
    modulus_poly = sp.Poly(modulus, S, domain=sp.QQ)
    numerator_poly = sp.Poly(sp.expand(numerator), S, domain=sp.QQ).rem(modulus_poly)
    denominator_poly = sp.Poly(sp.expand(denominator), S, domain=sp.QQ).rem(modulus_poly)
    denominator_inverse = sp.invert(denominator_poly, modulus_poly)
    return sp.Poly(
        numerator_poly * denominator_inverse,
        S,
        domain=sp.QQ,
    ).rem(modulus_poly).as_expr()


def check_scheme_reconstruction(P: sp.Expr, a: sp.Expr) -> None:
    """Check both sides of the reconstruction in QQ[S]/(P(a+S))."""
    translated = sp.expand(P.subs(S, S + a))
    G = sp.expand(translated - P.subs(S, a))
    poly = sp.Poly(G, S, domain=sp.QQ)
    degree = poly.degree()
    coeff = {k: poly.coeff_monomial(S**k) for k in range(1, degree + 1)}
    g1 = coeff[1]
    g3 = coeff[3]
    assert g1 != 0 and g3 != 0 and coeff[degree] != 0

    pi = sp.Integer(1)
    b_target = sp.Integer(0)
    c_target = sp.cancel(-2 * P.subs(S, a) / g1)
    E = sp.expand(G - g1 * (b_target * S**2 + c_target) / 2)
    assert sp.expand(E - translated) == 0

    E_poly = sp.Poly(E, S, domain=sp.QQ)
    derivative_poly = E_poly.diff()
    assert sp.gcd(E_poly, derivative_poly).degree() == 0

    # Bezout supplies the unit used by the scheme-theoretic reconstruction.
    derivative_inverse = sp.invert(derivative_poly, E_poly).as_expr()
    assert quotient_reduce(sp.diff(E, S) * derivative_inverse - 1, E) == 0

    beta = sp.cancel((sp.diff(G, S) / g1 - 1 - pi * S**2) / S)
    d = quotient_reduce(sp.diff(E, S) / g1, E)
    d_inverse = quotient_reduce(1 / d, E)
    Q = quotient_reduce(b_target - beta, E)
    t_bar = d_inverse
    x_bar = quotient_reduce(S * d_inverse, E)
    y_bar = quotient_reduce(Q - pi * S, E)
    q_bar = quotient_reduce(pi * d, E)
    z_bar = quotient_reduce(
        (q_bar - (g1 / g3) * y_bar**2 * (1 + 3 * t_bar)) / t_bar**2,
        E,
    )

    # Reconstruct the rational chart globally in the quotient algebra.
    assert quotient_reduce(1 + x_bar * y_bar - t_bar, E) == 0
    assert quotient_reduce(
        t_bar**2 * z_bar
        + (g1 / g3) * y_bar**2 * (1 + 3 * t_bar)
        - q_bar,
        E,
    ) == 0
    assert quotient_reduce(x_bar / t_bar - S, E) == 0
    assert quotient_reduce(y_bar + x_bar * q_bar - Q, E) == 0

    # Substitution into all three map coordinates returns the target in R.
    outputs = [
        quotient_reduce(
            component.subs({x: x_bar, y: y_bar, z: z_bar}),
            E,
        )
        for component in quadratic_gauge_map(G)
    ]
    assert quotient_reduce(outputs[0] - pi, E) == 0
    assert quotient_reduce(outputs[1] - b_target, E) == 0
    assert quotient_reduce(outputs[2] - c_target, E) == 0


def check_polynomial_to_fiber_transfer() -> None:
    examples = [
        S**3 - S - 1,
        2 * S**4 - S**3 - S**2 + S + 1,
        (S**3 - 19) * (S**2 + S + 1),
    ]
    for P in examples:
        p_poly = sp.Poly(P, S, domain=sp.QQ)
        degree = p_poly.degree()
        assert sp.gcd(p_poly, p_poly.diff()).degree() == 0
        a = next(
            sp.Integer(k)
            for k in range(-10, 11)
            if sp.diff(P, S).subs(S, k) != 0
            and sp.diff(P, S, 3).subs(S, k) != 0
        )
        G = sp.expand(P.subs(S, S + a) - P.subs(S, a))
        g_poly = sp.Poly(G, S, domain=sp.QQ)
        g1 = g_poly.coeff_monomial(S)
        assert g1 != 0
        assert g_poly.coeff_monomial(S**3) != 0
        assert g_poly.degree() == degree

        mapping = quadratic_gauge_map(G)
        assert sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det()) == -2

        target_c = sp.cancel(-2 * P.subs(S, a) / g1)
        inverse = sp.expand(G - g1 * target_c / 2)
        assert sp.expand(inverse - P.subs(S, S + a)) == 0
        check_scheme_reconstruction(P, a)


def check_minimal_hasse_map() -> None:
    P5 = sp.expand((S**3 - 19) * (S**2 + S + 1))
    G = sp.expand(P5 - P5.subs(S, 0))
    assert P5 == S**5 + S**4 + S**3 - 19 * S**2 - 19 * S - 19
    assert G == S**5 + S**4 + S**3 - 19 * S**2 - 19 * S

    normalized = quadratic_gauge_map(G)
    assert sp.factor(sp.Matrix(normalized).jacobian((x, y, z)).det()) == -2

    t = 1 + x * y
    q = t**2 * z - 19 * y**2 * (1 + 3 * t)
    displayed = (
        t * q,
        19 * y
        - 3 * x * q
        + 38 * t * q
        - 4 * t**2 * x**2 * q**4
        - 5 * t**2 * x**3 * q**5,
        19 * x * (5 - 3 * t)
        + x**3 * z
        + 2 * (x * q) ** 4
        + 3 * (x * q) ** 5,
    )

    # The integral map is exactly diag(1,19,19) after the normalized gauge.
    assert sp.expand(displayed[0] - normalized[0]) == 0
    assert sp.expand(displayed[1] - 19 * normalized[1]) == 0
    assert sp.expand(displayed[2] - 19 * normalized[2]) == 0
    assert sp.factor(sp.Matrix(displayed).jacobian((x, y, z)).det()) == -722
    assert -722 == 19**2 * (-2)

    # The normalized target (1,0,-2) becomes the displayed target (1,0,-38).
    normalized_target = (sp.Integer(1), sp.Integer(0), sp.Integer(-2))
    displayed_target = (
        normalized_target[0],
        19 * normalized_target[1],
        19 * normalized_target[2],
    )
    assert displayed_target == (1, 0, -38)

    g1 = sp.Poly(G, S).coeff_monomial(S)
    normalized_inverse = sp.expand(
        G - g1 * (normalized_target[1] * S**2 + normalized_target[2]) / 2
    )
    assert normalized_inverse == P5
    assert sp.gcd(sp.Poly(P5, S), sp.Poly(sp.diff(P5, S), S)).degree() == 0
    check_scheme_reconstruction(P5, sp.Integer(0))


def check_infinite_family() -> None:
    G = (
        S**5
        - sp.Rational(3, 2) * S**4
        + sp.Rational(3, 2) * S**3
        - sp.Rational(5, 4) * S**2
        + sp.Rational(9, 16) * S
    )
    mapping = quadratic_gauge_map(G)
    assert sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det()) == -2

    B = sp.Rational(32, 9) * ell
    C = (8 * ell + 1) / 3
    inverse = sp.expand(G - sp.Rational(9, 32) * (B * S**2 + C))
    expected = sp.expand(
        ((S - sp.Rational(1, 2)) ** 3 - ell)
        * (S**2 + sp.Rational(3, 4))
    )
    assert sp.expand(inverse - expected) == 0


if __name__ == "__main__":
    check_polynomial_to_fiber_transfer()
    check_minimal_hasse_map()
    check_infinite_family()
    print("PASS: squarefree polynomials in degrees 3, 4, and 5 transfer")
    print("PASS: quotient-ring reconstruction composes in both directions")
    print("PASS: the explicit Hasse map is diag(1,19,19) of the normalized gauge")
    print("PASS: determinant -722, target -38, and inverse P5 agree exactly")
    print("PASS: one fixed quadratic-gauge map carries the infinite family")
