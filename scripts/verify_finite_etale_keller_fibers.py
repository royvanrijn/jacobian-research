#!/usr/bin/env python3
"""Exact checks for finite etale algebras as Keller fibers."""
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


def check_minimal_hasse_map() -> None:
    t = 1 + x * y
    q = t**2 * z - 19 * y**2 * (1 + 3 * t)
    mapping = (
        t * q,
        19 * y - 3 * x * q + 38 * t * q - 4 * t**2 * x**2 * q**4 - 5 * t**2 * x**3 * q**5,
        19 * x * (5 - 3 * t) + x**3 * z + 2 * (x * q) ** 4 + 3 * (x * q) ** 5,
    )
    assert sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det()) == -722

    P5 = sp.expand((S**3 - 19) * (S**2 + S + 1))
    G = sp.expand(P5 - P5.subs(S, 0))
    assert P5 == S**5 + S**4 + S**3 - 19 * S**2 - 19 * S - 19
    assert G == S**5 + S**4 + S**3 - 19 * S**2 - 19 * S
    # The displayed map scales the normalized B,C coordinates by 19.
    # Thus C=-38 corresponds to normalized C=-2.
    assert sp.expand(G + sp.Rational(1, 2) * (-38) - P5) == 0
    assert sp.gcd(sp.Poly(P5, S), sp.Poly(sp.diff(P5, S), S)).degree() == 0


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
    print("PASS: squarefree polynomials of degrees 3, 4, and 5 transfer to determinant-minus-two maps")
    print("PASS: the explicit minimal Hasse map has determinant -722 and inverse P5")
    print("PASS: one fixed quadratic-gauge map carries the infinite intersective family")
