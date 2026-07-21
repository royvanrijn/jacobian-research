#!/usr/bin/env python3
"""Exact regressions for GENERALIZED_CANCELLATION_MECHANISM.md.

Finite ranges in this script are regression evidence, not all-degree proofs.
"""
from __future__ import annotations

import sympy as sp

from master_cancellation import integrate_unit_interval, phi


def check_two_weight_jacobian() -> None:
    s, y, B, C = sp.symbols("s y B C", nonzero=True)
    # Check the factorized coordinate identities on a nonmonomial f.  This
    # avoids making a huge expanded determinant the test oracle.
    f = y**2 + y + 1
    for a, b in ((1, 2), (1, 3), (2, 2)):
        exponent = a + b - 2
        A = 1 / (1 - s * f)
        P = A**a * B
        Q = y + s * P
        middle = sp.det(sp.Matrix([P, Q]).jacobian([y, B]))
        assert sp.cancel(middle + A**a) == 0

        D = 1 - s * f.subs(y, Q - s * P)
        assert sp.cancel(D - 1 / A) == 0

        source_factor = A**b * A**-2 * middle
        resolvent_factor = C * D**exponent
        assert sp.cancel(source_factor * resolvent_factor + C) == 0

        # The forced monomial gives a constant; adding a neighboring power
        # immediately reintroduces A-dependence.
        assert sp.cancel(A**exponent * D**exponent - 1) == 0
        perturbed = sp.cancel(A**exponent * (D**exponent + D ** (exponent + 1)))
        assert sp.cancel(sp.diff(perturbed, s)) != 0


def generalized_phi(
    f: sp.Expr,
    exponent: int,
    a: int,
    A: sp.Symbol,
    y: sp.Symbol,
    G: sp.Expr,
) -> sp.Expr:
    u = sp.Dummy("u")
    argument = y + (A - 1) * A ** (a - 1) * G * (1 - u) / f
    bracket = A - (A - 1) * u * sp.cancel(f.subs(y, argument) / f)
    return sp.cancel(integrate_unit_interval(sp.expand(bracket**exponent), u))


def check_resolvent_degree_and_critical_power() -> None:
    T, P, Q, R, C = sp.symbols("T P Q R C")
    t = sp.Dummy("t")
    f_at_line = (Q - P * t) ** 2 + (Q - P * t) + 1
    D = 1 - T * ((Q - P * T) ** 2 + (Q - P * T) + 1)
    for exponent in (1, 2, 3):
        integrand = (1 - t * f_at_line) ** exponent
        primitive = sp.integrate(sp.expand(integrand), t)
        psi = sp.expand(C * (primitive.subs(t, T) - primitive.subs(t, 0)) - R)
        assert sp.Poly(psi, T).degree() == exponent * 3 + 1
        assert sp.cancel(sp.diff(psi, T) - C * D**exponent) == 0


def check_weight_obstruction_and_c24_reduction() -> None:
    A, y, H = sp.symbols("A y H")
    f = y**3
    G = y**4 * H
    for exponent in (1, 2, 3):
        obstructed = generalized_phi(f, exponent, 2, A, y, G)
        assert sp.cancel(obstructed.subs(A, 0) - sp.Rational(1, exponent + 1)) == 0

        recovered = generalized_phi(f, exponent, 1, A, y, G)
        expected = phi(3, exponent, A, H)
        assert sp.cancel(recovered - expected) == 0


def spectral_polynomial(n: int, q: sp.Symbol) -> sp.Expr:
    u = sp.Dummy("u")
    return sp.factor(integrate_unit_interval(u * (1 - q * (1 - u)) ** n, u))


def check_e1_spectrum() -> None:
    y, alpha, q = sp.symbols("y alpha q")
    for n in range(1, 9):
        f = (y - alpha) ** n
        v = q * (y - alpha)
        lhs = sum(
            (-v) ** j
            * sp.diff(f, y, j)
            / (sp.factorial(j) * (j + 1) * (j + 2))
            for j in range(n + 1)
        )
        assert sp.factor(lhs - f * spectral_polynomial(n, q)) == 0

    polynomials = [sp.Poly(spectral_polynomial(n, q), q) for n in range(1, 13)]
    for index, left in enumerate(polynomials):
        for right in polynomials[:index]:
            assert sp.gcd(left, right).degree() == 0


def main() -> None:
    check_two_weight_jacobian()
    print("PASS: representative two-weight Jacobians are constant")
    check_resolvent_degree_and_critical_power()
    print("PASS: generalized resolvent degree and critical power")
    check_weight_obstruction_and_c24_reduction()
    print("PASS: weight obstruction and exact C24 operator reduction")
    check_e1_spectrum()
    print("PASS: e=1 spectral identity and pairwise gcds through degree 12")


if __name__ == "__main__":
    main()
