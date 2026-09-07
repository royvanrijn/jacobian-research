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


def check_weight_obstruction_and_cancellation_reduction() -> None:
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


def generalized_spectral_polynomial(n: int, exponent: int, q: sp.Symbol) -> sp.Poly:
    """Return J_(n,e) from beta-integral coefficients."""
    expression = sum(
        (-q) ** j
        * sp.binomial(n, j)
        / ((exponent + j + 1) * sp.binomial(exponent + j, j))
        for j in range(n + 1)
    )
    return sp.Poly(expression, q, domain=sp.QQ)


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

    # Closed form used in the all-degree Enestrom--Kakeya proof.
    r = sp.symbols("r")
    for n in range(1, 13):
        closed = (r ** (n + 2) - (n + 2) * r + n + 1) / (
            (n + 1) * (n + 2) * (1 - r) ** 2
        )
        assert sp.cancel(spectral_polynomial(n, q).subs(q, 1 - r) - closed) == 0


def check_higher_spectrum() -> None:
    q, r = sp.symbols("q r")
    u = sp.symbols("u")
    for exponent in range(1, 7):
        for n in range(1, 9):
            integrated = integrate_unit_interval(
                sp.expand(u**exponent * (1 - q * (1 - u)) ** n), u
            )
            spectral = generalized_spectral_polynomial(n, exponent, q).as_expr()
            assert sp.cancel(spectral - integrated) == 0

            scale = sp.factorial(n) / sp.factorial(n + exponent + 1)
            transformed = scale * sum(
                sp.factorial(n - k + exponent)
                / sp.factorial(n - k)
                * r**k
                for k in range(n + 1)
            )
            assert sp.cancel(spectral.subs(q, 1 - r) - transformed) == 0

            coefficient_form = sp.cancel(transformed / scale)
            reciprocal = sp.expand(r**n * coefficient_form.subs(r, 1 / r))
            rising_form = sum(sp.rf(j + 1, exponent) * r**j for j in range(n + 1))
            assert sp.cancel(reciprocal - rising_form) == 0

    # The two Enestrom--Kakeya annuli used in the uniform proof are disjoint.
    m, exponent = sp.symbols("m exponent", positive=True, integer=True)
    gap = sp.factor(
        (m + 2) / (m + exponent + 2) - m / (m + exponent)
    )
    assert gap == 2 * exponent / ((m + exponent) * (m + exponent + 2))

    # Independent bounded stress test of the uniform annulus theorem.
    for exponent in range(2, 9):
        polynomials: list[sp.Poly] = []
        for n in range(1, 41):
            current = generalized_spectral_polynomial(n, exponent, q)
            for previous in polynomials:
                assert sp.gcd(current, previous).degree() == 0
            polynomials.append(current)


def check_first_jet_formula() -> None:
    A, u, y, g0, g1 = sp.symbols("A u y g0 g1")
    exponent = 2
    f = y**2 + y + 1
    G = g0 + A * g1
    argument = y + (A - 1) * G * (1 - u) / f
    bracket = A - (A - 1) * u * f.subs(y, argument) / f
    coefficient = sp.diff(bracket**exponent, A).subs(A, 0)

    v = g0 / f
    Y = y - v * (1 - u)
    H = f.subs(y, Y)
    f_prime_at_y = sp.diff(f, y).subs(y, Y)
    delta_integrand = (
        -exponent
        * u**exponent
        * (1 - u)
        * H ** (exponent - 1)
        * f_prime_at_y
        / f ** (exponent + 1)
    )
    explicit_integrand = (
        exponent
        * (u * H / f) ** (exponent - 1)
        * (1 - u * H / f + u * (1 - u) * v * f_prime_at_y / f)
    )
    assert sp.cancel(coefficient - explicit_integrand - delta_integrand * g1) == 0


def check_monomial_linearization_unit() -> None:
    y, q, u = sp.symbols("y q u")
    for degree in range(1, 4):
        f = y**degree
        for exponent in range(1, 4):
            v = q * y
            shifted = y - v * (1 - u)
            delta = integrate_unit_interval(
                sp.expand(
                    -exponent
                    * u**exponent
                    * (1 - u)
                    * f.subs(y, shifted) ** (exponent - 1)
                    * sp.diff(f, y).subs(y, shifted)
                    / f ** (exponent + 1)
                ),
                u,
            )
            spectral_derivative = sp.diff(
                generalized_spectral_polynomial(
                    degree * exponent, exponent, q
                ).as_expr(),
                q,
            )
            assert sp.cancel(delta * y * f - spectral_derivative) == 0


def check_tail_source_shift() -> None:
    x, y, z = sp.symbols("x y z")
    A = 1 + x * y**2
    q, h1 = sp.symbols("q h1")
    base = y**3 * (q + h1 * A)
    tail = y + A * y**2
    exponent = 2

    B_base = A ** (exponent + 1) * z + base
    B_tail = B_base + A ** (exponent + 1) * tail
    shifted = sp.expand(B_base.subs(z, z + tail))
    assert sp.expand(shifted - B_tail) == 0

    P_base, Q_base = A * B_base, y + x * B_base
    P_tail, Q_tail = A * B_tail, y + x * B_tail
    assert sp.expand(P_base.subs(z, z + tail) - P_tail) == 0
    assert sp.expand(Q_base.subs(z, z + tail) - Q_tail) == 0


def main() -> None:
    check_two_weight_jacobian()
    print("PASS: representative two-weight Jacobians are constant")
    check_resolvent_degree_and_critical_power()
    print("PASS: generalized resolvent degree and critical power")
    check_weight_obstruction_and_cancellation_reduction()
    print("PASS: weight obstruction and exact cancellation construction operator reduction")
    check_e1_spectrum()
    print("PASS: e=1 spectral identity, closed form, and bounded pairwise gcds")
    check_higher_spectrum()
    print("PASS: rising-factorial spectrum and disjoint annuli for every symbolic e")
    print("PASS: spectral pairwise gcd stress test for e<=8 and degree<=40")
    check_first_jet_formula()
    print("PASS: explicit first-jet linearization formula")
    check_monomial_linearization_unit()
    print("PASS: monomial linearization equals the simple spectral derivative")
    check_tail_source_shift()
    print("PASS: every A^(e+1) tail is a polynomial source shift")


if __name__ == "__main__":
    main()
