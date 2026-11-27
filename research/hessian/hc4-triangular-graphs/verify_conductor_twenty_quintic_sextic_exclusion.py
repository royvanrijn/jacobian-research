#!/usr/bin/env python3
"""Exact parameter integrals and initial-shape tests for conductor twenty (6,5)."""

import sympy as s


def equal(actual, expected):
    assert s.expand(actual - expected) == 0, (actual, expected)


def main():
    z, a, b, c, g, D, E, F = s.symbols("z a b c g D E F")
    da, db, dc, dg = s.symbols("da db dc dg")
    variables = (a, b, c, g)
    rates = (da, db, dc, dg)
    q = s.Rational

    def differential(value):
        return s.expand(sum(s.diff(value, x) * rate for x, rate in zip(variables, rates)))

    def radial(form):
        t = s.Symbol("t")
        scaling = {x: t * x for x in variables}
        rate_values = dict(zip(rates, variables))
        path = s.Poly(form.subs(scaling, simultaneous=True).subs(rate_values), t)
        return s.expand(sum(coefficient / (index + 1) for (index,), coefficient in path.terms()))

    U = z**5 + a * z**3 + b * z**2 + c * z
    d, e, f0, dd, de, df = s.symbols("d e f0 dd de df")
    general_V = z**6 + d * z**4 + e * z**3 + f0 * z**2 + g * z
    general_dV = dd * z**4 + de * z**3 + df * z**2 + dg * z
    general = s.Poly(U.diff(z) * general_dV - differential(U) * general_V.diff(z), z)
    equal(general.nth(8), 5 * dd - 6 * da)
    equal(general.nth(7), 5 * de - 6 * db)

    d_value, e_value = 6 * a / 5 + D, 6 * b / 5 + E
    f_value = 4 * D * a / 5 + 3 * a**2 / 25 + 6 * c / 5 + F
    top = {d: d_value, e: e_value, dd: 6 * da / 5, de: 6 * db / 5}
    equal(general.nth(6).subs(top, simultaneous=True), 5 * (df - differential(f_value)))
    V = general_V.subs({d: d_value, e: e_value, f0: f_value}, simultaneous=True)
    theta_V = (25 * dg - 20 * D * db - 15 * E * da - 6 * a * db - 6 * b * da) / 30
    theta_U = (
        20 * D * a * da - 100 * D * dc - 75 * E * db - 50 * F * da
        + 12 * a**2 * da - 30 * a * dc - 30 * c * da - 30 * b * db
    ) / 125
    residual = s.Poly(
        U.diff(z) * differential(V) - V.diff(z) * differential(U)
        - theta_U * U.diff(z) - theta_V * V.diff(z),
        z,
    )
    assert residual.degree() == 3
    I7 = radial(residual.nth(3))
    I8 = radial(residual.nth(2))
    I9 = radial(residual.nth(1) + 2 * I7 * da / 5)
    G10 = radial(residual.nth(0) + (I7 * db + I8 * da) / 5)
    equal(differential(I7), residual.nth(3))
    equal(differential(I8), residual.nth(2))
    equal(differential(I9), residual.nth(1) + 2 * I7 * da / 5)
    equal(differential(G10), residual.nth(0) + (I7 * db + I8 * da) / 5)

    t = s.Symbol("t")
    weights = {a: 2, b: 3, c: 4, g: 5, D: 2, E: 3, F: 4}
    scaling = {variable: variable * t**weight for variable, weight in weights.items()}
    for expression, weight in ((I7, 7), (I8, 8), (I9, 9), (G10, 10)):
        equal(expression.subs(scaling, simultaneous=True), t**weight * expression)

    initial = {D: 0, E: 0, F: 0}
    equal(I7.subs(initial), (18 * a**2 * b - 25 * a * g - 30 * b * c) / 25)
    equal(I8.subs(initial), -(9 * a**4 - 60 * a**2 * c - 90 * a * b**2
                                + 125 * b * g + 75 * c**2) / 125)
    equal(I9.subs(initial), (6 * a**3 * b - 25 * a**2 * g + 60 * a * b * c
                               + 20 * b**3 - 125 * c * g) / 125)
    equal(G10.subs(initial), -(108 * a**5 - 2700 * a**2 * b**2
                                  - 4500 * a * c**2 + 15625 * g**2) / 37500)

    vanishing_a = {a: 0, b: 0, c: 0}
    U0, V0 = (p.subs(initial).subs(vanishing_a) for p in (U, V))
    equal(G10.subs(initial).subs(vanishing_a), -q(5, 12) * g**2)
    equal(s.resultant(U0.diff(z), V0.diff(z), z), 3125 * g**4)

    normalized = {a: 1, g: b * (18 - 30 * c) / 25}
    U0, V0 = (p.subs(initial).subs(normalized) for p in (U, V))
    initial_G = G10.subs(initial).subs(normalized)
    immersion = s.resultant(U0.diff(z), V0.diff(z), z)
    equal(I8.subs(initial).subs(normalized), 3 * (50 * b**2 * c - 25 * c**2 + 20 * c - 3) / 125)
    equal(I9.subs(initial).subs(normalized), 2 * b * (10 * b**2 + 75 * c**2 - 6) / 125)

    square_relation = b**2 - (q(3, 5) - q(15, 2) * c**2)
    for vanishing_b, equation, scale, derivative in (
        (True, 75 * c**2 - 60 * c + 9, 6 * (50 * c - 9) / 3125,
         1296 * (365 * c - 72) / 3125),
        (False, 375 * c**3 + 25 * c**2 - 50 * c + 3,
         91 * (125 * c**2 - 25 * c - 3) / 6250,
         -2548 * (24200 * c**2 - 7555 * c + 132) / 3125),
    ):
        def reduce_shape(expression):
            if vanishing_b:
                expression = expression.subs(b, 0)
            else:
                expression = s.rem(expression, square_relation, b)
            return s.rem(expression, equation, c)

        equal(reduce_shape(initial_G), scale)
        equal(reduce_shape(immersion), derivative)
        assert s.gcd(scale, equation) == s.gcd(derivative, equation) == 1

    print("Verified the quintic-sextic parameter integrals and all homogeneous limit tests.")


if __name__ == "__main__":
    main()
