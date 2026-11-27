#!/usr/bin/env python3
"""Replay the quartic-septic integrals and homogeneous limit tests."""

import sympy as s


def equal(actual, expected):
    assert s.expand(actual - expected) == 0, (actual, expected)


def main():
    z, a, b, e, c, D, F, G = s.symbols("z a b e c D F G")
    da, db, de = s.symbols("da db de")
    q = s.Rational
    variables = (a, b, e)
    rates = (da, db, de)

    def differential(value):
        return s.expand(sum(s.diff(value, x) * rate for x, rate in zip(variables, rates)))

    def radial(form):
        t = s.Symbol("t")
        scaling = {x: t * x for x in variables}
        rate_values = dict(zip(rates, variables))
        path = s.Poly(form.subs(scaling, simultaneous=True).subs(rate_values), t)
        return s.expand(sum(coefficient / (index + 1) for (index,), coefficient in path.terms()))

    U = z**4 + a * z**2 + b * z
    d, f0, g0, dc, dd, df, dg = s.symbols("d f0 g0 dc dd df dg")
    generic_V = z**7 + c * z**6 + d * z**5 + e * z**3 + f0 * z**2 + g0 * z
    generic_dV = dc * z**6 + dd * z**5 + de * z**3 + df * z**2 + dg * z
    generic = s.Poly(U.diff(z) * generic_dV - differential(U) * generic_V.diff(z), z)
    equal(generic.nth(9), 4 * dc)
    equal(generic.nth(8), 4 * dd - 7 * da)

    f = (
        F - q(15, 14) * D * c * a + q(5, 4) * D * b - q(27, 16) * c * a**2
        - q(7, 16) * a * b - q(9, 7) * c**2 * b + q(6, 7) * c * e
    )
    g = (
        G - q(25, 28) * D**2 * a - q(5, 4) * D * a**2 - q(15, 14) * D * c * b
        + q(5, 7) * D * e - q(35, 64) * a**3 - q(15, 8) * c * a * b
        + q(3, 4) * a * e - q(35, 32) * b**2
    )
    V = generic_V.subs({d: q(7, 4) * a + D, f0: f, g0: g}, simultaneous=True)
    form_U2 = -(6 * c * da + 7 * db) / 8
    form_V = (16 * de - 20 * D * da - 21 * a * da - 24 * c * db) / 28
    form_U = (
        6 * a**2 * c * da + 7 * a**2 * db + 4 * a * differential(f) - 4 * f * da
        + 2 * b * de - 6 * e * db
    ) / 8
    residual = s.Poly(
        U.diff(z) * differential(V) - V.diff(z) * differential(U)
        - form_U2 * (U**2).diff(z) - form_V * V.diff(z) - form_U * U.diff(z),
        z,
    )
    assert residual.degree() == 2
    I8 = radial(residual.nth(2))
    I9 = radial(residual.nth(1))
    G10 = radial(residual.nth(0) + I8 * da / 4)
    equal(differential(I8), residual.nth(2))
    equal(differential(I9), residual.nth(1))
    equal(differential(G10), residual.nth(0) + I8 * da / 4)

    initial = {c: 0, D: 0, F: 0, G: 0}
    equal(
        I8.subs(initial),
        -(1225 * a**4 - 2688 * a**2 * e - 1176 * a * b**2 + 1536 * e**2) / 1792,
    )
    equal(I9.subs(initial), -b * (49 * a**3 - 48 * a * e - 14 * b**2) / 64)
    equal(
        G10.subs(initial),
        -(833 * a**5 - 2240 * a**3 * e + 2352 * a**2 * b**2 + 1536 * a * e**2 - 2688 * b**2 * e)
        / 7168,
    )
    U0 = U.subs(a, 1)
    V0 = V.subs(initial).subs(a, 1)
    gamma = G10.subs(initial).subs(a, 1)
    derivative_resultant = s.resultant(U0.diff(z), V0.diff(z), z)
    for vanishing_b, shape, action, tangent, action_norm, tangent_norm in (
        (
            True, 1536 * e**2 - 2688 * e + 1225,
            -(8 * e - 7) / 128, -245 * (48 * e - 49) / 128,
            q(49, 256), q(44118375, 64),
        ),
        (
            False, 1536 * e**2 + 1344 * e - 2891,
            13 * (160 * e - 161) / 512, 15925 * (12576 * e - 12607) / 128,
            q(140777, 512), -q(8238885939375, 32),
        ),
    ):
        def reduce_shape(expression):
            if vanishing_b:
                expression = expression.subs(b, 0)
            else:
                expression = s.rem(expression, b**2 - (49 - 48 * e) / 14, b)
            return s.rem(expression, shape, e)

        equal(reduce_shape(gamma), action)
        equal(reduce_shape(derivative_resultant), tangent)
        equal(s.resultant(action, shape, e), action_norm)
        equal(s.resultant(tangent, shape, e), tangent_norm)
        assert action_norm != 0 and tangent_norm != 0

    print("Verified the quartic-septic parameter integrals and the two homogeneous limit tests.")


if __name__ == "__main__":
    main()
