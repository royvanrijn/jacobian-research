#!/usr/bin/env python3
"""Replay the exact quartic-quintic identities used in the conductor-12 proof."""

import sympy as sp


q = sp.Rational


def equal(left, right):
    assert sp.expand(left - right) == 0, (left, right)


def parameter_equation():
    z, a, b, f, d, e, D, E = sp.symbols("z a b f d e D E")
    da, db, dd, de, df = sp.symbols("da db dd de df")
    U = z**4 + a * z**2 + b * z
    V = z**5 + d * z**3 + e * z**2 + f * z
    rates = ((a, da), (b, db), (d, dd), (e, de), (f, df))

    def derivative(expression):
        return sum(expression.diff(variable) * rate for variable, rate in rates)

    bracket = sp.Poly(U.diff(z) * derivative(V) - derivative(U) * V.diff(z), z)
    assert bracket.degree() == 6
    equal(bracket.nth(6), 4 * dd - 5 * da)
    equal(bracket.nth(5), 4 * de - 5 * db)
    normalized = {d: 5 * a / 4 + D, e: 5 * b / 4 + E, dd: 5 * da / 4, de: 5 * db / 4}
    V = sp.expand(V.subs(normalized, simultaneous=True))
    HU = -(5 * a / 16 + 3 * D / 4) * db - (5 * b / 16 + E / 2) * da
    HV = 4 * df / 5 - (a / 4 + 3 * D / 5) * da
    A = (
        (a + 12 * D / 5) * f + 5 * b**2 / 8 + 2 * E * b
        - 5 * a**3 / 16 - 3 * D * a**2 / 2 - 9 * D**2 * a / 5
    )
    B = (
        (b + 8 * E / 5) * f - (5 * a**2 / 8 + 3 * D * a / 2) * b
        - 3 * E * a**2 / 4 - 6 * D * E * a / 5
    )
    G = (
        -2 * f**2 / 5 + (5 * a / 32 + 3 * D / 8) * b**2 + A * a / 4
        + 5 * a**4 / 256 + D * a**3 / 8 + 9 * D**2 * a**2 / 40
    )
    residual = bracket.as_expr().subs(normalized, simultaneous=True) - HU * U.diff(z) - HV * V.diff(z)
    equal(residual, -derivative(A) * z**2 - derivative(B) * z + derivative(G) - a * derivative(A) / 4)

    scale = sp.Symbol("scale", nonzero=True)
    scaling = {
        variable: variable * scale**weight
        for variable, weight in ((a, 2), (b, 3), (f, 4), (D, 2), (E, 3))
    }
    for invariant, weight in ((A, 6), (B, 7), (G, 8)):
        equal(invariant.subs(scaling, simultaneous=True), scale**weight * invariant)

    A0, B0, G0 = (sp.expand(expression.subs({D: 0, E: 0})) for expression in (A, B, G))
    first = {f: 5 * a**2 / 8}
    curve = 2 * b**2 + a**3
    for expression in (A0, B0):
        equal(sp.rem(expression.subs(first), curve, b), 0)
    equal(sp.rem(G0.subs(first), curve, b), -55 * a**4 / 256)
    cases = (
        ({b: 0, f: 5 * a**2 / 16}, -5 * a**4 / 256),
        ({a: 0, b: 0}, -2 * f**2 / 5),
    )
    for substitution, expected in cases:
        equal(A0.subs(substitution, simultaneous=True), 0)
        equal(B0.subs(substitution, simultaneous=True), 0)
        equal(G0.subs(substitution, simultaneous=True), expected)

    b0, delta = sp.symbols("b0 delta")
    ramified = {D: 0, E: -5 * b0 / 8, b: b0 + delta}
    equal(A.subs(ramified, simultaneous=True) + 5 * b0**2 / 8, a * f + 5 * delta**2 / 8 - 5 * a**3 / 16)
    equal(B.subs(ramified, simultaneous=True), delta * f - a**2 * (5 * delta / 8 + 5 * b0 / 32))
    base_G = -2 * f**2 / 5 + 5 * a * ((b0 + delta)**2 - b0**2) / 32 + 5 * a**4 / 256
    equal(G.subs(ramified, simultaneous=True), base_G + a * (a * f + 5 * delta**2 / 8 - 5 * a**3 / 16) / 4)


def immersion_and_implicit_relations():
    z, beta, X, Y, output = sp.symbols("z beta X Y output")
    shapes = (
        (1, beta, q(5, 8), beta * (80 * X + 41) / 1024),
        (1, 0, q(5, 16), 0),
        (0, 0, 1, 0),
    )
    for a, b, f, linear in shapes:
        U = z**4 + a * z**2 + b * z
        V = z**5 + q(5, 4) * a * z**3 + q(5, 4) * b * z**2 + f * z
        first = U.diff(z)
        second = q(4, 5) * (V.diff(z) - q(5, 4) * z * U.diff(z))
        resultant = sp.resultant(first, second, z)
        if b == beta:
            equal(resultant, -beta**2 * (2 * beta - 3) * (2 * beta + 3))
            equal(sp.rem(resultant, 2 * beta**2 + 1, beta), -11 / sp.Integer(2))
        else:
            assert resultant != 0
        implicit = sp.Poly(sp.resultant(U - X, V - Y, z), Y)
        assert implicit.degree() == 4 and implicit.LC() == 1
        depressed = sp.Poly(implicit.as_expr().subs(Y, output - implicit.nth(3) / 4), output)
        equal(depressed.nth(3), 0)
        coefficient = depressed.nth(1)
        if b == beta:
            coefficient = sp.rem(coefficient, 2 * beta**2 + 1, beta)
        equal(coefficient, linear)


def ramified_gap_jet():
    z, r, h, v, f = sp.symbols("z r h v f")
    U = z**4 - 4 * r**3 * z
    V = z**5 - q(5, 2) * r**3 * z**2 + f * z
    endpoint = r + h
    equations = [
        sp.cancel((polynomial.subs(z, endpoint + v) - polynomial.subs(z, endpoint - v)) / (2 * v))
        for polynomial in (U, V)
    ]
    equal(equations[0].diff(h).subs({h: 0, v: 0}), 12 * r**2)
    jets = [sp.series(eq.subs(h, -v**2 / (3 * r)), v, 0, 4).removeO() for eq in equations]
    equal(jets[0], 0)
    equal(jets[1], f + 5 * r**2 * v**2)


def main():
    parameter_equation()
    immersion_and_implicit_relations()
    ramified_gap_jet()
    print("Verified the quartic-quintic parameter, immersion, and specialization identities.")


if __name__ == "__main__":
    main()
