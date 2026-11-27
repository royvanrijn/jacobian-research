#!/usr/bin/env python3
"""Exact coefficient and critical-point tests for conductor twenty (12,8,5)."""

import sympy as sp


def equal(actual, expected):
    assert sp.cancel(actual - expected) == 0, (actual, expected)


def main():
    z, a, b, c, ell = sp.symbols("z a b c ell")
    delta, X, rate = sp.symbols("Delta X rate", nonzero=True)
    m = z**3 + a * z**2 + b * z + c
    L = (
        z**4 + 3 * a * z**3 / 2 + (3 * b / 2 + 3 * a**2 / 8) * z**2
        + (3 * c / 2 + 3 * a * b / 4 - a**3 / 16) * z + ell
    )
    defect = sp.Poly(m**3 - z * L**2, z)
    assert defect.degree() <= 5
    remainder = sp.Poly(sp.rem(m**2, L, z), z)
    c_value = (4 * a * b - a**3) / 8
    equal(remainder.nth(3), (c - c_value) / 2)
    ell_value = (5 * a**4 - 24 * a**2 * b + 80 * a * c + 16 * b**2) / 64
    equal(remainder.nth(2), ell_value - ell)

    def reduced(expression):
        return sp.expand(expression.subs(ell, ell_value).subs(c, c_value)
                         .subs(b, (a**2 - delta) / 4))

    equal(reduced(remainder.nth(1)), a * delta**2 / 128)
    equal(reduced(remainder.nth(0)), delta**2 * (2 * a**2 + delta) / 512)
    equal(reduced(defect.nth(5)), delta**2 / 64)

    # The level k is proportional to a**16 * X**6 * (2 + X).
    level = a**16 * X**6 * (2 + X)
    primitive = a**19 * X**8
    a_rate = -a * rate * (6 / X + 1 / (2 + X)) / 16
    equal(level.diff(a) * a_rate + level.diff(X) * rate, 0)
    equal(
        (primitive.diff(a) * a_rate + primitive.diff(X) * rate) / primitive,
        (28 - 5 * X) * rate / (16 * X * (2 + X)),
    )

    L0 = reduced(L).subs({a: 1, delta: -2})
    M0 = sp.Rational(3, 16) * reduced(m).subs({a: 1, delta: -2})
    quotient, rest = sp.div(M0**2, L0, z)
    equal(rest, sp.Rational(9, 8192) * z)
    W0 = L0**2 + M0
    V0 = L0**3 + 3 * L0 * M0 / 2 + 3 * quotient / 8
    T0 = sp.Poly(V0**2 - W0**3, z)
    assert T0.degree() == 5
    equal(T0.LC(), sp.Rational(27, 524288))
    assert sp.Poly(sp.gcd(W0.diff(z), T0.diff(z).as_expr()), z).monic().as_expr() == 1

    print("Verified the (12,8,5) coefficient reduction and the fixed critical-point test.")


if __name__ == "__main__":
    main()
