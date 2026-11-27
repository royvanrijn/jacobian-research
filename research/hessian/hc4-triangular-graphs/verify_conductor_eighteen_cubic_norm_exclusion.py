#!/usr/bin/env python3
"""Exact cubic-norm replay for conductor eighteen and sequence (10,6,3)."""

import sympy as sp


def equal(actual, expected):
    assert sp.expand(actual - expected) == 0, (actual, expected)


def main():
    z, U, Y, a, h, p, q = sp.symbols("z U Y a h p q")
    cubic = z**3 + a * z + h - U
    candidate = q * z + p * (z**2 + 2 * a / 3)
    norm = sp.Poly(sp.resultant(cubic, Y - candidate, z), Y)
    assert norm.degree() == 3 and norm.LC() == 1 and norm.nth(2) == 0
    equal(norm.nth(1), a * q**2 - 3 * (U - h) * p * q - a**2 * p**2 / 3)
    equal(
        norm.nth(0),
        -(U - h)**2 * p**3 + (U - h) * a * p**2 * q - (U - h) * q**3
        - 2 * a**3 * p**3 / 27 - 2 * a**2 * p * q**2 / 3,
    )

    p0, p1, p2, q0, q1, q2 = sp.symbols("p0 p1 p2 q0 q1 q2")
    terms = {p: p2 * U**2 + p1 * U + p0, q: U**3 + q2 * U**2 + q1 * U + q0}
    A = sp.Poly(norm.nth(1).subs(terms, simultaneous=True), U)
    B = sp.Poly(norm.nth(0).subs(terms, simultaneous=True), U)

    equal(B.nth(9), h - 3 * q2)
    top = {q2: h / 3}
    p1_value = 2 * a * h / 9 + 2 * h * p2 / 3
    equal(A.nth(5).subs(top, simultaneous=True), -3 * (p1 - p1_value))
    top[p1] = p1_value
    q0_value = h * (-4 * a**3 - 12 * a**2 * p2 + 8 * h**2 + 27 * q1) / 81
    equal(B.nth(7).subs(top, simultaneous=True), -3 * (q0 - q0_value))
    top[q0] = q0_value
    equal(B.nth(1).subs(top, simultaneous=True).subs(h, 0), 0)

    p0_value = (4 * a**4 + 12 * a**3 * p2 - 17 * a * h**2 - 27 * a * q1 - 15 * h**2 * p2) / 81
    equal(A.nth(3).subs(top, simultaneous=True), 2 * h * (p0 - p0_value))
    top[p0] = p0_value

    def normalized(coefficient):
        return sp.expand(coefficient.subs(top, simultaneous=True).subs(h, 1))

    q1_value = 13 * a**3 / 27 + a**2 * p2 - sp.Rational(5, 9)
    equal(normalized(B.nth(5)), 8 * (q1 - q1_value) / 27)
    equal(normalized(A.nth(1)).subs(q1, q1_value), 32 * (5 * a + 6 * p2) / 2187)
    reduced_B3 = normalized(B.nth(3)).subs(q1, q1_value)
    equal(
        reduced_B3,
        32 * (-15 * a**3 - 28 * a**2 * p2 + 3 * a * p2**2 + 18 * p2**3 + 4) / 6561,
    )
    equal(reduced_B3.subs(p2, -5 * a / 6), sp.Rational(128, 6561))

    print("Verified the cubic norm and coefficient contradiction for sequence (10,6,3).")


if __name__ == "__main__":
    main()
