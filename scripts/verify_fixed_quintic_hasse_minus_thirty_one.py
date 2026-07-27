#!/usr/bin/env python3
"""Exact certificate for the Q(sqrt(-31)) Hasse fiber in the fixed quintic map."""

from __future__ import annotations

import sympy as sp


def main() -> None:
    x = sp.symbols("x")
    pi = sp.Integer(5)
    b_target = sp.Rational(-144, 5)
    c_target = sp.Rational(-188, 3125)

    q = x**2 - 8 * x + 47
    h = x**3 + 8 * x**2 + 12 * x + 8
    normalized = sp.expand(q * h)
    expected_normalized = (
        x**5
        - 5 * x**3
        - 2 * b_target * pi * x**2
        + 4 * pi**3 * x
        - 2 * c_target * pi**5
    )
    assert sp.expand(normalized - expected_normalized) == 0

    disc_q = sp.discriminant(q, x)
    disc_h = sp.discriminant(h, x)
    assert disc_q == -31 * 2**2
    assert disc_h == -31 * 8**2
    assert sp.resultant(q, h, x) == 406503

    # q has nonsquare discriminant.  The cubic is irreducible modulo 5.
    assert all(
        (r**3 + 8 * r**2 + 12 * r + 8) % 5
        for r in range(5)
    )

    # These are the only primes dividing the two polynomial discriminants.
    # At 2, -31 == 1 (mod 8), so q splits over Q_2.  At 31, h has an
    # ordinary Hensel root at 15.
    assert (-31) % 8 == 1
    assert int(h.subs(x, 15)) % 31 == 0
    assert int(sp.diff(h, x).subs(x, 15)) % 31 != 0

    s = sp.symbols("s")
    inverse = (
        pi**5 * s**5
        - 5 * pi * s**3
        - 2 * b_target * s**2
        + 4 * s
        - 2 * c_target
    )
    scaled_factorization = sp.expand(
        pi**-5
        * (pi**4 * s**2 - 8 * pi**2 * s + 47)
        * (
            pi**6 * s**3
            + 8 * pi**4 * s**2
            + 12 * pi**2 * s
            + 8
        )
    )
    assert sp.expand(inverse - scaled_factorization) == 0

    target = [sp.Integer(1), pi, b_target, c_target]
    denominator_lcm = sp.ilcm(*[sp.denom(value) for value in target])
    projective = [int(denominator_lcm * value) for value in target]
    common_divisor = sp.igcd(*projective)
    projective = [value // common_divisor for value in projective]
    assert projective == [3125, 15625, -90000, -188]

    print("PASS: normalized and unnormalized target factorizations are exact")
    print("PASS: both factors are irreducible and have common resolvent Q(sqrt(-31))")
    print("PASS: exact local witnesses cover the exceptional primes 2 and 31")
    print("PASS: the complete fiber is everywhere locally soluble with no Q-point")
    print("PASS: primitive projective target =", projective)


if __name__ == "__main__":
    main()
