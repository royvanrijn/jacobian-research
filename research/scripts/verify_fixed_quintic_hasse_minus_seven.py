#!/usr/bin/env python3
"""Exact certificate for the Q(sqrt(-7)) Hasse fiber in the fixed quintic map."""

from __future__ import annotations

import sympy as sp


def main() -> None:
    x = sp.symbols("x")
    pi = sp.Integer(-7)
    b_target = sp.Rational(387, 14)
    c_target = sp.Rational(400, 2401)

    q = x**2 - 4 * x + 32
    h = x**3 + 4 * x**2 - 21 * x + 175
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
    assert disc_q == -7 * 4**2
    assert disc_h == -7 * (5 * 79) ** 2
    assert sp.resultant(q, h, x) == 27477

    # q is irreducible because its discriminant is nonsquare.  The primitive
    # cubic is irreducible because x^3+x+1 has no root over F_2.
    assert all((r**3 + r + 1) % 2 for r in range(2))

    # Exact local witnesses at every prime dividing either factor's
    # polynomial discriminant.  At 2, -7 == 1 (mod 8), hence -7 is a square
    # in Q_2 and q splits.  The other witnesses are ordinary Hensel lifts.
    assert (-7) % 8 == 1
    witnesses = {
        5: (h, 0),
        7: (h, 3),
        79: (q, 31),
    }
    for prime, (poly, residue) in witnesses.items():
        assert int(poly.subs(x, residue)) % prime == 0
        assert int(sp.diff(poly, x).subs(x, residue)) % prime != 0

    # Check the unnormalized inverse equation as well:
    # E(S) = pi^-5 * normalized(pi^2 S).
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
        * (pi**4 * s**2 - 4 * pi**2 * s + 32)
        * (
            pi**6 * s**3
            + 4 * pi**4 * s**2
            - 21 * pi**2 * s
            + 175
        )
    )
    assert sp.expand(inverse - scaled_factorization) == 0

    primitive_target = [
        1,
        pi,
        b_target,
        c_target,
    ]
    denominator_lcm = sp.ilcm(
        *[sp.denom(value) for value in primitive_target]
    )
    projective = [
        int(denominator_lcm * value) for value in primitive_target
    ]
    common_divisor = sp.igcd(*projective)
    projective = [value // common_divisor for value in projective]
    assert projective == [4802, -33614, 132741, 800]

    print("PASS: normalized and unnormalized target factorizations are exact")
    print("PASS: both factors are irreducible and have common resolvent Q(sqrt(-7))")
    print("PASS: exact local witnesses cover 2, 5, 7, and 79")
    print("PASS: the complete fiber is everywhere locally soluble with no Q-point")
    print("PASS: primitive projective target =", projective)


if __name__ == "__main__":
    main()
