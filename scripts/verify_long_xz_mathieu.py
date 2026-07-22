#!/usr/bin/env python3
"""Exact, dependency-free checks for Long's (xz) and SU(2) witnesses.

Source: Christopher D. Long, "Counterexamples to the (xz)-Conjecture and
the Mathieu Conjecture for SU(2)", arXiv:2607.19012v1 (21 July 2026).

The (xz) calculation is self-contained.  This script verifies the SU(2)
algebraic substitution and monomial identities on the right side of the
Mueger--Tuset integration formula quoted by Long.  The companion
verify_long_su2_haar.py supplies the independent Haar-measure proof.
"""

from fractions import Fraction
from math import comb, factorial


Laurent = dict[tuple[int, int, int], Fraction]  # powers of x, z1, z2


def add(*polys: Laurent) -> Laurent:
    out: Laurent = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] = out.get(monomial, Fraction(0)) + coefficient
            if not out[monomial]:
                del out[monomial]
    return out


def multiply(left: Laurent, right: Laurent) -> Laurent:
    out: Laurent = {}
    for (ax, az1, az2), ac in left.items():
        for (bx, bz1, bz2), bc in right.items():
            monomial = (ax + bx, az1 + bz1, az2 + bz2)
            out[monomial] = out.get(monomial, Fraction(0)) + ac * bc
    return {monomial: coefficient for monomial, coefficient in out.items() if coefficient}


def power(poly: Laurent, exponent: int) -> Laurent:
    out: Laurent = {(0, 0, 0): Fraction(1)}
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = multiply(out, base)
        base = multiply(base, base)
        n //= 2
    return out


def integral_constant_term(poly: Laurent) -> Fraction:
    """Apply integral_0^1 CT_(z1,z2), term by term."""
    total = Fraction(0)
    for (x_degree, z1_degree, z2_degree), coefficient in poly.items():
        if z1_degree == z2_degree == 0:
            total += coefficient * Fraction(1, x_degree + 1)
    return total


ONE: Laurent = {(0, 0, 0): Fraction(1)}
X: Laurent = {(1, 0, 0): Fraction(1)}
Z1: Laurent = {(0, 1, 0): Fraction(1)}
Z1_INV: Laurent = {(0, -1, 0): Fraction(1)}
Z2: Laurent = {(0, 0, 1): Fraction(1)}
Z2_INV: Laurent = {(0, 0, -1): Fraction(1)}

# f=(1-z^-1)((1-x)+xz)
f = multiply(
    add(ONE, {(0, -1, 0): Fraction(-1)}),
    add(ONE, {(1, 0, 0): Fraction(-1)}, multiply(X, Z1)),
)


def beta_integral(n: int, k: int) -> Fraction:
    """Integral of x^k(1-x)^(n-k), expanded exactly."""
    return sum(
        Fraction((-1) ** j * comb(n - k, j), k + j + 1)
        for j in range(n - k + 1)
    )


def mueger_tuset_monomial_rhs(r: int, s: int, t: int, u: int) -> Fraction:
    """Algebraic beta-map torus/Beta integral for a^r b^s c^t d^u."""
    if r != u or s != t:
        return Fraction(0)
    return Fraction((-1) ** t * factorial(r) * factorial(s), factorial(r + s + 1))


def main() -> None:
    # Proof-oriented beta/binomial identity.
    for n in range(1, 21):
        for k in range(n + 1):
            assert comb(n, k) * beta_integral(n, k) == Fraction(1, n + 1)

    # Direct Laurent-polynomial regression of the two moments.
    f_power = ONE
    for n in range(1, 16):
        f_power = multiply(f_power, f)
        assert integral_constant_term(f_power) == 0
        assert integral_constant_term(multiply(Z1_INV, f_power)) == Fraction(
            (-1) ** (n - 1), n + 1
        )

    # Mueger--Tuset beta substitution, in Long's coordinate order (a,b,c,d):
    # ((1-x)z2, xz1, -z1^-1, z2^-1).
    a = add(Z2, multiply({(0, 0, 1): Fraction(-1)}, X))
    b = multiply(X, Z1)
    c = {(0, -1, 0): Fraction(-1)}
    d = Z2_INV
    F_beta = multiply(add(ONE, c), add(multiply(a, d), b))
    G_beta = {(0, -1, 0): Fraction(1)}  # -c
    assert F_beta == f
    assert G_beta == Z1_INV

    # Check the quoted formula's monomial right side in a useful exact box.
    for r in range(5):
        for s in range(5):
            for t in range(5):
                for u in range(5):
                    beta_image = multiply(
                        multiply(power(a, r), power(b, s)),
                        multiply(power(c, t), power(d, u)),
                    )
                    assert integral_constant_term(beta_image) == mueger_tuset_monomial_rhs(
                        r, s, t, u
                    )

    print("PASS Long xz: beta/binomial identity n=1..20")
    print("PASS Long xz: exact Laurent moments n=1..15")
    print("PASS Long SU(2): beta substitution and monomial RHS in degrees 0..4")
    print("PASS Long SU(2): combine with verify_long_su2_haar.py for the full Haar proof")


if __name__ == "__main__":
    main()
