#!/usr/bin/env python3
"""Quick local-geometry screen for the proposed A4 two-mask divisors."""

import sympy as sp


a, b, z = sp.symbols("a b z")

B = (
    a**3
    - 3 * a * b**2
    + 2 * b**3
    - 9 * a * b
    + 9 * b**2
    - 27 * a
    + 27 * b
    + 27
)
rho = b**2 + 3 * b + 9
sigma = (
    2 * a**3 * b
    + 3 * a**3
    - 3 * a**2 * b**2
    - 9 * a**2 * b
    - 27 * a**2
    + b**4
    + 6 * b**3
    + 27 * b**2
    + 54 * b
    + 81
)
p3 = (
    8 * a**2 * b
    + 12 * a**2
    - 16 * a * b**2
    - 48 * a * b
    - 81 * a
    + 8 * b**3
    + 36 * b**2
    + 108 * b
    + 108
)


def reduced_groebner(expressions):
    basis = sp.groebner(
        expressions,
        a,
        b,
        order="lex",
        domain=sp.QQ,
    )
    return tuple(sp.factor(poly.as_expr()) for poly in basis.polys)


def tangent_determinant(first, second):
    return (
        sp.diff(first, a) * sp.diff(second, b)
        - sp.diff(first, b) * sp.diff(second, a)
    )


# B and rho are smooth, but sigma has two conjugate singular points.
assert reduced_groebner([B, sp.diff(B, a), sp.diff(B, b)]) == (1,)
assert reduced_groebner(
    [rho, sp.diff(rho, a), sp.diff(rho, b)]
) == (1,)
assert reduced_groebner(
    [sigma, sp.diff(sigma, a), sp.diff(sigma, b)]
) == (a**2, rho)

# Every pair fails transversality at the same two points, and all three
# components meet there.
for first, second in ((B, rho), (B, sigma), (rho, sigma)):
    assert reduced_groebner([
        first,
        second,
        tangent_determinant(first, second),
    ]) == (a**2, rho)
assert reduced_groebner([B, rho, sigma]) == (a**3, rho)

# Use z=rho and c=2b+3.  Since c^2+27=4z, c is a unit at z=0.
c = 2 * b + 3
assert sp.expand(c**2 + 27 - 4 * rho) == 0


def rewrite_with_z(expression):
    return sp.factor(
        sp.rem(expression, b**2 + 3 * b + 9 - z, b)
    )


assert sp.expand(
    rewrite_with_z(B) - (a**3 - 3 * a * z + c * z)
) == 0
assert sp.expand(
    rewrite_with_z(sigma) - (c * a**3 - 3 * a**2 * z + z**2)
) == 0
assert sp.expand(
    rewrite_with_z(p3)
    - (4 * c * a**2 - 16 * a * z + 63 * a + 4 * c * z)
) == 0

# At the common cluster (a,z)=(0,0), the reduced orders are
#
# ord(B)=1, ord(rho)=1, ord(sigma)=2, ord(p3)=1.
#
# Hence B^2 has order two and rho*sigma has order three, while the T^3
# numerator coefficient supplies only one order of cancellation.
print("PASS: B and rho are smooth; sigma is singular at (a^2,rho)")
print("PASS: all three divisor pairs are nontransverse at the same cluster")
print("PASS: the triple-intersection ideal is (a^3,rho)")
print("PASS: local forms give orders B^2=2 and rho*sigma=3")
print("PASS: the p3 numerator coefficient has only order one there")
print("NOTE: the simple normal-crossings two-mask viability test fails")
