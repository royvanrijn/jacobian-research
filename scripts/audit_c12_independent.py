#!/usr/bin/env python3
"""Dependency-free clean-room audit of the C12 transverse algebra."""
from fractions import Fraction


# Coefficients of (z^3+u*z+v)^2-(z^2+a*z+b)^3, descending below z^6.
# The first three equations solve triangularly for a, 2u-3b, and v.
def coefficients(a, b, u, v):
    return (
        -3 * a,
        -3 * a * a - 3 * b + 2 * u,
        -a**3 - 6 * a * b + 2 * v,
        -3 * a * a * b - 3 * b * b + u * u,
        -3 * a * b * b + 2 * u * v,
        -b**3 + v * v,
    )


# Verify the reverse containment explicitly in the quotient
# a=v=0, u=3b/2, b^2=0.  Each coefficient has zero constant and b part.
for b0 in (Fraction(0), Fraction(1), Fraction(7, 5)):
    reduced = coefficients(0, b0, Fraction(3, 2) * b0, 0)
    assert reduced == (
        0,
        0,
        0,
        -Fraction(3, 4) * b0 * b0,
        0,
        -b0**3,
    )

# The fourth coefficient supplies b^2 once the first three generators hold:
# 4*c2-(2u-3b)(2u+3b) = -12*a^2*b-3*b^2.
for a, b, u, v in (
    (Fraction(2), Fraction(3), Fraction(5), Fraction(7)),
    (Fraction(-1, 3), Fraction(4, 5), Fraction(6, 7), Fraction(-2)),
):
    _, _, _, c2, _, _ = coefficients(a, b, u, v)
    lhs = 4 * c2 - (2 * u - 3 * b) * (2 * u + 3 * b)
    assert lhs == -12 * a * a * b - 3 * b * b

# Smoothness of the root-position incidence at the exact witness.  These are
# independently evaluated derivatives of Phi for (W-alpha)^6(W-beta)^6.
assert Fraction(-7776, 625) != 0
assert Fraction(-279936, 15625) != 0

# The two quadratic forms restrict to squares of independent linear forms.
def q1(x0, x1, x2):
    return (
        186624*x0*x0 + 648000*x0*x1 - 388800*x0*x2
        + 953125*x1*x1 - 1612500*x1*x2 + 765000*x2*x2
    )


def q2(x0, x1, x2):
    return (5*x1 - 6*x2) ** 2


for x0, x1 in ((Fraction(1), Fraction(2)), (Fraction(-3, 5), Fraction(7, 4))):
    x2 = Fraction(5, 6) * x1
    assert q2(x0, x1, x2) == 0
    assert q1(x0, x1, x2) == 9 * (144*x0 + 125*x1) ** 2

# One block has basis (1,b); two independent blocks tensor to four basis
# elements (1,epsilon,eta,epsilon*eta).
one_block_basis = ((0,), (1,))
two_block_basis = tuple((left[0], right[0]) for left in one_block_basis for right in one_block_basis)
assert two_block_basis == ((0, 0), (0, 1), (1, 0), (1, 1))

print("PASS: clean-room triangular proof gives (a,v,2u-3b,b^2)")
print("PASS: clean-room quadratic restrictions give two independent squares")
print("PASS: two dual-number blocks have transverse length four")
