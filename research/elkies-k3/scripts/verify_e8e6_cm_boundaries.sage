#!/usr/bin/env sage
"""Verify the CM boundary divisors in the q=60 E8+E6 chart.

The pinned lattice neighbor predicts a generic II* fiber at infinity and an
IV* fiber at T=0.  The two exact CM closures enhance the latter to III* and
II*, respectively.  This script checks that these are precisely the
coefficient hyperplanes b0=0 and a0=b0=0 in the six-coefficient ambient
Weierstrass family.
"""

from sage.all import PolynomialRing, QQ


R = PolynomialRing(QQ, names=("a0", "a1", "b0", "b1", "b2", "b3", "T"))
a0, a1, b0, b1, b2, b3, T = R.gens()


def ord_T(poly):
    """T-adic order of a nonzero polynomial."""
    t_index = R.ngens() - 1
    return min(exponent[t_index] for exponent in poly.exponents())


def infinity_orders(a4, a6, discriminant):
    """Orders at infinity for a degree-(8,12,24) elliptic K3 model."""
    return (
        8 - a4.degree(T),
        12 - a6.degree(T),
        24 - discriminant.degree(T),
    )


a4 = T**3 * (a0 + a1 * T)
a6 = T**4 * (b0 + b1 * T + b2 * T**2 + b3 * T**3)
discriminant = -16 * (4 * a4**3 + 27 * a6**2)

# Generic q=60 child: IV* (E6) at zero and II* (E8) at infinity.
assert (ord_T(a4), ord_T(a6), ord_T(discriminant)) == (3, 4, 8)
assert infinity_orders(a4, a6, discriminant) == (4, 5, 10)

# Delta=-24 closure: E6 -> E7.  Generically on b0=0, a0 is nonzero.
a4_24 = a4.subs(b0=0)
a6_24 = a6.subs(b0=0)
discriminant_24 = discriminant.subs(b0=0)
assert (ord_T(a4_24), ord_T(a6_24), ord_T(discriminant_24)) == (3, 5, 9)

# Delta=-3 closure: E6 -> E8.  Generically on a0=b0=0, b1 is nonzero.
a4_3 = a4.subs(a0=0, b0=0)
a6_3 = a6.subs(a0=0, b0=0)
discriminant_3 = discriminant.subs(a0=0, b0=0)
assert (ord_T(a4_3), ord_T(a6_3), ord_T(discriminant_3)) == (4, 5, 10)

# The standard discriminant-3 Inose anchor lies on the second boundary:
#   Y^2 = X^3 + T^5*(T-1)^2.
anchor_a4 = a4.subs(a0=0, a1=0, b0=0, b1=1, b2=-2, b3=1)
anchor_a6 = a6.subs(a0=0, a1=0, b0=0, b1=1, b2=-2, b3=1)
anchor_discriminant = discriminant.subs(
    a0=0, a1=0, b0=0, b1=1, b2=-2, b3=1
)
assert anchor_a4 == 0
assert anchor_a6 == T**5 * (T - 1)**2
assert anchor_discriminant == -432 * T**10 * (T - 1)**4

print(
    "E8E6BOUNDARY|generic=II*+IV*|ord0=3,4,8|ordinf=4,5,10"
)
print("E8E6BOUNDARY|Delta=-24|condition=b0=0|ord0=3,5,9|root=E7")
print("E8E6BOUNDARY|Delta=-3|condition=a0=b0=0|ord0=4,5,10|root=E8")
print(
    "E8E6BOUNDARY|Delta=-3-anchor=T^5*(T-1)^2|"
    "fibers=II*+II*+IV|status=PASS"
)
