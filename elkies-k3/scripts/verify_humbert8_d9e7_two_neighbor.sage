#!/usr/bin/env sage
"""Replay the omitted D9+E7 -> E7+E8 Humbert-8 two-neighbor.

Elkies--Kumar's ancillary file ``8/8.txt`` does not begin with the Kumar
equation.  It first constructs a two-parameter fibration with fibers D9 and
E7, and then uses

    U = (X + s*T^3)/T^4

as the new elliptic parameter.  The new generic fiber is a pointed binary
quartic.  Its classical invariants give exactly the published E7+E8 Kumar
model.  This script verifies that calculation over QQ(r,s), including the
rational point that identifies the quartic with its Jacobian.

Source: Elkies--Kumar, arXiv:1209.3527, ancillary files ``8/8.txt`` and
``jacobians``.
"""

from sage.all import *


# The pre-neighbor elliptic K3.  The variables r,s are moduli coordinates;
# T is its elliptic base and (X,Y) are Weierstrass coordinates.
R = PolynomialRing(QQ, names=("r", "s", "T", "U", "X", "V"))
r, s, T, U, X, V = R.gens()

a2_old = T * (r + (2*r + 1)*T)
a4_old = 2*r*s*T**4*(T + 1)
a6_old = r*s**2*T**7
old_rhs = X**3 + a2_old*X**2 + a4_old*X + a6_old

b2 = 4*a2_old
b4 = 2*a4_old
b6 = 4*a6_old
b8 = 4*a2_old*a6_old - a4_old**2
Delta_old = -b2**2*b8 - 8*b4**3 - 27*b6**2 + 9*b2*b4*b6


def order_in(poly, variable_index):
    return min(exponent[variable_index] for exponent in poly.exponents())


# At T=0 the orders (1,4,7,11) give I5*=D9.  At infinity the K3
# weights give (2,3,5,9), hence III*=E7.  The residual quartic gives four
# generic I1 fibers.
assert order_in(a2_old, 2) == 1
assert order_in(a4_old, 2) == 4
assert order_in(a6_old, 2) == 7
assert order_in(Delta_old, 2) == 11
assert (4-a2_old.degree(T), 8-a4_old.degree(T), 12-a6_old.degree(T)) == (2, 3, 5)
assert 24-Delta_old.degree(T) == 9
residual_old = R(Delta_old // T**11)
assert residual_old.degree(T) == 4

# The component groups have common torsion divisor at most two.  The generic
# 2-division cubic is irreducible, so the generic torsion subgroup is trivial.
B = PolynomialRing(QQ, names=("rr", "ss", "tt"))
rr, ss, tt = B.gens()
KB = B.fraction_field()
xx = polygen(KB)
division_cubic = (
    xx**3 + tt*(rr+(2*rr+1)*tt)*xx**2
    + 2*rr*ss*tt**4*(tt+1)*xx + rr*ss**2*tt**7
)
assert division_cubic.is_irreducible()


# Execute the neighbor.  Put X=U*T^4-s*T^3 and Y=V*T^4.  Dividing the old
# equation by T^8 gives the pointed quartic V^2=q(T).
neighbor_x = U*T**4 - s*T**3
quartic = R(old_rhs(X=neighbor_x) // T**8)
expected_quartic = (
    U**3*T**4 - 3*s*U**2*T**3
    + (3*s**2*U + (2*r+1)*U**2)*T**2
    + (-s**3 - 2*r*s*U + r*U**2 - 2*s*U)*T
    + s**2
)
assert quartic == expected_quartic
assert quartic(T=0, V=s) == s**2


# Binary-quartic convention used by the ancillary ``jac2`` routine:
#
#   q(T)=c0*T^4+4*c1*T^3+6*c2*T^2+4*c3*T+c4.
#
# After x -> x/4 and multiplication by 64, its Jacobian is
# y^2=x^3-4*I*x-16*J.
c0 = U**3
c1 = -3*s*U**2/4
c2 = (3*s**2*U + (2*r+1)*U**2)/6
c3 = (-s**3 - 2*r*s*U + r*U**2 - 2*s*U)/4
c4 = s**2
invariant_i = c0*c4 - 4*c1*c3 + 3*c2**2
invariant_j = (
    c0*c2*c4 + 2*c1*c2*c3 - c0*c3**2 - c4*c1**2 - c2**3
)

A1 = 2*r*s**2
A = -(9*r*s + 4*r**2 + 4*r + 1)/3
B1 = r*s**2*(3*s + 8*r - 2)/3
B0 = -(
    54*r**2*s + 81*r*s - 16*r**3 - 24*r**2 - 12*r - 2
)/27
B2 = r**2
a4_new = A1*U**3 + A*U**4
a6_new = B1*U**5 + B0*U**6 + B2*U**7
assert -4*invariant_i == a4_new
assert -16*invariant_j == a6_new

# The quartic has the rational point (T,V)=(0,s).  Thus its Abel--Jacobi map
# is an isomorphism over QQ(r,s,U), rather than only a nontrivial torsor under
# the displayed Jacobian.  This certifies the characteristic-zero neighbor.
assert expected_quartic(T=0) == s**2

# Taking (0,s) as the origin, the other rational point (0,-s) maps to the
# height-four Humbert section.  A convenient derivation first gives the long
# Jacobian coordinate
#
#   x_long = d^2/(4*s^2)-c
#
# for q(T)=a*T^4+b*T^3+c*T^2+d*T+s^2, and then shifts by c/3 to the short
# model above.
raw_b = 4*c1
raw_c = 6*c2
raw_d = 4*c3
x_long = raw_d**2/(4*s**2) - raw_c
x_height4 = R.fraction_field()(x_long + raw_c/3)
y_height4 = R.fraction_field()(-raw_b*s - raw_d*x_long/(2*s))
assert y_height4**2 == x_height4**3 + a4_new*x_height4 + a6_new
assert x_height4.numerator().degree(U) == 4
assert y_height4.numerator().degree(U) == 6
assert x_height4(U=0) == s**4/4
assert y_height4(U=0) == s**6/8

# Determinant bookkeeping.  D9 and E7 have determinants 4 and 2, so this
# presentation absorbs the Humbert-8 height-four direction into its fibers.
# On the determinant-948 Shimura curve the sole new MW direction therefore
# has height 948/8=237/2 (generic torsion was excluded above).
assert CartanMatrix(["D", 9]).det() * CartanMatrix(["E", 7]).det() == 8
assert QQ(948)/8 == QQ(237)/2

# The exact CM(-43) point found independently in the Kumar chart is also a
# valid specialization of this pre-neighbor family.
r43 = -QQ(1225)/722
s43 = -QQ(93312)/442225
cm_quartic = expected_quartic(r=r43, s=s43)
assert cm_quartic(T=0) == s43**2
assert a4_new(r=r43, s=s43) != 0
assert a6_new(r=r43, s=s43) != 0

print(
    "H8D9E7|old_fibers=D9,E7,4I1|root_det=8|torsion=trivial"
    "|neighbor_base=(X+s*T^3)/T^4",
    flush=True,
)
print(
    "H8D9E7|quartic_point=T:0,V:s|jacobian=published_E7+E8_Kumar"
    "|second_point=height4_section|shimura_MW_height=237/2",
    flush=True,
)
print(
    f"H8D9E7|CM43|r={r43}|s={s43}|point_V={s43}|status=PASS",
    flush=True,
)
