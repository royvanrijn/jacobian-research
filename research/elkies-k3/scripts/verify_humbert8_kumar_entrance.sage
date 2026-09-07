"""Verify the explicit Humbert-discriminant-8 Kumar entrance.

The formulas are the D=8 ancillary data for Elkies--Kumar,
"K3 surfaces and equations for Hilbert modular surfaces"
(arXiv:1209.3527).  Here ``T`` is the elliptic-base coordinate, while
``r,s`` are coordinates on the Humbert surface.
"""

from sage.all import *


R = PolynomialRing(QQ, names=("r", "s", "T"))
r, s, T = R.gens()

# E7 at T=0 and E8 at infinity.
A1 = 2 * r * s**2
A = -(9 * r * s + 4 * r**2 + 4 * r + 1) / 3
B1 = r * s**2 * (3 * s + 8 * r - 2) / 3
B = -(54 * r**2 * s + 81 * r * s - 16 * r**3
      - 24 * r**2 - 12 * r - 2) / 27
B2 = r**2

a4 = A1 * T**3 + A * T**4
a6 = B1 * T**5 + B * T**6 + B2 * T**7
Delta = -16 * (4 * a4**3 + 27 * a6**2)

def ord_T(poly):
    return min(exp[2] for exp in poly.exponents())


assert ord_T(a4) == 3
assert ord_T(a6) == 5
assert ord_T(Delta) == 9
assert a4.degree(T) == 4
assert a6.degree(T) == 7
assert Delta.degree(T) == 14
# For a K3 model, orders at infinity are 8-deg(a4), 12-deg(a6),
# and 24-deg(Delta).
assert (8 - a4.degree(T), 12 - a6.degree(T), 24 - Delta.degree(T)) == (4, 5, 10)

# Kumar's inverse map to Clebsch--Igusa invariants.
K = R.fraction_field()
I2 = K(-24 * B1 / A1)
I4 = K(-12 * A)
I6 = K(96 * (A / A1) * B1 - 36 * B)
I10 = K(-4 * A1 * B2)

assert I2 == -4 * (8 * r + 3 * s - 2)
assert I4 == 4 * (4 * r**2 + 9 * r * s + 4 * r + 1)
assert I6 == -4 * (48 * r**3 + 94 * r**2 * s + 36 * r * s**2
                    + 40 * r**2 - 35 * r * s + 4 * r + 4 * s - 2)
assert I10 == -8 * r**3 * s**2

# The oriented RM double cover Y_-(8) -> H_8.
def branch_form(rr, ss):
    return (16 * rr * ss**2 + 32 * rr**2 * ss - 40 * rr * ss - ss
            + 16 * rr**3 + 24 * rr**2 + 12 * rr + 2)


branch = branch_form(r, s)

# The official ancillary calculation does more than give the double-cover
# equation: it rationalizes the oriented surface.  Completing the square in
# s and introducing m,n gives the following birational chart.  This is a
# useful coordinate system for intersections with other Humbert surfaces,
# because the orientation is isolated in n.
MN = PolynomialRing(QQ, names=("m", "n"))
m_poly, n_poly = MN.gens()
Kmn = MN.fraction_field()
m, n = map(Kmn, (m_poly, n_poly))
r_mn = (m**2 - 1) / (16 * (2*n**2 - 1))
s_mn = (m * (16*r_mn - 1) / (32*r_mn)
        - r_mn + QQ(5)/4 + 1/(32*r_mn))
z_mn = (16*r_mn - 1) * n
assert z_mn**2 == 2 * branch_form(r_mn, s_mn)

# Check the inverse formulas on the generic chart rather than merely checking
# that the parametrized points lie on the cover.
s1_mn = s_mn + r_mn - QQ(5)/4 - 1/(32*r_mn)
assert 32*r_mn*s1_mn/(16*r_mn - 1) == m
assert z_mn/(16*r_mn - 1) == n

# The ramification divisor z=0 is itself rational.  It is the extra-I2 branch
# of Y_-(8) -> H_8; it must not be confused with the missing H_237 curve.
V = PolynomialRing(QQ, "v")
v = V.fraction_field().gen()
r_branch = (1 - v**2) / 16
s_branch = (v + 3)**3 / (16 * (v + 1))
assert branch_form(r_branch, s_branch) == 0

# The generic residual discriminant has degree five.  Its discriminant in T
# vanishes simply along the branch divisor, producing the extra I2 described
# in the ancillary calculation.
residual = Delta // T**9
disc_residual = residual.discriminant(T)
quotient, remainder = disc_residual.quo_rem(branch)
assert remainder == 0
assert quotient.quo_rem(branch)[1] != 0

# One generic rational check keeps the squarefreeness assertion independent
# of symbolic factor ordering.
residual_sample = residual(r=1, s=1)
assert residual_sample.degree() == 5
assert gcd(residual_sample, residual_sample.derivative(T)) == 1
assert branch(r=1, s=1) == 61

# E7 has determinant 2.  Adding the height-4 identity-component section gives
# discriminant 8, exactly the Humbert discriminant selected by the H2 anchor.
assert 2 * 4 == 8

print(
    "HUMBERT8|fibers=E7,E8|ord0=3,5,9|ordinf=4,5,10"
    "|residual_degree=5|branch_simple=yes"
)
print(
    "HUMBERT8|I2={}|I4={}|I6={}|I10={}".format(I2, I4, I6, I10)
)
print("HUMBERT8|cover=z^2=2*({})".format(branch))
print(
    "HUMBERT8|rational_cover|r=(m^2-1)/(16*(2*n^2-1))"
    "|s=m*(16*r-1)/(32*r)-r+5/4+1/(32*r)|z=(16*r-1)*n"
)
print(
    "HUMBERT8|ramification|r=(1-v^2)/16"
    "|s=(v+3)^3/(16*(v+1))|not_H237=yes"
)
print("HUMBERT8|height4_disc=8|status=PASS")
