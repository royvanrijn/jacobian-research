#!/usr/bin/env sage
"""Exact deconstruction of the two elliptic quotients of X_0^6(79)/w_474.

The genus-two model is Elkies's

    u^2 = 16*t^6 - 19*t^4 + 88*t^2 - 48.

This script identifies its three visible involutions with Atkin--Lehner
classes, checks the two elliptic quotient labels, and applies the
polarization criterion that selects the symmetric Kumar MW lattice H2.
"""

from pathlib import Path

from sage.all import (
    EllipticCurve,
    Matrix,
    PolynomialRing,
    QQ,
    QuaternionAlgebra,
    ZZ,
    diagonal_matrix,
)


BASE = Path(__file__).resolve().parents[1]


def hall_class(n):
    """Represent a Hall divisor of 474 by its prime-support vector."""
    return frozenset(p for p in (2, 3, 79) if n % p == 0)


def hall_product(m, n):
    """Multiply Atkin--Lehner labels modulo squares."""
    support = hall_class(m).symmetric_difference(hall_class(n))
    out = 1
    for p in support:
        out *= p
    return out


def subgroup(*generators):
    elements = {1}
    for generator in generators:
        elements |= {hall_product(x, generator) for x in tuple(elements)}
    return elements


# Padurariu--Saia independently record the quotient by <w_474> as
# (GenusAtMost2 commit 6cc368fe37aa67187783118f18d149b2b1fd6230)
#
#     y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576.
#
# Its exact coordinate change from Elkies's model is unusually small:
#
#     x = 2/t,       y = 6*u/t^3,
#     t = 2/x,       u = 4*y/(3*x^3).
#
# This also fixes the involution conventions without relying on a numerical
# genus-two isomorphism.  Clearing t^6 gives the polynomial identity below.
R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()
elkies_sextic = 16 * T**6 - 19 * T**4 + 88 * T**2 - 48
ps_sextic = lambda value: (
    -27 * value**6 + 198 * value**4 - 171 * value**2 + 576
)
assert K(T**6 * ps_sextic(2 / T)) == K(36 * elkies_sextic)

# Elkies's non-CM point.
t = QQ(14) / 13
u = QQ(64 * 251) / 13**3
f = 16 * t**6 - 19 * t**4 + 88 * t**2 - 48
assert u**2 == f
x_ps = 2 / t
y_ps = 6 * u / t**3
assert x_ps == QQ(13) / 7
assert y_ps == QQ(12048) / 343
assert y_ps**2 == ps_sextic(x_ps)

# The rational CM loci visible in Elkies's model.  In the infinity chart
# s=1/t, v=u/t^3 the equation is
#
#     v^2 = 16 - 19*s^2 + 88*s^4 - 48*s^6.
#
# Thus the two points at infinity are (s,v)=(0,+/-4).  The actions become
# alpha:(s,v)->(-s,-v), beta:(s,v)->(-s,v), h:(s,v)->(s,-v),
# so beta=w_3 fixes each infinity point.  The four affine points with
# t=+/-2, u=+/-32 form one free orbit under alpha,beta,h.
cm_affine = {(tt, uu) for tt in (QQ(-2), QQ(2)) for uu in (QQ(-32), QQ(32))}
assert all(uu**2 == 16 * tt**6 - 19 * tt**4 + 88 * tt**2 - 48 for tt, uu in cm_affine)
alpha = lambda point: (-point[0], point[1])
beta = lambda point: (-point[0], -point[1])
hyperelliptic = lambda point: (point[0], -point[1])
seed_cm = (QQ(2), QQ(32))
assert {seed_cm, alpha(seed_cm), beta(seed_cm), hyperelliptic(seed_cm)} == cm_affine
cm_infinity = {(QQ(0), QQ(-4)), (QQ(0), QQ(4))}
beta_infinity = lambda point: (-point[0], point[1])
assert {beta_infinity(point) for point in cm_infinity} == cm_infinity
assert all(beta_infinity(point) == point for point in cm_infinity)

# Under the Padurariu--Saia coordinate, the points at infinity of Elkies's
# model become (x,y)=(0,+/-24), while the t=+/-2 CM orbit becomes
# (x,y)=(+/-1,+/-24).  The three involutions become
#
#     alpha: (x,y) -> (-x,-y),
#     beta:  (x,y) -> (-x, y),
#     h:     (x,y) -> ( x,-y).
ps_cm_infinity = {(QQ(0), QQ(-24)), (QQ(0), QQ(24))}
ps_cm_affine = {
    (2 / tt, 6 * uu / tt**3)
    for tt, uu in cm_affine
}
assert ps_cm_affine == {
    (QQ(xx), QQ(yy))
    for xx in (-1, 1)
    for yy in (-24, 24)
}
assert all(yy**2 == ps_sextic(xx) for xx, yy in ps_cm_infinity | ps_cm_affine)

# alpha:(t,u)->(-t,u).  With x=t^2 and y=u, scaling X=16*x,
# Y=16*y gives this integral Weierstrass model.
E_alpha = EllipticCurve(QQ, [0, -19, 0, 1408, -12288])
P_alpha = E_alpha(16 * t**2, 16 * u)
E_alpha_min = E_alpha.global_minimal_model()
Q_alpha = E_alpha.isomorphism_to(E_alpha_min)(P_alpha)

# beta:(t,u)->(-t,-u).  Put x=t^2 and v=t*u, then z=1/x and
# w=v/x^2.  The cubic is w^2=-48*z^3+88*z^2-19*z+16;
# X=-48*z and Y=-48*w give the model below.
E_beta = EllipticCurve(QQ, [0, 88, 0, 912, 36864])
x = t**2
v = t * u
z = 1 / x
w = v / x**2
P_beta = E_beta(-48 * z, -48 * w)
E_beta_min = E_beta.global_minimal_model()
Q_beta = E_beta.isomorphism_to(E_beta_min)(P_beta)

assert E_alpha.conductor() == 474
assert E_beta.conductor() == 474
assert E_alpha_min.cremona_label() == "474a1"
assert E_beta_min.cremona_label() == "474b1"
assert Q_alpha == -3 * E_alpha_min.gens(proof=False)[0]
assert Q_beta == 9 * E_beta_min.gens(proof=False)[0]

# Padurariu--Saia's exact quotient table gives
#   474a1 = X_0^6(79)/<w_6,w_79>,
#   474b1 = X_0^6(79)/<w_3,w_158>.
# Their intersection is Elkies's genus-two quotient by w_474.
W_base = subgroup(474)
W_alpha = subgroup(6, 79)
W_beta = subgroup(3, 158)
W_hyperelliptic = subgroup(2, 237)
assert W_base == {1, 474}
assert W_alpha & W_beta == W_base
assert W_alpha == subgroup(474, 6)
assert W_beta == subgroup(474, 3)
assert W_hyperelliptic == subgroup(474, 2)

# Hence alpha is w_6=w_79, beta is w_3=w_158, and
# alpha*beta is the hyperelliptic involution w_2=w_237.
assert hall_product(6, 3) == 2

# The symmetric Kumar height lattice has determinant 474.  The extra
# polarization involution is the 474=2*237 factorization: (-474,2) and
# (-474,237) are precisely the quaternion algebra ramified at 2 and 3.
H2 = Matrix(QQ, [[4, 0], [0, QQ(237) / 2]])
assert H2.det() == 474
assert QuaternionAlgebra(QQ, -474, 2).ramified_primes() == [2, 3]
assert QuaternionAlgebra(QQ, -474, 237).ramified_primes() == [2, 3]

# There is also an intrinsic lattice check of the label.  In the pinned H2
# frame, changing the sign of the height-4 section is an integral isometry.
# On the cyclic discriminant group Z/948 it acts by 475, which is -1 on the
# 2-primary part and +1 on the 3- and 79-primary parts: exactly the w_2
# action.  Multiplication by the global -1 gives the equivalent w_237 action.
frame_path = BASE / "data/fibrations/kumar_e7e8_mw2_frame_2.txt"
frame = Matrix(
    ZZ,
    [
        [ZZ(value) for value in line.split()]
        for line in frame_path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ],
)
flip_height4 = diagonal_matrix(ZZ, [1] * 15 + [-1, 1])
assert flip_height4.transpose() * frame * flip_height4 == frame
smith, left, _ = frame.smith_form()
assert list(smith.diagonal()) == [1] * 16 + [948]
disc_action = left * flip_height4.inverse().transpose() * left.inverse()
disc_multiplier = ZZ(disc_action[16, 16]) % 948
assert disc_multiplier == 475
assert (disc_multiplier % 4, disc_multiplier % 3, disc_multiplier % 79) == (
    3,
    1,
    1,
)

print("X0679|point|t={}|u={}|verified=yes".format(t, u))
print(
    "X0679|padurariu_saia|x=2/t|y=6*u/t^3|non_cm=({}, {})|"
    "cm_infinity=(0,+-24)|cm_affine=(+-1,+-24)|verified=yes".format(
        x_ps, y_ps
    )
)
print(
    "X0679|alpha=(-t,u)|label={}|ainvs={}|point={}|multiple=-3".format(
        E_alpha_min.cremona_label(), E_alpha_min.ainvs(), Q_alpha
    )
)
print(
    "X0679|beta=(-t,-u)|label={}|ainvs={}|point={}|multiple=9".format(
        E_beta_min.cremona_label(), E_beta_min.ainvs(), Q_beta
    )
)
print(
    "X0679|atkin_lehner|alpha=w6=w79|beta=w3=w158|"
    "hyperelliptic=w2=w237"
)
print(
    "X0679|cm_geometry|infinity_points=(s=0,v=+-4)|fixed_by=w3|"
    "affine_orbit=(t=+-2,u=+-32)|orbit_size=4"
)
print(
    "X0679|kumar_anchor|height={}|det=474|disc_action={}|"
    "factor_involution=w2=w237".format(
        H2.list(), disc_multiplier
    )
)
print("X0679|status=PASS")
