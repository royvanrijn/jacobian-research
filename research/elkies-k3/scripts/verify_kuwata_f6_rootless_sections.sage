#!/usr/bin/env sage
"""Exact checks for a rootless Kuwata F^(6) K3 and four line sections.

The checked equation is

  y^2 = x^3 + 909792*t^4*x
        + 1224440064*t^12 + 308629440*t^6 - 25509168.

It is the F^(6) surface attached (in the convention of Kuwata's F_n
construction) to the 2-isogenous non-CM pair

  E: y^2 = x^3 + 54*x - 189,
  F: y^2 = x^3 - 351*x - 1890.

Kuwata's rank formula therefore gives geometric MW rank 17.  This script
does not attempt to reprove that theorem or assert rank 17 over QQ(t): its
purpose is to check the equation, semistability/rootlessness, and explicit
line sections produced by the cubic-surface F^(3)-twist construction.

The four sections are R_3,...,R_6 in the Kumar--Kuwata normalization after
an exact Legendre change of variables.  Two are QQ(t)-rational; two are over
QQ(sqrt(-3))(t).
"""

QQx.<x> = PolynomialRing(QQ)
K.<r> = NumberField(x^2 + 3)
Ru.<u> = PolynomialRing(K)


def require(label, condition):
    if not condition:
        raise RuntimeError("FAILED: {}".format(label))
    print("PASS: {}".format(label))


# The original rational 2-isogeny, retained as an exact provenance check.
# E0 -> F0 is (x,y) |-> (x+1+1/x, y*(1-1/x^2)).
RT.<T> = PolynomialRing(QQ)
FT = FractionField(RT)
T = FT(T)
source_rhs = T^3 + T^2 + T
image_x = T + 1 + 1/T
image_y = (1 - 1/T^2)
require(
    "the displayed degree-two map sends E0 to F0",
    image_y^2 * source_rhs
    == image_x^3 - 2 * image_x^2 - 3 * image_x,
)

# Short models used by the F^(6) equation.
a, b = K(54), K(-189)
c, d = K(-351), K(-1890)
delta_E = -16 * (4 * a^3 + 27 * b^2)
delta_F = -16 * (4 * c^3 + 27 * d^2)
require(
    "E0 is the displayed short model E after x=9*x0+3, y=27*y0",
    27^2 * source_rhs == (9*T + 3)^3 + a*(9*T + 3) + b,
)
require(
    "F0 is the displayed short model F after x=9*x0-6, y=27*y0",
    27^2 * (T^3 - 2*T^2 - 3*T)
    == (9*T - 6)^3 + c*(9*T - 6) + d,
)
require("Delta(E) = -25509168", delta_E == -25509168)
require("Delta(F) = 1224440064", delta_F == 1224440064)
require(
    "j(E)=2048/3 is nonintegral",
    1728 * 4*a^3 / (4*a^3 + 27*b^2) == QQ(2048)/3,
)
require(
    "j(F)=35152/9 is nonintegral",
    1728 * 4*c^3 / (4*c^3 + 27*d^2) == QQ(35152)/9,
)

A = -48 * a * c
B = delta_F * u^12 + 864 * b * d * u^6 + delta_E
require("F^(6) x coefficient", A == 909792)
require("F^(6) constant coefficient", B == 1224440064*u^12 + 308629440*u^6 - 25509168)

# Rootlessness: the finite discriminant is squarefree, and the standard
# s=1/u chart is smooth at infinity.  For y^2=x^3+A*u^4*x+B, the cubic
# discriminant differs from 4(Au^4)^3+27B^2 only by the nonzero factor -16.
finite_discriminant = -16 * (4 * (A*u^4)^3 + 27 * B^2)
require("finite discriminant has degree 24", finite_discriminant.degree() == 24)
require(
    "all 24 finite singular fibers are simple",
    finite_discriminant.gcd(finite_discriminant.derivative()).degree() == 0,
)

# In the infinity chart x=u^4*x_inf, y=u^6*y_inf, s=1/u, the equation is
# y_inf^2=x_inf^3+A*x_inf+delta_F+O(s^6), whose discriminant is nonzero.
infinity_cubic_discriminant = -16 * (4 * A^3 + 27 * delta_F^2)
require("fiber at infinity is smooth", infinity_cubic_discriminant != 0)

# Exact sections.  The first and third are rational; the other two are their
# nontrivial 2-torsion-field companions.  They come from four lines on the
# cubic surface of the quadratic twist of F^(3).
sections = [
    (
        3888*u^4 + 756*u^2 + 324,
        244944*u^6 + 69984*u^4 + 40824*u^2 + 2916,
        "R3 (rational)",
    ),
    (
        3888*u^4 + (1134*r - 378)*u^2 - 162*r - 486,
        244944*u^6 + (104976*r - 34992)*u^4
        + (-20412*r - 61236)*u^2 - 5832*r + 8748,
        "R4 (over QQ(sqrt(-3)))",
    ),
    (
        -972*u^4 - 216*u^2 + 324,
        17496*u^6 - 17496*u^4 - 11664*u^2 + 2916,
        "R5 (rational)",
    ),
    (
        -972*u^4 + (-324*r + 108)*u^2 - 162*r - 486,
        17496*u^6 + (-26244*r + 8748)*u^4
        + (5832*r + 17496)*u^2 - 5832*r + 8748,
        "R6 (over QQ(sqrt(-3)))",
    ),
]

for section_x, section_y, label in sections:
    require(
        "{} lies on F^(6)".format(label),
        section_y^2 == section_x^3 + A*u^4*section_x + B,
    )

print("KUWATA_F6_ROOTLESS_SECTIONS: PASS")
