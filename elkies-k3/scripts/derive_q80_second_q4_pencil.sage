#!/usr/bin/env sage
"""Derive the second q80-to-rootless q=4 pencil.

This script reuses the exact first-child Jacobian constructed by
``derive_q80_first_q4_pencil.sage``.  Let U0=d-1 be its fixed multiplicative
fiber.  The transported A4 correction forces the numerator to agree through
first order with the moving nodal x-center at U0.  At the I5* fiber at
infinity, the D9 chamber selects the double root of the leading cubic.  With

    v = U-U0,
    x0 = nodal x-coordinate at U0,
    x1 = derivative of the critical x-center,
    alpha = double root of the infinity cubic,

the second coordinate is

    W = (X-alpha*v^3-x1*v-x0)/v^2.

On the unrestricted four-parameter ambient family, division by ``v^4`` gives
an exact cubic quotient (a binary quartic with its leading coefficient zero):
the residual vanishes identically.  The exact lattice transport identifies
the marked rank-19 member with the pinned D7+D5/MW5 frame.  At CM24 the formula
specializes exactly to the D7+E6+3A1 branch;
the alternative infinity root gives a different D12 chamber.
"""

from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


HERE = Path(__file__).resolve().parent
load(str(HERE / "derive_q80_first_q4_pencil.sage"))


# ``jacobian_a`` and ``jacobian_b`` are the first-child short Weierstrass
# coefficients in the base variable U over K=QQ(d,p,q,e).
child_A = jacobian_a
child_B = jacobian_b
u0 = d-1
x0 = -3*child_B(u0)/(2*child_A(u0))
assert x0**3+child_A(u0)*x0+child_B(u0) == 0
assert 3*x0**2+child_A(u0) == 0
x1 = -child_A.derivative()(u0)/(6*x0)

leading_A = child_A[6]
leading_B = child_B[9]
alpha = -3*leading_B/(2*leading_A)
assert alpha**3+leading_A*alpha+leading_B == 0
assert 3*alpha**2+leading_A == 0
assert alpha == 3

KW = PolynomialRing(K, "W")
W = KW.gen()
local = PolynomialRing(KW, "v")
v = local.gen()


def shift_base(polynomial):
    return sum(KW(coefficient)*(v+u0)**degree for degree, coefficient in enumerate(polynomial.list()))


A_v = shift_base(child_A)
B_v = shift_base(child_B)
f = -alpha*v**3-x1*v-x0
X_substitution = v**2*W-f
second_curve, remainder = (
    X_substitution**3+A_v*X_substitution+B_v
).quo_rem(v**4)
assert second_curve.degree() <= 5

assert second_curve.degree() == 3
assert remainder == 0

def specialize_scalar(value):
    value = K(value)
    return QQ(value.numerator().subs(cm24))/QQ(value.denominator().subs(cm24))


cm_W = PolynomialRing(QQ, "W")
cm_local = PolynomialRing(cm_W, "v")
cm_v = cm_local.gen()
cm_quintic = cm_local(
    [cm_W([specialize_scalar(value) for value in KW(coefficient).list()]) for coefficient in second_curve.list()]
)
cm_remainder = cm_local(
    [cm_W([specialize_scalar(value) for value in KW(coefficient).list()]) for coefficient in remainder.list()]
)
assert cm_remainder == 0
cm_quartic = cm_quintic
assert cm_quartic.degree() == 3
expected_cm = (
    (9*cm_W.gen()**2+243*cm_W.gen()+QQ(6561)/4)*cm_v**3
    +(cm_W.gen()**3-QQ(2187)/4*cm_W.gen()-QQ(19683)/4)*cm_v**2
    +(QQ(81)/2*cm_W.gen()**2+QQ(6561)/8*cm_W.gen()+QQ(59049)/16)*cm_v
    -QQ(243)/2*cm_W.gen()**2-QQ(6561)/4*cm_W.gen()-QQ(177147)/32
)
assert cm_quartic == expected_cm
assert specialize_scalar(alpha) == 3
assert specialize_scalar(x0) == -QQ(81)/2
assert specialize_scalar(x1) == QQ(27)/2

q0, q1, q2, q3, q4 = [cm_quartic[index] for index in range(5)]
cm_i = 12*q4*q0-3*q3*q1+q2**2
cm_j = (
    72*q4*q2*q0+9*q3*q2*q1-27*q4*q1**2
    -27*q3**2*q0-2*q2**3
)
cm_A2 = -27*cm_i
cm_B2 = -27*cm_j
cm_Delta2 = 4*cm_A2**3+27*cm_B2**2


def valuation_at(polynomial, factor):
    valuation = 0
    while polynomial % factor == 0:
        polynomial //= factor
        valuation += 1
    return valuation


cm_iv_star = cm_W.gen()+QQ(27)/2
assert (
    valuation_at(cm_A2, cm_iv_star),
    valuation_at(cm_B2, cm_iv_star),
    valuation_at(cm_Delta2, cm_iv_star),
) == (3, 4, 8)
cm_factors = tuple(
    (str(factor.monic()), factor.degree(), ZZ(exponent))
    for factor, exponent in cm_Delta2.factor()
)
assert sorted((degree, exponent) for _, degree, exponent in cm_factors) == [
    (1, 1), (1, 2), (1, 8), (2, 2)
]
assert (8-cm_A2.degree(), 12-cm_B2.degree(), 24-cm_Delta2.degree()) == (2, 3, 9)

print(
    "Q80SECONDQ4PENCIL|U0=d-1|"
    "x0=-3*B1(U0)/(2*A1(U0))|x1=-A1'(U0)/(6*x0)|alpha=3",
    flush=True,
)
print(
    "Q80SECONDQ4PENCIL|v=U-d+1|"
    f"W=(X-alpha*v^3-x1*v-x0)/v^2|ambient_quotient_degree={second_curve.degree()}|"
    "ambient_remainder=0",
    flush=True,
)
print(
    "Q80SECONDQ4PENCIL|rank19_curve_degree=3|"
    "generic_fibers_from_pinned_frame=I3*,I1*,8I1|"
    "ADE=D7+D5|MW=5",
    flush=True,
)
print(
    f"Q80SECONDQ4PENCIL|cm24_alpha=3|cm24_x0=-81/2|cm24_x1=27/2|"
    f"cm24_Delta_factors={cm_factors}|"
    "cm24_finite_special_valuations=3,4,8|"
    "cm24_fibers=I3*,IV*,3I2,I1|ADE=D7+E6+3A1|geometric_MW=2",
    flush=True,
)
print(
    f"Q80SECONDQ4PENCIL|cm24_A={cm_A2}|cm24_B={cm_B2}",
    flush=True,
)
print("Q80SECONDQ4PENCIL|status=PASS", flush=True)
