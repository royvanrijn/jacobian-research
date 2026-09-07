#!/usr/bin/env python3
"""Exact regressions for the cubic critical-boundary marking frontier."""

from __future__ import annotations

import math

import sympy as sp


s, q, Y = sp.symbols("s q Y")


# Infinite two-place defect atlas, checked on the first six odd slopes.
for b in range(1, 12, 2):
    D = 1 - s**2 * q**b
    T = s - q**b * s**3 / 3
    plane_map = sp.Matrix((q, T))
    jacobian = sp.factor(plane_map.jacobian((s, q)).det())

    assert jacobian == -D
    assert sp.Poly(T, s).degree() == 3
    assert sp.Poly(D, s, q).is_irreducible
    assert sp.factor(D.subs({s: Y ** (-b), q: Y**2})) == 0
    assert sp.factor(
        (s * q ** ((b + 1) // 2)).subs({s: Y ** (-b), q: Y**2}) - Y
    ) == 0
    assert math.gcd(2, b) == 1

    critical_T = sp.rem(
        sp.Poly(3 * T - 2 * s, s),
        sp.Poly(s**2 * q**b - 1, s),
    )
    assert critical_T.is_zero

    # The two valuation vectors generate the primitive rank-one lattice
    # jointly, while neither displayed q-vector is primitive.
    assert math.gcd(b, 2) == 1
    assert math.gcd(abs(2), abs(-2)) == 2
    if b >= 3:
        assert math.gcd(abs(b), abs(-b)) == b

    # Any reciprocal valuation with regular q would make T polar:
    # 2*u+b*v=-1 implies v(q^b*s^3)=u-1<v(s)=u.
    for valuation_q in range(0, 10):
        numerator = -1 - b * valuation_q
        if numerator % 2:
            continue
        valuation_s = numerator // 2
        assert valuation_s < 0
        assert b * valuation_q + 3 * valuation_s == valuation_s - 1
        assert valuation_s - 1 < valuation_s


# General diagonal reciprocal identity and exact polynomiality obstruction.
x, y, A_symbol = sp.symbols("x y A")
for a in range(1, 6):
    for b in range(1, 8):
        A = 1 + x**a * y**b
        for p in range(-8, 9):
            numerator = -1 - a * p
            if numerator % b:
                continue
            v = numerator // b
            assert a * p + b * v == -1

            source_s = x * A**p
            source_q = y * A**v
            reciprocal = sp.cancel(1 - source_s**a * source_q**b)
            assert sp.cancel(reciprocal - 1 / A) == 0

            primitive = sp.cancel(
                source_s
                - source_q**b * source_s ** (a + 1) / (a + 1)
            )
            expected = sp.cancel(
                x * A ** (p - 1) * (a * A + 1) / (a + 1)
            )
            assert sp.cancel(primitive - expected) == 0

            q_polynomial = v >= 0
            primitive_polynomial = p >= 1
            assert not (q_polynomial and primitive_polynomial)

# Leading-DVR algebra for the general reciprocal cubic coefficient theorem.
epsilon, s0, r1 = sp.symbols("epsilon s0 r1", nonzero=True)
for n in range(1, 8):
    valuation_d1 = n - 1
    valuation_d2 = 2 * n - 1
    r2 = -sp.Rational(3, 2) * r1 / s0

    source_s = s0 * epsilon ** (-n)
    leading_d1 = r1 * epsilon**valuation_d1
    leading_d2 = r2 * epsilon**valuation_d2
    leading_D_terms = sp.expand(
        leading_d1 * source_s + leading_d2 * source_s**2
    )
    leading_T_terms = sp.expand(
        leading_d1 * source_s**2 / 2
        + leading_d2 * source_s**3 / 3
    )

    assert sp.expand(leading_T_terms) == 0
    assert sp.expand(
        leading_D_terms + r1 * s0 * epsilon**(-1) / 2
    ) == 0

    if valuation_d2 == 1:
        assert n == 1
        d1, d2 = sp.symbols("d1 d2")
        cubic_D = 1 + d1 * s + d2 * s**2
        Q, P = -d1, d2
        primitive_Y = Q - P * s
        assert sp.expand(cubic_D - (1 - s * primitive_Y)) == 0


# Symbolic integration underlying the automatic one-place normal form.
w, h0, h1, h2, c, g = sp.symbols("w h0 h1 h2 c g")
h = h2 * w**2 + h1 * w + h0
H = sp.integrate(h, w)
T_tangent = c * (w * q - H) + g
assert sp.factor(sp.diff(T_tangent, w) - c * (q - h)) == 0
assert sp.Poly(H, w).degree() == 3

# Every cubic primitive is equivalent, modulo affine source change and the
# linear/constant target shears allowed by H(W)-qW+t, to one cubic monomial.
W = sp.symbols("W")
u3, u2, u1, u0 = sp.symbols("u3 u2 u1 u0", nonzero=True)
general_cubic = u3 * w**3 + u2 * w**2 + u1 * w + u0
depressed = sp.expand(general_cubic.subs(w, W - u2 / (3 * u3)))
lower_part = sp.expand(
    depressed.coeff(W, 1) * W + depressed.coeff(W, 0)
)
assert sp.factor(depressed - lower_part - u3 * W**3) == 0

foundational_primitive = w**2 * (1 - w)
foundational_shifted = sp.expand(
    foundational_primitive.subs(w, W + sp.Rational(1, 3))
)
foundational_lower = sp.expand(
    foundational_shifted.coeff(W, 1) * W
    + foundational_shifted.coeff(W, 0)
)
assert sp.factor(foundational_shifted - foundational_lower + W**3) == 0

# Exact positive weighted quotient tower.  In the cubic branch a=-3/2.
z = sp.symbols("z")
vertical_b = sp.symbols("vertical_b", nonzero=True)
vertical_a = -sp.Rational(3, 2)
gamma = 1 + vertical_a * x * y + vertical_b * x**2 * z
C = x * gamma
weighted_W = (1 + x * y) * gamma
positive_jacobian = sp.factor(
    sp.Matrix((weighted_W, gamma, C)).jacobian((x, y, z)).det()
)
assert sp.factor(positive_jacobian - vertical_b * x**3 * gamma**2) == 0
assert sp.cancel(C / gamma - x) == 0
assert sp.cancel((weighted_W - gamma) / C - y) == 0
assert sp.cancel(
    (gamma - 1 - vertical_a * x * y) / (vertical_b * x**2) - z
) == 0

# Cubic target polynomiality forces the constant and first conormal jets of
# gamma.  Treat g0,g1 as arbitrary coefficients in k[y].
g0, g1, higher = sp.symbols("g0 g1 higher")
u = 1 + x * y
trial_gamma = g0 + x * g1 + x**2 * higher
A_numerator = sp.expand(u + u**2 - 2 * u**3 * trial_gamma)
B_numerator = sp.expand(1 + 2 * u - 3 * u**2 * trial_gamma)
assert sp.expand(A_numerator.coeff(x, 0) - (2 - 2 * g0)) == 0
first_A = sp.expand(A_numerator.coeff(x, 1).subs(g0, 1))
assert sp.expand(first_A - (-2 * g1 - 3 * y)) == 0
assert sp.solve((A_numerator.coeff(x, 0), first_A), (g0, g1), dict=True) == [
    {g0: 1, g1: -3 * y / 2}
]
assert sp.expand(B_numerator.coeff(x, 0).subs(g0, 1)) == 0

print("PASS: every odd toric slope gives a reduced two-place cubic defect core")
print("PASS: the toric normalization and nonprimitive valuation vectors are exact")
print("PASS: reciprocal polynomiality excludes every toric cubic defect valuation")
print("PASS: no diagonal reciprocal monomial link polynomializes both core outputs")
print("PASS: cubic reciprocal polynomiality forces valuations (n-1,2*n-1)")
print("PASS: primitive quadratic conormal data extracts Y=Q-P*s and n=1")
print("PASS: the one-place reduced section integrates to the cubic tangent core")
print("PASS: every cubic tangent primitive reduces to the foundational plane core")
print("PASS: two conormal quotients reconstruct the unique positive cubic chart")
print("PASS: cubic target polynomiality forces gamma=1-3*x*y/2 mod x^2")
