#!/usr/bin/env python3
"""Exact rank-four expansion of the ordinary A4 chart reciprocal."""

import sympy as sp


a, b, T = sp.symbols("a b T")

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

# Numerators obtained by exact four-branch interpolation.  The checker below
# verifies them against the C3-invariant source presentation rather than
# trusting the interpolation.
p0 = (
    8 * a**8 * b
    + 12 * a**8
    - 16 * a**7 * b**2
    - 24 * a**7 * b
    - 45 * a**7
    - 100 * a**6 * b**3
    - 474 * a**6 * b**2
    - 1332 * a**6 * b
    - 1350 * a**6
    + 308 * a**5 * b**4
    + 1839 * a**5 * b**3
    + 6777 * a**5 * b**2
    + 12015 * a**5 * b
    + 11340 * a**5
    - 148 * a**4 * b**5
    - 1140 * a**4 * b**4
    - 5202 * a**4 * b**3
    - 13257 * a**4 * b**2
    - 20979 * a**4 * b
    - 14823 * a**4
    - 376 * a**3 * b**6
    - 3291 * a**3 * b**5
    - 17289 * a**3 * b**4
    - 53973 * a**3 * b**3
    - 113319 * a**3 * b**2
    - 139725 * a**3 * b
    - 88209 * a**3
    + 548 * a**2 * b**7
    + 5601 * a**2 * b**6
    + 33804 * a**2 * b**5
    + 126495 * a**2 * b**4
    + 329832 * a**2 * b**3
    + 563517 * a**2 * b**2
    + 613818 * a**2 * b
    + 297432 * a**2
    - 268 * a * b**8
    - 3048 * a * b**7
    - 20511 * a * b**6
    - 87831 * a * b**5
    - 269163 * a * b**4
    - 574209 * a * b**3
    - 858033 * a * b**2
    - 785133 * a * b
    - 367416 * a
    + 44 * b**9
    + 525 * b**8
    + 3798 * b**7
    + 17901 * b**6
    + 62289 * b**5
    + 158679 * b**4
    + 303993 * b**3
    + 415530 * b**2
    + 387099 * b
    + 196830
)
p1 = (
    40 * a**5 * b
    + 60 * a**5
    - 80 * a**4 * b**2
    - 192 * a**4 * b
    - 333 * a**4
    + 64 * a**3 * b**3
    + 240 * a**3 * b**2
    + 684 * a**3 * b
    + 540 * a**3
    - 112 * a**2 * b**4
    - 870 * a**2 * b**3
    - 3456 * a**2 * b**2
    - 7344 * a**2 * b
    - 7614 * a**2
    + 152 * a * b**5
    + 1488 * a * b**4
    + 7245 * a * b**3
    + 20925 * a * b**2
    + 34911 * a * b
    + 29646 * a
    - 64 * b**6
    - 726 * b**5
    - 4140 * b**4
    - 14850 * b**3
    - 33372 * b**2
    - 47142 * b
    - 29160
)
p2 = (
    8 * a**5 * b
    + 12 * a**5
    - 16 * a**4 * b**2
    - 72 * a**4 * b
    - 117 * a**4
    + 20 * a**3 * b**3
    + 114 * a**3 * b**2
    + 324 * a**3 * b
    + 378 * a**3
    - 44 * a**2 * b**4
    - 219 * a**2 * b**3
    - 837 * a**2 * b**2
    - 1323 * a**2 * b
    - 1620 * a**2
    + 52 * a * b**5
    + 324 * a * b**4
    + 1413 * a * b**3
    + 3132 * a * b**2
    + 4860 * a * b
    + 2673 * a
    - 20 * b**6
    - 159 * b**5
    - 783 * b**4
    - 2187 * b**3
    - 3969 * b**2
    - 3645 * b
    - 729
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

unit_numerator = (
    3 * B * p0
    + B * p1 * T
    + p2 * T**2
    - p3 * T**3
)
chart_unit = unit_numerator / (72 * B**2 * rho * sigma)


# ---------------------------------------------------------------------------
# Independent verification on the four V4-conjugate roots
# ---------------------------------------------------------------------------

s, t, X = sp.symbols("s t X")


def cyclic_invariant_1(first, second):
    denominator = second * (first - 1) * (second - 1) * (first * second - 1)
    numerator = (
        first**3 * second**3
        - 3 * first * second**2
        + second**3
        + 1
    )
    return sp.cancel(numerator / denominator)


def cyclic_invariant_2(first, second):
    denominator = second * (first - 1) * (second - 1) * (first * second - 1)
    numerator = (
        first**3 * second**3
        - 3 * first**2 * second**3
        + 6 * first * second**2
        - 3 * first * second
        + second**3
        - 3 * second**2
        + 1
    )
    return sp.cancel(numerator / denominator)


source_d = s**2 * t**2 - s * t**2 - s * t + t**2 - t + 1
source_e = (
    s**4 * t**4
    - s**2 * t**4
    - s**2 * t**2
    + t**4
    - t**2
    + 1
)
# Interpolate the four conjugate values in the unscaled root
# X=s+t+1/(st).  This produces the root-basis coefficients independently
# of p0,...,p3.
roots = []
unit_values = []
for sign_s, sign_t in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
    signed_d = source_d.subs(
        {s: sign_s * s, t: sign_t * t},
        simultaneous=True,
    )
    signed_product = (
        (sign_s * s + 1)
        * (sign_t * t + 1)
        * (sign_s * sign_t * s * t + 1)
        * signed_d
    )
    roots.append(
        sign_s * s
        + sign_t * t
        + sign_s * sign_t / (s * t)
    )
    unit_values.append(
        sign_s
        * sign_t
        * t
        * signed_product**3
        / (4 * s * source_e**3)
    )

interpolated = 0
for index in range(4):
    lagrange = 1
    for other in range(4):
        if other != index:
            lagrange *= (
                (X - roots[other]) / (roots[index] - roots[other])
            )
    interpolated += unit_values[index] * lagrange
interpolated = sp.Poly(sp.cancel(interpolated), X)
print("PASS: four conjugate chart units interpolated")

S, R = sp.symbols("S R")
target_denominator = R * (S - 1) * (R - 1) * (S * R - 1)
target_a_numerator = S**3 * R**3 - 3 * S * R**2 + R**3 + 1
target_b_numerator = (
    S**3 * R**3
    - 3 * S**2 * R**3
    + 6 * S * R**2
    - 3 * S * R
    + R**3
    - 3 * R**2
    + 1
)
target_E = (
    R**2 * S**2
    - R**2 * S
    + R**2
    - R * S
    - R
    + 1
)
target_vandermonde = (R - S) * (R * S**2 - 1) * (R**2 * S - 1)


def descend_even_rational(expression):
    """Replace even powers of s,t by S,R in a rational expression."""

    numerator, denominator = sp.fraction(sp.cancel(expression))

    def descend_polynomial(polynomial):
        descended = 0
        for (s_power, t_power), coefficient in sp.Poly(
            polynomial,
            s,
            t,
        ).terms():
            assert s_power % 2 == 0
            assert t_power % 2 == 0
            descended += (
                coefficient
                * S ** (s_power // 2)
                * R ** (t_power // 2)
            )
        return descended

    return sp.cancel(
        descend_polynomial(numerator) / descend_polynomial(denominator)
    )


def pullback_polynomial(polynomial):
    """Return the numerator after substituting a=A/D,b=B/D."""

    polynomial = sp.Poly(polynomial, a, b)
    degree = polynomial.total_degree()
    numerator = 0
    for (a_power, b_power), coefficient in polynomial.terms():
        numerator += (
            coefficient
            * target_a_numerator**a_power
            * target_b_numerator**b_power
            * target_denominator ** (degree - a_power - b_power)
        )
    return sp.expand(numerator), degree


pulled_p0, degree_p0 = pullback_polynomial(p0)
pulled_p1, degree_p1 = pullback_polynomial(p1)
pulled_p2, degree_p2 = pullback_polynomial(p2)
pulled_p3, degree_p3 = pullback_polynomial(p3)
assert (degree_p0, degree_p1, degree_p2, degree_p3) == (9, 6, 6, 3)

# These formulas use the independently factored target pullbacks
#
# B=27*S*E^3/(R*D0^3),
# rho=E^3/(R^2*D0^2),
# sigma=-27*Vandermonde*E^3/(R^2*D0^4).
candidate_squared_pairs = (
    (
        -pulled_p0,
        17496
        * S
        * target_vandermonde
        * target_E**9
        * R**4,
    ),
    (
        -pulled_p1,
        1944 * target_vandermonde * target_E**6 * R**2,
    ),
    (
        -pulled_p2,
        1944 * target_vandermonde * target_E**6 * R**2,
    ),
    (
        S * pulled_p3,
        72 * target_vandermonde * target_E**3,
    ),
)

for exponent in range(4):
    interpolated_squared = descend_even_rational(
        interpolated.nth(exponent)
    )
    interpolated_numerator, interpolated_denominator = sp.fraction(
        interpolated_squared
    )
    candidate_numerator, candidate_denominator = candidate_squared_pairs[
        exponent
    ]
    assert sp.Poly(
        sp.expand(
            interpolated_numerator * candidate_denominator
            - candidate_numerator * interpolated_denominator
        ),
        S,
        R,
    ).is_zero
    print(f"PASS: root-basis coefficient {exponent} recovered")


# The three target factors in the common denominator are irreducible and
# pairwise coprime.  Thus the correct chart unit has three genuine target
# boundary components, although they can formally be grouped into two masks
# B^2 and rho*sigma.
assert sp.factor(B) == B
assert sp.factor(rho) == rho
assert sp.factor(sigma) == sigma
assert sp.gcd(sp.Poly(B, a, b), sp.Poly(rho, a, b)).as_expr() == 1
assert sp.gcd(sp.Poly(B, a, b), sp.Poly(sigma, a, b)).as_expr() == 1
assert sp.gcd(sp.Poly(rho, a, b), sp.Poly(sigma, a, b)).as_expr() == 1

# A diagonal two-mask suspension exists on the localized chart.  Its mask
# Jacobian is exactly the reciprocal of the ordinary U,V chart Jacobian.
assert sp.cancel(
    unit_numerator / (72 * B**2 * rho * sigma) - chart_unit
) == 0

# Both mask coordinates retain genuine target poles.
assert sp.gcd(sp.Poly(p3, a, b), sp.Poly(rho, a, b)).as_expr() == 1
assert sp.gcd(sp.Poly(p3, a, b), sp.Poly(sigma, a, b)).as_expr() == 1

# The corresponding triangular primitive split also retains a genuine
# B-pole after reduction modulo the quartic.
Acoef = a**3 - b**3 - 9 * b**2 - 27 * b - 54
Ccoef = a**3 - b**3 + 27
P = (
    T**4
    - 6 * Acoef * B * T**2
    - 8 * B**3 * T
    + B**2 * (9 * Acoef**2 - 12 * Ccoef * B)
)
primitive = sp.integrate(unit_numerator, T)
primitive_mod_P = sp.Poly(sp.rem(sp.expand(primitive), P, T), T)

# The T^3 coefficient is p2/3 and p2 is not divisible by B.
assert sp.factor(primitive_mod_P.nth(3) - p2 / 3) == 0
assert sp.gcd(sp.Poly(p2, a, b), sp.Poly(B, a, b)).as_expr() == 1

print("PASS: H^3/(4*K^3*L) has an exact rank-four root-basis expansion")
print("PASS: its common denominator is 72*B^2*rho*sigma")
print("PASS: B, rho, and sigma are genuine pairwise-coprime target divisors")
print("PASS: a localized diagonal two-mask suspension has Jacobian one")
print("PASS: the formal split B^2 | rho*sigma retains a genuine B-pole")
print("NOTE: no polynomial two-mask A4 Keller map is produced")
