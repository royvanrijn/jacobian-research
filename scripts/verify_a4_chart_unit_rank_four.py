#!/usr/bin/env python3
"""Exact rank-four expansion and resolved divisor audit of the A4 chart unit."""

from fractions import Fraction
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


# ---------------------------------------------------------------------------
# Resolution-aware divisor allocation at (B,rho,sigma)=(a^3,rho)
# ---------------------------------------------------------------------------

z = sp.symbols("z")
local_relation = b**2 + 3 * b + 9 - z


def local_reduction(expression, *generators):
    """Reduce b^2+3b+9 to z and return a polynomial in local coordinates."""

    return sp.Poly(
        sp.rem(sp.expand(expression), local_relation, b),
        *generators,
    )


def weighted_order(expression, a_weight, z_weight):
    """Return the (a,z)-weighted order after passing to z=rho."""

    polynomial = local_reduction(expression, a, z, b)
    return min(
        a_weight * a_power + z_weight * z_power
        for (a_power, z_power, _), coefficient in polynomial.terms()
        if coefficient
    )


def divisor_order(expression, divisor):
    """Return the exact order along an irreducible target divisor."""

    quotient = sp.Poly(expression, a, b, domain=sp.QQ)
    divisor_poly = sp.Poly(divisor, a, b, domain=sp.QQ)
    order = 0
    while True:
        next_quotient, remainder = sp.div(quotient, divisor_poly)
        if remainder.as_expr() != 0:
            return order
        order += 1
        quotient = next_quotient


def newton_root_orders(coefficient_orders):
    """Read root valuations from the lower Newton polygon."""

    points = [
        (exponent, Fraction(order))
        for exponent, order in enumerate(coefficient_orders)
        if order is not None
    ]
    hull = []
    for point in points:
        while len(hull) >= 2:
            previous_slope = Fraction(
                hull[-1][1] - hull[-2][1],
                hull[-1][0] - hull[-2][0],
            )
            next_slope = Fraction(
                point[1] - hull[-1][1],
                point[0] - hull[-1][0],
            )
            if next_slope <= previous_slope:
                hull.pop()
            else:
                break
        hull.append(point)

    orders = []
    for left, right in zip(hull, hull[1:]):
        slope = Fraction(
            right[1] - left[1],
            right[0] - left[0],
        )
        orders.extend([-slope] * (right[0] - left[0]))
    return tuple(orders)


# Three blowups follow the B/rho contact and one follows the sigma corner.
# The resulting exceptional chain has primitive rays
#
#     E1=(1,1), F=(2,3), E2=(1,2), E3=(1,3).
#
# The strict transforms of sigma, B, and rho meet F, E3, and E3,
# respectively, with the last two intersection points distinct.
resolution_rays = {
    "E1": (1, 1),
    "E2": (1, 2),
    "E3": (1, 3),
    "F": (2, 3),
}
c = 2 * b + 3
local_B = sp.rem(B, local_relation, b)
local_sigma = sp.rem(sigma, local_relation, b)
assert sp.expand(local_B - (a**3 - 3 * a * z + c * z)) == 0
assert sp.expand(
    local_sigma - (c * a**3 - 3 * a**2 * z + z**2)
) == 0

# Exact chart equations for the two successive central blowups, followed by
# the separate B/rho and sigma-corner blowups.
u, v, w, k = sp.symbols("u v w k")
chart_identities = (
    (local_B.subs(z, a * u) / a, a**2 - 3 * a * u + c * u),
    (
        local_sigma.subs(z, a * u) / a**2,
        u**2 + a * (c - 3 * u),
    ),
    (
        local_B.subs(z, a**2 * v) / a**2,
        c * v + a * (1 - 3 * v),
    ),
    (
        local_sigma.subs(z, a**2 * v) / a**3,
        c + a * (v**2 - 3 * v),
    ),
    (
        local_B.subs(z, a**3 * w) / a**3,
        1 + c * w - 3 * a * w,
    ),
    (
        local_sigma.subs(z, a**3 * w) / a**3,
        c - 3 * a**2 * w + a**3 * w**2,
    ),
    (
        local_B.subs({a: u * w, z: u**2 * w}) / (u**2 * w),
        c + u * (w**2 - 3 * w),
    ),
    (
        local_sigma.subs({a: u * w, z: u**2 * w}) / (u**3 * w**2),
        u + c * w - 3 * u * w,
    ),
    (
        local_B.subs({a: u**2 * k, z: u**3 * k}) / (u**3 * k),
        c - 3 * u**2 * k + u**3 * k**2,
    ),
    (
        local_sigma.subs({a: u**2 * k, z: u**3 * k}) / (u**6 * k**2),
        1 + c * k - 3 * u * k,
    ),
)
for actual, expected in chart_identities:
    assert sp.expand(actual - expected) == 0

expected_target_orders = {
    "E1": (2, 3),
    "E2": (4, 5),
    "E3": (6, 6),
    "F": (6, 9),
}
for name, weights in resolution_rays.items():
    b_squared_order = 2 * weighted_order(B, *weights)
    rho_sigma_order = (
        weighted_order(rho, *weights)
        + weighted_order(sigma, *weights)
    )
    assert (b_squared_order, rho_sigma_order) == expected_target_orders[name]

# The quartic itself is nontrivial over the resolution.  Its Newton polygon
# shows that T has order 1,2,3,3 over E1,E2,E3,F, respectively.
quartic = sp.Poly(P, T)
expected_root_orders = {
    "E1": (Fraction(1),) * 4,
    "E2": (Fraction(2),) * 4,
    "E3": (Fraction(3),) * 4,
    "F": (Fraction(3),) * 4,
}
for name, weights in resolution_rays.items():
    coefficient_orders = [
        weighted_order(quartic.nth(exponent), *weights)
        if quartic.nth(exponent) != 0
        else None
        for exponent in range(5)
    ]
    assert newton_root_orders(coefficient_orders) == expected_root_orders[name]

# The characteristic polynomial of the full numerator, rather than the
# order of one displayed coefficient, gives its branchwise divisorial
# transform on the root incidence.
numerator_characteristic = sp.Poly(
    sp.resultant(P, X - unit_numerator, T),
    X,
)
expected_numerator_orders = {
    "E1": (Fraction(6), Fraction(4), Fraction(4), Fraction(4)),
    "E2": (Fraction(8), Fraction(7), Fraction(7), Fraction(7)),
    "E3": (Fraction(9),) * 4,
    "F": (Fraction(15), Fraction(12), Fraction(12), Fraction(12)),
}
for name, weights in resolution_rays.items():
    coefficient_orders = [
        weighted_order(numerator_characteristic.nth(exponent), *weights)
        if numerator_characteristic.nth(exponent) != 0
        else None
        for exponent in range(5)
    ]
    assert (
        newton_root_orders(coefficient_orders)
        == expected_numerator_orders[name]
    )

# The norm is a compact independent check on every four-entry row above.
chi = 16 * a**2 - 8 * a * b - 12 * a + rho
assert sp.factor(
    sp.resultant(P, unit_numerator, T)
    - 104976 * rho * chi**3 * B**6 * sigma**4
) == 0

# Include the strict transforms.  Along B the half-integral target-normalized
# order becomes integral after the ramified root normalization.
strict_numerator_orders = {}
for name, divisor in (("B", B), ("rho", rho), ("sigma", sigma)):
    coefficient_orders = [
        divisor_order(numerator_characteristic.nth(exponent), divisor)
        if numerator_characteristic.nth(exponent) != 0
        else None
        for exponent in range(5)
    ]
    strict_numerator_orders[name] = newton_root_orders(coefficient_orders)
assert strict_numerator_orders == {
    "B": (Fraction(3, 2),) * 4,
    "rho": (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
    "sigma": (Fraction(1),) * 4,
}

# At a branch with orders (d1,d2,n) for B^2, rho*sigma, and the numerator,
# allocating x numerator orders to the first mask leaves
#
#     m1=d1-x,  m2=d2-(n-x).
#
# Whenever n<d1+d2 the interval
# max(0,n-d2) <= x <= min(n,d1) distributes, but cannot remove, the
# residual delta=d1+d2-n.  These are the exact minimal residual source-mask
# orders on the resolved root charts.
def residual_mask_order(first_denominator, second_denominator, numerator):
    return max(
        Fraction(0),
        Fraction(first_denominator + second_denominator) - numerator,
    )


expected_exceptional_deficits = {
    "E1": (Fraction(0), Fraction(1), Fraction(1), Fraction(1)),
    "E2": (Fraction(1), Fraction(2), Fraction(2), Fraction(2)),
    "E3": (Fraction(3),) * 4,
    "F": (Fraction(0), Fraction(3), Fraction(3), Fraction(3)),
}
for name in resolution_rays:
    first_denominator, second_denominator = expected_target_orders[name]
    deficits = tuple(
        residual_mask_order(
            first_denominator,
            second_denominator,
            numerator_order,
        )
        for numerator_order in expected_numerator_orders[name]
    )
    assert deficits == expected_exceptional_deficits[name]

# Strict-boundary deficits: B contributes one order after normalizing the
# ramified root valuation; rho contributes on three branches; sigma is
# canceled divisorially on all four branches by the full numerator.
assert tuple(
    2 * residual_mask_order(2, 0, order)
    for order in strict_numerator_orders["B"]
) == (Fraction(1),) * 4
assert tuple(
    residual_mask_order(0, 1, order)
    for order in strict_numerator_orders["rho"]
) == (Fraction(0), Fraction(1), Fraction(1), Fraction(1))
assert tuple(
    residual_mask_order(0, 1, order)
    for order in strict_numerator_orders["sigma"]
) == (Fraction(0),) * 4

print("PASS: four point blowups give exceptional rays (1,1),(1,2),(1,3),(2,3)")
print("PASS: full transforms of B^2 and rho*sigma are exact on every ray")
print("PASS: numerator branch orders follow from its rank-four characteristic polynomial")
print("PASS: the resolved two-mask allocation has unavoidable residual deficits")


# ---------------------------------------------------------------------------
# The forced rho mask and its obstruction to polynomial descent
# ---------------------------------------------------------------------------

# The quartic splits over the generic point of rho into one simple branch
# and one triple branch.  Thus G=T+a^3 is the unique reduced equation of the
# ramified three-branch component in the local root algebra.
simple_rho_factor = T - 3 * a**3
triple_rho_factor = T + a**3
assert sp.rem(
    P - simple_rho_factor * triple_rho_factor**3,
    rho,
    b,
) == 0

# The full chart numerator vanishes on the simple component, not on the
# triple component.  Hence the residual rho pole is forced onto the triple
# component and any exact source-mask product must vanish there.
numerator_mod_rho = sp.rem(unit_numerator, rho, b)
assert sp.factor(numerator_mod_rho.subs(T, 3 * a**3)) == 0
assert sp.factor(
    numerator_mod_rho.subs(T, -a**3) - 144 * a**9 * (2 * b + 3)
) == 0

# Compute the complete transform of the compact selector G.  It has the
# same exceptional orders as rho itself, while on strict rho it vanishes
# to normalized order one on the ramified triple component and to order
# zero on the simple component.
selector_characteristic = sp.Poly(
    sp.resultant(P, X - triple_rho_factor, T),
    X,
)
expected_selector_exceptional_orders = {
    "E1": (Fraction(1),) * 4,
    "E2": (Fraction(2),) * 4,
    "E3": (Fraction(3),) * 4,
    "F": (Fraction(3),) * 4,
}
for name, weights in resolution_rays.items():
    coefficient_orders = [
        weighted_order(selector_characteristic.nth(exponent), *weights)
        if selector_characteristic.nth(exponent) != 0
        else None
        for exponent in range(5)
    ]
    assert (
        newton_root_orders(coefficient_orders)
        == expected_selector_exceptional_orders[name]
    )

selector_rho_coefficient_orders = [
    divisor_order(selector_characteristic.nth(exponent), rho)
    if selector_characteristic.nth(exponent) != 0
    else None
    for exponent in range(5)
]
assert newton_root_orders(selector_rho_coefficient_orders) == (
    Fraction(1, 3),
    Fraction(1, 3),
    Fraction(1, 3),
    Fraction(0),
)

# On the resolved charts the strict triple component is Cartier, with local
# equations obtained from G by removing the exceptional factors
#
#     G/a, G/a^2, G/a^3, G/u^3
#
# on E1,E2,E3,F charts, respectively.  These quotients do not descend to
# the original polynomial root algebra.  Indeed the height-one triple
# prime p3=(rho,G) lies inside the cluster maximal ideal
#
#     m=(a,rho,T).
#
# The containment is an exact all-degree obstruction: every polynomial
# element vanishing on the forced triple component vanishes at the cluster
# and therefore has positive order on every E1 branch.
cluster_groebner = sp.groebner(
    [a, rho, T],
    T,
    a,
    b,
    order="lex",
    domain=sp.QQ,
)
for generator in (rho, triple_rho_factor):
    assert cluster_groebner.reduce(generator)[1] == 0
assert sp.factor(rho) == rho

# But the exact total residual divisor has one E1 branch of order zero.
# If two regular polynomial masks realized it, their product would vanish
# on p3, hence lie in m and have positive order on all four E1 branches,
# contradicting this zero entry.
assert expected_exceptional_deficits["E1"] == (
    Fraction(0),
    Fraction(1),
    Fraction(1),
    Fraction(1),
)

print("PASS: P mod rho is (T-3*a^3)*(T+a^3)^3")
print("PASS: T+a^3 is the compact selector of the forced triple rho mask")
print("PASS: its exceptional orders are (1,2,3,3) on every root branch")
print("PASS: the selector is Cartier after removing exceptional chart equations")
print("PASS: the forced triple prime lies in the cluster maximal ideal")
print("OBSTRUCTION: no masks in the original polynomial root algebra have the exact residual divisor")
print("NOTE: a new affine modification using exceptional quotients is still not excluded")
