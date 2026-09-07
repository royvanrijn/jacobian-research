#!/usr/bin/env python3
"""Exact audit of proportional moving-tangent sections for the Davenport pair."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import (  # noqa: E402
    A,
    T,
    Y,
    Z,
    conjugate_a,
    davenport_pair,
    reduce_a,
)


g, h = davenport_pair()
q, d, k, s, W = sp.symbols("q d k s W")
r = q + d


def tangent_quotient(polynomial: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    """Return the polynomial second divided difference at (q,q+d)."""
    return reduce_a(
        (
            polynomial.subs(variable, r)
            - polynomial.subs(variable, q)
            - sp.diff(polynomial, variable).subs(variable, q) * d
        )
        / d**2,
        T,
        q,
        d,
        k,
    )


moving = tangent_quotient(g, Y)
proportional = sp.Poly(moving.subs(d, k * q), T)
assert proportional.degree() == 2

# A polynomial T(q) can occur only when the quadratic discriminant in T is
# a square in K[q].  Its residual quadratic discriminant factors into two
# quintics in the proportionality constant k.
discriminant_T = reduce_a(
    sp.discriminant(proportional.as_expr(), T) / q**4,
    q,
    k,
)
discriminant_k = reduce_a(sp.discriminant(discriminant_T, q), k)
factor_one = (
    k**5
    + 7 * k**4
    + 21 * k**3
    + 35 * k**2
    + 35 * k
    + 21
)
factor_two = (
    (11 * A - 7) * k**5
    + (133 * A - 105) * k**4
    + (700 * A - 560) * k**3
    + (1988 * A - 1456) * k**2
    + (2940 * A - 1960) * k
    + 1764 * A
    - 1176
)
assert reduce_a(
    discriminant_k + sp.Rational(16, 49) * factor_one * factor_two,
    k,
) == 0

field_variable = sp.symbols("field_variable")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(field_variable**2 + field_variable + 2, field_variable),
    alias="alpha",
)
alpha = number_field.ext
factorization_one = sp.factor_list(
    sp.Poly(factor_one, k, domain=number_field)
)
assert len(factorization_one[1]) == 1
assert factorization_one[1][0][0].degree() == 5
assert factorization_one[1][0][1] == 1

split_factor_two = (
    (11 * A - 7)
    * (k + 1 - sp.Rational(3, 2) * A)
    * (k + sp.Rational(48, 23) - sp.Rational(9, 23) * A)
    * (k + 2 + sp.Rational(1, 2) * A)
    * (k + 4 + A) ** 2
)
assert reduce_a(factor_two - split_factor_two, k) == 0
factorization_two = sp.factor_list(
    sp.Poly(factor_two.subs(A, alpha), k, domain=number_field)
)
assert sorted(
    (factor.degree(), multiplicity)
    for factor, multiplicity in factorization_two[1]
) == [(1, 1), (1, 1), (1, 1), (1, 2)]

# Exactly one of the two T-roots is polynomial for each of the four roots of
# factor_two.  Direct substitution is a shorter and stronger certificate
# than recording the corresponding square roots of discriminant_T.
polynomial_branches = (
    (
        sp.Rational(3, 2) * A - 1,
        -sp.Rational(3, 44) * (5 * A + 2),
    ),
    (
        sp.Rational(1, 23) * (9 * A - 48),
        sp.Rational(9, 529) * (16 * A - 1),
    ),
    (-sp.Rational(1, 2) * A - 2, sp.Rational(1, 4) * A),
    (-A - 4, -1),
)
for proportionality, coefficient in polynomial_branches:
    assert reduce_a(
        moving.subs(
            {
                d: proportionality * q,
                T: coefficient * q**2,
            }
        ),
        q,
    ) == 0
    other_root = reduce_a(
        -proportional.nth(1).subs(k, proportionality)
        / proportional.nth(2).subs(k, proportionality)
        - coefficient * q**2,
        q,
    )
    other_denominator = sp.denom(sp.cancel(other_root))
    assert sp.Poly(other_denominator, q).degree() > 0

print("PASS: the proportional square-discriminant locus factors into two quintics")
print("PASS: the first quintic is irreducible over Q(sqrt(-7))")
print("PASS: the second quintic supplies four exact quadratic tangent sections")


# The rational branch T=-q^2 is fixed by conjugation.  It therefore gives
# simultaneous polynomial tangent marks for g and h over one s-line.
T_section = -s**2
g_first = s
g_second = -(A + 3) * s
h_first = s
h_second = (A - 2) * s


def tangent_identity(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    first: sp.Expr,
    second: sp.Expr,
) -> sp.Expr:
    return reduce_a(
        polynomial.subs({variable: second, T: T_section})
        - polynomial.subs({variable: first, T: T_section})
        - sp.diff(polynomial, variable).subs(
            {variable: first, T: T_section}
        )
        * (second - first),
        s,
    )


assert tangent_identity(g, Y, g_first, g_second) == 0
assert tangent_identity(h, Z, h_first, h_second) == 0
assert reduce_a(conjugate_a(g_second, s) - h_second, s) == 0


def primitive_data(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    first: sp.Expr,
    second: sp.Expr,
):
    difference = second - first
    primitive = reduce_a(
        polynomial.subs(
            {
                variable: first + difference * W,
                T: T_section,
            }
        )
        - polynomial.subs({variable: first, T: T_section})
        - difference
        * sp.diff(polynomial, variable).subs(
            {variable: first, T: T_section}
        )
        * W,
        s,
        W,
    )
    assert reduce_a(primitive.subs(W, 0), s) == 0
    assert reduce_a(sp.diff(primitive, W).subs(W, 0), s) == 0
    assert reduce_a(primitive.subs(W, 1), s) == 0
    for index in range(8):
        divided_coefficient = reduce_a(
            sp.Poly(primitive, W).nth(index) / s**6,
            s,
        )
        assert sp.denom(sp.cancel(divided_coefficient)) == 1
    reduced_primitive = reduce_a(primitive / s**6, s, W)
    c = reduce_a(-sp.diff(reduced_primitive, W).subs(W, 1), s)
    kappa = reduce_a(
        sp.diff(reduced_primitive, W, 2).subs(W, 1) / c,
        s,
    )
    return difference, primitive, reduced_primitive, c, kappa


dg, Hg, reduced_Hg, cg, kappa_g = primitive_data(
    g, Y, g_first, g_second
)
dh, Hh, reduced_Hh, ch, kappa_h = primitive_data(
    h, Z, h_first, h_second
)

assert reduce_a(cg - 28 * (A - 10) * (2 * s + 1), s) == 0
assert reduce_a(ch + 28 * (A + 11) * (2 * s + 1), s) == 0
assert reduce_a(
    kappa_g
    - (2 * A * s + A - 50 * s - 11) / (2 * (2 * s + 1)),
    s,
) == 0
assert reduce_a(
    kappa_h
    + (2 * A * s + A + 52 * s + 12) / (2 * (2 * s + 1)),
    s,
) == 0
assert reduce_a(conjugate_a(reduced_Hg, s, W) - reduced_Hh, s, W) == 0
assert reduce_a(conjugate_a(cg, s) - ch, s) == 0

endpoint_g = 2 * A * s + A - 42 * s - 7
endpoint_h = 2 * A * s + A + 44 * s + 8
endpoint_numerator_g = 2 * A * s + A - 46 * s - 9
endpoint_numerator_h = 2 * A * s + A + 48 * s + 10
assert reduce_a(
    (kappa_g + 2) * 2 * (2 * s + 1) - endpoint_g,
    s,
) == 0
assert reduce_a(
    (kappa_h + 2) * 2 * (2 * s + 1) + endpoint_h,
    s,
) == 0
assert reduce_a(
    endpoint_g * endpoint_h
    + 2 * (928 * s**2 + 326 * s + 29),
    s,
) == 0
assert reduce_a(sp.resultant(2 * s + 1, endpoint_g, s)) == 28
assert reduce_a(sp.resultant(2 * s + 1, endpoint_h, s)) == -28
assert reduce_a(
    sp.resultant(endpoint_numerator_g, endpoint_g, s)
) == -56
assert reduce_a(
    sp.resultant(endpoint_numerator_h, endpoint_h, s)
) == -56

# In the universal weighted ansatz det(F)=b*c.  Since c is a nonunit of
# K[s], no polynomial b can make this a nonzero constant.  The uniquely
# forced a=-(1+kappa)/(2+kappa) has the coprime endpoint denominator
# certified by the two resultants -56 above.
for c_value in (cg, ch):
    assert sp.Poly(c_value, s).degree() == 1

print("PASS: conjugate g/h tangent sections share T=-s^2")
print("PASS: division by s^6 extends both weighted primitives across s=0")
print("PASS: the remaining determinant and endpoint poles are explicit")
print("PASS: both poles are forced and coprime in the universal weighted ansatz")


# On C=1 the two weighted inverse pencils recover the original Davenport
# fibers over the same (s,u)-base.  The common factor s^6 is absorbed into
# the target embedding.
u = sp.symbols("u")
for polynomial, variable, first, difference, primitive in (
    (g, Y, g_first, dg, Hg),
    (h, Z, h_first, dh, Hh),
):
    base_value = polynomial.subs({variable: first, T: T_section})
    base_derivative = sp.diff(polynomial, variable).subs(
        {variable: first, T: T_section}
    )
    target_B = -difference * base_derivative / s**6
    target_cA = (base_value - u) / s**6
    inverse = reduce_a(
        primitive / s**6 - target_B * W + target_cA,
        s,
        W,
        u,
    )
    expected = reduce_a(
        (
            polynomial.subs(
                {
                    variable: first + difference * W,
                    T: T_section,
                }
            )
            - u
        )
        / s**6,
        s,
        W,
        u,
    )
    assert inverse == expected

print("PASS: both inverse pencils recover one common Davenport Sunada base")


# A second one of the four proportional branches also lies over T=-s^2
# after a K-rational rescaling q=mu*s.  Together with the first chart its
# admissible locus covers the entire s-line.
k_two = (9 * A - 48) / 23
tau_two = 9 * (16 * A - 1) / 529
mu = (5 + 2 * A) / 3
assert reduce_a(tau_two * mu**2 + 1, s) == 0
assert reduce_a(k_two * mu - (-A - 4), s) == 0

g_two_first = mu * s
g_two_second = (1 + k_two) * g_two_first
h_two_first = conjugate_a(mu, s) * s
h_two_second = (1 + conjugate_a(k_two, s)) * h_two_first
assert tangent_identity(g, Y, g_two_first, g_two_second) == 0
assert tangent_identity(h, Z, h_two_first, h_two_second) == 0

dg_two, _, reduced_Hg_two, cg_two, kappa_g_two = primitive_data(
    g, Y, g_two_first, g_two_second
)
dh_two, _, reduced_Hh_two, ch_two, kappa_h_two = primitive_data(
    h, Z, h_two_first, h_two_second
)
endpoint_g_two = reduce_a(
    sp.diff(reduced_Hg_two, W, 2).subs(W, 1) + 2 * cg_two,
    s,
)
endpoint_h_two = reduce_a(
    sp.diff(reduced_Hh_two, W, 2).subs(W, 1) + 2 * ch_two,
    s,
)
assert reduce_a(
    cg_two
    + sp.Rational(28, 81)
    * (34 * A * s - 297 * A - 564 * s - 54),
    s,
) == 0
assert reduce_a(
    ch_two
    - sp.Rational(28, 81)
    * (34 * A * s - 297 * A + 598 * s - 243),
    s,
) == 0
assert reduce_a(
    endpoint_g_two
    + sp.Rational(28, 81)
    * (3670 * A * s + 1161 * A + 6420 * s - 1026),
    s,
) == 0
assert reduce_a(
    endpoint_h_two
    - sp.Rational(28, 81)
    * (3670 * A * s + 1161 * A - 2750 * s + 2187),
    s,
) == 0

# Up to nonzero constants, these are the complete bad polynomials of the
# two simultaneous charts.  Their gcd is one already over Q[s].
bad_four = (2 * s + 1) * (928 * s**2 + 326 * s + 29)
bad_two = (
    (1516 * s**2 - 648 * s + 729)
    * (796300 * s**2 + 3240 * s + 88209)
)
assert reduce_a(
    (2 * s + 1) * endpoint_g * endpoint_h + 2 * bad_four,
    s,
) == 0
assert reduce_a(
    cg_two * ch_two * endpoint_g_two * endpoint_h_two
    - sp.Rational(28**4 * 12544, 81**4) * bad_two,
    s,
) == 0
assert sp.gcd(sp.Poly(bad_four, s), sp.Poly(bad_two, s)).degree() == 0

# The root coordinates on the two charts differ only by a constant
# translation.  The primitives then differ by an affine polynomial, exactly
# the harmless slope/intercept change in H(W)-slope*W+intercept.
transition_g = reduce_a(
    (g_first - g_two_first + dg * W) / dg_two,
    s,
    W,
)
transition_h = reduce_a(
    (h_first - h_two_first + dh * W) / dh_two,
    s,
    W,
)
assert reduce_a(
    transition_g - W - (3 * A + 5) / 21,
    W,
) == 0
assert reduce_a(
    transition_h - W - (2 - 3 * A) / 21,
    W,
) == 0
difference_g = reduce_a(
    reduced_Hg - reduced_Hg_two.subs(W, transition_g),
    s,
    W,
)
difference_h = reduce_a(
    reduced_Hh - reduced_Hh_two.subs(W, transition_h),
    s,
    W,
)
assert reduce_a(sp.diff(difference_g, W, 2), s, W) == 0
assert reduce_a(sp.diff(difference_h, W, 2), s, W) == 0
slope_shift_g = reduce_a(sp.diff(difference_g, W), s)
slope_shift_h = reduce_a(sp.diff(difference_h, W), s)
assert slope_shift_g != 0
assert slope_shift_h != 0

# If the weighted boundary coordinate C is retained, incidence slope is
# slope=B*C.  The nonzero affine shift above can only lift as
# B_four=B_two+slope_shift/C, which is not polynomial across C=0.
B_two_coordinate, B_four_coordinate, C_coordinate = sp.symbols(
    "B_two_coordinate B_four_coordinate C_coordinate"
)
for slope_shift in (slope_shift_g, slope_shift_h):
    forced_B_four = sp.cancel(
        (
            B_two_coordinate * C_coordinate + slope_shift
        )
        / C_coordinate
    )
    forced_denominator = sp.denom(forced_B_four)
    assert sp.Poly(forced_denominator, C_coordinate).degree() == 1
    assert forced_denominator.subs(C_coordinate, 0) == 0
    remainder_at_boundary = sp.rem(
        sp.Poly(
            B_two_coordinate * C_coordinate + slope_shift,
            C_coordinate,
        ),
        sp.Poly(C_coordinate, C_coordinate),
    ).as_expr()
    assert reduce_a(remainder_at_boundary - slope_shift, s) == 0

    product_shift = sp.Poly(
        (
            B_two_coordinate * C_coordinate + slope_shift
        ).subs(A, alpha),
        B_two_coordinate,
        C_coordinate,
        domain=number_field.frac_field(s),
    )
    _, product_factors = sp.factor_list(product_shift)
    assert len(product_factors) == 1
    assert product_factors[0][0].total_degree() == 2
    assert product_factors[0][1] == 1

print("PASS: a second simultaneous chart also lies over T=-s^2")
print("PASS: the two bad polynomials are coprime and cover the full s-line")
print("PASS: their marked-root incidences glue by constant root translations")
print("PASS: C-preserving weighted gluing has a forced nonzero 1/C slope pole")
print("PASS: irreducibility forbids every weighted-product coordinate lift")
print("PASS proportional Davenport tangent sections")
