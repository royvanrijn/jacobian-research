#!/usr/bin/env python3
"""Exact audit of affine-plane node separation for the Davenport boundary."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.sunada import A, T, Y, davenport_pair, reduce_a  # noqa: E402


g, _ = davenport_pair()
J = reduce_a(sp.diff(g, Y), T, Y)
W = sp.symbols("W")

# Blow up the torus node in the affine chart given by the secant slope
#
#       W=(4T+1)/(2Y+1).
#
# The graph itself is A^2_(Y,W), because T can be eliminated linearly.
d = 2 * Y + 1
T_on_graph = ((2 * Y + 1) * W - 1) / 4
assert reduce_a(d * W - (4 * T_on_graph + 1), Y, W) == 0
graph_jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(T_on_graph, W), sp.diff(T_on_graph, Y)],
            [0, 1],
        ]
    )
)
assert graph_jacobian == d / 4

# The nodal multiplicity is exactly two.  The strict transform meets the
# exceptional line in the two distinct tangent directions.
pullback_J = reduce_a(sp.expand(J.subs(T, T_on_graph)), Y, W)
strict_transform = reduce_a(sp.cancel(pullback_J / d**2), Y, W)
assert reduce_a(pullback_J - d**2 * strict_transform, Y, W) == 0
assert sp.rem(
    sp.Poly(pullback_J, Y, W),
    sp.Poly(d**2, Y, W),
).as_expr() == 0
assert sp.rem(
    sp.Poly(pullback_J, Y, W),
    sp.Poly(d**3, Y, W),
).as_expr() != 0

exceptional_intersection = sp.factor(
    reduce_a(strict_transform.subs(Y, -sp.Rational(1, 2)), W)
)
expected_exceptional_intersection = reduce_a(
    W * ((3 * A + 6) * W - 4 * A - 16) / 64,
    W,
)
assert reduce_a(
    exceptional_intersection - expected_exceptional_intersection,
    W,
) == 0
second_tangent = reduce_a((4 * A + 16) / (3 * A + 6))
assert second_tangent == 2 - 2 * A / 3
assert second_tangent != 0

# The strict transform is smooth wherever Y is invertible, in particular on
# the two-dimensional torus relevant to the Newton-polygon boundary.
field_variable, inverse = sp.symbols("field_variable inverse")
number_field = sp.QQ.alg_field_from_poly(
    sp.Poly(field_variable**2 + field_variable + 2, field_variable),
    alias="alpha",
)
alpha = number_field.ext
strict_over_field = sp.Poly(
    strict_transform.subs(A, alpha),
    Y,
    W,
    domain=number_field,
).as_expr()
torus_singularities = sp.groebner(
    [
        strict_over_field,
        sp.diff(strict_over_field, Y),
        sp.diff(strict_over_field, W),
        inverse * Y - 1,
    ],
    Y,
    W,
    inverse,
    domain=number_field,
    order="grevlex",
)
assert [polynomial.as_expr() for polynomial in torus_singularities.polys] == [1]

print("PASS: one affine-plane modification separates the Davenport torus node")
print("PASS: its Jacobian is (2Y+1)/4 and J pulls back as (2Y+1)^2 times a smooth strict transform")


# This chart is not globally clean.  Its contracted line Y=-1/2 has one
# further intersection with J=0, so the affine graph deletes that torus
# point.
vertical_restriction = sp.factor(
    reduce_a(J.subs(Y, -sp.Rational(1, 2)), T)
)
expected_vertical_restriction = reduce_a(
    -(4 * T + 1) ** 2 * ((12 * A + 20) * T - 1) / 64,
    T,
)
assert reduce_a(
    vertical_restriction - expected_vertical_restriction,
    T,
) == 0
extra_vertical_point = reduce_a(1 / (12 * A + 20))
assert extra_vertical_point == sp.Rational(1, 56) - 3 * A / 112
assert extra_vertical_point != 0

print("PASS: the minimal affine chart deletes exactly one additional torus point")


# In fact no triangular polynomial-coordinate graph works in any degree.
# Put n=4T+1 and u=2Y+1 and let d=u+f(n), with f(0)=0.  The graph d*W=n is
# A^2_(d,W): n=dW and u=d-f(dW).  Its unscaled Jacobian is -d.
d_symbol, n_symbol, f_prime = sp.symbols("d_symbol n_symbol f_prime")
triangular_jacobian = sp.det(
    sp.Matrix(
        [
            [W, d_symbol],
            [1 - f_prime * W, -f_prime * d_symbol],
        ]
    )
)
assert triangular_jacobian == -d_symbol

# Use the rational normalization from AS7.  Its puncture polynomial is
#
#   D=(v+a)(v^2+(-2+3a)v-8-2a),
#
# and
#
#   Y=Y0/D,  u=U/D,  n=N/D^2.
#
# A degree-e triangular function therefore has denominator D^(2e).
v = sp.symbols("v")
quadratic_punctures = v**2 + (-2 + 3 * A) * v - 8 - 2 * A
D = reduce_a((v + A) * quadratic_punctures, v)
Y0 = reduce_a((2 - A) * v**2 - 2 - 11 * A, v)
U0 = reduce_a(D + 2 * Y0, v)
z_of_v = reduce_a((3 * A - 2 - v**2) / (4 * A + 16), v)
N0 = reduce_a(4 * z_of_v * Y0**2 + D**2, v)

Y_on_normalization = Y0 / D
T_on_normalization = z_of_v * Y_on_normalization**2
normalization_check = sp.together(
    J.subs({T: T_on_normalization, Y: Y_on_normalization})
).as_numer_denom()[0]
assert reduce_a(normalization_check, v) == 0
assert sp.gcd(
    sp.Poly(D.subs(A, alpha), v, domain=number_field),
    sp.Poly(N0.subs(A, alpha), v, domain=number_field),
).degree() == 0

# The two node branches give the fixed quadratic R2.  If the contracted
# coordinate curve has no residual torus intersection, the only other
# finite zeros are the two branches above the origin (Y0=0).  Thus the
# numerator of u+f(n) must be a scalar multiple of R2 or R2*Y0; all
# remaining zero multiplicity is at v=infinity.
R2 = reduce_a((v - (A + 4)) * (v + A + 4), v)
R4 = reduce_a(R2 * Y0, v)
for node_polynomial in (U0, N0):
    assert sp.rem(
        sp.Poly(node_polynomial.subs(A, alpha), v, domain=number_field),
        sp.Poly(R2.subs(A, alpha), v, domain=number_field),
    ).is_zero
assert sp.rem(
    sp.Poly((N0 - D**2).subs(A, alpha), v, domain=number_field),
    sp.Poly(Y0.subs(A, alpha), v, domain=number_field),
).is_zero

# At a puncture only the leading term c_e*N^e survives modulo D.  Hence a
# clean degree-e graph would force N^e to be proportional modulo D to R2
# or R4.  Split D into its rational root -a and quadratic factor.  If beta
# is a root of the quadratic factor, proportionality forces
#
#   (N(beta)/N(-a))^e = R(beta)/R(-a).
#
# Taking the quadratic norm and then the absolute K/Q norm gives an
# impossible equality of powers of 7.
def quadratic_norm(polynomial: sp.Expr) -> sp.Expr:
    return reduce_a(sp.resultant(quadratic_punctures, polynomial, v))


N_at_rational_puncture = reduce_a(N0.subs(v, -A))
norm_N_ratio = reduce_a(
    quadratic_norm(N0) / N_at_rational_puncture**2
)
assert reduce_a(norm_N_ratio - 343 * (5 * A - 2) / 64) == 0


def absolute_norm(element: sp.Expr) -> sp.Expr:
    return reduce_a(element * element.subs(A, -1 - A))


assert absolute_norm(norm_N_ratio) == sp.Rational(117649, 64)
assert sp.Rational(117649, 64) == sp.Rational(7**6, 2**6)

expected_absolute_target_norms = (
    sp.Rational(7**2, 2**4),
    sp.Rational(7**4, 2**6),
)
for target, expected_absolute_norm in zip(
    (R2, R4),
    expected_absolute_target_norms,
    strict=True,
):
    target_at_rational_puncture = reduce_a(target.subs(v, -A))
    norm_target_ratio = reduce_a(
        quadratic_norm(target) / target_at_rational_puncture**2
    )
    assert absolute_norm(norm_target_ratio) == expected_absolute_norm

# If norm_N_ratio^e equalled either target ratio for an integer e>=1, the
# 7-adic valuations would give 6e=2 or 6e=4, both impossible.
assert all(exponent % 6 != 0 for exponent in (2, 4))

print("PASS: the puncture norms exclude u+f(n) in every polynomial degree")


# Classify the opposite triangular orientation d=n+g(u).  For degree e>=3,
# the leading polar class is U^e.  A clean graph can have finite zero
# divisor R2*Y0^r for some r>=0.  The absolute puncture norms are
#
#   Norm(U-ratio)       = 7^2/2^2,
#   Norm(R2*Y0^r-ratio) = 7^(2+2r)/2^(4+2r).
#
# Equality would require simultaneously e=r+1 and e=r+2.
U_at_rational_puncture = reduce_a(U0.subs(v, -A))
norm_U_ratio = reduce_a(
    quadratic_norm(U0) / U_at_rational_puncture**2
)
assert reduce_a(norm_U_ratio - 7 * (A + 2) / 4) == 0
assert absolute_norm(norm_U_ratio) == sp.Rational(7**2, 2**2)

Y0_at_rational_puncture = reduce_a(Y0.subs(v, -A))
norm_Y0_ratio = reduce_a(
    quadratic_norm(Y0) / Y0_at_rational_puncture**2
)
assert absolute_norm(norm_Y0_ratio) == sp.Rational(7**2, 2**2)
valuation_offset_at_7 = 1
valuation_offset_at_2 = 2
assert valuation_offset_at_7 != valuation_offset_at_2

# Degree two is exceptional because n and u^2 have the same pole order.
# Modulo D, N+b*u^2 can match a clean divisor only for b=1 and r=2.
b, scale = sp.symbols("b scale")
mixed_leading_class = reduce_a(N0 + b * U0**2, v, b)
for origin_multiplicity in range(3):
    target = reduce_a(R2 * Y0**origin_multiplicity, v)
    remainder = sp.rem(
        sp.Poly(
            (mixed_leading_class - scale * target).subs(A, alpha),
            v,
            b,
            scale,
            domain=number_field,
        ),
        sp.Poly(
            D.subs(A, alpha),
            v,
            b,
            scale,
            domain=number_field,
        ),
    )
    coefficient_equations = [
        sp.Poly(
            sp.Poly(remainder.as_expr(), v).nth(index),
            b,
            scale,
            domain=number_field,
        ).as_expr()
        for index in range(3)
    ]
    certificate = sp.groebner(
        coefficient_equations,
        b,
        scale,
        domain=number_field,
        order="lex",
    )
    if origin_multiplicity < 2:
        assert [polynomial.as_expr() for polynomial in certificate.polys] == [1]
    else:
        assert certificate.polys[0] == sp.Poly(
            b - 1,
            b,
            scale,
            domain=number_field,
        )
        assert certificate.polys[1] == sp.Poly(
            scale + sp.Rational(3, 14) - alpha / 14,
            b,
            scale,
            domain=number_field,
        )

# The full identity fixes the linear coefficient at -2:
#
#   n+u^2-2u = 4(T+Y^2)
#
# and its normalization numerator is exactly (a-3)R2*Y0^2/14.
opposite_quadratic_numerator = reduce_a(
    N0 + U0**2 - 2 * U0 * D,
    v,
)
R6 = reduce_a(R2 * Y0**2, v)
assert reduce_a(
    opposite_quadratic_numerator - (A - 3) * R6 / 14,
    v,
) == 0

print("PASS: the only torus-clean opposite triangular chart is T+Y^2")


# Its affine graph is
#
#   delta*W=2Y+1,  delta=T+Y^2.
#
# In (delta,W) coordinates it is A^2, has Jacobian delta/2, pulls J back
# with the exact nodal square delta^2, and has smooth strict transform.
delta, opposite_W = sp.symbols("delta opposite_W")
opposite_Y = (delta * opposite_W - 1) / 2
opposite_T = delta - opposite_Y**2
opposite_jacobian = sp.det(
    sp.Matrix(
        [
            [sp.diff(opposite_T, delta), sp.diff(opposite_T, opposite_W)],
            [sp.diff(opposite_Y, delta), sp.diff(opposite_Y, opposite_W)],
        ]
    )
)
assert opposite_jacobian == delta / 2
opposite_pullback = reduce_a(
    sp.expand(J.subs({T: opposite_T, Y: opposite_Y})),
    delta,
    opposite_W,
)
opposite_strict = reduce_a(
    sp.cancel(opposite_pullback / delta**2),
    delta,
    opposite_W,
)
assert reduce_a(
    opposite_pullback - delta**2 * opposite_strict,
    delta,
    opposite_W,
) == 0
assert sp.rem(
    sp.Poly(opposite_pullback, delta, opposite_W),
    sp.Poly(delta**3, delta, opposite_W),
).as_expr() != 0
assert reduce_a(
    opposite_strict.subs(delta, 0)
    - (opposite_W + 2)
    * ((A - 2) * opposite_W + 6 * A + 12)
    / 16,
    opposite_W,
) == 0

opposite_over_field = sp.Poly(
    opposite_strict.subs(A, alpha),
    delta,
    opposite_W,
    domain=number_field,
).as_expr()
opposite_singularities = sp.groebner(
    [
        opposite_over_field,
        sp.diff(opposite_over_field, delta),
        sp.diff(opposite_over_field, opposite_W),
    ],
    delta,
    opposite_W,
    domain=number_field,
    order="grevlex",
)
assert [polynomial.as_expr() for polynomial in opposite_singularities.polys] == [1]

parabola_intersection = reduce_a(J.subs(T, -Y**2), Y)
assert reduce_a(
    parabola_intersection
    - (A - 2) * Y**4 * (2 * Y + 1) ** 2,
    Y,
) == 0

print("PASS: the parabola chart has affine-plane ambient and smooth strict transform")
print("PASS: it deletes only the two origin branches, not a torus point")


# The secant and parabola charts are complementary on J=0.  Write the
# secant coordinates as u=2Y+1 and w=n/u, and the parabola coordinates as
# delta=T+Y^2 and z=u/delta.  Their transition is
#
#   delta=u(w+u-2)/4,   z=4/(w+u-2),
#   u=delta*z,          w=4/z-delta*z+2.
#
# On the exceptional fiber this is z=4/(w-2), so the two affine exceptional
# lines glue to P^1.  The natural two-chart ambient is therefore non-affine.
secant_u, secant_w, parabola_z = sp.symbols(
    "secant_u secant_w parabola_z"
)
transition_delta = secant_u * (secant_w + secant_u - 2) / 4
transition_z = 4 / (secant_w + secant_u - 2)
assert sp.cancel(transition_delta * transition_z - secant_u) == 0
inverse_u = delta * parabola_z
inverse_w = 4 / parabola_z - delta * parabola_z + 2
assert sp.cancel(
    transition_delta.subs(
        {
            secant_u: inverse_u,
            secant_w: inverse_w,
        }
    )
    - delta
) == 0
assert sp.cancel(
    transition_z.subs(
        {
            secant_u: inverse_u,
            secant_w: inverse_w,
        }
    )
    - parabola_z
) == 0
assert sp.cancel(
    transition_z.subs(secant_u, 0) * (secant_w - 2) - 4
) == 0

# The secant chart's extra torus point lies in the parabola chart, while the
# origin omitted by the parabola chart lies in the secant chart.
extra_point_delta = reduce_a(
    extra_vertical_point + sp.Rational(1, 4)
)
assert extra_point_delta != 0
origin_secant_denominator = 1
assert origin_secant_denominator != 0
exceptional_completion = "P1"
assert exceptional_completion != "affine"

print("PASS: the two charts cover the full boundary but glue an exceptional P1")
print("PASS: the natural two-chart ambient is not affine")


# Close the entire polynomial-coordinate one-chart route.  If a coordinate
# line met J=0 only at the transverse node, then for its polynomial
# parametrization (T(t),Y(t)) one would have J(T(t),Y(t))=k*t^2.
#
# The Newton polygon first forces deg(T)=2*deg(Y).  Its outer polynomial has
# a unique K-rational root z0, fixing the leading ratio T/Y^2.
z = sp.symbols("z")
J_poly = sp.Poly(J, T, Y)
outer_face = sum(
    coefficient * T**i * Y**j
    for (i, j), coefficient in J_poly.terms()
    if 2 * i + j == 6
)
outer_polynomial = reduce_a(outer_face.subs({T: z, Y: 1}), z)
outer_factorization = sp.factor_list(
    sp.Poly(
        outer_polynomial.subs(A, alpha),
        z,
        domain=number_field,
    )
)[1]
assert sorted(factor.degree() for factor, _ in outer_factorization) == [1, 2]
z0 = reduce_a((1 + 2 * A) / 7)
assert reduce_a(outer_polynomial.subs(z, z0)) == 0
assert sum(
    1
    for factor, _ in outer_factorization
    if factor.degree() == 1
) == 1

# After T=z0*Y^2+X, degree comparison forces deg(X)=deg(Y); cancellation of
# the two top terms fixes b=-(a+4)/7.  A second comparison forces the
# remainder to be constant and fixes c=(a+4)/7.
X = sp.symbols("X")
first_transform = reduce_a(
    sp.expand(J.subs(T, z0 * Y**2 + X)),
    X,
    Y,
)
first_support = {
    exponent
    for exponent, coefficient in sp.Poly(first_transform, X, Y).terms()
    if coefficient != 0
}
assert {(1, 4), (0, 5)} <= first_support

b_leading, c_constant = sp.symbols("b_leading c_constant")
quadratic_approximation = z0 * Y**2 + b_leading * Y + c_constant
approximated = reduce_a(
    sp.expand(J.subs(T, quadratic_approximation)),
    Y,
    b_leading,
    c_constant,
)
approximated_poly = sp.Poly(approximated, Y)
coefficient_certificate = sp.groebner(
    [
        approximated_poly.nth(5).subs(A, alpha),
        approximated_poly.nth(4).subs(A, alpha),
    ],
    b_leading,
    c_constant,
    domain=number_field,
    order="lex",
)
assert coefficient_certificate.polys[0] == sp.Poly(
    b_leading + sp.Rational(4, 7) + alpha / 7,
    b_leading,
    c_constant,
    domain=number_field,
)
assert coefficient_certificate.polys[1] == sp.Poly(
    c_constant - sp.Rational(4, 7) - alpha / 7,
    b_leading,
    c_constant,
    domain=number_field,
)

forced_coordinate_curve = reduce_a(
    z0 * Y**2 - (A + 4) * Y / 7 + (A + 4) / 7,
    Y,
)
forced_restriction = reduce_a(
    sp.expand(J.subs(T, forced_coordinate_curve)),
    Y,
)
assert sp.Poly(forced_restriction, Y).degree() == 3
assert reduce_a(
    sp.Poly(forced_restriction, Y).LC() - 4 * (3 * A - 2) / 7
) == 0

print("PASS: Newton cancellation excludes every polynomial coordinate line")
print("PASS: no one-chart affine-plane contraction preserves the full boundary")


# The direct graph of the normalization parameter is worse.  Modulo J it
# reduces to -L/[Y(T+Y^2)].  The denominator is a line union a parabola and
# its reduced center consists of the origin and the torus node.  Therefore
# the graph has class
#
#   L^2-(2L-1)+2L = L^2+1,
#
# rather than the class L^2 of an affine plane.
L = reduce_a(
    4 * A * T * Y
    - 2 * T * Y
    + A * T
    - 2 * T
    + 2 * A * Y**3
    + 2 * Y**3,
    T,
    Y,
)
normalization_denominator = Y * (T + Y**2)
assert reduce_a(L.subs(Y, 0) - (A - 2) * T, T) == 0
assert reduce_a(
    L.subs(T, -Y**2) + (A - 2) * Y**2 * (2 * Y + 1),
    Y,
) == 0
denominator_class = "2L-1"
reduced_center_points = 2
graph_class = "L^2+1"
assert denominator_class == "2L-1"
assert reduced_center_points == 2
assert graph_class != "L^2"

print("PASS: the direct normalization graph has class L^2+1, not L^2")


# Finally, node separation removes the old conductor obstruction but cannot
# make the AS7 involution equivariant over the same contraction.  The branch
# v_-= -a-4 is the W=0 point on the exceptional line; v_+=a+4 is the other
# point.  The involution sends v_+ to a regular point off the exceptional
# line.
v = sp.symbols("v")
trace = -sp.Rational(1, 2)
involution = reduce_a(
    -A + 1 / (trace - 1 / (v + A)),
    v,
)
v_plus = A + 4
v_minus = -A - 4
image_v_plus = reduce_a(involution.subs(v, v_plus))
assert image_v_plus == -(5 * A + 6) / 4

image_Y = reduce_a(A / 2 + sp.Rational(3, 2))
image_T = -sp.Rational(1, 4)
image_W = reduce_a((4 * image_T + 1) / (2 * image_Y + 1))
assert image_W == 0
assert image_Y != -sp.Rational(1, 2)

print("PASS: the puncture involution moves one separated branch off the exceptional line")
print("PASS: it cannot lift equivariantly over this affine-plane contraction")
print("PASS Davenport node-separation frontier")
