#!/usr/bin/env python3
"""Saturate the visible LND pencil on the normalized (2,3) slice.

The previously known commuting derivations D10 and D22 span the generic
Ga^2 distribution, but that pair is not primitive along a=0.  If s is the
global D10-slice and h=D22(s), then D22-h*D10 vanishes on a=0.  This script
constructs the quotient derivation

    E = (D22-h*D10)/a

in the coordinate ring, verifies that it is a homogeneous LND, and checks
that D10 and E no longer lose rank at the sample boundary point.  Since
E(s)=0 and k[X]=ker(D10)[s], their two kernels generate k[X]; consequently
the full Derksen algebra is k[X].

This is a saturation result for the visible LND pencil, not a classification
of all homogeneous LNDs.  In particular it does not compute ML(X) or prove
flexibility.  The final check also proves that the displayed degree-two LND
of the trinomial quotient cannot lift to any global derivation of X, even
after arbitrary regular corrections in the p- and s-directions.
"""

from __future__ import annotations

import sympy as sp


T = sp.Symbol("T")
a, c, z, p, d, y, x = sp.symbols("a c z p d y x")
variables = (a, c, z, p, d, y, x)
weights = dict(zip(variables, (4, 3, 2, -3, -4, -5, -6), strict=True))

A = a * T**2 + c * T + z
B = p * T**3 + d * T**2 + y * T + x
m = a * d + c * p
resultant = sp.expand(sp.resultant(A, B, T))


def apply_derivation(
    polynomial: sp.Expr, images: tuple[sp.Expr, ...]
) -> sp.Expr:
    return sp.expand(
        sum(
            sp.diff(polynomial, variable) * image
            for variable, image in zip(variables, images, strict=True)
        )
    )


def assert_homogeneous(polynomial: sp.Expr, expected_weight: int) -> None:
    """Check the effective coefficient Gm-weight of every nonzero monomial."""

    expanded = sp.Poly(sp.expand(polynomial), *variables, domain=sp.QQ)
    for exponents, coefficient in expanded.terms():
        if coefficient == 0:
            continue
        monomial_weight = sum(
            exponent * weights[variable]
            for variable, exponent in zip(variables, exponents, strict=True)
        )
        assert monomial_weight == expected_weight


# The two previously visible commuting LNDs.
D22 = (
    sp.Integer(0),
    sp.Integer(0),
    sp.Integer(0),
    a**2,
    -a * c,
    a * z - 2 * c**2,
    -2 * c * z,
)
D10 = (
    sp.Integer(0),
    -2 * a**2,
    -a * c,
    sp.Integer(0),
    2 * a * p,
    -2 * a * d + 5 * c * p,
    -a * y + 5 * p * z,
)

slice_coordinate = (
    14 * a * p * x
    - 14 * c * d**2
    - 7 * c * p * y
    - 35 * d * p * z
    + 9 * y
) / 10

U = 4 * a * z - c**2
K = 4 * a**2 * y - 4 * c + 8 * c**2 * p - 4 * a * p * z
J = (
    8 * a**3 * x
    - 4 * a**2 * c * y
    + 4 * c**2
    - 8 * c**3 * p
    + 20 * a * c * p * z
    - 8 * a * z
)
h = sp.expand((7 * J - 12 * U) / 40)

assert sp.rem(apply_derivation(slice_coordinate, D22) - h, m - 1, d) == 0
assert sp.rem(apply_derivation(slice_coordinate, D10) - 1, m - 1, d) == 0
for invariant in (a, h):
    assert apply_derivation(invariant, D22) == 0
    assert sp.rem(apply_derivation(invariant, D10), m - 1, d) == 0

# Delta=D22-h*D10 vanishes modulo (a,m-1).  For y and x, its ambient
# representatives have a residual multiple of m-1; subtracting those
# multiples makes the divisibility by a literal in the polynomial ring.
delta_images = tuple(
    sp.expand(left - h * right)
    for left, right in zip(D22, D10, strict=True)
)
boundary_corrections = (
    sp.Integer(0),
    sp.Integer(0),
    sp.Integer(0),
    sp.Integer(0),
    sp.Integer(0),
    c**2 * (7 * c * p + 2),
    c * z * (7 * c * p + 2),
)

E_list: list[sp.Expr] = []
coefficient_domain = sp.QQ[c, z, p, d, y, x]
for delta_image, correction in zip(
    delta_images, boundary_corrections, strict=True
):
    numerator = sp.expand(delta_image - correction * (m - 1))
    quotient, remainder = sp.div(
        sp.Poly(numerator, a, domain=coefficient_domain),
        sp.Poly(a, a, domain=coefficient_domain),
    )
    assert remainder.as_expr() == 0
    E_list.append(sp.expand(quotient.as_expr()))
E = tuple(E_list)

# Thus a*E=Delta in k[X].  The chosen representatives preserve m exactly and
# preserve resultant=1 and s=0 modulo m-1.
for E_image, delta_image in zip(E, delta_images, strict=True):
    assert sp.rem(a * E_image - delta_image, m - 1, d) == 0
assert apply_derivation(m, E) == 0
assert sp.rem(apply_derivation(resultant, E), m - 1, d) == 0
assert sp.rem(apply_derivation(slice_coordinate, E), m - 1, d) == 0

# Quotient-ring reduction independently verifies preservation of both defining
# equations and local nilpotence on every algebra generator.
quotient_groebner = sp.groebner(
    [m - 1, resultant - 1],
    x,
    y,
    z,
    d,
    p,
    c,
    a,
    order="grevlex",
    domain=sp.QQ,
)


def reduce_on_X(polynomial: sp.Expr) -> sp.Expr:
    return sp.expand(quotient_groebner.reduce(sp.expand(polynomial))[1])


assert reduce_on_X(apply_derivation(m, E)) == 0
assert reduce_on_X(apply_derivation(resultant, E)) == 0
assert reduce_on_X(apply_derivation(slice_coordinate, E)) == 0

nilpotence_orders: dict[sp.Symbol, int] = {}
for variable in variables:
    iterate = variable
    for order in range(1, 9):
        iterate = reduce_on_X(apply_derivation(iterate, E))
        if iterate == 0:
            nilpotence_orders[variable] = order
            break
    else:
        raise AssertionError(f"nilpotence bound failed on {variable}")
assert nilpotence_orders == {a: 1, c: 2, z: 3, p: 2, d: 3, y: 4, x: 5}

# The effective grading is inherited from relative coefficient scaling.
assert_homogeneous(m - 1, 0)
assert_homogeneous(resultant - 1, 0)
assert_homogeneous(h, 6)
assert_homogeneous(slice_coordinate, -5)
for variable, image in zip(variables, D10, strict=True):
    if image != 0:
        assert_homogeneous(image, weights[variable] + 5)
for variable, image in zip(variables, D22, strict=True):
    if image != 0:
        assert_homogeneous(image, weights[variable] + 11)
for variable, image in zip(variables, E, strict=True):
    if image != 0:
        assert_homogeneous(image, weights[variable] + 7)

# The old pair loses rank at this point of a=0, while its primitive saturation
# does not.  E is not divisible by a again: E(z)=-2/5 at this point.
boundary_point = {a: 0, c: 1, z: 0, p: 1, d: 0, y: 0, x: -1}
assert m.subs(boundary_point) == 1
assert resultant.subs(boundary_point) == 1
assert slice_coordinate.subs(boundary_point) == 0


def value_vector(images: tuple[sp.Expr, ...]) -> sp.Matrix:
    return sp.Matrix([image.subs(boundary_point) for image in images])


old_rank = sp.Matrix.hstack(value_vector(D10), value_vector(D22)).rank()
saturated_rank = sp.Matrix.hstack(value_vector(D10), value_vector(E)).rank()
assert old_rank == 1
assert saturated_rank == 2
assert E[2].subs(boundary_point) == sp.Rational(-2, 5)
assert E[3] == a

# The quotient LND L(a)=K^2, L(U)=192*a^2, L(K)=L(J)=0 has
# L(h)=-288*a^2/5.  On X_a, an arbitrary lift is allowed to have
# P=L(p) and S=L(s).  Formula (53) gives
#
#   L(c)=-576*a^2*p/5 + 2*h*P - 4*a*K^2*s - 2*a^2*S.
#
# On a=0, regularity of d=(1-cp)/a forces P=-80*c*d.  Regularity of
# z=(U+c^2)/(4a) would then force z=2*c^2*d.  Since c is invertible and d,z
# are independent coordinates on the boundary, no global derivation can
# induce this quotient action.  Notice that S disappears from both boundary
# conditions, so a vertical D10 correction cannot repair the obstruction.
P_boundary = sp.Symbol("P_boundary")
h_boundary = -sp.Rational(2, 5) * c**2
K_boundary = 4 * c
p_boundary = 1 / c
Lc_boundary = sp.expand(2 * h_boundary * P_boundary)
d_regularity = sp.factor(
    p_boundary * Lc_boundary
    + c * P_boundary
    + d * K_boundary**2
)
assert d_regularity == c * (P_boundary + 80 * c * d) / 5
forced_P = -80 * c * d
z_regularity = sp.factor(
    2 * c * Lc_boundary - 4 * K_boundary**2 * z
)
assert sp.expand(
    z_regularity.subs(P_boundary, forced_P)
    - 64 * c**2 * (2 * c**2 * d - z)
) == 0

# Derksen conclusion: D10(s)=1 gives k[X]=ker(D10)[s] by the slice theorem,
# while E is nonzero and E(s)=0.  Hence ker(D10) and ker(E) generate k[X].
print("PASS: the effective grading has weights (4,3,2,-3,-4,-5,-6)")
print("PASS: E=(D22-h*D10)/a is a global homogeneous LND of degree 7")
print("PASS: the primitive pair D10,E has rank two at the tested a=0 point")
print("PASS: ker(D10) and ker(E) generate k[X], so HD(X)=k[X]")
print("PASS: the degree-two quotient LND has no global derivation lift to X")
print("OPEN: this saturation does not compute ML(X) or decide flexibility")
