#!/usr/bin/env python3
"""Verify the PSL_2(11) action-level Keller inverse-Galois benchmark.

The symbolic half checks the corrected factorization, derivative, discriminant,
and boundary-unit ledger of Klein's two conjugate degree-eleven Shabat
polynomials.  The finite-group half constructs PSL_2(F_11) directly as
SL_2(F_11)/{+/-I}; it verifies the natural degree-twelve action, the two
degree-eleven A_5-coset actions, their Gassmann equivalence, and the
(2,3,11) Riemann--Hurwitz genera.  It also identifies the normalized
degree-five/six correspondence components as the A_4 and D_10 quotients,
checks both boundary-projection profiles, and verifies an exact X_0(11)
Weierstrass model and degree-twelve j-map.  It also proves that the
degree-five normalization has a K-point, reduces it to a sparse cubic with
j=-121, descends it to a conductor-121 Weierstrass model, and separates it
from X_0(11) by exact traces above 23.  For the genus-two component it
checks a canonical adjoint pencil, an exact even hyperelliptic model, and
the two elliptic quotients giving a (2,2)-split Jacobian.  The optional
Singular replay certifies the affine conductor, ordinary-node schemes,
normalized boundary-prime degrees, the genus-two canonical elimination, and
both quadratic Cremona transformations.  It then constructs every C5/C6 mask
as an exact normalization-module Riemann--Roch kernel.  The portable half also
checks the three full positive-genus boundary-unit lattices; optional PARI
and a separate exact residue-trace replay certify their Mordell--Weil inputs.
It finally computes the residual projection-exchange modules, constructs
effective divisor bases for the remaining mask classes, exhausts the
degree-five simple-pole divisors, certifies the descended C6 infinity
imbalance, and pins the first factor-rich normal support degrees.

The checker does not infer the monodromy group of the displayed polynomial
from its passport alone: that identification and the common regular cover are
the external Jones--Zvonkin input recorded in the accompanying note.
"""

import argparse
from collections import Counter
from itertools import product
from pathlib import Path
import shutil
import subprocess

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


# ---------------------------------------------------------------------------
# The conjugate degree-eleven Shabat polynomials over Q(sqrt(-11)).

x, y, u = sp.symbols("x y u")
number_field = sp.QQ.algebraic_field(sp.sqrt(-11))
r = number_field.ext

p1 = 2 * x + 11 - 3 * r
p2 = 2 * x**2 - (11 - 3 * r) * x - (22 + 6 * r)
p3 = x**2 + 11 * x + 55 + 9 * r
q1 = 2 * x + 5 + 3 * r
q2 = (
    2 * x**3
    + (15 - 3 * r) * x**2
    - (12 - 12 * r) * x
    + 56
    + 96 * r
)
q3 = 2 * x**3 - 18 * x**2 + (21 + 45 * r) * x - (175 + 279 * r)

normalizing_constant = 2**12 * 3**14
numerator = sp.Poly(p1**3 * p2**3 * p3, x, domain=number_field)


def polynomial_is_zero(expression, variable=x):
    """Test an exact polynomial identity over Q(sqrt(-11))."""

    return sp.Poly(sp.expand(expression), variable, domain=number_field).is_zero


# The printed Jones--Zvonkin display has a sign/parenthesis typo in its
# P-1 line.  Direct expansion gives the corrected identity below.
assert polynomial_is_zero(
    numerator.as_expr()
    - normalizing_constant
    - 2 * q1**2 * q2**2 * q3
)
assert polynomial_is_zero(
    numerator.diff().as_expr() - 11 * p1**2 * p2**2 * q1 * q2
)

# Every displayed boundary factor is irreducible over Q(sqrt(-11)).
boundary_factors = (p1, p2, p3, q1, q2, q3)
for factor in boundary_factors:
    _, factorization = sp.factor_list(factor, x, extension=r)
    assert len(factorization) == 1
    irreducible_factor, exponent = factorization[0]
    assert exponent == 1
    assert sp.degree(irreducible_factor, x) == sp.degree(factor, x)

# The discriminant is a square, as required by PSL_2(11) < A_11.
coefficient_ring = number_field.poly_ring(u)
raw_polynomial = (
    sp.Poly(p1, x, domain=coefficient_ring) ** 3
    * sp.Poly(p2, x, domain=coefficient_ring) ** 3
    * sp.Poly(p3, x, domain=coefficient_ring)
    - sp.Poly(normalizing_constant * u, x, domain=coefficient_ring)
)
normalized_discriminant = sp.cancel(
    sp.discriminant(raw_polynomial.as_expr(), x) / normalizing_constant**20
)
expected_discriminant = (
    -sp.Rational(11**11, 2**60 * 3**140) * u**6 * (u - 1) ** 4
)
assert sp.Poly(
    sp.together(normalized_discriminant - expected_discriminant),
    u,
    domain=number_field,
).is_zero
discriminant_square_root = (
    11**5 * r * u**3 * (u - 1) ** 2 / (2**30 * 3**70)
)
assert sp.Poly(
    sp.together(normalized_discriminant - discriminant_square_root**2),
    u,
    domain=number_field,
).is_zero

# The direct common-closure correspondence.  Exact factorization over
# Q(sqrt(-11)) gives bidegrees (5,5) and (6,6), matching the cross-action
# A_5 orbit sizes computed below.  The factors are recorded explicitly so
# the portable checker needs no factorization package at replay time.
half = sp.Rational(1, 2)
correspondence_five = (
    x**5
    + (-half * r + half) * x**4 * y
    - x**3 * y**2
    + x**2 * y**3
    + (-half * r - half) * x * y**4
    - y**5
    + (-half * r + sp.Rational(11, 2)) * x**4
    - 2 * r * x**3 * y
    - 2 * r * x * y**3
    + (-half * r - sp.Rational(11, 2)) * y**4
    + (16 * r + 11) * x**3
    + (sp.Rational(3, 2) * r + sp.Rational(99, 2)) * x**2 * y
    + (sp.Rational(3, 2) * r - sp.Rational(99, 2)) * x * y**2
    + (16 * r - 11) * y**3
    + (sp.Rational(165, 2) * r - sp.Rational(473, 2)) * x**2
    + 176 * r * x * y
    + (sp.Rational(165, 2) * r + sp.Rational(473, 2)) * y**2
    + (605 * r - 1430) * x
    + (605 * r + 1430) * y
    - 550 * r
)
correspondence_six = (
    x**6
    + (half * r - half) * x**5 * y
    + (-half * r - sp.Rational(3, 2)) * x**4 * y**2
    + 2 * x**3 * y**3
    + (half * r - sp.Rational(3, 2)) * x**2 * y**4
    + (-half * r - half) * x * y**5
    + y**6
    + (half * r + sp.Rational(11, 2)) * x**5
    + (sp.Rational(3, 2) * r - sp.Rational(11, 2)) * x**4 * y
    - 2 * r * x**3 * y**2
    + 2 * r * x**2 * y**3
    + (-sp.Rational(3, 2) * r - sp.Rational(11, 2)) * x * y**4
    + (-half * r + sp.Rational(11, 2)) * y**5
    + (sp.Rational(67, 2) * r - sp.Rational(77, 2)) * x**4
    + (-35 * r - 110) * x**3 * y
    + 198 * x**2 * y**2
    + (35 * r - 110) * x * y**3
    + (-sp.Rational(67, 2) * r - sp.Rational(77, 2)) * y**4
    + (66 * r - 583) * x**3
    + (-242 * r + 385) * x**2 * y
    + (242 * r + 385) * x * y**2
    + (-66 * r - 583) * y**3
    + (-902 * r - 3212) * x**2
    + 4048 * x * y
    + (902 * r - 3212) * y**2
    + (-5588 * r + 8712) * x
    + (5588 * r + 8712) * y
    + 43560
)

conjugate_numerator = (
    (2 * y + 11 + 3 * r) ** 3
    * (2 * y**2 - (11 + 3 * r) * y - (22 - 6 * r)) ** 3
    * (y**2 + 11 * y + 55 - 9 * r)
)
assert sp.Poly(
    numerator.as_expr()
    - conjugate_numerator
    - 64 * correspondence_five * correspondence_six,
    x,
    y,
    domain=number_field,
).is_zero
assert (
    sp.degree(correspondence_five, x),
    sp.degree(correspondence_five, y),
) == (5, 5)
assert (
    sp.degree(correspondence_six, x),
    sp.degree(correspondence_six, y),
) == (6, 6)

# Although the two factors have bidegrees (5,5) and (6,6) in P1 x P1,
# their affine equations have total degrees five and six.  Their top forms
# are squarefree, so their ordinary projective-plane closures are smooth at
# infinity.  The exact affine conductor calculation is optionally replayed
# in Singular below.
for correspondence, total_degree in (
    (correspondence_five, 5),
    (correspondence_six, 6),
):
    polynomial = sp.Poly(correspondence, x, y, domain=number_field)
    assert polynomial.total_degree() == total_degree
    top_form = sum(
        coefficient * x**monomial[0] * y**monomial[1]
        for monomial, coefficient in polynomial.terms()
        if sum(monomial) == total_degree
    )
    infinity_polynomial = sp.Poly(top_form.subs(y, 1), x, domain=number_field)
    assert infinity_polynomial.degree() == total_degree
    assert sp.gcd(infinity_polynomial, infinity_polynomial.diff()).degree() == 0


# ---------------------------------------------------------------------------
# The normalized quintic is an elliptic curve, not a nontrivial torsor.

# Its rational affine node is the intersection q1(x)=q1^-(y)=0.  The
# inverse image of that ordinary node on the normalization is a K-rational
# effective divisor of degree two.  The line at infinity cuts a K-rational
# effective divisor of degree five.  Since 3*2-5=1, Riemann--Roch on the
# genus-one normalization supplies a K-point.  The ordinary-node assertion
# is replayed independently by the optional Singular calculation below.
component_five_rational_node = {
    x: -(5 + 3 * r) / 2,
    y: -(5 - 3 * r) / 2,
}
assert number_field.from_sympy(
    sp.expand(correspondence_five.subs(component_five_rational_node))
) == number_field.zero
component_five_gradient_at_rational_node = (
    sp.diff(correspondence_five, x).subs(component_five_rational_node),
    sp.diff(correspondence_five, y).subs(component_five_rational_node),
)
assert all(
    number_field.from_sympy(sp.expand(value)) == number_field.zero
    for value in component_five_gradient_at_rational_node
)
component_five_tangent_discriminant = sp.expand(
    sp.diff(correspondence_five, x, y).subs(component_five_rational_node) ** 2
    - sp.diff(correspondence_five, x, 2).subs(component_five_rational_node)
    * sp.diff(correspondence_five, y, 2).subs(component_five_rational_node)
)
assert number_field.from_sympy(component_five_tangent_discriminant) != number_field.zero
assert sp.gcd(2, 5) == 1

# Two quadratic Cremona transformations reduce the five-nodal plane quintic
# to the sparse cubic below.  The first transformation is centered at the
# rational node and the degree-two p2(x),p1^-(y) node orbit.  The second is
# centered at the remaining degree-two node orbit and the K-rational point
# obtained by contracting the line through the first conjugate pair.  The
# optional Singular replay verifies both substitutions exactly.
#
# The displayed cubic is over K(t), t^2=2(11+r).  Its coefficients involve
# only this one of the two temporary quadratic node fields.
cremona_r = sp.sqrt(-11)
cremona_t = sp.sqrt(2 * (11 + cremona_r))
cremona_field = sp.QQ.algebraic_field(cremona_r, cremona_t)
cremona_u, cremona_v, cremona_w = sp.symbols(
    "cremona_u cremona_v cremona_w"
)
cremona_coefficients = {
    (2, 1, 0): (
        sp.Rational(264627, 8)
        - sp.Rational(29403, 8) * cremona_t
        + sp.Rational(8019, 8) * cremona_r
        + sp.Rational(2673, 4) * cremona_r * cremona_t
    ),
    (2, 0, 1): (
        sp.Rational(617463, 16)
        - sp.Rational(65043, 32) * cremona_r * cremona_t
        - sp.Rational(421443, 32) * cremona_t
        + sp.Rational(376893, 16) * cremona_r
    ),
    (1, 2, 0): (
        sp.Rational(264627, 8)
        - sp.Rational(2673, 4) * cremona_r * cremona_t
        + sp.Rational(29403, 8) * cremona_t
        + sp.Rational(8019, 8) * cremona_r
    ),
    (1, 1, 1): 78408 + 60588 * cremona_r,
    (1, 0, 2): (
        -sp.Rational(16335, 2)
        - sp.Rational(5049, 4) * cremona_r * cremona_t
        - sp.Rational(35937, 4) * cremona_t
        + sp.Rational(98307, 2) * cremona_r
    ),
    (0, 2, 1): (
        sp.Rational(617463, 16)
        + sp.Rational(421443, 32) * cremona_t
        + sp.Rational(65043, 32) * cremona_r * cremona_t
        + sp.Rational(376893, 16) * cremona_r
    ),
    (0, 1, 2): (
        -sp.Rational(16335, 2)
        + sp.Rational(35937, 4) * cremona_t
        + sp.Rational(5049, 4) * cremona_r * cremona_t
        + sp.Rational(98307, 2) * cremona_r
    ),
    (0, 0, 3): -28314 + 16434 * cremona_r,
}
component_five_cremona_cubic = sp.Poly(
    sum(
        coefficient
        * cremona_u**monomial[0]
        * cremona_v**monomial[1]
        * cremona_w**monomial[2]
        for monomial, coefficient in cremona_coefficients.items()
    ),
    cremona_u,
    cremona_v,
    cremona_w,
    domain=cremona_field,
)
assert component_five_cremona_cubic.total_degree() == 3
assert len(component_five_cremona_cubic.terms()) == 8


def sparse_ternary_cubic_invariants(cubic):
    """Return c4,c6,Delta for the sparse cubic used above.

    This is the Rodriguez-Villegas ternary-cubic-to-Weierstrass formula,
    specialized to cubics with no U^3 or V^3 term.  It is the same formula
    implemented in Singular's tropical.lib.
    """

    affine_cubic = sp.Poly(
        cubic.as_expr().subs(cremona_w, 1),
        cremona_u,
        cremona_v,
        domain=cremona_field,
    )

    def coefficient(u_degree, v_degree):
        return cremona_field.from_sympy(
            affine_cubic.coeff_monomial(
                cremona_u**u_degree * cremona_v**v_degree
            )
        )

    s1 = coefficient(0, 2)
    s0 = coefficient(1, 2)
    r2 = coefficient(0, 1)
    r1 = coefficient(1, 1)
    r0 = coefficient(2, 1)
    q3 = coefficient(0, 0)
    q2 = coefficient(1, 0)
    q1 = coefficient(2, 0)
    assert coefficient(0, 3) == cremona_field.zero
    assert coefficient(3, 0) == cremona_field.zero

    weierstrass_a1 = r1
    weierstrass_a2 = -(s0 * q2 + s1 * q1 + r0 * r2)
    weierstrass_a3 = -(
        s0 * r0 * q3 + s1 * r0 * q2 + s0 * r2 * q1
    )
    weierstrass_a4 = (
        (s0**2 * q1 + s1 * r0**2) * q3
        + (s1 * s0 * q1 + s0 * r0 * r2) * q2
        + s1 * r0 * r2 * q1
    )
    weierstrass_a6 = (
        (
            -s1 * s0 * r0**2 * q2
            - s1 * s0**2 * q1**2
            + (
                -s0**2 * r0 * r2
                + s1 * s0 * r0 * r1
                - s1**2 * r0**2
            )
            * q1
        )
        * q3
        - s1 * s0 * r0 * r2 * q1 * q2
    )
    b2 = weierstrass_a1**2 + 4 * weierstrass_a2
    b4 = 2 * weierstrass_a4 + weierstrass_a1 * weierstrass_a3
    b6 = weierstrass_a3**2 + 4 * weierstrass_a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    discriminant = (c4**3 - c6**2) / 1728
    return c4, c6, discriminant


(
    component_five_c4,
    component_five_c6,
    component_five_discriminant,
) = sparse_ternary_cubic_invariants(component_five_cremona_cubic)
component_five_j = (
    component_five_c4**3 / component_five_discriminant
)
assert component_five_j == cremona_field.convert(-121)

# Compare with E_121: v^2+uv=u^3+u^2-2u-7, whose invariants are
# c4=121, c6=5203, Delta=-11^4.  The temporary cubic is its twist by
# 2(11+r), up to the displayed square; descent to K is resolved below by
# the two reductions above 23.
component_five_twist_parameter = (
    121 * component_five_c6 / (5203 * component_five_c4)
)
expected_component_five_twist_parameter = cremona_field.from_sympy(
    (22 + 2 * cremona_r) * (4455 + 891 * cremona_r) ** 2
)
assert (
    component_five_twist_parameter
    == expected_component_five_twist_parameter
)


# ---------------------------------------------------------------------------
# The natural degree-twelve quotient as the classical modular curve X_0(11).

elliptic_x, elliptic_y = sp.symbols("elliptic_x elliptic_y")
elliptic_cubic = elliptic_x**3 - elliptic_x**2 - 10 * elliptic_x - 20
elliptic_equation = elliptic_y**2 + elliptic_y - elliptic_cubic

elliptic_basis_12 = (
    1,
    elliptic_x,
    elliptic_y,
    elliptic_x**2,
    elliptic_x * elliptic_y,
    elliptic_x**3,
    elliptic_x**2 * elliptic_y,
    elliptic_x**4,
    elliptic_x**3 * elliptic_y,
    elliptic_x**5,
    elliptic_x**4 * elliptic_y,
    elliptic_x**6,
)
elliptic_basis_13 = elliptic_basis_12 + (elliptic_x**5 * elliptic_y,)
elliptic_denominator_coefficients = (
    8278096,
    -9028165,
    -3062640,
    799092,
    977047,
    288374,
    15778,
    -11099,
    -3511,
    -482,
    -34,
    -1,
)
elliptic_numerator_coefficients = (
    122758525012,
    82027395739,
    5600100593,
    4120902300,
    19113285610,
    -12574009616,
    6795826151,
    -1879636243,
    262696343,
    -15351268,
    172103,
    -709,
    1,
)
elliptic_denominator = sp.expand(
    sum(
        coefficient * basis_element
        for coefficient, basis_element in zip(
            elliptic_denominator_coefficients, elliptic_basis_12, strict=True
        )
    )
)
elliptic_numerator = sp.expand(
    sum(
        coefficient * basis_element
        for coefficient, basis_element in zip(
            elliptic_numerator_coefficients, elliptic_basis_13, strict=True
        )
    )
)


def elliptic_function_norm(expression):
    """Norm under the involution y -> -1-y on the X_0(11) model."""

    polynomial = sp.Poly(expression, elliptic_y, domain=sp.QQ[elliptic_x])
    assert polynomial.degree() <= 1
    constant = polynomial.coeff_monomial(1)
    linear = polynomial.coeff_monomial(elliptic_y)
    return sp.Poly(
        sp.expand(constant * (constant - linear) - linear**2 * elliptic_cubic),
        elliptic_x,
        domain=sp.QQ,
    )


elliptic_zero_quartic = sp.Poly(
    elliptic_x**4
    - 52820 * elliptic_x**3
    + 1333262 * elliptic_x**2
    + 4971236 * elliptic_x
    + 9789217,
    elliptic_x,
    domain=sp.QQ,
)
elliptic_one_sextic = sp.Poly(
    elliptic_x**6
    - 288318 * elliptic_x**5
    + 141521931 * elliptic_x**4
    + 169928888 * elliptic_x**3
    + 8135691435 * elliptic_x**2
    + 30544230678 * elliptic_x
    + 28453700753,
    elliptic_x,
    domain=sp.QQ,
)
assert elliptic_function_norm(elliptic_denominator) == sp.Poly(
    (elliptic_x - 16) ** 12, elliptic_x, domain=sp.QQ
)
assert elliptic_function_norm(elliptic_numerator) == sp.Poly(
    -(elliptic_x - 16) * elliptic_zero_quartic.as_expr() ** 3,
    elliptic_x,
    domain=sp.QQ,
)
assert elliptic_function_norm(
    elliptic_numerator - 1728 * elliptic_denominator
) == sp.Poly(
    -(elliptic_x - 16) * elliptic_one_sextic.as_expr() ** 2,
    elliptic_x,
    domain=sp.QQ,
)
assert sp.gcd(elliptic_zero_quartic, elliptic_zero_quartic.diff()).degree() == 0
assert sp.gcd(elliptic_one_sextic, elliptic_one_sextic.diff()).degree() == 0
assert sp.gcd(elliptic_zero_quartic, elliptic_one_sextic).degree() == 0

elliptic_cusp_zero = {elliptic_x: 16, elliptic_y: -61}
elliptic_common_zero = {elliptic_x: 16, elliptic_y: 60}
assert elliptic_equation.subs(elliptic_cusp_zero) == 0
assert elliptic_equation.subs(elliptic_common_zero) == 0
assert elliptic_denominator.subs(elliptic_cusp_zero) == 0
assert elliptic_denominator.subs(elliptic_common_zero) == 0
elliptic_denominator_quotient, elliptic_denominator_remainder = sp.div(
    sp.Poly(elliptic_denominator, elliptic_x, elliptic_y, domain=sp.QQ),
    sp.Poly(elliptic_x - 16, elliptic_x, elliptic_y, domain=sp.QQ),
)
assert elliptic_denominator_remainder.is_zero
assert elliptic_denominator_quotient.as_expr().subs(elliptic_common_zero) != 0
assert elliptic_numerator.subs(elliptic_cusp_zero) != 0
assert elliptic_numerator.subs(elliptic_common_zero) == 0
assert (
    elliptic_numerator - 1728 * elliptic_denominator
).subs(elliptic_cusp_zero) != 0
assert (
    elliptic_numerator - 1728 * elliptic_denominator
).subs(elliptic_common_zero) == 0

# The reduced boundary divisors above j=0 and j=1728 are K-prime: their
# x-coordinate polynomials remain irreducible over K=Q(sqrt(-11)).  Denote
# them by Z (degree four) and W (degree six).  The two poles are
# P=(16,-61), of order eleven, and O, of order one.
for boundary_polynomial, boundary_degree in (
    (elliptic_zero_quartic, 4),
    (elliptic_one_sextic, 6),
):
    _, boundary_factorization = sp.factor_list(
        boundary_polynomial.as_expr(),
        elliptic_x,
        extension=r,
    )
    assert len(boundary_factorization) == 1
    irreducible_boundary_factor, boundary_exponent = boundary_factorization[0]
    assert boundary_exponent == 1
    assert sp.degree(irreducible_boundary_factor, elliptic_x) == boundary_degree

# P has order five.  The tangent at P meets E again at 3P=(5,5), while the
# tangent at 2P=(5,-6) meets E again at P.  The resulting Miller function
# has divisor 5P-5O.
elliptic_tangent_at_p = elliptic_y + 6 * elliptic_x - 35
elliptic_tangent_at_2p = elliptic_y + 5 * elliptic_x - 19
assert sp.expand(
    elliptic_equation.subs(elliptic_y, -6 * elliptic_x + 35)
    + (elliptic_x - 16) ** 2 * (elliptic_x - 5)
) == 0
assert sp.expand(
    elliptic_equation.subs(elliptic_y, -5 * elliptic_x + 19)
    + (elliptic_x - 16) * (elliptic_x - 5) ** 2
) == 0
elliptic_torsion_unit = (
    elliptic_tangent_at_p**2
    * elliptic_tangent_at_2p
    / (elliptic_x - 5) ** 2
)

# Primitive units for the reduced zero and one fibers.  Their divisors are
# Z-2P-2O and W-3P-3O.  The identities below show directly that their cubes
# and squares recover j and j-1728 after the torsion correction.
elliptic_zero_unit_numerator = (
    -elliptic_x**3
    - 2856 * elliptic_x**2
    + 242 * elliptic_x * elliptic_y
    - 10206 * elliptic_x
    + 10769 * elliptic_y
    + 20068
)
elliptic_one_unit_numerator = (
    501 * elliptic_x**4
    + elliptic_x**3 * elliptic_y
    + 191605 * elliptic_x**3
    - 18682 * elliptic_x**2 * elliptic_y
    + 1619779 * elliptic_x**2
    - 852403 * elliptic_x * elliptic_y
    - 3303386 * elliptic_x
    - 1070227 * elliptic_y
    - 10226604
)
elliptic_zero_unit = elliptic_zero_unit_numerator / (elliptic_x - 16) ** 2
elliptic_one_unit = elliptic_one_unit_numerator / (elliptic_x - 16) ** 3


def elliptic_function_is_zero(expression):
    """Test a rational-function identity on the X_0(11) model."""

    numerator = sp.together(expression).as_numer_denom()[0]
    rational_x_field = sp.QQ.frac_field(elliptic_x)
    remainder = sp.rem(
        sp.Poly(numerator, elliptic_y, domain=rational_x_field),
        sp.Poly(elliptic_equation, elliptic_y, domain=rational_x_field),
    )
    return remainder.is_zero


elliptic_j_function = elliptic_numerator / elliptic_denominator
assert elliptic_function_is_zero(
    elliptic_zero_unit**3 - elliptic_j_function * elliptic_torsion_unit
)
assert elliptic_function_is_zero(
    elliptic_one_unit**2
    + (elliptic_j_function - 1728) * elliptic_torsion_unit
)

# In the ordered K-prime boundary basis (Z,W,P,O), the full unit lattice is
# generated by the three primitive divisor rows below.  The quotient has
# Smith form Z plus Z/5: O supplies degree one and P-O supplies the torsion
# generator, so these rows are the exact kernel of Div_S(E)->Pic(E).
elliptic_boundary_unit_lattice = sp.Matrix(
    (
        (1, 0, -2, -2),
        (0, 1, -3, -3),
        (0, 0, 5, -5),
    )
)
elliptic_boundary_degrees = sp.Matrix((4, 6, 1, 1))
elliptic_boundary_torsion_classes = sp.Matrix((2, 3, 1, 0))
for lattice_row in elliptic_boundary_unit_lattice.tolist():
    row = sp.Matrix(1, 4, lattice_row)
    assert (row * elliptic_boundary_degrees)[0] == 0
    assert (row * elliptic_boundary_torsion_classes)[0] % 5 == 0
assert smith_normal_form(
    elliptic_boundary_unit_lattice,
    domain=ZZ,
) == sp.Matrix(((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 5, 0)))

# The long Weierstrass coefficients [0,-1,1,-10,-20] have discriminant
# -11^5.  The norm identities prove the passport
# 3^4 | 2^6 | 11 1 for j = numerator/denominator after cancellation of
# the common simple zero (16,60).
a1, a2, a3, a4, a6 = 0, -1, 1, -10, -20
b2 = a1**2 + 4 * a2
b4 = 2 * a4 + a1 * a3
b6 = a3**2 + 4 * a6
b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
elliptic_discriminant = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
assert elliptic_discriminant == -(11**5)


def number_field_coefficient_mod_prime(coefficient, square_root, prime):
    """Reduce a coefficient in Q(sqrt(-11)) at a split prime."""

    value = 0
    for rational in number_field.from_sympy(coefficient).to_sympy_list():
        value = (
            value * square_root
            + int(rational.p) * pow(int(rational.q), -1, prime)
        ) % prime
    return value


def polynomial_terms_mod_prime(polynomial, square_root, prime):
    """Return sparse bivariate terms after r -> square_root modulo prime."""

    return [
        (
            monomial,
            number_field_coefficient_mod_prime(
                coefficient, square_root, prime
            ),
        )
        for monomial, coefficient in polynomial.terms()
    ]


def evaluate_bivariate_terms(terms, x_value, y_value, prime):
    """Evaluate sparse bivariate terms over a prime field."""

    return sum(
        coefficient
        * pow(x_value, monomial[0], prime)
        * pow(y_value, monomial[1], prime)
        for monomial, coefficient in terms
    ) % prime


def ideal_intersection_mod_prime(first, second, prime):
    """Intersect two plane ideals by one-variable elimination."""

    auxiliary = sp.symbols("ideal_intersection_auxiliary")
    elimination_basis = sp.groebner(
        [auxiliary * generator for generator in first]
        + [(1 - auxiliary) * generator for generator in second],
        auxiliary,
        x,
        y,
        order="lex",
        modulus=prime,
    )
    return [
        polynomial.as_expr()
        for polynomial in elimination_basis.polys
        if not polynomial.as_expr().has(auxiliary)
    ]


def count_weierstrass_points_mod_prime(coefficients, prime):
    """Count a generalized Weierstrass model over F_p, including infinity."""

    curve_a1, curve_a2, curve_a3, curve_a4, curve_a6 = coefficients
    return 1 + sum(
        (
            y_value**2
            + curve_a1 * x_value * y_value
            + curve_a3 * y_value
            - x_value**3
            - curve_a2 * x_value**2
            - curve_a4 * x_value
            - curve_a6
        )
        % prime
        == 0
        for x_value in range(prime)
        for y_value in range(prime)
    )


# Exact good-reduction descent at the two primes of K above 23.  The nodal
# plane model has 20 affine and five infinite F_23-points.  Its three
# rational nodes have nonsplit tangents, so normalization removes one point
# at each node and has 22 points, trace 2.  The singular ideal remains the
# disjoint reduced length-five node support and the Hessian is a unit on it;
# hence no extra geometric singularity appears in this reduction.
descent_prime = 23
component_five_modular_counts = []
component_five_modular_traces = []
component_five_integral_polynomial = sp.Poly(
    2 * correspondence_five, x, y, domain=number_field
)
component_five_integral_top_form = sp.Poly(
    sum(
        coefficient * x**monomial[0] * y**monomial[1]
        for monomial, coefficient in component_five_integral_polynomial.terms()
        if sum(monomial) == 5
    ),
    x,
    y,
    domain=number_field,
)
for modular_square_root in (9, 14):
    assert (modular_square_root**2 + 11) % descent_prime == 0

    reduced_polynomials = []
    for polynomial in (
        component_five_integral_polynomial,
        component_five_integral_polynomial.diff(x),
        component_five_integral_polynomial.diff(y),
        component_five_integral_polynomial.diff(x).diff(x),
        component_five_integral_polynomial.diff(x).diff(y),
        component_five_integral_polynomial.diff(y).diff(y),
    ):
        reduced_polynomials.append(
            polynomial_terms_mod_prime(
                polynomial, modular_square_root, descent_prime
            )
        )
    (
        reduced_curve,
        reduced_x_derivative,
        reduced_y_derivative,
        reduced_xx_derivative,
        reduced_xy_derivative,
        reduced_yy_derivative,
    ) = reduced_polynomials
    reduced_top_form = polynomial_terms_mod_prime(
        component_five_integral_top_form,
        modular_square_root,
        descent_prime,
    )

    affine_points = [
        (x_value, y_value)
        for x_value in range(descent_prime)
        for y_value in range(descent_prime)
        if evaluate_bivariate_terms(
            reduced_curve, x_value, y_value, descent_prime
        )
        == 0
    ]
    infinite_points = [
        (1, y_value)
        for y_value in range(descent_prime)
        if evaluate_bivariate_terms(
            reduced_top_form, 1, y_value, descent_prime
        )
        == 0
    ]
    if (
        evaluate_bivariate_terms(reduced_top_form, 0, 1, descent_prime)
        == 0
    ):
        infinite_points.append((0, 1))
    rational_nodes = [
        point
        for point in affine_points
        if evaluate_bivariate_terms(
            reduced_x_derivative, *point, descent_prime
        )
        == 0
        and evaluate_bivariate_terms(
            reduced_y_derivative, *point, descent_prime
        )
        == 0
    ]
    tangent_branch_counts = []
    inverse_two = pow(2, -1, descent_prime)
    for node in rational_nodes:
        tangent_xx = (
            inverse_two
            * evaluate_bivariate_terms(
                reduced_xx_derivative, *node, descent_prime
            )
        ) % descent_prime
        tangent_xy = evaluate_bivariate_terms(
            reduced_xy_derivative, *node, descent_prime
        )
        tangent_yy = (
            inverse_two
            * evaluate_bivariate_terms(
                reduced_yy_derivative, *node, descent_prime
            )
        ) % descent_prime
        tangent_branch_counts.append(
            sum(
                (
                    tangent_xx
                    + tangent_xy * slope
                    + tangent_yy * slope**2
                )
                % descent_prime
                == 0
                for slope in range(descent_prime)
            )
            + int(tangent_yy == 0)
        )

    # Check the full geometric singular scheme against the exact support
    # V(p2(x),p1^-(y)) union V(p1(x),p2^-(y)) union V(q1(x),q1^-(y)).
    reduced_expression = sp.Poly(
        sum(
            coefficient * x**monomial[0] * y**monomial[1]
            for monomial, coefficient in reduced_curve
        ),
        x,
        y,
        modulus=descent_prime,
    ).as_expr()
    reduced_p1_x = 2 * x + 11 - 3 * modular_square_root
    reduced_p2_x = (
        2 * x**2
        - (11 - 3 * modular_square_root) * x
        - (22 + 6 * modular_square_root)
    )
    reduced_q1_x = 2 * x + 5 + 3 * modular_square_root
    reduced_p1_minus_y = 2 * y + 11 + 3 * modular_square_root
    reduced_p2_minus_y = (
        2 * y**2
        - (11 + 3 * modular_square_root) * y
        - (22 - 6 * modular_square_root)
    )
    reduced_q1_minus_y = 2 * y + 5 - 3 * modular_square_root
    node_support = ideal_intersection_mod_prime(
        ideal_intersection_mod_prime(
            [reduced_p2_x, reduced_p1_minus_y],
            [reduced_p1_x, reduced_p2_minus_y],
            descent_prime,
        ),
        [reduced_q1_x, reduced_q1_minus_y],
        descent_prime,
    )
    singular_basis = sp.groebner(
        [
            reduced_expression,
            sp.diff(reduced_expression, x),
            sp.diff(reduced_expression, y),
        ],
        x,
        y,
        order="grevlex",
        modulus=descent_prime,
    )
    support_basis = sp.groebner(
        node_support,
        x,
        y,
        order="grevlex",
        modulus=descent_prime,
    )
    assert all(
        singular_basis.reduce(generator)[1] == 0
        for generator in node_support
    )
    assert all(
        support_basis.reduce(generator.as_expr())[1] == 0
        for generator in singular_basis.polys
    )
    leading_monomials = [
        polynomial.LM(order=singular_basis.order).exponents
        for polynomial in singular_basis.polys
    ]
    standard_monomials = [
        (x_degree, y_degree)
        for x_degree in range(5)
        for y_degree in range(5)
        if not any(
            x_degree >= leading[0] and y_degree >= leading[1]
            for leading in leading_monomials
        )
    ]
    assert len(standard_monomials) == 5
    reduced_hessian = (
        sp.diff(reduced_expression, x, 2)
        * sp.diff(reduced_expression, y, 2)
        - sp.diff(reduced_expression, x, y) ** 2
    )
    hessian_basis = sp.groebner(
        node_support + [reduced_hessian],
        x,
        y,
        order="grevlex",
        modulus=descent_prime,
    )
    assert len(hessian_basis.polys) == 1
    assert hessian_basis.polys[0].as_expr() == 1

    infinity_polynomial_mod_prime = sp.Poly(
        sum(
            coefficient * x**monomial[0] * y**monomial[1]
            for monomial, coefficient in reduced_top_form
        ).subs(y, 1),
        x,
        modulus=descent_prime,
    )
    assert sp.gcd(
        infinity_polynomial_mod_prime,
        infinity_polynomial_mod_prime.diff(),
    ).degree() == 0
    assert len(affine_points) == 20
    assert len(infinite_points) == 5
    assert len(rational_nodes) == 3
    assert tangent_branch_counts == [0, 0, 0]
    normalized_point_count = (
        len(affine_points)
        + len(infinite_points)
        + sum(branches - 1 for branches in tangent_branch_counts)
    )
    component_five_modular_counts.append(normalized_point_count)
    component_five_modular_traces.append(
        descent_prime + 1 - normalized_point_count
    )

assert component_five_modular_counts == [22, 22]
assert component_five_modular_traces == [2, 2]

# The two possible temporary node extensions have square classes
# e_s=2(11-r), e_t=2(11+r).  Their residue characters at the two primes
# above 23 distinguish all four descent twists.  The observed trace pair
# (2,2) therefore selects the untwisted E_121 model over K.
twist_character_pairs = []
for twist_rational_part, twist_r_coefficient in (
    (1, 0),
    (22, -2),
    (22, 2),
    (528, 0),
):
    character_pair = []
    for modular_square_root in (9, 14):
        residue = (
            twist_rational_part
            + twist_r_coefficient * modular_square_root
        ) % descent_prime
        legendre_value = pow(residue, (descent_prime - 1) // 2, descent_prime)
        character_pair.append(1 if legendre_value == 1 else -1)
    twist_character_pairs.append(tuple(character_pair))
assert twist_character_pairs == [(1, 1), (1, -1), (-1, 1), (-1, -1)]

component_five_weierstrass_coefficients = (1, 1, 0, -2, -7)
component_five_weierstrass_count = count_weierstrass_points_mod_prime(
    component_five_weierstrass_coefficients, descent_prime
)
natural_weierstrass_count = count_weierstrass_points_mod_prime(
    (0, -1, 1, -10, -20), descent_prime
)
assert component_five_weierstrass_count == 22
assert natural_weierstrass_count == 25
assert descent_prime + 1 - component_five_weierstrass_count == 2
assert descent_prime + 1 - natural_weierstrass_count == -1

# Invariants of the selected small model, and its rational -11 twist
# 121a2.  The twist becomes isomorphic over K because r^2=-11.
(
    component_model_a1,
    component_model_a2,
    component_model_a3,
    component_model_a4,
    component_model_a6,
) = component_five_weierstrass_coefficients
component_model_b2 = component_model_a1**2 + 4 * component_model_a2
component_model_b4 = 2 * component_model_a4 + component_model_a1 * component_model_a3
component_model_b6 = component_model_a3**2 + 4 * component_model_a6
component_model_b8 = (
    component_model_a1**2 * component_model_a6
    + 4 * component_model_a2 * component_model_a6
    - component_model_a1 * component_model_a3 * component_model_a4
    + component_model_a2 * component_model_a3**2
    - component_model_a4**2
)
component_five_model_c4 = component_model_b2**2 - 24 * component_model_b4
component_five_model_c6 = (
    -component_model_b2**3
    + 36 * component_model_b2 * component_model_b4
    - 216 * component_model_b6
)
component_five_model_discriminant = (
    -component_model_b2**2 * component_model_b8
    - 8 * component_model_b4**3
    - 27 * component_model_b6**2
    + 9 * component_model_b2 * component_model_b4 * component_model_b6
)
assert (
    component_five_model_c4,
    component_five_model_c6,
    component_five_model_discriminant,
) == (121, 5203, -(11**4))
assert component_five_model_c4**3 // component_five_model_discriminant == -121
assert 14641 == (-11) ** 2 * component_five_model_c4
assert -6925193 == (-11) ** 3 * component_five_model_c6

# The K-prime punctures of the normalization above the seven prime boundary
# divisors p1,p2,p3,q1,q2,q3,infinity have the following degrees.  The
# optional Singular replay proves the affine decomposition and separates the
# two branches at each nodal boundary point; irreducibility of the top form
# gives the final degree-five prime at infinity.
component_five_top_form = sum(
    coefficient * x**monomial[0] * y**monomial[1]
    for monomial, coefficient in sp.Poly(
        correspondence_five, x, y, domain=number_field
    ).terms()
    if sum(monomial) == 5
)
_, component_five_infinity_factorization = sp.factor_list(
    component_five_top_form.subs(y, 1),
    x,
    extension=r,
)
assert len(component_five_infinity_factorization) == 1
assert component_five_infinity_factorization[0][1] == 1
assert sp.degree(component_five_infinity_factorization[0][0], x) == 5
component_five_boundary_blocks = {
    "p1": (1, 4),
    "p2": (4, 4, 2),
    "p3": (4, 2),
    "q1": (2, 3),
    "q2": (3, 6, 6),
    "q3": (6, 3),
    "infinity": (5,),
}
component_five_boundary_labels = (
    "p1_p1",
    "p1_p2_node",
    "p2_p1_node",
    "p2_p2",
    "p2_p3",
    "p3_p3",
    "p3_p2",
    "q1_q1_node",
    "q1_q2",
    "q2_q1",
    "q2_q2",
    "q2_q3",
    "q3_q2",
    "q3_q3",
    "infinity_5",
)
component_five_boundary_degrees = tuple(
    degree
    for block in component_five_boundary_blocks.values()
    for degree in block
)
assert len(component_five_boundary_labels) == 15
assert len(component_five_boundary_degrees) == 15
assert sum(component_five_boundary_degrees) == 55
assert sum(
    (component_five_boundary_blocks[label] for label in ("p1", "p2", "p3")),
    (),
) == (
    1,
    4,
    4,
    4,
    2,
    4,
    2,
)
assert sum(
    (component_five_boundary_blocks[label] for label in ("q1", "q2", "q3")),
    (),
) == (
    2,
    3,
    3,
    6,
    6,
    6,
    3,
)

# The elliptic model has trivial K-torsion: at the split good primes 3 and
# 5 its reductions have coprime orders 2 and 5.  The optional PARI replay
# proves that this curve and its -11 twist both have rational rank zero, so
# E_121(K)=0.  Consequently every K-rational degree-zero boundary divisor is
# principal, and the full unit lattice is simply the degree kernel.
assert count_weierstrass_points_mod_prime((0, -1, 1, -10, -20), 3) == 5
assert count_weierstrass_points_mod_prime((0, -1, 1, -10, -20), 5) == 5
assert count_weierstrass_points_mod_prime(
    component_five_weierstrass_coefficients, 3
) == 2
assert count_weierstrass_points_mod_prime(
    component_five_weierstrass_coefficients, 5
) == 5
component_five_boundary_degree_matrix = sp.Matrix(
    (component_five_boundary_degrees,)
)
component_five_boundary_unit_lattice = sp.zeros(
    len(component_five_boundary_degrees) - 1,
    len(component_five_boundary_degrees),
)
for row_index, boundary_degree in enumerate(
    component_five_boundary_degrees[1:]
):
    component_five_boundary_unit_lattice[row_index, 0] = -boundary_degree
    component_five_boundary_unit_lattice[row_index, row_index + 1] = 1
assert (
    component_five_boundary_unit_lattice
    * component_five_boundary_degree_matrix.T
) == sp.zeros(14, 1)
assert smith_normal_form(
    component_five_boundary_unit_lattice,
    domain=ZZ,
) == sp.Matrix.hstack(sp.eye(14), sp.zeros(14, 1))


# ---------------------------------------------------------------------------
# The genus-two normalization and its two elliptic quotients.

# Adjunction realizes the canonical system of an eight-nodal plane sextic
# as the cubics through its nodes.  These two cubics are a reduced basis of
# that pencil.  The three ideals below are the disjoint degree 4+1+3 node
# orbits from the exact conductor calculation.
component_six_adjoint_a = (
    4 * x**2 * y
    + (r - 1) * x * y**2
    - 2 * y**3
    + 4 * x**2
    + 2 * r * x * y
    + (r - 7) * y**2
    + (-8 * r - 44) * x
    + (20 * r + 44) * y
    - 44 * r
    + 220
)
component_six_adjoint_b = (
    x**3
    + 3 * x**2
    - sp.Rational(1, 2) * x * y**2
    + sp.Rational(1, 4) * (11 + 3 * r) * x * y
    + sp.Rational(1, 2) * (-55 + 33 * r) * x
    + sp.Rational(1, 4) * (1 + r) * y**3
    + sp.Rational(1, 4) * (1 + 3 * r) * y**2
    + (44 - 12 * r) * y
    - 198
    - 44 * r
)
component_six_node_ideals = (
    (
        p2,
        2 * y**2 - (11 + 3 * r) * y - (22 - 6 * r),
    ),
    (p1, 2 * y + 11 + 3 * r),
    (
        2 * y**3
        + (15 + 3 * r) * y**2
        - (12 + 12 * r) * y
        + 56
        - 96 * r,
        18 * x + (r + 1) * y**2 + (8 * r - 10) * y + 88 - 20 * r,
    ),
)
component_six_node_bases = tuple(
    sp.groebner(ideal, x, y, order="lex", domain=number_field)
    for ideal in component_six_node_ideals
)
for node_basis in component_six_node_bases:
    assert node_basis.reduce(component_six_adjoint_a)[1] == 0
    assert node_basis.reduce(component_six_adjoint_b)[1] == 0

# The restriction map from the ten plane cubics to the reduced length-eight
# node scheme has rank eight.  Hence its kernel is exactly the canonical
# pencil spanned by the two independent adjoints above.
plane_cubic_monomials = tuple(
    x**x_degree * y**y_degree
    for total_degree in range(4)
    for x_degree in range(total_degree + 1)
    for y_degree in (total_degree - x_degree,)
)
component_six_node_remainders = tuple(
    tuple(
        sp.Poly(
            node_basis.reduce(monomial)[1],
            x,
            y,
            domain=number_field,
        )
        for monomial in plane_cubic_monomials
    )
    for node_basis in component_six_node_bases
)
component_six_restriction_rows = sorted(
    {
        (orbit_index, monomial)
        for orbit_index, remainders in enumerate(component_six_node_remainders)
        for remainder in remainders
        for monomial, _ in remainder.terms()
    }
)
component_six_restriction_matrix = sp.Matrix(
    [
        [
            component_six_node_remainders[orbit_index][column].coeff_monomial(
                x**monomial[0] * y**monomial[1]
            )
            for column in range(len(plane_cubic_monomials))
        ]
        for orbit_index, monomial in component_six_restriction_rows
    ]
)
assert component_six_restriction_matrix.shape == (8, 10)
assert component_six_restriction_matrix.rank() == 8

# Put canonical_t=B/A.  Singular eliminates y from C6 and B-tA, removes
# the node factors, and verifies that the remaining quadratic discriminant
# has the following sextic square class.  The optional replay below proves
# that elimination identity; the portable half checks all subsequent exact
# changes of variable.
canonical_t, hyperelliptic_z = sp.symbols("canonical_t hyperelliptic_z")
component_six_canonical_sextic = (
    canonical_t**6
    + (sp.Rational(3, 4) * r - sp.Rational(69, 52)) * canonical_t**5
    + (-sp.Rational(345, 416) * r - sp.Rational(659, 416))
    * canonical_t**4
    + (sp.Rational(7, 104) * r + sp.Rational(357, 208))
    * canonical_t**3
    + (sp.Rational(489, 6656) * r - sp.Rational(2675, 6656))
    * canonical_t**2
    + (-sp.Rational(243, 13312) * r + sp.Rational(615, 13312))
    * canonical_t
    + sp.Rational(15, 13312) * r
    - sp.Rational(71, 53248)
)
component_six_discriminant_constant = 7907328 * r + 2635776
component_six_canonical_square_factor = (
    canonical_t**3
    + (sp.Rational(1, 40) * r - sp.Rational(7, 40)) * canonical_t**2
    + (sp.Rational(7, 160) * r - sp.Rational(19, 160)) * canonical_t
    - sp.Rational(1, 320) * r
    + sp.Rational(1, 160)
)

# Clearing the canonical coordinate by u=4*canonical_t gives a small
# integral sextic.  The two discriminant representatives differ by the
# square (3r/2)^2 in K.
component_six_hyperelliptic_sextic = (
    26 * u**6
    + (78 * r - 138) * u**5
    + (-345 * r - 659) * u**4
    + (112 * r + 2856) * u**3
    + (489 * r - 2675) * u**2
    + (-486 * r + 1230) * u
    + 120 * r
    - 142
)
component_six_hyperelliptic_polynomial = (
    -(3 * r + 1) * component_six_hyperelliptic_sextic
)
assert polynomial_is_zero(
    component_six_discriminant_constant * component_six_canonical_sextic
    - (3 * r / 2) ** 2
    * component_six_hyperelliptic_polynomial.subs(u, 4 * canonical_t),
    canonical_t,
)

# The canonical sextic has the K-rational involution
# u |-> (-6u+1-7r)/((1-r)u+6).  Moving its two fixed points to zero and
# infinity sends the involution to z |-> -z and produces the even model.
component_six_involution = (
    -6 * u + 1 - 7 * r
) / ((1 - r) * u + 6)
component_six_fixed_alpha = (3 - r) / 2
component_six_fixed_beta = (-5 - r) / 2
component_six_z_of_u = (
    (u - component_six_fixed_alpha) / (u - component_six_fixed_beta)
)
component_six_involution_check = sp.together(
    component_six_z_of_u.subs(u, component_six_involution)
    + component_six_z_of_u
).as_numer_denom()[0]
assert polynomial_is_zero(component_six_involution_check, u)
component_six_u_of_z = (
    component_six_fixed_beta * hyperelliptic_z - component_six_fixed_alpha
) / (hyperelliptic_z - 1)
component_six_even_sextic = (
    1
    - 11 * hyperelliptic_z**2
    - 77 * hyperelliptic_z**4
    - 121 * hyperelliptic_z**6
)
assert polynomial_is_zero(
    sp.cancel(
        (hyperelliptic_z - 1) ** 6
        * component_six_hyperelliptic_polynomial.subs(
            u, component_six_u_of_z
        )
    )
    - 512 * (1 + 3 * r) * component_six_even_sextic,
    hyperelliptic_z,
)

# Thus, after scaling the ordinate by (z-1)^3/16, the normalization is
#   Y^2=d(1-11z^2-77z^4-121z^6),  d=2(1+3r).
# Its two non-hyperelliptic quotient maps are
#   (z,Y) -> (X=z^2,Y) and (z,Y) -> (X=z^2,W=zY).
# The first target is a cubic and the second a quartic.  We compute their
# Jacobian invariants and compare twists with small rational models.
component_six_twist_d = 2 * (1 + 3 * r)


def generalized_weierstrass_invariants(coefficients):
    """Return c4,c6,Delta,j for long Weierstrass coefficients."""

    curve_a1, curve_a2, curve_a3, curve_a4, curve_a6 = coefficients
    curve_b2 = curve_a1**2 + 4 * curve_a2
    curve_b4 = 2 * curve_a4 + curve_a1 * curve_a3
    curve_b6 = curve_a3**2 + 4 * curve_a6
    curve_b8 = (
        curve_a1**2 * curve_a6
        + 4 * curve_a2 * curve_a6
        - curve_a1 * curve_a3 * curve_a4
        + curve_a2 * curve_a3**2
        - curve_a4**2
    )
    curve_c4 = curve_b2**2 - 24 * curve_b4
    curve_c6 = -curve_b2**3 + 36 * curve_b2 * curve_b4 - 216 * curve_b6
    curve_discriminant = (
        -curve_b2**2 * curve_b8
        - 8 * curve_b4**3
        - 27 * curve_b6**2
        + 9 * curve_b2 * curve_b4 * curve_b6
    )
    return (
        curve_c4,
        curve_c6,
        curve_discriminant,
        sp.cancel(sp.sympify(curve_c4) ** 3 / curve_discriminant),
    )


# For E_+, x=-121*d*X and y=-121*d*Y give the displayed long twist of
# y^2=x^3-77x^2+1331x+14641.
component_six_plus_coefficients = (
    0,
    -77 * component_six_twist_d,
    0,
    1331 * component_six_twist_d**2,
    14641 * component_six_twist_d**3,
)
component_six_plus_invariants = generalized_weierstrass_invariants(
    component_six_plus_coefficients
)
assert component_six_plus_invariants[3] == -sp.Rational(4096, 11)

# E_11: y^2+y=x^3-x^2 has invariants (16,-152,-11).  E_+ is its twist
# by 44*d, which is already a square in K.
component_six_eleven_model = (0, -1, 1, 0, 0)
component_six_eleven_invariants = generalized_weierstrass_invariants(
    component_six_eleven_model
)
assert component_six_eleven_invariants == (
    16,
    -152,
    -11,
    -sp.Rational(4096, 11),
)
component_six_plus_twist = 44 * component_six_twist_d
assert polynomial_is_zero(component_six_plus_twist - (22 + 6 * r) ** 2, u)
assert polynomial_is_zero(
    component_six_plus_invariants[0]
    - component_six_eleven_invariants[0] * component_six_plus_twist**2,
    u,
)
assert polynomial_is_zero(
    component_six_plus_invariants[1]
    - component_six_eleven_invariants[1] * component_six_plus_twist**3,
    u,
)

# For E_-, the binary quartic is
# d*(-121X^4-77X^3-11X^2+X).  Its classical invariants I,J give the
# Jacobian y^2=x^3-27Ix-27J.  The resulting curve is a square twist of the
# CM conductor-121 model y^2+y=x^3-x^2-7x+10.
component_six_minus_quartic = (-121, -77, -11, 1, 0)
quartic_a, quartic_b, quartic_c, quartic_d, quartic_e = (
    component_six_minus_quartic
)
component_six_minus_i = (
    12 * quartic_a * quartic_e
    - 3 * quartic_b * quartic_d
    + quartic_c**2
)
component_six_minus_j = (
    72 * quartic_a * quartic_c * quartic_e
    + 9 * quartic_b * quartic_c * quartic_d
    - 27 * quartic_a * quartic_d**2
    - 27 * quartic_b**2 * quartic_e
    - 2 * quartic_c**3
)
assert (component_six_minus_i, component_six_minus_j) == (352, 13552)
component_six_minus_coefficients = (
    0,
    0,
    0,
    -27 * component_six_minus_i * component_six_twist_d**2,
    -27 * component_six_minus_j * component_six_twist_d**3,
)
component_six_minus_invariants = generalized_weierstrass_invariants(
    component_six_minus_coefficients
)
assert component_six_minus_invariants[3] == -32768
component_six_cm_model = (0, -1, 1, -7, 10)
component_six_cm_invariants = generalized_weierstrass_invariants(
    component_six_cm_model
)
assert component_six_cm_invariants == (352, -6776, -1331, -32768)
component_six_minus_twist = -36 * component_six_twist_d
assert polynomial_is_zero(
    component_six_minus_twist - (6 * (3 - r)) ** 2, u
)
assert polynomial_is_zero(
    component_six_minus_invariants[0]
    - component_six_cm_invariants[0] * component_six_minus_twist**2,
    u,
)
assert polynomial_is_zero(
    component_six_minus_invariants[1]
    - component_six_cm_invariants[1] * component_six_minus_twist**3,
    u,
)


# ---------------------------------------------------------------------------
# The full boundary-unit lattice of the genus-two normalization.

# The twenty K-prime punctures occur in the following normalized blocks.
# At p1/p1 and q2/q2 the two geometric node branches form one K-prime;
# at p2/p2 they form two degree-four K-primes.  The top form has one
# rational and one degree-five factor.
component_six_top_form = sum(
    coefficient * x**monomial[0] * y**monomial[1]
    for monomial, coefficient in sp.Poly(
        correspondence_six, x, y, domain=number_field
    ).terms()
    if sum(monomial) == 6
)
_, component_six_infinity_factorization = sp.factor_list(
    component_six_top_form.subs(y, 1),
    x,
    extension=r,
)
assert sorted(
    sp.degree(factor, x)
    for factor, exponent in component_six_infinity_factorization
    if exponent == 1
) == [1, 5]

component_six_boundary_labels = (
    "p1_node",
    "p1_a",
    "p1_b",
    "p2_a",
    "p2_node_a",
    "p2_node_b",
    "p2_b",
    "p3_a",
    "p3_b",
    "q1_a",
    "q1_b",
    "q2_a",
    "q2_node",
    "q2_b",
    "q2_c",
    "q3_a",
    "q3_b",
    "q3_c",
    "infinity_1",
    "infinity_5",
)
component_six_boundary_degrees = sp.Matrix(
    (2, 2, 2, 2, 4, 4, 2, 2, 2, 3, 3, 3, 6, 3, 6, 3, 3, 6, 1, 5)
)
assert len(component_six_boundary_labels) == 20
assert sum(component_six_boundary_degrees) == 64

# Use T=(0,0) on E_+: y^2+y=x^3+2x^2+x.  It has order five.  On
# E_-: y^2+y=x^3+2x^2-6x+3 use the independent points
#   G0=(3,5), G2=(-8,(-1+11r)/2).
# Here G0 is the rational Heegner generator and G2 is the image of the
# generator (88,660) on the -11 twist.  The optional PARI replay proves the
# two rational ranks.  The exact boundary-image elimination is replayed by
# verify_psl2_11_c6_boundary_images.py.
component_six_plus_model = (0, 2, 1, 1, 0)
component_six_minus_model = (0, 2, 1, -6, 3)


def weierstrass_point_is_on_curve(point, coefficients):
    """Test a finite point on a long Weierstrass model over K."""

    point_x, point_y = point
    curve_a1, curve_a2, curve_a3, curve_a4, curve_a6 = coefficients
    return polynomial_is_zero(
        point_y**2
        + curve_a1 * point_x * point_y
        + curve_a3 * point_y
        - point_x**3
        - curve_a2 * point_x**2
        - curve_a4 * point_x
        - curve_a6,
        u,
    )


def weierstrass_add(left, right, coefficients):
    """Add two points on a long Weierstrass model; None denotes infinity."""

    if left is None:
        return right
    if right is None:
        return left
    x_left, y_left = left
    x_right, y_right = right
    curve_a1, curve_a2, curve_a3, curve_a4, curve_a6 = coefficients
    if polynomial_is_zero(x_left - x_right, u):
        if polynomial_is_zero(
            y_left + y_right + curve_a1 * x_left + curve_a3,
            u,
        ):
            return None
        denominator = 2 * y_left + curve_a1 * x_left + curve_a3
        slope = (
            3 * x_left**2
            + 2 * curve_a2 * x_left
            + curve_a4
            - curve_a1 * y_left
        ) / denominator
        intercept = (
            -x_left**3
            + curve_a4 * x_left
            + 2 * curve_a6
            - curve_a3 * y_left
        ) / denominator
    else:
        slope = (y_right - y_left) / (x_right - x_left)
        intercept = (
            y_left * x_right - y_right * x_left
        ) / (x_right - x_left)
    result_x = sp.cancel(
        slope**2
        + curve_a1 * slope
        - curve_a2
        - x_left
        - x_right
    )
    result_y = sp.cancel(
        -(slope + curve_a1) * result_x - intercept - curve_a3
    )
    assert weierstrass_point_is_on_curve(
        (result_x, result_y), coefficients
    )
    return result_x, result_y


def weierstrass_multiply(scalar, point, coefficients):
    """Multiply a point by an integer on a long Weierstrass model."""

    if scalar < 0:
        point_x, point_y = point
        curve_a1, _, curve_a3, _, _ = coefficients
        inverse = (
            point_x,
            -point_y - curve_a1 * point_x - curve_a3,
        )
        return weierstrass_multiply(-scalar, inverse, coefficients)
    answer = None
    summand = point
    while scalar:
        if scalar & 1:
            answer = weierstrass_add(answer, summand, coefficients)
        scalar >>= 1
        if scalar:
            summand = weierstrass_add(summand, summand, coefficients)
    return answer


component_six_plus_generator = (sp.Integer(0), sp.Integer(0))
assert weierstrass_multiply(
    5, component_six_plus_generator, component_six_plus_model
) is None
assert all(
    weierstrass_multiply(
        scalar, component_six_plus_generator, component_six_plus_model
    )
    is not None
    for scalar in range(1, 5)
)
component_six_minus_rational_generator = (sp.Integer(3), sp.Integer(5))
component_six_minus_anti_generator = (
    sp.Integer(-8),
    (-1 + 11 * r) / 2,
)
component_six_minus_boundary_trace = (
    (1 + r) / 2,
    1 - r,
)
for point in (
    component_six_minus_rational_generator,
    component_six_minus_anti_generator,
    component_six_minus_boundary_trace,
):
    assert weierstrass_point_is_on_curve(point, component_six_minus_model)
component_six_minus_generator_relation = weierstrass_add(
    component_six_minus_anti_generator,
    weierstrass_multiply(
        -1,
        component_six_minus_rational_generator,
        component_six_minus_model,
    ),
    component_six_minus_model,
)
assert all(
    polynomial_is_zero(left - right, u)
    for left, right in zip(
        component_six_minus_generator_relation,
        component_six_minus_boundary_trace,
        strict=True,
    )
)

# Each column is the class of one K-prime boundary divisor after pushforward
# to E_+ x E_-.  The first row is the coefficient of T modulo five; the
# next two are the coefficients of G0 and G2.  The q2-node column follows
# from div(q2(x)); all other columns are direct exact residue-field traces.
component_six_plus_classes = sp.Matrix(
    (0, 1, 0, 1, 4, 3, 4, 0, 4, 1, 0, 4, 1, 1, 2, 4, 0, 0, 2, 4)
)
component_six_minus_rational_classes = sp.Matrix(
    (4, -1, 3, -1, 0, 12, 1, 3, 1, 2, 4, 0, 2, 2, 14, 0, 4, 10, 1, 5)
)
component_six_minus_anti_classes = sp.Matrix(
    (0, 1, -1, -1, 0, 0, 1, 1, -1, -1, 1, -1, 0, 1, 0, 1, -1, 0, 0, 0)
)

# Principal fiber checks.  The p3-primes have ramification index three;
# the two degree-three q3-primes have index two.  All other listed affine
# boundary primes are unramified for the first projection.
component_six_boundary_class_rows = (
    component_six_plus_classes,
    component_six_minus_rational_classes,
    component_six_minus_anti_classes,
)
component_six_infinity_indices = (18, 19)
component_six_fiber_checks = (
    ((0, 1, 2), (1, 1, 1), 1),
    ((3, 4, 5, 6), (1, 1, 1, 1), 2),
    ((7, 8), (3, 3), 2),
    ((9, 10), (1, 1), 1),
    ((11, 12, 13, 14), (1, 1, 1, 1), 3),
    ((15, 16, 17), (2, 2, 1), 3),
)
for fiber_indices, multiplicities, polynomial_degree in component_six_fiber_checks:
    for row_index, class_row in enumerate(component_six_boundary_class_rows):
        zero_sum = sum(
            multiplicity * class_row[index]
            for index, multiplicity in zip(
                fiber_indices, multiplicities, strict=True
            )
        )
        pole_sum = polynomial_degree * sum(
            class_row[index] for index in component_six_infinity_indices
        )
        difference = zero_sum - pole_sum
        if row_index == 0:
            assert difference % 5 == 0
        else:
            assert difference == 0

# The (2,2)-isogeny kernel is the even F_2-subspace on the three pairs of
# Weierstrass roots.  The cubic in X=z^2 is irreducible over K and has square
# discriminant (176r)^2, so its cyclic cubic Galois group permutes the three
# nonzero kernel points transitively.  Hence the isogeny has no nonzero
# K-rational kernel point, and vanishing in both elliptic factors is exactly
# vanishing in Jac(C6)(K).
component_six_pair_cubic = sp.Poly(
    1 - 11 * x - 77 * x**2 - 121 * x**3,
    x,
    domain=number_field,
)
_, component_six_pair_factorization = sp.factor_list(
    component_six_pair_cubic.as_expr(), x, extension=r
)
assert len(component_six_pair_factorization) == 1
assert sp.degree(component_six_pair_factorization[0][0], x) == 3
assert component_six_pair_factorization[0][1] == 1
assert polynomial_is_zero(
    sp.discriminant(component_six_pair_cubic.as_expr(), x) - (176 * r) ** 2,
    u,
)

# The exact unit lattice is the simultaneous degree/Mordell--Weil kernel,
# with the E_+ row imposed modulo five.  The displayed row-Hermite basis has
# rank seventeen.  Its Smith form shows
#   Div_S / div(O(U)^*/K^*) = Z^3 + Z/5,
# exactly matching the three free constraints and the order-five class.
component_six_boundary_unit_lattice = sp.Matrix(
    (
        (20, -10, -25, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (10, -6, -13, 7, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (2, -2, -6, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (6, -4, -8, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (0, -1, -1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (6, -3, -8, 4, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (17, -9, -23, 12, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (-1, -1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        (7, -4, -10, 4, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
        (3, -3, -5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
        (17, -10, -23, 13, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
        (7, -5, -13, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
        (7, -5, -10, 5, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
        (19, -10, -26, 14, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
        (-1, -1, -2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
        (11, -6, -15, 8, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
        (5, -4, -9, 4, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
    )
)
assert component_six_boundary_unit_lattice.rank() == 17
for constraint_index, constraint in enumerate(
    (
        component_six_boundary_degrees,
        component_six_minus_rational_classes,
        component_six_minus_anti_classes,
        component_six_plus_classes,
    )
):
    values = component_six_boundary_unit_lattice * constraint
    if constraint_index < 3:
        assert values == sp.zeros(17, 1)
    else:
        assert all(value % 5 == 0 for value in values)
component_six_boundary_unit_smith = smith_normal_form(
    component_six_boundary_unit_lattice,
    domain=ZZ,
)
assert component_six_boundary_unit_smith == sp.Matrix.hstack(
    sp.diag(*([1] * 16 + [5])),
    sp.zeros(17, 3),
)


# ---------------------------------------------------------------------------
# Pullback/cokernel ledgers for both correspondence projections.

# Compactify the six finite colors by adjoining infinity.  This is the full
# unit lattice of either exceptional rational quotient X_+ or X_-.
exceptional_compact_boundary_degrees = sp.Matrix((1, 2, 2, 1, 3, 3, 1))
exceptional_compact_unit_lattice = sp.zeros(6, 7)
for color_index, color_degree in enumerate(
    exceptional_compact_boundary_degrees[:-1]
):
    exceptional_compact_unit_lattice[color_index, color_index] = 1
    exceptional_compact_unit_lattice[color_index, 6] = -color_degree
assert exceptional_compact_unit_lattice.rank() == 6

# The common triangle-base units and the derivative unit, now with their
# pole orders at infinity included.
exceptional_triangle_unit_lattice = sp.Matrix(
    (
        (3, 3, 1, 0, 0, 0, -11),
        (0, 0, 0, 2, 2, 1, -11),
    )
)
exceptional_derivative_unit = sp.Matrix(
    ((2, 2, 0, 1, 1, 0, -10),)
)


def compact_boundary_pullback(number_source_primes, fiber_rows):
    """Return the seven-row divisor pullback matrix for one projection."""

    assert len(fiber_rows) == 7
    matrix = sp.zeros(7, number_source_primes)
    for target_index, source_terms in enumerate(fiber_rows):
        for source_index, ramification_index in source_terms:
            assert matrix[target_index, source_index] == 0
            matrix[target_index, source_index] = ramification_index
    return matrix


def unit_sublattice_quotient(ambient_basis, sublattice):
    """Compute free rank and torsion of ambient/sub using exact coordinates."""

    coordinate_rows = []
    for sublattice_row in sublattice.tolist():
        solution = next(
            iter(sp.linsolve((ambient_basis.T, sp.Matrix(sublattice_row))))
        )
        assert all(coordinate.is_Integer for coordinate in solution)
        coordinate_rows.append(solution)
    coordinate_matrix = sp.Matrix(coordinate_rows)
    smith = smith_normal_form(coordinate_matrix, domain=ZZ)
    rank = coordinate_matrix.rank()
    diagonal = tuple(
        abs(int(smith[index, index]))
        for index in range(min(smith.shape))
        if smith[index, index]
    )
    assert len(diagonal) == rank
    torsion = tuple(invariant for invariant in diagonal if invariant != 1)
    return ambient_basis.rows - rank, torsion, coordinate_matrix


def assert_exact_common_projection_intersection(x_image, y_image):
    """Prove the two rank-six images meet in exactly the two base rows."""

    relation_matrix = x_image.T.row_join(-y_image.T)
    expected_integer_kernel = sp.Matrix(
        (
            (3, 3, 1, 0, 0, 0, 3, 3, 1, 0, 0, 0),
            (0, 0, 0, 2, 2, 1, 0, 0, 0, 2, 2, 1),
        )
    )
    assert relation_matrix.rank() == 10
    assert relation_matrix * expected_integer_kernel.T == sp.zeros(
        relation_matrix.rows,
        2,
    )
    # Columns 2 and 5 form the identity, so this rank-two kernel lattice is
    # saturated in Z^12.  Since the rational nullity is two, it is the full
    # integer kernel, not merely a finite-index sublattice of it.
    assert expected_integer_kernel[:, (2, 5)].det() == 1


# For C5 the labels record both x- and y-colors.  Node labels denote the
# single K-prime containing the two normalized tangent branches.
component_five_x_pullback = compact_boundary_pullback(
    15,
    (
        ((0, 1), (1, 1)),
        ((2, 1), (3, 1), (4, 1)),
        ((5, 1), (6, 3)),
        ((7, 1), (8, 1)),
        ((9, 1), (10, 1), (11, 1)),
        ((12, 2), (13, 1)),
        ((14, 1),),
    ),
)
component_five_y_pullback = compact_boundary_pullback(
    15,
    (
        ((0, 1), (2, 1)),
        ((1, 1), (3, 1), (6, 1)),
        ((4, 3), (5, 1)),
        ((7, 1), (9, 1)),
        ((8, 1), (10, 1), (12, 1)),
        ((11, 2), (13, 1)),
        ((14, 1),),
    ),
)
component_five_boundary_degree_column = sp.Matrix(
    component_five_boundary_degrees
)
for projection_pullback in (
    component_five_x_pullback,
    component_five_y_pullback,
):
    assert (
        projection_pullback * component_five_boundary_degree_column
        == 5 * exceptional_compact_boundary_degrees
    )
assert (
    exceptional_triangle_unit_lattice * component_five_x_pullback
    == exceptional_triangle_unit_lattice * component_five_y_pullback
)
component_five_x_unit_image = (
    exceptional_compact_unit_lattice * component_five_x_pullback
)
component_five_y_unit_image = (
    exceptional_compact_unit_lattice * component_five_y_pullback
)
component_five_two_projection_image = sp.Matrix.vstack(
    component_five_x_unit_image,
    component_five_y_unit_image,
)
assert component_five_x_unit_image.rank() == 6
assert component_five_y_unit_image.rank() == 6
assert component_five_two_projection_image.rank() == 10
assert_exact_common_projection_intersection(
    component_five_x_unit_image,
    component_five_y_unit_image,
)
component_five_projection_quotient = unit_sublattice_quotient(
    component_five_boundary_unit_lattice,
    component_five_two_projection_image,
)
assert component_five_projection_quotient[:2] == (4, ())

# For C6 the y-colors are read from the exact boundary ideals replayed by
# verify_psl2_11_c6_boundary_images.py.  The p3 and q3 coefficients are the
# ramification indices 3 and 2.
component_six_x_pullback = compact_boundary_pullback(
    20,
    (
        ((0, 1), (1, 1), (2, 1)),
        ((3, 1), (4, 1), (5, 1), (6, 1)),
        ((7, 3), (8, 3)),
        ((9, 1), (10, 1)),
        ((11, 1), (12, 1), (13, 1), (14, 1)),
        ((15, 2), (16, 2), (17, 1)),
        ((18, 1), (19, 1)),
    ),
)
component_six_y_pullback = compact_boundary_pullback(
    20,
    (
        ((0, 1), (3, 1), (7, 1)),
        ((1, 1), (4, 1), (5, 1), (8, 1)),
        ((2, 3), (6, 3)),
        ((13, 1), (16, 1)),
        ((9, 1), (12, 1), (14, 1), (15, 1)),
        ((10, 2), (11, 2), (17, 1)),
        ((18, 1), (19, 1)),
    ),
)
for projection_pullback in (
    component_six_x_pullback,
    component_six_y_pullback,
):
    assert (
        projection_pullback * component_six_boundary_degrees
        == 6 * exceptional_compact_boundary_degrees
    )
assert (
    exceptional_triangle_unit_lattice * component_six_x_pullback
    == exceptional_triangle_unit_lattice * component_six_y_pullback
)
component_six_x_unit_image = (
    exceptional_compact_unit_lattice * component_six_x_pullback
)
component_six_y_unit_image = (
    exceptional_compact_unit_lattice * component_six_y_pullback
)
component_six_two_projection_image = sp.Matrix.vstack(
    component_six_x_unit_image,
    component_six_y_unit_image,
)
assert component_six_x_unit_image.rank() == 6
assert component_six_y_unit_image.rank() == 6
assert component_six_two_projection_image.rank() == 10
assert_exact_common_projection_intersection(
    component_six_x_unit_image,
    component_six_y_unit_image,
)
component_six_projection_quotient = unit_sublattice_quotient(
    component_six_boundary_unit_lattice,
    component_six_two_projection_image,
)
assert component_six_projection_quotient[:2] == (7, ())

# The smallest two-output divisor ledger adds the two pulled-back derivative
# rows to the common rank-two triangle image.  It is primitive in both full
# source-unit lattices and has rank four: the two derivative classes are
# independent modulo the common base.  This passes the weak determinant
# ledger, but it cannot be a full character completion.  Even after starting
# from the whole rank-ten two-projection image, two extra mask rows leave at
# least ranks two and five on C5 and C6 respectively.
for ambient_lattice, x_pullback, y_pullback, expected_free_rank in (
    (
        component_five_boundary_unit_lattice,
        component_five_x_pullback,
        component_five_y_pullback,
        10,
    ),
    (
        component_six_boundary_unit_lattice,
        component_six_x_pullback,
        component_six_y_pullback,
        13,
    ),
):
    common_triangle_image = exceptional_triangle_unit_lattice * x_pullback
    assert (
        common_triangle_image
        == exceptional_triangle_unit_lattice * y_pullback
    )
    two_output_ledger = sp.Matrix.vstack(
        common_triangle_image,
        exceptional_derivative_unit * x_pullback,
        exceptional_derivative_unit * y_pullback,
    )
    assert two_output_ledger.rank() == 4
    assert unit_sublattice_quotient(
        ambient_lattice,
        two_output_ledger,
    )[:2] == (expected_free_rank, ())
assert component_five_boundary_unit_lattice.rows - 12 == 2
assert component_six_boundary_unit_lattice.rows - 12 == 5

# On X_0(11), the two triangle pullbacks leave one free primitive direction.
# Adding the Miller unit closes the rank but leaves the exact order-six
# saturation quotient generated by the missing square/cube roots.
elliptic_triangle_pullback_lattice = sp.Matrix(
    (
        (3, 0, -11, -1),
        (0, 2, -11, -1),
    )
)
elliptic_evident_unit_lattice = sp.Matrix.vstack(
    elliptic_triangle_pullback_lattice,
    sp.Matrix(((0, 0, 5, -5),)),
)
assert unit_sublattice_quotient(
    elliptic_boundary_unit_lattice,
    elliptic_triangle_pullback_lattice,
)[:2] == (1, ())
assert unit_sublattice_quotient(
    elliptic_boundary_unit_lattice,
    elliptic_evident_unit_lattice,
)[:2] == (0, (6,))


# ---------------------------------------------------------------------------
# Residual mask bases and the first factor-rich two-output support.


def boundary_permutation_matrix(permutation):
    """Return the row-action matrix for a permutation of boundary primes."""

    matrix = sp.zeros(len(permutation), len(permutation))
    for source_index, target_index in enumerate(permutation):
        matrix[source_index, target_index] = 1
    assert matrix * matrix == sp.eye(len(permutation))
    return matrix


def unit_lattice_coordinate_action(ambient_basis, boundary_action):
    """Express a boundary permutation as an action in a unit-lattice basis."""

    transformed_basis = ambient_basis * boundary_action
    free_rank, torsion, coordinates = unit_sublattice_quotient(
        ambient_basis,
        transformed_basis,
    )
    assert (free_rank, torsion) == (0, ())
    assert coordinates * coordinates == sp.eye(ambient_basis.rows)
    return coordinates


def nonnegative_weighted_compositions(weights, total):
    """Enumerate every nonnegative vector with the prescribed weight."""

    if len(weights) == 1:
        if total % weights[0] == 0:
            yield (total // weights[0],)
        return
    for first_coefficient in range(total // weights[0] + 1):
        for tail in nonnegative_weighted_compositions(
            weights[1:],
            total - first_coefficient * weights[0],
        ):
            yield (first_coefficient,) + tail


# These maps send coordinates in the displayed bases of L5 and L6 to
# coordinates in their free projection cokernels.  The exact kernel checks
# below prove that they are quotient maps, rather than relying on a chosen
# Smith transformation at replay time.
component_five_residual_coordinate_map = sp.Matrix(
    (
        (0, 1, 0, 0),
        (0, 1, 0, 0),
        (0, 1, -1, 0),
        (0, 0, 1, 0),
        (0, 2, -3, 0),
        (0, 0, 1, 0),
        (1, -2, 0, 1),
        (-1, 3, 0, -1),
        (-1, 3, 0, -1),
        (0, 0, 0, 1),
        (1, 0, 0, 0),
        (1, 0, 0, 0),
        (-2, 3, 0, 0),
        (0, 1, 0, 0),
    )
)
component_six_residual_coordinate_map = sp.Matrix(
    (
        (6, 6, 3, -6, -6, 0, 0),
        (4, 4, 2, -4, -4, 0, -1),
        (0, 0, 0, 0, 0, 0, 1),
        (2, 2, 1, -2, -2, 0, 0),
        (0, 0, 0, 0, 0, 0, 0),
        (2, 2, 1, -2, -2, 0, 0),
        (0, 0, 0, 1, 1, 0, 0),
        (0, 0, 0, 0, 0, 0, 0),
        (-3, -3, -2, 4, 4, 0, 0),
        (3, 4, 2, -4, -4, -1, 0),
        (6, 5, 3, -5, -5, 0, 0),
        (0, 0, 0, 0, 0, 1, 0),
        (3, 2, 1, -2, -2, 0, 0),
        (0, 1, 0, 0, 0, 0, 0),
        (0, 0, 1, 0, 0, 0, 0),
        (0, 0, 0, 1, 0, 0, 0),
        (0, 0, 0, 0, 1, 0, 0),
    )
)
assert (
    component_five_projection_quotient[2]
    * component_five_residual_coordinate_map
) == sp.zeros(12, 4)
assert (
    component_six_projection_quotient[2]
    * component_six_residual_coordinate_map
) == sp.zeros(12, 7)

# Sparse complements certify surjectivity of the quotient maps.  Since the
# two projection images are already primitive, their rank-ten coordinate
# lattices are exactly the kernels.
component_five_complement_indices = (3, 6, 7, 9)
component_six_complement_indices = (2, 5, 8, 9, 12, 15, 16)
component_five_complement_quotient = component_five_residual_coordinate_map[
    component_five_complement_indices,
    :,
]
component_six_complement_quotient = component_six_residual_coordinate_map[
    component_six_complement_indices,
    :,
]
assert abs(component_five_complement_quotient.det()) == 1
assert abs(component_six_complement_quotient.det()) == 1
assert unit_sublattice_quotient(
    component_five_boundary_unit_lattice,
    sp.Matrix.vstack(
        component_five_two_projection_image,
        component_five_boundary_unit_lattice[
            component_five_complement_indices,
            :,
        ],
    ),
)[:2] == (0, ())
assert unit_sublattice_quotient(
    component_six_boundary_unit_lattice,
    sp.Matrix.vstack(
        component_six_two_projection_image,
        component_six_boundary_unit_lattice[
            component_six_complement_indices,
            :,
        ],
    ),
)[:2] == (0, ())

# Swapping the two projections also conjugates r.  On K-prime labels it is
# the following involution.  The exact pullback matrices are exchanged.
component_five_projection_swap = boundary_permutation_matrix(
    (0, 2, 1, 3, 6, 5, 4, 7, 9, 8, 10, 12, 11, 13, 14)
)
component_six_projection_swap = boundary_permutation_matrix(
    (0, 3, 7, 1, 4, 5, 8, 2, 6, 13, 16, 15, 12, 9, 14, 11, 10, 17, 18, 19)
)
assert (
    component_five_x_pullback * component_five_projection_swap
    == component_five_y_pullback
)
assert (
    component_six_x_pullback * component_six_projection_swap
    == component_six_y_pullback
)
component_five_unit_swap = unit_lattice_coordinate_action(
    component_five_boundary_unit_lattice,
    component_five_projection_swap,
)
component_six_unit_swap = unit_lattice_coordinate_action(
    component_six_boundary_unit_lattice,
    component_six_projection_swap,
)
component_five_complement_swap_quotient = (
    component_five_unit_swap[component_five_complement_indices, :]
    * component_five_residual_coordinate_map
)
component_six_complement_swap_quotient = (
    component_six_unit_swap[component_six_complement_indices, :]
    * component_six_residual_coordinate_map
)
component_five_residual_swap = (
    component_five_complement_swap_quotient
    * component_five_complement_quotient.inv()
)
component_six_residual_swap = (
    component_six_complement_swap_quotient
    * component_six_complement_quotient.inv()
)
assert component_five_residual_swap == sp.eye(4)
assert component_six_residual_swap == sp.Matrix(
    (
        (1, 0, 0, 0, 0, 0, 0),
        (0, 1, 0, 0, 0, 0, 0),
        (0, 0, -1, 0, 2, 0, 0),
        (0, 0, 0, 1, 0, 0, 0),
        (0, 0, 0, 0, 1, 0, 0),
        (0, 0, -1, 0, 1, 1, 0),
        (0, 0, -1, 0, 1, 0, 1),
    )
)

# In the chosen C6 complement, the first five rows below are fixed and the
# last two are exchanged.  Its determinant one proves the integral module
# decomposition Q6 = Z^5_triv + Z[C2], not merely a rational eigensplitting.
component_six_permutation_basis = sp.Matrix(
    (
        (0, 0, 0, 0, 0, -1, 1),
        (0, 0, 0, 0, 1, 0, 0),
        (0, 0, 0, 1, 0, 0, 0),
        (0, 1, 0, 0, 0, 0, 0),
        (1, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, -1),
        (0, 0, 1, 0, -1, 0, -1),
    )
)
component_six_permutation_action = sp.diag(*([1] * 5), 0, 0)
component_six_permutation_action[5, 6] = 1
component_six_permutation_action[6, 5] = 1
assert abs(component_six_permutation_basis.det()) == 1
assert (
    component_six_permutation_basis * component_six_residual_swap
    == component_six_permutation_action * component_six_permutation_basis
)
assert (
    4 - (component_five_residual_swap + sp.eye(4)).rank(),
    4 - (component_five_residual_swap - sp.eye(4)).rank(),
) == (0, 4)
assert (
    7 - (component_six_residual_swap + sp.eye(7)).rank(),
    7 - (component_six_residual_swap - sp.eye(7)).rank(),
) == (1, 6)

# Effective regular-mask representatives: all finite-boundary coefficients
# are nonnegative and only the infinity primes carry poles.  The C5 basis is
# optimal in pole order.  Exhausting all degree-five effective zero divisors
# (simple pole at infinity_5) gives Smith index two, while three such masks
# plus the displayed double-pole mask give determinant one.
component_five_effective_masks = sp.Matrix(
    (
        (0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, -1),
        (0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, -1),
        (1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -1),
        (0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, -2),
    )
)
component_five_effective_coordinates = unit_sublattice_quotient(
    component_five_boundary_unit_lattice,
    component_five_effective_masks,
)[2]
component_five_effective_quotient = (
    component_five_effective_coordinates
    * component_five_residual_coordinate_map
)
assert component_five_effective_quotient.det() == -1
assert all(
    coefficient >= 0
    for coefficient in component_five_effective_masks[:, :14]
)
assert all(
    coefficient < 0
    for coefficient in component_five_effective_masks[:, 14]
)
component_five_simple_pole_masks = sp.Matrix(
    tuple(
        finite_divisor + (-1,)
        for finite_divisor in nonnegative_weighted_compositions(
            component_five_boundary_degrees[:-1],
            5,
        )
    )
)
assert component_five_simple_pole_masks.rows == 26
component_five_simple_pole_coordinates = unit_sublattice_quotient(
    component_five_boundary_unit_lattice,
    component_five_simple_pole_masks,
)[2]
component_five_simple_pole_quotient = (
    component_five_simple_pole_coordinates
    * component_five_residual_coordinate_map
)
component_five_simple_pole_smith = smith_normal_form(
    component_five_simple_pole_quotient,
    domain=ZZ,
)
assert tuple(
    abs(int(component_five_simple_pole_smith[index, index]))
    for index in range(4)
) == (1, 1, 1, 2)

# A fully effective C6 permutation-basis lift.  The first five quotient rows
# are fixed and the last two are exchanged; their particular divisors need
# only agree with that action modulo the two projection images.
component_six_effective_masks = sp.Matrix(
    (
        (0, 4, 8, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -7, -5),
        (8, 0, 0, 1, 0, 0, 5, 0, 4, 0, 0, 1, 0, 1, 0, 0, 0, 0, -7, -7),
        (5, 1, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, -4, -4),
        (5, 0, 0, 0, 0, 0, 5, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5, -5),
        (2, 1, 0, 1, 0, 1, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -4),
        (0, 5, 10, 0, 0, 0, 0, 4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -6, -7),
        (0, 0, 4, 5, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -6, -7),
    )
)
component_six_effective_coordinates = unit_sublattice_quotient(
    component_six_boundary_unit_lattice,
    component_six_effective_masks,
)[2]
component_six_effective_quotient = (
    component_six_effective_coordinates
    * component_six_residual_coordinate_map
)
assert component_six_effective_quotient == (
    component_six_permutation_basis * component_six_complement_quotient
)
assert component_six_effective_quotient.det() == -1
assert all(
    coefficient >= 0
    for coefficient in component_six_effective_masks[:, :18]
)
assert all(
    coefficient <= 0
    for coefficient in component_six_effective_masks[:, 18:]
)

# The difference of the two infinity coefficients vanishes on both
# projection images, so it descends to an intrinsic homomorphism on the
# rank-seven residual quotient.  In the effective permutation basis its
# values are (-2,0,0,0,0,1,1).  In particular the N1 and exchanged N6/N7
# classes cannot be represented by a divisor with equal infinity orders.
assert all(
    row[18] == row[19]
    for row in component_six_two_projection_image.tolist()
)
component_six_infinity_imbalance = tuple(
    int(mask[18] - mask[19])
    for mask in component_six_effective_masks.tolist()
)
assert component_six_infinity_imbalance == (-2, 0, 0, 0, 0, 1, 1)

# If every residual mask occupies a distinct nonconstant normal monomial and
# the zero section must preserve the old degree-eleven fields, the smallest
# two-normal supports have degrees two on C5 and three on C6.  C6 uses all
# five nonconstant quadratic-or-lower slots for its fixed classes and the
# exchanged cubic pair s^3/t^3.
component_five_normal_support = ((1, 0), (0, 1), (2, 0), (1, 1))
component_six_fixed_normal_support = (
    (1, 0),
    (0, 1),
    (2, 0),
    (1, 1),
    (0, 2),
)
component_six_exchanged_normal_support = ((3, 0), (0, 3))
number_nonconstant_two_variable_monomials = lambda degree: (
    (degree + 1) * (degree + 2) // 2 - 1
)
assert number_nonconstant_two_variable_monomials(1) < len(
    component_five_normal_support
) <= number_nonconstant_two_variable_monomials(2)
assert number_nonconstant_two_variable_monomials(2) < (
    len(component_six_fixed_normal_support) + 1
) <= number_nonconstant_two_variable_monomials(3)

# Good reduction at the split primes 3 and 5 proves that E_+(K)_tors is
# exactly <T> of order five and that E_-(K) has trivial torsion.
assert (
    count_weierstrass_points_mod_prime(component_six_plus_model, 3),
    count_weierstrass_points_mod_prime(component_six_plus_model, 5),
) == (5, 5)
assert (
    count_weierstrass_points_mod_prime(component_six_minus_model, 3),
    count_weierstrass_points_mod_prime(component_six_minus_model, 5),
) == (5, 9)


def singular_polynomial(expression):
    """Translate an exact SymPy polynomial over Q(sqrt(-11)) to Singular."""

    return (
        str(sp.Poly(expression, x, y, domain=number_field).as_expr())
        .replace("sqrt(11)*I", "r")
        .replace("**", "^")
    )


def singular_polynomial_coefficient_first(expression):
    """Translate with every rational coefficient before its monomial."""

    polynomial = sp.Poly(expression, x, y, domain=number_field)
    terms = []
    for (x_degree, y_degree), coefficient in polynomial.terms():
        term = (
            "("
            + str(coefficient)
            .replace("sqrt(11)*I", "r")
            .replace("**", "^")
            + ")"
        )
        if x_degree:
            term += f"*x^{x_degree}"
        if y_degree:
            term += f"*y^{y_degree}"
        terms.append(term)
    return "+".join(terms)


def singular_homogeneous_evaluation(expression, total_degree):
    """Translate F(Xo,Yo,Zo) to Singular without adding three variables."""

    polynomial = sp.Poly(expression, x, y, domain=number_field)
    terms = []
    for (x_degree, y_degree), coefficient in polynomial.terms():
        term = (
            "("
            + str(coefficient)
            .replace("sqrt(11)*I", "rr")
            .replace("**", "^")
            + ")"
        )
        if x_degree:
            term += f"*Xo^{x_degree}"
        if y_degree:
            term += f"*Yo^{y_degree}"
        z_degree = total_degree - x_degree - y_degree
        if z_degree:
            term += f"*Zo^{z_degree}"
        terms.append(term)
    return "+".join(terms)


def run_singular_normalization_check():
    """Replay normalization, canonical elimination, and Cremona reduction."""

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for --singular-normalization")

    programs = []
    for label, correspondence, expected_delta, support in (
        (
            "C5",
            correspondence_five,
            5,
            "ideal V=std(intersect(ideal(p2x,p1my),"
            "intersect(ideal(q1x,q1my),ideal(p1x,p2my))));",
        ),
        (
            "C6",
            correspondence_six,
            8,
            "poly ell=18*x+(r+1)*y2+(8*r-10)*y+88-20*r;"
            "ideal V=std(intersect(ideal(p2x,p2my),"
            "intersect(ideal(p1x,p1my),ideal(q2my,ell))));",
        ),
    ):
        programs.extend(
            (
                "ring R=(0,r),(x,y),dp; minpoly=r2+11;",
                f"poly f={singular_polynomial(correspondence)}; ideal I=f;",
                "ideal S=std(radical(I+jacob(I)));",
                "poly p1x=2*x+11-3*r;",
                "poly p2x=2*x2-(11-3*r)*x-(22+6*r);",
                "poly q1x=2*x+5+3*r;",
                "poly p1my=2*y+11+3*r;",
                "poly p2my=2*y2-(11+3*r)*y-(22-6*r);",
                "poly q1my=2*y+5-3*r;",
                "poly q2my=2*y3+(15+3*r)*y2-(12+12*r)*y+56-96*r;",
                support,
                f'if (vdim(S)!={expected_delta}) '
                f'{{print("FAIL {label} singular length"); exit(1);}}',
                f'if (size(reduce(S,V))!=0 || size(reduce(V,S))!=0) '
                f'{{print("FAIL {label} singular support"); exit(1);}}',
                "ideal C=std(normalConductor(I));",
                f'if (vdim(C)!={expected_delta}) '
                f'{{print("FAIL {label} conductor length"); exit(1);}}',
                f'if (size(reduce(C,S))!=0 || size(reduce(S,C))!=0) '
                f'{{print("FAIL {label} conductor support"); exit(1);}}',
                "poly h=diff(diff(f,x),x)*diff(diff(f,y),y)"
                "-diff(diff(f,x),y)^2;",
                "ideal T=std(S+ideal(h));",
                f'if (T[1]!=1) {{print("FAIL {label} nonnode"); exit(1);}}',
                f'print("PASS Singular normalization: {label} has '
                f'{expected_delta} affine nodes and reduced conductor");',
                "kill R;",
            )
        )

    # Decompose the six affine boundary fibers on the plane quintic and
    # then separate both tangent branches at the three nodal boundary
    # orbits.  Replacing the plane-node degrees 2,2,1 by the normalized
    # branch degrees 4,4,2 gives the fifteen K-prime degrees recorded above.
    programs.extend(
        (
            "ring R=(0,r),(x,y,m),lp; minpoly=r2+11;",
            f"poly f={singular_polynomial(correspondence_five)};",
            "proc countDegree(list L,int d)"
            "{int c=0;int i;for(i=1;i<=size(L);i++)"
            "{if(vdim(std(L[i]))==d){c++;}}return(c);}",
            "poly p1x=2*x+11-3*r;",
            "poly p2x=2*x2-(11-3*r)*x-(22+6*r);",
            "poly p3x=x2+11*x+55+9*r;",
            "poly q1x=2*x+5+3*r;",
            "poly q2x=2*x3+(15-3*r)*x2-(12-12*r)*x+56+96*r;",
            "poly q3x=2*x3-18*x2+(21+45*r)*x-(175+279*r);",
            "ideal J=std(radical(ideal(f,p1x,m)));list L=minAssGTZ(J);",
            "if(vdim(J)!=3 || size(L)!=2 || countDegree(L,1)!=1"
            " || countDegree(L,2)!=1)"
            '{print("FAIL C5 p1 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,p2x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=8 || size(L)!=3 || countDegree(L,2)!=2"
            " || countDegree(L,4)!=1)"
            '{print("FAIL C5 p2 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,p3x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=6 || size(L)!=2 || countDegree(L,2)!=1"
            " || countDegree(L,4)!=1)"
            '{print("FAIL C5 p3 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,q1x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=4 || size(L)!=2 || countDegree(L,1)!=1"
            " || countDegree(L,3)!=1)"
            '{print("FAIL C5 q1 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,q2x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=15 || size(L)!=3 || countDegree(L,3)!=1"
            " || countDegree(L,6)!=2)"
            '{print("FAIL C5 q2 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,q3x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=9 || size(L)!=2 || countDegree(L,3)!=1"
            " || countDegree(L,6)!=1)"
            '{print("FAIL C5 q3 boundary primes");exit(1);}',
            "poly tangent=diff(diff(f,x),x)+2*diff(diff(f,x),y)*m"
            "+diff(diff(f,y),y)*m2;",
            "ideal B=std(ideal(p2x,2*y+11+3*r,tangent));",
            "L=minAssGTZ(B);",
            "if(vdim(B)!=4 || size(L)!=1 || countDegree(L,4)!=1)"
            '{print("FAIL C5 first node branches");exit(1);}',
            "B=std(ideal(p1x,2*y2-(11+3*r)*y-(22-6*r),tangent));",
            "L=minAssGTZ(B);",
            "if(vdim(B)!=4 || size(L)!=1 || countDegree(L,4)!=1)"
            '{print("FAIL C5 second node branches");exit(1);}',
            "B=std(ideal(q1x,2*y+5-3*r,tangent));L=minAssGTZ(B);",
            "if(vdim(B)!=2 || size(L)!=1 || countDegree(L,2)!=1)"
            '{print("FAIL C5 rational-node branches");exit(1);}',
            'print("PASS Singular boundary decomposition: C5 has 15 K-prime punctures");',
            "kill R;",
        )
    )

    # The same calculation on C6 yields twenty normalized K-prime
    # punctures.  The radical plane fibers have degrees
    # 5,8,4,6,15,12; replacing the three nodal orbits by their tangent
    # branches changes these to normalized totals 6,12,4,6,18,12.
    programs.extend(
        (
            "ring R=(0,r),(x,y,m),lp; minpoly=r2+11;",
            f"poly f={singular_polynomial(correspondence_six)};",
            "proc countDegree(list L,int d)"
            "{int c=0;int i;for(i=1;i<=size(L);i++)"
            "{if(vdim(std(L[i]))==d){c++;}}return(c);}",
            "poly p1x=2*x+11-3*r;",
            "poly p2x=2*x2-(11-3*r)*x-(22+6*r);",
            "poly p3x=x2+11*x+55+9*r;",
            "poly q1x=2*x+5+3*r;",
            "poly q2x=2*x3+(15-3*r)*x2-(12-12*r)*x+56+96*r;",
            "poly q3x=2*x3-18*x2+(21+45*r)*x-(175+279*r);",
            "poly p1my=2*y+11+3*r;",
            "poly p2my=2*y2-(11+3*r)*y-(22-6*r);",
            "poly q2my=2*y3+(15+3*r)*y2-(12+12*r)*y+56-96*r;",
            "ideal J=std(radical(ideal(f,p1x,m)));list L=minAssGTZ(J);",
            "if(vdim(J)!=5 || size(L)!=3 || countDegree(L,1)!=1"
            " || countDegree(L,2)!=2)"
            '{print("FAIL C6 p1 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,p2x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=8 || size(L)!=3 || countDegree(L,2)!=2"
            " || countDegree(L,4)!=1)"
            '{print("FAIL C6 p2 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,p3x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=4 || size(L)!=2 || countDegree(L,2)!=2)"
            '{print("FAIL C6 p3 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,q1x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=6 || size(L)!=2 || countDegree(L,3)!=2)"
            '{print("FAIL C6 q1 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,q2x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=15 || size(L)!=4 || countDegree(L,3)!=3"
            " || countDegree(L,6)!=1)"
            '{print("FAIL C6 q2 boundary primes");exit(1);}',
            "J=std(radical(ideal(f,q3x,m)));L=minAssGTZ(J);",
            "if(vdim(J)!=12 || size(L)!=3 || countDegree(L,3)!=2"
            " || countDegree(L,6)!=1)"
            '{print("FAIL C6 q3 boundary primes");exit(1);}',
            "poly tangent=diff(diff(f,x),x)+2*diff(diff(f,x),y)*m"
            "+diff(diff(f,y),y)*m2;",
            "ideal N=std(ideal(p1x,p1my,tangent));L=minAssGTZ(N);",
            "if(vdim(N)!=2 || size(L)!=1 || countDegree(L,2)!=1)"
            '{print("FAIL C6 rational-node branches");exit(1);}',
            "N=std(ideal(p2x,p2my,tangent));L=minAssGTZ(N);",
            "if(vdim(N)!=8 || size(L)!=2 || countDegree(L,4)!=2)"
            '{print("FAIL C6 degree-four node branches");exit(1);}',
            "poly ell=18*x+(r+1)*y2+(8*r-10)*y+88-20*r;",
            "N=std(ideal(q2my,ell,tangent));L=minAssGTZ(N);",
            "if(vdim(N)!=6 || size(L)!=1 || countDegree(L,6)!=1)"
            '{print("FAIL C6 degree-six node branches");exit(1);}',
            'print("PASS Singular boundary decomposition: C6 has 20 K-prime punctures");',
            "kill R;",
        )
    )

    # Eliminate the plane sextic against the canonical pencil B-tA.  The
    # node contributions are the exact powers l^2*q^2*p^4.  What remains is
    # quadratic in x; its discriminant is a sextic times a cubic square.
    programs.extend(
        (
            "ring R=(0,r),(x,y,t),dp; minpoly=r2+11;",
            f"poly f={singular_polynomial(correspondence_six)};",
            "poly A="
            f"{singular_polynomial_coefficient_first(component_six_adjoint_a)};",
            "poly B="
            f"{singular_polynomial_coefficient_first(component_six_adjoint_b)};",
            "poly Res=resultant(f,B-t*A,y);",
            "poly ell=x+(-3*r+11)/2;",
            "poly q=x3+(-3*r+15)/2*x2+(6*r-6)*x+48*r+28;",
            "poly p=x2+(3*r-11)/2*x-3*r-11;",
            "poly nodes=ell^2*q^2*p^4;",
            "poly Q=Res/nodes;",
            "if(Res-Q*nodes!=0)"
            '{print("FAIL C6 canonical node factors"); exit(1);}',
            "poly qa=diff(diff(Q,x),x)/2;",
            "poly qb=diff(Q,x)-2*qa*x;",
            "poly qc=Q-qa*x^2-qb*x;",
            "if(diff(qa,x)!=0 || diff(qb,x)!=0 || diff(qc,x)!=0)"
            '{print("FAIL C6 canonical residual degree"); exit(1);}',
            "poly disc=qb^2-4*qa*qc;",
            "poly sext=t6+(3/4*r-69/52)*t5"
            "+(-345/416*r-659/416)*t4"
            "+(7/104*r+357/208)*t3"
            "+(489/6656*r-2675/6656)*t2"
            "+(-243/13312*r+615/13312)*t"
            "+15/13312*r-71/53248;",
            "poly square=t3+(1/40*r-7/40)*t2"
            "+(7/160*r-19/160)*t-1/320*r+1/160;",
            "poly expected=(7907328*r+2635776)*sext*square^2;",
            "if(disc!=expected)"
            '{print("FAIL C6 canonical discriminant"); exit(1);}',
            'print("PASS Singular canonical elimination: C6 has the pinned '
            'hyperelliptic sextic");',
            "kill R;",
        )
    )

    homogeneous_component_five = singular_homogeneous_evaluation(
        correspondence_five, 5
    )
    programs.extend(
        (
            'LIB "linalg.lib";',
            "ring R=(0,a),(A,B,C,D,E,F),dp;",
            "minpoly=a8-132*a6+6182*a4+251196*a2+833569;",
            "number rr=9*a7/1504624-1105*a5/1504624"
            "+4477*a3/136784+254411*a/136784;",
            "number ss=-9*a7/3009248+3*a6/72512+1105*a5/3009248"
            "-437*a4/72512-4477*a3/273568+2179*a2/6592"
            "-117627*a/273568+39827/6592;",
            "number tt=-9*a7/3009248-3*a6/72512+1105*a5/3009248"
            "+437*a4/72512-4477*a3/273568-2179*a2/6592"
            "-117627*a/273568-39827/6592;",
            "if(rr*rr+11!=0 || ss*ss-(22-2*rr)!=0"
            " || tt*tt-(22+2*rr)!=0)"
            '{print("FAIL primitive node field"); quit;}',
            "number x1=(11-3*rr+3*ss)/4;",
            "number x2=(11-3*rr-3*ss)/4;",
            "number yA=-(11+3*rr)/2;",
            "number x0=-(5+3*rr)/2;",
            "number y0=-(5-3*rr)/2;",
            "matrix M[3][3]=x1,x2,x0,yA,yA,y0,1,1,1;",
            "matrix Mi=inverse(M);",
            "poly Xo=x1*B*C+x2*C*A+x0*A*B;",
            "poly Yo=yA*B*C+yA*C*A+y0*A*B;",
            "poly Zo=B*C+C*A+A*B;",
            f"poly Graw={homogeneous_component_five};",
            "poly G=Graw/(A2*B2*C2);",
            "if(Graw-G*A2*B2*C2!=0)"
            '{print("FAIL first Cremona transform"); quit;}',
            "number xB=-(11-3*rr)/2;",
            "number y3=(11+3*rr+3*tt)/4;",
            "number y4=(11+3*rr-3*tt)/4;",
            "matrix P3[3][1]=xB,y3,1;",
            "matrix P4[3][1]=xB,y4,1;",
            "matrix q3=Mi*P3; matrix q4=Mi*P4;",
            "poly q31=q3[2,1]*q3[3,1];",
            "poly q32=q3[3,1]*q3[1,1];",
            "poly q33=q3[1,1]*q3[2,1];",
            "poly q41=q4[2,1]*q4[3,1];",
            "poly q42=q4[3,1]*q4[1,1];",
            "poly q43=q4[1,1]*q4[2,1];",
            "poly AA=q31*E*F+q41*F*D;",
            "poly BB=q32*E*F+q42*F*D;",
            "poly CC=q33*E*F+q43*F*D+D*E;",
            "poly Hraw=subst(subst(subst(G,A,AA),B,BB),C,CC);",
            "poly H=Hraw/(D2*E2*F);",
            "if(Hraw-H*D2*E2*F!=0)"
            '{print("FAIL second Cremona transform"); quit;}',
            "number c210=264627/8-29403/8*tt+8019/8*rr"
            "+2673/4*rr*tt;",
            "number c201=617463/16-65043/32*rr*tt-421443/32*tt"
            "+376893/16*rr;",
            "number c120=264627/8-2673/4*rr*tt+29403/8*tt"
            "+8019/8*rr;",
            "number c111=78408+60588*rr;",
            "number c102=-16335/2-5049/4*rr*tt-35937/4*tt"
            "+98307/2*rr;",
            "number c021=617463/16+421443/32*tt+65043/32*rr*tt"
            "+376893/16*rr;",
            "number c012=-16335/2+35937/4*tt+5049/4*rr*tt"
            "+98307/2*rr;",
            "number c003=-28314+16434*rr;",
            "poly Hp=c210*D2*E+c201*D2*F+c120*D*E2"
            "+c111*D*E*F+c102*D*F2+c021*E2*F+c012*E*F2+c003*F3;",
            "if(H-Hp!=0)"
            '{print("FAIL sparse cubic certificate"); quit;}',
            'print("PASS Singular Cremona reduction: C5 has sparse cubic j=-121");',
            "kill R;",
        )
    )

    program = (
        'LIB "normal.lib";\nLIB "primdec.lib";\n'
        + "\n".join(programs)
        + "\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=False,
    )
    singular_output = result.stdout + result.stderr
    required_singular_markers = (
        "PASS Singular normalization: C5",
        "PASS Singular normalization: C6",
        "PASS Singular boundary decomposition: C5",
        "PASS Singular boundary decomposition: C6",
        "PASS Singular canonical elimination: C6",
        "PASS Singular Cremona reduction: C5",
    )
    if (
        result.returncode != 0
        or "FAIL" in singular_output
        or "   ?" in singular_output
        or any(marker not in singular_output for marker in required_singular_markers)
    ):
        raise RuntimeError(
            "Singular normalization replay failed:\n"
            + result.stdout
            + result.stderr
        )
    for line in result.stdout.splitlines():
        if line.startswith("PASS Singular"):
            print(line)

    # Replay the actual Riemann--Roch representatives separately.  Singular
    # does not reliably return a nonzero status for errors raised inside a
    # procedure, so require every terminal marker and reject both explicit
    # failures and interpreter diagnostics.
    mask_script = Path(__file__).with_name(
        "verify_psl2_11_normalization_masks.sing"
    )
    mask_result = subprocess.run(
        [singular, "-q", str(mask_script)],
        text=True,
        capture_output=True,
        check=False,
    )
    mask_output = mask_result.stdout + mask_result.stderr
    required_mask_markers = (
        "PASS C5 explicit normalization masks M1--M4",
        "PASS C6 filtered Riemann--Roch dimensions 23, 29, and 41",
        "PASS C6 uniform-pole masks N2--N5",
        "PASS C6 exchanged masks N6--N7 with forced pole pair (6,7)",
        "PASS C6 asymmetric mask N1 with pole pair (7,5)",
        "PASS explicit normalization representatives for all C5/C6 masks",
    )
    if (
        mask_result.returncode != 0
        or "FAIL" in mask_output
        or "   ?" in mask_output
        or any(marker not in mask_output for marker in required_mask_markers)
    ):
        raise RuntimeError(
            "Singular normalization-mask replay failed:\n" + mask_output
        )
    for line in mask_result.stdout.splitlines():
        if line.startswith("PASS C5") or line.startswith("PASS C6"):
            print(line)
    print("PASS Singular normalization-module representatives: all C5/C6 masks")


def run_pari_mordell_weil_check():
    """Prove the rational ranks used after quadratic base change to K."""

    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError("PARI/GP is required for --pari-mordell-weil")
    program = r"""
EX=ellinit([0,-1,1,-10,-20]);
TX=elltwist(EX,-11);
R=ellrank(EX); if(R[1]!=0 || R[2]!=0,error("X0(11) rational rank"));
R=ellrank(TX); if(R[1]!=0 || R[2]!=0,error("X0(11) twist rank"));
if(elltors(EX)[1]!=5 || elltors(TX)[1]!=1,error("X0(11) torsion"));
E5=ellinit([1,1,0,-2,-7]);
T5=elltwist(E5,-11);
R=ellrank(E5); if(R[1]!=0 || R[2]!=0,error("C5 rational rank"));
R=ellrank(T5); if(R[1]!=0 || R[2]!=0,error("C5 twist rank"));
if(elltors(E5)[1]!=1 || elltors(T5)[1]!=1,error("C5 rational torsion"));
Ep=ellinit([0,2,1,1,0]);
Tp=elltwist(Ep,-11);
R=ellrank(Ep); if(R[1]!=0 || R[2]!=0,error("Eplus rational rank"));
R=ellrank(Tp); if(R[1]!=0 || R[2]!=0,error("Eplus twist rank"));
if(elltors(Ep)[1]!=5 || elltors(Tp)[1]!=1,error("Eplus torsion"));
Em=ellinit([0,2,1,-6,3]);
Tm=elltwist(Em,-11);
R=ellrank(Em); if(R[1]!=1 || R[2]!=1,error("Eminus rational rank"));
R=ellrank(Tm); if(R[1]!=1 || R[2]!=1,error("Eminus twist rank"));
if(elltors(Em)[1]!=1 || elltors(Tm)[1]!=1,error("Eminus torsion"));
if(ellheegner(Em)!=[3,5],error("Eminus Heegner generator"));
if(ellheegner(Tm)!=[88,660],error("Eminus twist Heegner generator"));
if(ellmul(Em,[3,5],3)!=[1/4,7/8],error("Eminus generator relation"));
print("PASS PARI Mordell-Weil: X0 and C5 rank 0; Eplus rank 0; Eminus rank 2 over K");
"""
    result = subprocess.run(
        [gp, "-fq"],
        input=program,
        text=True,
        capture_output=True,
        check=False,
    )
    output = result.stdout + result.stderr
    marker = "PASS PARI Mordell-Weil:"
    if result.returncode != 0 or marker not in output:
        raise RuntimeError("PARI Mordell-Weil replay failed:\n" + output)
    for line in result.stdout.splitlines():
        if line.startswith(marker):
            print(line)

# Boundary-unit pullback and derivative vectors in the ordered basis
# (p1,p2,p3,q1,q2,q3).
pullback_zero = sp.Matrix([3, 3, 1, 0, 0, 0])
pullback_one = sp.Matrix([0, 0, 0, 2, 2, 1])
derivative_vector = sp.Matrix([2, 2, 0, 1, 1, 0])
target_lattice = sp.Matrix.hstack(pullback_zero, pullback_one)
assert target_lattice.rank() == 2
assert sp.Matrix.hstack(target_lattice, derivative_vector).rank() == 3
assert len(boundary_factors) - target_lattice.rank() == 4


# ---------------------------------------------------------------------------
# PSL_2(F_11) as determinant-one matrices modulo the scalar pair {+/-I}.

prime = 11


def canonical(matrix):
    matrix = tuple(entry % prime for entry in matrix)
    negative = tuple((-entry) % prime for entry in matrix)
    return min(matrix, negative)


def multiply(left, right):
    return canonical(
        (
            left[0] * right[0] + left[1] * right[2],
            left[0] * right[1] + left[1] * right[3],
            left[2] * right[0] + left[3] * right[2],
            left[2] * right[1] + left[3] * right[3],
        )
    )


identity = canonical((1, 0, 0, 1))
group = {
    canonical(matrix)
    for matrix in product(range(prime), repeat=4)
    if (matrix[0] * matrix[3] - matrix[1] * matrix[2]) % prime == 1
}
assert len(group) == 660


def inverse(matrix):
    return canonical((matrix[3], -matrix[1], -matrix[2], matrix[0]))


def element_order(matrix):
    power = identity
    for order in range(1, 67):
        power = multiply(power, matrix)
        if power == identity:
            return order
    raise AssertionError("element order exceeded the PSL_2(11) exponent bound")


elements_by_order = {}
for element in group:
    elements_by_order.setdefault(element_order(element), []).append(element)
assert {order: len(elements) for order, elements in elements_by_order.items()} == {
    1: 1,
    2: 55,
    3: 110,
    5: 264,
    6: 110,
    11: 120,
}


def generated_subgroup(generators):
    subgroup = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = multiply(generator, current)
            if candidate not in subgroup:
                subgroup.add(candidate)
                frontier.append(candidate)
    return frozenset(subgroup)


def conjugate_subgroup(subgroup, element):
    element_inverse = inverse(element)
    return frozenset(
        multiply(multiply(element, member), element_inverse) for member in subgroup
    )


# A_5 is the finite (2,3,5) triangle group.  The order-sixty subgroups found
# from such pairs are therefore A_5 subgroups without relying on a group table.
a5_subgroups = set()
for involution in elements_by_order[2]:
    for order_three in elements_by_order[3]:
        if element_order(multiply(involution, order_three)) != 5:
            continue
        subgroup = generated_subgroup((involution, order_three))
        if len(subgroup) == 60:
            a5_subgroups.add(subgroup)
assert len(a5_subgroups) == 22
for subgroup in a5_subgroups:
    assert Counter(element_order(element) for element in subgroup) == Counter(
        {1: 1, 2: 15, 3: 20, 5: 24}
    )

unused_subgroups = set(a5_subgroups)
a5_conjugacy_orbits = []
while unused_subgroups:
    representative = next(iter(unused_subgroups))
    orbit = {conjugate_subgroup(representative, element) for element in group}
    assert orbit <= a5_subgroups
    a5_conjugacy_orbits.append(orbit)
    unused_subgroups -= orbit
assert sorted(len(orbit) for orbit in a5_conjugacy_orbits) == [11, 11]

a5_first = next(iter(a5_conjugacy_orbits[0]))
a5_second = next(iter(a5_conjugacy_orbits[1]))


def left_coset_action(subgroup):
    unseen = set(group)
    cosets = []
    coset_index = {}
    while unseen:
        representative = next(iter(unseen))
        coset = {multiply(representative, member) for member in subgroup}
        index = len(cosets)
        cosets.append(coset)
        for member in coset:
            coset_index[member] = index
        unseen -= coset

    def action(element):
        return tuple(
            coset_index[multiply(element, next(iter(coset)))] for coset in cosets
        )

    return cosets, action


first_cosets, first_action = left_coset_action(a5_first)
second_cosets, second_action = left_coset_action(a5_second)
assert len(first_cosets) == len(second_cosets) == 11


def natural_action(matrix):
    a, b, c, d = matrix
    images = []
    for point in range(prime):
        denominator = (c * point + d) % prime
        if denominator == 0:
            images.append(prime)
        else:
            images.append(
                ((a * point + b) * pow(denominator, -1, prime)) % prime
            )
    images.append(prime if c == 0 else (a * pow(c, -1, prime)) % prime)
    return tuple(images)


def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


# Gassmann equivalence is checked in its stronger all-powers form: the two
# actions have the same cycle partition for every element of the common group.
for element in group:
    assert cycle_type(first_action(element)) == cycle_type(second_action(element))

# The actions are nevertheless nonisomorphic because their stabilizers lie in
# different conjugacy classes.  Their cross-action subdegrees are 5 and 6.
assert all(
    conjugate_subgroup(a5_first, element) != a5_second for element in group
)
unseen_points = set(range(11))
cross_orbits = []
while unseen_points:
    point = next(iter(unseen_points))
    representative = next(iter(second_cosets[point]))
    orbit = {
        next(
            index
            for index, coset in enumerate(second_cosets)
            if multiply(element, representative) in coset
        )
        for element in a5_first
    }
    cross_orbits.append(orbit)
    unseen_points -= orbit
assert sorted(len(orbit) for orbit in cross_orbits) == [5, 6]


# Conjugacy classes and the two rigid (2,3,11) Nielsen orbits.
unseen_elements = set(group)
conjugacy_classes = []
while unseen_elements:
    representative = next(iter(unseen_elements))
    conjugacy_class = {
        multiply(multiply(element, representative), inverse(element))
        for element in group
    }
    conjugacy_classes.append(conjugacy_class)
    unseen_elements -= conjugacy_class
assert sorted(
    (element_order(next(iter(conjugacy_class))), len(conjugacy_class))
    for conjugacy_class in conjugacy_classes
) == [(1, 1), (2, 55), (3, 110), (5, 132), (5, 132), (6, 110), (11, 60), (11, 60)]

generating_triples = []
for involution in elements_by_order[2]:
    for order_three in elements_by_order[3]:
        product_element = multiply(involution, order_three)
        if element_order(product_element) != 11:
            continue
        if len(generated_subgroup((involution, order_three))) != len(group):
            continue
        generating_triples.append(
            (involution, order_three, inverse(product_element))
        )
assert len(generating_triples) == 1320

unused_triples = set(generating_triples)
triple_orbits = []
while unused_triples:
    triple = next(iter(unused_triples))
    orbit = {
        tuple(
            multiply(multiply(element, entry), inverse(element))
            for entry in triple
        )
        for element in group
    }
    triple_orbits.append(orbit)
    unused_triples -= orbit
assert sorted(len(orbit) for orbit in triple_orbits) == [660, 660]

order_eleven_classes = [
    conjugacy_class
    for conjugacy_class in conjugacy_classes
    if element_order(next(iter(conjugacy_class))) == 11
]
assert len(order_eleven_classes) == 2
assert {
    next(iter(orbit))[2] in order_eleven_classes[0] for orbit in triple_orbits
} == {False, True}

involution, order_three, order_eleven = next(iter(triple_orbits[0]))
natural_passport = (
    cycle_type(natural_action(order_three)),
    cycle_type(natural_action(involution)),
    cycle_type(natural_action(order_eleven)),
)
exceptional_passport = (
    cycle_type(first_action(order_three)),
    cycle_type(first_action(involution)),
    cycle_type(first_action(order_eleven)),
)
assert natural_passport == (
    (3, 3, 3, 3),
    (2, 2, 2, 2, 2, 2),
    (11, 1),
)
assert exceptional_passport == (
    (3, 3, 3, 1, 1),
    (2, 2, 2, 2, 1, 1, 1),
    (11,),
)


def permutation_cycles(permutation):
    unseen = set(range(len(permutation)))
    cycles = []
    while unseen:
        start = min(unseen)
        cycle = []
        current = start
        while current not in cycle:
            cycle.append(current)
            unseen.discard(current)
            current = permutation[current]
        cycles.append(tuple(cycle))
    return cycles


def boundary_projection_profile(source_permutation, target_permutation, projection):
    """Ramification partitions above every target inertia cycle."""

    source_cycles = permutation_cycles(source_permutation)
    target_cycles = permutation_cycles(target_permutation)
    target_cycle_index = {
        point: index
        for index, target_cycle in enumerate(target_cycles)
        for point in target_cycle
    }
    ramification = [[] for _ in target_cycles]
    for source_cycle in source_cycles:
        target_indices = {
            target_cycle_index[projection[source_point]]
            for source_point in source_cycle
        }
        assert len(target_indices) == 1
        target_index = next(iter(target_indices))
        target_length = len(target_cycles[target_index])
        assert len(source_cycle) % target_length == 0
        ramification[target_index].append(len(source_cycle) // target_length)
    return tuple(
        sorted(
            (len(target_cycles[index]), tuple(sorted(indices)))
            for index, indices in enumerate(ramification)
        )
    )


first_coset_lookup = {
    member: index for index, coset in enumerate(first_cosets) for member in coset
}
second_coset_lookup = {
    member: index for index, coset in enumerate(second_cosets) for member in coset
}
normalization_components = {}
for cross_orbit in cross_orbits:
    point = min(cross_orbit)
    cross_representative = next(iter(second_cosets[point]))
    intersection_subgroup = a5_first.intersection(
        conjugate_subgroup(a5_second, cross_representative)
    )
    component_cosets, component_action = left_coset_action(intersection_subgroup)
    projection_to_first = []
    projection_to_second = []
    for component_coset in component_cosets:
        representative = next(iter(component_coset))
        projection_to_first.append(first_coset_lookup[representative])
        projection_to_second.append(
            second_coset_lookup[multiply(representative, cross_representative)]
        )

    component_degree = len(component_cosets)
    component_passport = (
        cycle_type(component_action(order_three)),
        cycle_type(component_action(involution)),
        cycle_type(component_action(order_eleven)),
    )
    normalization_components[len(cross_orbit)] = {
        "subgroup": intersection_subgroup,
        "degree": component_degree,
        "action": component_action,
        "passport": component_passport,
        "projection_to_first": tuple(projection_to_first),
        "projection_to_second": tuple(projection_to_second),
    }

assert set(normalization_components) == {5, 6}
component_five = normalization_components[5]
component_six = normalization_components[6]
assert len(component_five["subgroup"]) == 12
assert Counter(element_order(element) for element in component_five["subgroup"]) == Counter(
    {1: 1, 2: 3, 3: 8}
)
assert len(component_six["subgroup"]) == 10
assert Counter(element_order(element) for element in component_six["subgroup"]) == Counter(
    {1: 1, 2: 5, 5: 4}
)
assert component_five["degree"] == 55
assert component_six["degree"] == 66
assert tuple(Counter(passport_row) for passport_row in component_five["passport"]) == (
    Counter({3: 17, 1: 4}),
    Counter({2: 26, 1: 3}),
    Counter({11: 5}),
)
assert tuple(Counter(passport_row) for passport_row in component_six["passport"]) == (
    Counter({3: 22}),
    Counter({2: 30, 1: 6}),
    Counter({11: 6}),
)

expected_projection_profiles = {
    5: (
        (
            (1, (1, 1, 3)),
            (1, (1, 1, 3)),
            (3, (1, 1, 1, 1, 1)),
            (3, (1, 1, 1, 1, 1)),
            (3, (1, 1, 1, 1, 1)),
        ),
        (
            (1, (1, 2, 2)),
            (1, (1, 2, 2)),
            (1, (1, 2, 2)),
            (2, (1, 1, 1, 1, 1)),
            (2, (1, 1, 1, 1, 1)),
            (2, (1, 1, 1, 1, 1)),
            (2, (1, 1, 1, 1, 1)),
        ),
        ((11, (1, 1, 1, 1, 1)),),
    ),
    6: (
        (
            (1, (3, 3)),
            (1, (3, 3)),
            (3, (1, 1, 1, 1, 1, 1)),
            (3, (1, 1, 1, 1, 1, 1)),
            (3, (1, 1, 1, 1, 1, 1)),
        ),
        (
            (1, (1, 1, 2, 2)),
            (1, (1, 1, 2, 2)),
            (1, (1, 1, 2, 2)),
            (2, (1, 1, 1, 1, 1, 1)),
            (2, (1, 1, 1, 1, 1, 1)),
            (2, (1, 1, 1, 1, 1, 1)),
            (2, (1, 1, 1, 1, 1, 1)),
        ),
        ((11, (1, 1, 1, 1, 1, 1)),),
    ),
}
for projection_degree, component in normalization_components.items():
    source_actions = tuple(
        component["action"](element)
        for element in (order_three, involution, order_eleven)
    )
    for target_action, projection in (
        (first_action, component["projection_to_first"]),
        (second_action, component["projection_to_second"]),
    ):
        profiles = tuple(
            boundary_projection_profile(
                source_action,
                target_action(element),
                projection,
            )
            for source_action, element in zip(
                source_actions,
                (order_three, involution, order_eleven),
                strict=True,
            )
        )
        assert profiles == expected_projection_profiles[projection_degree]


def riemann_hurwitz_genus(degree, passport):
    total_index = sum(degree - len(cycles) for cycles in passport)
    assert total_index % 2 == 0
    return 1 - degree + total_index // 2, total_index


natural_genus, natural_index = riemann_hurwitz_genus(12, natural_passport)
exceptional_genus, exceptional_index = riemann_hurwitz_genus(
    11, exceptional_passport
)
assert (natural_genus, natural_index) == (1, 24)
assert (exceptional_genus, exceptional_index) == (0, 20)

component_five_genus, component_five_index = riemann_hurwitz_genus(
    component_five["degree"], component_five["passport"]
)
component_six_genus, component_six_index = riemann_hurwitz_genus(
    component_six["degree"], component_six["passport"]
)
assert (component_five_genus, component_five_index) == (1, 110)
assert (component_six_genus, component_six_index) == (2, 134)

# The total-degree-five/six plane closures have arithmetic genera six and ten.
# Their exact affine node counts five and eight (the optional Singular replay)
# therefore give the same normalization genera as the group calculation.
assert (5 - 1) * (5 - 2) // 2 - 5 == component_five_genus
assert (6 - 1) * (6 - 2) // 2 - 8 == component_six_genus

print("PASS corrected degree-eleven Shabat factorization and derivative")
print("PASS discriminant square and rank-six/rank-two boundary-unit ledger")
print("PASS PSL_2(11): two nonconjugate Gassmann A_5-coset actions")
print("PASS exact degree-5/6 correspondence and cross-action subdegrees 5+6")
print("PASS two rigid (2,3,11) Nielsen orbits")
print("PASS Riemann-Hurwitz: exceptional degree 11 has genus 0")
print("PASS Riemann-Hurwitz: natural degree 12 has genus 1")
print("PASS normalized C5/A4 and C6/D10 quotients have genera 1 and 2")
print("PASS normalized C5 has a K-point and elliptic j-invariant -121")
print("PASS C5 trace 2 versus X_0(11) trace -1 above 23: not isogenous")
print("PASS normalized C6 has an even hyperelliptic model and split Jacobian")
print("PASS Jac(C6) is (2,2)-isogenous to the conductor-11 and CM-121 curves")
print("PASS exact positive-genus boundary-unit lattices have ranks 3, 14, and 17")
print("PASS both correspondence pullback images have rank 10 and primitive cokernel")
print("PASS two-output derivative ledger is primitive but not a full character completion")
print("PASS residual swap modules are Z^4 and Z^5 plus one exchanged pair")
print("PASS effective mask bases force minimal normal-support degrees 2 and 3")
print("PASS C6 infinity imbalance descends with values (-2,0,0,0,0,1,1)")
print("PASS exact X_0(11) model and degree-12 j-map passport")

if __name__ == "__main__":
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument(
        "--singular-normalization",
        action="store_true",
        help=(
            "replay affine normalization, conductor, node support, the C6 "
            "canonical elimination, the C5 Cremona reduction, and explicit "
            "normalization-module representatives for every C5/C6 mask"
        ),
    )
    argument_parser.add_argument(
        "--pari-mordell-weil",
        action="store_true",
        help=(
            "certify the rational ranks and Heegner generators whose "
            "quadratic rank decompositions give the K-Mordell-Weil ranks"
        ),
    )
    arguments = argument_parser.parse_args()
    if arguments.singular_normalization:
        run_singular_normalization_check()
    if arguments.pari_mordell_weil:
        run_pari_mordell_weil_check()
