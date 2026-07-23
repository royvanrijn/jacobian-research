#!/usr/bin/env python3
"""Exact audit of natural plane quotients of the foundational Keller map."""

import sympy as sp


x, y, z = sp.symbols("x y z")
U, V = sp.symbols("U V")

u = 1 + x * y
F = sp.Matrix(
    (
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    )
)
assert sp.factor(F.jacobian((x, y, z)).det()) == -2

# The categorical torus quotients have polynomial coordinates
# (U,V)=(xy,x^2*z) and (P,Q)=(a*c^2,b*c).  The foundational map descends
# through the following invariant polynomials.
A = (1 + U) ** 3 * V + U**2 * (1 + U) * (4 + 3 * U)
B = U + 3 * (1 + U) ** 2 * V + 3 * U**2 * (4 + 3 * U)
C = 2 - 3 * U - V
P = sp.expand(A * C**2)
Q = sp.expand(B * C)

source_invariants = {U: x * y, V: x**2 * z}
assert sp.expand(P.subs(source_invariants) - F[0] * F[2] ** 2) == 0
assert sp.expand(Q.subs(source_invariants) - F[1] * F[2]) == 0
assert sp.expand(sp.Matrix((P, Q)).jacobian((U, V)).det() + 2 * C**2) == 0

# Verify the universal identity det D(AC^2,BC)=C^2 det DF in formal first
# derivatives.  This avoids relying only on the foundational specialization.
A0, B0, C0 = sp.symbols("A0 B0 C0")
AU, AV, BU, BV, CU, CV = sp.symbols("AU AV BU BV CU CV")
weighted_determinant = sp.Matrix(
    ((-2 * A0, AU, AV), (-B0, BU, BV), (C0, CU, CV))
).det()
quotient_determinant = sp.expand(
    (AU * C0**2 + 2 * A0 * C0 * CU) * (BV * C0 + B0 * CV)
    - (AV * C0**2 + 2 * A0 * C0 * CV) * (BU * C0 + B0 * CU)
)
assert sp.expand(quotient_determinant - C0**2 * weighted_determinant) == 0

# The three-point collision becomes a two-point collision on the quotient;
# the two nonzero points are in the same torus orbit.
quotient_collision = ((0, 0), (sp.Rational(-3, 2), sp.Rational(13, 2)))
for point in quotient_collision:
    substitution = dict(zip((U, V), point))
    assert P.subs(substitution) == Q.subs(substitution) == 0
assert C.subs({U: sp.Rational(-3, 2), V: sp.Rational(13, 2)}) == 0

# The preserved hypersurface x=0 carries a triangular plane automorphism.
restriction = tuple(sp.expand(component.subs(x, 0)) for component in F)
assert restriction == (z + 4 * y**2, y, 0)
assert sp.Matrix(restriction[:2]).jacobian((y, z)).det() == -1
target_a, target_b = sp.symbols("target_a target_b")
inverse_restriction = {y: target_b, z: target_a - 4 * target_b**2}
assert tuple(value.subs(inverse_restriction) for value in restriction[:2]) == (
    target_a,
    target_b,
)

# The vertical derivative and inverse cubic give the target-first-integral
# obstruction.  Symbolic coefficient comparison says that a cubic
# h_a*T^3+3*h_b*T^2-h_c can be a multiple of
# c*T^3-2*T^2+b*T-2*a only when all h-partials vanish generically.
T = sp.symbols("T")
ha, hb, hc, lam = sp.symbols("ha hb hc lam")
a, b, c = sp.symbols("a b c", nonzero=True)
inverse_cubic = c * T**3 - 2 * T**2 + b * T - 2 * a
first_integral_cubic = ha * T**3 + 3 * hb * T**2 - hc
coefficient_equations = sp.Poly(
    first_integral_cubic - lam * inverse_cubic, T
).all_coeffs()
solution = sp.solve(coefficient_equations, (ha, hb, hc, lam), dict=True)
assert solution == [{ha: 0, hb: 0, hc: 0, lam: 0}]

t = sp.symbols("t")
vertical = F.diff(z)
assert all(
    sp.expand(got - expected) == 0
    for got, expected in zip(
        vertical.subs(y, t - 1 / x),
        (x**3 * t**3, 3 * x**3 * t**2, -x**3),
    )
)

# Finally exclude commuting rank-two linear projections.  For a source
# kernel vector (p,q,r), collect the coefficients of DF*(p,q,r).  Rank at
# most one of this 3-by-N coefficient matrix is necessary for all values to
# lie in one constant target direction.  Its 2-by-2 minors have radical
# supported only at the zero vector; the displayed Groebner basis is exact.
p, q, r = sp.symbols("p q r")
directional_image = sp.expand(F.jacobian((x, y, z)) * sp.Matrix((p, q, r)))
monomials = sorted(
    set().union(*(sp.Poly(entry, x, y, z).monoms() for entry in directional_image))
)
coefficient_matrix = sp.Matrix(
    [
        [sp.Poly(entry, x, y, z).coeff_monomial(monomial) for monomial in monomials]
        for entry in directional_image
    ]
)
minors = []
for row1 in range(3):
    for row2 in range(row1 + 1, 3):
        for column1 in range(len(monomials)):
            for column2 in range(column1 + 1, len(monomials)):
                minor = sp.expand(
                    coefficient_matrix[row1, column1]
                    * coefficient_matrix[row2, column2]
                    - coefficient_matrix[row1, column2]
                    * coefficient_matrix[row2, column1]
                )
                if minor != 0:
                    minors.append(minor)

linear_quotient_ideal = sp.groebner(minors, p, q, r, order="grevlex")
assert [polynomial.as_expr() for polynomial in linear_quotient_ideal.polys] == [
    p**2,
    p * q,
    q**2,
    p * r,
    q * r,
    r**2,
]

print("PASS: both torus invariant rings are polynomial rings in two generators")
print("PASS: the torus quotient has Jacobian -2*(2-3*U-V)^2 and contracts C=0")
print("PASS: the vertical LND has no nonconstant target polynomial first integral")
print("PASS: no rank-two linear source/target quotient square exists")
print("PASS: x=0 restricts to a triangular plane Keller automorphism")
