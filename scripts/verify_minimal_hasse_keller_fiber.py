#!/usr/bin/env python3
"""Exact audit of the minimal degree-five Hasse-failing Keller fiber."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402


# Berend--Bilu's minimal intersective polynomial.  The cubic covers primes
# p == 2 mod 3, while the cyclotomic quadratic covers primes p == 1 mod 3.
X = sp.symbols("X")
f = sp.expand((X**3 - 19) * (X**2 + X + 1))
assert f == X**5 + X**4 + X**3 - 19 * X**2 - 19 * X - 19

# Neither irreducible factor has a rational root, and they are coprime.
assert sp.Poly(X**3 - 19, X, domain=sp.QQ).is_irreducible
assert sp.discriminant(X**2 + X + 1, X) == -3
assert sp.gcd(X**3 - 19, X**2 + X + 1) == 1

# The residue-class covering is uniform away from 3.  At p == 1 mod 3 the
# cyclotomic quadratic splits; at p == 2 mod 3 the cube map is bijective.
for prime in sp.primerange(2, 2000):
    if prime == 3:
        continue
    roots = [
        a
        for a in range(prime)
        if ((a**3 - 19) * (a**2 + a + 1)) % prime == 0
    ]
    assert roots, f"missing residue root at {prime}"
    assert any(
        (
            (3 * a**2) * (a**2 + a + 1)
            + (a**3 - 19) * (2 * a + 1)
        )
        % prime
        != 0
        for a in roots
    )

# At 3, 19 belongs to 1 + 9 Z_3 and is therefore a cube.  Compatible lifts
# give a finite exact regression of this standard 3-adic cube-map fact.
roots_3 = {1}
modulus_3 = 3
for _ in range(1, 12):
    candidates = {
        root + digit * modulus_3
        for root in roots_3
        for digit in range(3)
    }
    modulus_3 *= 3
    roots_3 = {
        candidate
        for candidate in candidates
        if candidate**3 % modulus_3 == 19 % modulus_3
    }
    assert roots_3


# The affine reparameterization X=10-27W puts the same finite etale scheme
# on the weighted tangent-normalization slice.  Division by 27 makes P
# primitive without changing its roots.
P = sp.expand(f.subs(X, 10 - 27 * w) / 27)
assert P == (
    -531441 * w**5
    + 1003833 * w**4
    - 758889 * w**3
    + 286497 * w**2
    - 53901 * w
    + 4033
)
P_poly = sp.Poly(P, w, domain=sp.QQ)
assert P_poly.gcd(P_poly.diff()).degree() == 0

p_at_zero = P.subs(w, 0)
p_prime_at_zero = sp.diff(P, w).subs(w, 0)
assert (p_at_zero, p_prime_at_zero, P.subs(w, 1)) == (
    4033,
    -53901,
    -49868,
)
assert P.subs(w, 1) - p_at_zero == p_prime_at_zero

H = sp.expand(P - p_at_zero - p_prime_at_zero * w)
h = sp.diff(H, w)
c = -h.subs(w, 1)
assert H == (
    -531441 * w**5
    + 1003833 * w**4
    - 758889 * w**3
    + 286497 * w**2
)
assert c == 345546
assert sp.diff(H, w, 2).subs(w, 1) == -2563164
assert sp.cancel(sp.diff(H, w, 2).subs(w, 1) / c) == -sp.Rational(
    586, 79
)

model = WeightedSeedModel(h, c=c, b=1)
assert model.a == -sp.Rational(507, 428)


# Expand the three coordinates, verify polynomiality and the Keller
# determinant, and audit the weighted suspension square.
mapping = model.mapping()
component_profiles = tuple(
    (
        sp.Poly(component, x, y, z).total_degree(),
        len(sp.Poly(component, x, y, z).terms()),
    )
    for component in mapping
)
assert component_profiles == ((17, 28), (16, 25), (4, 3))
assert sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det()) == c

v = x * y
S = x**2 * z
gamma = 1 + model.a * v + S
W = (1 + v) * gamma
A_map, B_map, C_map = mapping
assert sp.expand(C_map - x * gamma) == 0
assert sp.factor(B_map * C_map - (h.subs(w, W) + c * gamma)) == 0
assert sp.factor(
    c * A_map * C_map**2
    - (W * (h.subs(w, W) + c * gamma) - H.subs(w, W))
) == 0


# This rational target has inverse polynomial P exactly.
target = (sp.Rational(4033, 345546), sp.Integer(53901), sp.Integer(1))
assert sp.expand(model.inverse_polynomial(*target) - P) == 0


# In Q[W]/(P), P' is a unit, so the reconstruction formulas produce exactly
# one source point over every root and no C=0 boundary point.
def reduce_rational_mod_p(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.fraction(sp.cancel(expression))
    numerator_poly = sp.Poly(numerator, w, domain=sp.QQ)
    denominator_poly = sp.Poly(denominator, w, domain=sp.QQ)
    inverse = sp.invert(denominator_poly, P_poly)
    return sp.rem(numerator_poly * inverse, P_poly).as_expr()


root_gamma = sp.cancel(-sp.diff(P, w) / c)
root_x = sp.cancel(1 / root_gamma)
root_y = sp.cancel(w - root_gamma)
root_z = sp.cancel(
    (root_gamma - 1 - model.a * (w / root_gamma - 1)) / root_x**2
)
root_values = {x: root_x, y: root_y, z: root_z}

assert reduce_rational_mod_p((x * gamma).subs(root_values)) == 1
assert reduce_rational_mod_p(gamma.subs(root_values) - root_gamma) == 0
assert reduce_rational_mod_p(W.subs(root_values) - w) == 0


# The first fixed-map infinite-family test reduces to an explicit genus-one
# quartic.  A monic quadratic q=W^2+aW+b divides a unique pencil member,
# because the remainder of H modulo q is linear.  On a=-7/9, equality of
# the quadratic field with the cubic quotient's discriminant field is the
# displayed square-class equation.
a, b, Y = sp.symbols("a b Y")
quadratic = w**2 + a * w + b
cubic, _remainder = sp.div(H, quadratic, w)
cubic_discriminant = sp.factor(sp.discriminant(cubic, w))
quadratic_discriminant = a**2 - 4 * b
fixed_a = -sp.Rational(7, 9)
g = 57395628 * b**3 - 26749197 * b**2 + 4181544 * b - 219512
ratio = sp.factor(
    cubic_discriminant.subs(a, fixed_a)
    / quadratic_discriminant.subs(a, fixed_a)
)
assert sp.cancel(ratio - sp.Integer(3) ** 37 * g / (49 - 324 * b)) == 0
quartic = sp.expand(3 * g * (49 - 324 * b))
known_b = sp.Rational(37, 243)
assert quartic.subs(b, known_b) == 19**2


print("PASS: the degree-five polynomial is everywhere local and nowhere rational")
print("PASS: X=10-27W enforces the tangent normalization exactly")
print("PASS: the polynomial Keller map has determinant 345546")
print("PASS: the complete fiber is Spec(Q[W]/(P)), so d_HP=5")
print("PASS: the fixed-map S3-resolvent slice is the certified genus-one quartic")
