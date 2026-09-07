#!/usr/bin/env python3
"""Exact audit of the everywhere-local, nowhere-rational Keller fiber."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402


# The first three factors form an elementary intersective polynomial.  The
# last factor enforces the weighted tangent normalization without adding a
# rational or real root.
Q = sp.expand((w**2 - 2) * (w**2 - 17) * (w**2 - 34))
R = 157 * w**2 - 267 * w + 399
P = sp.expand(Q * R)

assert Q == w**6 - 53 * w**4 + 680 * w**2 - 1156
assert P == (
    157 * w**8
    - 267 * w**7
    - 7922 * w**6
    + 14151 * w**5
    + 85613 * w**4
    - 181560 * w**3
    + 89828 * w**2
    + 308652 * w
    - 461244
)


# No factor has a rational root.  A rational square root of an integer is an
# integer, while the final quadratic is positive definite.
for radicand in (2, 17, 34):
    assert math.isqrt(radicand) ** 2 != radicand
assert sp.discriminant(R, w) == -179283 < 0

# The eight geometric roots are distinct.  In particular every local root
# used below is a regular reconstruction point.
P_poly = sp.Poly(P, w, domain=sp.QQ)
assert P_poly.gcd(P_poly.diff()).degree() == 0


# Local covering away from 2 and 17: if the square classes of 2 and 17 are
# both nontrivial, their product 34 is trivial.  This four-case truth table is
# the only uniform residue argument needed.
for square_class_2 in (-1, 1):
    for square_class_17 in (-1, 1):
        assert 1 in (
            square_class_2,
            square_class_17,
            square_class_2 * square_class_17,
        )

# The exceptional completions are covered by sqrt(17) in Q_2 and sqrt(2) in
# Q_17.  The congruences are the standard Hensel starting data.
assert 17 % 8 == 1
assert 6**2 % 17 == 2


def hensel_square_root(radicand: int, prime: int, exponent: int) -> int:
    """Lift one square root through p^exponent by exhaustive Hensel steps."""
    roots = [
        residue
        for residue in range(prime)
        if residue**2 % prime == radicand % prime
    ]
    assert roots
    root = roots[0]
    modulus = prime
    for _ in range(1, exponent):
        candidates = tuple(root + digit * modulus for digit in range(prime))
        modulus *= prime
        root = next(
            candidate
            for candidate in candidates
            if candidate**2 % modulus == radicand % modulus
        )
    return root


# At 2 the usual derivative form of Hensel's lemma is replaced by the unit
# criterion u == 1 mod 8; explicit compatible lifts give a finite regression.
roots_2 = {1}
modulus_2 = 2
for _ in range(1, 16):
    candidates = {
        root + digit * modulus_2
        for root in roots_2
        for digit in (0, 1)
    }
    modulus_2 *= 2
    roots_2 = {
        candidate
        for candidate in candidates
        if candidate**2 % modulus_2 == 17 % modulus_2
    }
    assert roots_2
assert all(root**2 % modulus_2 == 17 % modulus_2 for root in roots_2)

root_17 = hensel_square_root(2, 17, 5)
assert root_17**2 % 17**5 == 2

# Direct prime-field checks through a substantial finite range guard the
# square-class proof and its exceptional-prime split against transcription
# errors.  Every selected odd root is simple and therefore lifts to Q_p.
for prime in sp.primerange(3, 2000):
    radicands = (2,) if prime == 17 else (2, 17, 34)
    witnesses = [
        residue
        for radicand in radicands
        for residue in range(prime)
        if residue * residue % prime == radicand % prime
        and 2 * residue % prime != 0
    ]
    assert witnesses, f"missing simple local root at {prime}"


# Tangent normalization.  The identity P(1)-P(0)=P'(0) is exactly what is
# needed after removing the affine part of P.
p_at_zero = P.subs(w, 0)
p_prime_at_zero = sp.diff(P, w).subs(w, 0)
assert (p_at_zero, p_prime_at_zero, P.subs(w, 1)) == (
    -461244,
    308652,
    -152592,
)
assert P.subs(w, 1) - p_at_zero == p_prime_at_zero

H = sp.expand(P - p_at_zero - p_prime_at_zero * w)
assert H == (
    157 * w**8
    - 267 * w**7
    - 7922 * w**6
    + 14151 * w**5
    + 85613 * w**4
    - 181560 * w**3
    + 89828 * w**2
)

h = sp.diff(H, w)
c = -h.subs(w, 1)
assert (H.subs(w, 0), h.subs(w, 0), H.subs(w, 1), h.subs(w, 1)) == (
    0,
    0,
    0,
    38,
)
assert c == -38
assert sp.diff(H, w, 2).subs(w, 1) == 160590
assert sp.cancel(sp.diff(H, w, 2).subs(w, 1) / c) == -sp.Rational(80295, 19)

model = WeightedSeedModel(h, c=c, b=1)
assert model.a == -sp.Rational(80276, 80257)


# Construct and expand the three polynomial coordinates.  This independently
# exercises all exact divisions in the weighted formulas before checking the
# constant determinant.
mapping = model.mapping()
component_profiles = tuple(
    (
        sp.Poly(component, x, y, z).total_degree(),
        len(sp.Poly(component, x, y, z).terms()),
    )
    for component in mapping
)
assert component_profiles == ((32, 82), (31, 76), (4, 3))
assert sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det()) == -38


# Audit the suspension square directly for this seed.  These identities imply
# that every point above C=1 is represented by a root of the inverse pencil.
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


# The integral target has inverse polynomial P exactly.
target = (sp.Integer(12138), sp.Integer(-308652), sp.Integer(1))
assert sp.expand(model.inverse_polynomial(*target) - P) == 0


# In Q[W]/(P), P' is a unit.  The reconstruction formulas therefore produce
# one source point over every root field.  The following quotient-ring audit
# checks W=(1+xy)gamma and the definitions of gamma and C exactly.
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


print("PASS: P has no rational root and has a root over R and every Q_p")
print("PASS: the normalized degree-eight seed is weighted-admissible")
print("PASS: the explicit polynomial map has determinant -38")
print("PASS: the complete C=1 fiber is Spec(Q[W]/(P))")
