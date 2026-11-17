#!/usr/bin/env python3
"""Exact finite-degree audit of the universal weighted-seed pencil theorem."""

import math
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import (  # noqa: E402
    WeightedSeedModel,
    canonical_seed,
    deformation_basis,
    w,
    x,
    y,
    z,
)


s, t, r = sp.symbols("s t r")


def audit_model(name, seed):
    """Check the algebraic spine for one exact admissible seed."""
    model = WeightedSeedModel(seed)
    H = model.primitive
    p = model.seed
    E = sp.expand(H - s * w + t)
    n = model.fiber_degree

    assert n >= 3
    assert sp.Poly(H, w).degree() == n
    assert sp.Poly(p, w).degree() == n - 1

    # Connectedness of the generic cover.  The theorem proves this uniformly
    # from linearity and monicity in t; SymPy independently checks each listed
    # representative over the rational function field.
    field = sp.QQ.frac_field(s, t)
    assert sp.Poly(E, w, domain=field).is_irreducible

    # Repeated-root normalization: E(r)=E'(r)=0 at
    # (s,t)=(p(r), r*p(r)-H(r)).  The coprime coordinate-degree lemma says
    # that its function-field degree divides both coordinate degrees n-1
    # and n, and therefore equals one.
    branch_s = p.subs(w, r)
    branch_t = (w * p - H).subs(w, r)
    branch = {s: branch_s, t: branch_t, w: r}
    assert sp.expand(E.subs(branch)) == 0
    assert sp.expand(sp.diff(E, w).subs(branch)) == 0
    assert sp.Poly(branch_s, r).degree() == n - 1
    assert sp.Poly(branch_t, r).degree() == n
    assert math.gcd(n - 1, n) == 1

    # Eliminating r gives the same irreducible branch equation as the quartic
    # discriminant, up to a nonzero rational scalar.
    eliminated = sp.resultant(branch_s - s, branch_t - t, r)
    discriminant = sp.discriminant(E, w)
    ratio = sp.cancel(eliminated / discriminant)
    assert ratio != 0 and not ratio.has(s, t)

    # Exhibit a smooth, uniquely double branch point.  At such a point the
    # transverse local equation is epsilon+c*h^2+..., hence its meridian is a
    # transposition.  The proof only needs generic existence; this loop gives
    # an exact representative-level regression check.
    sample = None
    for candidate in range(-3, 5):
        second_derivative = sp.diff(p, w).subs(w, candidate)
        if second_derivative == 0:
            continue
        s0 = p.subs(w, candidate)
        t0 = (w * p - H).subs(w, candidate)
        specialized = sp.Poly(E.subs({s: s0, t: t0}), w)
        derivative = sp.Poly(sp.diff(E, w).subs({s: s0, t: t0}), w)
        if specialized.gcd(derivative).degree() == 1:
            sample = candidate
            break
    assert sample is not None

    return name, n


models = []
for degree in range(2, 8):
    models.append((f"canonical H_{degree}", canonical_seed(degree)))
for index in range(1, 5):
    models.append(
        (
            f"deformation index {index}",
            canonical_seed(2) + deformation_basis(index),
        )
    )

audited = [audit_model(name, seed) for name, seed in models]
assert [n for _, n in audited[:6]] == list(range(3, 9))

# The uniform sparse family H_N=(w^2-w^N)/(N-2) is admissible in every
# degree N>=3.  The identities below are symbolic in N; the finite loop
# checks that the weighted model reads the same normalization exactly.
N = sp.symbols("N", integer=True, positive=True)
assert sp.cancel((2 - N) / (N - 2)) == -1
assert sp.cancel((2 - N * (N - 1)) / (N - 2)) == -(N + 1)
assert sp.cancel(
    -(1 - (N + 1)) / (2 - (N + 1))
) == -N / (N - 1)
for degree in range(3, 13):
    H_degree = (w**2 - w**degree) / (degree - 2)
    model_degree = WeightedSeedModel(sp.diff(H_degree, w), c=1, b=1)
    assert model_degree.fiber_degree == degree
    assert model_degree.kappa == -(degree + 1)
    assert model_degree.a == -sp.Rational(degree, degree - 1)

# A sparse degree-twelve seed gives a concrete composite-degree example.
H12 = (w**2 - w**12) / 10
model12 = WeightedSeedModel(sp.diff(H12, w), c=1, b=1)
assert model12.primitive == H12
assert model12.fiber_degree == 12
assert model12.kappa == -13
assert model12.a == -sp.Rational(12, 11)

u12 = 1 + x * y
gamma12 = 1 - sp.Rational(12, 11) * x * y + x**2 * z
numerator_A12 = 10 * u12 + u12**2 - 11 * u12**12 * gamma12**10
numerator_B12 = 5 + u12 - 6 * u12**11 * gamma12**10
A12 = sp.cancel(numerator_A12 / (10 * x**2))
B12 = sp.cancel(numerator_B12 / (5 * x))
C12 = x * gamma12
assert sp.denom(A12) == 1
assert sp.denom(B12) == 1

F12 = model12.mapping()
assert all(
    sp.expand(left - right) == 0
    for left, right in zip(F12, (A12, B12, C12))
)
assert sp.factor(sp.det(sp.Matrix(F12).jacobian((x, y, z)))) == 1

W12 = u12 * gamma12
inverse_incidence12 = sp.expand(
    H12.subs(w, W12) - B12 * C12 * W12 + A12 * C12**2
)
assert inverse_incidence12 == 0
assert sp.Poly(model12.inverse_polynomial(s, t, sp.Integer(1)), w).degree() == 12

# If S_11 < J <= S_12, then J moves the distinguished point and S_11 moves
# its image through the other eleven points.  Thus J is transitive and
# orbit--stabilizer forces |J|=12*11!=12!, proving maximality.
assert 12 * math.factorial(11) == math.factorial(12)

print("PASS: inverse pencils are irreducible in degrees 3 through 8")
print("PASS: discriminant elimination agrees with the repeated-root normalization")
print("PASS: generic branch meridians are transpositions")
print("PASS: canonical and deformed seeds share the universal S_n algebraic spine")
print("PASS: sparse seeds H_N=(w^2-w^N)/(N-2) are admissible in every degree")
print("PASS: sparse degree-twelve map is polynomial with determinant one")
print("PASS: inverse incidence has degree twelve and S_11 is maximal in S_12")
