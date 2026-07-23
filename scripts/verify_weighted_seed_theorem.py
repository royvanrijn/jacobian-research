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

print("PASS: inverse pencils are irreducible in degrees 3 through 8")
print("PASS: discriminant elimination agrees with the repeated-root normalization")
print("PASS: generic branch meridians are transpositions")
print("PASS: canonical and deformed seeds share the universal S_n algebraic spine")
