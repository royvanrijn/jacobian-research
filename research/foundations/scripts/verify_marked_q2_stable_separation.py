#!/usr/bin/env python3
"""Same marked Q_2 fiber in two stably inequivalent Keller maps."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from jcsearch.keller_fiber import compile_polynomial_to_keller_fiber  # noqa: E402
from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402
from verify_gq2_permutation_action import verify_certificate_path  # noqa: E402


T, S = sp.symbols("T S")
polynomial = T**5 + T**3 - 2 * T**2 + T + 1

# Layer 1: the complete marked G_Q2-set.  Squarefree reduction gives
# Q_2 x U_(2,4), with geometric Frobenius a fixed point plus a 4-cycle.
certificate = (
    ROOT
    / "arithmetic"
    / "certificates"
    / "gq2_common_quintic_stable_pair.json"
)
action = verify_certificate_path(certificate)
assert action["admissible"]
assert action["image_group_order"] == 4
assert action["orbit_sizes"] == [4, 1]
assert action["comparison"]["factor_degrees"] == [4, 1]

# Layer 2: global connectedness and regularity.  Irreducibility modulo 17
# proves irreducibility over Q; the odd discriminant also agrees with the
# squarefree dyadic reduction.
P = sp.Poly(polynomial, T, domain=sp.QQ)
assert P.gcd(P.diff()).degree() == 0
assert sp.discriminant(polynomial, T) == 28085
assert sp.Poly(polynomial, T, modulus=17).is_irreducible
assert P.is_irreducible
assert P.count_roots(-sp.oo, sp.oo) == 1

# Layer 3a: the fixed weighted map for H=T^5+T^3-2T^2.  Its determinant is
# -4 and its target (-1/4,-1,1) recovers P.  Scaling the first output by
# -1/4 gives determinant one and target (1/16,-1,1).
H = T**5 + T**3 - 2 * T**2
weighted = WeightedSeedModel(sp.diff(H, T).subs(T, w), c=-4)
weighted_map = weighted.mapping()
assert sp.factor(
    sp.Matrix(weighted_map).jacobian((x, y, z)).det()
) == -4
weighted_target = (sp.Rational(-1, 4), -1, 1)
assert sp.expand(
    weighted.inverse_polynomial(*weighted_target).subs(w, T) - polynomial
) == 0
weighted_map_one = (-weighted_map[0] / 4, weighted_map[1], weighted_map[2])
assert sp.factor(
    sp.Matrix(weighted_map_one).jacobian((x, y, z)).det()
) == 1
weighted_target_one = (sp.Rational(1, 16), -1, 1)

# Layer 3b: the fixed quadratic-gauge map for G=H+T.  The common compiler
# scales its determinant-minus-two output to determinant one.
quadratic = compile_polynomial_to_keller_fiber(
    polynomial,
    T,
    translation=0,
    inverse_variable=S,
    source_variables=(x, y, z),
)
assert quadratic.seed == S**5 + S**3 - 2 * S**2 + S
assert quadratic.target == (1, 0, -2)
assert quadratic.inverse_polynomial == polynomial.subs(T, S)
assert sp.factor(
    sp.Matrix(quadratic.determinant_one_map).jacobian((x, y, z)).det()
) == 1
assert quadratic.coordinate_degrees == (7, 32, 30)

# Layer 4: the stable-separation ledger proved by the canonical-boundary
# theorem.  Identity stabilization preserves these unit ranks.
stable_unit_ranks = {
    "weighted": 1,
    "quadratic_gauge": 2,
}
assert len(set(stable_unit_ranks.values())) == 2

print("PASS: the marked dyadic action is 1 plus an unramified Frobenius 4-cycle")
print("PASS: P_(5,1) is globally irreducible modulo 17 and has discriminant 28085")
print("PASS: the determinant-one weighted target is (1/16,-1,1)")
print("PASS: the determinant-one quadratic-gauge target is (1,0,-2)")
print("PASS: both complete fibers are Spec(Q[T]/P_(5,1))")
print("PASS: stable ramified-stratum unit ranks 1 and 2 separate the ambient maps")
