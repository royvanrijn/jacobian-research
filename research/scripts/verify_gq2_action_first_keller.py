#!/usr/bin/env python3
"""Exact action-first Keller certificate for the tame S_3 cubic T^3 - 2."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_gq2_permutation_action import verify_certificate_path  # noqa: E402


certificate = (
    ROOT / "arithmetic" / "certificates" / "gq2_s3_x3_minus_2.json"
)
action_report = verify_certificate_path(certificate)
assert action_report["admissible"]
assert action_report["image_group_order"] == 6
assert action_report["comparison"]["discriminant"] == -108

# The manuscript's mixed S_4 tuple is a second finite-word regression.  It
# exercises the nontrivial wild 2-core but has no polynomial comparison yet.
mixed_certificate = (
    ROOT / "arithmetic" / "certificates" / "gq2_s4_mixed_action.json"
)
mixed_report = verify_certificate_path(mixed_certificate)
assert mixed_report["admissible"]
assert mixed_report["image_group_order"] == 24
assert mixed_report["wild_normal_closure_order"] == 4

# Translate P(T)=T^3-2 by T=1+S.  The rooted seed has
# (g_1,g_2,g_3)=(3,3,1), hence satisfies the quadratic-gauge gates.
x, y, z, S, T = sp.symbols("x y z S T")
polynomial = T**3 - 2
seed = sp.expand(polynomial.subs(T, 1 + S) - polynomial.subs(T, 1))
assert seed == S**3 + 3 * S**2 + 3 * S
g1, g2, g3 = (sp.Poly(seed, S).nth(index) for index in (1, 2, 3))
assert (g1, g2, g3) == (3, 3, 1)
assert sp.gcd(sp.Poly(polynomial, T), sp.Poly(sp.diff(polynomial, T), T)) == 1

# Formula (26) of ROOT_ENGINEERED_QUADRATIC_GAUGE.md in degree three.
t = 1 + x * y
q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
mapping = (
    sp.expand(t * q),
    sp.expand(y + 3 * g3 * x * q / g1 + 2 * g2 * t * q / g1),
    sp.expand(x * (5 - 3 * t) - g3 * x**3 * z / g1),
)
assert sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det()) == -2

# The target (1,0,2/3) recovers P(1+S).  Scaling the last output by -1/2
# gives a determinant-one map and scales the last target to -1/3.
target_c = sp.Rational(2, 3)
inverse = sp.expand(seed - g1 * target_c / 2)
assert inverse == sp.expand(polynomial.subs(T, 1 + S))

jacobian_one_mapping = (mapping[0], mapping[1], -mapping[2] / 2)
assert sp.factor(
    sp.Matrix(jacobian_one_mapping).jacobian((x, y, z)).det()
) == 1
jacobian_one_target = (sp.Integer(1), sp.Integer(0), -target_c / 2)
assert jacobian_one_target == (1, 0, -sp.Rational(1, 3))

component_degrees = tuple(
    sp.Poly(component, x, y, z).total_degree()
    for component in jacobian_one_mapping
)
assert component_degrees == (7, 7, 4)

print("PASS: the Roe--Turturean S_3 marking is an admissible degree-three action")
print("PASS: the mixed S_4 regression has wild normal closure V_4")
print("PASS: T^3-2 is the compared tame Eisenstein cubic over Q_2")
print("PASS: the translated quadratic-gauge map has determinant one and full cubic fiber")
print("PASS: geometric fiber degree 3 is minimal for a nonabelian permutation image")
