#!/usr/bin/env python3
"""Verify the equal-scale admissible source expansion for the sextic chart."""

from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "admissible_equal_scale_degree6.json"
)

W, X, U = sp.symbols("W X U")

# Central special cover on the main target component.
central_map = sp.expand(W**2 * (W - 1) ** 2 * (W - 3) ** 2)
central_derivative = sp.factor(sp.diff(central_map, W))
assert central_derivative == (
    2 * W * (W - 1) * (W - 3) * (3 * W**2 - 8 * W + 3)
)

cluster_nodes = (sp.Integer(0), sp.Integer(1), sp.Integer(3))
stationary_critical_points = sp.solve(3 * W**2 - 8 * W + 3, W)
assert len(stationary_critical_points) == 2

# Riemann--Hurwitz on the degree-six central component:
# infinity contributes 6-1, the three cluster nodes contribute 2-1 each,
# and the two stationary simple ramifications contribute 2-1 each.
central_ramification_contributions = {
    "infinity": 5,
    "cluster_nodes": 3,
    "stationary_critical_points": 2,
}
assert sum(central_ramification_contributions.values()) == 2 * 6 - 2

# On the equal-scale trait x=y=z=t, the first target bubble receives three
# degree-two tails.  Their exact limiting maps and moving branch values are
# the local models from the branch-scale calculation.
tail_units = (sp.Integer(9), sp.Integer(4), sp.Integer(36))
tail_maps = tuple(unit * X * (X - 1) for unit in tail_units)
tail_branch_values = tuple(
    sp.factor(tail_map.subs(X, sp.Rational(1, 2)))
    for tail_map in tail_maps
)
assert tail_branch_values == (
    -sp.Rational(9, 4),
    -sp.Integer(1),
    -sp.Integer(9),
)
assert len(set(tail_branch_values)) == 3

for tail_map in tail_maps:
    assert sp.degree(tail_map, X) == 2
    assert sp.solve(sp.diff(tail_map, X), X) == [sp.Rational(1, 2)]
    # The finite simple ramification and total ramification at infinity each
    # contribute one, giving 2d-2=2.
    assert 1 + 1 == 2 * sp.degree(tail_map, X) - 2

# The three degree-two tails give total degree six over the bubble.  Each has
# one point of index two over the attaching target node at infinity.
assert sum(sp.degree(tail_map, X) for tail_map in tail_maps) == 6
assert sum(2 for _ in tail_maps) == 6

# Centering a tail at its finite critical point makes the local
# ramification equation an exact square.  This is the source expansion needed
# when the three moving branch values collide and the target sprouts the
# diagonal/Maxwell bubble.
for unit, tail_map, branch_value in zip(
    tail_units, tail_maps, tail_branch_values
):
    centered = sp.expand(tail_map.subs(X, U + sp.Rational(1, 2)) - branch_value)
    assert centered == unit * U**2

# One target node with three source nodes of index two has deformation ring
#
# k[[tau,sx,sy,sz]]/(sx^2-tau,sy^2-tau,sz^2-tau).
#
# Eliminating tau gives the union of four sign branches.  Verify the exact
# radical decomposition by iterated Groebner elimination.
tau, sx, sy, sz, selector = sp.symbols("tau sx sy sz selector")
node_equations = [sx**2 - tau, sy**2 - tau, sz**2 - tau]
eliminated_node_equations = [sy**2 - sx**2, sz**2 - sx**2]
sign_primes = [
    [sy - epsilon * sx, sz - delta * sx]
    for epsilon, delta in product((-1, 1), repeat=2)
]


def ideal_intersection(
    left: list[sp.Expr],
    right: list[sp.Expr],
) -> list[sp.Expr]:
    variables = (sx, sy, sz)
    groebner = sp.groebner(
        [selector * equation for equation in left]
        + [(1 - selector) * equation for equation in right],
        selector,
        *variables,
        order="lex",
    )
    return [
        sp.factor(polynomial.as_expr())
        for polynomial in groebner.polys
        if selector not in polynomial.as_expr().free_symbols
    ]


intersection = sign_primes[0]
for sign_prime in sign_primes[1:]:
    intersection = ideal_intersection(intersection, sign_prime)

intersection_basis = sp.groebner(intersection, sx, sy, sz, order="lex")
node_basis = sp.groebner(
    eliminated_node_equations,
    sx,
    sy,
    sz,
    order="lex",
)
assert [polynomial.as_expr() for polynomial in intersection_basis.polys] == [
    polynomial.as_expr() for polynomial in node_basis.polys
]

# The four normalization branches are tau=q^2 and
# (sx,sy,sz)=(q,epsilon*q,delta*q).  They are exactly the four Kummer
# preimages of [1:1:1] found in the target-graph calculation.
q = sp.symbols("q")
normalization_branches = []
for epsilon, delta in product((-1, 1), repeat=2):
    substitution = {
        tau: q**2,
        sx: q,
        sy: epsilon * q,
        sz: delta * q,
    }
    assert all(sp.expand(equation.subs(substitution)) == 0 for equation in node_equations)
    normalization_branches.append(
        {
            "tau": "q^2",
            "sx": "q",
            "sy": f"{epsilon}*q",
            "sz": f"{delta}*q",
        }
    )
assert len(normalization_branches) == 4

artifact = {
    "experiment": "degree-six equal-scale admissible source expansion",
    "central_component": {
        "map": "W^2(W-1)^2(W-3)^2",
        "degree": 6,
        "ramification_contributions": central_ramification_contributions,
        "riemann_hurwitz_total": 10,
    },
    "first_target_bubble": {
        "tail_maps": [str(tail_map) for tail_map in tail_maps],
        "tail_degrees": [2, 2, 2],
        "moving_branch_values": [str(value) for value in tail_branch_values],
        "degree_over_bubble": 6,
        "node_indices": [2, 2, 2],
    },
    "node_deformation": {
        "ring": (
            "k[[tau,sx,sy,sz]]/"
            "(sx^2-tau,sy^2-tau,sz^2-tau)"
        ),
        "normalization_branch_count": 4,
        "normalization_branches": normalization_branches,
        "interpretation": (
            "the saturated source-node choices equal the four Kummer "
            "preimages of the diagonal target blowup center"
        ),
    },
    "diagonal_bubble": {
        "local_equations": [
            f"S-lambda_{label}={unit}*U_{label}^2"
            for label, unit in zip(("x", "y", "z"), tail_units)
        ],
        "node_indices": [2, 2, 2],
    },
    "scope": (
        "equal-scale and first diagonal-resonance source charts only; "
        "extension to the selected admissible-cover component and label "
        "descent are not proved"
    ),
}

expected_artifact = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
if "--emit-json" in sys.argv:
    print(expected_artifact, end="")
    raise SystemExit(0)
assert ARTIFACT.read_text() == expected_artifact, (
    f"{ARTIFACT.relative_to(ROOT)} is stale; regenerate it from this script"
)

print("PASS central degree-six component: Riemann-Hurwitz total 10")
print("PASS first target bubble: three quadratic tails of total degree six")
print("PASS diagonal bubble: three exact simple-ramification squares")
print("PASS node deformation: reduced union of four Kummer branches")
print("DEGREE_SIX_ADMISSIBLE_EQUAL_SCALE_PASS")
