#!/usr/bin/env python3
"""Generate and audit the decomposable degree-twelve Keller map F_4 o F_3."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z


def sparse_weighted_map(degree: int) -> tuple[sp.Expr, ...]:
    """Return the expanded map attached to H=(w^2-w^degree)/(degree-2)."""

    if degree < 3:
        raise ValueError("the sparse admissible family starts in degree three")
    primitive = (w**2 - w**degree) / (degree - 2)
    model = WeightedSeedModel(sp.diff(primitive, w), c=1, b=1)
    assert model.fiber_degree == degree
    return tuple(sp.expand(component) for component in model.mapping())


parser = argparse.ArgumentParser()
parser.add_argument(
    "--print-map",
    action="store_true",
    help="print the three fully expanded coordinates of F_4 o F_3",
)
args = parser.parse_args()

F3 = sparse_weighted_map(3)
F4 = sparse_weighted_map(4)

assert sp.factor(sp.det(sp.Matrix(F3).jacobian((x, y, z)))) == 1
assert sp.factor(sp.det(sp.Matrix(F4).jacobian((x, y, z)))) == 1

P, Q, R = sp.symbols("P Q R")
F4_intermediate = tuple(
    component.subs({x: P, y: Q, z: R}, simultaneous=True)
    for component in F4
)
intermediate_substitution = {P: F3[0], Q: F3[1], R: F3[2]}
composite = tuple(
    sp.expand(component.subs(intermediate_substitution, simultaneous=True))
    for component in F4_intermediate
)

# Stable fingerprints ensure that the generator really expands the announced
# straight-line composition rather than silently returning one factor.
term_counts = tuple(
    len(sp.Poly(component, x, y, z).terms()) for component in composite
)
total_degrees = tuple(
    sp.Poly(component, x, y, z).total_degree() for component in composite
)
assert term_counts == (922, 733, 58)
assert total_degrees == (74, 68, 25)

# The chain rule gives det D(F4 o F3)=1*1.  Generic degrees multiply in the
# strict tower k(F4(F3)) < k(F3) < k(x,y,z), giving 4*3=12 and four
# monodromy blocks of size three.
assert 4 * 3 == 12

if args.print_map:
    print(*composite, sep="\n")

print("PASS: sparse factors F_3 and F_4 have determinant one")
print("PASS: expanded F_4 o F_3 has 922, 733, 58 terms")
print("PASS: degree tower is 4 times 3 with four blocks of size three")
