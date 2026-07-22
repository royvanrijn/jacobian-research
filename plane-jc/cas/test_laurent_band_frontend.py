#!/usr/bin/env python3
"""Regress the polygon-to-band-to-de Rham front end on both audited cases."""

import sympy as sp

from laurent_band_frontend import (
    audited_72_108_laurent_case,
    bracket_layers,
    compile_laurent_leading_block,
)
from newton_derham_compiler import compile_weighted_wronskian


t = sp.symbols("t")
compiled_cases = [
    compile_laurent_leading_block(audited_72_108_laurent_case(case))
    for case in (1, 2)
]
case1, case2 = compiled_cases

assert case1.source.chart.exponent_determinant == -1

assert len(case1.source.P_polygon.lattice_points()) == 61
assert len(case1.source.Q_polygon.lattice_points()) == 125
assert len(case2.source.P_polygon.lattice_points()) == 25
assert len(case2.source.Q_polygon.lattice_points()) == 47

P1 = dict(case1.P_bands)
Q1 = dict(case1.Q_bands)
P2 = dict(case2.P_bands)
Q2 = dict(case2.Q_bands)
assert tuple(P1) == tuple(range(2, -9, -1))
assert tuple(Q1) == tuple(range(3, -13, -1))
assert tuple(P2) == (2, 1, 0)
assert tuple(Q2) == (3, 2, 1, 0)
assert P1[2] == P2[2] == tuple(range(1, 9))
assert Q1[3] == Q2[3] == tuple(range(2, 13))
assert case1.rhs_tz_exponent == case2.rhs_tz_exponent == (2, 4)

# Both polygons compile to the same normalized first block and then to the
# same six de Rham obstruction coordinates.
for result in compiled_cases:
    block = result.weighted_wronskian
    assert block.covering_exponent == 2
    assert block.primitive_weight == 3
    assert block.A == t + sum(sp.Symbol(f"a{k}") * t**k for k in range(2, 8)) + t**8
    assert block.R == t**2
    assert block.primitive_exponents == tuple(range(2, 13))
    de_rham = compile_weighted_wronskian(block)
    assert de_rham.genus == 3
    assert de_rham.compact_dimension == 6
    assert de_rham.certificate.tail_basis_certified

# Compile all common upper bracket bands and recover J4,...,J0 exactly.
A, B, C, D, E, F, G = [sp.Function(name)(t) for name in "ABCDEFG"]
layers = bracket_layers(
    t,
    {2: A, 1: B, 0: C},
    {3: D, 2: E, 1: F, 0: G},
    sp.S.NegativeOne,
)
expected = {
    4: 2 * A * sp.diff(D, t) - 3 * sp.diff(A, t) * D,
    3: 2 * (A * sp.diff(E, t) - sp.diff(A, t) * E)
    + B * sp.diff(D, t)
    - 3 * sp.diff(B, t) * D,
    2: 2 * A * sp.diff(F, t)
    - sp.diff(A, t) * F
    + B * sp.diff(E, t)
    - 2 * sp.diff(B, t) * E
    - 3 * sp.diff(C, t) * D,
    1: 2 * A * sp.diff(G, t)
    + B * sp.diff(F, t)
    - sp.diff(B, t) * F
    - 2 * sp.diff(C, t) * E,
    0: B * sp.diff(G, t) - sp.diff(C, t) * F,
}
assert tuple(layers) == (4, 3, 2, 1, 0)
assert all(sp.expand(layers[layer] - expression) == 0 for layer, expression in expected.items())

print("PASS: both audited polygons compile to their exact Laurent band supports")
print("PASS: polygon top bands compile to the certified genus-three first block")
print("PASS: the front end reproduces all five upper bracket layers J4,...,J0")
