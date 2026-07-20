#!/usr/bin/env python3
"""Smooth quotient normalizations of the maximal exceptional components."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    collision_precedes,
    component_decomposition_count,
    maximal_two_three_partitions,
    maximal_two_three_phi,
)


def full_contact_partitions(total, maximum=None):
    if total == 0:
        yield ()
        return
    maximum = min(total, maximum or total)
    for first in range(maximum, 1, -1):
        for tail in full_contact_partitions(total - first, first):
            yield (first,) + tail


# Stable endpoint-coordinate proof.  On Phi=0 and D!=0, x=Q(0) and y=R(0)
# are units: if either vanished, Phi would force Q(1)R(1)=0, making both
# endpoint derivatives of M vanish and hence D=0.
x, u, X, y, v, Y = sp.symbols("x u X y v Y")
universal_phi = sp.expand(
    X**2 * Y**3 - x**2 * y**3 - 2 * x * u * y**3 - 3 * x**2 * y**2 * v
)
assert sp.diff(universal_phi, u) == -2 * x * y**3
assert sp.diff(universal_phi, v) == -3 * x**2 * y**2

# For a>=3, u is independent; for b>=3, v is independent.  The seven cases
# with neither stable coordinate are checked by saturating their singular
# ideals only by D and weighted admissibility.  Root collisions are retained.
small_cases = ((0, 1), (0, 2), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2))
for double_count, triple_count in small_cases:
    model = maximal_two_three_phi(
        double_count,
        triple_count,
        prefix=f"normal_{double_count}_{triple_count}",
    )
    W = model.variable
    coordinates = model.double_coordinates + model.triple_coordinates
    D = sp.expand(
        sp.diff(model.M, W).subs(W, 1) - sp.diff(model.M, W).subs(W, 0)
    )
    weighted_open = sp.expand(sp.diff(model.M, W, 2).subs(W, 1) - 2 * D)
    gate = sp.symbols(f"normal_gate_{double_count}_{triple_count}")
    singular_equations = (
        model.phi,
        *(sp.diff(model.phi, coordinate) for coordinate in coordinates),
        1 - gate * D * weighted_open,
    )
    basis = sp.groebner(singular_equations, gate, *coordinates, order="lex")
    assert any(polynomial.as_expr() == 1 for polynomial in basis.polys)

# Generic exact 2/3 polynomials determine Q and R uniquely.  At collisions,
# normalization fibers are counted by local solutions 2*i+3*j=m with global
# Q,R degrees fixed.
for degree in range(3, 15):
    for maximal in maximal_two_three_partitions(degree):
        double_count = maximal.count(2)
        triple_count = maximal.count(3)
        assert component_decomposition_count(
            maximal, double_count, triple_count
        ) == 1
        for collision in full_contact_partitions(degree):
            if collision_precedes(maximal, collision):
                assert component_decomposition_count(
                    collision, double_count, triple_count
                ) >= 1

# The first simple self-identification example: at two multiplicity-six roots,
# the (a,b)=(3,2) normalization can assign Q^2 and R^3 to the two roots in
# either order.
assert component_decomposition_count((6, 6), 3, 2) == 2

print("PASS: stable maximal Phi hypersurfaces are smooth on the admissible open")
print("PASS: all seven endpoint-rank singular ideals saturate to the unit ideal")
print("PASS: exact 2/3 polynomials give generic degree-one quotient maps")
print("PASS: collision fibers are counted by the local equations 2*i+3*j=m")
print("PASS: every maximal seed component has an explicit smooth normalization")
