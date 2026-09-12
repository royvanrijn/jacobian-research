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
    refinement_branch_orbits,
    refinement_branches,
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
                count = component_decomposition_count(
                    collision, double_count, triple_count
                )
                assert count >= 1
                assert count == len(refinement_branches(maximal, collision))

# The first geometric multiple-branch fiber occurs in degree twelve.  The two
# multiplicity-six entries name two distinct generic roots, so exchanging
# their Q^2 and R^3 roles gives two distinct normalization points.
for degree in range(3, 12):
    for maximal in maximal_two_three_partitions(degree):
        for collision in full_contact_partitions(degree):
            if collision_precedes(maximal, collision):
                assert len(refinement_branches(maximal, collision)) == 1

degree_twelve = (3, 3, 2, 2, 2)
branches_twelve = refinement_branches(degree_twelve, (6, 6))
assert branches_twelve == (((3, 0), (0, 2)), ((0, 2), (3, 0)))
assert component_decomposition_count((6, 6), 3, 2) == 2

# This degree-twelve fiber is realized on the admissible open, not merely in
# the projective collision boundary.
W = sp.symbols("witness_W")
r, s = -sp.Rational(6, 5), sp.Integer(1)
M = sp.expand((W - r) ** 6 * (W - s) ** 6)
D = sp.expand(sp.diff(M, W).subs(W, 1) - sp.diff(M, W).subs(W, 0))
phi = sp.expand(M.subs(W, 1) - M.subs(W, 0) - sp.diff(M, W).subs(W, 0))
weighted_open = sp.expand(sp.diff(M, W, 2).subs(W, 1) - 2 * D)
assert phi == 0
assert D == sp.Rational(46656, 15625)
assert weighted_open == -sp.Rational(93312, 15625)
normalization_pairs = (
    (sp.expand((W - r) ** 3), sp.expand((W - s) ** 2)),
    (sp.expand((W - s) ** 3), sp.expand((W - r) ** 2)),
)
assert normalization_pairs[0] != normalization_pairs[1]
assert all(sp.expand(Q**2 * R**3 - M) == 0 for Q, R in normalization_pairs)

# If one additionally forgets the identities of equal-multiplicity collision
# roots, degree fourteen is the first ambiguity of abstract allocation types.
# This orbit count is deliberately not used as a geometric fiber count.
for degree in range(3, 14):
    for maximal in maximal_two_three_partitions(degree):
        for collision in full_contact_partitions(degree):
            if collision_precedes(maximal, collision):
                assert len(refinement_branch_orbits(maximal, collision)) == 1

degree_fourteen = (3, 3, 2, 2, 2, 2)
assert len(refinement_branch_orbits(degree_fourteen, (8, 6))) == 2

print("PASS: stable maximal Phi hypersurfaces are smooth on the admissible open")
print("PASS: all seven endpoint-rank singular ideals saturate to the unit ideal")
print("PASS: exact 2/3 polynomials give generic degree-one quotient maps")
print("PASS: collision fibers are counted by the local equations 2*i+3*j=m")
print("PASS: the first generic geometric two-branch collision is degree 12")
print("PASS: the first equal-multiplicity orbit-type ambiguity is degree 14")
print("PASS: every maximal seed component has an explicit smooth normalization")
