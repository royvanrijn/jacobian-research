#!/usr/bin/env python3
"""Audit the tangent-core starting point for a structural BCW bypass.

The foundational suspension square has five natural structural relations for
``gamma, W, C, s, t``.  Once these quantities are independent variables, all
five relations have degree at most three.  This is the right factored search
seed, but it is not itself a Keller stabilization: on the collision boundary
the five residual differentials have rank three.  Consequently no three
additional coordinate functions can complete all five residuals to an
eight-dimensional Keller map.

The audit separates these two statements exactly.  Future searches must use
the relations transiently in a chain of stable source/target shears; they may
not retain all five residuals as output coordinates.
"""

from itertools import combinations

import sympy as sp


x, y, z = sp.symbols("x y z")
A, B, C = sp.symbols("A B C")
W, gamma = sp.symbols("W gamma")

u = 1 + x * y
source_gamma = 1 - sp.Rational(3, 2) * x * y - sp.Rational(1, 2) * x**2 * z
source_W = u * source_gamma
source_C = x * source_gamma

H = W**2 * (1 - W)
core_s = sp.diff(H, W) + gamma
core_t = sp.expand(W * core_s - H)

# The weighted normalization of the foundational map.  Reordering and
# rescaling as (C,2B,A) gives the identity linear part used by the BCW route.
foundational_A = sp.expand(u**3 * z + y**2 * u * (4 + 3 * x * y))
foundational_B = sp.expand(
    (y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y)) / 2
)
foundational_C = sp.expand((2 * x - 3 * x**2 * y - x**3 * z) / 2)
normalized = sp.Matrix([foundational_C, 2 * foundational_B, foundational_A])
assert normalized.jacobian((x, y, z)).subs({x: 0, y: 0, z: 0}) == sp.eye(3)
assert sp.factor(normalized.jacobian((x, y, z)).det()) == 1

# The suspension square commutes exactly.
source_substitution = {W: source_W, gamma: source_gamma}
assert sp.expand(foundational_C - source_C) == 0
assert sp.factor(foundational_B * foundational_C - core_s.subs(source_substitution)) == 0
assert sp.factor(
    foundational_A * foundational_C**2 - core_t.subs(source_substitution)
) == 0

# The finite core and every structural relation are already quadratic-cubic.
assert sp.Poly(core_s, W, gamma).total_degree() == 2
assert sp.Poly(core_t, W, gamma).total_degree() == 3

residuals = sp.Matrix(
    [
        gamma - source_gamma,
        W - u * gamma,
        C - x * gamma,
        B * C - core_s,
        A * C**2 - core_t,
    ]
)
variables = (x, y, z, A, B, C, W, gamma)
residual_degrees = tuple(sp.Poly(entry, *variables).total_degree() for entry in residuals)
assert residual_degrees == (3, 3, 2, 2, 3)

# All three known collision points satisfy the same five structural equations.
collision_points = [
    (sp.Integer(0), sp.Integer(0), -sp.Rational(1, 4)),
    (sp.Integer(1), -sp.Rational(3, 2), sp.Rational(13, 2)),
    (-sp.Integer(1), sp.Rational(3, 2), sp.Rational(13, 2)),
]
target = {A: -sp.Rational(1, 4), B: 0, C: 0}
structural_lifts = []
for point in collision_points:
    source = dict(zip((x, y, z), point))
    lift = {
        **source,
        **target,
        W: sp.expand(source_W.subs(source)),
        gamma: sp.expand(source_gamma.subs(source)),
    }
    assert residuals.subs(lift) == sp.zeros(5, 1)
    structural_lifts.append(lift)
assert [(lift[W], lift[gamma], lift[C]) for lift in structural_lifts] == [
    (1, 1, 0),
    (0, 0, 0),
    (0, 0, 0),
]

# Obstruction to the tempting eight-dimensional residual-coordinate map.
# If five of eight output rows are these residuals, their row rank must be
# five everywhere: three further rows can raise rank by at most three.  On
# every lifted collision point the residual rank is only three.
residual_jacobian = residuals.jacobian(variables)
collision_ranks = tuple(residual_jacobian.subs(lift).rank() for lift in structural_lifts)
assert collision_ranks == (3, 3, 3)

# Stronger pruning for a coupled-shear search.  A retained subset of k
# residual coordinates must have row rank k at every collision lift, since
# the other 8-k output rows can increase the full rank by at most 8-k.
independent_subsets: dict[int, list[tuple[int, ...]]] = {}
for subset_size in range(1, 6):
    independent_subsets[subset_size] = []
    for indices in combinations(range(5), subset_size):
        ranks = tuple(
            residual_jacobian[list(indices), :].subs(lift).rank()
            for lift in structural_lifts
        )
        if all(rank == subset_size for rank in ranks):
            independent_subsets[subset_size].append(indices)

assert independent_subsets[5] == []
assert independent_subsets[4] == []
assert independent_subsets[3] == [(0, 1, 2), (0, 2, 3)]

# In particular r_t cannot be a terminal coordinate at all: its differential
# vanishes at the two lifts with (W,gamma,C)=(0,0,0).
rt_gradient_ranks = tuple(
    residual_jacobian[[4], :].subs(lift).rank() for lift in structural_lifts
)
assert rt_gradient_ranks == (1, 0, 0)

# As a bounded regression, none of the 56 completions by three input
# coordinates has constant nonzero determinant.  The rank argument above is
# stronger and also excludes arbitrary nonlinear choices of the other three
# outputs while the five residual rows are retained.
constant_coordinate_completions = []
for indices in combinations(range(len(variables)), 3):
    outputs = sp.Matrix([variables[index] for index in indices] + list(residuals))
    determinant = sp.factor(outputs.jacobian(variables).det())
    if determinant != 0 and not determinant.free_symbols:
        constant_coordinate_completions.append((indices, determinant))
assert constant_coordinate_completions == []

print("PASS tangent core: the foundational suspension square commutes exactly")
print("PASS tangent core: structural residual degrees are (3,3,2,2,3)")
print("PASS tangent core: all three collision points lift to the cubic incidence")
print("PASS tangent core: residual Jacobian ranks are (3,3,3) on the collision")
print("PASS tangent core: no Keller map can retain all five residuals as outputs")
print("PASS tangent core: no four residuals are jointly admissible terminal rows")
print("PASS tangent core: only the gamma-W-C and gamma-C-s triples survive rank pruning")
print("PASS tangent core: the t residual must be transient in every Keller trace")
