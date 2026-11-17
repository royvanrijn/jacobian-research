#!/usr/bin/env python3
"""Verify the finite algebra in the even quartic--sextic HC4 closure.

The conceptual inputs are:

* Gordan--Noether/de Bondt classification of singular Hessians in the
  indicated low dimensions;
* HC4DCK for a constant common kernel direction; and
* HC4HQ1 for the sextic-free boundary.

This script checks the remaining internal algebra:

* bihomogeneity leaves (dual degree, cone degree)=(2,1) as the only
  possible nonconstant binary cone;
* that cone is a square of a bilinear form X^T M U; and
* the dual-degree-four part of the spatial z^4 determinant face is exactly
  48*c^4*det(M)^2*(X^T M U)^4, even with every compatible lower block and
  the source-only sextic block present; and
* the rank-one residual ternary Hessian has the required homogeneous
  scaling covariance.
"""

from __future__ import annotations

from itertools import permutations

import sympy as sp


# 1. If the highest nonlinear dual degree is d and the primitive
# projective cone [a:b] has homogeneous source degree delta, the scalar
# coefficient has degree 4-d-d*delta.  A moving cone has delta>=1.
admissible = [
    (dual_degree, cone_degree, 4 - dual_degree - dual_degree * cone_degree)
    for dual_degree in range(2, 5)
    for cone_degree in range(5)
    if 4 - dual_degree - dual_degree * cone_degree >= 0
]
moving_admissible = [
    row for row in admissible if row[1] >= 1
]
assert moving_admissible == [(2, 1, 0)]


# 2. The unique moving possibility is c*(X^T M U)^2.  Its full Hessian
# determinant is the claimed relative invariant.
x, y, t, m, z, epsilon, scalar = sp.symbols(
    "x y t m z epsilon scalar"
)
p, q, r, w = sp.symbols("p q r w")
variables = (t, m, x, y)
bilinear = (p * x + q * y) * t + (r * x + w * y) * m
moving_quartic = scalar * bilinear**2
moving_hessian = sp.hessian(moving_quartic, variables)
moving_determinant = sp.factor(moving_hessian.det())
expected_moving_determinant = (
    48
    * scalar**4
    * (p * w - q * r) ** 2
    * bilinear**4
)
assert sp.expand(moving_determinant - expected_moving_determinant) == 0


# 3. Verify that no compatible lower quartic block, quadratic block, or
# source-only sextic block changes the z^4, dual-degree-four coefficient.
# Order variables as U=(t,m), X=(x,y).  Besides the moving square, an
# at-most-dual-linear quartic has UU block zero, UX block of dual degree
# zero, and XX block of dual degree at most one.  The rank-two sextic is
# source-only, hence occupies only the XX block.
def symmetric_two(entries: tuple[sp.Expr, sp.Expr, sp.Expr]) -> sp.Matrix:
    return sp.Matrix(
        [[entries[0], entries[1]], [entries[1], entries[2]]]
    )


b_entries = sp.symbols("b0:4")
c0_entries = sp.symbols("c0:3")
ct_entries = sp.symbols("ct0:3")
cm_entries = sp.symbols("cm0:3")
lower_ux = sp.Matrix(2, 2, b_entries)
lower_xx = (
    symmetric_two(c0_entries)
    + t * symmetric_two(ct_entries)
    + m * symmetric_two(cm_entries)
)
lower_quartic_hessian = sp.zeros(4)
lower_quartic_hessian[:2, 2:] = lower_ux
lower_quartic_hessian[2:, :2] = lower_ux.T
lower_quartic_hessian[2:, 2:] = lower_xx
quartic_hessian = moving_hessian + lower_quartic_hessian

sextic_entries = sp.symbols("s0:3")
sextic_hessian = sp.zeros(4)
sextic_hessian[2:, 2:] = symmetric_two(sextic_entries)

h_entries = sp.symbols("h0:10")
quadratic_hessian = sp.Matrix(
    [
        [h_entries[0], h_entries[1], h_entries[2], h_entries[3]],
        [h_entries[1], h_entries[4], h_entries[5], h_entries[6]],
        [h_entries[2], h_entries[5], h_entries[7], h_entries[8]],
        [h_entries[3], h_entries[6], h_entries[8], h_entries[9]],
    ]
)

spatial_pencil = (
    quadratic_hessian
    + z * quartic_hessian
    + z**2 * sextic_hessian
).subs({t: epsilon * t, m: epsilon * m})


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


target_face = 0
for permutation in permutations(range(4)):
    determinant_term: sp.Expr = permutation_sign(permutation)
    for row, column in enumerate(permutation):
        determinant_term *= spatial_pencil[row, column]
    target_face += sp.Poly(
        determinant_term, z, epsilon
    ).coeff_monomial(z**4 * epsilon**4)

assert sp.expand(target_face - expected_moving_determinant) == 0


# 4. On the rank-one sextic boundary, write the quotient coordinate as s
# and its constant three-plane kernel as (u1,u2,u3).  Every quartic
# monomial satisfies Hess_U h(lambda*s,lambda*U)=lambda^2 Hess_U h(s,U).
s_base, u1, u2, u3, lam = sp.symbols("s_base u1 u2 u3 lambda")
rank_one_variables = (s_base, u1, u2, u3)
kernel_variables = (u1, u2, u3)


def weak_compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, length - 1):
            yield (first, *tail)


for exponents in weak_compositions(4, 4):
    monomial = sp.prod(
        variable**exponent
        for variable, exponent in zip(
            rank_one_variables, exponents, strict=True
        )
    )
    kernel_hessian = sp.hessian(monomial, kernel_variables)
    scaled_kernel_hessian = kernel_hessian.subs(
        {
            variable: lam * variable
            for variable in rank_one_variables
        }
    )
    assert all(
        sp.expand(entry) == 0
        for entry in scaled_kernel_hessian - lam**2 * kernel_hessian
    )


print("admissible (dual degree, cone degree, scalar degree):", admissible)
print("unique moving cone: dual degree 2, source degree 1")
print(
    "moving determinant:",
    "48*c^4*det(M)^2*(X^T M U)^4",
)
print("generic z^4 / dual-degree-4 face: no lower-layer contribution")
print("rank-one kernel Hessian scaling: 35 quartic monomials")
print("all even quartic--sextic HC4 closure checks passed")
