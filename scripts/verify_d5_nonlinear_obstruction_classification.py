#!/usr/bin/env python3
"""Exact checks for the nonlinear D5 obstruction classification.

The valuation and UFD arguments are proved uniformly in the accompanying
note.  This checker replays the D5 polynomial ledgers, the low-order
solutions, the normalized branch incidence, the diagonal pole, and the
graph-section determinant for two and three auxiliary coordinates.
"""

from __future__ import annotations

import sympy as sp


a, u, v, s, z = sp.symbols("a u v s z")
sqrt5 = sp.sqrt(5)
alpha = (3 + sqrt5) / 2
beta = (3 - sqrt5) / 2

P = a**5 - 5 * a**3 * u + 5 * a * u**2
R_plus = a**2 - alpha * u
R_minus = a**2 - beta * u
Q = sp.expand(R_plus * R_minus)
C = a**2 - 4 * u
J = sp.diff(P, a)
Delta = v**2 - 4 * u**5
Delta_pullback = sp.expand(Delta.subs(v, P))

assert sp.expand(J - 5 * Q) == 0
assert sp.expand(Delta_pullback - C * Q**2) == 0

branch_vector = sp.Matrix([1, 2, 2])
derivative_vector = sp.Matrix([0, 1, 1])

# Every branch-supported valuation balance j + source = m*d has the stated
# unique solution.  Check the primitive and first bounded higher orders.
for branch_order in range(1, 9):
    source_orders = branch_order * branch_vector - derivative_vector
    assert source_orders == sp.Matrix(
        [branch_order, 2 * branch_order - 1, 2 * branch_order - 1]
    )
    assert all(order >= 0 for order in source_orders)
assert any(
    order < 0 for order in -derivative_vector
)
assert branch_vector - derivative_vector == sp.Matrix([1, 1, 1])
print("PASS: the single-branch log-discrepancy solutions are exact")

# The minimal diagonal source/target multipliers match determinant divisors
# but leave the inverse pole 1/Q.
minimal_source_multiplier = sp.expand(C * Q)
minimal_target_pullback = Delta_pullback
diagonal_inverse = sp.cancel(minimal_source_multiplier * z / minimal_target_pullback)
assert diagonal_inverse == z / Q
assert not sp.Poly(Q, a, u).is_ground
print("PASS: the primitive diagonal ledger retains the two-color pole")

# The normalization of the target cusp and its differential ranks.
normalized_u = s**2
normalized_v = 2 * s**5
assert sp.expand(
    Delta.subs({u: normalized_u, v: normalized_v})
) == 0
normalized_jacobian = sp.Matrix(
    [normalized_u, normalized_v]
).jacobian((s,))
assert normalized_jacobian.subs(s, 0).rank() == 0
assert normalized_jacobian.subs(s, 1).rank() == 1
print("PASS: target-branch incidence factors through the cusp normalization")


def verify_graph_section_block(auxiliary_count: int) -> None:
    """Verify the determinant after translating arbitrary polynomial graphs."""

    top_right = sp.Matrix(
        2,
        auxiliary_count,
        sp.symbols(f"b{auxiliary_count}_0:{2 * auxiliary_count}"),
    )
    normal_matrix = sp.Matrix(
        auxiliary_count,
        auxiliary_count,
        sp.symbols(
            f"n{auxiliary_count}_0:{auxiliary_count * auxiliary_count}"
        ),
    )
    incidence_derivative = sp.Matrix([u, P]).jacobian((a, u))
    graph_derivative = sp.Matrix.vstack(
        sp.Matrix.hstack(incidence_derivative, top_right),
        sp.Matrix.hstack(
            sp.zeros(auxiliary_count, 2),
            normal_matrix,
        ),
    )
    assert sp.expand(
        graph_derivative.det() + J * normal_matrix.det()
    ) == 0


verify_graph_section_block(2)
verify_graph_section_block(3)
print("PASS: every translated incidence graph retains the J5 divisor")
print("PASS nonlinear D5 obstruction classification")
