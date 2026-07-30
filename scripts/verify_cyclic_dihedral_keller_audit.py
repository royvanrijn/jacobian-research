#!/usr/bin/env python3
"""Exact regressions for the cyclic/dihedral inverse-Galois Keller audit.

The all-degree proofs in the note use P_n(x+y,xy)=x^n+y^n.  This checker
replays the polynomial, derivative, branch-pullback, and discriminant
identities through degree twelve and verifies the stated low-degree cards.
"""

import sympy as sp


A, U, V, X, Y, Z = sp.symbols("A U V X Y Z")


def dickson_power_sums(limit):
    """Return P_0,...,P_limit with P_m=A P_(m-1)-U P_(m-2)."""

    values = [sp.Integer(2), A]
    for _ in range(2, limit + 1):
        values.append(sp.expand(A * values[-1] - U * values[-2]))
    return values


power_sums = dickson_power_sums(12)

for n in range(2, 13):
    polynomial = power_sums[n]
    derivative = sp.diff(polynomial, A)

    # Uniform invariant and derivative identities in the splitting variables.
    assert sp.expand(
        polynomial.subs({A: X + Y, U: X * Y}) - X**n - Y**n
    ) == 0
    derivative_identity = sp.cancel(
        derivative.subs({A: X + Y, U: X * Y})
        - n * (X**n - Y**n) / (X - Y)
    )
    assert derivative_identity == 0

    # Pullback of the reduced target discriminant.
    branch_pullback = sp.factor(
        polynomial**2
        - 4 * U**n
        - (A**2 - 4 * U) * (derivative / n) ** 2
    )
    assert branch_pullback == 0

    discriminant = sp.factor(sp.discriminant(polynomial - V, A))
    if n % 2:
        expected = n**n * (4 * U**n - V**2) ** ((n - 1) // 2)
    else:
        m = n // 2
        expected = (
            n**n
            * (2 * U**m - V) ** (m - 1)
            * (2 * U**m + V) ** m
        )
    assert sp.factor(discriminant - expected) == 0

    source_unit_rank = 1 + n // 2
    target_unit_rank = 1 if n % 2 else 2
    assert source_unit_rank >= target_unit_rank

# Low-degree cards and the Jacobian of the derivative-unit suspension.
P3, P4, P5 = power_sums[3], power_sums[4], power_sums[5]
assert P3 == A**3 - 3 * A * U
assert P4 == A**4 - 4 * A**2 * U + 2 * U**2
assert P5 == A**5 - 5 * A**3 * U + 5 * A * U**2

J3 = sp.diff(P3, A)
chart = sp.Matrix([U, P3, Z / J3])
chart_jacobian = sp.factor(chart.jacobian((A, U, Z)).det())
assert chart_jacobian == -1

# A polynomial triangular pole clearing retaining the first two outputs has
# determinant -J_n*A_coeff; a nonzero constant would require 1/J_n.
A_coeff, B_coeff = sp.symbols("A_coeff B_coeff")
formal_triangular_determinant = sp.factor(-J3 * A_coeff)
assert formal_triangular_determinant == -3 * A_coeff * (A**2 - U)

# The rational n=5 Dickson row cannot have arithmetic group D5: the natural
# D5 is contained in A5, while the generic discriminant has constant square
# class 5.
disc5 = sp.factor(sp.discriminant(P5 - V, A))
assert disc5 == 5**5 * (4 * U**5 - V**2) ** 2
assert sp.factorint(5**5)[5] % 2 == 1

print("PASS: Dickson invariant and derivative identities for 2 <= n <= 12")
print("PASS: reduced-branch pullback is (A^2-4U)*(J_n/n)^2")
print("PASS: odd/even discriminant formulas for 2 <= n <= 12")
print("PASS: the derivative-unit suspension has Jacobian -1")
print("PASS: the cyclic regular-action and rational D5 audit flags are exact")
