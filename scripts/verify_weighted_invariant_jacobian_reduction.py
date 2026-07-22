#!/usr/bin/env python3
"""Exact checks of the weighted invariant-coordinate Jacobian formula."""

import sympy as sp


x, y, z = sp.symbols("x y z", nonzero=True)
u_symbol, v_symbol = sp.symbols("u v")

A_polynomials = (
    1 + 2 * u_symbol + 3 * v_symbol + u_symbol * v_symbol,
    2 - u_symbol + v_symbol**2,
    3 + u_symbol**2 - 2 * v_symbol + u_symbol * v_symbol,
)


def verify_case(k: int, weights: tuple[int, int, int]) -> None:
    u = x * y
    v = x**k * z
    substitution = {u_symbol: u, v_symbol: v}
    components = tuple(
        sp.cancel(x**weight * polynomial.subs(substitution))
        for weight, polynomial in zip(weights, A_polynomials)
    )
    direct = sp.factor(sp.Matrix(components).jacobian((x, y, z)).det())
    reduced_matrix = sp.Matrix(
        [
            (
                weight * polynomial,
                sp.diff(polynomial, u_symbol),
                sp.diff(polynomial, v_symbol),
            )
            for weight, polynomial in zip(weights, A_polynomials)
        ]
    )
    reduced = sp.factor(
        x ** (sum(weights) + k)
        * reduced_matrix.det().subs(substitution)
    )
    assert sp.cancel(direct - reduced) == 0


test_cases = (
    (1, (-1, 0, 0)),
    (2, (-2, -1, 1)),
    (3, (-3, 1, -1)),
    (4, (-2, -1, -1)),
)
for test_case in test_cases:
    verify_case(*test_case)

assert all(sum(weights) == -k for k, weights in test_cases)

# The coordinate-change determinant used in the proof is x^(k+1).
for k in range(1, 6):
    u = x * y
    v = x**k * z
    coordinate_jacobian = sp.Matrix((x, u, v)).jacobian((x, y, z)).det()
    assert sp.factor(coordinate_jacobian - x ** (k + 1)) == 0

print("PASS: weighted Jacobian reduction for k=1,2,3,4")
print("PASS: balanced output weights remove the residual x-power")
print("PASS: invariant-coordinate Jacobian is x^(k+1) for k=1,...,5")
