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

# The target invariant quotient has Jacobian -Lambda^k times the reduced
# three-row determinant.
A, B, Lambda = A_polynomials
for k in range(1, 6):
    P = B * Lambda
    Q = A * Lambda**k
    quotient_jacobian = sp.Matrix((P, Q)).jacobian(
        (u_symbol, v_symbol)
    ).det()
    reduced_matrix = sp.Matrix(
        (
            (-k * A, sp.diff(A, u_symbol), sp.diff(A, v_symbol)),
            (-B, sp.diff(B, u_symbol), sp.diff(B, v_symbol)),
            (
                Lambda,
                sp.diff(Lambda, u_symbol),
                sp.diff(Lambda, v_symbol),
            ),
        )
    )
    assert sp.factor(
        quotient_jacobian + Lambda**k * reduced_matrix.det()
    ) == 0

# Every admissible k=2 seed has a quotient factorization whose two plane
# Jacobians contribute one copy of gamma each.  Treat H and its derivatives
# as independent symbolic functions to keep the check degree-free.
a0, b0, c = sp.symbols("a0 b0 c", nonzero=True)
gamma = 1 + a0 * u_symbol + b0 * v_symbol
W = (1 + u_symbol) * gamma
W_symbol, gamma_symbol = sp.symbols("W gamma")
H = sp.Function("H")
P_seed = sp.diff(H(W_symbol), W_symbol) + c * gamma_symbol
Q_seed = (
    W_symbol * P_seed - H(W_symbol)
) / c
theta_jacobian = sp.factor(
    sp.Matrix((W, gamma)).jacobian((u_symbol, v_symbol)).det()
)
psi_jacobian = sp.simplify(
    sp.Matrix((P_seed, Q_seed)).jacobian(
        (W_symbol, gamma_symbol)
    ).det()
)
assert sp.factor(theta_jacobian - b0 * gamma) == 0
assert sp.simplify(psi_jacobian + c * gamma_symbol) == 0

# The restrictions to x=0 used in the stabilizer table have nonzero linear
# y- and z-coefficients for every admissible seed.  Expand H' and q through
# the required order at W=1; higher seed coefficients only affect y^2 in A.
epsilon, y_boundary, z_boundary = sp.symbols("epsilon y_boundary z_boundary")
kappa, p_second = sp.symbols("kappa p_second")
a0_admissible = -(1 + kappa) / (2 + kappa)
gamma_boundary = (
    1
    + a0_admissible * epsilon * y_boundary
    + b0 * epsilon**2 * z_boundary
)
W_boundary = (1 + epsilon * y_boundary) * gamma_boundary
delta_W = sp.expand(W_boundary - 1)
p_taylor = -c + c * kappa * delta_W + p_second * delta_W**2 / 2
q_second = kappa + p_second / c
q_taylor = -1 + kappa * delta_W + q_second * delta_W**2 / 2
B_numerator = c + p_taylor / gamma_boundary
A_numerator = (
    1 + epsilon * y_boundary + q_taylor / gamma_boundary**2
)
B_boundary = sp.simplify(
    sp.diff(B_numerator, epsilon).subs(epsilon, 0)
)
A_z_coefficient = sp.simplify(
    sp.diff(
        sp.diff(A_numerator, epsilon, 2).subs(epsilon, 0) / 2,
        z_boundary,
    )
)
assert sp.simplify(B_boundary + c * y_boundary / (2 + kappa)) == 0
assert sp.simplify(A_z_coefficient - b0 * (2 + kappa)) == 0

print("PASS: weighted Jacobian reduction for k=1,2,3,4")
print("PASS: balanced output weights remove the residual x-power")
print("PASS: invariant-coordinate Jacobian is x^(k+1) for k=1,...,5")
print("PASS: target quotient Jacobian is -Lambda^k times the reduced determinant")
print("PASS: admissible-seed quotient splits -b0*c*gamma^2 into two simple factors")
print("PASS: admissible boundary restrictions preserve the mu_2 stabilizer lines")
