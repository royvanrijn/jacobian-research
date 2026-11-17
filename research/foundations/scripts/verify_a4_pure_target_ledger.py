#!/usr/bin/env python3
"""Exact checks for the pure-target A4 determinant-ledger lift."""

import itertools

import sympy as sp


U, V, W, z = sp.symbols("U V W z")

H = (
    8 * U**3
    - 6 * U * V**2
    - 18 * U * V
    - 54 * U
    - 2 * V**3
    - 9 * V**2
    - 27 * V
    - 27
)
K = 4 * U**2 + 4 * U * V + 6 * U + V**2 + 3 * V + 9
M = U**2 + 2 * V**2 + 6 * V + 18
L = (
    U**3
    - 3 * U * V**2
    - 9 * U * V
    - 27 * U
    + 2 * V**3
    + 9 * V**2
    + 27 * V
    + 27
)
N1 = sp.expand(M * K)
N2 = (
    8 * U**3 * V
    + 12 * U**2 * V**2
    + 36 * U**2 * V
    + 108 * U**2
    + 6 * U * V**3
    + 36 * U * V**2
    + 108 * U * V
    + 162 * U
    + V**4
    + 9 * V**3
    + 27 * V**2
    + 54 * V
)

P = W * N1
Q = W * N2
R = W * H
S = W * L * z / 4


# ---------------------------------------------------------------------------
# 1. Pure-target determinant identity
# ---------------------------------------------------------------------------

target_B = (
    P**3
    - 3 * P * Q**2
    + 2 * Q**3
    - 9 * P * Q * R
    + 9 * Q**2 * R
    - 27 * P * R**2
    + 27 * Q * R**2
    + 27 * R**3
)
lift = sp.Matrix([P, Q, R, S])
lift_jacobian = sp.factor(lift.jacobian((U, V, W, z)).det())

assert sp.factor(target_B - W**3 * K**3 * L**2) == 0
assert lift_jacobian == W**3 * K**3 * L**2
assert sp.factor(lift_jacobian - target_B) == 0
assert sp.cancel(lift_jacobian / target_B) == 1

# The new source coordinate is reconstructed generically.
reconstructed_z = sp.cancel(4 * S / (W * L))
assert reconstructed_z == z


# ---------------------------------------------------------------------------
# 2. Double-incidence encoding of the three-puncture residual cubic
# ---------------------------------------------------------------------------

r, v_plane = sp.symbols("r v_plane")
scaled_L = sp.factor(L.subs({U: v_plane + 3 * r, V: v_plane}) / 27)
expected_scaled_L = r * (r - 1) * v_plane + r**3 - 3 * r + 1
assert sp.expand(scaled_L - expected_scaled_L) == 0

u_incidence = 3 - r**2 - (r - 1) * v_plane
v_incidence = r * v_plane + r**2 + r - 2
assert sp.expand(1 - r * u_incidence - scaled_L) == 0
assert sp.expand(1 - (r - 1) * v_incidence + scaled_L) == 0


# ---------------------------------------------------------------------------
# 3. Every pair of cone-coordinate differentials has a W factor
# ---------------------------------------------------------------------------

cone_coordinates = [P, Q, R]
variables = [U, V, W]
for first, second in itertools.combinations(cone_coordinates, 2):
    jacobian_rows = sp.Matrix([first, second]).jacobian(variables)
    for first_column, second_column in itertools.combinations(range(3), 2):
        minor = sp.factor(
            jacobian_rows[:, [first_column, second_column]].det()
        )
        assert sp.rem(minor, W, W) == 0


# ---------------------------------------------------------------------------
# 4. The naive target scaling quotient is not polynomial
# ---------------------------------------------------------------------------

naive_fourth_coordinate = sp.cancel(S / target_B)
assert sp.factor(
    naive_fourth_coordinate - z / (4 * W**2 * K**3 * L)
) == 0


print("PASS: the fourfold A4 lift has determinant equal to target B")
print("PASS: the logarithmic target volume pulls back to ordinary volume")
print("PASS: the added source coordinate is generically reconstructed")
print("PASS: L has two polynomial incidence defects with opposite signs")
print("PASS: retaining any two cone outputs leaves a W divisor")
print("PASS: the naive multiplicative target factorization is nonpolynomial")
