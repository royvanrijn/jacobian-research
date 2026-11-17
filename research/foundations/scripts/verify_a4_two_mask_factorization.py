#!/usr/bin/env python3
"""Exact two-mask target factorization checks for the A4 lift."""

import itertools

import sympy as sp


U, V, W, z1, z2 = sp.symbols("U V W z1 z2")

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


def target_cubic(first, second, third):
    return sp.expand(
        first**3
        - 3 * first * second**2
        + 2 * second**3
        - 9 * first * second * third
        + 9 * second**2 * third
        - 27 * first * third**2
        + 27 * second * third**2
        + 27 * third**3
    )


def target_matrix_data(first, second, third):
    difference = sp.expand(first - second)
    first_incidence = sp.expand(
        27 * third**2
        - difference**2
        - 3 * (difference - 3 * third) * second
    )
    second_incidence = sp.expand(
        3 * difference * second
        + difference**2
        + 3 * difference * third
        - 18 * third**2
    )
    return difference, first_incidence, second_incidence


# ---------------------------------------------------------------------------
# 1. The target cubic has a two-puncture determinantal representation
# ---------------------------------------------------------------------------

P, Q, R, A, B = sp.symbols("P Q R A B")
X, U_incidence, V_incidence = target_matrix_data(P, Q, R)
target_B = target_cubic(P, Q, R)

assert sp.expand(target_B - (27 * R**3 - X * U_incidence)) == 0
assert sp.expand(
    target_B - ((X - 3 * R) * V_incidence - 27 * R**3)
) == 0

target_vertical_outputs = sp.Matrix([
    27 * R**3 * A + X * B,
    U_incidence * A + B,
])
vertical_jacobian = sp.factor(
    target_vertical_outputs.jacobian((A, B)).det()
)
assert sp.expand(vertical_jacobian - target_B) == 0

S, T = sp.symbols("S T")
inverse_A = sp.cancel((S - X * T) / target_B)
inverse_B = sp.cancel(
    (-U_incidence * S + 27 * R**3 * T) / target_B
)
round_trip = target_vertical_outputs.subs(
    {A: inverse_A, B: inverse_B},
    simultaneous=True,
)
assert sp.factor(round_trip[0] - S) == 0
assert sp.factor(round_trip[1] - T) == 0


# ---------------------------------------------------------------------------
# 2. All coordinate permutations of the pure lift plus identity fail
# ---------------------------------------------------------------------------

outputs = [
    W * N1,
    W * N2,
    W * H,
    W * L * z1 / 4,
    z2,
]
variables = (U, V, W, z1, z2)
survivors = []

for permutation in itertools.permutations(range(5)):
    first, second, third, fourth, fifth = [
        outputs[index] for index in permutation
    ]
    difference, incidence, _ = target_matrix_data(first, second, third)
    denominator = target_cubic(first, second, third)
    numerator_1 = sp.expand(fourth - difference * fifth)
    numerator_2 = sp.expand(
        -incidence * fourth + 27 * third**3 * fifth
    )

    _, remainder_1 = sp.div(
        numerator_1,
        denominator,
        *variables,
    )
    if remainder_1 != 0:
        continue

    _, remainder_2 = sp.div(
        numerator_2,
        denominator,
        *variables,
    )
    if remainder_2 == 0:
        survivors.append(permutation)

assert survivors == []


# ---------------------------------------------------------------------------
# 3. The six simplest singular-line incidence charts fail
# ---------------------------------------------------------------------------

cone_outputs = outputs[:3]
singular_line_survivors = []
for repeated_index in range(3):
    remaining_indices = [
        index for index in range(3)
        if index != repeated_index
    ]
    for mask_order in (
        remaining_indices,
        list(reversed(remaining_indices)),
    ):
        first = cone_outputs[repeated_index]
        second = first + outputs[3]
        third = outputs[4]
        fourth = cone_outputs[mask_order[0]]
        fifth = cone_outputs[mask_order[1]]

        difference, incidence, _ = target_matrix_data(
            first,
            second,
            third,
        )
        denominator = target_cubic(first, second, third)
        numerator_1 = sp.expand(fourth - difference * fifth)
        numerator_2 = sp.expand(
            -incidence * fourth + 27 * third**3 * fifth
        )
        _, remainder_1 = sp.div(
            numerator_1,
            denominator,
            *variables,
        )
        if remainder_1 != 0:
            continue
        _, remainder_2 = sp.div(
            numerator_2,
            denominator,
            *variables,
        )
        if remainder_2 == 0:
            singular_line_survivors.append(
                (repeated_index, tuple(mask_order))
            )

assert singular_line_survivors == []


# ---------------------------------------------------------------------------
# 4. Polynomial normalization parametrization of the target cubic
# ---------------------------------------------------------------------------

lam, parameter = sp.symbols("lam parameter")
normalized_R = lam * parameter * (parameter - 1)
normalized_Q = lam * (-parameter**3 + 3 * parameter - 1)
normalized_P = lam * (
    2 * parameter**3
    - 3 * parameter**2
    + 3 * parameter
    - 1
)
assert sp.factor(
    target_cubic(normalized_P, normalized_Q, normalized_R)
) == 0


# ---------------------------------------------------------------------------
# 5. The one-primitive constant-term divisibility gate
# ---------------------------------------------------------------------------

base_pullback = target_cubic(outputs[0], outputs[1], outputs[2])
assert sp.factor(base_pullback - W**3 * K**3 * L**2) == 0
primitive_coefficient = W * L / 4
_, primitive_remainder = sp.div(
    primitive_coefficient,
    base_pullback,
    U,
    V,
    W,
)
assert primitive_remainder != 0


print("PASS: the target cubic has two puncture-adapted Bezout identities")
print("PASS: the two-mask target blowdown has Jacobian B and rational inverse")
print("PASS: every one-primitive triangular rechart fails divisibility")
print("PASS: all 120 coordinate assignments fail two-mask factorization")
print("PASS: all six singular-line incidence recharts fail factorization")
print("PASS: the smooth target boundary has a polynomial normalization chart")
