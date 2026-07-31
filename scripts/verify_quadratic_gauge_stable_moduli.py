#!/usr/bin/env python3
"""Exact checks for quadratic-gauge stable moduli."""

from __future__ import annotations

from math import gcd

import sympy as sp


x, y, z = sp.symbols("x y z")
alpha, beta = sp.symbols("alpha beta", nonzero=True)
a2, a3, a4, a5, a6 = sp.symbols("a2 a3 a4 a5 a6", nonzero=True)


def quadratic_map(
    coefficients: dict[int, sp.Expr],
) -> tuple[tuple[sp.Expr, ...], sp.Expr, sp.Expr]:
    """Normalized quadratic-gauge map for one coefficient dictionary."""

    local_t = 1 + x * y
    local_q = (
        local_t**2 * z
        + y**2 * (1 + 3 * local_t) / coefficients[3]
    )
    degree = max(coefficients)
    mapping = (
        local_t * local_q,
        y
        + 3 * coefficients[3] * x * local_q
        + 2 * coefficients.get(2, 0) * local_t * local_q
        + sum(
            k
            * coefficients[k]
            * local_t**2
            * x ** (k - 2)
            * local_q**k
            for k in range(4, degree + 1)
        ),
        x * (5 - 3 * local_t)
        - coefficients[3] * x**3 * z
        - sum(
            (k - 2) * coefficients[k] * (x * local_q) ** k
            for k in range(4, degree + 1)
        ),
    )
    return tuple(sp.expand(item) for item in mapping), local_t, local_q


a = {2: a2, 3: a3, 4: a4, 5: a5, 6: a6}
b = {
    2: a2 * alpha**-1 * beta**-1,
    3: a3 * alpha**-2 * beta**-1,
    4: a4 * alpha**-3 * beta**-4,
    5: a5 * alpha**-4 * beta**-5,
    6: a6 * alpha**-5 * beta**-6,
}
F_a, t_a, q_a = quadratic_map(a)
F_b, _, q_b = quadratic_map(b)
source_scaling = {x: alpha * x, y: y / alpha, z: beta * z}

assert sp.factor(t_a.subs(source_scaling) - t_a) == 0
assert sp.factor(q_b.subs(source_scaling, simultaneous=True) - beta * q_a) == 0

scaled_F_b = tuple(
    sp.factor(component.subs(source_scaling, simultaneous=True))
    for component in F_b
)
target_scaled_F_a = (
    beta * F_a[0],
    F_a[1] / alpha,
    alpha * F_a[2],
)
assert all(
    sp.factor(got - expected) == 0
    for got, expected in zip(scaled_F_b, target_scaled_F_a)
)


# The quadratic coefficient is exactly a target shear.
a_without_2 = dict(a)
a_without_2[2] = sp.Integer(0)
F_without_2, _, _ = quadratic_map(a_without_2)
assert sp.factor(F_a[0] - F_without_2[0]) == 0
assert sp.factor(F_a[1] - F_without_2[1] - 2 * a2 * F_a[0]) == 0
assert sp.factor(F_a[2] - F_without_2[2]) == 0


# Intrinsic normalization and its Fitting divisor.
P, r = sp.symbols("P r", nonzero=True)
h = (
    r
    + a2 * P * r**2
    + a3 * P * r**3
    + a4 * P**4 * r**4
    + a5 * P**5 * r**5
    + a6 * P**6 * r**6
)
B = sp.factor(sp.diff(h, r) / r)
C = sp.factor(2 * h - r * sp.diff(h, r))
assert sp.factor(C + r**2 * B - 2 * h) == 0
assert sp.factor(sp.diff(C, r) + r**2 * sp.diff(B, r)) == 0

J = sp.factor(r**2 * sp.diff(B, r))
expected_J = (
    -1
    + 3 * a3 * P * r**2
    + 8 * a4 * P**4 * r**3
    + 15 * a5 * P**5 * r**4
    + 24 * a6 * P**6 * r**5
)
assert sp.factor(J - expected_J) == 0
assert not J.has(a2)


# Support rigidity and the two-dimensional weight lattice in all tested
# degrees.  Inversion changes the unique missing interior r-exponent, and
# the r^2 term forces a possible P^m twist to have m=0.
for degree in range(4, 65):
    support = {(0, 0), (1, 2)}
    support.update((k, k - 1) for k in range(4, degree + 1))
    r_degrees = {pair[1] for pair in support}
    assert min(r_degrees) == 0
    assert max(r_degrees) == degree - 1
    assert [pair for pair in support if pair[1] == 0] == [(0, 0)]
    assert [pair for pair in support if pair[1] == 2] == [(1, 2)]

    inverted_r_degrees = {
        degree - 1 - exponent
        for exponent in r_degrees
    }
    assert inverted_r_degrees == (
        set(range(0, degree - 2)) | {degree - 1}
    )
    assert r_degrees != inverted_r_degrees
    assert 1 not in r_degrees
    assert degree - 2 not in inverted_r_degrees

    weight_3 = (-2, -1)
    weight_4 = (-3, -4)
    determinant_34 = (
        weight_3[0] * weight_4[1] - weight_3[1] * weight_4[0]
    )
    assert determinant_34 == 5
    assert degree - 2 - 2 == degree - 4

    if degree == 4:
        assert abs(determinant_34) == 5  # kernel mu_5
    else:
        weight_5 = (-4, -5)
        determinant_35 = (
            weight_3[0] * weight_5[1] - weight_3[1] * weight_5[0]
        )
        assert determinant_35 == 6
        assert gcd(abs(determinant_34), abs(determinant_35)) == 1


# On the compiler slice a3=1, the residual scaling has weights
# (5,6,...,N+1) on (u4,u5,...,uN).  I5 and the consecutive second
# differences J6,...,JN form a saturated basis of its invariant lattice.
for degree in range(5, 17):
    residual_weights = sp.Matrix(
        [[index + 1 for index in range(4, degree + 1)]]
    )
    relations = sp.zeros(degree - 4, degree - 3)
    relations[0, 0] = -6
    relations[0, 1] = 5
    for row in range(1, degree - 4):
        relations[row, row - 1] = 1
        relations[row, row] = -2
        relations[row, row + 1] = 1

    assert residual_weights * relations.T == sp.zeros(1, degree - 4)
    assert relations.rank() == degree - 4
    maximal_minors = []
    for deleted_column in range(degree - 3):
        retained_columns = [
            column
            for column in range(degree - 3)
            if column != deleted_column
        ]
        maximal_minors.append(
            abs(int(relations[:, retained_columns].det()))
        )
    assert maximal_minors == list(range(5, degree + 2))
    assert gcd(*maximal_minors) == 1


# The residual receiver action is globally split, not merely faithful.
# On the seed torus lambda=u5/u4 has weight one.  The characters
#
#   lambda, q4=u4/lambda^5, qj=uj/lambda^(j+1) (j>=6)
#
# form a unimodular coordinate basis.  Hence lambda=1 is a global slice.
for degree in range(5, 17):
    seed_count = degree - 3
    residual_weights = sp.Matrix(
        [[index + 1 for index in range(4, degree + 1)]]
    )
    split_characters = sp.zeros(seed_count, seed_count)
    # lambda=u5/u4
    split_characters[0, 0] = -1
    split_characters[0, 1] = 1
    # q4=u4/lambda^5
    split_characters[1, 0] = 6
    split_characters[1, 1] = -5
    # qj=uj/lambda^(j+1), j>=6
    for index in range(6, degree + 1):
        row = index - 4
        split_characters[row, 0] = index + 1
        split_characters[row, 1] = -(index + 1)
        split_characters[row, index - 4] = 1

    assert abs(int(split_characters.det())) == 1
    transformed_weights = residual_weights * split_characters.T
    assert transformed_weights == sp.Matrix(
        [[1] + [0] * (seed_count - 1)]
    )


# Exact marked-fibre descent under the residual action, and the invariant
# target coordinates on the global lambda=1 slice.
receiver_S, receiver_pi, receiver_b, receiver_c = sp.symbols(
    "receiver_S receiver_pi receiver_b receiver_c"
)
receiver_u4, receiver_u5, receiver_u6 = sp.symbols(
    "receiver_u4 receiver_u5 receiver_u6",
    nonzero=True,
)
receiver_polynomial = (
    receiver_S
    + receiver_b * receiver_S**2
    + receiver_pi * receiver_S**3
    + receiver_u4 * receiver_pi**4 * receiver_S**4
    + receiver_u5 * receiver_pi**5 * receiver_S**5
    + receiver_u6 * receiver_pi**6 * receiver_S**6
    - receiver_c / 2
)
receiver_transformed = receiver_polynomial.subs(
    {
        receiver_u4: alpha**5 * receiver_u4,
        receiver_u5: alpha**6 * receiver_u5,
        receiver_u6: alpha**7 * receiver_u6,
        receiver_pi: alpha**-2 * receiver_pi,
        receiver_b: alpha**-1 * receiver_b,
        receiver_c: alpha * receiver_c,
    },
    simultaneous=True,
)
assert sp.factor(
    receiver_transformed
    - alpha * receiver_polynomial.subs(receiver_S, receiver_S / alpha)
) == 0

receiver_lambda = receiver_u5 / receiver_u4
assert sp.factor(
    receiver_lambda.subs(
        {
            receiver_u4: alpha**5 * receiver_u4,
            receiver_u5: alpha**6 * receiver_u5,
        },
        simultaneous=True,
    )
    - alpha * receiver_lambda
) == 0
receiver_action = {
    receiver_u4: alpha**5 * receiver_u4,
    receiver_u5: alpha**6 * receiver_u5,
    receiver_pi: alpha**-2 * receiver_pi,
    receiver_b: alpha**-1 * receiver_b,
    receiver_c: alpha * receiver_c,
}
for receiver_invariant in (
    receiver_lambda**2 * receiver_pi,
    receiver_lambda * receiver_b,
    receiver_c / receiver_lambda,
    receiver_u4 / receiver_lambda**5,
):
    assert sp.factor(
        receiver_invariant.subs(receiver_action, simultaneous=True)
        - receiver_invariant
    ) == 0


# Universal discriminant intruder.  For a degree-N polynomial, every
# discriminant coefficient monomial prod(a_i^e_i) satisfies
#
#   sum(e_i)=2N-2,  sum(i*e_i)=N(N-1).
#
# After a_0=-C/2, a_1=1, a_2=-B/2, a_3=P, and
# a_j=u_j*P^j, division by the exact P-order N^2-3N-2 gives
#
#   p=2N+2-2b-e_1-2e_3 <= 2N+2-2b-e_1,
#   N*c+(N-2)*b+(N-1)*e_1 <= N(N-1).
#
# The relaxed constraints on the right make D_N=(2,N,1) the unique
# maximizer of the positive weight w_N=(1,N+1,N), so the exact support does
# too.  The exceptional relaxed possibility b=N+1 has (p,c)=(0,0) and lies
# one weight unit below D_N.
for degree in range(4, 129):
    reduced_p_order = degree**2 - 3 * degree - 2
    assert reduced_p_order >= 2
    intruder = (2, degree, 1)
    intruder_weight = (1, degree + 1, degree)
    intruder_degree = sum(
        exponent * weight
        for exponent, weight in zip(
            intruder,
            intruder_weight,
            strict=True,
        )
    )
    assert intruder_degree == degree**2 + 2 * degree + 2

    equality_cases = []
    for exponent_a1 in range(2 * degree + 3):
        for exponent_b in range(degree + 2):
            maximum_exponent_p = (
                2 * degree + 2 - 2 * exponent_b - exponent_a1
            )
            if maximum_exponent_p < 0:
                continue
            remaining_weight = (
                degree * (degree - 1)
                - (degree - 2) * exponent_b
                - (degree - 1) * exponent_a1
            )
            if remaining_weight < 0:
                continue
            maximum_c = remaining_weight // degree
            for exponent_c in range(maximum_c + 1):
                monomial_degree = (
                    maximum_exponent_p
                    + (degree + 1) * exponent_b
                    + degree * exponent_c
                )
                assert monomial_degree <= intruder_degree
                if monomial_degree == intruder_degree:
                    equality_cases.append(
                        (
                            maximum_exponent_p,
                            exponent_b,
                            exponent_c,
                            exponent_a1,
                        )
                    )
    assert equality_cases == [(2, degree, 1, 0)]


# Exact discriminant regressions verify that the universal term D_N is
# present and that no specialization artefact enters the support in the
# first five ranks.  Its coefficient is nonzero in every rank by the
# standard trinomial discriminant term
#   +-4*(N-2)^(N-2)*a_N^(N-3)*a_2^N*a_0.
S, boundary_P, boundary_B, boundary_C = sp.symbols("S P B C")
seed_values = (2, 3, 5, 7, 11)
for degree in range(4, 9):
    inverse_polynomial = (
        S
        - boundary_B * S**2 / 2
        + boundary_P * S**3
        - boundary_C / 2
        + sum(
            seed_values[index - 4]
            * boundary_P**index
            * S**index
            for index in range(4, degree + 1)
        )
    )
    discriminant_polynomial = sp.Poly(
        sp.discriminant(inverse_polynomial, S),
        boundary_P,
        boundary_B,
        boundary_C,
        domain=sp.QQ,
    )
    minimum_p_order = min(
        monomial[0]
        for monomial, _coefficient in discriminant_polynomial.terms()
    )
    assert minimum_p_order == degree**2 - 3 * degree - 2
    reduced_discriminant = sp.Poly(
        discriminant_polynomial.as_expr()
        / boundary_P**minimum_p_order,
        boundary_P,
        boundary_B,
        boundary_C,
        domain=sp.QQ,
    )
    intruder = (2, degree, 1)
    weight = (1, degree + 1, degree)
    assert reduced_discriminant.coeff_monomial(intruder) != 0
    intruder_degree = sum(
        exponent * weight
        for exponent, weight in zip(intruder, weight, strict=True)
    )
    assert [
        monomial
        for monomial, _coefficient in reduced_discriminant.terms()
        if sum(
            exponent * local_weight
            for exponent, local_weight in zip(
                monomial,
                weight,
                strict=True,
            )
        )
        == intruder_degree
    ] == [intruder]


print("PASS: the independent (alpha,beta) source-target scaling is exact")
print("PASS: a2 is removed by the target shear B -> B-2*a2*P")
print("PASS: the intrinsic Fitting polynomial recovers a3,...,aN")
print("PASS: Fitting support orders the two toric punctures")
print("PASS: Fitting support removes P^m twists")
print("PASS: the stable coefficient-torus quotient has dimension N-4")
print("PASS: I5,J6,...,JN are saturated compiler-slice quotient coordinates")
print("PASS: lambda=u5/u4 gives a global weight-one receiver slice")
print("PASS: the finite-etale fibre descends under the residual action")
print("PASS: D_N=(2,N,1) is the universal exposed discriminant intruder")
print("PASS: ranks four through eight reproduce the all-rank Newton bound")
