#!/usr/bin/env python3
"""Verify the first discriminant obstruction for quadratic two-shear times."""

from __future__ import annotations

import itertools

import sympy as sp


C0, C1, C2, C3 = sp.symbols("C0 C1 C2 C3")
variables = (C0, C1, C2, C3)


def monomial_basis(degree: int) -> list[sp.Expr]:
    return [
        sp.prod(variable**exponent for variable, exponent in zip(variables, exponents))
        for exponents in itertools.product(range(degree + 1), repeat=4)
        if sum(exponents) == degree
    ]


def d_plus(polynomial: sp.Expr) -> sp.Expr:
    return sp.expand(
        3 * C0 * sp.diff(polynomial, C1)
        + 2 * C1 * sp.diff(polynomial, C2)
        + C2 * sp.diff(polynomial, C3)
    )


def d_minus(polynomial: sp.Expr) -> sp.Expr:
    return sp.expand(
        C1 * sp.diff(polynomial, C0)
        + 2 * C2 * sp.diff(polynomial, C1)
        + 3 * C3 * sp.diff(polynomial, C2)
    )


def h_weight(polynomial: sp.Expr) -> sp.Expr:
    return sp.expand(
        3 * C0 * sp.diff(polynomial, C0)
        + C1 * sp.diff(polynomial, C1)
        - C2 * sp.diff(polynomial, C2)
        - 3 * C3 * sp.diff(polynomial, C3)
    )


def coefficient_vector(polynomial: sp.Expr, basis: list[sp.Expr]) -> sp.Matrix:
    poly = sp.Poly(sp.expand(polynomial), *variables)
    return sp.Matrix(
        [poly.coeff_monomial(monomial) for monomial in basis]
    )


def derivation_matrix(
    derivation,
    basis: list[sp.Expr],
) -> sp.Matrix:
    return sp.Matrix.hstack(
        *[
            coefficient_vector(derivation(monomial), basis)
            for monomial in basis
        ]
    )


basis2 = monomial_basis(2)
d_plus2 = derivation_matrix(d_plus, basis2)
d_minus2 = derivation_matrix(d_minus, basis2)
linear_cancellation2 = d_plus2.row_join(d_minus2)
kernel2 = linear_cancellation2.nullspace()
assert len(basis2) == 10
assert len(kernel2) == 10

parameters = sp.symbols(f"u0:{len(kernel2)}")
kernel_vector = sum(
    (
        parameter * vector
        for parameter, vector in zip(parameters, kernel2)
    ),
    sp.zeros(2 * len(basis2), 1),
)
p2 = sp.expand(
    sum(
        kernel_vector[index] * basis2[index]
        for index in range(len(basis2))
    )
)
q2 = sp.expand(
    sum(
        kernel_vector[len(basis2) + index] * basis2[index]
        for index in range(len(basis2))
    )
)
assert sp.expand(d_plus(p2) + d_minus(q2)) == 0

basis4 = monomial_basis(4)
d_plus4 = derivation_matrix(d_plus, basis4)
d_minus4 = derivation_matrix(d_minus, basis4)
linear_cancellation4 = d_plus4.row_join(d_minus4)
invariant_functionals4 = linear_cancellation4.T.nullspace()
assert len(invariant_functionals4) == 1
invariant_functional4 = invariant_functionals4[0]

discriminant = (
    C1**2 * C2**2
    - 4 * C0 * C2**3
    - 4 * C1**3 * C3
    - 27 * C0**2 * C3**2
    + 18 * C0 * C1 * C2 * C3
)
discriminant_vector = coefficient_vector(discriminant, basis4)
normalization = (invariant_functional4.T * discriminant_vector)[0]
assert normalization != 0
invariant_functional4 /= normalization

interaction4 = sp.expand(
    -q2 * h_weight(p2)
    + d_plus(p2) * d_minus(q2)
    - d_minus(p2) * d_plus(q2)
)
obstruction4 = sp.factor(
    (
        invariant_functional4.T
        * coefficient_vector(interaction4, basis4)
    )[0]
)

obstruction_matrix, _ = sp.linear_eq_to_matrix(
    [sp.diff(obstruction4, parameter) / 2 for parameter in parameters],
    parameters,
)
obstruction_rank = obstruction_matrix.rank()


def homogeneous_part(polynomial: sp.Expr, degree: int) -> sp.Expr:
    expanded = sp.Poly(sp.expand(polynomial), *variables)
    return sp.expand(
        sum(
            coefficient
            * sp.prod(
                variable**exponent
                for variable, exponent in zip(variables, monomial)
            )
            for monomial, coefficient in expanded.terms()
            if sum(monomial) == degree
        )
    )


def jacobian_minus_one(p_time: sp.Expr, q_time: sp.Expr) -> sp.Expr:
    x_p = sp.expand(
        d_plus(p_time)
        - q_time * h_weight(p_time)
        - q_time**2 * d_minus(p_time)
    )
    x_q = sp.expand(
        d_plus(q_time)
        - q_time * h_weight(q_time)
        - q_time**2 * d_minus(q_time)
    )
    return sp.expand(
        x_p
        + d_minus(q_time)
        + x_p * d_minus(q_time)
        - d_minus(p_time) * x_q
    )


def zero_free_parameters(solution: sp.Matrix, free_parameters: sp.Matrix) -> sp.Matrix:
    substitutions = {
        parameter: 0
        for parameter in free_parameters
    }
    return solution.subs(substitutions)


def first_recursive_obstruction(
    initial_vector: sp.Matrix,
    maximum_degree: int = 8,
) -> tuple[int | None, sp.Expr]:
    p_time = sp.expand(
        sum(
            initial_vector[index] * basis2[index]
            for index in range(len(basis2))
        )
    )
    q_time = sp.expand(
        sum(
            initial_vector[len(basis2) + index] * basis2[index]
            for index in range(len(basis2))
        )
    )
    assert homogeneous_part(jacobian_minus_one(p_time, q_time), 2) == 0

    for degree in range(3, maximum_degree + 1):
        remainder = homogeneous_part(
            jacobian_minus_one(p_time, q_time),
            degree,
        )
        if remainder == 0:
            continue
        basis = monomial_basis(degree)
        linear_map = derivation_matrix(d_plus, basis).row_join(
            derivation_matrix(d_minus, basis)
        )
        right_hand_side = -coefficient_vector(remainder, basis)
        try:
            solution, free_parameters = linear_map.gauss_jordan_solve(
                right_hand_side
            )
        except ValueError:
            invariant_functionals = linear_map.T.nullspace()
            assert len(invariant_functionals) == 1
            obstruction = sp.factor(
                (
                    invariant_functionals[0].T
                    * coefficient_vector(remainder, basis)
                )[0]
            )
            return degree, obstruction
        solution = zero_free_parameters(solution, free_parameters)
        split = len(basis)
        p_time += sp.expand(
            sum(solution[index] * basis[index] for index in range(split))
        )
        q_time += sp.expand(
            sum(
                solution[split + index] * basis[index]
                for index in range(split)
            )
        )

    return None, sp.Integer(0)


recursive_results = []
for index, kernel_basis_vector in enumerate(kernel2):
    if (
        d_plus(
            sum(
                kernel_basis_vector[position] * basis2[position]
                for position in range(len(basis2))
            )
        )
        == 0
        or d_minus(
            sum(
                kernel_basis_vector[len(basis2) + position] * basis2[position]
                for position in range(len(basis2))
            )
        )
        == 0
    ):
        continue
    recursive_results.append(
        (index, first_recursive_obstruction(kernel_basis_vector))
    )

print(f"quadratic cancellation-kernel dimension: {len(kernel2)}")
print(f"degree-four invariant-cokernel dimension: {len(invariant_functionals4)}")
print(f"first obstruction: {obstruction4}")
print(f"quadratic-form rank: {obstruction_rank}")
print(f"coupled kernel-basis recursion through degree eight: {recursive_results}")
