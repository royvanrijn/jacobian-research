#!/usr/bin/env python3
"""Exact low-degree audit of the failed all-k transfer block conductor-ribbon norm map.

The calculations prove the stated facts only for k=1,2,3.  In particular,
they do not establish flatness over the full completed base for arbitrary k.
"""
from __future__ import annotations

import itertools
import math

import sympy as sp


Z = sp.Symbol("Z")


def reduce_by(basis: sp.GroebnerBasis, expression: sp.Expr) -> sp.Expr:
    return sp.expand(basis.reduce(sp.expand(expression))[1])


def reduce_z_coefficients(basis: sp.GroebnerBasis, expression: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), Z)
    return sp.expand(sum(
        reduce_by(basis, coefficient) * Z**degree[0]
        for degree, coefficient in polynomial.terms()
    ))


def boolean_norm_data(k: int):
    roots = sp.symbols(f"r0:{k}")
    epsilons = sp.symbols(f"e0:{k}")
    boolean = sp.groebner(
        [epsilon**2 for epsilon in epsilons],
        *roots,
        *epsilons,
        order="grevlex",
        domain=sp.QQ,
    )
    U = sp.Integer(1)
    V = sp.Integer(1)
    for root, epsilon in zip(roots, epsilons):
        q = Z - root
        U *= q**3 + sp.Rational(3, 2) * epsilon * q
        V *= q**2 + epsilon
    U = sp.Poly(reduce_z_coefficients(boolean, U), Z)
    V = sp.Poly(reduce_z_coefficients(boolean, V), Z)
    assert all(
        reduce_by(boolean, coefficient) == 0
        for coefficient in sp.Poly(U.as_expr() ** 2 - V.as_expr() ** 3, Z).all_coeffs()
    )
    return roots, epsilons, boolean, U, V


def extracted_square_root(
    k: int,
    variables: tuple[sp.Symbol, ...],
    boolean: sp.GroebnerBasis,
    V: sp.Poly,
):
    coefficients: list[sp.Expr] = []
    S = Z**k
    for index in range(1, k + 1):
        degree = 2 * k - index
        residual = sp.Poly(V.as_expr() - S**2, Z).coeff_monomial(Z**degree)
        coefficient = reduce_by(boolean, residual / 2)
        coefficients.append(coefficient)
        S += coefficient * Z ** (k - index)
    remainder = sp.Poly(reduce_z_coefficients(boolean, V.as_expr() - S**2), Z)
    assert remainder.degree() < k or remainder.is_zero
    transverse = [
        reduce_by(boolean, remainder.coeff_monomial(Z**degree))
        for degree in range(k)
    ]
    return coefficients, transverse


def standard_monomials(basis: sp.GroebnerBasis, variable_count: int):
    leading = [polynomial.LM(order=basis.order).exponents for polynomial in basis.polys]
    # The audited quotients have nilpotency bounded by 2k.  Increase the box
    # until no standard monomial touches its upper face.
    for bound in range(2, 2 * variable_count + 4):
        standard = [
            exponent
            for exponent in itertools.product(range(bound), repeat=variable_count)
            if not any(
                all(exponent[index] >= monomial[index] for index in range(variable_count))
                for monomial in leading
            )
        ]
        if standard and all(max(exponent) < bound - 1 for exponent in standard):
            return standard
    raise AssertionError("standard-monomial box did not stabilize")


def quotient_linear_algebra(basis: sp.GroebnerBasis, variables, standard):
    index = {exponent: position for position, exponent in enumerate(standard)}

    def expression(exponent):
        return sp.prod(variable**power for variable, power in zip(variables, exponent))

    monomials = [expression(exponent) for exponent in standard]

    def vector(value):
        remainder = sp.Poly(reduce_by(basis, value), *variables)
        answer = sp.zeros(len(standard), 1)
        for exponent, coefficient in remainder.terms():
            answer[index[exponent]] = coefficient
        return answer

    return monomials, vector


def column_basis(vectors: list[sp.Matrix]) -> list[sp.Matrix]:
    if not vectors:
        return []
    return sp.Matrix.hstack(*vectors).columnspace()


def generated_subalgebra_dimension(generators, monomials, vector) -> int:
    basis_vectors = [vector(1)]
    while True:
        candidates = list(basis_vectors)
        for generator in generators:
            for current in basis_vectors:
                representative = sum(
                    current[index] * monomial for index, monomial in enumerate(monomials)
                )
                candidates.append(vector(generator * representative))
        enlarged = column_basis(candidates)
        if len(enlarged) == len(basis_vectors):
            return len(enlarged)
        basis_vectors = enlarged


def invariant_dimension(k, roots, epsilons, monomials, vector) -> int:
    size = len(monomials)
    equations = []
    for index in range(k - 1):
        substitution = {
            roots[index]: roots[index + 1],
            roots[index + 1]: roots[index],
            epsilons[index]: epsilons[index + 1],
            epsilons[index + 1]: epsilons[index],
        }
        action = sp.Matrix.hstack(
            *[vector(monomial.xreplace(substitution)) for monomial in monomials]
        )
        equations.append(action - sp.eye(size))
    if not equations:
        return size
    return size - sp.Matrix.vstack(*equations).rank()


def factorization_fiber_length(k: int) -> int:
    transverse = sp.symbols(f"a0:{k}")
    u = sp.symbols(f"u0:{3 * k}")
    V = Z ** (2 * k) + sum(transverse[index] * Z**index for index in range(k))
    U = Z ** (3 * k) + sum(u[index] * Z**index for index in range(3 * k))
    difference = sp.Poly(sp.expand(U**2 - V**3), Z)
    solution = {}
    for degree, variable in zip(range(6 * k - 1, 3 * k - 1, -1), reversed(u)):
        equation = sp.expand(difference.coeff_monomial(Z**degree).subs(solution))
        coefficient = sp.diff(equation, variable)
        assert coefficient != 0 and not coefficient.has(variable)
        solution[variable] = sp.expand(-(equation - coefficient * variable) / coefficient)
    equations = [
        sp.expand(difference.coeff_monomial(Z**degree).subs(solution))
        for degree in range(3 * k)
    ]
    basis = sp.groebner(equations, *transverse, order="grevlex", domain=sp.QQ)
    standard = standard_monomials(basis, k)
    return len(standard)


def audit(k: int):
    roots, epsilons, boolean, _U, V = boolean_norm_data(k)
    variables = roots + epsilons
    base, transverse = extracted_square_root(k, variables, boolean, V)

    if k == 2:
        natural_p = -(roots[0] + roots[1])
        natural_q = roots[0] * roots[1]
        assert reduce_by(boolean, base[0] - natural_p) == 0
        assert reduce_by(
            boolean,
            base[1] - natural_q - (epsilons[0] + epsilons[1]) / 2,
        ) == 0

    collision_ideal = [epsilon**2 for epsilon in epsilons] + base
    collision = sp.groebner(
        collision_ideal,
        *variables,
        order="grevlex",
        domain=sp.QQ,
    )
    standard = standard_monomials(collision, len(variables))
    monomials, vector = quotient_linear_algebra(collision, variables, standard)
    ordered_length = len(standard)
    assert ordered_length == math.factorial(k) * 2**k

    norm_dimension = generated_subalgebra_dimension(transverse, monomials, vector)
    fixed_dimension = invariant_dimension(k, roots, epsilons, monomials, vector)
    factor_length = factorization_fiber_length(k)
    assert fixed_dimension == factor_length == 2**k
    return ordered_length, norm_dimension, fixed_dimension, factor_length


def main() -> None:
    results = {}
    for k in (1, 2, 3):
        result = audit(k)
        results[k] = result
        ordered_length, norm_dimension, invariant_dimension_value, factor_length = result
        print(
            f"PASS k={k}: ordered collision length={ordered_length}, "
            f"norm={norm_dimension}, invariant={invariant_dimension_value}, "
            f"factorization={factor_length}"
        )
    assert results[1][1:] == (2, 2, 2)
    assert results[2][1:] == (3, 4, 4)

    t = sp.Symbol("t")
    adversarial_V = Z**4 + t * Z
    adversarial_U = Z**6 + sp.Rational(3, 2) * t * Z**3 + sp.Rational(3, 8) * t**2
    difference = sp.rem(
        sp.Poly(sp.expand(adversarial_U**2 - adversarial_V**3), t),
        sp.Poly(t**3, t),
    ).as_expr()
    assert difference == 0
    print("PASS: k=2 norm base differs from the forgetful symmetric-product base")
    print("PASS: k=2 norm misses one invariant/factorization direction")
    print("PASS: the explicit K[t]/(t^3) factorization deformation is exact")


if __name__ == "__main__":
    main()
