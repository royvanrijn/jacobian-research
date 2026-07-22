#!/usr/bin/env python3
"""Rank-compress a BCW quadratic-cubic map before homogenization.

The input is a map K=X+Q+C over QQ together with a finite collision.  If the
coefficient matrix of the cubic vector C has row rank k, this module computes
C=Bc with c consisting of k independent component polynomials and constructs

    V(X,Y,T)=(X,Y,T)+(T Q(X)+T^2 B Y, -c(X), 0).

Thus V has n+k+1 variables, rather than the 2n+1 variables in the usual
nilpotent doubling.  The functions here are deliberately independent of the
particular 16-variable shared-factor trace.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class RankFactorization:
    """A rational factorization C=Bc of a cubic component vector."""

    monomials: tuple[tuple[int, ...], ...]
    basis_components: tuple[int, ...]
    B: sp.Matrix
    c: tuple[sp.Poly, ...]


@dataclass(frozen=True)
class EssentialQuotient:
    """Constant-kernel quotient of a cubic-homogeneous map."""

    ambient_variables: tuple[sp.Symbol, ...]
    ambient_h: tuple[sp.Poly, ...]
    kernel: sp.Matrix
    B: sp.Matrix
    C: sp.Matrix
    quotient_variables: tuple[sp.Symbol, ...]
    quotient_h: tuple[sp.Poly, ...]


def monomial_expression(exponents: tuple[int, ...], variables: list[sp.Symbol]) -> sp.Expr:
    return sp.prod(variable**exponent for variable, exponent in zip(variables, exponents))


def homogeneous_part(expression: sp.Expr, variables: list[sp.Symbol], degree: int) -> sp.Poly:
    poly = sp.Poly(expression, *variables, domain=sp.QQ)
    selected = sum(
        coefficient * monomial_expression(exponents, variables)
        for exponents, coefficient in poly.terms()
        if sum(exponents) == degree
    )
    return sp.Poly(selected, *variables, domain=sp.QQ)


def extract_quadratic_cubic(
    expressions: list[sp.Expr], variables: list[sp.Symbol]
) -> tuple[list[sp.Poly], list[sp.Poly]]:
    """Extract Q and C and certify K=X+Q+C with identity linear part."""
    if len(expressions) != len(variables):
        raise ValueError("K must be a square polynomial map")
    quadratic = [homogeneous_part(expression, variables, 2) for expression in expressions]
    cubic = [homogeneous_part(expression, variables, 3) for expression in expressions]
    for index, (expression, q_part, c_part) in enumerate(zip(expressions, quadratic, cubic)):
        remainder = sp.Poly(
            sp.expand(expression - variables[index] - q_part.as_expr() - c_part.as_expr()),
            *variables,
            domain=sp.QQ,
        )
        if not remainder.is_zero:
            raise ValueError(f"component {index} is not X_i+Q_i+C_i")
    return quadratic, cubic


def factor_cubic_output(cubic: list[sp.Poly]) -> RankFactorization:
    """Compute C=Bc over QQ using independent component rows as c."""
    if not cubic:
        return RankFactorization((), (), sp.zeros(0, 0), ())
    variables = list(cubic[0].gens)
    monomials = tuple(
        sorted(
            {
                exponents
                for poly in cubic
                for exponents, coefficient in poly.terms()
                if coefficient
            }
        )
    )
    coefficient_matrix = sp.Matrix(
        [[poly.coeff_monomial(exponents) for exponents in monomials] for poly in cubic]
    )
    # Pivot columns of A^T are precisely independent rows of A.
    basis_components = tuple(coefficient_matrix.T.rref()[1])
    basis_matrix = coefficient_matrix[list(basis_components), :]
    k = len(basis_components)
    rows: list[list[sp.Expr]] = []
    for row_index in range(coefficient_matrix.rows):
        if k:
            solution, parameters = basis_matrix.T.gauss_jordan_solve(
                coefficient_matrix[row_index, :].T
            )
            assert parameters.rows == 0
            rows.append([sp.cancel(value) for value in solution])
        else:
            rows.append([])
    B = sp.Matrix(rows) if k else sp.zeros(len(cubic), 0)
    assert B * basis_matrix == coefficient_matrix
    c = tuple(cubic[index] for index in basis_components)
    assert all(
        sp.Poly(
            cubic[index].as_expr()
            - sum(B[index, j] * c[j].as_expr() for j in range(k)),
            *variables,
            domain=sp.QQ,
        ).is_zero
        for index in range(len(cubic))
    )
    return RankFactorization(monomials, basis_components, B, c)


def verify_parametric_factorization(
    variables: list[sp.Symbol],
    quadratic: list[sp.Poly],
    cubic: list[sp.Poly],
    factorization: RankFactorization,
) -> None:
    """Verify id+tN=P_t o (E_t x id) o A_t component by component."""
    n = len(variables)
    k = len(factorization.c)
    y = list(sp.symbols(f"rank_y0:{k}"))
    t = sp.Symbol("rank_t")
    B = factorization.B
    c = factorization.c

    e_first = [
        variables[i] + t * quadratic[i].as_expr() + t**2 * cubic[i].as_expr()
        for i in range(n)
    ]
    a_second = [y[j] - t * c[j].as_expr() for j in range(k)]
    right_first = [
        e_first[i] + t * sum(B[i, j] * a_second[j] for j in range(k))
        for i in range(n)
    ]
    left_first = [
        variables[i]
        + t
        * (
            quadratic[i].as_expr()
            + sum(B[i, j] * y[j] for j in range(k))
        )
        for i in range(n)
    ]
    assert all(sp.expand(left - right) == 0 for left, right in zip(left_first, right_first))

    # Verify E_t(X)=t^{-1}K(tX) as a polynomial identity.
    scaling = {variable: t * variable for variable in variables}
    scaled_k = [
        sp.cancel(
            (
                variables[i]
                + quadratic[i].as_expr()
                + cubic[i].as_expr()
            ).subs(scaling)
            / t
        )
        for i in range(n)
    ]
    assert all(sp.expand(left - right) == 0 for left, right in zip(e_first, scaled_k))

    # The upper-left Jacobian block of the final homogeneous map has Schur
    # complement DE_t.  Since JC=B Jc, this also checks the determinant bridge.
    jq = sp.Matrix([poly.as_expr() for poly in quadratic]).jacobian(variables)
    jc = sp.Matrix([poly.as_expr() for poly in c]).jacobian(variables)
    de = sp.Matrix(e_first).jacobian(variables)
    assert de == sp.eye(n) + t * jq + t**2 * B * jc


def rank_compressed_homogeneous_map(
    variables: list[sp.Symbol],
    quadratic: list[sp.Poly],
    factorization: RankFactorization,
) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Poly, ...]]:
    """Construct the homogeneous part of the rank-compressed BCW map."""

    n = len(variables)
    k = len(factorization.c)
    y = list(sp.symbols(f"essential_y0:{k}"))
    t = sp.Symbol("essential_t")
    all_variables = tuple(variables + y + [t])
    B = factorization.B
    first = [
        sp.Poly(
            t * quadratic[i].as_expr()
            + t**2 * sum(B[i, j] * y[j] for j in range(k)),
            *all_variables,
            domain=sp.QQ,
        )
        for i in range(n)
    ]
    second = [
        sp.Poly(-poly.as_expr(), *all_variables, domain=sp.QQ)
        for poly in factorization.c
    ]
    zero = sp.Poly(0, *all_variables, domain=sp.QQ)
    homogeneous = tuple(first + second + [zero])
    assert all(
        sum(exponents) == 3
        for poly in homogeneous
        for exponents, coefficient in poly.terms()
        if coefficient
    )
    return all_variables, homogeneous


def constant_jacobian_kernel(
    components: tuple[sp.Poly, ...], variables: tuple[sp.Symbol, ...]
) -> sp.Matrix:
    """Compute the constant intersection of the kernels of ``JH(x)``."""

    jacobian = sp.Matrix([poly.as_expr() for poly in components]).jacobian(variables)
    monomials = sorted(
        {
            exponents
            for entry in jacobian
            for exponents, coefficient in sp.Poly(entry, *variables, domain=sp.QQ).terms()
            if coefficient
        }
    )
    coefficient_matrix = sp.Matrix(
        [
            [
                sp.Poly(jacobian[i, j], *variables, domain=sp.QQ).coeff_monomial(exponents)
                for j in range(len(variables))
            ]
            for i in range(len(components))
            for exponents in monomials
        ]
    )
    basis = coefficient_matrix.nullspace()
    return sp.Matrix.hstack(*basis) if basis else sp.zeros(len(variables), 0)


def constant_kernel_quotient(
    variables: tuple[sp.Symbol, ...], components: tuple[sp.Poly, ...]
) -> EssentialQuotient:
    """Construct ``B,C`` and the essential-input quotient of ``I+H``."""

    dimension = len(variables)
    kernel = constant_jacobian_kernel(components, variables)
    quotient_dimension = dimension - kernel.cols
    if kernel.cols:
        B = sp.Matrix.hstack(*kernel.T.nullspace()).T
    else:
        B = sp.eye(dimension)
    assert B.shape == (quotient_dimension, dimension)
    assert B * kernel == sp.zeros(quotient_dimension, kernel.cols)
    pivots = list(B.rref()[1])
    C = sp.zeros(dimension, quotient_dimension)
    pivot_inverse = B[:, pivots].inv()
    for local_index, ambient_index in enumerate(pivots):
        C[ambient_index, :] = pivot_inverse[local_index, :]
    assert B * C == sp.eye(quotient_dimension)

    q = tuple(sp.symbols(f"essential_q0:{quotient_dimension}"))
    substitution = dict(zip(variables, list(C * sp.Matrix(q))))
    descended = B * sp.Matrix([poly.as_expr() for poly in components]).subs(
        substitution, simultaneous=True
    )
    quotient_h = tuple(
        sp.Poly(sp.expand(value), *q, domain=sp.QQ) for value in descended
    )
    return EssentialQuotient(
        ambient_variables=variables,
        ambient_h=components,
        kernel=kernel,
        B=B,
        C=C,
        quotient_variables=q,
        quotient_h=quotient_h,
    )


def jacobian_coefficient_matrices(
    components: tuple[sp.Poly, ...], variables: tuple[sp.Symbol, ...]
) -> tuple[sp.Matrix, ...]:
    """Coefficient matrices of the polynomial matrix ``JH``."""

    jacobian = sp.Matrix([poly.as_expr() for poly in components]).jacobian(variables)
    monomials = sorted(
        {
            exponents
            for entry in jacobian
            for exponents, coefficient in sp.Poly(entry, *variables, domain=sp.QQ).terms()
            if coefficient
        }
    )
    return tuple(
        sp.Matrix(
            [
                [
                    sp.Poly(jacobian[i, j], *variables, domain=sp.QQ).coeff_monomial(exponents)
                    for j in range(len(variables))
                ]
                for i in range(len(components))
            ]
        )
        for exponents in monomials
    )


def cyclic_invariant_row_module_dimensions(
    components: tuple[sp.Poly, ...],
    variables: tuple[sp.Symbol, ...],
    prime: int = 1_000_003,
) -> tuple[int, ...]:
    """Good-prime dimensions of coordinate-generated cyclic row modules.

    This deterministic finite-field diagnostic is suitable inside a search;
    shortlisted candidates still require an exact characteristic-zero module
    audit such as ``audit_bcw_22_linear_quotients.py``.
    """

    rational_matrices = jacobian_coefficient_matrices(components, variables)
    dimension = len(variables)
    matrices = [
        tuple(
            tuple(
                int(value.p) * pow(int(value.q), -1, prime) % prime
                for value in matrix.row(i)
            )
            for i in range(dimension)
        )
        for matrix in rational_matrices
    ]
    dimensions: list[int] = []
    for seed_index in range(dimension):
        pivots: dict[int, tuple[int, ...]] = {}
        seed = tuple(1 if index == seed_index else 0 for index in range(dimension))
        pending = [seed]

        def add(row: tuple[int, ...]) -> tuple[int, ...] | None:
            reduced = list(row)
            for pivot in sorted(pivots):
                if reduced[pivot]:
                    factor = reduced[pivot]
                    reduced = [
                        (left - factor * right) % prime
                        for left, right in zip(reduced, pivots[pivot])
                    ]
            pivot = next((index for index, value in enumerate(reduced) if value), None)
            if pivot is None:
                return None
            inverse = pow(reduced[pivot], -1, prime)
            normalized = tuple(value * inverse % prime for value in reduced)
            pivots[pivot] = normalized
            return normalized

        while pending:
            row = pending.pop()
            normalized = add(row)
            if normalized is None:
                continue
            pending.extend(
                tuple(
                    sum(normalized[k] * matrix[k][j] for k in range(dimension)) % prime
                    for j in range(dimension)
                )
                for matrix in matrices
            )
        dimensions.append(len(pivots))
    return tuple(sorted(set(dimensions)))
