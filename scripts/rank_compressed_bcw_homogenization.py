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
