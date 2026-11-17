#!/usr/bin/env python3
"""Exact regressions for optimal weighted-bridge moment coordinates.

For g=u*h(g), the bridge gives

    [u^m]g = E(Q P^m)/m!.

This script constructs those coefficients by univariate Lagrange inversion,
reverts g without using the identity being tested, and verifies h=z/eta.
The written argument in WEIGHTED_GAUSSIAN_BRIDGE.md is the all-order proof.
Besides full series reversion, this script checks that the N-3 moments
M_3,...,M_(N-1) are triangular coordinates on the normalized degree-N seed
slice, verifies the variable-scale N-2 coordinate version, and checks the
Toeplitz--Hessenberg determinant formula for the reciprocal series.
"""

from __future__ import annotations

import sympy as sp


u, z = sp.symbols("u z")


def truncate(expr: sp.Expr, variable: sp.Symbol, order: int) -> sp.Expr:
    return sp.expand(sp.series(expr, variable, 0, order).removeO())


def mixed_generating_series(h: sp.Expr, degree: int) -> sp.Expr:
    """Return g through u^(degree+1) from the exact mixed-moment formula."""
    result = sp.Integer(0)
    for m in range(1, degree + 2):
        coefficient = sp.expand(h**m).coeff(z, m - 1) / sp.Integer(m)
        result += coefficient * u**m
    return sp.expand(result)


def normalized_moment(h: sp.Expr, m: int) -> sp.Expr:
    """Return mu_m=M_m/m!=[u^m]g by Lagrange inversion."""
    return sp.expand(h**m).coeff(z, m - 1) / sp.Integer(m)


def recover_triangular_coefficients(
    moments: dict[int, sp.Expr], last_degree: int
) -> dict[int, sp.Expr]:
    """Recover a_1,...,a_last from mu_2,... by triangular elimination."""
    recovered: dict[int, sp.Expr] = {1: sp.Integer(0)}
    for degree in range(2, last_degree + 1):
        trial = sp.symbols(f"trial_a_{degree}")
        partial_h = 1 + sum(
            recovered[index] * z**index for index in range(1, degree)
        ) + trial * z**degree
        equation = sp.expand(
            normalized_moment(partial_h, degree + 1) - moments[degree + 1]
        )
        solution = sp.solve(sp.Eq(equation, 0), trial)
        assert len(solution) == 1
        recovered[degree] = sp.factor(solution[0])
    return recovered


def endpoint_completion(
    coefficients: dict[int, sp.Expr], degree: int, derivative_at_one: sp.Expr
) -> dict[int, sp.Expr]:
    """Use h(1)=1 and h'(1)=derivative_at_one for the last two terms."""
    first_sum = sum(coefficients[index] for index in range(1, degree - 1))
    derivative_sum = sum(
        index * coefficients[index] for index in range(1, degree - 1)
    )
    total = -first_sum
    weighted = derivative_at_one - derivative_sum
    completed = dict(coefficients)
    completed[degree - 1] = sp.expand(degree * total - weighted)
    completed[degree] = sp.expand(weighted - (degree - 1) * total)
    return completed


def toeplitz_hessenberg(values: list[sp.Expr], size: int) -> sp.Matrix:
    """Return T_size(e), where values[n]=e_n and e_0=1."""
    return sp.Matrix(
        size,
        size,
        lambda row, column: (
            values[row - column + 1]
            if column <= row
            else sp.Integer(1)
            if column == row + 1
            else sp.Integer(0)
        ),
    )


def compositional_inverse(g: sp.Expr, order: int) -> sp.Expr:
    """Revert g=u+O(u^2) through z^(order-1), coefficient by coefficient."""
    assert sp.expand(g).coeff(u, 1) == 1
    eta = z
    for n in range(2, order):
        trial_coefficient = sp.symbols(f"c{n}")
        trial = eta + trial_coefficient * z**n
        composed = truncate(g.subs(u, trial), z, n + 1)
        equation = sp.expand(composed - z).coeff(z, n)
        solution = sp.solve(sp.Eq(equation, 0), trial_coefficient)
        assert len(solution) == 1
        eta = sp.expand(trial.subs(trial_coefficient, solution[0]))
    return eta


def recover(h: sp.Expr, degree: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    g = mixed_generating_series(h, degree)
    eta = compositional_inverse(g, degree + 2)
    recovered = truncate(z / eta, z, degree + 1)
    return g, eta, recovered


def main() -> None:
    # A universal degree-four audit: M_1,...,M_5 recover every coefficient.
    a1, a2, a3, a4 = sp.symbols("a1 a2 a3 a4")
    symbolic_h = 1 + a1 * z + a2 * z**2 + a3 * z**3 + a4 * z**4
    _, symbolic_eta, symbolic_recovered = recover(symbolic_h, 4)
    assert sp.expand(symbolic_recovered - symbolic_h) == 0
    assert truncate(symbolic_eta * symbolic_h - z, z, 6) == 0
    print("PASS Gaussian fingerprint: symbolic degree four recovered from M_1,...,M_5")

    # A concrete degree-five weighted seed with rational scale.
    seed = z**4 * (1 - z)
    normalized_seed = sp.expand(-seed / sp.diff(seed, z).subs(z, 1))
    scale = sp.Rational(2, 3)
    assert sp.diff(normalized_seed, z).subs(z, 1) == -1
    concrete_h = 1 + scale * normalized_seed
    concrete_g, concrete_eta, concrete_recovered = recover(concrete_h, 5)
    assert sp.expand(concrete_recovered - concrete_h) == 0
    assert truncate(concrete_g.subs(u, concrete_eta) - z, z, 7) == 0
    recovered_scale = -sp.diff(concrete_recovered, z).subs(z, 1)
    recovered_seed = sp.cancel((concrete_recovered - 1) / recovered_scale)
    assert recovered_scale == scale
    assert sp.expand(recovered_seed - normalized_seed) == 0
    print("PASS Gaussian fingerprint: weighted degree five recovered from M_1,...,M_6")

    # The normalized degree-N slice has free coordinates c_2,...,c_(N-2).
    # Its N-3 optimal moments have triangular Jacobian lambda^(N-3), and the
    # two endpoint equations recover c_(N-1),c_N.  Check several universal
    # degrees with symbolic free coordinates and symbolic nonzero scale.
    scale_symbol = sp.symbols("lambda", nonzero=True)
    degree_six_h = None
    for degree in range(4, 9):
        free = sp.symbols(f"c2:{degree - 1}")
        seed_coefficients = {
            index: free[index - 2] for index in range(2, degree - 1)
        }
        completed_seed = endpoint_completion(
            {1: sp.Integer(0), **seed_coefficients}, degree, sp.Integer(-1)
        )
        normalized_h = 1 + scale_symbol * sum(
            completed_seed[index] * z**index for index in range(2, degree + 1)
        )
        minimal_moments = {
            moment_index: normalized_moment(normalized_h, moment_index)
            for moment_index in range(3, degree)
        }
        jacobian = sp.Matrix(list(minimal_moments.values())).jacobian(free)
        assert sp.factor(jacobian.det() - scale_symbol ** (degree - 3)) == 0

        recovered_scaled = recover_triangular_coefficients(
            minimal_moments, degree - 2
        )
        recovered_seed = {
            index: sp.cancel(recovered_scaled[index] / scale_symbol)
            for index in range(1, degree - 1)
        }
        recovered_seed = endpoint_completion(
            recovered_seed, degree, sp.Integer(-1)
        )
        for index in range(2, degree + 1):
            assert sp.expand(
                recovered_seed[index] - completed_seed[index]
            ) == 0
        if degree == 6:
            degree_six_h = normalized_h
    print("PASS optimal fingerprint: M_3,...,M_(N-1) are N-3 uniform coordinates")

    # When the scale varies, M_3,...,M_N recover a_2,...,a_(N-1), h(1)=1
    # recovers a_N, and -h'(1) recovers lambda.
    degree = 6
    assert degree_six_h is not None
    variable_scale_moments = {
        moment_index: normalized_moment(degree_six_h, moment_index)
        for moment_index in range(3, degree + 1)
    }
    recovered_variable = recover_triangular_coefficients(
        variable_scale_moments, degree - 1
    )
    recovered_variable[degree] = -sum(
        recovered_variable[index] for index in range(1, degree)
    )
    recovered_h = 1 + sum(
        recovered_variable[index] * z**index
        for index in range(1, degree + 1)
    )
    recovered_scale = -sp.diff(recovered_h, z).subs(z, 1)
    assert sp.factor(recovered_scale - scale_symbol) == 0
    assert sp.expand(recovered_h - degree_six_h) == 0
    print("PASS optimal fingerprint: variable scale uses exactly N-2 moments")

    # If eta/z=e=sum e_n z^n, then a_n=[z^n](1/e)=(-1)^n det T_n(e).
    b1, b2, b3, b4 = sp.symbols("b1 b2 b3 b4")
    determinant_h = 1 + b1 * z + b2 * z**2 + b3 * z**3 + b4 * z**4
    reciprocal = sp.series(1 / determinant_h, z, 0, 8).removeO().expand()
    reciprocal_coefficients = [
        reciprocal.coeff(z, index) for index in range(8)
    ]
    for index in range(1, 8):
        determinant_coefficient = sp.expand(
            (-1) ** index
            * toeplitz_hessenberg(reciprocal_coefficients, index).det()
        )
        expected = determinant_h.coeff(z, index)
        assert sp.factor(determinant_coefficient - expected) == 0
    print("PASS moment equations: Toeplitz-Hessenberg determinants recover h")


if __name__ == "__main__":
    main()
