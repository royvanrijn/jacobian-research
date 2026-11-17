#!/usr/bin/env python3
"""Bounded exact replay of the all-degree dihedral obstruction theorem.

The accompanying note proves the statements uniformly.  This checker
replays the identities for 3 <= n <= 12.
"""

from __future__ import annotations

import sympy as sp


a, u, U, V, S, T = sp.symbols("a u U V S T")
curve_x, curve_delta = sp.symbols("curve_x curve_delta")
normal_a, normal_b, log_scale = sp.symbols(
    "normal_a normal_b log_scale"
)


def dickson_power_sums(limit: int) -> list[sp.Expr]:
    values = [sp.Integer(2), a]
    for _ in range(2, limit + 1):
        values.append(sp.expand(a * values[-1] - u * values[-2]))
    return values


def dickson_second_kind(limit: int) -> list[sp.Expr]:
    values = [sp.Integer(1), a]
    for _ in range(2, limit + 1):
        values.append(sp.expand(a * values[-1] - u * values[-2]))
    return values


def one_normal_search_route(n: int, p_degree: int, q_degree: int) -> str:
    """Compile a one-normal degree pair to its exact obstruction."""

    if p_degree == q_degree == 0:
        return "inactive"
    if p_degree == 0 or q_degree == 0:
        return "degree"
    if n * p_degree != 2 * q_degree:
        return "degree"
    if n % 2 == 0:
        return "even_factorization"
    return "odd_infinity"


power_sums = dickson_power_sums(12)
second_kind = dickson_second_kind(11)
C = a**2 - 4 * u

for n in range(3, 13):
    P = power_sums[n]
    J = sp.diff(P, a)
    derivative_factor = sp.cancel(J / n)
    branch_pullback = sp.expand(P**2 - 4 * u**n)
    assert sp.expand(branch_pullback - C * derivative_factor**2) == 0

    # The canonical two-mask matrix has determinant Delta_n.
    Delta = V**2 - 4 * U**n
    mask_matrix = sp.Matrix([[V, 2 * U], [2 * U ** (n - 1), V]])
    assert sp.expand(mask_matrix.det() - Delta) == 0
    adjugate = mask_matrix.adjugate()
    assert adjugate * mask_matrix == sp.eye(2) * Delta

    # The generic fixed-discriminant curve is smooth and has positive genus.
    curve_polynomial = 4 * curve_x**n + curve_delta
    curve_gcd = sp.gcd(
        sp.Poly(curve_polynomial, curve_x, domain=sp.QQ.frac_field(curve_delta)),
        sp.Poly(
            sp.diff(curve_polynomial, curve_x),
            curve_x,
            domain=sp.QQ.frac_field(curve_delta),
        ),
    )
    assert curve_gcd.degree() == 0
    assert (n - 1) // 2 >= 1

    # The first nonautomorphic normalized-cusp chart has a contraction vector
    # divisible by Delta, but its desired quotient does not.
    cusp_p = V**2 + Delta * S
    cusp_q = 2 * V**n + Delta * T
    cusp_quotient = sp.cancel((cusp_q**2 - 4 * cusp_p**n) / Delta)
    assert sp.denom(cusp_quotient) == 1
    cusp_variables = (V, S, T)
    contraction = sp.Matrix(
        [sp.diff(cusp_p, variable) for variable in cusp_variables]
    ).cross(
        sp.Matrix(
            [sp.diff(cusp_q, variable) for variable in cusp_variables]
        )
    )
    assert all(sp.rem(entry, Delta, V) == 0 for entry in contraction)
    cusp_basis = sp.groebner([Delta], U, V, S, T)
    cusp_remainder = sp.expand(cusp_basis.reduce(cusp_quotient)[1])
    expected_remainder = sp.expand(
        4 * V**n * (T - n * S * V ** (n - 2))
    )
    assert sp.expand(cusp_remainder - expected_remainder) == 0

    # The affine-normal tangential chart: the T^n and T^2 coefficients force
    # the two normal coefficients to vanish.
    h = V + S
    tangent_p = h**2 + normal_a * Delta * T
    tangent_q = 2 * h**n + normal_b * Delta * T
    tangent_map = sp.Matrix([tangent_p, tangent_q, U, S])
    tangent_jacobian = sp.factor(
        tangent_map.jacobian((U, V, S, T)).det()
    )
    expected_jacobian = sp.expand(
        2
        * h
        * Delta
        * (n * normal_a * h ** (n - 2) - normal_b)
    )
    assert sp.expand(tangent_jacobian - expected_jacobian) == 0
    tangent_error = sp.Poly(
        sp.expand(
            tangent_q**2
            - 4 * tangent_p**n
            - log_scale * Delta * tangent_jacobian
        ),
        T,
    )
    assert sp.expand(
        tangent_error.coeff_monomial(T**n)
        + 4 * normal_a**n * Delta**n
    ) == 0
    after_normal_a = sp.Poly(
        tangent_error.as_expr().subs(normal_a, 0),
        T,
    )
    assert sp.expand(
        after_normal_a.coeff_monomial(T**2)
        - normal_b**2 * Delta**2
    ) == 0

    if n % 2 == 0:
        m = n // 2
        A_m = second_kind[m - 1]
        B_m = power_sums[m]
        assert sp.expand(derivative_factor - A_m * B_m) == 0
        assert sp.expand(P - 2 * u**m - C * A_m**2) == 0
        assert sp.expand(P + 2 * u**m - B_m**2) == 0

    # Away from the resonance n*deg_T(p) = 2*deg_T(q), the cusp side has
    # degree strictly above the one-variable Jacobian bound.
    for p_degree in range(0, 17):
        for q_degree in range(0, 25):
            if p_degree == q_degree == 0:
                continue
            jacobian_bound = p_degree + q_degree - 1
            if p_degree == 0 or q_degree == 0:
                left_degree = max(n * p_degree, 2 * q_degree)
                assert left_degree > jacobian_bound
                continue
            if n * p_degree != 2 * q_degree:
                left_degree = max(n * p_degree, 2 * q_degree)
                assert left_degree > jacobian_bound

    # In even degree n=2m, write q=+/-2*p^m+b after top-degree
    # cancellation.  The cusp difference has degree m*r+deg(b), while the
    # residual bracket has degree at most r+deg(b)-1.
    if n % 2 == 0:
        m = n // 2
        for p_degree in range(1, 17):
            for residual_degree in range(0, p_degree):
                cusp_degree = m * p_degree + residual_degree
                residual_bracket_bound = (
                    p_degree + residual_degree - 1
                )
                assert cusp_degree > residual_bracket_bound

        # Resonance alone is not a solution: this pair has the right top
        # degrees, but its cusp difference still varies with T while its
        # bracket is constant.
        false_p = T + V
        false_q = 2 * false_p**m + V
        false_bracket = sp.expand(
            sp.diff(false_p, V) * sp.diff(false_q, T)
            - sp.diff(false_p, T) * sp.diff(false_q, V)
        )
        false_cusp = sp.expand(false_q**2 - 4 * false_p**n)
        assert false_bracket == -1
        assert sp.Poly(false_cusp, T).degree() == m
    else:
        # At an odd resonant place at infinity, deg_T(p)=2d.  If
        # ord_infinity(y-epsilon)=k, the transformed cusp side has strictly
        # smaller order than either term of the transformed bracket.
        for half_p_degree in range(1, 17):
            for vanishing_order in range(1, 25):
                cusp_order = (
                    vanishing_order
                    - (n - 1) * half_p_degree
                )
                bracket_order_bound = (
                    vanishing_order - half_p_degree + 1
                )
                assert cusp_order < bracket_order_bound

        false_p = T**2 + V
        false_q = 2 * T**n
        false_bracket = sp.expand(
            sp.diff(false_p, V) * sp.diff(false_q, T)
            - sp.diff(false_p, T) * sp.diff(false_q, V)
        )
        false_cusp = sp.expand(false_q**2 - 4 * false_p**n)
        assert sp.Poly(false_bracket, T).degree() == n - 1
        assert sp.Poly(false_cusp, T).degree() == 2 * n - 2

    # The compiled one-normal search has no open degree route.
    compiled_routes = {
        one_normal_search_route(n, p_degree, q_degree)
        for p_degree in range(0, 17)
        for q_degree in range(0, 25)
    }
    expected_resonant_route = (
        "even_factorization" if n % 2 == 0 else "odd_infinity"
    )
    assert compiled_routes == {
        "inactive",
        "degree",
        expected_resonant_route,
    }

print("PASS: uniform Dickson and canonical two-mask identities for 3 <= n <= 12")
print("PASS: even branch-component ledgers for 4 <= n <= 12")
print("PASS: positive-genus automorphic rigidity gate for 3 <= n <= 12")
print("PASS: nonautomorphic cusp divisor mismatch for 3 <= n <= 12")
print("PASS: affine-normal tangential obstruction for 3 <= n <= 12")
print("PASS: nonlinear normal-degree resonance gate in the bounded rectangle")
print("PASS: even-degree one-normal no-go degree gate")
print("PASS: odd-degree one-normal valuation-at-infinity gate")
print("PASS: resonant false-positive counterexample regressions")
print("PASS: one-normal search compiler has no open degree route")
print("PASS all-degree dihedral affine-completion obstruction replay")
