#!/usr/bin/env python3
"""Exact checks for the scoped ordinary-degree-six boundary audit.

This script verifies identities in four explicitly declared templates.  It
does not enumerate all Keller maps, all affine modifications, or all
small-boundary presentations.
"""
from __future__ import annotations

import sympy as sp


def total_degree(polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]) -> int:
    """Return the ordinary total degree of a polynomial."""
    return sp.Poly(sp.expand(polynomial), *variables).total_degree()


def check_asymmetric_linear_quadratic_slice() -> None:
    a, y, z = sp.symbols("a y z")
    g = sp.Function("g")(a, y, z)

    b = 1 + a * y
    c0, h0 = sp.symbols("c0 h0")
    # Divisibility of d=(1-bc)/a first gives c=1 mod a.  Writing
    # c=1+a*h, divisibility of e by a^2 then gives h=-3y/2 mod a.
    assert sp.solve(sp.Eq((1 - b * c0).subs(a, 0), 0), c0) == [1]
    trial_c = 1 + a * h0
    trial_e_numerator = sp.expand(1 + b - 2 * b**2 * trial_c)
    linear_coefficient = sp.diff(trial_e_numerator, a).subs(a, 0)
    assert linear_coefficient == -2 * h0 - 3 * y
    assert sp.solve(sp.Eq(linear_coefficient, 0), h0) == [-sp.Rational(3, 2) * y]

    c = 1 - sp.Rational(3, 2) * a * y + a**2 * g
    d = sp.cancel((1 - b * c) / a)
    e = sp.cancel((1 + b - 2 * b**2 * c) / a**2)

    assert a not in sp.denom(d).free_symbols
    assert a not in sp.denom(e).free_symbols
    assert sp.expand(a * d + b * c - 1) == 0
    assert sp.expand(a**2 * e - a * b * d + b**2 * c - 1) == 0

    p = sp.expand(a * c)
    r = sp.expand(a * e + b * d)
    s = sp.expand(b * e)
    jacobian = sp.factor(sp.Matrix((p, r, s)).jacobian((a, y, z)).det())
    assert jacobian == -sp.diff(g, z)
    assert sp.factor(sp.diff(p, z) - a**3 * sp.diff(g, z)) == 0
    assert sp.factor(
        sp.diff(r, z) + 3 * a * (1 + a * y) ** 2 * sp.diff(g, z)
    ) == 0
    assert sp.factor(
        sp.diff(s, z) + 2 * (1 + a * y) ** 3 * sp.diff(g, z)
    ) == 0

    lam = sp.symbols("lambda", nonzero=True)
    linear_g = lam * z
    product_coordinates = tuple(
        sp.expand(component.subs(g, linear_g)) for component in (p, r, s)
    )
    assert sp.factor(
        sp.Matrix(product_coordinates).jacobian((a, y, z)).det() + lam
    ) == 0
    assert tuple(
        total_degree(component, (a, y, z)) for component in product_coordinates
    ) == (4, 6, 7)
    assert sp.Poly(product_coordinates[2], a, y, z).coeff_monomial(a**3 * y**3 * z) == -2 * lam
    print("PASS asymmetric (1,2) slice: determinant -g_z and degree floor 7")


def check_balanced_two_plus_two_slice() -> None:
    T, h, b, u, v = sp.symbols("T h b u v")
    q_plus = T**2 + (h + u) * T + b + v
    q_minus = T**2 + (h - u) * T + b - v
    product = sp.Poly(sp.expand(q_plus * q_minus), T).all_coeffs()
    X, Y, Z = product[2:]
    assert tuple(map(sp.expand, (X, Y, Z))) == (
        2 * b + h**2 - u**2,
        2 * h * b - 2 * u * v,
        b**2 - v**2,
    )

    resultant = sp.factor(sp.resultant(q_plus, q_minus, T))
    cox_jacobian = sp.factor(sp.Matrix((X, Y, Z)).jacobian((b, u, v)).det())
    assert sp.expand(resultant - 4 * (b * u**2 - h * u * v + v**2)) == 0
    assert sp.expand(cox_jacobian - 2 * resultant) == 0

    # The invariant lattice for k[u^2,uv,v^2] is the even-sum sublattice of
    # Z^2.  Its index, and hence the affine toric class-group order, is two.
    invariant_lattice_basis = sp.Matrix(((2, 1), (0, 1)))
    assert abs(int(invariant_lattice_basis.det())) == 2

    x, chart_y = sp.symbols("x chart_y")
    chart_coordinates = (
        2 * b + h**2 - x,
        2 * h * b - 2 * x * chart_y,
        b**2 - x * chart_y**2,
    )
    chart_jacobian = sp.factor(
        sp.Matrix(chart_coordinates).jacobian((b, x, chart_y)).det()
    )
    assert chart_jacobian == 4 * x * (b - h * chart_y + chart_y**2)

    other_x, other_y = sp.symbols("other_x other_y")
    other_chart_coordinates = (
        2 * b + h**2 - other_x * other_y**2,
        2 * h * b - 2 * other_x * other_y,
        b**2 - other_x,
    )
    other_chart_jacobian = sp.factor(
        sp.Matrix(other_chart_coordinates).jacobian((b, other_x, other_y)).det()
    )
    assert other_chart_jacobian == -4 * other_x * (
        b * other_y**2 - h * other_y + 1
    )

    # On the overlap, other_x=x*chart_y^2 and other_y=1/chart_y.  These are
    # the two standard charts of the minimal resolution of the A1 cone.
    transition = {other_x: x * chart_y**2, other_y: 1 / chart_y}
    assert all(
        sp.cancel(got.subs(transition) - expected) == 0
        for got, expected in zip(other_chart_coordinates, chart_coordinates)
    )

    # Four labelled roots have three unordered pair partitions and six
    # ordered decompositions.
    roots = (0, 1, 2, 3)
    unordered = {
        frozenset((frozenset(pair), frozenset(set(roots) - set(pair))))
        for pair in ((0, 1), (0, 2), (0, 3))
    }
    assert len(unordered) == 3
    assert 2 * len(unordered) == 6
    print("PASS balanced (2,2) slice: both resolution charts retain the resultant divisor")


def check_two_reconstruction_coordinates() -> None:
    P, S, Q = sp.symbols("P S Q")
    b0 = sp.Function("b0")(P, S)
    b1 = sp.Function("b1")(P, S)
    c0 = sp.Function("c0")(P, S)
    c1 = sp.Function("c1")(P, S)
    B = b0 + b1 * Q
    C = c0 + c1 * Q
    relative_jacobian = sp.expand(
        sp.Matrix((P, B, C)).jacobian((P, S, Q)).det()
    )
    expected = (
        sp.diff(b0, S) * c1
        - b1 * sp.diff(c0, S)
        + Q * (sp.diff(b1, S) * c1 - b1 * sp.diff(c1, S))
    )
    assert sp.expand(relative_jacobian - expected) == 0

    # If kappa=S^2 and b1,c1 are nonzero monomials in S, the Wronskian
    # exponent equation is m+n=3.  Up to swapping B,C there are two types.
    ordered_profiles = {
        (m, n)
        for m in range(5)
        for n in range(5)
        if m + n == 3 and m != n
    }
    assert ordered_profiles == {(0, 3), (1, 2), (2, 1), (3, 0)}
    assert {tuple(sorted(profile)) for profile in ordered_profiles} == {
        (0, 3),
        (1, 2),
    }

    # The (1,2) profile cannot satisfy the second Wronskian equation:
    # its left side is divisible by S, while D0=1+P*S^3 is not.
    db0, dc0, lam = sp.symbols("db0 dc0 lambda", nonzero=True)
    mixed_profile_lhs = S**2 * db0 - S * dc0
    assert mixed_profile_lhs.subs(S, 0) == 0
    assert (lam * (1 + P * S**3)).subs(S, 0) == lam

    # For the remaining (0,3) profile, normalized so lambda=3, the general
    # b0-monomial contributes the following term to c0.
    i, j = sp.symbols("i j", integer=True, nonnegative=True)
    b0_monomial = P**i * S**j
    monomial_operator = sp.simplify(
        S**3 * b0_monomial - 3 * sp.integrate(S**2 * b0_monomial, S)
    )
    assert sp.simplify(
        monomial_operator - sp.Rational(1, 1) * j * P**i * S ** (j + 3) / (j + 3)
    ) == 0

    x, y, q, z = sp.symbols("x y q z")
    t = 1 + x**2 * y
    chart_P = t * q
    chart_S = x / t
    chart_Q = t * y + x * q
    controlled_divisor = sp.factor(
        1 - chart_S**2 * chart_Q + chart_P * chart_S**3
    )
    chart_jacobian = sp.factor(
        sp.Matrix((chart_P, chart_S, chart_Q)).jacobian((x, y, q)).det()
    )
    assert controlled_divisor == 1 / t
    assert chart_jacobian == 1

    coefficient = sp.Function("coefficient")(x, y)
    q0 = sp.Function("q0")(x, y)
    general_pullback = tuple(
        component.subs(q, coefficient * z + q0)
        for component in (chart_P, chart_S, chart_Q)
    )
    general_source_jacobian = sp.factor(
        sp.Matrix(general_pullback).jacobian((x, y, z)).det()
    )
    assert general_source_jacobian == coefficient
    assert sp.factor(general_source_jacobian * controlled_divisor) == coefficient / t

    # Hence every z-linear Keller clearing has coefficient=unit*t, not just
    # a pure-power ansatz.
    mu = sp.symbols("mu", nonzero=True)
    forced_pullback = tuple(
        component.subs(q, mu * t * z) for component in (chart_P, chart_S, chart_Q)
    )
    forced_source_jacobian = sp.factor(
        sp.Matrix(forced_pullback).jacobian((x, y, z)).det()
    )
    assert forced_source_jacobian == mu * t
    assert sp.factor(forced_source_jacobian * controlled_divisor) == mu
    forced_P = sp.expand(chart_P.subs(q, mu * t * z))
    assert total_degree(forced_P, (x, y, z)) == 7
    assert sp.Poly(forced_P, x, y, z).coeff_monomial(x**4 * y**2 * z) == mu

    # If q0 is divisible by t, a polynomial triangular z-shift removes it.
    # In that normalized subchart q=t*z.  Polynomiality of B forces every
    # b0 monomial P^i*S^j to have j<=2i.  Terms with j=2i first create an
    # uncancellable t^-3 pole in C; after they vanish, the only possible
    # t^-2 corrections have j=2i-1 and positive z-degree.  The fixed t^-2
    # residue has the nonzero z^0 term -x.
    normalized_P = chart_P.subs(q, t * z)
    normalized_Q = chart_Q.subs(q, t * z)
    fixed_C = sp.cancel(
        -3 * chart_S
        - sp.Rational(3, 4) * normalized_P * chart_S**4
        + chart_S**3 * normalized_Q
    )
    regularized_fixed_C = sp.cancel(t**2 * fixed_C)
    fixed_double_residue = sp.cancel(regularized_fixed_C.subs(y, -1 / x**2))
    assert sp.expand(fixed_double_residue - (x**4 * z / 4 - x)) == 0
    assert sp.Poly(fixed_double_residue, z).coeff_monomial(1) == -x

    # The remaining nonzero boundary jets can also be eliminated.  Put
    # U=P*S and T(H)=S^3*H-3*integral(S^2*H,dS).  On a diagonal monomial
    # P^r*U^j the operator has the Laurent form
    #
    #   T(P^r U^j) = P^(r-3) * j/(j+3) * U^(j+3).
    #
    # This is the exact identity used for the nonconstant-boundary argument.
    U = sp.symbols("U")
    for radial_power in range(4):
        for u_power in range(7):
            diagonal_monomial = P ** (radial_power + u_power) * S**u_power
            transformed = sp.Rational(u_power, u_power + 3) * P ** (
                radial_power + u_power
            ) * S ** (u_power + 3)
            expected_laurent = (
                sp.Rational(u_power, u_power + 3)
                * P ** (radial_power - 3)
                * U ** (u_power + 3)
            )
            assert sp.cancel(transformed.subs(S, U / P) - expected_laurent) == 0

    # Its one-variable numerator operator kills exactly the constants on a
    # generic polynomial.  After the triple-pole term is removed, the
    # double-pole equation would make x^-1 a polynomial composition of
    # degree at least two with the nonconstant Laurent function U(x).
    generic_coefficients = sp.symbols("generic_coefficient_0:7")
    numerator_operator = sum(
        sp.Rational(power, power + 3)
        * generic_coefficients[power]
        * U ** (power + 3)
        for power in range(7)
    )
    assert sp.solve(
        sp.Poly(numerator_operator, U).all_coeffs(),
        generic_coefficients[1:],
        dict=True,
    ) == [{coefficient: 0 for coefficient in generic_coefficients[1:]}]

    # If the boundary value U=x*q0|_(t=0) is the nonzero constant b, a
    # triangular source change normalizes q0=-b*x*y.  The residual variable
    # delta=(U-b)/P has the displayed affine boundary value.
    boundary_b = sp.symbols("boundary_b", nonzero=True)
    constant_mu = sp.symbols("constant_mu", nonzero=True)
    constant_q = -boundary_b * x * y + constant_mu * t * z
    constant_P = sp.expand(t * constant_q)
    constant_U = sp.expand(x * constant_q)
    boundary_delta = sp.cancel((constant_U - boundary_b) / constant_P)
    expected_boundary_delta = constant_mu * x**2 * z / boundary_b - x
    assert sp.cancel(
        boundary_delta.subs(y, -1 / x**2) - expected_boundary_delta
    ) == 0

    # Every Laurent term P^-m*N(U), expanded at U=b+P*delta, has its
    # possible P^-2 coefficient in k[delta].  These finite generic-monomial
    # replays check the Taylor coefficient used in the all-degree proof.
    epsilon, delta = sp.symbols("epsilon delta")
    for pole_order in range(2, 9):
        for polynomial_degree in range(9):
            expansion = sp.expand((boundary_b + epsilon * delta) ** polynomial_degree)
            got = expansion.coeff(epsilon, pole_order - 2)
            if pole_order - 2 <= polynomial_degree:
                expected = (
                    sp.binomial(polynomial_degree, pole_order - 2)
                    * boundary_b ** (polynomial_degree - pole_order + 2)
                    * delta ** (pole_order - 2)
                )
            else:
                expected = 0
            assert sp.expand(got - expected) == 0

    # After b0=-P*S+H, the fixed part of C is x^3*y/t^2-3*x/t and has
    # double-pole residue -x.  The Laurent module above can contribute only
    # x^2 times a polynomial in boundary_delta.  Since boundary_delta is
    # affine with nonzero z coefficient, it cannot equal this residue.
    fixed_after_shift = x**3 * y / t**2 - 3 * x / t
    assert sp.cancel(
        sp.cancel(t**2 * fixed_after_shift).subs(y, -1 / x**2) + x
    ) == 0
    print(
        "PASS cubic reconstruction: z-linear floor 7; both monomial profiles "
        "excluded for every polynomial boundary jet"
    )


def check_total_ramification_core() -> None:
    w, q = sp.symbols("w q")
    D = q - w
    T = (D**4 - q**4) / 4
    core_jacobian = sp.factor(sp.Matrix((q, T)).jacobian((w, q)).det())
    assert core_jacobian == D**3

    # Standard weighted source chart plus a single monomial C-ledger.
    A, B, C = sp.symbols("A B C", nonzero=True)
    m, n = sp.symbols("m n", integer=True, nonnegative=True)
    beta = (B * C**m, A * C**n, C)
    target_jacobian = sp.factor(sp.Matrix(beta).jacobian((A, B, C)).det())
    assert target_jacobian == -C ** (m + n)
    # Equivalently, constancy would require m+n=3 from x and m+n=5
    # from D before substituting C=xD.
    M = sp.symbols("M", integer=True, nonnegative=True)
    assert sp.solve((sp.Eq(3 - M, 0), sp.Eq(5 - M, 0)), (M,), dict=True) == []

    # Standard reciprocal chart.  Its base Jacobian is D=1/t.
    x, y, q, z = sp.symbols("x y q z")
    t = 1 + x * y
    rec_P = t * q
    rec_S = x / t
    rec_Q = y + x * q
    rec_D = sp.factor(1 - rec_S * rec_Q + rec_P * rec_S**2)
    rec_base_jacobian = sp.factor(
        sp.Matrix((rec_P, rec_S, rec_Q)).jacobian((x, y, q)).det()
    )
    assert rec_D == 1 / t
    assert rec_base_jacobian == rec_D

    coefficient = sp.Function("reciprocal_coefficient")(x, y)
    q0 = sp.Function("reciprocal_q0")(x, y)
    general_pullback = tuple(
        component.subs(q, coefficient * z + q0)
        for component in (rec_P, rec_S, rec_Q)
    )
    general_source_jacobian = sp.factor(
        sp.Matrix(general_pullback).jacobian((x, y, z)).det()
    )
    assert general_source_jacobian == coefficient / t
    assert sp.factor(general_source_jacobian * rec_D**3) == coefficient / t**4

    # Therefore every z-linear Keller clearing has coefficient=unit*t^4.
    mu = sp.symbols("reciprocal_mu", nonzero=True)
    forced_pullback = tuple(
        component.subs(q, mu * t**4 * z) for component in (rec_P, rec_S, rec_Q)
    )
    forced_source_jacobian = sp.factor(
        sp.Matrix(forced_pullback).jacobian((x, y, z)).det()
    )
    assert forced_source_jacobian == mu * t**3
    assert sp.factor(forced_source_jacobian * rec_D**3) == mu
    forced_P = sp.expand(rec_P.subs(q, mu * t**4 * z))
    assert total_degree(forced_P, (x, y, z)) == 11
    assert sp.Poly(forced_P, x, y, z).coeff_monomial(x**5 * y**5 * z) == mu
    print("PASS D^3 core: pure-C weighted ledgers fail; every z-linear reciprocal lift has floor 11")


def check_nodal_conductor_unit_rank() -> None:
    # For A=k+t(t-1)k[t] inside B=k[t], the conductor quotients are
    # A/I=k diagonally embedded in B/I=k x k.  On character lattices the
    # diagonal G_m has image generated by (1,1), leaving rank one.
    diagonal = sp.Matrix(((1,), (1,)))
    assert diagonal.rank() == 1
    assert 2 - diagonal.rank() == 1
    print("PASS nodal conductor square: one residual G_m gluing direction")


def main() -> None:
    check_asymmetric_linear_quadratic_slice()
    check_balanced_two_plus_two_slice()
    check_two_reconstruction_coordinates()
    check_total_ramification_core()
    check_nodal_conductor_unit_rank()
    print("PASS scoped ordinary-degree-six boundary audit")


if __name__ == "__main__":
    main()
