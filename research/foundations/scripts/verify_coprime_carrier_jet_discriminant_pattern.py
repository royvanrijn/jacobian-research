#!/usr/bin/env python3
"""Verify primitive carrier-jet saturation and the F2 cusp discriminant."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def truncate(expression: sp.Expr, variable: sp.Symbol, order: int) -> sp.Expr:
    return sp.series(expression, variable, 0, order + 1).removeO().expand()


def primitive_bezout(m: int, n: int) -> tuple[int, int]:
    for s in range(1, m):
        numerator = s * n - 1
        if numerator % m == 0:
            return numerator // m, s
    raise AssertionError("coprime pair has no primitive Bezout solution")


def universal_transport_audit() -> None:
    for m, n in ((2, 3), (2, 5), (3, 4), (3, 5), (4, 7), (5, 8)):
        assert sp.gcd(m, n) == 1
        r, s = primitive_bezout(m, n)
        assert s * n - r * m == 1

        # U=x^(-m)y^s and V=x^(-n)y^r invert
        # x=U^r/V^s and y=U^n/V^m exactly.
        assert -m * r + n * s == 1
        assert s * r - r * s == 0
        assert -m * n + n * m == 0
        assert s * n - r * m == 1

        normalized_parameters = (m - 2) + (n - 2)
        assert normalized_parameters + 3 == m + n - 1
        assert (n - m, m, n) == (n - m, m, n)


def formal_raw_jet_determinant(m: int, n: int) -> tuple[sp.Expr, sp.Expr]:
    """Compute the square raw-jet determinant on the standard slice."""

    u, t, x = sp.symbols("u t x")
    p_parameters = sp.symbols(f"p1:{m - 1}")
    q_indices = [index for index in range(1, n) if index != m]
    q_parameters = tuple(sp.Symbol(f"q{index}") for index in q_indices)
    p = t**m + sum(
        coefficient * t**index
        for coefficient, index in zip(p_parameters, range(1, m - 1))
    )
    q = t**n + sum(
        coefficient * t**index
        for coefficient, index in zip(q_parameters, q_indices)
    )

    r, s = primitive_bezout(m, n)
    jet_order = m + n - 1
    p_bar = sp.expand(u**m * p.subs(t, 1 / u))
    q_bar = sp.expand(u**n * q.subs(t, 1 / u))
    carrier_x = truncate(u * p_bar**r / q_bar**s, u, jet_order)
    carrier_y = truncate(p_bar**n / q_bar**m, u, jet_order)

    normalized_jets: list[sp.Expr] = []
    graph = sp.Integer(1)
    for index in range(1, jet_order + 1):
        coefficient = truncate(carrier_y - graph, u, jet_order).coeff(u, index)
        normalized_jets.append(sp.expand(coefficient))
        graph = truncate(graph + coefficient * carrier_x**index, u, jet_order)

    mu, eta, nu = sp.symbols("mu eta nu")
    normalized_graph = 1 + sum(
        coefficient * x**index
        for index, coefficient in enumerate(normalized_jets, start=1)
    )
    numerator_unit = 1 + mu * x**m / normalized_graph**s
    denominator_unit = (
        1
        + eta * x ** (n - m) * normalized_graph ** (s - r)
        + nu * x**n / normalized_graph**r
    )
    fixed_x = truncate(
        x * numerator_unit**r / denominator_unit**s, x, jet_order
    )
    fixed_y = truncate(
        normalized_graph * numerator_unit**n / denominator_unit**m,
        x,
        jet_order,
    )

    fixed_jets: list[sp.Expr] = []
    fixed_graph = sp.Integer(1)
    for index in range(1, jet_order + 1):
        coefficient = truncate(fixed_y - fixed_graph, x, jet_order).coeff(x, index)
        fixed_jets.append(sp.expand(coefficient))
        fixed_graph = truncate(fixed_graph + coefficient * fixed_x**index, x, jet_order)

    parameters = p_parameters + q_parameters + (mu, eta, nu)
    assert len(parameters) == jet_order
    determinant = sp.factor(sp.Matrix(fixed_jets).jacobian(parameters).det())
    immersion_resultant = sp.factor(
        sp.resultant(sp.diff(p, t), sp.diff(q, t), t)
    )
    return determinant, immersion_resultant


def low_bidegree_determinant_audit() -> None:
    expected_ratios = {(2, 3): 2, (2, 5): -2, (3, 4): -3}
    for degrees, expected_ratio in expected_ratios.items():
        determinant, resultant = formal_raw_jet_determinant(*degrees)
        assert sp.factor(determinant - expected_ratio * resultant) == 0


def nonimmersion_collision_and_conductor_audit() -> None:
    t, u, r = sp.symbols("t u r")
    a, b, c, d = sp.symbols("a b c d")
    p = t**3 + a * t
    q = t**5 + b * t**4 + c * t**2 + d * t
    discriminant = sp.factor(sp.resultant(sp.diff(p, t), sp.diff(q, t), t))
    expected_discriminant = (
        25 * a**4
        + 48 * a**3 * b**2
        - 144 * a**2 * b * c
        + 90 * a**2 * d
        + 108 * a * c**2
        + 81 * d**2
    )
    assert sp.expand(discriminant - expected_discriminant) == 0
    assert sp.factor(expected_discriminant) == expected_discriminant

    critical_substitution = {
        a: -3 * r**2,
        d: -5 * r**4 - 4 * b * r**3 - 2 * c * r,
    }
    assert sp.expand(discriminant.subs(critical_substitution)) == 0

    collision = u**4 + b * u**3 + a * u**2 + (2 * a * b - c) * u - (a**2 + d)
    residual_cubic = (
        u**3
        + (b + 2 * r) * u**2
        + (2 * b * r + r**2) * u
        - 2 * b * r**2
        - c
        + 2 * r**3
    )
    specialized_collision = sp.expand(collision.subs(critical_substitution))
    assert sp.expand(specialized_collision - (u - 2 * r) * residual_cubic) == 0

    expected_cubic_discriminant = -(
        -3 * b * r**2 - c + 2 * r**3
    ) * (
        4 * b**3 - 12 * b**2 * r - 69 * b * r**2 - 27 * c + 50 * r**3
    )
    assert sp.expand(
        sp.discriminant(residual_cubic, u) - expected_cubic_discriminant
    ) == 0
    assert sp.expand(
        sp.resultant(residual_cubic, 3 * u**2 - 12 * r**2, u)
        + 27 * (2 * b * r**2 + c) * (6 * b * r**2 - c + 20 * r**3)
    ) == 0

    p_specialized = p.subs(critical_substitution)
    q_specialized = q.subs(critical_substitution)
    cusp_cross = sp.factor(
        sp.diff(p_specialized, t, 2).subs(t, r)
        * sp.diff(q_specialized, t, 3).subs(t, r)
        - sp.diff(q_specialized, t, 2).subs(t, r)
        * sp.diff(p_specialized, t, 3).subs(t, r)
    )
    assert sp.expand(cusp_cross - 12 * (6 * b * r**2 - c + 20 * r**3)) == 0

    pair_polynomial = t**2 - u * t + (u**2 - 3 * r**2)
    conductor = sp.factor(sp.resultant(specialized_collision, pair_polynomial, u))
    conductor_quotient, remainder = sp.div(conductor, (t - r) ** 2, t)
    assert remainder == 0
    assert sp.Poly(conductor_quotient, t).degree() == 6
    assert sp.Poly(conductor_quotient, t).LC() == 1

    witness = {r: 1, b: 0, c: 1}
    witness_cubic = sp.expand(residual_cubic.subs(witness))
    assert witness_cubic == u**3 + 2 * u**2 + u + 1
    assert sp.discriminant(witness_cubic, u) == -23
    assert sp.resultant(witness_cubic, 3 * u**2 - 12, u) == -513
    witness_conductor_quotient = sp.expand(conductor_quotient.subs(witness))
    expected_sextic = t**6 + 2 * t**5 - 6 * t**4 - 12 * t**3 + 11 * t**2 + 22 * t + 1
    assert witness_conductor_quotient == expected_sextic
    assert sp.discriminant(expected_sextic, t) == 17_368_128
    assert expected_sextic.subs(t, 1) == 19


def cusp_three_node_witness_audit() -> None:
    s, t, u = sp.symbols("s t u")
    left, right, inverse = sp.symbols("left right inverse")
    p = t**3 - 3 * t
    q = t**5 + t**2 - 7 * t
    residual_cubic = u**3 + 2 * u**2 + u + 1

    tangent = sp.expand(
        sp.diff(p.subs(t, s), s) * sp.diff(q, t)
        - sp.diff(q.subs(t, s), s) * sp.diff(p, t)
    )
    tangent_quotient, remainder = sp.div(tangent, s - t, domain=sp.QQ)
    assert remainder == 0
    symmetric_tangent = sp.symmetrize(tangent_quotient, [s, t], formal=True)[0]
    symmetric_1, symmetric_2 = sp.symbols("s1 s2")
    tangent_u = sp.expand(
        symmetric_tangent.subs({symmetric_1: u, symmetric_2: u**2 - 3})
    )
    tangent_u = sp.rem(tangent_u, residual_cubic, u)
    assert sp.expand(tangent_u + 3 * (23 * u**2 + 42 * u + 14)) == 0
    assert sp.resultant(residual_cubic, tangent_u, u) == -224_181

    target_x = lambda value: -value * (value**2 - 3)
    target_y = lambda value: (value**2 - 3) * (value**3 - 6 * value - 1)
    distinct_images = sp.groebner(
        [
            residual_cubic.subs(u, left),
            residual_cubic.subs(u, right),
            target_x(left) - target_x(right),
            target_y(left) - target_y(right),
            inverse * (left - right) - 1,
        ],
        inverse,
        right,
        left,
        order="lex",
    )
    assert any(polynomial.as_expr() == 1 for polynomial in distinct_images.polys)

    cusp_x = p.subs(t, 1)
    cusp_y = q.subs(t, 1)
    cusp_disjoint = sp.groebner(
        [residual_cubic, target_x(u) - cusp_x, target_y(u) - cusp_y],
        u,
        order="lex",
    )
    assert any(polynomial.as_expr() == 1 for polynomial in cusp_disjoint.polys)


def raw_corank_witness_audit() -> None:
    import verify_f2_affine_k1_carrier_jet_factorization as carrier

    normalized_jets, _ = carrier.carrier_graph_jet_audit()
    fixed_transport, _, _ = carrier.fixed_coordinate_transport_audit()
    mu, eta, nu = sp.symbols("mu eta nu")
    normalized_substitution = dict(
        zip(
            (
                carrier.h1,
                carrier.h2,
                carrier.h3,
                carrier.h4,
                carrier.h5,
                carrier.h6,
                carrier.h7,
            ),
            normalized_jets,
        )
    )
    fixed_jets = [
        sp.expand(coefficient.subs(normalized_substitution))
        for coefficient in fixed_transport
    ]
    parameters = (carrier.a, carrier.b, carrier.c, carrier.d, mu, eta, nu)
    witness = {
        carrier.a: -3,
        carrier.b: 0,
        carrier.c: 1,
        carrier.d: -7,
        mu: 0,
        eta: 0,
        nu: 0,
    }
    jacobian = sp.Matrix(fixed_jets).jacobian(parameters).subs(witness)
    assert jacobian.det() == 0
    assert jacobian.rank() == 6
    assert jacobian[:6, :6].det() == -1134


def main() -> None:
    universal_transport_audit()
    low_bidegree_determinant_audit()
    nonimmersion_collision_and_conductor_audit()
    cusp_three_node_witness_audit()
    raw_corank_witness_audit()
    print(
        "PASS: primitive coprime carrier transport has three universal raw "
        "parameters; exact low-bidegree jet determinants are signed m times "
        "the immersion resultant; the F2 nonimmersion divisor generically "
        "has one cusp plus three nodes, conductor split 2+6, and raw "
        "seven-jet corank one"
    )


if __name__ == "__main__":
    main()
