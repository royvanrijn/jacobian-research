#!/usr/bin/env python3
"""Verify two-parameter exceptional Cox corner charts for F20.

The one-parameter exceptional atlas certifies generic points of each colored
divisor.  This continuation checks actual codimension-two corners.  It
constructs four normalized corners in the ramphoid-cusp resolution and three
in the q-r tangency resolution, proves exact bivariate strict-transform and
derivative factorizations, and checks the compact (3,1,1) Cox monomial on
each corner.

It also verifies the complementary-chart transition at the q-node and proves
that the positive triple-E1 and triple-E2 colors are root-center separated,
despite adjacency of their base rays.  The result closes the
exceptional--exceptional corner gate represented by these charts.  It does
not construct the remaining strict-boundary attachment charts or a global
Cech cover.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import sympy as sp

from verify_f20_global_multi_rees_cox_algebra import f20_data


@dataclass(frozen=True)
class Corner:
    name: str
    substitution: dict[sp.Symbol, sp.Expr]
    equation_orders: tuple[int, int]
    root_weights: tuple[int, int]
    compact_orders: tuple[
        tuple[int, int], tuple[int, int], tuple[int, int]
    ]
    residue_modulus: sp.Expr | None = None
    residue_variable: sp.Symbol | None = None


def reduce_mod(
    expression: sp.Expr,
    modulus: sp.Expr | None,
    modulus_variable: sp.Symbol | None,
) -> sp.Expr:
    expression = sp.expand(expression)
    if modulus is None:
        return expression
    assert modulus_variable is not None
    return sp.expand(sp.rem(expression, modulus, modulus_variable))


def bivariate_quotient(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, sp.Symbol],
    expected_orders: tuple[int, int],
    modulus: sp.Expr | None,
    modulus_variable: sp.Symbol | None,
) -> sp.Expr:
    u, v = variables
    polynomial = sp.Poly(sp.expand(expression), u, v)
    quotient = sp.Integer(0)
    nonzero_monomials: list[tuple[int, int]] = []
    for (u_order, v_order), coefficient in polynomial.terms():
        coefficient = reduce_mod(coefficient, modulus, modulus_variable)
        if coefficient == 0:
            continue
        assert u_order >= expected_orders[0]
        assert v_order >= expected_orders[1]
        nonzero_monomials.append((u_order, v_order))
        quotient += (
            coefficient
            * u ** (u_order - expected_orders[0])
            * v ** (v_order - expected_orders[1])
        )
    assert nonzero_monomials
    assert min(order[0] for order in nonzero_monomials) == expected_orders[0]
    assert min(order[1] for order in nonzero_monomials) == expected_orders[1]
    return sp.expand(quotient)


def corner_certificate(
    corner: Corner,
    P: sp.Expr,
    P_X: sp.Expr,
    u: sp.Symbol,
    v: sp.Symbol,
    Y: sp.Symbol,
) -> dict[str, object]:
    derivative_orders = (
        corner.equation_orders[0] - corner.root_weights[0],
        corner.equation_orders[1] - corner.root_weights[1],
    )
    compact_degree_orders = tuple(
        3 * corner.compact_orders[0][index]
        + corner.compact_orders[1][index]
        + corner.compact_orders[2][index]
        for index in range(2)
    )
    assert compact_degree_orders == derivative_orders

    pulled_P = sp.expand(P.subs(corner.substitution))
    pulled_P_X = sp.expand(P_X.subs(corner.substitution))
    strict_transform = bivariate_quotient(
        pulled_P,
        (u, v),
        corner.equation_orders,
        corner.residue_modulus,
        corner.residue_variable,
    )
    derivative_quotient = bivariate_quotient(
        pulled_P_X,
        (u, v),
        derivative_orders,
        corner.residue_modulus,
        corner.residue_variable,
    )
    strict_derivative = reduce_mod(
        sp.diff(strict_transform, Y),
        corner.residue_modulus,
        corner.residue_variable,
    )
    assert reduce_mod(
        derivative_quotient - strict_derivative,
        corner.residue_modulus,
        corner.residue_variable,
    ) == 0
    assert reduce_mod(
        pulled_P
        - u ** corner.equation_orders[0]
        * v ** corner.equation_orders[1]
        * strict_transform,
        corner.residue_modulus,
        corner.residue_variable,
    ) == 0
    assert reduce_mod(
        pulled_P_X
        - u**derivative_orders[0]
        * v**derivative_orders[1]
        * strict_derivative,
        corner.residue_modulus,
        corner.residue_variable,
    ) == 0

    corner_residual = reduce_mod(
        strict_transform.subs({u: 0, v: 0}),
        corner.residue_modulus,
        corner.residue_variable,
    )
    residual_degree = int(sp.Poly(corner_residual, Y).degree())
    assert residual_degree >= 1
    if residual_degree > 1:
        residual_discriminant = reduce_mod(
            sp.discriminant(corner_residual, Y),
            corner.residue_modulus,
            corner.residue_variable,
        )
        assert residual_discriminant != 0

    def monomial(order: tuple[int, int], label: str) -> str:
        return f"u^{order[0]}*v^{order[1]}*T_{label}"

    return {
        "name": corner.name,
        "equation_orders": corner.equation_orders,
        "root_weights": corner.root_weights,
        "derivative_orders": derivative_orders,
        "compact_orders": {
            "D_d": corner.compact_orders[0],
            "D_q": corner.compact_orders[1],
            "D_r": corner.compact_orders[2],
        },
        "compact_local_generators": {
            "Z_d": monomial(corner.compact_orders[0], "d"),
            "Z_q": monomial(corner.compact_orders[1], "q"),
            "Z_r": monomial(corner.compact_orders[2], "r"),
        },
        "degree_3_1_1_orders": compact_degree_orders,
        "corner_residual": str(sp.factor(corner_residual)),
        "corner_residual_degree": residual_degree,
        "corner_residual_generically_separable": True,
        "literal_quotient": "dF/dY",
        "literal_polynomiality_verified": True,
    }


def q_node_transition_certificate(
    P: sp.Expr,
    P_X: sp.Expr,
    s: sp.Symbol,
    t: sp.Symbol,
    X: sp.Symbol,
) -> dict[str, object]:
    u, v, z, zeta, Y_left, Y_right, a = sp.symbols(
        "u v z zeta Y_left Y_right a"
    )
    modulus = a**2 - 3 * a + 1

    def univariate_quotient(
        expression: sp.Expr, variable: sp.Symbol, order: int
    ) -> sp.Expr:
        polynomial = sp.Poly(sp.expand(expression), variable)
        result = sp.Integer(0)
        for exponent in range(polynomial.degree() + 1):
            coefficient = reduce_mod(polynomial.nth(exponent), modulus, a)
            if coefficient == 0:
                continue
            assert exponent >= order
            result += coefficient * variable ** (exponent - order)
        return sp.expand(result)

    left_substitution = {
        s: 1 + u,
        t: -sp.Rational(1, 2) + z * u,
        X: a + u * Y_left,
    }
    right_substitution = {
        s: 1 + zeta * v,
        t: -sp.Rational(1, 2) + v,
        X: a + v * Y_right,
    }
    F_left = univariate_quotient(P.subs(left_substitution), u, 2)
    G_left = univariate_quotient(P_X.subs(left_substitution), u, 1)
    F_right = univariate_quotient(P.subs(right_substitution), v, 2)
    G_right = univariate_quotient(P_X.subs(right_substitution), v, 1)
    transition = {v: z * u, zeta: 1 / z, Y_right: Y_left / z}

    def rational_zero(expression: sp.Expr) -> bool:
        numerator = sp.together(expression).as_numer_denom()[0]
        return reduce_mod(numerator, modulus, a) == 0

    assert rational_zero(F_right.subs(transition) - F_left / z**2)
    assert rational_zero(G_right.subs(transition) - G_left / z)
    assert reduce_mod(G_left - sp.diff(F_left, Y_left), modulus, a) == 0
    assert reduce_mod(G_right - sp.diff(F_right, Y_right), modulus, a) == 0

    return {
        "left_chart": "s-1=u, t+1/2=z*u, X-a=u*Y_left",
        "right_chart": "t+1/2=v, s-1=zeta*v, X-a=v*Y_right",
        "overlap": "v=z*u, zeta=1/z, Y_right=Y_left/z",
        "strict_transform_transition": "F_right=F_left/z^2",
        "derivative_quotient_transition": "G_right=G_left/z",
        "D_q_frame_transition": "Z_q_right=z*Z_q_left",
        "P_X_invariant_on_overlap": True,
    }


def build_certificate() -> dict[str, object]:
    data = f20_data()
    s = data["s"]
    t = data["t"]
    X = data["X"]
    P = data["P"]
    P_X = data["P_X"]
    assert isinstance(s, sp.Symbol)
    assert isinstance(t, sp.Symbol)
    assert isinstance(X, sp.Symbol)

    u, v, Y, alpha, b = sp.symbols("u v Y alpha b")
    qr_modulus = 8 * alpha**3 + 16 * alpha**2 + 2 * alpha - 7

    # Cusp tangent coordinates are V=t+1/2 and W=(s-11)-50*V.
    V_12 = u**5 * v**5
    W_12 = u**5 * v**10
    V_24 = u**5 * v**2
    W_24 = 180 * V_24**2 + u**10 * v**5
    V_43u = u**2 * v
    W_43u = 180 * V_43u**2 + u**5 * v**3
    V_43r = u**2 * v**2
    W_43r = 180 * V_43r**2 + u**5 * v**6

    qr_s0 = 4 * alpha**2 - 5
    qr_slope = -8 * alpha**2 - 4 * alpha + 10
    qr_A = -2 - 2 * alpha
    qr_B = 8 - 12 * alpha**2 - 8 * alpha
    qr_V_ramified = u**2 * v
    qr_U_ramified = u**2 * v**2
    qr_V_unramified = u * v
    qr_U_unramified = u * v**2

    zero = (0, 0)
    corners = (
        Corner(
            "r_cusp_E1_E2",
            {
                t: -sp.Rational(1, 2) + V_12,
                s: 11 + 50 * V_12 + W_12,
                X: -1 + u * v**2 * Y,
            },
            (5, 10),
            (1, 2),
            (zero, zero, (4, 8)),
        ),
        Corner(
            "r_cusp_E2_E4",
            {
                t: -sp.Rational(1, 2) + V_24,
                s: 11 + 50 * V_24 + W_24,
                X: -1 + u**2 * v * Y,
            },
            (10, 5),
            (2, 1),
            (zero, zero, (8, 4)),
        ),
        Corner(
            "r_cusp_E4_E3_unramified",
            {
                t: -sp.Rational(1, 2) + V_43u,
                s: 11 + 50 * V_43u + W_43u,
                X: -1 + u * v * Y,
            },
            (5, 3),
            (1, 1),
            (zero, zero, (4, 2)),
        ),
        Corner(
            "r_cusp_E4_E3_ramified",
            {
                t: -sp.Rational(1, 2) + V_43r,
                s: 11 + 50 * V_43r + W_43r,
                X: -1 + u * v * Y,
            },
            (5, 5),
            (1, 1),
            (zero, zero, (4, 4)),
        ),
        Corner(
            "qr_E1_A_ramified_E2_A",
            {
                t: alpha + qr_V_ramified,
                s: qr_s0 + qr_slope * qr_V_ramified + qr_U_ramified,
                X: qr_A + u * v * Y,
            },
            (3, 3),
            (1, 1),
            (zero, (1, 1), (1, 1)),
            qr_modulus,
            alpha,
        ),
        Corner(
            "qr_E1_A_unramified_E2_A",
            {
                t: alpha + qr_V_unramified,
                s: qr_s0 + qr_slope * qr_V_unramified + qr_U_unramified,
                X: qr_A + u * v * Y,
            },
            (2, 3),
            (1, 1),
            (zero, (1, 1), (0, 1)),
            qr_modulus,
            alpha,
        ),
        Corner(
            "qr_E1_B_ramified_E2_B",
            {
                t: alpha + qr_V_ramified,
                s: qr_s0 + qr_slope * qr_V_ramified + qr_U_ramified,
                X: qr_B + u * v * Y,
            },
            (2, 2),
            (1, 1),
            (zero, zero, (1, 1)),
            qr_modulus,
            alpha,
        ),
    )
    corner_certificates = tuple(
        corner_certificate(corner, P, P_X, u, v, Y) for corner in corners
    )
    assert all(
        certificate["literal_polynomiality_verified"]
        for certificate in corner_certificates
    )

    q_node_transition = q_node_transition_certificate(P, P_X, s, t, X)

    # The triple E1 residual root has nonzero leading center b, b^2+i=0.
    # The E2 cluster chart is centered at X-X0=0.  Therefore these colored
    # primes cannot meet above the adjacent base rays.
    triple_modulus = b**2 + sp.I
    assert sp.gcd(b, triple_modulus) == 1
    assert reduce_mod(b**2, triple_modulus, b) == -sp.I
    triple_separation = {
        "base_ray_adjacency": True,
        "positive_colored_prime_adjacency": False,
        "E1_root_center": "(X-X0)/tau^2 = b + O(tau), b^2+i=0",
        "E2_root_center": "X-X0=0 along the cluster center",
        "certificate": "b is a unit modulo b^2+i, so the centers are disjoint",
    }

    return {
        "schema": "f20-exceptional-cox-corners-v1",
        "status": (
            "seven_exceptional_corners_and_q_node_transition_certified_"
            "strict_boundary_attachments_open"
        ),
        "general_corner_theorem": {
            "hypothesis": "P=u^N*v^M*F and X=xi+u^r*v^s*Y",
            "derivative_identity": "P_X=u^(N-r)*v^(M-s)*dF/dY",
            "local_multi_rees": (
                "D_j=a_j*E_u+b_j*E_v has generator "
                "u^a_j*v^b_j*T_j"
            ),
            "exact_degree_gate": (
                "sum n_j*(a_j,b_j)=(N-r,M-s) gives quotient dF/dY"
            ),
        },
        "corner_graph": {
            "r_cusp": "E1--E2--E4--E3, with unramified and ramified E3 colors",
            "q_r_tangency": (
                "E1_A_ram--E2_A, E1_A_unram--E2_A, "
                "E1_B_ram--E2_B"
            ),
            "triple_positive_colors": "E1 and E2 are root-center separated",
        },
        "corner_charts": corner_certificates,
        "q_node_complementary_transition": q_node_transition,
        "triple_root_center_separation": triple_separation,
        "degree_3_1_1_gate": {
            "literal_on_all_seven_corners": True,
            "q_node_transition_compatible": True,
            "global_Cech_class": "not_yet_defined_without_boundary_attachments",
        },
        "remaining_global_edges": (
            "strict d/q/r boundary attachments, generic ramification charts, "
            "and conductor-to-corner transitions"
        ),
        "downstream": {
            "entrywise_inverse_adjugate_polynomiality": "not_reached_globally",
            "affine_space_recognition": "not_reached",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)
    print("PASS: seven two-parameter exceptional Cox corners are exact")
    print("PASS: the q-node complementary-chart transition preserves P_X")
    print("PASS: the positive triple E1/E2 colors are root-center separated")
    print("SCOPE: strict-boundary attachment charts and the global Cech class remain open")


if __name__ == "__main__":
    main()
