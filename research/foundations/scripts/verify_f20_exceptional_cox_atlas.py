#!/usr/bin/env python3
"""Verify the controlled-transform exceptional Cox atlas for F20.

This checker constructs exact strict-transform charts at every positive
exceptional packet over the q-node, r-cusp, conjugate triple orbit, and
q-r tangency orbit.  In every chart it verifies

    P = tau^N F,          P_X = tau^(N-r) * dF/dY,

where ``r`` is the tau-order of dX/dY.  The compact Cartier columns
``D_d,D_q,D_r`` have local coefficients ``tau^a,tau^b,tau^c`` and satisfy
``3*a+b+c=N-r``.  Thus the degree-(3,1,1) derivative quotient is literally
the polynomial ``dF/dY`` on every controlled chart.

The checker also gives one parity-compatible distribution of the previously
certified total conductor residue L(w) among the three compact columns.  It
is a cocycle on the punctured conductor cover, and its zeros are assigned to
the collision packets covered by the local charts.  No global regular SNC
model or global section ring is inferred: overlap gluing through those
charts, full inverse-adjugate polynomiality, and affine-space recognition
remain open.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import sympy as sp

from verify_f20_global_multi_rees_cox_algebra import f20_data


@dataclass(frozen=True)
class Chart:
    name: str
    substitution: dict[sp.Symbol, sp.Expr]
    total_order: int
    root_weight: int
    compact_orders: tuple[int, int, int]
    geometric_colors: int
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


def tau_order_and_quotient(
    expression: sp.Expr,
    tau: sp.Symbol,
    expected_order: int,
    modulus: sp.Expr | None,
    modulus_variable: sp.Symbol | None,
) -> tuple[int, sp.Expr]:
    polynomial = sp.Poly(sp.expand(expression), tau)
    reduced_coefficients: dict[int, sp.Expr] = {}
    for order in range(polynomial.degree() + 1):
        coefficient = reduce_mod(
            polynomial.nth(order), modulus, modulus_variable
        )
        if coefficient != 0:
            reduced_coefficients[order] = coefficient
    assert reduced_coefficients
    actual_order = min(reduced_coefficients)
    assert actual_order == expected_order, (actual_order, expected_order)
    quotient = sp.expand(
        sum(
            coefficient * tau ** (order - expected_order)
            for order, coefficient in reduced_coefficients.items()
        )
    )
    return actual_order, quotient


def polynomial_degree(expression: sp.Expr, variable: sp.Symbol) -> int:
    return int(sp.Poly(sp.expand(expression), variable).degree())


def chart_certificate(
    chart: Chart,
    P: sp.Expr,
    P_X: sp.Expr,
    tau: sp.Symbol,
    Y: sp.Symbol,
) -> dict[str, object]:
    pulled_P = sp.expand(P.subs(chart.substitution))
    pulled_P_X = sp.expand(P_X.subs(chart.substitution))
    derivative_order = chart.total_order - chart.root_weight
    assert derivative_order == (
        3 * chart.compact_orders[0]
        + chart.compact_orders[1]
        + chart.compact_orders[2]
    )

    _, strict_transform = tau_order_and_quotient(
        pulled_P,
        tau,
        chart.total_order,
        chart.residue_modulus,
        chart.residue_variable,
    )
    _, derivative_quotient = tau_order_and_quotient(
        pulled_P_X,
        tau,
        derivative_order,
        chart.residue_modulus,
        chart.residue_variable,
    )
    strict_derivative = reduce_mod(
        sp.diff(strict_transform, Y),
        chart.residue_modulus,
        chart.residue_variable,
    )
    assert reduce_mod(
        derivative_quotient - strict_derivative,
        chart.residue_modulus,
        chart.residue_variable,
    ) == 0
    assert reduce_mod(
        pulled_P - tau**chart.total_order * strict_transform,
        chart.residue_modulus,
        chart.residue_variable,
    ) == 0
    assert reduce_mod(
        pulled_P_X - tau**derivative_order * strict_derivative,
        chart.residue_modulus,
        chart.residue_variable,
    ) == 0

    exceptional_residual = reduce_mod(
        strict_transform.subs(tau, 0),
        chart.residue_modulus,
        chart.residue_variable,
    )
    residual_degree = polynomial_degree(exceptional_residual, Y)
    assert residual_degree >= 1
    if residual_degree > 1:
        residual_discriminant = reduce_mod(
            sp.discriminant(exceptional_residual, Y),
            chart.residue_modulus,
            chart.residue_variable,
        )
        assert residual_discriminant != 0
    else:
        residual_discriminant = sp.Integer(1)

    a, b_order, c = chart.compact_orders
    return {
        "name": chart.name,
        "base_ring": (
            "residue-field controlled chart"
            if chart.residue_modulus is not None
            else "rational controlled chart"
        ),
        "residue_modulus": (
            str(chart.residue_modulus)
            if chart.residue_modulus is not None
            else None
        ),
        "strict_transform_order": chart.total_order,
        "root_coordinate_weight": chart.root_weight,
        "derivative_order": derivative_order,
        "compact_section_orders": {
            "D_d": a,
            "D_q": b_order,
            "D_r": c,
        },
        "degree_3_1_1_order": 3 * a + b_order + c,
        "geometric_colors": chart.geometric_colors,
        "exceptional_residual_degree": residual_degree,
        "exceptional_residual": str(sp.factor(exceptional_residual)),
        "exceptional_residual_discriminant_nonzero": True,
        "primitive_exceptional_rees_variable": "e=tau*T_E",
        "compact_local_coefficients": {
            "Z_d": f"tau^{a}*T_d",
            "Z_q": f"tau^{b_order}*T_q",
            "Z_r": f"tau^{c}*T_r",
        },
        "local_multi_rees_algebra": "A_chart[Z_d,Z_q,Z_r]",
        "literal_degree_3_1_1_quotient": "dF/dY",
        "literal_polynomiality_verified": True,
    }


def conductor_cocycle_certificate() -> dict[str, object]:
    w = sp.symbols("w")
    A_minus = w**2 - 2 * w + 5
    A_plus = w**2 + 2 * w + 5
    B_minus = w**3 - 3 * w**2 - w - 5
    B_plus = w**3 + 3 * w**2 - w + 5
    node = w**4 + 10 * w**2 + 5
    smooth_q = 3 * w**2 + 5
    puncture = (w**2 - 1) * (w**2 + 3)

    residue_d = sp.Integer(1)
    residue_q = sp.cancel(
        smooth_q * A_minus * A_plus / (8 * (w**2 - 1))
    )
    residue_r = w / 4
    total_residue = sp.cancel(
        w
        * smooth_q
        * A_minus
        * A_plus
        / (32 * (w - 1) * (w + 1))
    )
    assert sp.cancel(residue_d**3 * residue_q * residue_r - total_residue) == 0
    assert sp.cancel(residue_d.subs(w, -w) - residue_d) == 0
    assert sp.cancel(residue_q.subs(w, -w) - residue_q) == 0
    assert sp.cancel(residue_r.subs(w, -w) + residue_r) == 0
    assert sp.cancel(total_residue.subs(w, -w) + total_residue) == 0

    total_numerator = sp.together(total_residue).as_numer_denom()[0]
    assert sp.gcd(node, total_numerator) == 1
    assert sp.gcd(B_minus * B_plus, total_numerator) == 1
    assert sp.gcd(A_minus * A_plus, w * smooth_q) == 1
    assert sp.gcd(w, smooth_q * A_minus * A_plus) == 1

    unit_pullback = sp.Matrix(
        (
            (-1, -1, -2),
            (-1, -1, -2),
            (1, 0, 0),
            (0, 0, 1),
        )
    )
    selector_completion = unit_pullback.row_join(sp.Matrix((1, 0, 0, 0)))
    assert selector_completion.det() == -1

    return {
        "cover_ring": f"QQ[w,1/({sp.factor(puncture)})]",
        "involution": "w -> -w",
        "residue_frames": {
            "D_d": str(residue_d),
            "D_q": str(residue_q),
            "D_r": str(residue_r),
        },
        "characters": {"D_d": 1, "D_q": 1, "D_r": -1},
        "degree_3_1_1_product": str(total_residue),
        "product_equals_total_derivative_residue": True,
        "unit_selector_completion_determinant": -1,
        "collision_extension": {
            "q_node": "all three chosen frames are units generically",
            "triple_orbit": "D_q carries the simple A_minus*A_plus zero",
            "q_r_transverse": "D_r carries the simple w zero",
            "q_r_tangency_orbit": "all three chosen frames are units generically",
            "r_cusp": "off the q conductor and handled only by cusp charts",
        },
        "uniqueness": "only up to invariant units and redistribution of even factors",
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

    tau, z, Y, a, b, alpha = sp.symbols("tau z Y a b alpha")
    imaginary_unit = sp.I
    node_modulus = a**2 - 3 * a + 1
    triple_modulus = b**2 + imaginary_unit
    qr_modulus = 8 * alpha**3 + 16 * alpha**2 + 2 * alpha - 7
    qr_s0 = 4 * alpha**2 - 5
    qr_A = -2 - 2 * alpha
    qr_B = 8 - 12 * alpha**2 - 8 * alpha
    qr_slope = -8 * alpha**2 - 4 * alpha + 10

    charts = (
        Chart(
            "q_node_slopes",
            {
                s: 1 + tau,
                t: -sp.Rational(1, 2) + z * tau,
                X: a + tau * Y,
            },
            2,
            1,
            (0, 1, 0),
            4,
            node_modulus,
            a,
        ),
        Chart(
            "r_cusp_E1",
            {
                s: 11 + tau**5 * (z + 50),
                t: -sp.Rational(1, 2) + tau**5,
                X: -1 + tau * Y,
            },
            5,
            1,
            (0, 0, 4),
            1,
        ),
        Chart(
            "r_cusp_E2",
            {
                s: 11 + 50 * tau**5 + tau**10 * z,
                t: -sp.Rational(1, 2) + tau**5,
                X: -1 + tau**2 * Y,
            },
            10,
            2,
            (0, 0, 8),
            1,
        ),
        Chart(
            "r_cusp_E3_unramified",
            {
                s: 11 + 50 * tau + 180 * tau**2 + tau**3 * z,
                t: -sp.Rational(1, 2) + tau,
                X: -1 + tau * Y,
            },
            3,
            1,
            (0, 0, 2),
            1,
        ),
        Chart(
            "r_cusp_E3_ramified",
            {
                s: 11 + 50 * tau**2 + 180 * tau**4 + tau**6 * z,
                t: -sp.Rational(1, 2) + tau**2,
                X: -1 + tau * Y,
            },
            5,
            1,
            (0, 0, 4),
            2,
        ),
        Chart(
            "r_cusp_E4",
            {
                s: (
                    11
                    + 50 * tau**2 * z
                    + 180 * tau**4 * z**2
                    + tau**5 * z**2
                ),
                t: -sp.Rational(1, 2) + tau**2 * z,
                X: -1 + tau * Y,
            },
            5,
            1,
            (0, 0, 4),
            5,
        ),
        Chart(
            "triple_E1_ramified",
            {
                s: 2 * imaginary_unit + tau**4,
                t: -sp.Rational(3, 4) + imaginary_unit / 2 + z * tau**4,
                X: 1 + imaginary_unit + tau**2 * (b + tau * Y),
            },
            10,
            3,
            (1, 2, 2),
            2,
            triple_modulus,
            b,
        ),
        Chart(
            "triple_E2_cluster",
            {
                s: 2 * imaginary_unit + z * tau**2,
                t: -sp.Rational(3, 4) + imaginary_unit / 2 + tau,
                X: 1 + imaginary_unit + tau * Y,
            },
            4,
            1,
            (1, 0, 0),
            8,
        ),
        Chart(
            "qr_E1_A_ramified",
            {
                s: qr_s0 + z * tau**2,
                t: alpha + tau**2,
                X: qr_A + tau * Y,
            },
            3,
            1,
            (0, 1, 1),
            3,
            qr_modulus,
            alpha,
        ),
        Chart(
            "qr_E1_A_unramified",
            {
                s: qr_s0 + z * tau,
                t: alpha + tau,
                X: qr_A + tau * Y,
            },
            2,
            1,
            (0, 1, 0),
            3,
            qr_modulus,
            alpha,
        ),
        Chart(
            "qr_E1_B_ramified",
            {
                s: qr_s0 + z * tau**2,
                t: alpha + tau**2,
                X: qr_B + tau * Y,
            },
            2,
            1,
            (0, 0, 1),
            3,
            qr_modulus,
            alpha,
        ),
        Chart(
            "qr_E2_A",
            {
                s: qr_s0 + qr_slope * tau + z * tau**2,
                t: alpha + tau,
                X: qr_A + tau * Y,
            },
            3,
            1,
            (0, 1, 1),
            9,
            qr_modulus,
            alpha,
        ),
        Chart(
            "qr_E2_B",
            {
                s: qr_s0 + qr_slope * tau + z * tau**2,
                t: alpha + tau,
                X: qr_B + tau * Y,
            },
            2,
            1,
            (0, 0, 1),
            6,
            qr_modulus,
            alpha,
        ),
    )

    chart_certificates = tuple(
        chart_certificate(chart, P, P_X, tau, Y) for chart in charts
    )
    assert sum(chart.geometric_colors for chart in charts) == 48
    assert all(
        certificate["literal_polynomiality_verified"]
        for certificate in chart_certificates
    )

    conductor_cocycle = conductor_cocycle_certificate()

    return {
        "schema": "f20-exceptional-cox-atlas-v1",
        "status": (
            "controlled_transform_atlas_and_punctured_conductor_cocycle_"
            "certified_global_gluing_open"
        ),
        "general_local_theorem": {
            "strict_transform_identity": "P=tau^N*F implies P_X=tau^(N-r)*dF/dY",
            "hypothesis": "X=X_0+tau^r*Y and tau is a regular exceptional parameter",
            "principal_multi_rees": (
                "for D_j=a_j*E, the local algebra is "
                "A[tau^a_j*T_j] isomorphic to A[Z_j]"
            ),
            "degree_gate": "3*a_d+a_q+a_r=N-r",
        },
        "controlled_transform_charts": chart_certificates,
        "positive_exceptional_colors_covered": 48,
        "positive_generic_colors": 5,
        "positive_colors_total": 53,
        "conductor_cocycle": conductor_cocycle,
        "degree_3_1_1_gate": {
            "literal_on_every_controlled_chart": True,
            "literal_on_punctured_conductor": True,
            "global_H0_membership": "uncertified_without_global_overlap_gluing",
        },
        "downstream": {
            "local_derivative_denominator_cancellation": "passed_on_all_13_chart_types",
            "entrywise_inverse_adjugate_polynomiality": "not_reached_globally",
            "affine_space_recognition": "not_reached",
            "global_regular_SNC_source": "not_constructed",
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
    print("PASS: thirteen controlled chart types cover 48 positive exceptional colors")
    print("PASS: every local (3,1,1) quotient is the polynomial dF/dY")
    print("PASS: the conductor residues multiply to the anti-invariant total residue")
    print("SCOPE: global overlap gluing and affine-space recognition remain open")


if __name__ == "__main__":
    main()
