#!/usr/bin/env python3
"""Verify the fixed double point at infinity on the Q-residual branch."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import sympy as sp

from research_two_pair_sic_bidegree33_t0_Q_residual import (
    LEADING_ARTIFACT,
)
from research_two_pair_sic_bidegree33_t0_stratum_leading import ROOT
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (
    substitute,
)
from verify_two_pair_sic_bidegree33_corrected_boundary import (
    t0_open_localized_export,
)


PROJECTIVE_PROBE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_t0_stratum_Q_residual_"
        "projective_probe_s1_5_ell_7_u_2_exact.json"
    )
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def fibre_homogenization(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, sp.Symbol],
    homogenizing_variable: sp.Symbol,
) -> sp.Expr:
    polynomial = sp.Poly(expression, *variables)
    degree = polynomial.total_degree()
    return sp.expand(
        sum(
            coefficient
            * variables[0] ** monomial[0]
            * variables[1] ** monomial[1]
            * homogenizing_variable ** (degree - sum(monomial))
            for monomial, coefficient in polynomial.terms()
        )
    )


def main() -> None:
    arguments = parse_arguments()
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    export = t0_open_localized_export(
        singular,
        tuple(range(2, 6)),
        0,
        arguments.timeout,
    )
    moments = dict(
        zip(range(3, 6), export["polynomials"][:-1], strict=True)
    )
    q_replacement = (("s2", "(s1^2*u-(13/3)*u)"),)
    adapted_replacement = (
        ("t1", "(s1*u-ell)"),
        ("t2", "(a*u^2)"),
    )
    s1, ell, a, u, s3, s6, s5, h, x, y, r = sp.symbols(
        "s1 ell a u s3 s6 s5 h x y r"
    )
    symbols = {
        symbol.name: symbol
        for symbol in (s1, ell, a, u, s3, s6, s5)
    }
    parsed: dict[int, sp.Expr] = {}
    leading: dict[int, sp.Expr] = {}
    local_initial: dict[int, sp.Expr] = {}
    for order in (4, 5):
        expression = substitute(
            substitute(moments[order], q_replacement),
            adapted_replacement,
        )
        parsed[order] = sp.sympify(
            expression.replace("^", "**"),
            locals=symbols,
        )
        polynomial = sp.Poly(parsed[order], s6, s5)
        fibre_degree = polynomial.total_degree()
        leading[order] = sp.factor(
            sum(
                coefficient * s6**monomial[0] * s5**monomial[1]
                for monomial, coefficient in polynomial.terms()
                if sum(monomial) == fibre_degree
            )
        )
        homogenized = fibre_homogenization(
            parsed[order],
            (s6, s5),
            h,
        )
        local = sp.expand(
            homogenized.subs(
                {s5: 1, s6: 6 * s1 * u + y, h: x},
                simultaneous=True,
            )
        )
        local_polynomial = sp.Poly(local, x, y)
        minimum_degree = min(
            sum(monomial)
            for monomial, _ in local_polynomial.terms()
        )
        assert minimum_degree == 1
        local_initial[order] = sp.expand(
            sum(
                coefficient * x**monomial[0] * y**monomial[1]
                for monomial, coefficient in local_polynomial.terms()
                if sum(monomial) == minimum_degree
            )
        )

    infinity_linear = 6 * s1 * s5 * u - s6
    second_linear = 1092 * ell * s5 + 930 * s1 * s5 * u - 155 * s6
    second_quadratic = (
        180 * s1**2 * s5**2 * u**2
        - 60 * s1 * s5 * s6 * u
        + 196 * s5**2 * u**2
        + 5 * s6**2
    )
    assert sp.cancel(
        leading[4]
        / (
            sp.Rational(1244160, 49)
            * u
            * infinity_linear
            * second_linear
        )
    ) == 1
    assert sp.cancel(
        leading[5]
        / (
            -sp.Rational(1679616000, 49)
            * infinity_linear
            * second_quadratic
        )
    ) == 1

    extra_resultant = sp.factor(
        sp.resultant(
            second_linear.subs({s6: r, s5: 1}),
            second_quadratic.subs({s6: r, s5: 1}),
            r,
        )
    )
    inherited_factor = 6084 * ell**2 + 4805 * u**2
    assert sp.expand(extra_resultant - 980 * inherited_factor) == 0, (
        extra_resultant,
        inherited_factor,
    )

    tangent_matrix = sp.Matrix(
        [
            [
                local_initial[order].coeff(x),
                local_initial[order].coeff(y),
            ]
            for order in (4, 5)
        ]
    )
    tangent_determinant = sp.factor(tangent_matrix.det())
    tangent_primitive = sp.Poly(
        sp.cancel(tangent_determinant / u**5),
        s1,
        ell,
        a,
        u,
        s3,
        domain=sp.QQ,
    ).primitive()[1]
    residual_border_factor = sp.factor(tangent_primitive.as_expr())

    leading_payload = json.loads(
        LEADING_ARTIFACT.read_text(encoding="utf-8")
    )
    border_expression = substitute(
        leading_payload["leading_coefficient_lcm"],
        adapted_replacement,
    )
    parsed_border = sp.factor(
        sp.sympify(
            border_expression.replace("^", "**"),
            locals=symbols,
        )
    )
    border_quotient = sp.cancel(
        parsed_border
        / (
            u**10
            * inherited_factor**2
            * residual_border_factor**2
        )
    )
    assert border_quotient.is_Rational and border_quotient != 0
    tangent_quotient = sp.cancel(
        tangent_determinant / (u**5 * residual_border_factor)
    )
    assert tangent_quotient.is_Rational and tangent_quotient != 0

    projective_probe = json.loads(
        PROJECTIVE_PROBE.read_text(encoding="utf-8")
    )
    projective = projective_probe["projective_probe"]
    assert projective["homogeneous_degrees"] == [2, 3]
    assert projective["projective_gcd_degree"] == 0
    assert projective["infinity_gcd_degree"] == 1
    assert projective["infinity_gcd"] == "-s6+60*s5"
    assert projective["inferred_infinity_intersection_length"] == 2

    payload = {
        "format": "two-pair-sic-bidegree33-Q-residual-infinity-v1",
        "status": (
            "exact characteristic-zero verification of the fixed "
            "double point at infinity on the dense Q-residual open; "
            "this identifies a smaller resultant construction but is "
            "not by itself a global exclusion"
        ),
        "adapted_coordinates": {
            "ell": "s1*u-t1",
            "a": "t2/u^2",
        },
        "homogeneous_degrees": {"mu4": 2, "mu5": 3},
        "common_infinity_linear": "6*s1*s5*u-s6",
        "leading_mu4": sp.sstr(leading[4]),
        "leading_mu5": sp.sstr(leading[5]),
        "extra_infinity_resultant": sp.sstr(extra_resultant),
        "inherited_Q_cap_J_factor": sp.sstr(inherited_factor),
        "residual_border_factor": sp.sstr(residual_border_factor),
        "tangent_determinant_scalar_after_u5_D": sp.sstr(
            tangent_quotient
        ),
        "border_scalar_after_u10_JQ2_D2": sp.sstr(border_quotient),
        "dense_open_conclusion": (
            "On u*J_Q != 0 and D=0, mu4^h and mu5^h have the "
            "single infinity support [s6:s5:h]=[6*s1*u:1:0]; "
            "their tangent determinant vanishes and the exact closed "
            "probe shows intersection length two."
        ),
        "closed_probe": str(PROJECTIVE_PROBE.relative_to(ROOT)),
        "reproduction_command": " ".join(sys.argv),
    }
    if arguments.output is not None:
        output = arguments.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
