#!/usr/bin/env python3
"""Explore monomial orders for the transported degree-30 ``{2,3}`` pair."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import (  # noqa: E402
    W,
    Z,
    canonical_linear_lift,
    canonical_residuals,
)
from jcsearch.ritt_complex import dickson_vertex_factors  # noqa: E402


def serialize(expression: sp.Expr) -> str:
    return str(sp.together(expression).as_numer_denom()[0]).replace(
        "**", "^"
    )


def transformed_problem(
    model_base: bool = False,
) -> tuple[
    tuple[sp.Symbol, ...],
    tuple[sp.Symbol, ...],
    list[sp.Expr],
    sp.Expr,
]:
    inner_coefficients = sp.symbols("sync30t_b1:10")
    outer_coefficients = sp.symbols("sync30t_a1:3")
    source_parameters = inner_coefficients + outer_coefficients
    inner = W**10 + sum(
        inner_coefficients[power - 1] * W**power
        for power in range(1, 10)
    )
    outer = (
        Z**3
        + outer_coefficients[1] * Z**2
        + outer_coefficients[0] * Z
    )
    polynomial = sp.expand(outer.subs(Z, inner))
    residuals = canonical_residuals(
        polynomial,
        2,
        15,
        parameters=source_parameters,
        factor_output=False,
        minimum_coefficient_degree=2,
    )
    defect = (
        canonical_linear_lift(polynomial, 2, 15)
        - sp.Poly(polynomial, W).nth(1)
    )

    translation, parameter = sp.symbols("sync30t_t sync30t_z")
    quintic_coefficients = sp.symbols("sync30t_e1:5")
    quintic = W**5 + sum(
        quintic_coefficients[power - 1] * W**power
        for power in range(1, 5)
    )
    collision_outer, collision_inner = dickson_vertex_factors(
        (3, 2), W, translation, parameter
    )
    transported_inner = sp.expand(collision_inner.subs(W, quintic))
    coefficient_map = {
        inner_coefficients[power - 1]:
        sp.Poly(transported_inner, W).nth(power)
        for power in range(1, 10)
    }
    coefficient_map.update(
        {
            outer_coefficients[power - 1]:
            sp.Poly(collision_outer, W).nth(power)
            for power in range(1, 3)
        }
    )

    b5, b6, b7, b8, b9 = inner_coefficients[4:]
    e4 = b9 / 2
    e3 = (b8 - e4**2) / 2
    e2 = (b7 - 2 * e4 * e3) / 2
    e1 = (b6 - 2 * e4 * e2 - e3**2) / 2
    translation_inverse = (
        b5 - 2 * e4 * e1 - 2 * e3 * e2
    ) / 2
    parameter_inverse = (
        3 * translation_inverse**2 - outer_coefficients[1]
    ) / 6
    inverse = {
        quintic_coefficients[3]: e4,
        quintic_coefficients[2]: e3,
        quintic_coefficients[1]: e2,
        quintic_coefficients[0]: e1,
        translation: translation_inverse,
        parameter: parameter_inverse,
    }

    solved_variables = (
        inner_coefficients[0],
        inner_coefficients[1],
        inner_coefficients[2],
        inner_coefficients[3],
        outer_coefficients[0],
    )
    normal_coordinates = sp.symbols("sync30t_x1:6")
    if model_base:
        normal_change = dict(coefficient_map)
        normal_change.update(
            {
                variable: coefficient_map[variable] + normal
                for variable, normal in zip(
                    solved_variables, normal_coordinates
                )
            }
        )
        base_coordinates = (
            quintic_coefficients + (translation, parameter)
        )
    else:
        graph_values = {
            variable: sp.factor(coefficient_map[variable].subs(inverse))
            for variable in solved_variables
        }
        normal_change = {
            variable: graph_values[variable] + normal
            for variable, normal in zip(
                solved_variables, normal_coordinates
            )
        }
        base_coordinates = (
            b5, b6, b7, b8, b9, outer_coefficients[1]
        )
    transformed_residuals = [
        sp.together(
            residual.subs(normal_change, simultaneous=True)
        ).as_numer_denom()[0]
        for residual in residuals
    ]
    transformed_defect = sp.together(
        defect.subs(normal_change, simultaneous=True)
    ).as_numer_denom()[0]
    zero_normal = {normal: 0 for normal in normal_coordinates}
    assert all(
        sp.expand(residual.subs(zero_normal)) == 0
        for residual in transformed_residuals
    )
    assert sp.expand(transformed_defect.subs(zero_normal)) == 0
    return (
        normal_coordinates,
        base_coordinates,
        transformed_residuals,
        transformed_defect,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--order", default="(dp(5),dp(6))")
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="std")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--model-base", action="store_true")
    args = parser.parse_args()
    normals, bases, residuals, defect = transformed_problem(args.model_base)
    characteristic = args.prime if args.prime else 0
    variables = normals + bases
    singular = shutil.which("Singular")
    assert singular is not None
    program = (
        f'ring q={characteristic},({",".join(map(str, variables))}),'
        f"{args.order};\n"
        f'ideal I={",".join(serialize(item) for item in residuals)};\n'
        f"ideal G={args.algorithm}(I);\n"
        'print("TRANSPORTED_SYNC30_23");\n'
        f"print(reduce({serialize(defect)},G)==0);\n"
        "print(size(G));\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=args.timeout,
    )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, file=sys.stderr, end="")


if __name__ == "__main__":
    main()
