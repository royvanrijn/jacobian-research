#!/usr/bin/env python3
"""Explore the power-chart certificate for degree-30 ``{2,5}``."""

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


def serialize(expression: sp.Expr) -> str:
    return str(sp.together(expression).as_numer_denom()[0]).replace(
        "**", "^"
    )


def transformed_problem() -> tuple[
    tuple[sp.Symbol, ...],
    tuple[sp.Symbol, ...],
    list[sp.Expr],
    sp.Expr,
]:
    inner_coefficients = sp.symbols("sync30p_b1:6")
    outer_coefficients = sp.symbols("sync30p_a1:5")
    source_parameters = inner_coefficients + outer_coefficients
    inner = W**6 + sum(
        inner_coefficients[power - 1] * W**power
        for power in range(1, 6)
    )
    outer = Z**5 + sum(
        outer_coefficients[power - 1] * Z**power
        for power in range(1, 5)
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

    y, z = sp.symbols("sync30p_y sync30p_z")
    translation = sp.symbols("sync30p_t")
    w0, w1 = sp.symbols("sync30p_w0 sync30p_w1")
    e1, e2 = sp.symbols("sync30p_e1 sync30p_e2")
    common_right = W**3 + e2 * W**2 + e1 * W
    power_w = z**2 + w1 * z + w0
    unshifted_outer = z * power_w**2
    shifted_inner = sp.expand(
        (y + translation) ** 2 - translation**2
    )
    shifted_outer = sp.expand(
        unshifted_outer.subs(z, z + translation**2)
        - unshifted_outer.subs(z, translation**2)
    )
    model_inner = sp.expand(shifted_inner.subs(y, common_right))
    coefficient_map = {
        inner_coefficients[power - 1]:
        sp.Poly(model_inner, W).nth(power)
        for power in range(1, 6)
    }
    coefficient_map.update(
        {
            outer_coefficients[power - 1]:
            sp.Poly(shifted_outer, z).nth(power)
            for power in range(1, 5)
        }
    )

    solved_variables = (
        inner_coefficients[0],
        inner_coefficients[1],
        outer_coefficients[0],
        outer_coefficients[1],
    )
    normal_coordinates = sp.symbols("sync30p_x1:5")
    coordinate_change = dict(coefficient_map)
    coordinate_change.update(
        {
            variable: coefficient_map[variable] + normal
            for variable, normal in zip(
                solved_variables, normal_coordinates
            )
        }
    )
    transformed_residuals = [
        sp.together(
            residual.subs(coordinate_change, simultaneous=True)
        ).as_numer_denom()[0]
        for residual in residuals
    ]
    transformed_defect = sp.together(
        defect.subs(coordinate_change, simultaneous=True)
    ).as_numer_denom()[0]
    zero_normal = {normal: 0 for normal in normal_coordinates}
    assert all(
        sp.expand(residual.subs(zero_normal)) == 0
        for residual in transformed_residuals
    )
    assert sp.expand(transformed_defect.subs(zero_normal)) == 0
    base_coordinates = (e1, e2, translation, w0, w1)
    return (
        normal_coordinates,
        base_coordinates,
        transformed_residuals,
        transformed_defect,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--order", default="dp")
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="std")
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    normals, bases, residuals, defect = transformed_problem()
    characteristic = args.prime if args.prime else 0
    variables = normals + bases
    singular = shutil.which("Singular")
    assert singular is not None
    program = (
        f'ring q={characteristic},({",".join(map(str, variables))}),'
        f"{args.order};\n"
        f'ideal I={",".join(serialize(item) for item in residuals)};\n'
        f"ideal G={args.algorithm}(I);\n"
        'print("TRANSPORTED_SYNC30_25");\n'
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
