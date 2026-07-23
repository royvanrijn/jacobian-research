#!/usr/bin/env python3
"""Explore the transported degree-42 ``{2,7}`` synchronization pair.

The primitive degree-14 power collision is composed on the right with a
generic cubic.  Leading source-chart coefficients recover the six model
parameters, leaving five refinement-normal coordinates.
"""

from __future__ import annotations

import argparse
import itertools
import os
import re
import signal
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
    canonical_reconstruction,
)


def serialize(expression: sp.Expr) -> str:
    serialized = str(expression).replace("**", "^")
    # Singular treats ``(poly/number)^n`` as ``poly ^ number``.  Write every
    # numeric denominator as multiplication by a field constant while
    # preserving the compact arithmetic circuit.
    return re.sub(r"/([0-9]+)", r"*(1/\1)", serialized)


def solve_triangular(
    expression: sp.Expr,
    variable: sp.Symbol,
) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), variable)
    assert polynomial.degree() == 1
    return sp.cancel(-polynomial.nth(0) / polynomial.nth(1))


def raw_canonical_problem(
    polynomial: sp.Expr,
    outer_degree: int,
    inner_degree: int,
) -> tuple[list[sp.Expr], sp.Expr]:
    """Return canonical residual circuits without primitive expansion."""

    source, composition, reconstruction, used_degrees = (
        canonical_reconstruction(
            polynomial,
            outer_degree,
            inner_degree,
        )
    )
    residuals = [
        composition.nth(coefficient_degree).subs(
            reconstruction,
            simultaneous=True,
        )
        - source.nth(coefficient_degree)
        for coefficient_degree in range(2, outer_degree * inner_degree)
        if coefficient_degree not in used_degrees
    ]
    defect = (
        composition.nth(1).subs(reconstruction, simultaneous=True)
        - source.nth(1)
    )
    return residuals, defect


def transformed_problem() -> tuple[
    tuple[sp.Symbol, ...],
    tuple[sp.Symbol, ...],
    list[sp.Expr],
    sp.Expr,
]:
    inner_coefficients = sp.symbols("sync42p_b1:6")
    outer_coefficients = sp.symbols("sync42p_a1:7")
    source_parameters = inner_coefficients + outer_coefficients
    inner = W**6 + sum(
        inner_coefficients[power - 1] * W**power
        for power in range(1, 6)
    )
    outer = Z**7 + sum(
        outer_coefficients[power - 1] * Z**power
        for power in range(1, 7)
    )
    polynomial = sp.expand(outer.subs(Z, inner))
    residuals, defect = raw_canonical_problem(
        polynomial,
        2,
        21,
    )

    y, z = sp.symbols("sync42p_y sync42p_z")
    translation = sp.symbols("sync42p_t")
    w0, w1, w2 = sp.symbols("sync42p_w0 sync42p_w1 sync42p_w2")
    e1, e2 = sp.symbols("sync42p_e1 sync42p_e2")
    common_right = W**3 + e2 * W**2 + e1 * W
    power_w = z**3 + w2 * z**2 + w1 * z + w0
    unshifted_outer = z * power_w**2
    shifted_inner = sp.expand(
        (y + translation) ** 2 - translation**2
    )
    shifted_outer = sp.expand(
        unshifted_outer.subs(z, z + translation**2)
        - unshifted_outer.subs(z, translation**2)
    )
    core_source = shifted_outer.subs(z, shifted_inner)
    inner_constant = translation * power_w.subs(
        z,
        translation**2,
    )
    target_inner = (
        (y + translation)
        * power_w.subs(z, (y + translation) ** 2)
        - inner_constant
    )
    target_outer = Z**2 + 2 * inner_constant * Z
    assert sp.expand(
        core_source - target_outer.subs(Z, target_inner)
    ) == 0
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
            for power in range(1, 7)
        }
    )

    b3, b4, b5 = inner_coefficients[2:]
    a4, a5, a6 = outer_coefficients[3:]
    recovered_e2 = b5 / 2
    recovered_e1 = (b4 - recovered_e2**2) / 2
    recovered_translation = (
        b3 - 2 * recovered_e2 * recovered_e1
    ) / 2
    recovered_w2 = solve_triangular(
        coefficient_map[a6].subs(
            translation,
            recovered_translation,
        )
        - a6,
        w2,
    )
    recovered_w1 = solve_triangular(
        coefficient_map[a5].subs(
            {
                translation: recovered_translation,
                w2: recovered_w2,
            },
            simultaneous=True,
        )
        - a5,
        w1,
    )
    recovered_w0 = solve_triangular(
        coefficient_map[a4].subs(
            {
                translation: recovered_translation,
                w2: recovered_w2,
                w1: recovered_w1,
            },
            simultaneous=True,
        )
        - a4,
        w0,
    )
    inverse = {
        e1: recovered_e1,
        e2: recovered_e2,
        translation: recovered_translation,
        w0: recovered_w0,
        w1: recovered_w1,
        w2: recovered_w2,
    }
    assert all(
        sp.expand(
            recovered.subs(coefficient_map, simultaneous=True) - variable
        ) == 0
        for variable, recovered in inverse.items()
    )

    base_source_variables = (b3, b4, b5, a4, a5, a6)
    solved_variables = tuple(
        variable
        for variable in source_parameters
        if variable not in base_source_variables
    )
    assert len(solved_variables) == 5
    normal_coordinates = sp.symbols("sync42p_x1:6")
    coordinate_change = dict(coefficient_map)
    coordinate_change.update(
        {
            variable: coefficient_map[variable] + normal
            for variable, normal in zip(
                solved_variables,
                normal_coordinates,
            )
        }
    )
    # Keep the substituted polynomials as arithmetic circuits.  Calling
    # together/expand here causes degree-42 expression swell before Singular
    # sees the much smaller five-normal-coordinate problem.
    transformed_residuals = [
        residual.subs(coordinate_change, simultaneous=True)
        for residual in residuals
    ]
    transformed_defect = defect.subs(
        coordinate_change,
        simultaneous=True,
    )
    base_coordinates = (
        e1,
        e2,
        translation,
        w0,
        w1,
        w2,
    )
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
    parser.add_argument(
        "--algorithm",
        choices=("std", "slimgb", "modstd"),
        default="slimgb",
    )
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument(
        "--normal-order",
        type=int,
        default=0,
        help=(
            "work modulo the (normal-order + 1)-st power of the normal "
            "maximal ideal; zero requests the full ideal"
        ),
    )
    args = parser.parse_args()
    assert args.normal_order >= 0
    normals, bases, residuals, defect = transformed_problem()
    print(
        "PREPARED: degree-42 pair {2,7}; "
        f"{len(normals)} normal, {len(bases)} base, "
        f"{len(residuals)} residual coordinates"
    )
    if args.prepare_only:
        return

    characteristic = args.prime if args.prime else 0
    variables = normals + bases
    singular = shutil.which("Singular")
    assert singular is not None
    basis_command = (
        "modStd(I)" if args.algorithm == "modstd"
        else f"{args.algorithm}(I)"
    )
    truncation_generators = []
    if args.normal_order:
        truncation_generators = [
            sp.prod(monomial)
            for monomial in itertools.combinations_with_replacement(
                normals,
                args.normal_order + 1,
            )
        ]
    ideal_generators = residuals + truncation_generators
    program = (
        ('LIB "modstd.lib";\n' if args.algorithm == "modstd" else "")
        + f'ring q={characteristic},({",".join(map(str, variables))}),'
        f"{args.order};\n"
        f'ideal I={",".join(serialize(item) for item in ideal_generators)};\n'
        f"ideal G={basis_command};\n"
        'print("TRANSPORTED_SYNC42_27");\n'
        f"print(reduce({serialize(defect)},G)==0);\n"
        "print(size(G));\n"
    )
    process = subprocess.Popen(
        [singular, "-q"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(program, timeout=args.timeout)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        raise
    compact = stdout.split()
    marker = compact.index("TRANSPORTED_SYNC42_27")
    assert compact[marker + 1] == "1", stdout + stderr
    field = f"GF({characteristic})" if characteristic else "QQ"
    neighborhood = (
        f", normal order {args.normal_order}"
        if args.normal_order
        else ""
    )
    print(
        "PASS: degree-42 pair {2,7} synchronizes in transported "
        f"power coordinates over {field}{neighborhood}; basis size "
        f"{compact[marker + 2]}"
    )


if __name__ == "__main__":
    main()
