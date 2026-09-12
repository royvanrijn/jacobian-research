#!/usr/bin/env python3
"""Verify the transported degree-30 ``{2,5}`` pair over ``QQ``."""

from __future__ import annotations

import argparse
import os
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
    canonical_linear_lift,
    canonical_residuals,
)
from exact_membership_cache import (  # noqa: E402
    cached_basis_size,
    program_sha256,
    write_exact_membership_cache,
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
    b3, b4, b5 = inner_coefficients[2:]
    a3, a4 = outer_coefficients[2:]
    recovered_e2 = b5 / 2
    recovered_e1 = (b4 - recovered_e2**2) / 2
    recovered_translation = (
        b3 - 2 * recovered_e2 * recovered_e1
    ) / 2
    recovered_w1 = (
        a4 - 5 * recovered_translation**2
    ) / 2
    recovered_w0 = (
        a3
        - 10 * recovered_translation**4
        - 8 * recovered_w1 * recovered_translation**2
        - recovered_w1**2
    ) / 2
    inverse = {
        e1: recovered_e1,
        e2: recovered_e2,
        translation: recovered_translation,
        w0: recovered_w0,
        w1: recovered_w1,
    }
    assert all(
        sp.expand(
            recovered.subs(coefficient_map, simultaneous=True) - variable
        ) == 0
        for variable, recovered in inverse.items()
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
    parser.add_argument("--prime", type=int, default=0)
    parser.add_argument("--order", default="dp")
    parser.add_argument(
        "--algorithm",
        choices=("std", "slimgb", "modstd"),
        default="modstd",
    )
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--fingerprint-only", action="store_true")
    args = parser.parse_args()
    normals, bases, residuals, defect = transformed_problem()
    characteristic = args.prime if args.prime else 0
    variables = normals + bases
    basis_command = (
        "modStd(I)" if args.algorithm == "modstd"
        else f"{args.algorithm}(I)"
    )
    program = (
        ('LIB "modstd.lib";\n' if args.algorithm == "modstd" else "")
        +
        f'ring q={characteristic},({",".join(map(str, variables))}),'
        f"{args.order};\n"
        f'ideal I={",".join(serialize(item) for item in residuals)};\n'
        f"ideal G={basis_command};\n"
        'print("TRANSPORTED_SYNC30_25");\n'
        f"print(reduce({serialize(defect)},G)==0);\n"
        "print(size(G));\n"
    )
    cache_path = (
        ROOT / "certificates/degree30_hessian/transported_25.json"
    )
    cache_metadata = {
        "pair": "{2,5}",
        "characteristic": 0,
        "algorithm": "modStd",
        "order": "dp",
        "coordinate_mode": "power",
    }
    cache_eligible = (
        characteristic == 0
        and args.algorithm == "modstd"
        and args.order == "dp"
    )
    if args.fingerprint_only:
        print(program_sha256(program))
        return
    cached = (
        None
        if args.refresh or not cache_eligible
        else cached_basis_size(cache_path, program, cache_metadata)
    )
    if cached is not None:
        print(
            "CACHED PASS: degree-30 pair {2,5} synchronizes in exact "
            f"transported power coordinates; basis size {cached}"
        )
        return
    singular = shutil.which("Singular")
    assert singular is not None
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
    marker = compact.index("TRANSPORTED_SYNC30_25")
    assert compact[marker + 1] == "1", stdout + stderr
    basis_size = int(compact[marker + 2])
    if cache_eligible:
        version = subprocess.run(
            [singular, "--version"],
            text=True,
            capture_output=True,
            check=True,
        ).stdout.splitlines()[0]
        write_exact_membership_cache(
            cache_path, program, cache_metadata, basis_size, version
        )
    print(
        "PASS: degree-30 pair {2,5} synchronizes in exact transported "
        f"power coordinates; basis size {basis_size}"
    )


if __name__ == "__main__":
    main()
