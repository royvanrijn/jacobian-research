#!/usr/bin/env python3
"""Exact translated full-support quartic-plane audit for all cubic symbols.

Let psi_generic be the deterministic order-four lift used by the original
ten-stratum endpoint audit.  Let psi_plus and psi_minus be the sum and
alternating sum of the fixed 24 primitive kernel-basis tensors.  This
checker audits

    phi_h + psi_generic + u*psi_plus + v*psi_minus

over Q[u,v,x,y,z] for all nine nonzero ternary-cubic orbit representatives
and the zero symbol.
"""

from __future__ import annotations

import math
import multiprocessing
import shutil
import subprocess
from functools import cache

import sympy as sp

import verify_cubic_symbol_double_saturation as cubic_audit
from verify_cubic_symbol_dense_quartic_plane_saturation import (
    dense_directions,
)


u, v = sp.symbols("u v")


@cache
def generic_offset_coordinates() -> tuple[sp.Rational, ...]:
    """Coordinates of generic_quartic_tensor in the primitive basis."""

    _, _, _, raw_basis = cubic_audit.quartic_constraint_data()
    primitive_scales = []
    for vector in raw_basis:
        denominators = [
            int(sp.denom(entry)) for entry in vector if entry != 0
        ]
        primitive_scales.append(math.lcm(*denominators))

    raw_combination = sum(
        (
            (index + 1) * vector
            for index, vector in enumerate(raw_basis)
        ),
        sp.zeros(raw_basis[0].rows, 1),
    )
    combination_denominators = [
        int(sp.denom(entry)) for entry in raw_combination if entry != 0
    ]
    common_scale = math.lcm(*combination_denominators)
    coordinates = tuple(
        sp.Rational(common_scale * (index + 1), primitive_scale)
        for index, primitive_scale in enumerate(primitive_scales)
    )
    assert all(coordinate != 0 for coordinate in coordinates)
    return coordinates


@cache
def affine_plane_tensor() -> dict[tuple[int, int, int], sp.Expr]:
    """Return psi_generic+u*psi_plus+v*psi_minus."""

    plus, minus = dense_directions()
    offset = cubic_audit.generic_quartic_tensor()
    tensor = {
        triple: sp.expand(offset[triple] + u * plus[triple] + v * minus[triple])
        for triple in offset
    }

    # Verify the offset decomposition in the same primitive basis.  Hence
    # the generic coefficient lambda_i+u+(-1)^i*v is nonzero for every i.
    basis = cubic_audit.quartic_kernel_basis_tensors()
    coordinates = generic_offset_coordinates()
    reconstructed_offset = {
        triple: sp.expand(
            sum(
                coordinate * basis[index][triple]
                for index, coordinate in enumerate(coordinates)
            )
        )
        for triple in offset
    }
    assert reconstructed_offset == offset
    return tensor


def affine_subspace_program(cubic: sp.Expr) -> str:
    """Build the exact Singular program for the translated parameter plane."""

    program = cubic_audit.singular_program(
        cubic,
        affine_plane_tensor(),
    ).replace(
        "ring coefficient_ring=0,(x,y,z),dp;",
        "ring coefficient_ring=0,(u,v,x,y,z),dp;",
    )
    anchor = (
        'print("EXT2_VECTOR_DIMENSION="'
        "+string(vdim(support_ext2)));"
    )
    diagnostics = (
        'print("EXT2_MULTIPLICITY="+string(mult(support_ext2)));'
        "ideal ext2_fitting=fitting(support_ext2,0);"
        "ideal ext2_support=std(radical(ext2_fitting));"
        "ideal parameter_space=std(ideal(x,y,z));"
        "ideal first_support_difference=simplify("
        "reduce(ext2_support,parameter_space),2);"
        "ideal second_support_difference=simplify("
        "reduce(parameter_space,ext2_support),2);"
        'print("PARAMETER_SPACE_DIFFERENCE="'
        "+string(size(first_support_difference)"
        "+size(second_support_difference)));"
        "module pruned_ext2_presentation=std(prune(support_ext2));"
        "module central_ext2_presentation=std(prune("
        "subst(subst(support_ext2,u,0),v,0)));"
        "module first_presentation_difference=simplify("
        "reduce(pruned_ext2_presentation,central_ext2_presentation),2);"
        "module second_presentation_difference=simplify("
        "reduce(central_ext2_presentation,pruned_ext2_presentation),2);"
        'print("PRUNED_PRESENTATION_DIFFERENCE="'
        "+string(size(first_presentation_difference)"
        "+size(second_presentation_difference)));"
        'print("PRUNED_PRESENTATION_RANK="'
        "+string(nrows(pruned_ext2_presentation)));"
    )
    assert anchor in program
    return program.replace(anchor, anchor + diagnostics)


def audit_stratum(name: str) -> tuple[str, tuple[int, ...]]:
    """Run the translated polynomial-plane audit for one cubic symbol."""

    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for this checker"
    result = subprocess.run(
        [singular, "-q"],
        input=affine_subspace_program(cubic_audit.CUBIC_STRATA[name]),
        text=True,
        capture_output=True,
        check=True,
        timeout=600,
    )
    wanted = (
        "SATURATION_GENERATORS",
        "EXT2_MULTIPLICITY",
        "PARAMETER_SPACE_DIFFERENCE",
        "PRUNED_PRESENTATION_DIFFERENCE",
        "PRUNED_PRESENTATION_RANK",
    )
    values: dict[str, int] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            if key in wanted:
                values[key] = int(value)
    assert set(values) == set(wanted), result.stdout + result.stderr
    return name, tuple(values[key] for key in wanted)


def main() -> None:
    # Check full support once before launching the exact Singular workers.
    coordinates = generic_offset_coordinates()
    assert len(coordinates) == 24
    assert all(coordinate != 0 for coordinate in coordinates)

    names = sorted(cubic_audit.CUBIC_STRATA)
    context = multiprocessing.get_context("fork")
    with context.Pool(processes=4) as pool:
        for name, result in pool.imap_unordered(audit_stratum, names):
            assert result == (0, 6, 0, 0, 3), (name, result)
            print(
                f"PASS: {name}: translated dense quartic plane has "
                "saturated cotangent presentation and central Ext block"
            )

    print(
        "PASS: all ten cubic-symbol strata retain a pure relative "
        "length-six support defect on the translated quartic plane"
    )
    print(
        "PASS: the generic affine-plane tensor has nonzero coordinates in "
        "all 24 fixed quartic-kernel basis directions"
    )


if __name__ == "__main__":
    main()
