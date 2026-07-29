#!/usr/bin/env python3
"""Research finite-jet module extensions in the degree-42 Ritt flag.

This script extracts standard-monomial bases and multiplication matrices for

    I_6 < I_boundary < K

after adding the q-th power of the completed maximal ideal.  It then treats
the three kernels as modules over ``Q[tau,zeta]`` and tests whether

    0 -> I_boundary/I_6 -> K/I_6 -> K/I_boundary -> 0

splits on the selected finite jet.

The calculation is exploratory until its output is promoted to a pinned
artifact and verifier assertions.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import canonical_residuals  # noqa: E402
from explore_degree42_ritt_spectator_universality import (  # noqa: E402
    ALL_CUTS,
    DEGREE,
    W,
    build_chart,
    serialize_ideal,
)
from verify_degree42_ritt_relative_cotangent_cone import (  # noqa: E402
    dickson_normal_map,
    serialize_polynomial,
)


def singular_data(order: int) -> tuple[str, tuple[sp.Symbol, ...]]:
    """Return Singular output for standard bases, maps, and base actions."""

    parameters, factor_variables, polynomial = build_chart()
    base_cuts = {2, 14}
    requested_cuts = tuple(cut for cut in ALL_CUTS if cut not in base_cuts)
    residuals = {
        cut: canonical_residuals(
            polynomial,
            cut,
            DEGREE // cut,
            parameters=parameters,
            factor_output=False,
            minimum_coefficient_degree=1,
        )
        for cut in requested_cuts
    }
    endpoint = residuals[3] + residuals[21]
    thick = endpoint + residuals[7]
    boundary = endpoint + residuals[6] + residuals[7]

    normals, base_coordinates, images = dickson_normal_map(factor_variables)
    local_variables = normals + base_coordinates
    map_images = ",".join(
        serialize_polynomial(images[parameter]) for parameter in parameters
    )
    maximal = ",".join(map(str, local_variables))

    def print_reductions(
        source_basis: str,
        target_basis: str,
        marker: str,
        multiplier: str = "",
    ) -> str:
        return f"""
print("{marker}");
module sourceBasis={source_basis};
for (int i=1; i<=ncols(sourceBasis); i++)
{{
  print(reduce(({multiplier})*sourceBasis[i],{target_basis}));
}}
"""

    program = f"""
ring source=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={serialize_ideal(thick)};
ideal IBsource={serialize_ideal(boundary)};
ring q=0,({",".join(map(str, local_variables))}),(dp({len(normals)}),dp(2));
map phi=source,{map_images};
option(redSB);
ideal IT=phi(ITsource);
ideal IB=phi(IBsource);
ideal K={",".join(map(str, normals))};
ideal maximalIdeal={maximal};
ideal maximalPower=maximalIdeal^{order};
ideal G6=std(IT+maximalPower);
ideal GB=std(IB+maximalPower);
ideal GK=std(K+maximalPower);
module B6=kbase(G6);
module BB=kbase(GB);
module BK=kbase(GK);
print("BASIS_6");
print(B6);
print("BASIS_BOUNDARY");
print(BB);
print("BASIS_K");
print(BK);
{print_reductions("B6", "GB", "MAP_6_BOUNDARY", "1")}
{print_reductions("B6", "GK", "MAP_6_K", "1")}
{print_reductions("BB", "GK", "MAP_BOUNDARY_K", "1")}
{print_reductions("B6", "G6", "TAU_6", str(base_coordinates[0]))}
{print_reductions("B6", "G6", "ZETA_6", str(base_coordinates[1]))}
{print_reductions("BB", "GB", "TAU_BOUNDARY", str(base_coordinates[0]))}
{print_reductions("BB", "GB", "ZETA_BOUNDARY", str(base_coordinates[1]))}
"""
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=7200,
    )
    return result.stdout, local_variables


def parse_sections(output: str) -> dict[str, list[str]]:
    """Parse marker-delimited Singular polynomials."""

    markers = (
        "BASIS_6",
        "BASIS_BOUNDARY",
        "BASIS_K",
        "MAP_6_BOUNDARY",
        "MAP_6_K",
        "MAP_BOUNDARY_K",
        "TAU_6",
        "ZETA_6",
        "TAU_BOUNDARY",
        "ZETA_BOUNDARY",
    )
    sections: dict[str, list[str]] = {marker: [] for marker in markers}
    active: str | None = None
    for raw_line in output.splitlines():
        line = raw_line.strip()
        if line in sections:
            active = line
            continue
        if active is None or not line:
            continue
        sections[active].append(line.rstrip(","))
    return sections


def sympy_expression(text: str, symbols: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Parse Singular's polynomial display into SymPy."""

    names = {str(symbol): symbol for symbol in symbols}
    normalized = re.sub(r"\^", "**", text)
    return sp.expand(sp.sympify(normalized, locals=names))


def coordinate_matrix(
    expressions: list[str],
    basis_text: list[str],
    symbols: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    """Write displayed normal forms as columns in a standard-monomial basis."""

    basis = tuple(sympy_expression(text, symbols) for text in basis_text)
    basis_monomials = {
        sp.Poly(monomial, *symbols).monoms()[0]: index
        for index, monomial in enumerate(basis)
    }
    matrix = sp.zeros(len(basis), len(expressions))
    for column, text in enumerate(expressions):
        polynomial = sp.Poly(sympy_expression(text, symbols), *symbols)
        for monomial, coefficient in polynomial.terms():
            matrix[basis_monomials[monomial], column] = coefficient
    return matrix


def basis_matrix(nullspace: list[sp.Matrix], ambient: int) -> sp.Matrix:
    """Stack a nullspace basis as columns, including the zero case."""

    if not nullspace:
        return sp.zeros(ambient, 0)
    return sp.Matrix.hstack(*nullspace)


def restricted_action(
    action: sp.Matrix,
    subspace: sp.Matrix,
) -> sp.Matrix:
    """Restrict an ambient action to an invariant column subspace."""

    if subspace.cols == 0:
        return sp.zeros(0, 0)
    columns = []
    for column in range(subspace.cols):
        solution = subspace.gauss_jordan_solve(
            action * subspace[:, column]
        )[0]
        columns.append(solution)
    return sp.Matrix.hstack(*columns)


def induced_map(
    ambient_map: sp.Matrix,
    source_subspace: sp.Matrix,
    target_subspace: sp.Matrix,
) -> sp.Matrix:
    """Express an ambient map between invariant subspaces."""

    if source_subspace.cols == 0:
        return sp.zeros(target_subspace.cols, 0)
    columns = []
    for column in range(source_subspace.cols):
        image = ambient_map * source_subspace[:, column]
        solution = target_subspace.gauss_jordan_solve(image)[0]
        columns.append(solution)
    return sp.Matrix.hstack(*columns)


def splitting_solution(
    projection: sp.Matrix,
    total_actions: tuple[sp.Matrix, ...],
    quotient_actions: tuple[sp.Matrix, ...],
) -> tuple[bool, sp.Matrix | None]:
    """Solve for a base-linear section of a finite-dimensional extension."""

    total_dimension = projection.cols
    quotient_dimension = projection.rows
    variables = sp.symbols(
        f"s0:{total_dimension * quotient_dimension}"
    )
    section = sp.Matrix(
        total_dimension,
        quotient_dimension,
        variables,
    )
    equations = list(projection * section - sp.eye(quotient_dimension))
    for total_action, quotient_action in zip(
        total_actions, quotient_actions
    ):
        equations.extend(total_action * section - section * quotient_action)
    coefficient_matrix, right_hand_side = sp.linear_eq_to_matrix(
        equations, variables
    )
    try:
        solution, parameters = coefficient_matrix.gauss_jordan_solve(
            right_hand_side
        )
    except ValueError:
        return False, None
    zero_parameters = {
        parameter: sp.Integer(0) for parameter in parameters
    }
    values = solution.subs(zero_parameters)
    return True, sp.Matrix(
        total_dimension,
        quotient_dimension,
        list(values),
    )


def audit(order: int) -> dict[str, object]:
    """Compute the finite-jet modules and their extension."""

    output, symbols = singular_data(order)
    sections = parse_sections(output)
    basis_6 = sections["BASIS_6"]
    basis_boundary = sections["BASIS_BOUNDARY"]
    basis_k = sections["BASIS_K"]
    map_6_boundary = coordinate_matrix(
        sections["MAP_6_BOUNDARY"], basis_boundary, symbols
    )
    map_6_k = coordinate_matrix(sections["MAP_6_K"], basis_k, symbols)
    map_boundary_k = coordinate_matrix(
        sections["MAP_BOUNDARY_K"], basis_k, symbols
    )
    tau_6 = coordinate_matrix(sections["TAU_6"], basis_6, symbols)
    zeta_6 = coordinate_matrix(sections["ZETA_6"], basis_6, symbols)
    tau_boundary = coordinate_matrix(
        sections["TAU_BOUNDARY"], basis_boundary, symbols
    )
    zeta_boundary = coordinate_matrix(
        sections["ZETA_BOUNDARY"], basis_boundary, symbols
    )

    sector_basis = basis_matrix(
        map_6_boundary.nullspace(), len(basis_6)
    )
    total_basis = basis_matrix(map_6_k.nullspace(), len(basis_6))
    spectator_basis = basis_matrix(
        map_boundary_k.nullspace(), len(basis_boundary)
    )
    projection = induced_map(
        map_6_boundary, total_basis, spectator_basis
    )

    sector_actions = (
        restricted_action(tau_6, sector_basis),
        restricted_action(zeta_6, sector_basis),
    )
    total_actions = (
        restricted_action(tau_6, total_basis),
        restricted_action(zeta_6, total_basis),
    )
    spectator_actions = (
        restricted_action(tau_boundary, spectator_basis),
        restricted_action(zeta_boundary, spectator_basis),
    )
    splits, section = splitting_solution(
        projection, total_actions, spectator_actions
    )

    assert map_boundary_k * map_6_boundary == map_6_k
    assert projection.rank() == spectator_basis.cols
    assert sector_basis.cols + spectator_basis.cols == total_basis.cols
    assert all(
        action * sector_basis == sector_basis * restricted
        for action, restricted in zip(
            (tau_6, zeta_6), sector_actions
        )
    )
    return {
        "order": order,
        "ring_lengths": {
            "A6": len(basis_6),
            "A_boundary": len(basis_boundary),
            "B": len(basis_k),
        },
        "module_dimensions": {
            "sector": sector_basis.cols,
            "spectator": spectator_basis.cols,
            "total": total_basis.cols,
        },
        "sector_actions": {
            "tau": sector_actions[0].tolist(),
            "zeta": sector_actions[1].tolist(),
        },
        "spectator_actions": {
            "tau": spectator_actions[0].tolist(),
            "zeta": spectator_actions[1].tolist(),
        },
        "total_actions": {
            "tau": total_actions[0].tolist(),
            "zeta": total_actions[1].tolist(),
        },
        "projection": projection.tolist(),
        "splits_over_Q_tau_zeta": splits,
        "one_section": section.tolist() if section is not None else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=3)
    arguments = parser.parse_args()
    print(audit(arguments.order))


if __name__ == "__main__":
    main()
