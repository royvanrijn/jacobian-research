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
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from math import gcd

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import canonical_residuals  # noqa: E402
from explore_degree42_ritt_spectator_universality import (  # noqa: E402
    ALL_CUTS,
    DEGREE,
    WORD,
    W,
    build_chart,
)
from verify_degree42_ritt_relative_cotangent_cone import (  # noqa: E402
    dickson_normal_map,
    serialize_polynomial,
)

CACHE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_cellular_source_ideals.json"
)


def factor_variables_from_word():
    """Reconstruct the stable source-ring variable order without composition."""

    factor_variables = tuple(
        sp.symbols(f"x{position}_1:{degree}")
        for position, degree in enumerate(WORD)
    )
    parameters = tuple(
        variable
        for variables in factor_variables
        for variable in variables
    )
    return parameters, factor_variables


def source_ideal_data():
    """Load or construct the expensive ordinary residual equations."""

    if CACHE.is_file():
        cached = json.loads(CACHE.read_text())
        parameters, factor_variables = factor_variables_from_word()
        assert [str(parameter) for parameter in parameters] == cached[
            "parameters"
        ]
        return (
            parameters,
            factor_variables,
            cached["thick"],
            cached["boundary"],
        )

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
    cached = {
        "schema": "degree42-ritt-cellular-source-ideals.v1",
        "parameters": [str(parameter) for parameter in parameters],
        "thick": [serialize_polynomial(equation) for equation in thick],
        "boundary": [
            serialize_polynomial(equation) for equation in boundary
        ],
        "construction": (
            "ordinary polynomial residuals including degree one on the "
            "normalized 2 o 7 o 3 chart"
        ),
    }
    CACHE.write_text(json.dumps(cached, indent=2) + "\n")
    return parameters, factor_variables, cached["thick"], cached["boundary"]


def singular_data(
    order: int,
    filtration: str,
    normal_order: int | None = None,
) -> tuple[str, tuple[sp.Symbol, ...]]:
    """Return Singular output for standard bases, maps, and base actions."""

    parameters, factor_variables, thick, boundary = source_ideal_data()

    normals, base_coordinates, images = dickson_normal_map(factor_variables)
    local_variables = normals + base_coordinates
    map_images = ",".join(
        serialize_polynomial(images[parameter]) for parameter in parameters
    )
    if filtration == "maximal":
        filtration_expression = (
            f"ideal filtrationIdeal={','.join(map(str, local_variables))};\n"
            f"ideal filtrationPower=filtrationIdeal^{order};"
        )
    elif filtration == "base":
        if normal_order is None:
            raise ValueError("base filtration requires a normal cutoff")
        filtration_expression = (
            f"ideal baseIdeal={','.join(map(str, base_coordinates))};\n"
            f"ideal normalIdeal={','.join(map(str, normals))};\n"
            f"ideal filtrationPower=baseIdeal^{order}"
            f"+normalIdeal^{normal_order};"
        )
    else:
        raise ValueError(f"unknown filtration {filtration}")

    def print_reductions(
        source_basis: str,
        target_basis: str,
        marker: str,
        multiplier: str = "",
    ) -> str:
        local_name = "basis_" + marker.lower()
        return f"""
print("{marker}");
module {local_name}={source_basis};
for (int i=1; i<=ncols({local_name}); i++)
{{
  print(reduce(({multiplier})*{local_name}[i],{target_basis}));
}}
"""

    program = f"""
ring source=0,({",".join(map(str, parameters))}),dp;
ideal ITsource={",".join(thick)};
ideal IBsource={",".join(boundary)};
ring q=0,({",".join(map(str, local_variables))}),(dp({len(normals)}),dp(2));
map phi=source,{map_images};
option(redSB);
ideal IT=phi(ITsource);
ideal IB=phi(IBsource);
ideal K={",".join(map(str, normals))};
{filtration_expression}
ideal G6=std(IT+filtrationPower);
ideal GB=std(IB+filtrationPower);
ideal GK=std(K+filtrationPower);
module B6=kbase(G6);
module BB=kbase(GB);
module BK=kbase(GK);
print("BASIS_6");
for (int basisIndex6=1; basisIndex6<=ncols(B6); basisIndex6++)
{{
  print(B6[basisIndex6]);
}}
print("BASIS_BOUNDARY");
for (int basisIndexBoundary=1; basisIndexBoundary<=ncols(BB); basisIndexBoundary++)
{{
  print(BB[basisIndexBoundary]);
}}
print("BASIS_K");
for (int basisIndexK=1; basisIndexK<=ncols(BK); basisIndexK++)
{{
  print(BK[basisIndexK]);
}}
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
        if active is None or not line or line.startswith("//"):
            continue
        if line.startswith("[") and line.endswith("]"):
            line = line[1:-1]
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


def vector_space_section(projection: sp.Matrix) -> sp.Matrix:
    """Choose a rational section without imposing module linearity."""

    columns = []
    for index in range(projection.rows):
        target = sp.eye(projection.rows)[:, index]
        solution, parameters = projection.gauss_jordan_solve(target)
        solution = solution.subs(
            {parameter: sp.Integer(0) for parameter in parameters}
        )
        columns.append(solution)
    return sp.Matrix.hstack(*columns)


def primitive_integer_vector(vector: sp.Matrix) -> sp.Matrix:
    """Clear denominators and primitive-normalize a rational column."""

    rationals = [sp.Rational(value) for value in vector]
    denominator_lcm = sp.ilcm(*(value.q for value in rationals))
    integers = [int(value * denominator_lcm) for value in rationals]
    divisor = 0
    for value in integers:
        divisor = gcd(divisor, abs(value))
    divisor = divisor or 1
    integers = [value // divisor for value in integers]
    first = next((value for value in integers if value), 1)
    if first < 0:
        integers = [-value for value in integers]
    return sp.Matrix(integers)


def extension_cocycle_certificate(
    inclusion: sp.Matrix,
    projection: sp.Matrix,
    sector_actions: tuple[sp.Matrix, ...],
    total_actions: tuple[sp.Matrix, ...],
    spectator_actions: tuple[sp.Matrix, ...],
    action_names: tuple[str, ...] = ("tau", "zeta"),
) -> dict[str, object]:
    """Return an adapted-block cocycle and a nonsplitting functional."""

    section = vector_space_section(projection)
    adapted_basis = inclusion.row_join(section)
    assert adapted_basis.det() != 0
    adapted_inverse = adapted_basis.inv()
    sector_dimension = inclusion.cols
    spectator_dimension = section.cols

    couplings = []
    for sector_action, total_action, spectator_action in zip(
        sector_actions, total_actions, spectator_actions
    ):
        adapted_action = adapted_inverse * total_action * adapted_basis
        assert (
            adapted_action[:sector_dimension, :sector_dimension]
            == sector_action
        )
        assert (
            adapted_action[
                sector_dimension:,
                sector_dimension:,
            ]
            == spectator_action
        )
        assert (
            adapted_action[
                sector_dimension:,
                :sector_dimension,
            ]
            == sp.zeros(spectator_dimension, sector_dimension)
        )
        couplings.append(
            adapted_action[
                :sector_dimension,
                sector_dimension:,
            ]
        )

    h_variables = sp.symbols(
        f"h0:{sector_dimension * spectator_dimension}"
    )
    correction = sp.Matrix(
        sector_dimension,
        spectator_dimension,
        h_variables,
    )
    coboundary_expressions = []
    for sector_action, spectator_action in zip(
        sector_actions, spectator_actions
    ):
        coboundary_expressions.extend(
            sector_action * correction - correction * spectator_action
        )
    coboundary, _ = sp.linear_eq_to_matrix(
        coboundary_expressions, h_variables
    )
    cocycle = sp.Matrix(
        [
            value
            for coupling in couplings
            for value in coupling
        ]
    )
    augmented_rank = coboundary.row_join(-cocycle).rank()

    functional_terms = None
    witness_value = None
    if augmented_rank > coboundary.rank():
        witness = None
        for functional in coboundary.T.nullspace():
            value = (functional.T * cocycle)[0]
            if value:
                witness = primitive_integer_vector(functional)
                witness_value = (witness.T * cocycle)[0]
                break
        assert witness is not None
        assert (witness.T * coboundary) == sp.zeros(
            1, coboundary.cols
        )
        assert witness_value

        coordinate_names = tuple(
            f"{base}[{row},{column}]"
            for base in action_names
            for row in range(sector_dimension)
            for column in range(spectator_dimension)
        )
        functional_terms = {
            coordinate: int(value)
            for coordinate, value in zip(coordinate_names, witness)
            if value
        }
    return {
        "adapted_coupling": {
            name: coupling.tolist()
            for name, coupling in zip(action_names, couplings)
        },
        "coboundary_rank": coboundary.rank(),
        "augmented_rank": augmented_rank,
        "obstruction_functional": functional_terms,
        "obstruction_value": (
            str(sp.Rational(witness_value))
            if witness_value is not None
            else None
        ),
    }


def audit(
    order: int,
    filtration: str = "maximal",
    normal_order: int | None = None,
) -> dict[str, object]:
    """Compute the finite-jet modules and their extension."""

    output, symbols = singular_data(order, filtration, normal_order)
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
    inclusion = induced_map(
        sp.eye(len(basis_6)), sector_basis, total_basis
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
    assert projection * inclusion == sp.zeros(
        spectator_basis.cols, sector_basis.cols
    )
    sequence_homology = {
        "left_kernel": sector_basis.cols - inclusion.rank(),
        "middle": (
            total_basis.cols
            - projection.rank()
            - inclusion.rank()
        ),
        "right_cokernel": spectator_basis.cols - projection.rank(),
    }
    sequence_is_short_exact = sequence_homology == {
        "left_kernel": 0,
        "middle": 0,
        "right_cokernel": 0,
    }
    if sequence_is_short_exact:
        splits, section = splitting_solution(
            projection, total_actions, spectator_actions
        )
        extension_certificate = extension_cocycle_certificate(
            inclusion,
            projection,
            sector_actions,
            total_actions,
            spectator_actions,
        )
    else:
        splits = None
        section = None
        extension_certificate = None

    assert map_boundary_k * map_6_boundary == map_6_k
    assert all(
        action * sector_basis == sector_basis * restricted
        for action, restricted in zip(
            (tau_6, zeta_6), sector_actions
        )
    )
    return {
        "order": order,
        "filtration": filtration,
        "normal_order": normal_order,
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
        "sequence_homology": sequence_homology,
        "sequence_is_short_exact": sequence_is_short_exact,
        "splits_over_Q_tau_zeta": splits,
        "one_section": section.tolist() if section is not None else None,
        "extension_certificate": extension_certificate,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=3)
    parser.add_argument(
        "--filtration",
        choices=("maximal", "base"),
        default="maximal",
    )
    parser.add_argument(
        "--normal-order",
        type=int,
        help="normal-ideal cutoff used with --filtration base",
    )
    parser.add_argument(
        "--rebuild-source",
        action="store_true",
        help="recompute and replace the cached source residual ideals",
    )
    arguments = parser.parse_args()
    if arguments.rebuild_source and CACHE.is_file():
        CACHE.unlink()
    print(
        audit(
            arguments.order,
            arguments.filtration,
            arguments.normal_order,
        )
    )


if __name__ == "__main__":
    main()
