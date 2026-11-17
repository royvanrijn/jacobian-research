#!/usr/bin/env python3
"""Compute modular normal Kuranishi cones along the two coordinate strata.

After propagating the reduced-coordinate consequences, the 41 quadratic
obstructions leave 24 nonzero equations in 33 active variables.  The two
known coordinate P6s are therefore not components.  This script fixes a
generic point on either coordinate space, quotients its seven tangent
directions, projects the quadratic term to the Jacobian cokernel, and sends
the resulting small homogeneous ideal to Singular.
"""

from __future__ import annotations

import argparse
import subprocess

import sympy as sp
from sympy.polys.domains import GF, QQ
from sympy.polys.matrices.sdm import (
    sdm_irref,
    sdm_nullspace_from_rref,
)

from explore_rank_two_odd_mixed_quantization import essential_problem


RADICAL_COORDINATES = (5, 6, 7, 19, 34)
COORDINATE_SPACES = {
    "L1": (0, 1, 2, 3, 8, 9, 10),
    "L2": (0, 1, 9, 10, 23, 24, 35),
}


def reduce_value(value, field):
    return (
        field(int(value.numerator))
        / field(int(value.denominator))
    )


def reduced_quadrics(field):
    """Return the 24 surviving quadrics after five reduced vanishings."""

    _, _, _, _, parameter_monomials, equations = essential_problem()
    radical = set(RADICAL_COORDINATES)
    result = []
    for equation in equations.values():
        reduced = {
            parameter_monomials[column]: reduce_value(coefficient, field)
            for column, coefficient in equation.items()
            if not set(parameter_monomials[column]).intersection(radical)
        }
        reduced = {
            monomial: coefficient
            for monomial, coefficient in reduced.items()
            if coefficient
        }
        if reduced:
            result.append(reduced)
    assert len(result) == 24
    return result


def jacobian_rows(quadrics, variables, base_point, field):
    variable_index = {
        variable: index for index, variable in enumerate(variables)
    }
    rows = {}
    for row, equation in enumerate(quadrics):
        entries = {}
        for (left, right), coefficient in equation.items():
            if left == right:
                value = (
                    field(2)
                    * coefficient
                    * base_point.get(left, field.zero)
                )
                if value and left in variable_index:
                    entries[variable_index[left]] = (
                        entries.get(variable_index[left], field.zero)
                        + value
                    )
            else:
                if left in variable_index:
                    value = (
                        coefficient
                        * base_point.get(right, field.zero)
                    )
                    if value:
                        entries[variable_index[left]] = (
                            entries.get(variable_index[left], field.zero)
                            + value
                        )
                if right in variable_index:
                    value = (
                        coefficient
                        * base_point.get(left, field.zero)
                    )
                    if value:
                        entries[variable_index[right]] = (
                            entries.get(variable_index[right], field.zero)
                            + value
                        )
        entries = {
            column: value
            for column, value in entries.items()
            if value
        }
        if entries:
            rows[row] = entries
    return rows


def nullspace(rows, column_count, field):
    reduced, pivots, nonzero = sdm_irref(rows)
    basis, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        column_count,
        pivots,
        nonzero,
    )
    return len(pivots), basis


def evaluate_quadrics(quadrics, vector, variables, field):
    coordinates = {
        variable: vector.get(index, field.zero)
        for index, variable in enumerate(variables)
    }
    result = {}
    for row, equation in enumerate(quadrics):
        value = sum(
            (
                coefficient
                * coordinates.get(left, field.zero)
                * coordinates.get(right, field.zero)
                for (left, right), coefficient in equation.items()
            ),
            field.zero,
        )
        if value:
            result[row] = value
    return result


def add_vectors(left, right, field):
    return {
        index: value
        for index in set(left).union(right)
        if (
            value := left.get(index, field.zero)
            + right.get(index, field.zero)
        )
    }


def project(vector, functionals, field):
    return {
        functional_index: value
        for functional_index, functional in enumerate(functionals)
        if (
            value := sum(
                (
                    coefficient
                    * vector.get(row, field.zero)
                    for row, coefficient in functional.items()
                ),
                field.zero,
            )
        )
    }


def normal_cone(space, field):
    quadrics = reduced_quadrics(field)
    coordinate_space = COORDINATE_SPACES[space]
    active = [
        index
        for index in range(38)
        if index not in RADICAL_COORDINATES
    ]
    normal_variables = [
        index for index in active if index not in coordinate_space
    ]
    # Distinct nonzero values avoid accidental special subloci.
    base_point = {
        variable: field(position + 1)
        for position, variable in enumerate(coordinate_space)
    }
    jacobian = jacobian_rows(
        quadrics,
        normal_variables,
        base_point,
        field,
    )
    jacobian_rank, tangent_basis = nullspace(
        jacobian,
        len(normal_variables),
        field,
    )

    transpose = {}
    for row, entries in jacobian.items():
        for column, coefficient in entries.items():
            transpose.setdefault(column, {})[row] = coefficient
    _, cokernel = nullspace(transpose, len(quadrics), field)

    diagonal = [
        evaluate_quadrics(
            quadrics,
            vector,
            normal_variables,
            field,
        )
        for vector in tangent_basis
    ]
    projected_equations = [dict() for _ in cokernel]
    for left, left_vector in enumerate(tangent_basis):
        projected = project(diagonal[left], cokernel, field)
        for equation, coefficient in projected.items():
            projected_equations[equation][(left, left)] = coefficient
        for right in range(left + 1, len(tangent_basis)):
            mixed_value = evaluate_quadrics(
                quadrics,
                add_vectors(left_vector, tangent_basis[right], field),
                normal_variables,
                field,
            )
            cross = {
                row: value
                for row in set(mixed_value).union(
                    diagonal[left],
                    diagonal[right],
                )
                if (
                    value := mixed_value.get(row, field.zero)
                    - diagonal[left].get(row, field.zero)
                    - diagonal[right].get(row, field.zero)
                )
            }
            projected = project(cross, cokernel, field)
            for equation, coefficient in projected.items():
                projected_equations[equation][(left, right)] = coefficient
    projected_equations = [
        equation for equation in projected_equations if equation
    ]
    return {
        "field": field,
        "quadrics": quadrics,
        "active_variables": active,
        "normal_variable_indices": normal_variables,
        "base_point": base_point,
        "jacobian_rows": jacobian,
        "jacobian_rank": jacobian_rank,
        "normal_variables": len(normal_variables),
        "tangent_basis": tangent_basis,
        "cokernel_dimension": len(cokernel),
        "equations": projected_equations,
    }


def singular_invariants(data, prime, primary):
    variable_count = len(data["tangent_basis"])
    names = [f"u{index}" for index in range(variable_count)]
    polynomials = []
    for equation in data["equations"]:
        terms = []
        for (left, right), coefficient in equation.items():
            value = int(coefficient) % prime
            if value:
                terms.append(f"{value}*u{left}*u{right}")
        polynomials.append("+".join(terms))
    commands = [
        f"ring r={prime},({','.join(names)}),dp;",
        f"ideal I={','.join(polynomials)};",
        "option(redSB);",
        "timer=1;",
        "ideal G=slimgb(I);",
        '"GB_SIZE",size(G);',
        '"AFFINE_DIM",dim(G);',
        '"DEGREE",mult(G);',
        '"SECONDS",timer;',
        'for (int j=1; j<=size(G); j++)',
        "{",
        '  "GENERATOR",j,G[j];',
        '  "FACTORIZATION",factorize(G[j]);',
        "}",
    ]
    if primary:
        commands.extend(
            [
                'LIB "primdec.lib";',
                "list P=minAssGTZ(I);",
                '"MINIMAL_PRIMES",size(P);',
                "for (int i=1; i<=size(P); i++)",
                "{",
                '  "PRIME",i,"SIZE",size(P[i]),"DIM",dim(std(P[i])),'
                '"DEGREE",mult(std(P[i]));',
                "}",
            ]
        )
    result = subprocess.run(
        ["Singular", "-q"],
        input="\n".join(commands),
        text=True,
        capture_output=True,
        check=True,
        timeout=1200,
    )
    return result.stdout.strip()


def exact_invariants(data):
    """Return a reduced exact basis and its factorization over QQ."""

    monomials = sorted(
        set().union(*(set(equation) for equation in data["equations"]))
    )
    monomial_index = {
        monomial: index for index, monomial in enumerate(monomials)
    }
    rows = {
        row: {
            monomial_index[monomial]: coefficient
            for monomial, coefficient in equation.items()
        }
        for row, equation in enumerate(data["equations"])
    }
    reduced, pivots, _ = sdm_irref(rows)
    variables = sp.symbols(f"u0:{len(data['tangent_basis'])}")
    expressions = []
    for row in range(len(pivots)):
        expression = sp.S.Zero
        for column, coefficient in reduced[row].items():
            left, right = monomials[column]
            expression += (
                sp.Rational(
                    int(coefficient.numerator),
                    int(coefficient.denominator),
                )
                * variables[left]
                * variables[right]
            )
        expressions.append(sp.factor(expression))
    return len(pivots), expressions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("space", choices=sorted(COORDINATE_SPACES))
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--primary", action="store_true")
    parser.add_argument("--exact", action="store_true")
    arguments = parser.parse_args()

    field = QQ if arguments.exact else GF(arguments.prime)
    data = normal_cone(arguments.space, field)
    expected = {
        "L1": (12, 14, 12),
        "L2": (8, 18, 16),
    }
    assert (
        data["jacobian_rank"],
        len(data["tangent_basis"]),
        data["cokernel_dimension"],
    ) == expected[arguments.space]
    print(
        f"{arguments.space}: normal Jacobian rank "
        f"{data['jacobian_rank']}; tangent variables "
        f"{len(data['tangent_basis'])}; cokernel dimension "
        f"{data['cokernel_dimension']}; nonzero projected coordinates "
        f"{len(data['equations'])}"
    )
    if arguments.exact:
        rank, expressions = exact_invariants(data)
        print(f"EXACT_EQUATION_RANK {rank}")
        for index, expression in enumerate(expressions, 1):
            print(f"EXACT_GENERATOR {index}: {expression}")
    else:
        print(
            singular_invariants(
                data,
                arguments.prime,
                arguments.primary,
            )
        )


if __name__ == "__main__":
    main()
