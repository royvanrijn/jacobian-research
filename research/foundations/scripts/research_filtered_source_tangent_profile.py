#!/usr/bin/env python3
"""Modular first-order source-contact profiles for the maps F_4,F_5,F_6.

For a bounded coefficient tangent G, the source-only infinitesimal
trivializer is unique:

    V = adj(D F_N) G.

This script computes, modulo a good prime, the dimension of the subspace for
which deg(V)<=b at every cutoff b.  It therefore records exactly when the raw
bounded-box tangent space is exhausted by source vector fields of increasing
degree.  This is a filtered first-order calculation, not a stable-moduli
quotient and not a higher-Artin algebraization certificate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from research_quartic_coefficient_kuranishi import (  # noqa: E402
    VARIABLES,
    add_scalar,
    add_term,
    multiply,
    sparse_polynomial,
    tangent_echelon,
)
from verify_all_degree_coefficient_tangents import (  # noqa: E402
    explicit_seed,
    mapping_from_primitive,
    polynomial_degree,
)
from jcsearch.weighted import w  # noqa: E402


EXPECTED_TANGENT_DIMENSIONS = {4: 58, 5: 88, 6: 123}


def source_trivializer_rows(
    relations: list[dict[int, int]],
    monomials: list[tuple[int, int, int]],
    adjugate_terms: list[
        list[dict[tuple[int, int, int], int]]
    ],
    prime: int,
) -> dict[
    tuple[int, tuple[int, int, int]], dict[int, int]
]:
    """Return coefficient rows of V=adj(DF)G across a tangent basis."""
    monomial_count = len(monomials)
    rows: dict[
        tuple[int, tuple[int, int, int]], dict[int, int]
    ] = {}
    for tangent_index, relation in enumerate(relations):
        direction = [{} for _ in range(3)]
        for coefficient_index, coefficient in relation.items():
            component = coefficient_index // monomial_count
            exponent = monomials[coefficient_index % monomial_count]
            add_term(
                direction[component],
                exponent,
                coefficient,
                prime,
            )
        for output_component in range(3):
            entry: dict[tuple[int, int, int], int] = {}
            for input_component in range(3):
                product = multiply(
                    adjugate_terms[output_component][input_component],
                    direction[input_component],
                    prime,
                )
                for exponent, coefficient in product.items():
                    add_term(entry, exponent, coefficient, prime)
            for exponent, coefficient in entry.items():
                row = rows.setdefault(
                    (output_component, exponent), {}
                )
                add_scalar(row, tangent_index, coefficient, prime)
    return rows


def insert_row(
    original: dict[int, int],
    pivots: dict[int, dict[int, int]],
    prime: int,
) -> bool:
    """Insert one row into a finite-field echelon basis."""
    row = dict(original)
    while row:
        lead = max(row)
        coefficient = row[lead]
        pivot = pivots.get(lead)
        if pivot is None:
            inverse = pow(coefficient, -1, prime)
            pivots[lead] = {
                index: value * inverse % prime
                for index, value in row.items()
            }
            return True
        for index, value in pivot.items():
            add_scalar(
                row, index, -coefficient * value, prime
            )
    return False


def filtered_profile(
    tangent_dimension: int,
    rows: dict[
        tuple[int, tuple[int, int, int]], dict[int, int]
    ],
    prime: int,
) -> list[dict[str, int]]:
    """Compute dim{G:deg(adj(DF)G)<=b} for every relevant cutoff."""
    rows_by_degree: dict[int, list[dict[int, int]]] = {}
    for (_component, exponent), row in rows.items():
        rows_by_degree.setdefault(sum(exponent), []).append(row)
    maximum_degree = max(rows_by_degree, default=0)
    pivots: dict[int, dict[int, int]] = {}
    rank_above_cutoff = 0
    dimensions = {
        maximum_degree: tangent_dimension,
    }
    for degree in range(maximum_degree, 0, -1):
        for row in rows_by_degree.get(degree, []):
            if insert_row(row, pivots, prime):
                rank_above_cutoff += 1
        dimensions[degree - 1] = (
            tangent_dimension - rank_above_cutoff
        )

    breakpoints: list[dict[str, int]] = []
    previous_dimension: int | None = None
    for cutoff in range(maximum_degree + 1):
        dimension = dimensions[cutoff]
        if (
            dimension != previous_dimension
            or cutoff in (0, maximum_degree)
        ):
            breakpoints.append(
                {
                    "source_degree_cutoff": cutoff,
                    "filtered_tangent_dimension": dimension,
                    "remaining_quotient_dimension": (
                        tangent_dimension - dimension
                    ),
                }
            )
            previous_dimension = dimension
    return breakpoints


def visible_seed_source_rows(
    degree: int,
    primitive: sp.Expr,
    adjugate: sp.Matrix,
    prime: int,
) -> dict[
    tuple[int, tuple[int, int, int]], dict[int, int]
]:
    """Coefficient rows of the unique source gauges of seed directions."""
    parameter = sp.symbols("filtered_seed_parameter")
    rows: dict[
        tuple[int, tuple[int, int, int]], dict[int, int]
    ] = {}
    for seed_index in range(degree - 3):
        family = mapping_from_primitive(
            primitive
            + parameter
            * w ** (seed_index + 2)
            * (w - 1) ** 2
        )
        direction = sp.Matrix(
            [
                sp.diff(component, parameter).subs(parameter, 0)
                for component in family
            ]
        )
        trivializer = adjugate * direction
        for component in range(3):
            polynomial = sparse_polynomial(
                sp.expand(trivializer[component]), prime
            )
            for exponent, coefficient in polynomial.items():
                row = rows.setdefault(
                    (component, exponent), {}
                )
                add_scalar(
                    row, seed_index, coefficient, prime
                )
    return rows


def audit(degree: int, prime: int) -> dict[str, object]:
    primitive = explicit_seed(degree)
    mapping = mapping_from_primitive(primitive)
    coefficient_degree = max(
        polynomial_degree(component) for component in mapping
    )
    jacobian = mapping.jacobian(VARIABLES)
    assert sp.factor(jacobian.det()) == 1
    (
        monomials,
        _free_columns,
        relations,
        _linear_pivots,
        _pivot_combinations,
        adjugate_terms,
    ) = tangent_echelon(
        jacobian.adjugate(), coefficient_degree, prime
    )
    expected = EXPECTED_TANGENT_DIMENSIONS.get(degree)
    if expected is not None:
        assert len(relations) == expected
    rows = source_trivializer_rows(
        relations, monomials, adjugate_terms, prime
    )
    profile = filtered_profile(len(relations), rows, prime)
    seed_rows = visible_seed_source_rows(
        degree, primitive, jacobian.adjugate(), prime
    )
    seed_profile = filtered_profile(
        degree - 3, seed_rows, prime
    )
    return {
        "N": degree,
        "coefficient_degree": coefficient_degree,
        "prime": prime,
        "tangent_dimension": len(relations),
        "source_trivializer_maximum_degree": profile[-1][
            "source_degree_cutoff"
        ],
        "profile_breakpoints": profile,
        "visible_seed_profile_breakpoints": seed_profile,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "degrees", nargs="*", type=int, default=(4, 5, 6)
    )
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--json-output", type=Path)
    args = parser.parse_args()
    if not sp.isprime(args.prime):
        parser.error("--prime must be prime")
    summaries = [
        audit(degree, args.prime) for degree in args.degrees
    ]
    output = {
        "status": "modular filtered first-order computation",
        "interpretation": (
            "The quotient is source-only and degree-filtered. "
            "The unrestricted first-order quotient is zero at the final "
            "cutoff; no stable-moduli claim is made."
        ),
        "summaries": summaries,
    }
    print(json.dumps(output, indent=2, sort_keys=True))
    if args.json_output is not None:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps(output, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
