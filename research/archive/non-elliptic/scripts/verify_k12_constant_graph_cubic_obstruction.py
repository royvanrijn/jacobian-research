#!/usr/bin/env python3
"""Exact cubic-completion obstructions for five full K12 graph families.

For z4,z9,z10,z11,z12, build all target monomials through degree three in
the other ten retained raw outputs.  Sparse finite-field elimination at the
literal parameter point selects one full-column minor and one augmented row.
The selected determinants are then reconstructed exactly over the complete
rational parameter ring.  Both are nonzero constants in every family.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp

from audit_k12_coordinate_pair_frontier import build_k12
from audit_k12_z8_cubic_completion import (
    PRIME,
    encode,
    fixed_minor_certificate,
    high_degree_group,
    multiply,
)
from search_k12_cubic_graph_bilinear_completions import linear_graph_families
from verify_k12_cubic_graph_bilinear_obstruction import (
    exact_graph_restriction,
    original_component,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_constant_graph_cubic_obstruction.json"
)
DESIRED_DEGREE = 3
MAXIMUM_TARGET_DEGREE = 3

# Zero-based source pivot -> selected retained component.
SELECTED_BAD = {
    3: 2,
    8: 1,
    9: 2,
    10: 0,
    11: 0,
}

EXPECTED = {
    3: (-sp.Integer(2**32 * 3**74 * 23), sp.Rational(2, 3)),
    8: (sp.Integer(2**8 * 3**35 * 7), sp.Integer(-1)),
    9: (sp.Integer(2**23 * 3**85 * 23), sp.Integer(7)),
    10: (-sp.Integer(2**24 * 3**84 * 7 * 23), -sp.Rational(1, 2)),
    11: (-sp.Integer(2**39 * 3**84 * 7 * 23), -sp.Rational(1, 2)),
}


def build_completion_columns(
    encoded: list[dict[tuple[int, ...], sp.Rational]],
    bad_component: int,
) -> tuple[list[dict], dict, dict[int, int], int]:
    available = [
        component
        for index, component in enumerate(encoded)
        if index != bad_component
    ]
    assert len(available) == 10
    columns = []
    nonzero_by_degree = {}
    full_basis_count = 0
    for degree in range(1, MAXIMUM_TARGET_DEGREE + 1):
        nonzero = 0
        for indices in itertools.combinations_with_replacement(
            range(10), degree
        ):
            full_basis_count += 1
            product = available[indices[0]]
            for index in indices[1:]:
                product = multiply(product, available[index])
            column = high_degree_group(product)
            if column:
                nonzero += 1
                columns.append(column)
        nonzero_by_degree[degree] = nonzero
    target = high_degree_group(encoded[bad_component])
    return columns, target, nonzero_by_degree, full_basis_count


def certify_family(
    k12: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    family,
) -> dict[str, object]:
    pivot = family.pivot
    parameters = family.parameters
    graph, retained_variables, restricted = exact_graph_restriction(
        k12, variables, family
    )
    domain = sp.QQ.frac_field(*parameters)
    component_degrees = [
        sp.Poly(component, *retained_variables, domain=domain).total_degree()
        for component in restricted
    ]
    bad = [
        index
        for index, degree in enumerate(component_degrees)
        if degree > DESIRED_DEGREE
    ]
    selected_bad = SELECTED_BAD[pivot]
    assert selected_bad in bad

    combined_variables = retained_variables + parameters
    encoded = [
        encode(component, combined_variables) for component in restricted
    ]
    columns, target, nonzero_by_degree, full_basis_count = (
        build_completion_columns(encoded, selected_bad)
    )
    all_rows = set(target)
    for column in columns:
        all_rows.update(column)

    zero = (sp.Integer(0),) * len(parameters)
    certificate, determinant = fixed_minor_certificate(
        columns, target, parameters, zero
    )
    ratio = sp.sympify(certificate["augmented_to_column_ratio"])
    expected_determinant, expected_ratio = EXPECTED[pivot]
    assert determinant == expected_determinant
    assert ratio == expected_ratio
    assert not determinant.free_symbols
    assert not sp.sympify(certificate["augmented_determinant"]).free_symbols
    assert full_basis_count == 285
    assert nonzero_by_degree[3] == 220

    return {
        "source_pivot": pivot + 1,
        "parameter_count": len(parameters),
        "normalized_linear_target_coefficients": [
            str(coefficient) for coefficient in family.coefficients
        ],
        "source_graph_degree": sp.Poly(
            graph, *retained_variables, domain=domain
        ).total_degree(),
        "source_graph_term_count": len(
            sp.Poly(graph, *retained_variables, domain=domain).terms()
        ),
        "restricted_component_degrees": component_degrees,
        "bad_original_components": [
            original_component(index, pivot) for index in bad
        ],
        "selected_obstructed_original_component": original_component(
            selected_bad, pivot
        ),
        "full_target_basis_count": full_basis_count,
        "nonzero_high_degree_columns_by_target_degree": {
            str(degree): count
            for degree, count in nonzero_by_degree.items()
        },
        "nonzero_high_degree_column_count": len(columns),
        "high_degree_row_count": len(all_rows),
        "matrix_shape": [len(all_rows), len(columns)],
        "selection_prime": PRIME,
        "certificate": certificate,
    }


def main() -> None:
    variables, k12 = build_k12()
    families = {
        family.pivot: family
        for family in linear_graph_families(k12, variables)
    }
    records = []
    for pivot in SELECTED_BAD:
        record = certify_family(k12, variables, families[pivot])
        records.append(record)
        print(
            "PASS exact cubic graph obstruction: "
            f"z{pivot + 1}, matrix {tuple(record['matrix_shape'])}, "
            "ratio "
            f"{record['certificate']['augmented_to_column_ratio']}"
        )

    artifact = {
        "format": "k12-constant-graph-cubic-obstruction-v1",
        "status": "exact bounded obstruction over Q",
        "desired_restricted_degree": DESIRED_DEGREE,
        "maximum_target_completion_degree": MAXIMUM_TARGET_DEGREE,
        "families": records,
        "proof_logic": (
            "For each family, modular sparse elimination at the literal "
            "parameter point selects a full-column row set for all nonzero "
            "linear, quadratic, and cubic completion columns and one extra "
            "augmented row. Exact reconstruction over the full rational "
            "parameter ring gives a nonzero constant column determinant and "
            "a nonzero constant augmented determinant. Thus the selected "
            "defect is outside the complete target-degree-at-most-three span "
            "at every parameter value."
        ),
        "conclusion": (
            "None of the five full normalized linear graph-coordinate "
            "families with source pivot z4,z9,z10,z11,z12 admits a one-stage "
            "target completion of degree at most three that restores source "
            "degree at most three."
        ),
        "scope": (
            "The full cubic graph families z5,z6,z7,z8, target degree at "
            "least four, nonlinear target coordinates, and ordered target "
            "stages remain open. This is not a dimension-eleven lower bound."
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print("PASS five full graph families have exact cubic obstructions")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
