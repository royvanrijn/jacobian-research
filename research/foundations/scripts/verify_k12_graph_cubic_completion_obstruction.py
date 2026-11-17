#!/usr/bin/env python3
"""Exact target-degree-three obstruction for all nine K12 graph families.

The five families z4,z9,z10,z11,z12 have constant full-column and augmented
minors.  The z8 family also acquires a constant minor after cubic target
columns are included.  The remaining z5,z6,z7 families are covered by exact
determinant opens and closed-stratum column relations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

from audit_k12_coordinate_pair_frontier import build_k12
from audit_k12_z8_cubic_completion import encode, fixed_minor_certificate
from search_k12_cubic_graph_bilinear_completions import linear_graph_families
from verify_k12_constant_graph_cubic_obstruction import (
    DESIRED_DEGREE,
    MAXIMUM_TARGET_DEGREE,
    SELECTED_BAD as CONSTANT_SELECTED_BAD,
    build_completion_columns,
    certify_family as certify_constant_family,
)
from verify_k12_cubic_graph_bilinear_obstruction import (
    exact_graph_restriction,
    original_component,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_graph_cubic_completion_obstruction.json"
)

STRATIFIED_SELECTED_BAD = {
    4: 1,
    5: 1,
    6: 2,
    7: 3,
}


def restrict_parameter_prefix(grouped: dict, count: int) -> dict:
    """Set an initial parameter block to zero and remove its exponents."""

    result = {}
    for source_exponents, parameter_polynomial in grouped.items():
        restricted = {
            parameter_exponents[count:]: coefficient
            for parameter_exponents, coefficient in parameter_polynomial.items()
            if not any(parameter_exponents[:count])
        }
        if restricted:
            result[source_exponents] = restricted
    return result


def add_grouped_scaled(left: dict, right: dict, scale: sp.Expr) -> dict:
    result = {
        row: dict(parameter_polynomial)
        for row, parameter_polynomial in left.items()
    }
    for row, parameter_polynomial in right.items():
        target = result.setdefault(row, {})
        for exponents, coefficient in parameter_polynomial.items():
            value = target.get(exponents, sp.Integer(0)) + scale * coefficient
            if value:
                target[exponents] = value
            else:
                target.pop(exponents, None)
        if not target:
            result.pop(row, None)
    return result


def multiply_parameter(
    grouped: dict,
    parameter_index: int,
    scale: sp.Expr = sp.Integer(1),
) -> dict:
    return {
        row: {
            tuple(
                exponent + (1 if index == parameter_index else 0)
                for index, exponent in enumerate(parameter_exponents)
            ): scale * coefficient
            for parameter_exponents, coefficient in parameter_polynomial.items()
        }
        for row, parameter_polynomial in grouped.items()
    }


def assert_expression(actual: sp.Expr, expected: sp.Expr) -> None:
    assert sp.expand(actual - expected) == 0, (actual, expected)


def checked_certificate(
    columns: list[dict],
    target: dict,
    parameters: tuple[sp.Symbol, ...],
    point: tuple[sp.Rational, ...],
    expected_determinant: sp.Expr,
    expected_ratio: sp.Expr,
    original_columns: list[int] | None = None,
) -> dict[str, object]:
    certificate, determinant = fixed_minor_certificate(
        columns, target, parameters, point
    )
    ratio = sp.sympify(certificate["augmented_to_column_ratio"])
    assert_expression(determinant, expected_determinant)
    assert_expression(ratio, expected_ratio)
    if original_columns is not None:
        certificate["selected_original_columns"] = original_columns
    return certificate


def common_record(
    family,
    graph: sp.Expr,
    retained_variables: tuple[sp.Symbol, ...],
    restricted: list[sp.Expr],
    selected_bad: int,
    columns: list[dict],
    target: dict,
    nonzero_by_degree: dict[int, int],
    full_basis_count: int,
) -> dict[str, object]:
    parameters = family.parameters
    domain = sp.QQ.frac_field(*parameters)
    component_degrees = [
        sp.Poly(component, *retained_variables, domain=domain).total_degree()
        for component in restricted
    ]
    bad = [index for index, degree in enumerate(component_degrees) if degree > 3]
    assert selected_bad in bad
    all_rows = set(target)
    for column in columns:
        all_rows.update(column)
    assert full_basis_count == 285
    assert nonzero_by_degree[3] == 220
    return {
        "source_pivot": family.pivot + 1,
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
            original_component(index, family.pivot) for index in bad
        ],
        "selected_obstructed_original_component": original_component(
            selected_bad, family.pivot
        ),
        "full_target_basis_count": full_basis_count,
        "nonzero_high_degree_columns_by_target_degree": {
            str(degree): count for degree, count in nonzero_by_degree.items()
        },
        "nonzero_high_degree_column_count": len(columns),
        "high_degree_row_count": len(all_rows),
        "matrix_shape": [len(all_rows), len(columns)],
    }


def build_family_data(k12, variables, family, selected_bad):
    graph, retained_variables, restricted = exact_graph_restriction(
        k12, variables, family
    )
    encoded = [
        encode(component, retained_variables + family.parameters)
        for component in restricted
    ]
    columns, target, nonzero_by_degree, full_basis_count = (
        build_completion_columns(encoded, selected_bad)
    )
    record = common_record(
        family,
        graph,
        retained_variables,
        restricted,
        selected_bad,
        columns,
        target,
        nonzero_by_degree,
        full_basis_count,
    )
    return columns, target, record


def certify_z5(k12, variables, family) -> dict[str, object]:
    parameters = family.parameters
    columns, target, record = build_family_data(
        k12, variables, family, STRATIFIED_SELECTED_BAD[4]
    )
    generic = [
        checked_certificate(
            columns,
            target,
            parameters,
            (sp.Integer(1),) + (sp.Integer(0),) * 6,
            -2**9 * 3**47 * 7 * parameters[0] ** 43,
            sp.Integer(9),
        ),
        checked_certificate(
            columns,
            target,
            parameters,
            (sp.Integer(0), sp.Integer(1)) + (sp.Integer(0),) * 5,
            -2**11 * 3**47 * 7 * parameters[1] ** 42,
            sp.Integer(9),
        ),
    ]

    closed_columns = [restrict_parameter_prefix(column, 2) for column in columns]
    closed_target = restrict_parameter_prefix(target, 2)
    relation = add_grouped_scaled(
        closed_columns[46], closed_columns[2], -sp.Rational(1, 3)
    )
    relation = add_grouped_scaled(
        relation, multiply_parameter(closed_columns[47], 0), -sp.Integer(2)
    )
    assert not relation
    original_columns = [index for index in range(len(columns)) if index != 46]
    basis = [closed_columns[index] for index in original_columns]
    remaining = parameters[2:]
    closed = [
        checked_certificate(
            basis,
            closed_target,
            remaining,
            (sp.Integer(0),) * 5,
            -2**15 * 3**74 * 5 * 7 * (7 * remaining[2] + 3) ** 8,
            sp.Integer(9),
            original_columns,
        ),
        checked_certificate(
            basis,
            closed_target,
            remaining,
            (
                sp.Integer(0),
                sp.Integer(0),
                -sp.Rational(3, 7),
                sp.Integer(0),
                sp.Integer(0),
            ),
            -2**10 * 3**48 * 7 * remaining[2] ** 40,
            sp.Integer(9),
            original_columns,
        ),
    ]
    record.update(
        {
            "generic_open_certificates": generic,
            "closed_stratum_certificate": {
                "vanishing_parameters": [str(value) for value in parameters[:2]],
                "remaining_parameters": [str(value) for value in remaining],
                "omitted_column_relation": (
                    "column46=(1/3)*column2+2*tau2*column47"
                ),
                "cover_certificates": closed,
            },
        }
    )
    return record


def certify_z6(k12, variables, family) -> dict[str, object]:
    parameters = family.parameters
    columns, target, record = build_family_data(
        k12, variables, family, STRATIFIED_SELECTED_BAD[5]
    )
    unit = lambda index: tuple(
        sp.Integer(1) if position == index else sp.Integer(0)
        for position in range(8)
    )
    exceptional = 3 * parameters[2] + 7 * parameters[5]
    generic = [
        checked_certificate(
            columns,
            target,
            parameters,
            unit(0),
            2**9 * 3**72 * 7 * parameters[0] ** 51,
            sp.Integer(-1),
        ),
        checked_certificate(
            columns,
            target,
            parameters,
            unit(1),
            -2**14 * 3**72 * parameters[1] ** 50,
            sp.Integer(-1),
        ),
        checked_certificate(
            columns,
            target,
            parameters,
            unit(2),
            -2**10
            * 3**75
            * 7
            * parameters[2] ** 45
            * (3 * parameters[2] - 7 * parameters[5])
            * exceptional**3
            * (6 * parameters[2] + 7 * parameters[5]),
            sp.Integer(-1),
        ),
        checked_certificate(
            columns,
            target,
            parameters,
            (
                sp.Integer(0),
                sp.Integer(0),
                sp.Integer(7),
                sp.Integer(0),
                sp.Integer(0),
                -sp.Integer(3),
                sp.Integer(0),
                sp.Integer(0),
            ),
            2**10 * 3**74 * 7 * parameters[2] ** 2 * parameters[5] ** 56,
            sp.Integer(-1),
        ),
    ]

    closed_columns = [restrict_parameter_prefix(column, 3) for column in columns]
    closed_target = restrict_parameter_prefix(target, 3)
    relation = add_grouped_scaled(
        closed_columns[45], closed_columns[1], -sp.Rational(1, 3)
    )
    assert not relation
    original_columns = [index for index in range(len(columns)) if index != 45]
    basis = [closed_columns[index] for index in original_columns]
    remaining = parameters[3:]
    closed = checked_certificate(
        basis,
        closed_target,
        remaining,
        (sp.Integer(0),) * 5,
        -2**9 * 3**59 * 7,
        sp.Integer(-1),
        original_columns,
    )
    record.update(
        {
            "generic_open_certificates": generic,
            "closed_stratum_certificate": {
                "vanishing_parameters": [str(value) for value in parameters[:3]],
                "remaining_parameters": [str(value) for value in remaining],
                "omitted_column_relation": "column45=(1/3)*column1",
                "cover_certificates": [closed],
            },
        }
    )
    return record


def certify_z7(k12, variables, family) -> dict[str, object]:
    parameters = family.parameters
    columns, target, record = build_family_data(
        k12, variables, family, STRATIFIED_SELECTED_BAD[6]
    )
    unit = lambda index: tuple(
        sp.Integer(1) if position == index else sp.Integer(0)
        for position in range(9)
    )
    generic = [
        checked_certificate(
            columns,
            target,
            parameters,
            unit(0),
            2**22 * 3**88 * parameters[0] ** 7,
            sp.Integer(7),
        ),
        checked_certificate(
            columns,
            target,
            parameters,
            unit(1),
            2**22 * 3**88 * parameters[1] ** 7,
            sp.Integer(7),
        ),
        checked_certificate(
            columns,
            target,
            parameters,
            unit(2),
            -2**24 * 3**87 * parameters[2] ** 7,
            sp.Integer(7),
        ),
    ]

    closed_columns = [restrict_parameter_prefix(column, 3) for column in columns]
    closed_target = restrict_parameter_prefix(target, 3)
    relation = add_grouped_scaled(
        closed_columns[0], multiply_parameter(closed_columns[44], 0), -1
    )
    relation = add_grouped_scaled(
        relation, multiply_parameter(closed_columns[45], 1), 2
    )
    assert not relation
    original_columns = list(range(1, len(columns)))
    basis = closed_columns[1:]
    remaining = parameters[3:]
    closed = checked_certificate(
        basis,
        closed_target,
        remaining,
        (sp.Integer(0),) * 6,
        -2**22 * 3**88,
        sp.Integer(7),
        original_columns,
    )
    record.update(
        {
            "generic_open_certificates": generic,
            "closed_stratum_certificate": {
                "vanishing_parameters": [str(value) for value in parameters[:3]],
                "remaining_parameters": [str(value) for value in remaining],
                "omitted_column_relation": (
                    "column0=tau3*column44-2*tau4*column45"
                ),
                "cover_certificates": [closed],
            },
        }
    )
    return record


def certify_z8(k12, variables, family) -> dict[str, object]:
    parameters = family.parameters
    columns, target, record = build_family_data(
        k12, variables, family, STRATIFIED_SELECTED_BAD[7]
    )
    certificate = checked_certificate(
        columns,
        target,
        parameters,
        (sp.Integer(0),) * 8,
        2**23 * 3**144,
        sp.Integer(1),
    )
    record["constant_minor_certificate"] = certificate
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    variables, k12 = build_k12()
    families = {
        family.pivot: family
        for family in linear_graph_families(k12, variables)
    }

    constant_records = []
    for pivot in CONSTANT_SELECTED_BAD:
        record = certify_constant_family(k12, variables, families[pivot])
        constant_records.append(record)
        print(
            "PASS constant cubic graph obstruction: "
            f"z{pivot + 1}, matrix {tuple(record['matrix_shape'])}"
        )

    stratified_records = []
    certifiers = {4: certify_z5, 5: certify_z6, 6: certify_z7, 7: certify_z8}
    for pivot, certifier in certifiers.items():
        record = certifier(k12, variables, families[pivot])
        stratified_records.append(record)
        certificate_count = len(record.get("generic_open_certificates", []))
        certificate_count += 1 if "constant_minor_certificate" in record else 0
        print(
            "PASS stratified cubic graph obstruction: "
            f"z{pivot + 1}, matrix {tuple(record['matrix_shape'])}, "
            f"generic certificates={certificate_count}"
        )

    artifact = {
        "format": "k12-graph-cubic-completion-obstruction-v1",
        "status": "exact bounded obstruction over Q",
        "desired_restricted_degree": DESIRED_DEGREE,
        "maximum_target_completion_degree": MAXIMUM_TARGET_DEGREE,
        "constant_minor_families": constant_records,
        "stratified_families": stratified_records,
        "proof_logic": (
            "The full completion basis consists of all nonconstant target "
            "monomials through degree three. The families z4,z8,z9,z10,z11,z12 "
            "have constant nonzero full-column and augmented minors. For z5,z6,z7, "
            "monomial determinant opens cover the complement of a coordinate "
            "stratum; exact column relations and further nonzero augmented minors "
            "cover each closed stratum."
        ),
        "conclusion": (
            "None of the nine normalized linear graph-coordinate families with "
            "source pivot z4,...,z12 admits a one-stage target completion of "
            "degree at most three that restores source degree at most three."
        ),
        "scope": (
            "This is a theorem only for linear target graph coordinates and "
            "one-stage target completion of degree at most three. Target degree "
            "at least four on the full cubic families, nonlinear target "
            "coordinates, and ordered target stages remain open. It is not a "
            "dimension-eleven lower bound."
        ),
    }
    serialized = json.dumps(artifact, indent=2) + "\n"
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.exists(), f"missing {OUTPUT.relative_to(ROOT)}"
        assert OUTPUT.read_text() == serialized, (
            f"{OUTPUT.relative_to(ROOT)} is stale; regenerate with --write"
        )
    digest = hashlib.sha256(serialized.encode()).hexdigest()
    print("PASS all nine graph families have exact cubic-completion obstructions")
    print(f"PASS checked {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
