#!/usr/bin/env python3
"""Exact bilinear-completion obstructions for all cubic K12 graph families.

Five families have constant full-column and augmented minors.  The remaining
four are covered by exact determinant opens and lower-rank closed-stratum
column relations.  Together with the existing quadratic-graph theorem on the
closed z8 stratum, this closes one-stage bilinear completion for all nine
full normalized linear graph-coordinate families.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from audit_k12_coordinate_pair_frontier import build_k12
from audit_k12_parameterized_completion import completion_matrix
from search_k12_cubic_graph_bilinear_completions import (
    GraphFamily,
    linear_graph_families,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_cubic_graph_bilinear_obstruction.json"
)
DESIRED_DEGREE = 3

# Zero-based retained-output indices.  Each selected output alone obstructs
# completion of the whole graph restriction.
CONSTANT_SELECTED_BAD = {
    3: 2,   # pivot z4, retained K3
    8: 1,   # pivot z9, retained K2
    9: 2,   # pivot z10, retained K3
    10: 0,  # pivot z11, retained K1
    11: 0,  # pivot z12, retained K1
}

STRATIFIED_SELECTED_BAD = {
    4: 1,  # pivot z5, retained K2
    5: 1,  # pivot z6, retained K2
    6: 2,  # pivot z7, retained K3
    7: 3,  # pivot z8, retained K4
}

EXPECTED_CERTIFICATES = {
    3: (sp.Rational(3**23, 2), -sp.Rational(3**24, 2), sp.Integer(-3)),
    8: (-sp.Rational(3**7, 2), sp.Rational(3**9, 2), sp.Integer(-9)),
    9: (sp.Rational(3**25, 2), -sp.Rational(3**25, 2), sp.Integer(-1)),
    10: (-sp.Integer(3**24), sp.Rational(3**24, 2), -sp.Rational(1, 2)),
    11: (-sp.Integer(16 * 3**25), sp.Integer(8 * 3**25), -sp.Rational(1, 2)),
}


def exact_graph_restriction(
    k12: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    family: GraphFamily,
) -> tuple[sp.Expr, tuple[sp.Symbol, ...], list[sp.Expr]]:
    nonlinear = [
        sp.expand(component - variable)
        for component, variable in zip(k12, variables)
    ]
    collision_image = (
        sp.Integer(0),
        sp.Integer(0),
        -sp.Rational(1, 4),
        *([sp.Integer(0)] * 9),
    )
    coefficients = family.coefficients
    pivot = family.pivot
    slice_value = sum(
        coefficient * value
        for coefficient, value in zip(coefficients, collision_image)
    )
    graph = sp.expand(
        slice_value
        - sum(
            coefficients[index] * variables[index]
            for index in range(12)
            if index != pivot
        )
        - sum(
            coefficients[index] * nonlinear[index]
            for index in range(12)
        )
    )
    assert variables[pivot] not in graph.free_symbols
    retained_variables = tuple(
        variable
        for index, variable in enumerate(variables)
        if index != pivot
    )
    restricted = [
        sp.expand(component.subs(variables[pivot], graph))
        for index, component in enumerate(k12)
        if index != pivot
    ]
    return graph, retained_variables, restricted


def original_component(retained_component: int, pivot: int) -> int:
    return retained_component + 1 if retained_component < pivot else retained_component + 2


def determinant_with_constant_peeling(
    matrix: sp.Matrix,
    parameters: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, tuple[int, int]]:
    """Expand singleton constant rows/columns, then use the polynomial domain."""

    rows = list(range(matrix.rows))
    columns = list(range(matrix.cols))
    factor = sp.Integer(1)
    while True:
        selected = None
        for row_position, row in enumerate(rows):
            nonzero = [
                (column_position, column)
                for column_position, column in enumerate(columns)
                if matrix[row, column] != 0
            ]
            if len(nonzero) == 1:
                column_position, column = nonzero[0]
                value = matrix[row, column]
                if not value.free_symbols:
                    selected = (
                        row_position,
                        column_position,
                        row,
                        column,
                        value,
                    )
                    break
        if selected is None:
            for column_position, column in enumerate(columns):
                nonzero = [
                    (row_position, row)
                    for row_position, row in enumerate(rows)
                    if matrix[row, column] != 0
                ]
                if len(nonzero) == 1:
                    row_position, row = nonzero[0]
                    value = matrix[row, column]
                    if not value.free_symbols:
                        selected = (
                            row_position,
                            column_position,
                            row,
                            column,
                            value,
                        )
                        break
        if selected is None:
            break
        row_position, column_position, _, _, value = selected
        factor *= (-1) ** (row_position + column_position) * value
        rows.pop(row_position)
        columns.pop(column_position)

    core = matrix.extract(rows, columns)
    if not core.rows:
        return sp.factor(factor), core.shape
    ring = sp.QQ.poly_ring(*parameters)
    domain_matrix = DomainMatrix.from_Matrix(
        core, fmt="sparse"
    ).convert_to(ring)
    core_determinant = ring.to_sympy(domain_matrix.det())
    return sp.factor(factor * core_determinant), core.shape


def selected_minor_certificate(
    matrix: sp.Matrix,
    target: sp.Matrix,
    parameters: tuple[sp.Symbol, ...],
    point_values: tuple[sp.Rational, ...],
    selected_columns: list[int] | None = None,
) -> dict[str, object]:
    substitution = dict(zip(parameters, point_values))
    columns = selected_columns or list(range(matrix.cols))
    selected_matrix = matrix[:, columns]
    augmented = selected_matrix.row_join(target)
    column_rows = list(selected_matrix.subs(substitution).T.rref()[1])
    augmented_rows = list(augmented.subs(substitution).T.rref()[1])
    assert len(column_rows) == selected_matrix.cols
    assert len(augmented_rows) == selected_matrix.cols + 1
    determinant, column_core_shape = determinant_with_constant_peeling(
        selected_matrix[column_rows, :], parameters
    )
    augmented_determinant, augmented_core_shape = (
        determinant_with_constant_peeling(
            augmented[augmented_rows, :], parameters
        )
    )
    ratio = sp.factor(augmented_determinant / determinant)
    assert determinant.subs(substitution) != 0
    assert augmented_determinant.subs(substitution) != 0
    assert ratio != 0
    return {
        "selection_point": [str(value) for value in point_values],
        "selected_columns": columns,
        "column_rows": column_rows,
        "augmented_rows": augmented_rows,
        "column_core_shape_after_constant_peeling": list(column_core_shape),
        "augmented_core_shape_after_constant_peeling": list(
            augmented_core_shape
        ),
        "column_determinant": str(determinant),
        "augmented_determinant": str(augmented_determinant),
        "augmented_to_column_ratio": str(ratio),
    }


def certify_family(
    k12: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    family: GraphFamily,
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
    selected_bad = CONSTANT_SELECTED_BAD[pivot]
    assert selected_bad in bad
    matrix, target, row_monomials, full_basis_count = completion_matrix(
        restricted,
        retained_variables,
        parameters,
        selected_bad,
        maximum_target_degree=2,
    )
    augmented = matrix.row_join(target)
    zero = {parameter: sp.Integer(0) for parameter in parameters}
    column_rows = list(matrix.subs(zero).T.rref()[1])
    augmented_rows = list(augmented.subs(zero).T.rref()[1])
    assert len(column_rows) == matrix.cols
    assert len(augmented_rows) == matrix.cols + 1

    determinant = sp.factor(
        matrix[column_rows, :].det(method="domain-ge")
    )
    augmented_determinant = sp.factor(
        augmented[augmented_rows, :].det(method="domain-ge")
    )
    ratio = sp.factor(augmented_determinant / determinant)
    assert (determinant, augmented_determinant, ratio) == (
        EXPECTED_CERTIFICATES[pivot]
    )
    assert not determinant.free_symbols
    assert not augmented_determinant.free_symbols
    assert determinant != 0 and augmented_determinant != 0

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
        "nonzero_high_degree_columns": matrix.cols,
        "high_degree_row_count": matrix.rows,
        "matrix_shape": list(matrix.shape),
        "column_rows": column_rows,
        "column_row_monomials": [list(row_monomials[row]) for row in column_rows],
        "augmented_rows": augmented_rows,
        "augmented_row_monomials": [
            list(row_monomials[row]) for row in augmented_rows
        ],
        "column_determinant": str(determinant),
        "augmented_determinant": str(augmented_determinant),
        "augmented_to_column_ratio": str(ratio),
    }


GENERIC_POINT_INDICES = {
    4: [(0,), (1,)],
    5: [(0,), (1,), (2,), (2, 5)],
    6: [(0,), (1,), (2,)],
    7: [(0,), (1,), (2,)],
}


def generic_points(
    pivot: int,
    parameter_count: int,
) -> list[tuple[sp.Rational, ...]]:
    points = []
    for support in GENERIC_POINT_INDICES[pivot]:
        point = [sp.Integer(0)] * parameter_count
        if pivot == 5 and support == (2, 5):
            point[2] = sp.Integer(7)
            point[5] = sp.Integer(-3)
        else:
            for index in support:
                point[index] = sp.Integer(1)
        points.append(tuple(point))
    return points


def expected_generic_determinants(
    pivot: int,
    parameters: tuple[sp.Symbol, ...],
) -> list[tuple[sp.Expr, sp.Expr]]:
    if pivot == 4:
        return [
            (3**11 * parameters[0] ** 23, -3**13 * parameters[0] ** 23),
            (-8 * 3**11 * parameters[1] ** 23, -8 * 3**13 * parameters[1] ** 23),
        ]
    if pivot == 5:
        exceptional_factor = 3 * parameters[2] + 7 * parameters[5]
        return [
            (-sp.Rational(3**17, 8) * parameters[0] ** 17,
             sp.Rational(3**19, 8) * parameters[0] ** 17),
            (-sp.Rational(3**18, 4) * parameters[1] ** 17,
             sp.Rational(3**20, 4) * parameters[1] ** 17),
            (sp.Rational(3**20, 4) * parameters[2] ** 10 * exceptional_factor**3,
             -sp.Rational(3**22, 4) * parameters[2] ** 10 * exceptional_factor**3),
            (sp.Rational(3**19, 4) * parameters[2] * parameters[5] ** 13,
             -sp.Rational(3**21, 4) * parameters[2] * parameters[5] ** 13),
        ]
    if pivot == 6:
        return [
            (3**25 * parameters[0] ** 5, -3**25 * parameters[0] ** 5),
            (-3**26 * parameters[1] ** 5, -3**26 * parameters[1] ** 5),
            (-2 * 3**25 * parameters[2] ** 5, 2 * 3**25 * parameters[2] ** 5),
        ]
    if pivot == 7:
        return [
            (-sp.Rational(3**16, 8) * parameters[0] ** 21,
             -sp.Rational(3**16, 8) * parameters[0] ** 21),
            (sp.Rational(3**23, 4) * parameters[1] ** 16,
             sp.Rational(3**23, 4) * parameters[1] ** 16),
            (sp.Rational(3**23, 2) * parameters[2] ** 20,
             sp.Rational(3**23, 2) * parameters[2] ** 20),
        ]
    raise AssertionError(f"unexpected stratified pivot {pivot + 1}")


def closed_stratum_certificate(
    matrix: sp.Matrix,
    target: sp.Matrix,
    parameters: tuple[sp.Symbol, ...],
    vanishing_parameter_count: int,
    cover_points: list[tuple[sp.Rational, ...]],
) -> dict[str, object]:
    vanishing_parameters = parameters[:vanishing_parameter_count]
    remaining_parameters = parameters[vanishing_parameter_count:]
    substitution = {
        parameter: sp.Integer(0) for parameter in vanishing_parameters
    }
    restricted_matrix = matrix.subs(substitution)
    restricted_target = target.subs(substitution)
    reference = {
        parameter: sp.Integer(0) for parameter in remaining_parameters
    }
    basis_columns = list(restricted_matrix.subs(reference).rref()[1])
    basis = restricted_matrix[:, basis_columns]
    basis_rows = list(basis.subs(reference).T.rref()[1])
    assert len(basis_rows) == len(basis_columns)
    omitted_columns = sorted(set(range(matrix.cols)) - set(basis_columns))
    square_basis = basis[basis_rows, :]
    relations = []
    for omitted in omitted_columns:
        coefficients = sp.simplify(
            square_basis.inv() * restricted_matrix[basis_rows, omitted]
        )
        residual = restricted_matrix[:, omitted] - basis * coefficients
        assert all(sp.factor(value) == 0 for value in residual)
        relations.append(
            {
                "omitted_column": omitted,
                "basis_combination": [
                    {
                        "column": basis_columns[index],
                        "coefficient": str(sp.factor(value)),
                    }
                    for index, value in enumerate(coefficients)
                    if value
                ],
            }
        )
    certificates = [
        selected_minor_certificate(
            restricted_matrix,
            restricted_target,
            remaining_parameters,
            point,
            selected_columns=basis_columns,
        )
        for point in cover_points
    ]
    return {
        "vanishing_parameters": [str(value) for value in vanishing_parameters],
        "remaining_parameters": [str(value) for value in remaining_parameters],
        "basis_columns": basis_columns,
        "omitted_column_relations": relations,
        "cover_certificates": certificates,
    }


def certify_stratified_family(
    k12: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    family: GraphFamily,
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
    bad = [index for index, degree in enumerate(component_degrees) if degree > 3]
    selected_bad = STRATIFIED_SELECTED_BAD[pivot]
    assert selected_bad in bad
    matrix, target, _, full_basis_count = completion_matrix(
        restricted,
        retained_variables,
        parameters,
        selected_bad,
        maximum_target_degree=2,
    )
    points = generic_points(pivot, len(parameters))
    generic_certificates = [
        selected_minor_certificate(matrix, target, parameters, point)
        for point in points
    ]
    expected = expected_generic_determinants(pivot, parameters)
    assert [
        (
            sp.sympify(certificate["column_determinant"]),
            sp.sympify(certificate["augmented_determinant"]),
        )
        for certificate in generic_certificates
    ] == expected

    closed = None
    closed_dependency = None
    if pivot == 4:
        closed = closed_stratum_certificate(
            matrix,
            target,
            parameters,
            2,
            [
                (sp.Integer(0),) * 5,
                (
                    sp.Integer(0),
                    sp.Integer(0),
                    -sp.Rational(3, 7),
                    sp.Integer(0),
                    sp.Integer(0),
                ),
            ],
        )
        assert closed["omitted_column_relations"] == [
            {
                "omitted_column": 46,
                "basis_combination": [
                    {"column": 2, "coefficient": "1/3"},
                    {"column": 47, "coefficient": "2*tau2"},
                ],
            }
        ]
        expected_closed = [
            (24407490807 * (7 * parameters[4] + 3) ** 6,
             -219667417263 * (7 * parameters[4] + 3) ** 6),
            (-44641044 * parameters[4] ** 14,
             401769396 * parameters[4] ** 14),
        ]
    elif pivot == 5:
        closed = closed_stratum_certificate(
            matrix,
            target,
            parameters,
            3,
            [(sp.Integer(0),) * 5],
        )
        assert closed["omitted_column_relations"] == [
            {
                "omitted_column": 45,
                "basis_combination": [
                    {"column": 1, "coefficient": "1/3"}
                ],
            }
        ]
        expected_closed = [
            (-sp.Rational(3**18, 8), sp.Rational(3**20, 8))
        ]
    elif pivot == 6:
        closed = closed_stratum_certificate(
            matrix,
            target,
            parameters,
            3,
            [(sp.Integer(0),) * 6],
        )
        assert closed["omitted_column_relations"] == [
            {
                "omitted_column": 0,
                "basis_combination": [
                    {"column": 44, "coefficient": "tau3"},
                    {"column": 45, "coefficient": "-2*tau4"},
                ],
            }
        ]
        expected_closed = [(-sp.Integer(3**25), -sp.Integer(3**25))]
    else:
        closed_dependency = (
            "When tau0=tau1=tau2=0 this is exactly the five-parameter "
            "quadratic z8 graph family excluded over its full parameter "
            "space by BCR5."
        )
        expected_closed = []

    if closed is not None:
        actual_closed = [
            (
                sp.sympify(certificate["column_determinant"]),
                sp.sympify(certificate["augmented_determinant"]),
            )
            for certificate in closed["cover_certificates"]
        ]
        assert actual_closed == expected_closed

    return {
        "source_pivot": pivot + 1,
        "parameter_count": len(parameters),
        "normalized_linear_target_coefficients": [
            str(coefficient) for coefficient in family.coefficients
        ],
        "source_graph_degree": sp.Poly(
            graph, *retained_variables, domain=domain
        ).total_degree(),
        "restricted_component_degrees": component_degrees,
        "bad_original_components": [
            original_component(index, pivot) for index in bad
        ],
        "selected_obstructed_original_component": original_component(
            selected_bad, pivot
        ),
        "full_target_basis_count": full_basis_count,
        "matrix_shape": list(matrix.shape),
        "generic_open_certificates": generic_certificates,
        "closed_stratum_certificate": closed,
        "closed_stratum_dependency": closed_dependency,
    }


def main() -> None:
    variables, k12 = build_k12()
    families = {
        family.pivot: family
        for family in linear_graph_families(k12, variables)
    }
    records = []
    for pivot in CONSTANT_SELECTED_BAD:
        record = certify_family(k12, variables, families[pivot])
        records.append(record)
        print(
            "PASS cubic graph bilinear obstruction: "
            f"z{pivot + 1}, matrix {tuple(record['matrix_shape'])}, "
            f"ratio {record['augmented_to_column_ratio']}"
        )

    stratified_records = []
    for pivot in STRATIFIED_SELECTED_BAD:
        record = certify_stratified_family(k12, variables, families[pivot])
        stratified_records.append(record)
        print(
            "PASS stratified cubic graph bilinear obstruction: "
            f"z{pivot + 1}, matrix {tuple(record['matrix_shape'])}, "
            f"generic opens={len(record['generic_open_certificates'])}"
        )

    artifact = {
        "format": "k12-cubic-graph-bilinear-obstruction-v1",
        "status": "exact bounded obstruction over Q",
        "desired_restricted_degree": DESIRED_DEGREE,
        "maximum_target_completion_degree": 2,
        "constant_minor_families": records,
        "stratified_families": stratified_records,
        "proof_logic": (
            "Five families have constant nonzero column and augmented minors. "
            "For z5,z6,z7, monomial determinant opens cover the complement of a "
            "coordinate stratum; on that stratum exact column relations reduce "
            "the full completion span to a constant-rank basis whose augmented "
            "minors remain nonzero. For z8, three monomial opens cover the "
            "complement of tau0=tau1=tau2=0, and BCR5 excludes the resulting "
            "complete quadratic graph family on the closed stratum."
        ),
        "conclusion": (
            "None of the nine normalized linear graph-coordinate families with "
            "source pivot z4,...,z12 admits a one-stage target completion of "
            "degree at most two that restores source degree at most three."
        ),
        "scope": (
            "This is a theorem only for linear target graph coordinates and "
            "one-stage degree-at-most-two target completion. Target degree at "
            "least three on the full cubic families, nonlinear target "
            "coordinates, and ordered target stages remain open. It is not a "
            "dimension-eleven lower bound."
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print("PASS all nine full cubic graph families have exact bilinear obstructions")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
