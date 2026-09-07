#!/usr/bin/env python3
"""Exact target-completion audit over quadratic K12 graph families.

The linear graph-coordinate classification leaves six natural families,
with source pivots z7,...,z12, whose graph corrections are quadratic.
For every family this script obstructs arbitrary linear/quadratic target
completion of one bad output.  For the five single-defect families
z7,z9,z10,z11,z12 it also obstructs arbitrary cubic target completion.

Each obstruction is a fixed-minor certificate.  The selected column minor
has determinant Delta(a), and the corresponding augmented determinant is a
nonzero constant multiple of Delta(a).  The Delta(a) generate the unit ideal
in the graph parameters, so the certificates cover the entire family over
Qbar, hence over Q.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp

from audit_k12_coordinate_pair_frontier import build_k12


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_parameterized_completion_frontier.json"
)
DESIRED_DEGREE = 3


def quadratic_graph_family(
    k12: list[sp.Expr],
    z: tuple[sp.Symbol, ...],
    pivot: int,
) -> tuple[
    tuple[sp.Symbol, ...],
    tuple[sp.Expr, ...],
    tuple[sp.Symbol, ...],
    list[sp.Expr],
]:
    """Return the normalized general graph family with quadratic correction."""

    nonlinear = [
        sp.expand(component - variable)
        for component, variable in zip(k12, z)
    ]
    polynomials = [
        sp.Poly(component, *z, domain=sp.QQ)
        for component in nonlinear
    ]
    monomials = sorted(
        set().union(
            *[
                {
                    exponents
                    for exponents, coefficient in polynomial.terms()
                    if coefficient
                }
                for polynomial in polynomials
            ]
        )
    )
    constraints = [
        [
            polynomial.coeff_monomial(exponents)
            for polynomial in polynomials
        ]
        for exponents in monomials
        if exponents[pivot] or sum(exponents) == 3
    ]
    pivot_row = [sp.Integer(0)] * 12
    pivot_row[pivot] = sp.Integer(1)
    matrix = sp.Matrix.vstack(
        sp.Matrix(constraints), sp.Matrix([pivot_row])
    )
    target = sp.Matrix(
        [sp.Integer(0)] * len(constraints) + [sp.Integer(1)]
    )
    coefficients = tuple(next(iter(sp.linsolve((matrix, target)))))
    parameters = tuple(
        sorted(
            set().union(
                *[
                    coefficient.free_symbols
                    for coefficient in coefficients
                ]
            ),
            key=str,
        )
    )
    collision_image = (
        sp.Integer(0),
        sp.Integer(0),
        -sp.Rational(1, 4),
        *([sp.Integer(0)] * 9),
    )
    slice_value = sum(
        coefficients[index] * collision_image[index]
        for index in range(12)
    )
    graph = sp.expand(
        slice_value
        - sum(
            coefficients[index] * z[index]
            for index in range(12)
            if index != pivot
        )
        - sum(
            coefficients[index] * nonlinear[index]
            for index in range(12)
        )
    )
    retained_variables = tuple(
        variable for index, variable in enumerate(z) if index != pivot
    )
    restricted = [
        sp.expand(component.subs(z[pivot], graph))
        for index, component in enumerate(k12)
        if index != pivot
    ]
    return parameters, coefficients, retained_variables, restricted


def completion_matrix(
    restricted: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    bad_component: int,
    maximum_target_degree: int,
) -> tuple[
    sp.MutableSparseMatrix,
    sp.Matrix,
    list[tuple[int, ...]],
    int,
]:
    """Build all high-degree columns for target degree <= maximum."""

    domain = sp.QQ.frac_field(*parameters)
    available = [
        component
        for index, component in enumerate(restricted)
        if index != bad_component
    ]
    columns: list[dict[tuple[int, ...], sp.Expr]] = []
    row_monomials: set[tuple[int, ...]] = set()
    full_basis_count = 0
    for degree in range(1, maximum_target_degree + 1):
        for indices in itertools.combinations_with_replacement(
            range(len(available)), degree
        ):
            full_basis_count += 1
            polynomial = sp.Poly(
                sp.expand(
                    sp.prod(available[index] for index in indices)
                ),
                *variables,
                domain=domain,
            )
            column = {
                exponents: coefficient
                for exponents, coefficient in polynomial.terms()
                if coefficient and sum(exponents) > DESIRED_DEGREE
            }
            if column:
                columns.append(column)
                row_monomials.update(column)
    target_polynomial = sp.Poly(
        restricted[bad_component], *variables, domain=domain
    )
    target_dictionary = {
        exponents: coefficient
        for exponents, coefficient in target_polynomial.terms()
        if coefficient and sum(exponents) > DESIRED_DEGREE
    }
    row_monomials.update(target_dictionary)
    ordered_rows = sorted(row_monomials)
    row_lookup = {
        exponents: row for row, exponents in enumerate(ordered_rows)
    }
    entries: dict[tuple[int, int], sp.Expr] = {}
    target = sp.zeros(len(ordered_rows), 1)
    for exponents, coefficient in target_dictionary.items():
        target[row_lookup[exponents], 0] = coefficient
    for column_index, column in enumerate(columns):
        for exponents, coefficient in column.items():
            entries[(row_lookup[exponents], column_index)] = coefficient
    matrix = sp.MutableSparseMatrix(
        len(ordered_rows), len(columns), entries
    )
    return matrix, target, ordered_rows, full_basis_count


def fixed_minor_certificate(
    matrix: sp.MutableSparseMatrix,
    target: sp.Matrix,
    parameters: tuple[sp.Symbol, ...],
    point: tuple[sp.Rational, ...],
) -> dict[str, object]:
    substitution = dict(zip(parameters, point))
    specialized = matrix.subs(substitution)
    augmented = matrix.row_join(target)
    specialized_augmented = augmented.subs(substitution)
    column_rows = list(specialized.T.rref()[1])
    augmented_rows = list(specialized_augmented.T.rref()[1])
    assert len(column_rows) == matrix.cols
    assert len(augmented_rows) == matrix.cols + 1
    determinant = sp.factor(
        matrix[column_rows, :].det(method="domain-ge")
    )
    augmented_determinant = sp.factor(
        augmented[augmented_rows, :].det(method="domain-ge")
    )
    ratio = sp.factor(augmented_determinant / determinant)
    assert determinant.subs(substitution) != 0
    assert ratio != 0 and not ratio.free_symbols
    assert sp.expand(augmented_determinant - ratio * determinant) == 0
    return {
        "selection_point": [str(value) for value in point],
        "column_rows": column_rows,
        "augmented_rows": augmented_rows,
        "column_determinant": str(determinant),
        "augmented_determinant": str(augmented_determinant),
        "augmented_to_column_ratio": str(ratio),
    }


def audit_family(
    k12: list[sp.Expr],
    z: tuple[sp.Symbol, ...],
    pivot: int,
    maximum_target_degree: int,
) -> dict[str, object]:
    parameters, coefficients, variables, restricted = (
        quadratic_graph_family(k12, z, pivot)
    )
    domain = sp.QQ.frac_field(*parameters)
    component_degrees = [
        sp.Poly(component, *variables, domain=domain).total_degree()
        for component in restricted
    ]
    bad = [
        index
        for index, degree in enumerate(component_degrees)
        if degree > DESIRED_DEGREE
    ]
    expected_bad = {
        7: [3],
        8: [2, 3, 4],
        9: [2],
        10: [3],
        11: [1],
        12: [1],
    }
    bad_original = [
        index + 1 if index < pivot else index + 2 for index in bad
    ]
    assert bad_original == expected_bad[pivot + 1]
    selected_bad = bad[0]
    matrix, target, row_monomials, full_basis_count = completion_matrix(
        restricted,
        variables,
        parameters,
        selected_bad,
        maximum_target_degree,
    )

    zero = (sp.Integer(0),) * len(parameters)
    points = [zero]
    if pivot == 7 and maximum_target_degree == 2:
        # Four opens cover the exceptional locus of the zero-point minor.
        points.extend(
            [
                (
                    sp.Integer(0),
                    sp.Integer(0),
                    sp.Integer(1),
                    sp.Integer(0),
                    -sp.Rational(1, 7),
                ),
                (
                    sp.Integer(1),
                    -sp.Rational(1, 54),
                    sp.Integer(0),
                    sp.Integer(0),
                    sp.Integer(1),
                ),
                (
                    sp.Integer(1),
                    sp.Integer(0),
                    sp.Integer(0),
                    -sp.Rational(1, 9),
                    sp.Integer(1),
                ),
            ]
        )
    certificates = [
        fixed_minor_certificate(matrix, target, parameters, point)
        for point in points
    ]
    determinants = [
        sp.sympify(certificate["column_determinant"])
        for certificate in certificates
    ]
    cover_basis = sp.groebner(
        determinants, *parameters, order="grevlex"
    )
    cover_is_unit = cover_basis.contains(sp.Integer(1))
    assert cover_is_unit
    return {
        "source_pivot": pivot + 1,
        "maximum_target_degree": maximum_target_degree,
        "parameter_count": len(parameters),
        "linear_target_coefficients": [
            str(coefficient) for coefficient in coefficients
        ],
        "bad_original_components": bad_original,
        "selected_obstructed_component": (
            selected_bad + 1
            if selected_bad < pivot
            else selected_bad + 2
        ),
        "full_target_basis_count": full_basis_count,
        "nonzero_high_degree_columns": matrix.cols,
        "high_degree_row_count": matrix.rows,
        "matrix_shape": list(matrix.shape),
        "certificates": certificates,
        "cover_groebner_basis": [
            str(polynomial.as_expr()) for polynomial in cover_basis.polys
        ],
        "cover_ideal_is_unit": cover_is_unit,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    z, k12 = build_k12()
    quadratic_audits = []
    for pivot in range(6, 12):
        record = audit_family(k12, z, pivot, 2)
        quadratic_audits.append(record)
        print(
            "PASS quadratic completion: "
            f"pivot z{pivot + 1}, "
            f"{record['nonzero_high_degree_columns']} columns"
        )

    cubic_audits = []
    for pivot in (6, 8, 9, 10, 11):
        record = audit_family(k12, z, pivot, 3)
        cubic_audits.append(record)
        print(
            "PASS cubic completion: "
            f"pivot z{pivot + 1}, "
            f"{record['nonzero_high_degree_columns']} columns"
        )

    artifact = {
        "format": "k12-parameterized-completion-frontier-v1",
        "status": "exact bounded obstruction over Q",
        "quadratic_graph_families_quadratic_target_completion": (
            quadratic_audits
        ),
        "single_defect_families_cubic_target_completion": cubic_audits,
        "proof_logic": (
            "For every certificate, a selected full-column minor has "
            "determinant Delta(a), while the corresponding augmented "
            "determinant is a nonzero constant multiple of Delta(a). The "
            "Delta(a) for each family generate the unit ideal in the graph "
            "parameters. Hence at every parameter point at least one full "
            "column minor and its nonzero augmented minor coexist, proving "
            "that the high-degree target defect is outside the completion "
            "span over Qbar and therefore over Q."
        ),
        "conclusion": (
            "No quadratic graph-coordinate family with source pivot "
            "z7,...,z12 admits a one-stage target completion of degree at "
            "most two that restores degree at most three. For the five "
            "single-defect families z7,z9,z10,z11,z12, the same holds for "
            "target degree at most three."
        ),
        "scope": (
            "This does not cover graph families with cubic correction, "
            "target degree at least three for the z8 family, target degree "
            "at least four for the single-defect families, "
            "nonlinear target coordinates, or ordered multi-stage target "
            "automorphisms. It is not a dimension-eleven lower bound."
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
    print(
        "PASS K11: all six quadratic graph families have exact "
        "quadratic-completion cover certificates"
    )
    print(
        "PASS K11: all five single-defect families have exact "
        "cubic-completion constant-minor certificates"
    )
    print(f"PASS checked {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
