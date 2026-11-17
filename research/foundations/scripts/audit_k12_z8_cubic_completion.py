#!/usr/bin/env python3
"""Exact cubic target-completion obstruction for the z8 graph family.

The full fraction-field matrix is avoided.  Each restricted output is kept
as a sparse polynomial in two exponent blocks: eleven source variables and
five graph parameters.  Modular evaluation selects three fixed minors, and
only those minors are reconstructed and evaluated exactly over Q[a].
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import sympy as sp

from audit_k12_coordinate_pair_frontier import build_k12
from audit_k12_parameterized_completion import quadratic_graph_family


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_z8_cubic_completion_frontier.json"
)
PRIME = 1_000_003
SOURCE_COUNT = 11
PARAMETER_COUNT = 5
Exponent = tuple[int, ...]
CombinedPolynomial = dict[Exponent, sp.Rational]
ParameterPolynomial = dict[Exponent, sp.Rational]
GroupedPolynomial = dict[Exponent, ParameterPolynomial]


def encode(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> CombinedPolynomial:
    return {
        exponents: coefficient
        for exponents, coefficient in sp.Poly(
            expression, *variables, domain=sp.QQ
        ).terms()
        if coefficient
    }


def multiply(
    left: CombinedPolynomial,
    right: CombinedPolynomial,
) -> CombinedPolynomial:
    result: CombinedPolynomial = {}
    for left_exponents, left_coefficient in left.items():
        for right_exponents, right_coefficient in right.items():
            exponents = tuple(
                left_value + right_value
                for left_value, right_value in zip(
                    left_exponents, right_exponents
                )
            )
            coefficient = (
                result.get(exponents, sp.Integer(0))
                + left_coefficient * right_coefficient
            )
            if coefficient:
                result[exponents] = coefficient
            else:
                result.pop(exponents, None)
    return result


def high_degree_group(
    polynomial: CombinedPolynomial,
) -> GroupedPolynomial:
    result: GroupedPolynomial = {}
    for exponents, coefficient in polynomial.items():
        source_exponents = exponents[:SOURCE_COUNT]
        if sum(source_exponents) <= 3:
            continue
        parameter_exponents = exponents[SOURCE_COUNT:]
        parameter_polynomial = result.setdefault(source_exponents, {})
        parameter_polynomial[parameter_exponents] = (
            parameter_polynomial.get(
                parameter_exponents, sp.Integer(0)
            )
            + coefficient
        )
    return result


def residue(value: sp.Rational) -> int:
    return (
        int(value.p) % PRIME
    ) * pow(int(value.q) % PRIME, PRIME - 2, PRIME) % PRIME


def evaluate_parameter_polynomial(
    polynomial: ParameterPolynomial,
    point: tuple[sp.Rational, ...],
) -> int:
    point_mod = [
        residue(sp.Rational(value)) for value in point
    ]
    return sum(
        residue(coefficient)
        * sp.prod(
            pow(point_mod[index], exponent, PRIME)
            for index, exponent in enumerate(exponents)
        )
        for exponents, coefficient in polynomial.items()
    ) % PRIME


def add_column(
    pivots: dict[Exponent, dict[Exponent, int]],
    source: GroupedPolynomial,
    point: tuple[sp.Rational, ...],
) -> Exponent | None:
    column = {
        row: value
        for row, polynomial in source.items()
        if (
            value := evaluate_parameter_polynomial(polynomial, point)
        )
    }
    while column:
        pivot = min(column)
        coefficient = column[pivot]
        if pivot not in pivots:
            inverse = pow(coefficient, PRIME - 2, PRIME)
            pivots[pivot] = {
                row: value * inverse % PRIME
                for row, value in column.items()
            }
            return pivot
        existing = pivots[pivot]
        for row, value in existing.items():
            updated = (
                column.get(row, 0) - coefficient * value
            ) % PRIME
            if updated:
                column[row] = updated
            else:
                column.pop(row, None)
    return None


def parameter_expression(
    polynomial: ParameterPolynomial,
    parameters: tuple[sp.Symbol, ...],
) -> sp.Expr:
    return sum(
        coefficient
        * sp.prod(
            parameters[index] ** exponent
            for index, exponent in enumerate(exponents)
        )
        for exponents, coefficient in polynomial.items()
    )


def fixed_minor_certificate(
    columns: list[GroupedPolynomial],
    target: GroupedPolynomial,
    parameters: tuple[sp.Symbol, ...],
    point: tuple[sp.Rational, ...],
) -> tuple[dict[str, object], sp.Expr]:
    pivots: dict[Exponent, dict[Exponent, int]] = {}
    column_rows = [
        add_column(pivots, column, point) for column in columns
    ]
    augmented_row = add_column(pivots, target, point)
    assert None not in column_rows
    assert augmented_row is not None
    selected_rows = list(column_rows)
    augmented_rows = selected_rows + [augmented_row]
    size = len(columns)
    matrix_entries = {
        (row_index, column_index): parameter_expression(
            column[row], parameters
        )
        for column_index, column in enumerate(columns)
        for row_index, row in enumerate(selected_rows)
        if row in column
    }
    augmented_entries = {
        (row_index, column_index): parameter_expression(
            column[row], parameters
        )
        for column_index, column in enumerate(columns)
        for row_index, row in enumerate(augmented_rows)
        if row in column
    }
    for row_index, row in enumerate(augmented_rows):
        if row in target:
            augmented_entries[(row_index, size)] = parameter_expression(
                target[row], parameters
            )
    matrix = sp.MutableSparseMatrix(size, size, matrix_entries)
    augmented = sp.MutableSparseMatrix(
        size + 1, size + 1, augmented_entries
    )
    determinant = sp.factor(matrix.det(method="domain-ge"))
    augmented_determinant = sp.factor(
        augmented.det(method="domain-ge")
    )
    ratio = sp.factor(augmented_determinant / determinant)
    assert determinant.subs(dict(zip(parameters, point))) != 0
    assert ratio != 0 and not ratio.free_symbols
    assert (
        sp.expand(
            augmented_determinant
            - ratio * determinant
        )
        == 0
    )
    return (
        {
            "selection_point": [str(value) for value in point],
            "column_rows": [list(row) for row in selected_rows],
            "augmented_rows": [list(row) for row in augmented_rows],
            "column_determinant": str(determinant),
            "augmented_determinant": str(augmented_determinant),
            "augmented_to_column_ratio": str(ratio),
        },
        determinant,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    z, k12 = build_k12()
    parameters, coefficients, source_variables, restricted = (
        quadratic_graph_family(k12, z, 7)
    )
    assert len(parameters) == PARAMETER_COUNT
    assert [str(value) for value in coefficients] == [
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        "tau0",
        "1",
        "tau1",
        "tau2",
        "tau3",
        "tau4",
    ]
    combined_variables = source_variables + parameters
    encoded = [
        encode(component, combined_variables)
        for component in restricted
    ]
    # The first bad retained component is original component 2.
    bad_component = 1
    available = [
        component
        for index, component in enumerate(encoded)
        if index != bad_component
    ]
    columns: list[GroupedPolynomial] = []
    full_basis_count = 0
    for degree in (1, 2, 3):
        for indices in itertools.combinations_with_replacement(
            range(10), degree
        ):
            full_basis_count += 1
            if degree == 1:
                product = available[indices[0]]
            elif degree == 2:
                product = multiply(
                    available[indices[0]],
                    available[indices[1]],
                )
            else:
                product = multiply(
                    multiply(
                        available[indices[0]],
                        available[indices[1]],
                    ),
                    available[indices[2]],
                )
            column = high_degree_group(product)
            if column:
                columns.append(column)
    target = high_degree_group(encoded[bad_component])
    all_rows = set(target)
    for column in columns:
        all_rows.update(column)
    assert full_basis_count == 285
    assert len(columns) == 277
    assert len(all_rows) == 54_977

    points = [
        (sp.Integer(0),) * 5,
        (
            sp.Integer(0),
            sp.Integer(0),
            sp.Integer(1),
            sp.Integer(0),
            -sp.Integer(1),
        ),
        (
            sp.Integer(1),
            sp.Integer(0),
            sp.Integer(0),
            -sp.Rational(1, 18),
            sp.Integer(1),
        ),
    ]
    certificates: list[dict[str, object]] = []
    determinants: list[sp.Expr] = []
    for point in points:
        certificate, determinant = fixed_minor_certificate(
            columns, target, parameters, point
        )
        certificates.append(certificate)
        determinants.append(determinant)
    assert all(
        certificate["augmented_to_column_ratio"] == "9/7"
        for certificate in certificates
    )
    cover = sp.groebner(
        determinants, *parameters, order="grevlex"
    )
    assert cover.contains(sp.Integer(1))

    artifact = {
        "format": "k12-z8-cubic-completion-frontier-v1",
        "status": "exact bounded obstruction over Q",
        "source_pivot": 8,
        "linear_target_coefficients": [
            str(value) for value in coefficients
        ],
        "parameter_count": len(parameters),
        "bad_original_components": [2, 3, 4],
        "selected_obstructed_component": 2,
        "maximum_target_degree": 3,
        "full_target_basis_count": full_basis_count,
        "nonzero_high_degree_columns": len(columns),
        "high_degree_row_count": len(all_rows),
        "selection_prime": PRIME,
        "certificates": certificates,
        "cover_groebner_basis": [
            str(polynomial.as_expr()) for polynomial in cover.polys
        ],
        "cover_ideal_is_unit": True,
        "proof_logic": (
            "Sparse modular elimination selects three fixed full-column "
            "minors. Each corresponding exact augmented determinant is "
            "9/7 times its column determinant. The three column "
            "determinants generate the unit ideal in Q[tau0,...,tau4], so "
            "at every parameter point one certificate has full column rank "
            "and strictly larger augmented rank."
        ),
        "conclusion": (
            "The five-parameter quadratic graph family with source pivot "
            "z8 admits no one-stage target completion of degree at most "
            "three that restores degree at most three."
        ),
        "scope": (
            "Together with BCR5 this closes cubic one-stage completion for "
            "all six quadratic graph-coordinate families. It does not cover "
            "target degree at least four, cubic graph corrections, nonlinear "
            "target coordinates, or ordered multi-stage automorphisms."
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
        "PASS z8: sparse 54977-by-277 cubic completion system assembled"
    )
    print(
        "PASS z8: three exact augmented minors have constant ratio 9/7"
    )
    print(
        "PASS z8: determinant-open cover has Groebner basis [1]"
    )
    print(f"PASS checked {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
