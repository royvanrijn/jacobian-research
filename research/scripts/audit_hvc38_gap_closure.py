#!/usr/bin/env python3
"""Close bounded nonlinear-pivot and coordinated left-right HVC38 gaps.

This script extends the dimension-38 cross-construction frontier in two
directions:

1. completed nonlinear pivots for Traboulsi's d coordinate and the local
   K12 z8 coordinate, with target degree bounded by eight;
2. coordinated quadratic source shears and elementary quadratic target
   shears of K12, first linearly and then as genuine finite triangular
   automorphisms.

Everything is exact over QQ except the explicitly identified good-prime
rank calculations.  A good-prime rank lower bound, paired with an explicit
characteristic-zero kernel of the same codimension, gives the corresponding
exact rational rank.  These are bounded obstruction theorems, not global
minimality results.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import sympy as sp

from audit_hvc38_cross_construction_frontier import (
    homogeneous_coefficients,
    local_f12,
    traboulsi_g11,
)
from audit_k12_coordinate_pair_frontier import (
    PRIME,
    SparseColumnSpace,
    as_sparse_mod,
    multiply,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hvc38_gap_closure.json"
)
MAXIMUM_TARGET_DEGREE = 8

SparsePolynomial = dict[tuple[int, ...], int]
SparseVector = list[SparsePolynomial]
ColumnLabel = tuple[object, ...]


def symbolic_homogeneous_coefficients(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    degree: int,
) -> dict[tuple[int, ...], sp.Expr]:
    return {
        powers: coefficient
        for powers, coefficient in sp.Poly(expression, *variables).terms()
        if sum(powers) == degree
    }


def target_monomials(
    generators: list[SparsePolynomial],
    degree: int,
    dimension: int,
):
    one = {(0,) * dimension: 1}

    def recurse(
        remaining: int,
        start: int,
        current: SparsePolynomial,
    ):
        if remaining == 0:
            yield current
            return
        for index in range(start, len(generators)):
            yield from recurse(
                remaining - 1,
                index,
                multiply(current, generators[index]),
            )

    yield from recurse(degree, 0, one)


def bounded_pullback_filtration(
    variables: tuple[sp.Symbol, ...],
    mapping: list[sp.Expr],
    pivot: int,
    low_output_indices: tuple[int, ...],
    required_quadratic: tuple[int, ...],
) -> dict[str, object]:
    """Audit low source-degree pullbacks in a bounded target algebra."""

    allowed = [
        index
        for index, component in enumerate(mapping)
        if index != pivot and sp.diff(component, variables[pivot]) == 0
    ]
    assert set(low_output_indices).issubset(allowed)
    generators = [as_sparse_mod(mapping[index], variables) for index in allowed]
    space = SparseColumnSpace()
    cumulative_columns = 0
    degree_records = []
    for degree in range(1, MAXIMUM_TARGET_DEGREE + 1):
        degree_columns = 0
        for polynomial in target_monomials(
            generators, degree, len(variables)
        ):
            degree_columns += 1
            space.add(
                {
                    powers: coefficient
                    for powers, coefficient in polynomial.items()
                    if sum(powers) >= 3
                }
            )
        cumulative_columns += degree_columns
        nullity = cumulative_columns - space.rank
        assert nullity == len(low_output_indices)
        degree_records.append(
            {
                "target_degree": degree,
                "new_columns": degree_columns,
                "cumulative_columns": cumulative_columns,
                "rank_mod_prime": space.rank,
                "nullity": nullity,
            }
        )

    low_quadratics = [
        homogeneous_coefficients(mapping[index], variables, 2)
        for index in low_output_indices
    ]
    monomials = sorted(
        {required_quadratic}
        | set().union(*(row for row in low_quadratics))
    )
    coefficient_matrix = sp.Matrix(
        [
            [row.get(monomial, 0) for row in low_quadratics]
            for monomial in monomials
        ]
    )
    target = sp.zeros(len(monomials), 1)
    target[monomials.index(required_quadratic), 0] = 1
    rank = coefficient_matrix.rank()
    augmented_rank = coefficient_matrix.row_join(target).rank()
    assert rank == len(low_output_indices)
    assert augmented_rank == rank + 1

    return {
        "pivot": str(variables[pivot]),
        "allowed_target_outputs": [index + 1 for index in allowed],
        "maximum_target_degree": MAXIMUM_TARGET_DEGREE,
        "good_prime": PRIME,
        "degree_records": degree_records,
        "exact_low_source_degree_output_basis": [
            index + 1 for index in low_output_indices
        ],
        "quadratic_span_rank": rank,
        "quadratic_span_augmented_rank": augmented_rank,
    }


def nonlinear_pivot_audits() -> dict[str, object]:
    g_variables, g_mapping, _ = traboulsi_g11()
    k_variables, k_mapping, _ = local_f12()

    # Public pivot G7=d-xy.  If g=Y7+P and R=P(G_allowed), then on
    # g=0 one has d=xy-R and G4=a+(xy)^2-R^2.
    d_audit = bounded_pullback_filtration(
        g_variables,
        g_mapping,
        pivot=6,
        low_output_indices=(7, 9, 10),
        required_quadratic=(1, 1) + (0,) * 9,
    )
    x, y, *_ = g_variables
    a, d, auxiliary = sp.symbols("a d auxiliary")
    assert sp.expand(
        (a - d**2 + 2 * d * x * y).subs(d, x * y - auxiliary)
        - (a + (x * y) ** 2 - auxiliary**2)
    ) == 0
    d_audit["square_identity"] = "G4=a+(x*y)^2-R^2"
    d_audit["required_quadratic_part"] = "+/- x*y"

    # Local pivot K8=z8+z1*z2.  If g=Y8+P and R=P(K_allowed), then on
    # g=0 one has z8=-z1*z2-R and K4=z4+(z1*z2)^2-R^2.
    z8_audit = bounded_pullback_filtration(
        k_variables,
        k_mapping,
        pivot=7,
        low_output_indices=(6, 8, 9, 10, 11),
        required_quadratic=(1, 1) + (0,) * 10,
    )
    z1, z2 = k_variables[:2]
    z4, z8, auxiliary = sp.symbols("z4 z8 auxiliary")
    assert sp.expand(
        (z4 - 2 * z1 * z2 * z8 - z8**2).subs(
            z8, -z1 * z2 - auxiliary
        )
        - (z4 + (z1 * z2) ** 2 - auxiliary**2)
    ) == 0
    z8_audit["square_identity"] = "K4=z4+(z1*z2)^2-R^2"
    z8_audit["required_quadratic_part"] = "+/- z1*z2"

    return {
        "public_d_pivot": d_audit,
        "local_z8_pivot": z8_audit,
        "conclusion": (
            "For target degree at most eight, neither completed nonlinear "
            "pivot can have source degree at most three.  The square "
            "identity forces the completion pullback to have degree at "
            "most two and the stated missing quadratic part."
        ),
    }


def source_target_columns(
    variables: tuple[sp.Symbol, ...],
    mapping: list[sp.Expr],
) -> tuple[
    tuple[int, ...],
    list[list[sp.Expr]],
    list[SparseVector],
    list[ColumnLabel],
]:
    shear_indices = (3, 6, 8, 9, 10)
    base_indices = [
        index
        for index in range(len(variables))
        if index not in shear_indices
    ]
    quadratics = [
        variables[left] * variables[right]
        for left, right in itertools.combinations_with_replacement(
            base_indices, 2
        )
    ]
    exact_columns: list[list[sp.Expr]] = []
    modular_columns: list[SparseVector] = []
    labels: list[ColumnLabel] = []

    for index in shear_indices:
        for quadratic in quadratics:
            column = [
                sp.expand(
                    -sp.diff(component, variables[index]) * quadratic
                )
                for component in mapping
            ]
            exact_columns.append(column)
            modular_columns.append(
                [as_sparse_mod(component, variables) for component in column]
            )
            labels.append(("source", index, quadratic))

    modular_mapping = [
        as_sparse_mod(component, variables) for component in mapping
    ]
    for output in range(len(mapping)):
        other_outputs = [
            index for index in range(len(mapping)) if index != output
        ]
        for left, right in itertools.combinations_with_replacement(
            other_outputs, 2
        ):
            column = [sp.Integer(0)] * len(mapping)
            column[output] = sp.expand(mapping[left] * mapping[right])
            modular_column: SparseVector = [{} for _ in mapping]
            modular_column[output] = multiply(
                modular_mapping[left], modular_mapping[right]
            )
            exact_columns.append(column)
            modular_columns.append(modular_column)
            labels.append(("target", output, left, right))

    assert len(exact_columns) == 932
    return shear_indices, exact_columns, modular_columns, labels


def rational_reconstruction(value: int) -> sp.Rational:
    for denominator in range(1, 25):
        numerator = value * denominator % PRIME
        if numerator > PRIME // 2:
            numerator -= PRIME
        if (
            abs(numerator) <= 24
            and math.gcd(abs(numerator), denominator) == 1
        ):
            return sp.Rational(numerator, denominator)
    raise AssertionError(f"failed rational reconstruction of {value}")


def tracked_high_degree_relations(
    modular_columns: list[SparseVector],
) -> tuple[int, list[dict[int, sp.Rational]]]:
    pivots: dict[tuple[int, ...], dict[tuple[int, ...], int]] = {}
    pivot_representations: dict[tuple[int, ...], dict[int, int]] = {}
    modular_relations: list[dict[int, int]] = []

    for column_index, vector in enumerate(modular_columns):
        column = {
            (output,) + powers: coefficient
            for output, polynomial in enumerate(vector)
            for powers, coefficient in polynomial.items()
            if sum(powers) >= 4
        }
        representation = {column_index: 1}
        while column:
            pivot = min(column)
            coefficient = column[pivot]
            if pivot not in pivots:
                inverse = pow(coefficient, PRIME - 2, PRIME)
                pivots[pivot] = {
                    key: value * inverse % PRIME
                    for key, value in column.items()
                }
                pivot_representations[pivot] = {
                    key: value * inverse % PRIME
                    for key, value in representation.items()
                }
                break
            for key, value in pivots[pivot].items():
                updated = (
                    column.get(key, 0) - coefficient * value
                ) % PRIME
                if updated:
                    column[key] = updated
                else:
                    column.pop(key, None)
            for key, value in pivot_representations[pivot].items():
                updated = (
                    representation.get(key, 0) - coefficient * value
                ) % PRIME
                if updated:
                    representation[key] = updated
                else:
                    representation.pop(key, None)
        else:
            modular_relations.append(representation)

    assert len(pivots) == 896
    assert len(modular_relations) == 36
    relations = [
        {
            index: rational_reconstruction(coefficient)
            for index, coefficient in relation.items()
        }
        for relation in modular_relations
    ]
    assert max(len(relation) for relation in relations) == 5
    return len(pivots), relations


def correction_maps(
    variables: tuple[sp.Symbol, ...],
    exact_columns: list[list[sp.Expr]],
    relations: list[dict[int, sp.Rational]],
) -> list[list[sp.Expr]]:
    corrections = []
    for relation in relations:
        correction = [
            sp.expand(
                sum(
                    coefficient * exact_columns[column][output]
                    for column, coefficient in relation.items()
                )
            )
            for output in range(12)
        ]
        assert all(
            all(
                sum(powers) <= 3
                for powers, coefficient in sp.Poly(
                    component, *variables
                ).terms()
                if coefficient
            )
            for component in correction
        )
        corrections.append(correction)
    return corrections


def cubic_matrices(
    variables: tuple[sp.Symbol, ...],
    mapping: list[sp.Expr],
    corrections: list[list[sp.Expr]],
) -> tuple[sp.Matrix, list[sp.Matrix], list[tuple[int, ...]]]:
    base = [
        homogeneous_coefficients(component, variables, 3)
        for component in mapping
    ]
    correction_rows = [
        [
            homogeneous_coefficients(component, variables, 3)
            for component in correction
        ]
        for correction in corrections
    ]
    monomials = sorted(
        set().union(
            *(row for row in base),
            *(
                row
                for correction in correction_rows
                for row in correction
            ),
        )
    )
    base_matrix = sp.Matrix(
        [[row.get(monomial, 0) for monomial in monomials] for row in base]
    )
    correction_matrices = [
        sp.Matrix(
            [
                [row.get(monomial, 0) for monomial in monomials]
                for row in correction
            ]
        )
        for correction in correction_rows
    ]
    return base_matrix, correction_matrices, monomials


def linearized_rank_obstruction(
    variables: tuple[sp.Symbol, ...],
    mapping: list[sp.Expr],
    corrections: list[list[sp.Expr]],
) -> dict[str, object]:
    base, matrices, monomials = cubic_matrices(
        variables, mapping, corrections
    )
    invariant_columns = [
        column
        for column in range(len(monomials))
        if all(
            matrix[row, column] == 0
            for matrix in matrices
            for row in range(6)
        )
    ]
    block = base[:6, invariant_columns]
    pivot_positions = list(block.rref()[1])[:5]
    pivot_columns = [
        invariant_columns[position] for position in pivot_positions
    ]
    pivot_rows = list(block[:, pivot_positions].T.rref()[1])[:5]
    pivot = base.extract(pivot_rows, pivot_columns)
    assert pivot.det() == 12

    parameters = sp.symbols(f"u0:{len(matrices)}")
    family = base + sum(
        (
            parameter * matrix
            for parameter, matrix in zip(
                parameters, matrices, strict=True
            )
        ),
        sp.zeros(*base.shape),
    )
    inverse = pivot.inv()
    constant_schur_entries = []
    for row in range(family.rows):
        if row in pivot_rows:
            continue
        for column in range(family.cols):
            if column in pivot_columns:
                continue
            entry = sp.factor(
                family[row, column]
                - (
                    family.extract([row], pivot_columns)
                    * inverse
                    * family.extract(pivot_rows, [column])
                )[0]
            )
            if entry and not entry.free_symbols:
                constant_schur_entries.append(entry)
    assert -sp.Integer(2) in constant_schur_entries
    return {
        "fixed_minor_size": 5,
        "fixed_minor_determinant": str(pivot.det()),
        "constant_schur_obstruction": "-2",
        "conclusion": (
            "Every linearized source-target combination that cancels all "
            "terms above degree three still has cubic-output rank at least "
            "six."
        ),
    }


def integrate_relation(
    variables: tuple[sp.Symbol, ...],
    mapping: list[sp.Expr],
    shear_indices: tuple[int, ...],
    labels: list[ColumnLabel],
    relation: dict[int, sp.Rational],
    parameter: sp.Symbol,
) -> tuple[
    list[sp.Expr],
    int,
    int,
    int,
    dict[int, sp.Expr],
] | None:
    target_terms = [
        (coefficient, labels[column])
        for column, coefficient in relation.items()
        if labels[column][0] == "target"
    ]
    if not target_terms:
        return None
    assert len(target_terms) == 1
    target_coefficient, target_label = target_terms[0]
    assert target_coefficient == 1
    _, output, left, right = target_label
    shifts = {index: sp.Integer(0) for index in shear_indices}
    for column, coefficient in relation.items():
        label = labels[column]
        if label[0] == "source":
            _, index, quadratic = label
            shifts[index] += coefficient * quadratic
    substitutions = {
        variables[index]: variables[index] - parameter * shifts[index]
        for index in shear_indices
    }
    precomposed = [
        sp.expand(component.subs(substitutions, simultaneous=True))
        for component in mapping
    ]
    transformed = list(precomposed)
    transformed[output] = sp.expand(
        transformed[output]
        + parameter * precomposed[left] * precomposed[right]
    )
    return transformed, output, left, right, shifts


def high_degree_coefficients(
    variables: tuple[sp.Symbol, ...],
    mapping: list[sp.Expr],
) -> list[sp.Expr]:
    return [
        sp.factor(coefficient)
        for component in mapping
        for powers, coefficient in sp.Poly(component, *variables).terms()
        if sum(powers) >= 4 and coefficient
    ]


def finite_left_right_obstruction(
    variables: tuple[sp.Symbol, ...],
    mapping: list[sp.Expr],
    shear_indices: tuple[int, ...],
    labels: list[ColumnLabel],
    relations: list[dict[int, sp.Rational]],
) -> dict[str, object]:
    parameter = sp.Symbol("lambda")
    target_relations = []
    integrable = []
    nonintegrable = []
    for relation_index, relation in enumerate(relations):
        result = integrate_relation(
            variables,
            mapping,
            shear_indices,
            labels,
            relation,
            parameter,
        )
        if result is None:
            continue
        target_relations.append(relation_index)
        transformed, output, left, right, shifts = result
        high = high_degree_coefficients(variables, transformed)
        if high:
            gcd = sp.Poly(high[0], parameter)
            for coefficient in high[1:]:
                gcd = sp.gcd(gcd, sp.Poly(coefficient, parameter))
            assert sp.monic(gcd).as_expr() == parameter**2
            nonintegrable.append(relation_index)
        else:
            integrable.append(
                (
                    relation_index,
                    output,
                    left,
                    right,
                    shifts,
                )
            )
    assert len(target_relations) == 22
    assert [row[0] for row in integrable] == [
        14,
        15,
        17,
        18,
        21,
        22,
        23,
        24,
        25,
        27,
        28,
        30,
        31,
        32,
        33,
        34,
        35,
    ]
    assert nonintegrable == [16, 19, 20, 26, 29]

    parameters = sp.symbols(f"t0:{len(integrable)}")
    combined_shifts = {
        index: sp.Integer(0) for index in shear_indices
    }
    for family_parameter, (_, _, _, _, shifts) in zip(
        parameters, integrable, strict=True
    ):
        for index in shear_indices:
            combined_shifts[index] += family_parameter * shifts[index]
    substitutions = {
        variables[index]: variables[index] - combined_shifts[index]
        for index in shear_indices
    }
    precomposed = [
        sp.expand(component.subs(substitutions, simultaneous=True))
        for component in mapping
    ]
    combined = list(precomposed)
    for family_parameter, (_, output, left, right, _) in zip(
        parameters, integrable, strict=True
    ):
        assert output in (0, 1, 2)
        assert left not in (0, 1, 2) and right not in (0, 1, 2)
        combined[output] = sp.expand(
            combined[output]
            + family_parameter * precomposed[left] * precomposed[right]
        )

    high = list(
        dict.fromkeys(high_degree_coefficients(variables, combined))
    )
    assert len(high) == 48
    assert {
        sp.Poly(equation, *parameters).total_degree()
        for equation in high
    } == {2, 3}

    base_rows = [
        homogeneous_coefficients(component, variables, 3)
        for component in mapping
    ]
    base_monomials = sorted(set().union(*(row for row in base_rows)))
    base_matrix = sp.Matrix(
        [
            [row.get(monomial, 0) for monomial in base_monomials]
            for row in base_rows
        ]
    )
    pivot_positions = list(base_matrix[:6, :].rref()[1])[:6]
    assert base_matrix[:6, pivot_positions].det() == 6

    combined_rows = [
        symbolic_homogeneous_coefficients(component, variables, 3)
        for component in combined
    ]
    monomials = sorted(
        set(base_monomials) | set().union(*(row for row in combined_rows))
    )
    combined_matrix = sp.Matrix(
        [
            [row.get(monomial, 0) for monomial in monomials]
            for row in combined_rows
        ]
    )
    pivot_columns = [
        monomials.index(base_monomials[position])
        for position in pivot_positions
    ]
    minor = sp.factor(combined_matrix[:6, pivot_columns].det())
    groebner = sp.groebner(
        high + [minor], *parameters, order="grevlex"
    )
    assert len(groebner.polys) == 1
    assert groebner.polys[0].as_expr() == 1

    return {
        "target_bearing_linearized_relations": len(target_relations),
        "integrable_one_parameter_relations": [
            row[0] for row in integrable
        ],
        "nonintegrable_relation_indices": nonintegrable,
        "nonintegrable_common_parameter_factor": "lambda^2",
        "combined_family_parameters": len(parameters),
        "combined_degree_three_equations": len(high),
        "combined_degree_three_equation_degrees": [2, 3],
        "selected_cubic_minor_at_origin": "6",
        "selected_cubic_minor": str(minor),
        "degree_three_plus_rank_drop_groebner_basis": ["1"],
        "conclusion": (
            "On the exact degree-three locus of the combined seventeen-"
            "parameter triangular source-target family, the selected "
            "six-by-six cubic minor never vanishes.  Hence cubic-output "
            "rank cannot drop below six."
        ),
    }


def coordinated_left_right_audit() -> dict[str, object]:
    variables, mapping, _ = local_f12()
    (
        shear_indices,
        exact_columns,
        modular_columns,
        labels,
    ) = source_target_columns(variables, mapping)
    high_rank, relations = tracked_high_degree_relations(modular_columns)
    corrections = correction_maps(variables, exact_columns, relations)
    linearized = linearized_rank_obstruction(
        variables, mapping, corrections
    )
    finite = finite_left_right_obstruction(
        variables, mapping, shear_indices, labels, relations
    )
    return {
        "source_shear_coordinates": [
            str(variables[index]) for index in shear_indices
        ],
        "source_quadratic_columns": 140,
        "elementary_target_quadratic_columns": 792,
        "total_linearized_columns": len(exact_columns),
        "good_prime": PRIME,
        "high_degree_rank_mod_prime": high_rank,
        "high_degree_kernel_dimension": len(relations),
        "maximum_sparse_relation_support": max(
            len(relation) for relation in relations
        ),
        "linearized_cubic_rank_obstruction": linearized,
        "finite_triangular_family": finite,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    pivots = nonlinear_pivot_audits()
    left_right = coordinated_left_right_audit()
    artifact = {
        "format": "hvc38-gap-closure-v1",
        "field": "QQ",
        "nonlinear_pivots": pivots,
        "coordinated_source_target": left_right,
        "consequence": (
            "No dimension-36 construction is obtained.  The two nonlinear "
            "pivot-completion families are excluded through target degree "
            "eight, and the exact degree-three locus of the stated "
            "seventeen-parameter triangular left-right family has "
            "cubic-output rank at least six."
        ),
        "scope": (
            "The pivot statements are bounded by target degree eight.  The "
            "left-right theorem covers the explicitly derived quadratic "
            "source block and elementary quadratic target directions.  It "
            "does not exclude higher-degree completions, different source "
            "blocks, nonlinear target generators, non-nested state "
            "realizations, or non-doubling symmetric lifts."
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
        "PASS HVC38 pivots: public d and local z8 are obstructed "
        "through target degree 8"
    )
    print(
        "PASS HVC38 left-right: 36 exact high-degree kernel directions, "
        "17 integrable finite directions"
    )
    print(
        "PASS HVC38 left-right: combined degree-three/rank-drop "
        "Groebner basis is [1]"
    )
    print(f"PASS checked {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
