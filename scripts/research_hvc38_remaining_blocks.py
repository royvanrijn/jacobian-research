#!/usr/bin/env python3
"""Screen every maximal jointly affine source block of the local K12 map.

For a jointly affine block J, linearize the coordinated action consisting of

    z_j -> z_j - q_j(z_{not in J})                 (j in J)

with quadratic q_j, followed by every elementary quadratic target shear.
The screen records the ranks of the coefficient matrices in source degrees
at least four and at least three over the good prime 1,000,003.

This is an exploratory calculation.  It identifies blocks worth lifting over
QQ; by itself it is not a characteristic-zero obstruction theorem.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import shutil
import subprocess
from pathlib import Path

import sympy as sp

from audit_hvc38_cross_construction_frontier import local_f12
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
    / "hvc38_remaining_blocks_screen.json"
)

Exponent = tuple[int, ...]
SparsePolynomial = dict[Exponent, int]
SparseVector = list[SparsePolynomial]
SparseRelation = dict[int, sp.Rational]
ColumnLabel = tuple[object, ...]

# Zero-based coordinates.  These are all maximal subsets J for which every
# K12 component has total J-degree at most one.
BLOCKS = (
    (3, 5, 8, 9, 10),
    (3, 5, 8, 9, 11),
    (3, 6, 8, 9, 10),
    (3, 6, 8, 9, 11),
    (4, 6, 8, 9, 11),
    (2, 4, 6, 8, 9, 10),
)


def jointly_affine(
    mapping: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    block: tuple[int, ...],
) -> bool:
    for component in mapping:
        for powers, coefficient in sp.Poly(component, *variables).terms():
            if coefficient and sum(powers[index] for index in block) > 1:
                return False
    return True


def source_columns(
    mapping: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    block: tuple[int, ...],
) -> list[SparseVector]:
    base = [index for index in range(len(variables)) if index not in block]
    quadratics = [
        variables[left] * variables[right]
        for left, right in itertools.combinations_with_replacement(base, 2)
    ]
    return [
        [
            as_sparse_mod(
                sp.expand(-sp.diff(component, variables[index]) * quadratic),
                variables,
            )
            for component in mapping
        ]
        for index in block
        for quadratic in quadratics
    ]


def target_columns(
    mapping: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> list[SparseVector]:
    sparse_mapping = [
        as_sparse_mod(component, variables) for component in mapping
    ]
    columns: list[SparseVector] = []
    for output in range(len(mapping)):
        other_outputs = [
            index for index in range(len(mapping)) if index != output
        ]
        for left, right in itertools.combinations_with_replacement(
            other_outputs, 2
        ):
            column: SparseVector = [{} for _ in mapping]
            column[output] = multiply(
                sparse_mapping[left], sparse_mapping[right]
            )
            columns.append(column)
    assert len(columns) == 792
    return columns


def filtered_column(
    vector: SparseVector,
    cutoff: int,
) -> SparsePolynomial:
    return {
        (output,) + powers: coefficient
        for output, polynomial in enumerate(vector)
        for powers, coefficient in polynomial.items()
        if sum(powers) >= cutoff
    }


def rank(columns: list[SparseVector], cutoff: int) -> int:
    space = SparseColumnSpace()
    for column in columns:
        space.add(filtered_column(column, cutoff))
    return space.rank


def tracked_relations(
    columns: list[SparseVector],
    cutoff: int,
) -> tuple[int, list[SparseRelation]]:
    pivots: dict[tuple[int, ...], SparsePolynomial] = {}
    pivot_representations: dict[tuple[int, ...], dict[int, int]] = {}
    modular_relations: list[dict[int, int]] = []
    for column_index, vector in enumerate(columns):
        column = filtered_column(vector, cutoff)
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

    def reconstruct(value: int) -> sp.Rational:
        for denominator in range(1, 49):
            numerator = value * denominator % PRIME
            if numerator > PRIME // 2:
                numerator -= PRIME
            if (
                abs(numerator) <= 48
                and math.gcd(abs(numerator), denominator) == 1
            ):
                return sp.Rational(numerator, denominator)
        raise AssertionError(f"failed rational reconstruction of {value}")

    return len(pivots), [
        {
            index: reconstruct(coefficient)
            for index, coefficient in relation.items()
        }
        for relation in modular_relations
    ]


def exact_source_columns(
    mapping: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    block: tuple[int, ...],
) -> list[list[sp.Expr]]:
    base = [index for index in range(len(variables)) if index not in block]
    quadratics = [
        variables[left] * variables[right]
        for left, right in itertools.combinations_with_replacement(base, 2)
    ]
    return [
        [
            sp.expand(-sp.diff(component, variables[index]) * quadratic)
            for component in mapping
        ]
        for index in block
        for quadratic in quadratics
    ]


def exact_target_columns(mapping: list[sp.Expr]) -> list[list[sp.Expr]]:
    columns: list[list[sp.Expr]] = []
    for output in range(len(mapping)):
        others = [index for index in range(len(mapping)) if index != output]
        for left, right in itertools.combinations_with_replacement(others, 2):
            column = [sp.Integer(0)] * len(mapping)
            column[output] = sp.expand(mapping[left] * mapping[right])
            columns.append(column)
    return columns


def column_labels(
    variables: tuple[sp.Symbol, ...],
    block: tuple[int, ...],
) -> list[ColumnLabel]:
    base = [index for index in range(len(variables)) if index not in block]
    labels: list[ColumnLabel] = [
        ("source", index, variables[left] * variables[right])
        for index in block
        for left, right in itertools.combinations_with_replacement(base, 2)
    ]
    for output in range(len(variables)):
        others = [index for index in range(len(variables)) if index != output]
        labels.extend(
            ("target", output, left, right)
            for left, right in itertools.combinations_with_replacement(
                others, 2
            )
        )
    return labels


def exact_corrections(
    columns: list[list[sp.Expr]],
    relations: list[SparseRelation],
    variables: tuple[sp.Symbol, ...],
) -> list[list[sp.Expr]]:
    corrections = []
    for relation in relations:
        correction = [
            sp.expand(
                sum(
                    coefficient * columns[index][output]
                    for index, coefficient in relation.items()
                )
            )
            for output in range(len(variables))
        ]
        assert all(
            sp.Poly(component, *variables).total_degree() <= 3
            for component in correction
        )
        corrections.append(correction)
    return corrections


def cubic_rows(
    mapping: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> list[dict[Exponent, sp.Expr]]:
    return [
        {
            powers: coefficient
            for powers, coefficient in sp.Poly(
                component, *variables
            ).terms()
            if sum(powers) == 3
        }
        for component in mapping
    ]


def linearized_rank_witness(
    mapping: list[sp.Expr],
    corrections: list[list[sp.Expr]],
    variables: tuple[sp.Symbol, ...],
) -> dict[str, object]:
    base_rows = cubic_rows(mapping, variables)
    correction_rows = [
        cubic_rows(correction, variables) for correction in corrections
    ]
    monomials = sorted(
        set().union(
            *(row for row in base_rows),
            *(
                row
                for correction in correction_rows
                for row in correction
            ),
        )
    )
    base = sp.Matrix(
        [[row.get(monomial, 0) for monomial in monomials] for row in base_rows]
    )
    matrices = [
        sp.Matrix(
            [
                [row.get(monomial, 0) for monomial in monomials]
                for row in correction
            ]
        )
        for correction in correction_rows
    ]

    parameters = sp.symbols(f"u0:{len(matrices)}")
    family = base + sum(
        (
            parameter * matrix
            for parameter, matrix in zip(parameters, matrices, strict=True)
        ),
        sp.zeros(*base.shape),
    )

    first_six_rank = None
    best_rank = 0
    best_four_pivot: tuple[list[int], list[int], sp.Matrix] | None = None
    for selected_rows in itertools.combinations(range(base.rows), 6):
        invariant_columns = [
            column
            for column in range(len(monomials))
            if all(
                matrix[row, column] == 0
                for matrix in matrices
                for row in selected_rows
            )
        ]
        invariant_block = base.extract(selected_rows, invariant_columns)
        invariant_rank = invariant_block.rank()
        best_rank = max(best_rank, invariant_rank)
        if selected_rows == tuple(range(6)):
            first_six_rank = invariant_rank
        if invariant_rank == 4 and best_four_pivot is None:
            positions4 = list(invariant_block.rref()[1])[:4]
            columns4 = [
                invariant_columns[position] for position in positions4
            ]
            relative_rows4 = list(
                invariant_block[:, positions4].T.rref()[1]
            )[:4]
            rows4 = [
                selected_rows[position] for position in relative_rows4
            ]
            pivot4 = base.extract(rows4, columns4)
            assert pivot4.det()
            best_four_pivot = (rows4, columns4, pivot4)
        if invariant_rank >= 6:
            positions = list(invariant_block.rref()[1])[:6]
            determinant = invariant_block[:, positions].det()
            assert determinant
            return {
                "invariant_rank_on_first_six_outputs": first_six_rank,
                "selected_output_rows": [row + 1 for row in selected_rows],
                "fixed_minor_size": 6,
                "fixed_minor_determinant": str(determinant),
            }
        if invariant_rank < 5:
            continue

        positions = list(invariant_block.rref()[1])[:5]
        pivot_columns = [
            invariant_columns[position] for position in positions
        ]
        relative_pivot_rows = list(
            invariant_block[:, positions].T.rref()[1]
        )[:5]
        pivot_rows = [
            selected_rows[position] for position in relative_pivot_rows
        ]
        pivot = base.extract(pivot_rows, pivot_columns)
        assert pivot.det()
        inverse = pivot.inv()
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
                    return {
                        "invariant_rank_on_first_six_outputs": first_six_rank,
                        "selected_output_rows": [
                            selected_rows[position] + 1
                            for position in relative_pivot_rows
                        ],
                        "fixed_minor_size": 5,
                        "fixed_minor_determinant": str(pivot.det()),
                        "constant_schur_obstruction": str(entry),
                    }

    if best_four_pivot is not None:
        pivot_rows, pivot_columns, pivot = best_four_pivot
        inverse = pivot.inv()
        remaining_rows = [
            row for row in range(family.rows) if row not in pivot_rows
        ]
        remaining_columns = [
            column
            for column in range(family.cols)
            if column not in pivot_columns
        ]
        schur: dict[tuple[int, int], sp.Expr] = {}
        for row in remaining_rows:
            for column in remaining_columns:
                schur[row, column] = sp.factor(
                    family[row, column]
                    - (
                        family.extract([row], pivot_columns)
                        * inverse
                        * family.extract(pivot_rows, [column])
                    )[0]
                )
        constant_entries = [
            (row, column, entry)
            for (row, column), entry in schur.items()
            if entry and not entry.free_symbols
        ]
        for left_index, (r1, c1, _) in enumerate(constant_entries):
            for r2, c2, _ in constant_entries[left_index + 1 :]:
                if r1 == r2 or c1 == c2:
                    continue
                determinant = sp.factor(
                    schur[r1, c1] * schur[r2, c2]
                    - schur[r1, c2] * schur[r2, c1]
                )
                if determinant and not determinant.free_symbols:
                    return {
                        "invariant_rank_on_first_six_outputs": first_six_rank,
                        "maximum_invariant_rank_on_six_outputs": best_rank,
                        "fixed_minor_size": 4,
                        "fixed_minor_determinant": str(pivot.det()),
                        "constant_schur_minor_size": 2,
                        "constant_schur_minor_determinant": str(determinant),
                    }
    return {
        "invariant_rank_on_first_six_outputs": first_six_rank,
        "maximum_invariant_rank_on_six_outputs": best_rank,
        "witness": "no constant Schur entry found",
    }


def high_coefficients(
    mapping: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> list[sp.Expr]:
    return [
        sp.factor(coefficient)
        for component in mapping
        for powers, coefficient in sp.Poly(component, *variables).terms()
        if sum(powers) >= 4 and coefficient
    ]


def one_parameter_integrability(
    mapping: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    block: tuple[int, ...],
    labels: list[ColumnLabel],
    relations: list[SparseRelation],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    parameter = sp.Symbol("lambda")
    triangular_families: list[dict[str, object]] = []
    source_only = []
    target_integrable = []
    target_nonintegrable = []
    multiple_target_terms = []
    for relation_index, relation in enumerate(relations):
        shifts = {index: sp.Integer(0) for index in block}
        target_terms = []
        for column, coefficient in relation.items():
            label = labels[column]
            if label[0] == "source":
                _, index, quadratic = label
                shifts[index] += coefficient * quadratic
            else:
                target_terms.append((coefficient, label))
        substitutions = {
            variables[index]: variables[index] - parameter * shifts[index]
            for index in block
        }
        precomposed = [
            sp.expand(component.subs(substitutions, simultaneous=True))
            for component in mapping
        ]
        transformed = list(precomposed)
        if not target_terms:
            source_only.append(relation_index)
            assert not high_coefficients(transformed, variables)
            triangular_families.append(
                {
                    "relation_index": relation_index,
                    "kind": "source_only",
                    "shifts": shifts,
                }
            )
            continue
        target_outputs = {
            int(label[1]) for _, label in target_terms
        }
        target_factors = {
            int(label[position])
            for _, label in target_terms
            for position in (2, 3)
        }
        if (
            len(target_outputs) != 1
            or target_outputs & target_factors
        ):
            multiple_target_terms.append(
                {
                    "relation_index": relation_index,
                    "target_terms": [
                        {
                            "coefficient": str(coefficient),
                            "output": int(label[1]) + 1,
                            "left": int(label[2]) + 1,
                            "right": int(label[3]) + 1,
                        }
                        for coefficient, label in target_terms
                    ],
                }
            )
            continue
        output = next(iter(target_outputs))
        family = {
            "relation_index": relation_index,
            "kind": "source_target",
            "target_output": output,
            "target_terms": [
                (
                    coefficient,
                    int(label[2]),
                    int(label[3]),
                )
                for coefficient, label in target_terms
            ],
            "shifts": shifts,
        }
        triangular_families.append(family)
        for target_coefficient, target_label in target_terms:
            _, _, left, right = target_label
            transformed[output] = sp.expand(
                transformed[output]
                + parameter
                * target_coefficient
                * precomposed[left]
                * precomposed[right]
            )
        high = high_coefficients(transformed, variables)
        if high:
            gcd = sp.Poly(high[0], parameter)
            for coefficient in high[1:]:
                gcd = sp.gcd(gcd, sp.Poly(coefficient, parameter))
            target_nonintegrable.append(
                {
                    "relation_index": relation_index,
                    "common_parameter_factor": str(
                        sp.monic(gcd).as_expr()
                    ),
                }
            )
        else:
            target_integrable.append(relation_index)
    changed_outputs = sorted(
        {
            int(family["target_output"])
            for family in triangular_families
            if family["kind"] == "source_target"
        }
    )
    factor_outputs = sorted(
        {
            int(term[position])
            for family in triangular_families
            if family["kind"] == "source_target"
            for term in family["target_terms"]
            for position in (1, 2)
        }
    )
    return (
        {
            "source_only_relations": source_only,
            "target_integrable_relations": target_integrable,
            "target_nonintegrable_relations": target_nonintegrable,
            "multiple_target_term_relations": multiple_target_terms,
            "total_individually_integrable_parameters": (
                len(source_only) + len(target_integrable)
            ),
            "total_triangular_kernel_parameters": len(
                triangular_families
            ),
            "changed_target_outputs": [
                index + 1 for index in changed_outputs
            ],
            "target_factor_outputs": [
                index + 1 for index in factor_outputs
            ],
            "target_factors_avoid_changed_outputs": not (
                set(changed_outputs) & set(factor_outputs)
            ),
        },
        triangular_families,
    )


def combined_family_frontier(
    mapping: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    block: tuple[int, ...],
    families: list[dict[str, object]],
    run_groebner: bool,
) -> dict[str, object]:
    parameters = sp.symbols(f"t0:{len(families)}")
    shifts = {index: sp.Integer(0) for index in block}
    for parameter, family in zip(parameters, families, strict=True):
        family_shifts = family["shifts"]
        assert isinstance(family_shifts, dict)
        for index in block:
            shifts[index] += parameter * family_shifts[index]
    precomposed = [
        sp.expand(
            component.subs(
                {
                    variables[index]: variables[index] - shifts[index]
                    for index in block
                },
                simultaneous=True,
            )
        )
        for component in mapping
    ]
    combined = list(precomposed)
    changed_outputs = {
        int(family["target_output"])
        for family in families
        if family["kind"] == "source_target"
    }
    factor_outputs = {
        int(term[position])
        for family in families
        if family["kind"] == "source_target"
        for term in family["target_terms"]
        for position in (1, 2)
    }
    assert not (changed_outputs & factor_outputs)
    for parameter, family in zip(parameters, families, strict=True):
        if family["kind"] == "source_only":
            continue
        output = int(family["target_output"])
        for coefficient, left, right in family["target_terms"]:
            combined[output] = sp.expand(
                combined[output]
                + parameter
                * coefficient
                * precomposed[left]
                * precomposed[right]
            )

    high = list(dict.fromkeys(high_coefficients(combined, variables)))
    base_rows = cubic_rows(mapping, variables)
    base_monomials = sorted(set().union(*(row for row in base_rows)))
    base_matrix = sp.Matrix(
        [
            [row.get(monomial, 0) for monomial in base_monomials]
            for row in base_rows
        ]
    )
    pivot_positions = list(base_matrix[:6, :].rref()[1])[:6]
    assert len(pivot_positions) == 6
    base_minor = base_matrix[:6, pivot_positions].det()
    combined_rows = cubic_rows(combined, variables)
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
    initial_minor_terms = len(sp.Poly(minor, *parameters).terms())
    if initial_minor_terms > 8:
        active_base_columns = [
            column
            for column in range(base_matrix.cols)
            if any(base_matrix[row, column] for row in range(6))
        ]
        mandatory = sorted(
            {
                support[0]
                for row in range(6)
                if len(
                    support := [
                        column
                        for column in active_base_columns
                        if base_matrix[row, column]
                    ]
                )
                == 1
            }
        )
        optional = [
            column
            for column in active_base_columns
            if column not in mandatory
        ]
        best = (
            initial_minor_terms,
            sp.Poly(minor, *parameters).total_degree(),
            pivot_positions,
            minor,
        )
        for extra in itertools.combinations(
            optional, 6 - len(mandatory)
        ):
            candidate_positions = sorted(mandatory + list(extra))
            if not base_matrix[:6, candidate_positions].det():
                continue
            candidate_columns = [
                monomials.index(base_monomials[position])
                for position in candidate_positions
            ]
            candidate_minor = sp.factor(
                combined_matrix[:6, candidate_columns].det()
            )
            polynomial = sp.Poly(candidate_minor, *parameters)
            score = (
                len(polynomial.terms()),
                polynomial.total_degree(),
                candidate_positions,
                candidate_minor,
            )
            if score[:2] < best[:2]:
                best = score
        _, _, pivot_positions, minor = best
    active_high_parameters = sorted(
        {
            str(symbol)
            for equation in high
            for symbol in equation.free_symbols
            if symbol in parameters
        }
    )
    pure_power_constraints: dict[str, int] = {}
    for equation in high:
        polynomial = sp.Poly(equation, *parameters, domain=sp.QQ)
        terms = polynomial.terms()
        if len(terms) != 1:
            continue
        powers, coefficient = terms[0]
        support = [
            index for index, exponent in enumerate(powers) if exponent
        ]
        if len(support) == 1 and coefficient:
            index = support[0]
            pure_power_constraints[str(parameters[index])] = min(
                powers[index],
                pure_power_constraints.get(
                    str(parameters[index]), powers[index]
                ),
            )
    result = {
        "parameters": len(parameters),
        "high_degree_equations": len(high),
        "high_degree_equation_parameter_degrees": sorted(
            {
                sp.Poly(equation, *parameters).total_degree()
                for equation in high
            }
        ),
        "parameters_occurring_in_high_degree_equations": len(
            active_high_parameters
        ),
        "pure_power_parameter_constraints": pure_power_constraints,
        "selected_cubic_minor_at_origin": str(base_minor),
        "selected_cubic_minor": str(minor),
        "selected_cubic_minor_total_degree": sp.Poly(
            minor, *parameters
        ).total_degree(),
    }
    if run_groebner and minor.free_symbols:
        singular = shutil.which("Singular")
        assert singular is not None, "the full-kernel closure requires Singular"
        variables_string = ",".join(str(parameter) for parameter in parameters)
        integral_generators = [
            sp.Poly(equation, *parameters, domain=sp.QQ)
            .clear_denoms(convert=True)[1]
            .as_expr()
            for equation in high + [minor]
        ]
        generators = ",\n".join(
            sp.sstr(equation).replace("**", "^")
            for equation in integral_generators
        )
        program = f"""
ring coefficient_ring = 0,({variables_string}),dp;
option(redSB);
ideal rank_drop = {generators};
ideal certificate = std(rank_drop);
if (reduce(1,certificate)==0)
{{
  print("UNIT_IDEAL");
}}
else
{{
  print("NONUNIT_IDEAL");
}}
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=300,
        )
        assert "UNIT_IDEAL" in completed.stdout, (
            completed.stdout + completed.stderr + "\n" + program
        )
        assert "NONUNIT_IDEAL" not in completed.stdout, completed.stdout
        version = subprocess.run(
            [singular, "--version"],
            text=True,
            capture_output=True,
            check=True,
        ).stdout.splitlines()[0]
        result["degree_three_plus_rank_drop_groebner_basis"] = ["1"]
        result["degree_three_plus_rank_drop_unit_ideal"] = True
        result["groebner_engine"] = version
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--groebner-block",
        type=int,
        action="append",
        default=[],
        help="one-based BLOCKS entry on which to run the rank-drop Groebner test",
    )
    arguments = parser.parse_args()
    groebner_blocks = set(
        arguments.groebner_block or range(1, len(BLOCKS) + 1)
    )
    variables, mapping, _ = local_f12()
    for block_number, block in enumerate(BLOCKS, start=1):
        assert jointly_affine(mapping, variables, block)
        for candidate in range(len(variables)):
            if candidate not in block:
                assert not jointly_affine(
                    mapping, variables, tuple(sorted(block + (candidate,)))
                )

    targets = target_columns(mapping, variables)
    exact_targets = exact_target_columns(mapping)
    records = []
    for block_number, block in enumerate(BLOCKS, start=1):
        sources = source_columns(mapping, variables, block)
        columns = sources + targets
        rank_ge4 = rank(columns, 4)
        rank_ge3 = rank(columns, 3)
        tracked_rank, relations = tracked_relations(columns, 4)
        assert tracked_rank == rank_ge4
        exact_sources = exact_source_columns(mapping, variables, block)
        corrections = exact_corrections(
            exact_sources + exact_targets, relations, variables
        )
        witness = linearized_rank_witness(mapping, corrections, variables)
        labels = column_labels(variables, block)
        assert len(labels) == len(columns)
        integrability, triangular_families = one_parameter_integrability(
            mapping, variables, block, labels, relations
        )
        combined_frontier = combined_family_frontier(
            mapping,
            variables,
            block,
            triangular_families,
            block_number in groebner_blocks,
        )
        records.append(
            {
                "source_coordinates": [
                    str(variables[index]) for index in block
                ],
                "source_columns": len(sources),
                "target_columns": len(targets),
                "total_columns": len(columns),
                "rank_source_degree_at_least_4_mod_prime": rank_ge4,
                "kernel_dimension_source_degree_at_least_4": (
                    len(columns) - rank_ge4
                ),
                "rank_source_degree_at_least_3_mod_prime": rank_ge3,
                "kernel_dimension_source_degree_at_least_3": (
                    len(columns) - rank_ge3
                ),
                "linearized_cubic_correction_dimension": rank_ge3 - rank_ge4,
                "maximum_sparse_relation_support": max(
                    len(relation) for relation in relations
                ),
                "linearized_rank_witness": witness,
                "one_parameter_integrability": integrability,
                "combined_full_kernel_family": combined_frontier,
            }
        )
        print(
            tuple(index + 1 for index in block),
            f"columns={len(columns)}",
            f"rank>=4={rank_ge4}",
            f"null>=4={len(columns) - rank_ge4}",
            f"rank>=3={rank_ge3}",
            f"null>=3={len(columns) - rank_ge3}",
            f"witness={witness}",
            "kernel-parameters="
            f"{integrability['total_triangular_kernel_parameters']}",
            "source-only="
            f"{len(integrability['source_only_relations'])}",
            "target-integrable="
            f"{len(integrability['target_integrable_relations'])}",
            "target-nonintegrable="
            f"{len(integrability['target_nonintegrable_relations'])}",
        )

    artifact = {
        "format": "hvc38-remaining-blocks-screen-v1",
        "status": "experiment",
        "good_prime": PRIME,
        "records": records,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
