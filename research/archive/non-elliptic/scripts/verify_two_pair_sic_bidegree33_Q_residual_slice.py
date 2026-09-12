#!/usr/bin/env python3
"""Verify the exact one-parameter (3,3) Q-border residual slice.

The cubic factor of the ``mu_6`` norm is only a candidate survivor.  This
script specializes the exact multiplication matrices for ``mu_6``,
``mu_7``, and ``mu_8`` to every branch above that factor.  A nonzero
four-by-four minor proves that their images span the complete length-four
fiber algebra and hence that the corresponding corrected moment ideal is
the unit ideal.

If an actual one-dimensional common quotient is found instead, the script
continues with the higher corrected-moment checks and exports its exact
coordinates.  Thus the two outcomes are kept logically distinct.
"""

from __future__ import annotations

from itertools import combinations
import argparse
import json
from pathlib import Path

import sympy as sp

from research_two_pair_sic_bidegree33_t0_Q_residual import (
    ExtensionArithmetic,
    determinant_over_extension,
    evaluate_fraction_field_element_at_algebraic_root,
    specialize_extension_element_to_number_field,
)
from verify_two_pair_sic_bidegree33_corrected_boundary import ROOT


CHECKPOINT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_t0_stratum_Q_residual_"
        "slice_s1_0_ell_0_matrices678_exact.json"
    )
)
NORM_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_t0_stratum_Q_residual_"
        "slice_s1_0_ell_0_exact.json"
    )
)
PENCIL_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_t0_stratum_Q_residual_"
        "slice_s1_0_ell_0_pencil_exact.json"
    )
)
HIGHER_ARTIFACTS = {
    order: (
        ROOT
        / "artifacts"
        / "generated-results"
        / (
            "two_pair_sic_bidegree33_t0_stratum_Q_residual_"
            f"slice_s1_0_ell_0_cubic_mu{order}_exact.json"
        )
    )
    for order in (9, 10, 11, 12, 14)
}
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_Q_cubic_exceptional_factor.json"
)
REPRODUCTION_COMMAND = (
    "scripts/verify_two_pair_sic_bidegree33_Q_residual_slice.py"
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT)
    parser.add_argument("--factor-degree", type=int, default=3)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    return parser.parse_args()


def parse_matrix(serialized, arithmetic, locals_map):
    return [
        [
            arithmetic.make(sp.sympify(entry, locals=locals_map))
            for entry in row
        ]
        for row in serialized
    ]


def specialize_matrix(matrix, source, target):
    return [
        [
            specialize_extension_element_to_number_field(
                entry,
                source,
                target,
            )
            for entry in row
        ]
        for row in matrix
    ]


def dot(arithmetic, left, right):
    answer = arithmetic.zero
    for left_entry, right_entry in zip(left, right, strict=True):
        answer = arithmetic.add(
            answer,
            arithmetic.mul(left_entry, right_entry),
        )
    return answer


def first_rank_three_columns(arithmetic, columns):
    for selected_columns in combinations(range(len(columns)), 3):
        selected = [columns[index] for index in selected_columns]
        for selected_rows in combinations(range(4), 3):
            minor = [
                [selected[column][row] for column in range(3)]
                for row in selected_rows
            ]
            determinant = determinant_over_extension(arithmetic, minor)
            if not determinant.is_zero:
                return selected_columns, selected
    raise AssertionError("the common image has rank below three")


def first_rank_four_columns(arithmetic, columns):
    for selected_columns in combinations(range(len(columns)), 4):
        selected = [columns[index] for index in selected_columns]
        matrix = [
            [selected[column][row] for column in range(4)]
            for row in range(4)
        ]
        determinant = determinant_over_extension(arithmetic, matrix)
        if not determinant.is_zero:
            return selected_columns, determinant
    return None


def left_kernel_cofactor(arithmetic, columns):
    values = []
    for deleted_row in range(4):
        retained_rows = [row for row in range(4) if row != deleted_row]
        minor = [
            [columns[column][row] for column in range(3)]
            for row in retained_rows
        ]
        determinant = determinant_over_extension(arithmetic, minor)
        values.append(
            determinant
            if deleted_row % 2 == 0
            else arithmetic.neg(determinant)
        )
    assert all(dot(arithmetic, values, column).is_zero for column in columns)
    return values


def main() -> None:
    arguments = parse_arguments()
    checkpoint_path = arguments.checkpoint
    if not checkpoint_path.is_absolute():
        checkpoint_path = ROOT / checkpoint_path
    payload = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    checkpoint = payload["matrix_checkpoint"]

    parameter = sp.symbols(checkpoint["base_parameter"])
    root = sp.symbols(checkpoint["extension_symbol"])
    base = sp.QQ.frac_field(parameter)
    locals_map = {str(parameter): parameter, str(root): root}
    arithmetic = ExtensionArithmetic(
        base,
        root,
        sp.sympify(checkpoint["modulus"], locals=locals_map),
    )
    matrices = [
        parse_matrix(checkpoint[name], arithmetic, locals_map)
        for name in ("M6", "M7", "M8")
    ]

    norm_payload = json.loads(NORM_ARTIFACT.read_text(encoding="utf-8"))
    pencil_payload = json.loads(
        PENCIL_ARTIFACT.read_text(encoding="utf-8")
    )
    pencil_tests = pencil_payload["mu6_mu7_pencil"][
        "mu6_numerator_factor_tests"
    ]
    assert sorted(record["degree"] for record in pencil_tests) == [3, 100]
    degree_100_test = next(
        record for record in pencil_tests if record["degree"] == 100
    )
    assert degree_100_test["excluded_by_c1"] is True
    cubic_pencil_test = next(
        record for record in pencil_tests if record["degree"] == 3
    )
    assert cubic_pencil_test["excluded_by_c1"] is False
    records = []
    for record in norm_payload["mu6_norm"]["numerator_factors"]:
        factor = sp.Poly(
            sp.sympify(record["factor"], locals=locals_map),
            parameter,
            domain=sp.QQ,
        )
        if factor.degree() == arguments.factor_degree:
            records.append((record, factor))
    assert len(records) == 1
    factor_record, factor = records[0]
    assert sp.factor_list(factor.as_expr())[1] == [(factor.as_expr(), 1)]

    number_field = sp.QQ.alg_field_from_poly(
        factor,
        alias=f"theta_{factor.degree()}",
    )
    specialized_modulus = sp.Poly.from_dict(
        {
            (power,): evaluate_fraction_field_element_at_algebraic_root(
                arithmetic.base.convert(arithmetic.modulus.nth(power)),
                arithmetic.base,
                number_field,
            )
            for power in range(arithmetic.degree + 1)
            if arithmetic.modulus.nth(power) != 0
        },
        root,
        domain=number_field,
    )
    modulus_factorization = specialized_modulus.factor_list()[1]
    modulus_factor_degrees = sorted(
        polynomial.degree()
        for polynomial, multiplicity in modulus_factorization
        for _ in range(multiplicity)
    )
    # The Kummer quartic splits on the cubic norm component.  Work on one
    # exact linear branch; this lowers the counterexample field from the
    # apparent degree twelve tower to a cubic number field.
    assert modulus_factor_degrees == [1, 1, 1, 1]
    selected_branch = None
    branch_checks = []
    for branch_index, (branch_modulus, multiplicity) in enumerate(
        modulus_factorization
    ):
        assert branch_modulus.degree() == 1
        if branch_modulus.nth(0) == 0:
            # The Q-border chart uses u=s0^{-1}; its coordinate ring is
            # localized at u.  A branch supported at u=0 disappears after
            # that localization, including all of its nilpotent thickness.
            branch_checks.append(
                {
                    "branch_index": branch_index,
                    "branch_modulus": sp.sstr(branch_modulus.as_expr()),
                    "multiplicity": multiplicity,
                    "chart_admissible": False,
                    "reason": "supported_on_u=0_but_chart_localizes_at_u",
                }
            )
            continue
        candidate_field = ExtensionArithmetic(
            number_field,
            root,
            branch_modulus.as_expr(),
        )
        candidate_matrices = [
            specialize_matrix(matrix, arithmetic, candidate_field)
            for matrix in matrices
        ]
        candidate_columns = [
            [matrix[row][column] for row in range(4)]
            for matrix in candidate_matrices
            for column in range(4)
        ]
        rank_four = first_rank_four_columns(
            candidate_field,
            candidate_columns,
        )
        if rank_four is not None:
            rank_four_indices, rank_four_determinant = rank_four
            branch_checks.append(
                {
                    "branch_index": branch_index,
                    "branch_modulus": sp.sstr(branch_modulus.as_expr()),
                    "multiplicity": multiplicity,
                    "chart_admissible": True,
                    "common_image_rank": 4,
                    "unit_ideal": True,
                    "nonzero_minor_columns_zero_based_in_M6_M7_M8": list(
                        rank_four_indices
                    ),
                    "nonzero_minor_determinant": sp.sstr(
                        rank_four_determinant.as_expr()
                    ),
                }
            )
            continue
        try:
            candidate_indices, candidate_selected = (
                first_rank_three_columns(
                    candidate_field,
                    candidate_columns,
                )
            )
        except AssertionError:
            branch_checks.append(
                {
                    "branch_index": branch_index,
                    "branch_modulus": sp.sstr(branch_modulus.as_expr()),
                    "multiplicity": multiplicity,
                    "chart_admissible": True,
                    "common_image_rank": "at_most_2",
                    "unit_ideal": False,
                }
            )
            continue
        candidate_functional = left_kernel_cofactor(
            candidate_field,
            candidate_selected,
        )
        if candidate_functional[0].is_zero:
            branch_checks.append(
                {
                    "branch_index": branch_index,
                    "branch_modulus": sp.sstr(branch_modulus.as_expr()),
                    "multiplicity": multiplicity,
                    "chart_admissible": True,
                    "common_image_rank": 3,
                    "unit_ideal": False,
                    "one_dimensional_unital_quotient": False,
                }
            )
            continue
        assert all(
            dot(candidate_field, candidate_functional, column).is_zero
            for column in candidate_columns
        )
        selected_branch = (
            branch_modulus,
            candidate_field,
            candidate_columns,
            candidate_indices,
            candidate_functional,
        )
        branch_checks.append(
            {
                "branch_index": branch_index,
                "branch_modulus": sp.sstr(branch_modulus.as_expr()),
                "multiplicity": multiplicity,
                "chart_admissible": True,
                "common_image_rank": 3,
                "unit_ideal": False,
                "one_dimensional_unital_quotient": True,
            }
        )
        break
    if selected_branch is None:
        assert len(branch_checks) == len(modulus_factorization)
        assert all(
            not check["chart_admissible"] or check["unit_ideal"]
            for check in branch_checks
        ), branch_checks
        excluded_by_localization = all(
            not check["chart_admissible"] for check in branch_checks
        )
        assert excluded_by_localization
        result = {
            "format": (
                "two-pair-sic-bidegree33-Q-cubic-exception-v1"
            ),
            "status": (
                "exact characteristic-zero exclusion: every branch above "
                "the cubic mu_6-norm factor is either removed by the "
                "u=s0^-1 localization or has unit (mu_6,mu_7,mu_8) ideal"
            ),
            "coefficient_field": {
                "theta_minpoly": sp.sstr(factor.as_expr()),
                "theta_degree": factor.degree(),
                "u_minpoly_over_Q_theta": sp.sstr(
                    specialized_modulus.as_expr()
                ),
                "u_extension_factor_degrees": modulus_factor_degrees,
            },
            "branch_checks": branch_checks,
            "conclusion": {
                "cubic_mu6_norm_factor_excluded": True,
                "degree_100_mu6_norm_factor_excluded_by_mu7": True,
                "complete_one_parameter_slice_excluded": True,
                "corrected_moment_zero_point": False,
                "excluded_entirely_by_u_localization": (
                    excluded_by_localization
                ),
                "orders_used_on_admissible_branches": (
                    [] if excluded_by_localization else [6, 7, 8]
                ),
            },
            "reproduction_command": REPRODUCTION_COMMAND,
        }
        output = arguments.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(
            json.dumps(
                {
                    "cubic_mu6_norm_factor_excluded": True,
                    "complete_one_parameter_slice_excluded": True,
                    "branches": len(branch_checks),
                    "branch_outcomes": [
                        (
                            "outside_u_open"
                            if not check["chart_admissible"]
                            else check["common_image_rank"]
                        )
                        for check in branch_checks
                    ],
                },
                indent=2,
            )
        )
        return
    (
        branch_modulus,
        field,
        columns,
        selected_indices,
        functional,
    ) = selected_branch

    inverse_constant = field.inverse(functional[0])
    functional = [
        field.mul(entry, inverse_constant) for entry in functional
    ]
    assert functional[0] == field.one
    s5_value = functional[1]
    s6_value = functional[2]
    assert functional[3] == field.mul(s5_value, s5_value)

    higher_checks = {}
    for order, artifact_path in HIGHER_ARTIFACTS.items():
        higher = json.loads(artifact_path.read_text(encoding="utf-8"))
        vector = [
            arithmetic.make(sp.sympify(entry, locals=locals_map))
            for entry in higher["normal_form_vector"]
        ]
        specialized_vector = [
            specialize_extension_element_to_number_field(
                entry,
                arithmetic,
                field,
            )
            for entry in vector
        ]
        evaluation = dot(field, functional, specialized_vector)
        assert evaluation.is_zero
        higher_checks[str(order)] = {
            "normal_form_support": higher["normal_form_support"],
            "evaluation_zero": True,
        }

    pivot = specialize_extension_element_to_number_field(
        arithmetic.make(
            sp.sympify(checkpoint["pivot_s3"], locals=locals_map)
        ),
        arithmetic,
        field,
    )
    theta_value = field.constant(number_field.unit)
    u_value = field.generator
    s4 = field.mul(
        field.add(
            field.constant(-sp.Rational(325, 9)),
            field.mul(field.constant(-sp.Rational(2, 3)), theta_value),
        ),
        field.power(u_value, 3),
    )
    t3 = field.mul(
        field.neg(
            field.add(
                field.mul(field.constant(2), pivot),
                field.constant(-3),
            )
        ),
        u_value,
    )
    t2 = field.mul(theta_value, field.power(u_value, 2))
    s2 = field.mul(field.constant(-sp.Rational(13, 3)), u_value)
    t4_numerator = field.add(
        field.mul(
            field.constant(3),
            field.mul(s6_value, field.inverse(u_value)),
        ),
        field.mul(field.constant(45), field.mul(s2, s4)),
    )
    t4_numerator = field.add(
        t4_numerator,
        field.neg(field.mul(field.constant(30), field.mul(pivot, pivot))),
    )
    t4_numerator = field.add(
        t4_numerator,
        field.neg(field.mul(field.constant(42), field.mul(t2, t2))),
    )
    t4_numerator = field.add(t4_numerator, field.constant(-70))
    t4 = field.mul(field.constant(sp.Rational(1, 14)), t4_numerator)
    s0 = field.inverse(u_value)

    coordinates = {
        "s0": sp.sstr(s0.as_expr()),
        "s1": "0",
        "s2": sp.sstr(s2.as_expr()),
        "s3": sp.sstr(pivot.as_expr()),
        "s4": sp.sstr(s4.as_expr()),
        "s5": sp.sstr(s5_value.as_expr()),
        "s6": sp.sstr(s6_value.as_expr()),
        "t0": "1",
        "t1": "0",
        "t2": sp.sstr(t2.as_expr()),
        "t3": sp.sstr(t3.as_expr()),
        "t4": sp.sstr(t4.as_expr()),
    }
    result = {
        "format": "two-pair-sic-bidegree33-counterexample-candidate-v1",
        "status": (
            "exact characteristic-zero algebraic moment-zero point; "
            "semistable because its equivariant Sym^2 projection is the "
            "fixed non-null quadratic 2*X*T"
        ),
        "coefficient_field": {
            "theta_minpoly": sp.sstr(factor.as_expr()),
            "theta_degree": factor.degree(),
            "u_minpoly_over_Q_theta": sp.sstr(
                specialized_modulus.as_expr()
            ),
            "u_extension_factor_degrees": modulus_factor_degrees,
            "selected_u_linear_factor": sp.sstr(branch_modulus.as_expr()),
            "total_degree": factor.degree(),
        },
        "fiber_quotient": {
            "basis": ["1", "s5", "s6", "s5^2"],
            "common_ideal_rank": 3,
            "common_quotient_length": 1,
            "rank_three_columns_zero_based_in_M6_M7_M8": list(
                selected_indices
            ),
            "evaluation_functional": [
                sp.sstr(entry.as_expr()) for entry in functional
            ],
            "s5_square_consistency": True,
        },
        "normalized_irreducible_coordinates": coordinates,
        "moment_checks": {
            "orders_1_2": (
                "zero by the normalized non-null-quadratic chart and the "
                "three exact chart pivots"
            ),
            "orders_3_5": (
                "zero in the exact Q-residual field and length-four "
                "(mu4,mu5) quotient"
            ),
            "orders_6_8": {
                "evaluation_zero": True,
                "common_image_rank": 3,
            },
            "higher_orders": higher_checks,
        },
        "semistability": {
            "equivariant_projection": "Sym^2",
            "projection_value": "2*X*T",
            "projection_discriminant_nonzero": True,
            "outside_one_sided_nullcone": True,
        },
        "reproduction_command": REPRODUCTION_COMMAND,
    }
    if arguments.output is not None:
        output = arguments.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "field_total_degree": result["coefficient_field"][
                    "total_degree"
                ],
                "common_quotient_length": 1,
                "higher_orders_zero": list(higher_checks),
                "semistable": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
