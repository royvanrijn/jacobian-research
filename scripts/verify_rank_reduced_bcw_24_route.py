#!/usr/bin/env python3
"""Freeze the circuit-level 24D witness of generic Jacobian rank 17."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from rank_compressed_bcw_homogenization import (
    extract_quadratic_cubic,
    factor_cubic_output,
    iterated_constant_kernel_quotient,
    rank_compressed_homogeneous_map,
    verify_parametric_factorization,
)
from search_restricted_bcw_circuits import replay_encoded_plan
from verify_index_reduced_bcw_22_route import (
    encode_h,
    exact_rank_and_power_certificate,
    qtext,
    sparse_matrix_rows,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "rank_reduced_bcw_24_counterexample.json"
)


FROZEN_PLAN: dict[str, object] = {
    "circuit_atoms": ["v2r", "y2vb"],
    "monomial_plan": [
        {"component": 2, "first_factor": [[0, 1], [1, 1], [2, 1]], "second_factor": [[0, 2], [1, 2]]},
        {"component": 1, "first_factor": [[0, 1], [2, 1]], "second_factor": [[0, 1], [1, 1], [3, 1]]},
        {"component": 2, "first_factor": [[1, 2]], "second_factor": [[0, 1], [1, 1], [3, 1]]},
        {"component": 2, "first_factor": [[0, 1], [1, 1]], "second_factor": [[0, 1], [1, 1], [2, 1]]},
        {"component": 1, "first_factor": [[0, 1], [1, 1]], "second_factor": [[0, 1], [1, 1], [4, 1]]},
        {"component": 1, "first_factor": [[0, 1], [1, 1]], "second_factor": [[0, 1], [1, 2]]},
        {"component": 2, "first_factor": [[0, 1], [1, 1]], "second_factor": [[0, 1], [1, 1], [5, 1]]},
        {"component": 2, "first_factor": [[0, 1], [1, 1]], "second_factor": [[0, 1], [1, 1], [6, 1]]},
        {"component": 7, "first_factor": [[0, 1], [1, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 0, "first_factor": [[0, 1], [2, 1]], "second_factor": [[0, 2]]},
        {"component": 1, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 2, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [2, 1]]},
        {"component": 2, "first_factor": [[1, 1], [3, 1]], "second_factor": [[1, 1], [3, 1]]},
        {"component": 1, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [4, 1]]},
        {"component": 1, "first_factor": [[3, 2]], "second_factor": [[0, 1], [2, 1]]},
        {"component": 2, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [5, 1]]},
        {"component": 2, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [6, 1]]},
        {"component": 2, "first_factor": [[2, 1], [7, 1]], "second_factor": [[0, 1], [1, 1]]},
    ],
}


def main() -> None:
    state = replay_encoded_plan(FROZEN_PLAN)
    assert len(state.variables) == 20 and state.introduced == 17
    assert max(
        sp.Poly(expression, *state.variables, domain=sp.QQ).total_degree()
        for expression in state.expressions
    ) == 3

    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    assert len(factorization.c) == 9
    verify_parametric_factorization(
        state.variables, quadratic, cubic, factorization
    )
    ambient_variables, ambient_h = rank_compressed_homogeneous_map(
        state.variables, quadratic, factorization
    )
    assert len(ambient_variables) == 30
    quotient = iterated_constant_kernel_quotient(ambient_variables, ambient_h)
    assert [stage.kernel.cols for stage in quotient.stages] == [4, 2]
    assert len(quotient.quotient_variables) == 24

    quotient_points = []
    for point in state.collision_points:
        substitution = dict(zip(state.variables, point))
        quotient_points.append(
            sp.Matrix(
                list(point)
                + [poly.eval(substitution) for poly in factorization.c]
                + [sp.Integer(1)]
            )
        )
    for stage in quotient.stages:
        quotient_points = [stage.B * point for point in quotient_points]
    assert len({tuple(point) for point in quotient_points}) == 3
    images = []
    for point in quotient_points:
        substitution = dict(zip(quotient.quotient_variables, point))
        images.append(
            point
            + sp.Matrix(
                [
                    poly.as_expr().subs(substitution, simultaneous=True)
                    for poly in quotient.quotient_h
                ]
            )
        )
    assert images[0] == images[1] == images[2]

    certificate = exact_rank_and_power_certificate(
        quotient.quotient_h,
        expected_rank=17,
        expected_kernel_rank=7,
    )
    assert certificate["independent_generic_kernel_columns"] == 7
    assert certificate["specialized_rank_lower_bound"] == 17
    assert certificate["nonzero_entries_JH_power_17"] > 0
    assert certificate["nonzero_entries_JH_power_18"] == 0

    artifact = {
        "format": "rank-reduced-bcw-sparse-cubic-homogeneous-map-v1",
        "source": (
            "frozen two-atom polynomial-gate BCW circuit, shared-factor "
            "cleanup, rank-compressed homogenization, and fixed-point quotient"
        ),
        "dimension": 24,
        "linear_part": "identity",
        "jacobian_determinant": "1",
        "jacobian_certificate": (
            "direct polynomial identity (JH)^18=0 for the cubic-homogeneous "
            "correction, hence det(I+JH)=1"
        ),
        "circuit_plan": FROZEN_PLAN,
        "degree_reduced_dimension": len(state.variables),
        "rank_factorization": {
            "basis_components": list(factorization.basis_components),
            "B": [
                [qtext(value) for value in factorization.B.row(row)]
                for row in range(factorization.B.rows)
            ],
        },
        "ambient_homogeneous_dimension": len(ambient_variables),
        "quotient_factorization": {
            "kernel_dimensions": [
                stage.kernel.cols for stage in quotient.stages
            ],
            "stages": [
                {
                    "ambient_dimension": len(stage.ambient_variables),
                    "quotient_dimension": len(stage.quotient_variables),
                    "kernel_columns": sparse_matrix_rows(stage.kernel.T),
                    "B_rows": sparse_matrix_rows(stage.B),
                    "C_rows": sparse_matrix_rows(stage.C),
                }
                for stage in quotient.stages
            ],
        },
        "H": encode_h(quotient.quotient_h),
        "collision_points": [
            [qtext(value) for value in point] for point in quotient_points
        ],
        "common_image": [qtext(value) for value in images[0]],
        "statistics": {
            "polynomial_circuit_atoms": 2,
            "polynomial_circuit_gates": 3,
            "monomial_cleanup_steps": len(FROZEN_PLAN["monomial_plan"]),
            "introduced_variables": state.introduced,
            "cubic_output_rank": len(factorization.c),
            "ambient_constant_kernel_dimensions": [
                stage.kernel.cols for stage in quotient.stages
            ],
            "quotient_constant_kernel_dimension": 0,
            "generic_rank_JH_over_QQ_x": 17,
            "nilpotency_index_JH": 18,
            "specialized_jordan_type": [18, 1, 1, 1, 1, 1, 1],
            "exact_certificate": certificate,
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS rank-reduced BCW 24: two circuit atoms plus 18 cleanups reach dimension 20")
    print("PASS rank-reduced BCW 24: rank 9 homogenization has dimension 30")
    print("PASS rank-reduced BCW 24: successive constant kernels have dimensions 4 and 2")
    print("PASS rank-reduced BCW 24: quotient collision remains separated in dimension 24")
    print("PASS rank-reduced BCW 24: generic rank(JH) is exactly 17")
    print("PASS rank-reduced BCW 24: (JH)^17 is nonzero and (JH)^18 is identically zero")
    print(f"PASS rank-reduced BCW 24: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
