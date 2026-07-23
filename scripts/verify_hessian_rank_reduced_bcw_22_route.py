#!/usr/bin/env python3
"""Freeze the 22D cubic source whose quartic HN lift has Hessian rank 37."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from audit_fixed_rank_hessian_witness import (
    cotangent_hessian,
    deterministic_point,
    exact_singular_certificate,
    nilpotency_profile,
    specialization_profile,
)
from rank_compressed_bcw_homogenization import (
    extract_quadratic_cubic,
    factor_cubic_output,
    iterated_constant_kernel_quotient,
    rank_compressed_homogeneous_map,
    verify_parametric_factorization,
)
from restricted_rank_profiles import _fraction_polynomial
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
    / "hessian_rank_reduced_bcw_22_counterexample.json"
)


FROZEN_PLAN: dict[str, object] = {
    "circuit_atoms": ["qb", "x2s"],
    "monomial_plan": [
        {"component": 2, "first_factor": [[0, 1], [1, 2]], "second_factor": [[0, 2], [1, 1], [2, 1]]},
        {"component": 2, "first_factor": [[0, 1], [1, 2]], "second_factor": [[0, 1], [1, 2]]},
        {"component": 1, "first_factor": [[0, 1], [1, 2]], "second_factor": [[0, 2], [2, 1]]},
        {"component": 6, "first_factor": [[0, 1], [1, 2]], "second_factor": [[0, 2], [2, 1]]},
        {"component": 2, "first_factor": [[0, 1], [1, 1]], "second_factor": [[0, 1], [1, 1], [2, 1]]},
        {"component": 2, "first_factor": [[0, 1], [2, 1], [4, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 1, "first_factor": [[0, 1], [2, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 6, "first_factor": [[0, 1], [2, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 7, "first_factor": [[0, 1], [2, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 1, "first_factor": [[0, 1], [4, 1]], "second_factor": [[0, 1], [2, 1]]},
        {"component": 6, "first_factor": [[0, 1], [4, 1]], "second_factor": [[0, 1], [2, 1]]},
        {"component": 1, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 2, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [2, 1]]},
        {"component": 6, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 2, "first_factor": [[1, 1], [4, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 1, "first_factor": [[1, 1], [8, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 6, "first_factor": [[1, 1], [8, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 2, "first_factor": [[1, 1], [7, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 2, "first_factor": [[1, 2]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 2, "first_factor": [[3, 1], [4, 1]], "second_factor": [[0, 1], [2, 1]]},
    ],
}


def main() -> None:
    state = replay_encoded_plan(FROZEN_PLAN)
    assert len(state.variables) == 19
    assert state.introduced == 16
    assert len(state.monomial_plan) == 20
    assert max(
        sp.Poly(expression, *state.variables, domain=sp.QQ).total_degree()
        for expression in state.expressions
    ) == 3

    quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
    factorization = factor_cubic_output(cubic)
    assert len(factorization.c) == 8
    verify_parametric_factorization(
        state.variables, quadratic, cubic, factorization
    )
    ambient_variables, ambient_h = rank_compressed_homogeneous_map(
        state.variables, quadratic, factorization
    )
    assert len(ambient_variables) == 28
    quotient = iterated_constant_kernel_quotient(ambient_variables, ambient_h)
    assert [stage.kernel.cols for stage in quotient.stages] == [5, 1]
    assert len(quotient.quotient_variables) == 22

    quotient_points = []
    for point in state.collision_points:
        substitution = dict(zip(state.variables, point))
        ambient_point = sp.Matrix(
            list(point)
            + [poly.eval(substitution) for poly in factorization.c]
            + [sp.Integer(1)]
        )
        quotient_points.append(ambient_point)
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

    cubic_certificate = exact_rank_and_power_certificate(quotient.quotient_h)
    assert cubic_certificate["specialized_rank_lower_bound"] == 18
    assert cubic_certificate["independent_generic_kernel_columns"] == 4
    assert cubic_certificate["nonzero_entries_JH_power_17"] == 7
    assert cubic_certificate["nonzero_entries_JH_power_18"] == 0

    sparse_h = [_fraction_polynomial(poly) for poly in quotient.quotient_h]
    jacobian, upper_left, hessian = cotangent_hessian(sparse_h)
    profiles = [
        specialization_profile(jacobian, upper_left, hessian, 20_260_723 + offset)
        for offset in range(3)
    ]
    assert profiles == [(18, 16, 37)] * 3
    exact_point = deterministic_point(44, 1_000_003, 20_260_723)
    syzygies, specialized_rank, kernel_check = exact_singular_certificate(
        hessian,
        exact_point,
        expected_kernel_rank=7,
    )
    assert syzygies == 12
    assert specialized_rank == 37
    assert kernel_check
    transformed_power_ranks = nilpotency_profile(
        jacobian, upper_left, 20_260_723
    )
    assert transformed_power_ranks == [
        37,
        34,
        *range(32, -1, -1),
    ]

    artifact = {
        "format": "hessian-rank-reduced-bcw-sparse-cubic-homogeneous-map-v1",
        "source": (
            "frozen two-atom polynomial-gate BCW circuit, shared-factor "
            "cleanup, rank-compressed homogenization, and iterated "
            "constant-kernel quotient"
        ),
        "dimension": 22,
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
        "homogeneous_quartic_HN_lift": {
            "field": "QQ(I)",
            "dimension": 44,
            "variables": "u_0,...,u_21,v_0,...,v_21",
            "definition": (
                "P(u,v)=(1/4)*(u-I*v).H(u+I*v), with H encoded above"
            ),
            "degree": 4,
            "gradient_map": "(u,v)+grad(P) is Keller and noninjective",
            "hessian_nilpotency_certificate": (
                "the symmetric cotangent reduction applied to the exact "
                "nilpotent-Jacobian cubic source"
            ),
            "vanishing_counterexample_certificate": (
                "the standard HN/Vanishing equivalence applied to the "
                "three-point noninjective gradient map"
            ),
            "generic_hessian_rank": 37,
        },
        "statistics": {
            "polynomial_circuit_atoms": 2,
            "polynomial_circuit_gates": 4,
            "monomial_cleanup_steps": len(FROZEN_PLAN["monomial_plan"]),
            "introduced_variables": state.introduced,
            "cubic_output_rank": len(factorization.c),
            "ambient_constant_kernel_dimensions": [
                stage.kernel.cols for stage in quotient.stages
            ],
            "quotient_constant_kernel_dimension": 0,
            "generic_rank_JH_over_QQ_x": 18,
            "nilpotency_index_JH": 18,
            "cubic_exact_certificate": cubic_certificate,
            "cotangent_quartic_dimension": 44,
            "generic_cotangent_hessian_rank_over_QQ_xy": 37,
            "cotangent_hessian_exact_certificate": {
                "syzygy_generators": syzygies,
                "independent_generic_kernel_columns": 7,
                "specialized_rank_lower_bound": specialized_rank,
                "sampled_block_ranks": [18, 16, 37],
                "sample_count": len(profiles),
            },
            "transformed_hessian_sampled_power_ranks_mod_1000033": (
                transformed_power_ranks
            ),
            "transformed_hessian_sampled_nilpotency_index": 35,
            "transformed_hessian_sampled_index_status": (
                "diagnostic specialization only; it is not an exact "
                "polynomial-matrix nilpotency certificate"
            ),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS Hessian-rank BCW 22: exact replay reaches a 22D cubic collision")
    print("PASS Hessian-rank BCW 22: generic rank(JH)=18 and index(JH)=18")
    print("PASS Hessian-rank BCW 22: 12 syzygies contain 7 generically independent columns")
    print("PASS Hessian-rank BCW 22: exact specialized cotangent-Hessian rank is 37")
    print("PASS Hessian-rank BCW 22: transformed sampled Hessian index is 35")
    print("PASS Hessian-rank BCW 22: homogeneous quartic HN rank upper bound is 37")
    print(f"PASS Hessian-rank BCW 22: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
