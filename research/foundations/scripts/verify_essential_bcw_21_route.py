#!/usr/bin/env python3
"""Freeze and verify the 17 -> 24 -> 21 essential BCW route."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

from rank_compressed_bcw_homogenization import (
    constant_kernel_quotient,
    extract_quadratic_cubic,
    factor_cubic_output,
    rank_compressed_homogeneous_map,
    verify_parametric_factorization,
)
from search_essential_bcw import replay_collision_points
from search_rank_aware_bcw import State, initial_state, support
from verify_shared_bcw_33_route import apply_shared_step, dense_factor


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"

PLAN = [
    (2, [(0, 1), (1, 2)], [(0, 2), (1, 1), (2, 1)]),
    (2, [(0, 1), (1, 2)], [(0, 1), (1, 2)]),
    (1, [(0, 1), (1, 2)], [(0, 2), (2, 1)]),
    (2, [(0, 1), (2, 1)], [(0, 1), (1, 2)]),
    (2, [(1, 1), (3, 1)], [(0, 2), (2, 1)]),
    (1, [(0, 1), (1, 1)], [(0, 1), (1, 2)]),
    (0, [(0, 1), (2, 1)], [(0, 2)]),
    (1, [(0, 1), (2, 1)], [(0, 1), (1, 1)]),
    (2, [(1, 1), (3, 1)], [(0, 1), (1, 1)]),
    (4, [(0, 1), (2, 1)], [(0, 1), (1, 1)]),
    (1, [(0, 1), (3, 1)], [(0, 1), (2, 1)]),
    (1, [(1, 1), (5, 1)], [(0, 1), (1, 1)]),
    (1, [(1, 1), (8, 1)], [(0, 1), (1, 1)]),
    (2, [(0, 1), (7, 1)], [(0, 1), (2, 1)]),
    (2, [(1, 1), (4, 1)], [(0, 1), (1, 1)]),
    (2, [(1, 1), (6, 1)], [(0, 1), (1, 1)]),
    (2, [(1, 2)], [(0, 1), (1, 1)]),
]


def qtext(value: sp.Expr) -> str:
    value = sp.cancel(value)
    assert value.is_Rational
    return str(value)


def sparse_matrix_rows(matrix: sp.Matrix) -> list[list[list[object]]]:
    return [
        [[j, qtext(matrix[i, j])] for j in range(matrix.cols) if matrix[i, j]]
        for i in range(matrix.rows)
    ]


def encode_h(
    components: tuple[sp.Poly, ...], variables: tuple[sp.Symbol, ...]
) -> list[list[dict[str, object]]]:
    return [
        [
            {
                "coefficient": qtext(coefficient),
                "monomial": [[i, exponent] for i, exponent in enumerate(exponents) if exponent],
            }
            for exponents, coefficient in poly.terms()
            if coefficient
        ]
        for poly in components
    ]


state = initial_state()
steps = []
for component, first_support, second_support in PLAN:
    dimension = len(state.variables)
    first = dense_factor(first_support, dimension)
    second = dense_factor(second_support, dimension)
    removed = tuple(a + b for a, b in zip(first, second))
    polynomial = sp.Poly(state.expressions[component], *state.variables, domain=sp.QQ)
    coefficient = polynomial.coeff_monomial(removed)
    selected = (component, removed, coefficient, sum(removed))
    result = apply_shared_step(
        state.expressions,
        state.variables,
        state.registry,
        selected,
        (first, second),
        state.introduced,
    )
    assert result is not None
    steps.append(result[4])
    state = State(
        result[0], result[1], result[2], result[3],
        state.plan + ((component, support(first), support(second)),),
    )

assert len(state.variables) == 17 and state.introduced == 14
quadratic, cubic = extract_quadratic_cubic(state.expressions, state.variables)
factorization = factor_cubic_output(cubic)
assert len(factorization.c) == 6
verify_parametric_factorization(state.variables, quadratic, cubic, factorization)
ambient_variables, ambient_h = rank_compressed_homogeneous_map(
    state.variables, quadratic, factorization
)
assert len(ambient_variables) == 24
quotient = constant_kernel_quotient(ambient_variables, ambient_h)
assert quotient.kernel.cols == 3 and len(quotient.quotient_variables) == 21
assert constant_kernel_quotient(
    quotient.quotient_variables, quotient.quotient_h
).kernel.cols == 0

source_points = replay_collision_points(state)
ambient_points = []
for point in source_points:
    substitution = dict(zip(state.variables, point))
    ambient_points.append(
        sp.Matrix(point + [poly.eval(substitution) for poly in factorization.c] + [sp.Integer(1)])
    )
quotient_points = [quotient.B * point for point in ambient_points]
assert len({tuple(point) for point in quotient_points}) == 3
images = []
for point in quotient_points:
    substitution = dict(zip(quotient.quotient_variables, point))
    images.append(
        point + sp.Matrix([poly.as_expr().subs(substitution) for poly in quotient.quotient_h])
    )
assert images[0] == images[1] == images[2]

artifact = {
    "format": "essential-constant-kernel-bcw-sparse-cubic-homogeneous-map-v1",
    "source": "frozen 17-step shared-factor BCW trace, rank compression, and constant-kernel quotient",
    "dimension": 21,
    "linear_part": "identity",
    "jacobian_determinant": "1",
    "jacobian_certificate": (
        "17 determinant-one stable factor cancellations; rank-compressed BCW determinant bridge; "
        "three-dimensional constant-kernel triangular quotient"
    ),
    "degree_reduced_dimension": 17,
    "degree_lowering_steps": steps,
    "rank_factorization": {
        "basis_components": list(factorization.basis_components),
        "B": [[qtext(value) for value in factorization.B.row(i)] for i in range(factorization.B.rows)],
    },
    "ambient_homogeneous_dimension": 24,
    "quotient_factorization": {
        "kernel_dimension": quotient.kernel.cols,
        "kernel_columns": sparse_matrix_rows(quotient.kernel.T),
        "B_rows": sparse_matrix_rows(quotient.B),
        "C_rows": sparse_matrix_rows(quotient.C),
        "identities": ["BK=0", "BC=I_21", "H=HCB"],
    },
    "H": encode_h(quotient.quotient_h, quotient.quotient_variables),
    "collision_points": [[qtext(value) for value in point] for point in quotient_points],
    "common_image": [qtext(value) for value in images[0]],
    "statistics": {
        "target_cancellations": len(steps),
        "introduced_variables": state.introduced,
        "cubic_output_rank": len(factorization.c),
        "ambient_constant_kernel_dimension": quotient.kernel.cols,
        "quotient_constant_kernel_dimension": 0,
        "fixed_dimensional_GMC_target": 42,
    },
}
OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

print("PASS essential BCW 21: frozen 17-step trace reaches dimension 17")
print("PASS essential BCW 21: cubic-output rank 6 homogenizes in dimension 24")
print("PASS essential BCW 21: constant-kernel quotient has dimension 21")
print("PASS essential BCW 21: quotient has zero constant kernel")
print("PASS essential BCW 21: three-point collision descends and remains separated")
print(f"PASS essential BCW 21: wrote {OUTPUT.relative_to(ROOT)}")
