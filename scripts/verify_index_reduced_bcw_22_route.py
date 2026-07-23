#!/usr/bin/env python3
"""Freeze the 22-variable circuit-level witness of nilpotency index 18."""

from __future__ import annotations

import json
from pathlib import Path
import re
import shutil
import subprocess

import sympy as sp

from rank_compressed_bcw_homogenization import (
    extract_quadratic_cubic,
    factor_cubic_output,
    iterated_constant_kernel_quotient,
    rank_compressed_homogeneous_map,
    verify_parametric_factorization,
)
from search_restricted_bcw_circuits import replay_encoded_plan


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "index_reduced_bcw_22_counterexample.json"
)

# The circuit atom exposes q=x*y^2 and a=x*y, then applies the single
# multi-term target shear F_1 <- F_1-12 Q(2+3A).  The remaining entries are
# ordinary shared-factor cleanups in the post-circuit variable order.
FROZEN_PLAN: dict[str, object] = {
    "circuit_atoms": ["qb"],
    "monomial_plan": [
        {"component": 2, "first_factor": [[0, 1], [1, 2]], "second_factor": [[0, 2], [1, 1], [2, 1]]},
        {"component": 2, "first_factor": [[0, 1], [1, 2]], "second_factor": [[0, 1], [1, 2]]},
        {"component": 1, "first_factor": [[0, 1], [1, 2]], "second_factor": [[0, 2], [2, 1]]},
        {"component": 2, "first_factor": [[0, 1], [1, 1]], "second_factor": [[0, 1], [1, 1], [2, 1]]},
        {"component": 2, "first_factor": [[0, 1], [2, 1], [4, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 1, "first_factor": [[0, 1], [2, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 1, "first_factor": [[0, 1], [4, 1]], "second_factor": [[0, 1], [2, 1]]},
        {"component": 2, "first_factor": [[1, 1], [4, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 2, "first_factor": [[1, 2]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 5, "first_factor": [[0, 1], [2, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 1, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 2, "first_factor": [[1, 1], [3, 1]], "second_factor": [[0, 1], [2, 1]]},
        {"component": 2, "first_factor": [[2, 1], [3, 1]], "second_factor": [[0, 1], [4, 1]]},
        {"component": 1, "first_factor": [[1, 1], [6, 1]], "second_factor": [[0, 1], [1, 1]]},
        {"component": 0, "first_factor": [[0, 1], [2, 1]], "second_factor": [[0, 2]]},
        {"component": 2, "first_factor": [[1, 1], [5, 1]], "second_factor": [[0, 1], [1, 1]]},
    ],
}


def qtext(value: sp.Expr) -> str:
    value = sp.cancel(value)
    assert value.is_Rational
    return str(value)


def encode_h(components: tuple[sp.Poly, ...]) -> list[list[dict[str, object]]]:
    return [
        [
            {
                "coefficient": qtext(coefficient),
                "monomial": [
                    [index, exponent]
                    for index, exponent in enumerate(exponents)
                    if exponent
                ],
            }
            for exponents, coefficient in poly.terms()
            if coefficient
        ]
        for poly in components
    ]


def sparse_matrix_rows(matrix: sp.Matrix) -> list[list[list[object]]]:
    return [
        [[column, qtext(matrix[row, column])] for column in range(matrix.cols) if matrix[row, column]]
        for row in range(matrix.rows)
    ]


def singular_entry(poly: sp.Poly, variable: int) -> str:
    terms: list[str] = []
    for exponents, coefficient in poly.terms():
        exponent = exponents[variable]
        if not exponent:
            continue
        factors = [f"({qtext(coefficient)})", str(exponent)]
        for index, power in enumerate(exponents):
            power -= index == variable
            if power:
                factors.append(f"x{index}^{power}")
        terms.append("*".join(factors))
    return "+".join(terms).replace("+-", "-") or "0"


def exact_rank_and_power_certificate(
    components: tuple[sp.Poly, ...],
    *,
    expected_rank: int = 18,
    expected_kernel_rank: int = 4,
) -> dict[str, int]:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required for the exact rank/index certificate")
    dimension = len(components)
    entries = [
        singular_entry(poly, variable)
        for poly in components
        for variable in range(dimension)
    ]
    commands = [
        "ring r=0,(" + ",".join(f"x{i}" for i in range(dimension)) + "),dp;",
        f"matrix J[{dimension}][{dimension}]=" + ",".join(entries) + ";",
        "option(redSB);",
        "module K=syz(J);",
        "matrix Z=J*matrix(K);",
        "int product_zero=1;",
        "int i,j,k,nz17,nz18;",
        (
            f"for(i=1;i<={dimension};i++)"
            "{for(j=1;j<=ncols(K);j++){if(Z[i,j]!=0){product_zero=0;}}}"
        ),
        "matrix Js=J;",
        "module Ks=K;",
    ]
    point = [index * index + 3 * index + 5 for index in range(dimension)]
    for index, value in enumerate(point):
        commands.extend(
            [f"Js=subst(Js,x{index},{value});", f"Ks=subst(Ks,x{index},{value});"]
        )
    commands.extend(
        [
            "matrix P=J;",
            (
                "for(k=2;k<=18;k++){P=P*J;"
                f"if(k==17){{for(i=1;i<={dimension};i++)"
                f"{{for(j=1;j<={dimension};j++){{if(P[i,j]!=0){{nz17++;}}}}}}}}}}"
            ),
            (
                f"for(i=1;i<={dimension};i++)"
                f"{{for(j=1;j<={dimension};j++){{if(P[i,j]!=0){{nz18++;}}}}}}"
            ),
            'print("CERTIFICATE");',
            "ncols(K);",
            "product_zero;",
            "rank(Js);",
            "rank(Ks);",
            "nz17;",
            "nz18;",
            "quit;",
        ]
    )
    result = subprocess.run(
        [singular, "-q"],
        input="\n".join(commands) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if "CERTIFICATE" not in lines:
        raise RuntimeError(
            "Singular certificate did not complete:\n"
            + result.stdout
            + "\n"
            + result.stderr
        )
    marker = lines.index("CERTIFICATE")
    values = list(map(int, lines[marker + 1 : marker + 7]))
    syzygy_columns, product_zero, specialized_rank, syzygy_rank, nz17, nz18 = values
    assert syzygy_columns >= expected_kernel_rank
    assert product_zero == 1
    assert specialized_rank == expected_rank
    assert syzygy_rank == expected_kernel_rank
    assert nz17 > 0 and nz18 == 0
    return {
        "syzygy_generators": syzygy_columns,
        "independent_generic_kernel_columns": syzygy_rank,
        "specialized_rank_lower_bound": specialized_rank,
        "nonzero_entries_JH_power_17": nz17,
        "nonzero_entries_JH_power_18": nz18,
    }


def main() -> None:
    state = replay_encoded_plan(FROZEN_PLAN)
    assert len(state.variables) == 18
    assert state.introduced == 15
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
    assert len(ambient_variables) == 27
    quotient = iterated_constant_kernel_quotient(ambient_variables, ambient_h)
    assert [stage.kernel.cols for stage in quotient.stages] == [4, 1]
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

    certificate = exact_rank_and_power_certificate(quotient.quotient_h)
    # Four independent generic kernel columns give rank <= 18; the exact
    # specialization gives rank >= 18.  The polynomial power calculation is
    # a direct nilpotence/Keller certificate, independent of the BCW bridge.
    assert certificate["independent_generic_kernel_columns"] == 4
    assert certificate["specialized_rank_lower_bound"] == 18
    assert certificate["nonzero_entries_JH_power_17"] > 0
    assert certificate["nonzero_entries_JH_power_18"] == 0

    artifact = {
        "format": "index-reduced-bcw-sparse-cubic-homogeneous-map-v1",
        "source": (
            "frozen polynomial-gate BCW circuit, shared-factor cleanup, "
            "rank-compressed homogenization, and constant-kernel quotient"
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
        "statistics": {
            "polynomial_circuit_atoms": 1,
            "polynomial_circuit_gates": 2,
            "monomial_cleanup_steps": len(FROZEN_PLAN["monomial_plan"]),
            "introduced_variables": state.introduced,
            "cubic_output_rank": len(factorization.c),
            "ambient_constant_kernel_dimensions": [
                stage.kernel.cols for stage in quotient.stages
            ],
            "quotient_constant_kernel_dimension": 0,
            "generic_rank_JH_over_QQ_x": 18,
            "nilpotency_index_JH": 18,
            "specialized_jordan_type": [18, 2, 1, 1],
            "exact_certificate": certificate,
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS index-reduced BCW 22: polynomial circuit plus 16 cleanups reaches dimension 18")
    print("PASS index-reduced BCW 22: rank 8 homogenization has dimension 27")
    print("PASS index-reduced BCW 22: successive constant kernels have dimensions 4 and 1")
    print("PASS index-reduced BCW 22: quotient collision remains separated in dimension 22")
    print("PASS index-reduced BCW 22: generic rank(JH) is exactly 18")
    print("PASS index-reduced BCW 22: (JH)^17 is nonzero and (JH)^18 is identically zero")
    print(f"PASS index-reduced BCW 22: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
