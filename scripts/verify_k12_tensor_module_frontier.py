#!/usr/bin/env python3
"""Exact tensor-flattening and constant-module audit for K12 and G19."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

from audit_k12_coordinate_pair_frontier import build_k12
from verify_macfarlane_f12_reduction import homogeneous_part


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "k12_tensor_module_frontier.json"
)


def coefficient_matrix(
    components: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    polynomials = [sp.Poly(component, *variables, domain=sp.QQ) for component in components]
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
    return sp.Matrix(
        [
            [polynomial.coeff_monomial(monomial) for polynomial in polynomials]
            for monomial in monomials
        ]
    )


def directional_flattening(
    components: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> sp.Matrix:
    derivatives = [
        [
            sp.Poly(sp.diff(component, variable), *variables, domain=sp.QQ)
            for variable in variables
        ]
        for component in components
    ]
    monomials = sorted(
        set().union(
            *[
                {
                    exponents
                    for row in derivatives
                    for polynomial in row
                    for exponents, coefficient in polynomial.terms()
                    if coefficient
                }
            ]
        )
    )
    rows = []
    for output_row in derivatives:
        for monomial in monomials:
            rows.append(
                [
                    polynomial.coeff_monomial(monomial)
                    for polynomial in output_row
                ]
            )
    return sp.Matrix(rows)


def jacobian_coefficient_span_rank(
    components: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> int:
    jacobian = sp.Matrix(components).jacobian(variables)
    polynomials = [
        sp.Poly(jacobian[row, column], *variables, domain=sp.QQ)
        for row in range(jacobian.rows)
        for column in range(jacobian.cols)
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
    flattened_coefficient_matrices = sp.Matrix(
        [
            [polynomial.coeff_monomial(monomial) for polynomial in polynomials]
            for monomial in monomials
        ]
    )
    return flattened_coefficient_matrices.rank()


def encode_vector(vector: sp.Matrix) -> list[str]:
    return [str(sp.cancel(value)) for value in vector]


def tensor_record(
    name: str,
    components: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    homogeneous_degree: int,
) -> dict[str, object]:
    assert all(
        not component
        or sp.Poly(component, *variables, domain=sp.QQ).homogeneous_order()
        == homogeneous_degree
        for component in components
    )
    output_flattening = coefficient_matrix(components, variables)
    input_flattening = directional_flattening(components, variables)
    output_rank = output_flattening.rank()
    input_rank = input_flattening.rank()
    left_kernel = output_flattening.nullspace()
    right_kernel = input_flattening.nullspace()
    assert len(left_kernel) == len(components) - output_rank
    assert len(right_kernel) == len(variables) - input_rank
    return {
        "name": name,
        "input_dimension": len(variables),
        "output_dimension": len(components),
        "homogeneous_degree": homogeneous_degree,
        "nonzero_output_count": sum(bool(component) for component in components),
        "total_monomial_count": sum(
            len(sp.Poly(component, *variables, domain=sp.QQ).terms())
            for component in components
        ),
        "output_flattening_rank": output_rank,
        "input_directional_flattening_rank": input_rank,
        "common_left_kernel_dimension_of_jacobian_coefficients": len(left_kernel),
        "common_left_kernel_basis": [encode_vector(vector) for vector in left_kernel],
        "common_right_kernel_dimension_of_jacobian_coefficients": len(right_kernel),
        "common_right_kernel_basis": [encode_vector(vector) for vector in right_kernel],
        "common_column_module_dimension": output_rank,
        "common_row_module_dimension": input_rank,
        "jacobian_coefficient_matrix_span_dimension": jacobian_coefficient_span_rank(
            components, variables
        ),
        "partially_symmetric_cp_rank_lower_bound": max(output_rank, input_rank),
    }


def evaluate(components: list[sp.Expr], variables: tuple[sp.Symbol, ...], point: tuple[sp.Expr, ...]) -> tuple[sp.Expr, ...]:
    substitution = dict(zip(variables, point))
    return tuple(sp.cancel(component.subs(substitution)) for component in components)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    z, k12 = build_k12()
    nonlinear = [
        sp.expand(component - variable)
        for component, variable in zip(k12, z)
    ]
    quadratic = [homogeneous_part(component, z, 2) for component in nonlinear]
    cubic = [homogeneous_part(component, z, 3) for component in nonlinear]
    assert nonlinear == [sp.expand(q + c) for q, c in zip(quadratic, cubic)]

    k12_quadratic = tensor_record("K12 quadratic correction Q", quadratic, z, 2)
    k12_cubic = tensor_record("K12 cubic correction C", cubic, z, 3)
    assert k12_cubic["output_flattening_rank"] == 6
    assert k12_cubic["input_directional_flattening_rank"] == 12
    assert k12_cubic["common_right_kernel_dimension_of_jacobian_coefficients"] == 0

    w = tuple(sp.symbols("w1:7"))
    tau = sp.Symbol("tau")
    variables19 = z + w + (tau,)
    h19 = [
        sp.expand(
            tau * quadratic[index]
            + (tau**2 * w[index] if index < 6 else 0)
        )
        for index in range(12)
    ] + [-component for component in cubic[:6]] + [sp.Integer(0)]
    g19_record = tensor_record("G19 cubic homogeneous correction H", h19, variables19, 3)

    # The full cubic tensor uses every input direction.  Its only constant
    # left annihilator is the fixed tau output.
    assert g19_record["input_directional_flattening_rank"] == 19
    assert g19_record["output_flattening_rank"] == 18
    assert g19_record["common_right_kernel_dimension_of_jacobian_coefficients"] == 0
    assert g19_record["common_left_kernel_dimension_of_jacobian_coefficients"] == 1
    assert g19_record["common_left_kernel_basis"] == [
        ["0"] * 18 + ["1"]
    ]

    p12 = (
        sp.Integer(0), sp.Integer(0), -sp.Rational(1, 4),
        *([sp.Integer(0)] * 9),
    )
    q12 = (
        sp.Integer(1), -sp.Rational(3, 2), sp.Rational(13, 2),
        -sp.Rational(9, 4), sp.Integer(3), sp.Rational(3, 2),
        sp.Rational(99, 4), sp.Rational(3, 2), -sp.Rational(3, 4),
        -sp.Rational(45, 8), -sp.Rational(13, 2), sp.Rational(1, 2),
    )
    assert evaluate(k12, z, p12) == evaluate(k12, z, q12) == p12
    cubic_at_p = evaluate(cubic[:6], z, p12)
    cubic_at_q = evaluate(cubic[:6], z, q12)
    p19 = p12 + cubic_at_p + (sp.Integer(1),)
    q19 = q12 + cubic_at_q + (sp.Integer(1),)
    g19 = [
        sp.expand(variable + correction)
        for variable, correction in zip(variables19, h19)
    ]
    image19 = evaluate(g19, variables19, p19)
    assert p19 != q19 and evaluate(g19, variables19, q19) == image19

    # Recheck the companion/scaling determinant bridge used by BCR2.  This
    # audit deliberately does not expand the 19-by-19 determinant.
    e_tau = [
        sp.expand(z[index] + tau * quadratic[index] + tau**2 * cubic[index])
        for index in range(12)
    ]
    scaled_k12 = [
        sp.cancel(
            component.subs(
                {variable: tau * variable for variable in z},
                simultaneous=True,
            ) / tau
        )
        for component in k12
    ]
    assert e_tau == scaled_k12

    artifact = {
        "format": "k12-tensor-module-frontier-v1",
        "status": "exact rational linear-algebra theorem",
        "K12": {
            "quadratic_tensor": k12_quadratic,
            "cubic_tensor": k12_cubic,
            "collision_replayed": True,
        },
        "G19": {
            "cubic_tensor": g19_record,
            "collision_replayed": True,
            "determinant_certificate": (
                "exact companion cancellation plus E_tau(z)=tau^-1*K12(tau*z); "
                "the direct determinant-one K12 certificate is replayed by "
                "make verify-macfarlane-f12"
            ),
        },
        "consequences": {
            "K12_cubic_tensor_has_no_constant_source_direction": True,
            "G19_tensor_has_no_constant_source_direction": True,
            "G19_only_constant_output_annihilator": "the tau output",
            "linear_tensor_quotient_of_G19_is_obstructed": True,
            "pure_cube_summand_lower_bounds": {
                "K12_cubic_tensor": k12_cubic[
                    "partially_symmetric_cp_rank_lower_bound"
                ],
                "G19_cubic_tensor": g19_record[
                    "partially_symmetric_cp_rank_lower_bound"
                ],
            },
        },
        "scope": (
            "These flattening and common-kernel calculations exclude only constant "
            "linear source quotients and give lower bounds for decompositions "
            "H=sum u_s*ell_s^3. They do not exclude nonlinear graph coordinates, "
            "nonlinear row/column modules, Schur elimination, or a different tensor."
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
    print("PASS K12: cubic output rank 6 and full input-directional rank 12")
    print("PASS G19: output rank 18 and full input-directional rank 19")
    print("PASS G19: no constant right kernel; sole left kernel is tau output")
    print("PASS replayed both collisions and the exact determinant bridge")
    print(f"PASS checked {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
