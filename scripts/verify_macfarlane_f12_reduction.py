#!/usr/bin/env python3
"""Exact source-graph/target-square reduction of MacFarlane F13 to F12."""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

from audit_macfarlane_g20_dimension_reduction import (
    SOURCE_URL,
    build_maps,
    coefficient_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "macfarlane_f12_coordinate_pair_reduction.json"
)


def homogeneous_part(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    degree: int,
) -> sp.Expr:
    polynomial = sp.Poly(expression, *variables, domain=sp.QQ)
    return sp.expand(
        sum(
            coefficient
            * sp.prod(variable**exponent for variable, exponent in zip(variables, exponents))
            for exponents, coefficient in polynomial.terms()
            if sum(exponents) == degree
        )
    )


def exact_sparse_determinant(
    matrix: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
) -> tuple[dict[tuple[int, ...], Fraction], int]:
    """Laplace-expand a sparse polynomial matrix with exact Fraction arithmetic."""

    dimension = matrix.rows
    zero_exponents = (0,) * len(variables)

    def as_dict(expression: sp.Expr) -> dict[tuple[int, ...], Fraction]:
        return {
            exponents: Fraction(int(coefficient.p), int(coefficient.q))
            for exponents, coefficient in sp.Poly(
                expression, *variables, domain=sp.QQ
            ).terms()
            if coefficient
        }

    def add(
        left: dict[tuple[int, ...], Fraction],
        right: dict[tuple[int, ...], Fraction],
        sign: int,
    ) -> dict[tuple[int, ...], Fraction]:
        result = dict(left)
        for exponents, coefficient in right.items():
            updated = result.get(exponents, Fraction(0)) + sign * coefficient
            if updated:
                result[exponents] = updated
            else:
                result.pop(exponents, None)
        return result

    def multiply(
        left: dict[tuple[int, ...], Fraction],
        right: dict[tuple[int, ...], Fraction],
    ) -> dict[tuple[int, ...], Fraction]:
        result: dict[tuple[int, ...], Fraction] = {}
        for left_exponents, left_coefficient in left.items():
            for right_exponents, right_coefficient in right.items():
                exponents = tuple(
                    left_value + right_value
                    for left_value, right_value in zip(
                        left_exponents, right_exponents
                    )
                )
                result[exponents] = (
                    result.get(exponents, Fraction(0))
                    + left_coefficient * right_coefficient
                )
        return {
            exponents: coefficient
            for exponents, coefficient in result.items()
            if coefficient
        }

    entries = [
        [as_dict(matrix[row, column]) for column in range(dimension)]
        for row in range(dimension)
    ]
    row_order = sorted(
        range(dimension),
        key=lambda row: sum(bool(entry) for entry in entries[row]),
    )
    ordered = [entries[row] for row in row_order]
    inversions = sum(
        row_order[left] > row_order[right]
        for left in range(dimension)
        for right in range(left + 1, dimension)
    )
    permutation_sign = -1 if inversions % 2 else 1
    memo: dict[int, dict[tuple[int, ...], Fraction]] = {}

    def determinant(mask: int, row: int) -> dict[tuple[int, ...], Fraction]:
        if row == dimension:
            return {zero_exponents: Fraction(1)}
        if mask in memo:
            return memo[mask]
        result: dict[tuple[int, ...], Fraction] = {}
        position = 0
        for column in range(dimension):
            if mask & (1 << column):
                entry = ordered[row][column]
                if entry:
                    term = multiply(
                        entry,
                        determinant(mask ^ (1 << column), row + 1),
                    )
                    result = add(result, term, -1 if position % 2 else 1)
                position += 1
        memo[mask] = result
        return result

    value = determinant((1 << dimension) - 1, 0)
    return (
        {
            exponents: permutation_sign * coefficient
            for exponents, coefficient in value.items()
        },
        len(memo),
    )


def main() -> None:
    data = build_maps()
    x = data["x"]
    f13 = data["F13"]
    p13 = data["p13"]
    q13 = data["q13"]
    z = tuple(sp.symbols("z1:13"))
    s = sp.Symbol("s")

    # T(x)=(x1,...,x12,F13_13(x)) and
    # A(y)=(y1,y2,y3,y4-y8^2,y5,...,y13).
    assert f13[12] == x[12] + x[1] ** 2
    inverse_source = {x[index]: z[index] for index in range(12)}
    inverse_source[x[12]] = s - z[1] ** 2
    transformed = [
        sp.expand(component.subs(inverse_source, simultaneous=True))
        for component in f13
    ]
    transformed[3] = sp.expand(transformed[3] - transformed[7] ** 2)
    assert transformed[12] == s

    k12 = [sp.expand(component.subs(s, 0)) for component in transformed[:12]]
    expected_relative_defect = sp.zeros(13, 1)
    expected_relative_defect[3] = s * (2 * z[11] - z[0] ** 2)
    relative_difference = (
        sp.Matrix(transformed) - sp.Matrix(k12 + [s]) - expected_relative_defect
    )
    assert all(sp.expand(value) == 0 for value in relative_difference)

    # The two coordinate changes are triangular with determinant one.
    source_coordinates = list(z) + [s - z[1] ** 2]
    target_coordinates = list(sp.symbols("y1:14"))
    target_coordinates[3] -= target_coordinates[7] ** 2
    assert sp.Matrix(source_coordinates).jacobian(z + (s,)).det() == 1
    y = tuple(sp.symbols("y1:14"))
    assert sp.Matrix(target_coordinates).jacobian(y).det() == 1

    # Degree, tangent, collision, and a direct exact determinant calculation.
    jacobian12 = sp.Matrix(k12).jacobian(z)
    assert jacobian12.subs({variable: 0 for variable in z}) == sp.eye(12)
    component_degrees = [
        sp.Poly(component, *z, domain=sp.QQ).total_degree()
        for component in k12
    ]
    assert max(component_degrees) == 3
    nonlinear12 = [
        sp.expand(component - variable)
        for component, variable in zip(k12, z)
    ]
    nonlinear_degrees = {
        sum(exponents)
        for component in nonlinear12
        for exponents, coefficient in sp.Poly(
            component, *z, domain=sp.QQ
        ).terms()
        if coefficient
    }
    assert nonlinear_degrees == {2, 3}
    determinant12, determinant_memo_size = exact_sparse_determinant(jacobian12, z)
    assert determinant12 == {(0,) * 12: Fraction(1)}

    p12 = tuple(p13[index] for index in range(12))
    q12 = tuple(q13[index] for index in range(12))
    image12 = tuple(
        component.subs(dict(zip(z, p12))) for component in k12
    )
    assert p12 != q12
    assert image12 == tuple(
        component.subs(dict(zip(z, q12))) for component in k12
    )
    assert image12 == p12
    assert q13[12] + q13[1] ** 2 == 0

    quadratic = [homogeneous_part(value, z, 2) for value in nonlinear12]
    cubic = [homogeneous_part(value, z, 3) for value in nonlinear12]
    cubic_rank = coefficient_matrix(cubic, z).rank()
    assert cubic_rank == 6
    assert all(value == 0 for value in cubic[6:])
    assert coefficient_matrix(cubic[:6], z).rank() == 6

    # Rank-compressed homogeneous parent in 12+6+1=19 variables.
    w = tuple(sp.symbols("w1:7"))
    tau = sp.Symbol("tau")
    variables19 = z + w + (tau,)
    h19 = [
        sp.expand(
            tau * quadratic[index] + (tau**2 * w[index] if index < 6 else 0)
        )
        for index in range(12)
    ] + [-value for value in cubic[:6]] + [sp.Integer(0)]
    g19 = [
        sp.expand(variable + correction)
        for variable, correction in zip(variables19, h19)
    ]
    assert all(
        sp.Poly(value, *variables19, domain=sp.QQ).total_degree() == 3
        for value in h19
        if value
    )

    # Exact companion cancellation at tau=1 and the scaling identity
    # E_tau(z)=tau^{-1}K12(tau*z).  Together with det DK12=1 these prove
    # det DG19=1, including tau=0 by polynomial identity.
    e_tau = [
        sp.expand(
            z[index]
            + tau * quadratic[index]
            + tau**2 * cubic[index]
        )
        for index in range(12)
    ]
    scaled_k12 = [
        sp.cancel(
            component.subs(
                {variable: tau * variable for variable in z},
                simultaneous=True,
            )
            / tau
        )
        for component in k12
    ]
    assert e_tau == scaled_k12
    dehomogenized = [component.subs(tau, 1) for component in g19[:-1]]
    residual_companions = [
        sp.expand(w[index] - cubic[index]) for index in range(6)
    ]
    factorized = [
        sp.expand(
            k12[index]
            + (residual_companions[index] if index < 6 else 0)
        )
        for index in range(12)
    ] + residual_companions
    # M=A_B o (K12 x I_6) o S_C, with B the inclusion into the
    # first six output coordinates.
    assert dehomogenized == factorized

    cubic_at_p = tuple(value.subs(dict(zip(z, p12))) for value in cubic[:6])
    cubic_at_q = tuple(value.subs(dict(zip(z, q12))) for value in cubic[:6])
    p19 = p12 + cubic_at_p + (sp.Integer(1),)
    q19 = q12 + cubic_at_q + (sp.Integer(1),)
    image19_p = tuple(
        component.subs(dict(zip(variables19, p19))) for component in g19
    )
    image19_q = tuple(
        component.subs(dict(zip(variables19, q19))) for component in g19
    )
    assert p19 != q19 and image19_p == image19_q

    artifact = {
        "format": "macfarlane-f12-coordinate-pair-reduction-v1",
        "external_provenance": SOURCE_URL,
        "status": "exact theorem with independent standard-library replay",
        "source_target_reduction": {
            "source_coordinate": "s=x13+x2^2=F13_13",
            "target_shear": "y4 -> y4-y8^2",
            "relative_form": (
                "(K1,K2,K3,K4+s*(2*z12-z1^2),K5,...,K12,s)"
            ),
            "slice": "s=0",
        },
        "F12": {
            "dimension": 12,
            "component_degrees": component_degrees,
            "nonlinear_degrees": sorted(nonlinear_degrees),
            "determinant": "1",
            "determinant_memo_size": determinant_memo_size,
            "collision_p": [str(value) for value in p12],
            "collision_q": [str(value) for value in q12],
            "collision_image": [str(value) for value in image12],
            "cubic_output_rank": cubic_rank,
            "components": [str(value) for value in k12],
        },
        "G19": {
            "dimension": 19,
            "all_nonlinear_terms_homogeneous_degree": 3,
            "determinant": "1",
            "collision": True,
            "rank_compressed_count": "12+6+1=19",
            "variables": [str(value) for value in variables19],
            "nonlinear_components": [str(value) for value in h19],
            "collision_p": [str(value) for value in p19],
            "collision_q": [str(value) for value in q19],
            "collision_image": [str(value) for value in image19_p],
        },
        "consequences": {
            "direct_degree_three_upper_bound": 12,
            "cubic_homogeneous_upper_bound": 19,
            "homogeneous_quartic_HN_upper_bound": 38,
        },
        "scope": (
            "The construction is exact over Q. It proves upper bounds, not "
            "minimality, priority, formal verification, or external review."
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print("PASS F12: exact source-coordinate/target-square relative identity")
    print("PASS F12: degree <= 3, identity tangent, and exact rational collision")
    print(
        "PASS F12: direct sparse determinant is 1 "
        f"(memoized minors={determinant_memo_size})"
    )
    print("PASS F12: cubic-output rank is 6")
    print("PASS G19: cubic homogeneous parent, determinant bridge, and collision")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
