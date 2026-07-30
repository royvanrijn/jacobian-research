#!/usr/bin/env python3
"""Cross-audit the two dimension-38 quartic HN constructions.

The public Traboulsi construction starts from an 11-variable cubic Keller
collision whose cubic-output rank is seven.  The independent MacFarlane
reduction in this repository starts from twelve variables and rank six.
Both therefore have rank-compressed cost n+r=18.

This script checks two bounded attempts to lower that cost:

* quadratic pivot completions before deleting a graph coordinate of the
  public 11-variable lift;
* quadratic target completion of the nonlinear z8 pivot of the local map;
* coordinated quadratic source shears of the local 12-variable map.

All ranks are over QQ.  The searches are exact finite-dimensional
obstructions, not global lower bounds.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

from audit_macfarlane_g20_dimension_reduction import build_maps


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hvc38_cross_construction_frontier.json"
)
TRABOULSI_SOURCE = (
    "https://github.com/mtraboulsi689/"
    "Nested-Tail-components-for-Keller-Maps"
)
TRABOULSI_COMMIT = "73635c96034bb8364c036fbec2e366224e601b40"


def homogeneous_coefficients(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...] | list[sp.Symbol],
    degree: int,
) -> dict[tuple[int, ...], sp.Rational]:
    return {
        powers: coefficient
        for powers, coefficient in sp.Poly(
            expression, *variables, domain=sp.QQ
        ).terms()
        if sum(powers) == degree
    }


def coefficient_rank(
    expressions: list[sp.Expr],
    variables: tuple[sp.Symbol, ...] | list[sp.Symbol],
    degree: int,
) -> int:
    rows = [
        homogeneous_coefficients(expression, variables, degree)
        for expression in expressions
    ]
    monomials = sorted(set().union(*(row for row in rows)))
    return sp.Matrix(
        [[row.get(monomial, 0) for monomial in monomials] for row in rows]
    ).rank()


def traboulsi_g11() -> tuple[
    tuple[sp.Symbol, ...], list[sp.Expr], tuple[sp.Rational, ...]
]:
    x, y, z, a, b, c, d, q, s, h, k = sp.symbols(
        "x y z a b c d q s h k"
    )
    variables = (x, y, z, a, b, c, d, q, s, h, k)
    mapping = [
        -a * c
        - a * d * z
        - 3 * a * y**2
        - 2 * a * z
        - c * d**2
        + d**2 * z
        - d * s
        + 7 * d * y**2
        + s * x * y
        + 3 * x * y * z
        + 4 * y**2
        + z,
        -b * c
        - b * d * z
        - 3 * b * y**2
        - 2 * b * z
        - 3 * c * d * x
        - d * q
        + q * x * y
        + 12 * x * y**2
        + 3 * x * z
        + y,
        -h * k
        - h * x * z
        + k * x**2
        - 3 * x**2 * y
        + 2 * x,
        a - d**2 + 2 * d * x * y,
        b + 3 * x**2 * y,
        c + x * y * z + 3 * y**2 + 2 * z,
        d - x * y,
        b * z + 3 * c * x + q,
        s + a * z + c * x * y - x * y * z - 7 * y**2 + c * d - d * z,
        h - x**2,
        k + x * z,
    ]
    common_image = (-sp.Rational(1, 4),) + (sp.Rational(0),) * 10
    return variables, [sp.expand(value) for value in mapping], common_image


def local_f12() -> tuple[
    tuple[sp.Symbol, ...], list[sp.Expr], tuple[sp.Rational, ...]
]:
    data = build_maps()
    x = data["x"]
    f13 = data["F13"]
    z = tuple(sp.symbols("z1:13"))
    relative = sp.Symbol("relative")
    inverse_source = {x[index]: z[index] for index in range(12)}
    inverse_source[x[12]] = relative - z[1] ** 2
    transformed = [
        sp.expand(component.subs(inverse_source, simultaneous=True))
        for component in f13
    ]
    transformed[3] = sp.expand(transformed[3] - transformed[7] ** 2)
    mapping = [
        sp.expand(component.subs(relative, 0))
        for component in transformed[:12]
    ]
    common_image = tuple(data["p13"][:12])
    return z, mapping, common_image


def high_degree_system(
    base: list[sp.Expr],
    columns: list[list[sp.Expr]],
    variables: tuple[sp.Symbol, ...] | list[sp.Symbol],
    cutoff: int = 4,
) -> tuple[sp.MutableSparseMatrix, sp.MutableSparseMatrix]:
    base_polys = [sp.Poly(value, *variables, domain=sp.QQ) for value in base]
    column_polys = [
        [sp.Poly(value, *variables, domain=sp.QQ) for value in column]
        for column in columns
    ]
    rows: set[tuple[int, tuple[int, ...]]] = set()
    for output, polynomial in enumerate(base_polys):
        rows.update(
            (output, powers)
            for powers, coefficient in polynomial.terms()
            if coefficient and sum(powers) >= cutoff
        )
    for column in column_polys:
        for output, polynomial in enumerate(column):
            rows.update(
                (output, powers)
                for powers, coefficient in polynomial.terms()
                if coefficient and sum(powers) >= cutoff
            )
    ordered_rows = sorted(rows, key=str)
    matrix = sp.MutableSparseMatrix(
        len(ordered_rows), len(column_polys), {}
    )
    target = sp.MutableSparseMatrix(len(ordered_rows), 1, {})
    for row, (output, powers) in enumerate(ordered_rows):
        target[row, 0] = -base_polys[output].coeff_monomial(powers)
        for column, polynomials in enumerate(column_polys):
            coefficient = polynomials[output].coeff_monomial(powers)
            if coefficient:
                matrix[row, column] = coefficient
    return matrix, target


def quadratic_pivot_completion_audit(
    variables: tuple[sp.Symbol, ...],
    mapping: list[sp.Expr],
    common_image: tuple[sp.Rational, ...],
) -> list[dict[str, object]]:
    """Try g=Y_j+P(Y) with P linear/quadratic and x_j-affine."""

    records: list[dict[str, object]] = []
    for pivot in range(3, len(variables)):
        variable = variables[pivot]
        derivatives = [
            sp.expand(sp.diff(component, variable)) for component in mapping
        ]
        # The coordinate d occurs quadratically in other components, so a
        # completion of its pivot gives nonlinear coefficient equations and
        # lies outside this linear audit.
        if any(derivative.has(variable) for derivative in derivatives):
            continue
        allowed = [
            index
            for index, derivative in enumerate(derivatives)
            if derivative == 0 and index != pivot
        ]
        completions: list[sp.Expr] = []
        for degree in (1, 2):
            for indices in itertools.combinations_with_replacement(
                allowed, degree
            ):
                completions.append(
                    sp.expand(
                        sp.prod(mapping[index] for index in indices)
                        - sp.prod(common_image[index] for index in indices)
                    )
                )

        graph = sp.expand(
            common_image[pivot] - (mapping[pivot] - variable)
        )
        remaining = tuple(
            value
            for index, value in enumerate(variables)
            if index != pivot
        )
        base: list[sp.Expr] = []
        columns = [[] for _ in completions]
        for output, component in enumerate(mapping):
            if output == pivot:
                continue
            base.append(sp.expand(component.subs(variable, graph)))
            for column, completion in enumerate(completions):
                columns[column].append(
                    sp.expand(-derivatives[output] * completion)
                )
        matrix, target = high_degree_system(base, columns, remaining)
        rank = matrix.rank()
        augmented_rank = matrix.row_join(target).rank()
        assert augmented_rank == rank + 1
        records.append(
            {
                "pivot": str(variable),
                "allowed_target_outputs": len(allowed),
                "completion_columns": len(completions),
                "system_shape": list(matrix.shape),
                "rank": rank,
                "augmented_rank": augmented_rank,
            }
        )
    return records


def coordinated_source_shear_audit(
    variables: tuple[sp.Symbol, ...], mapping: list[sp.Expr]
) -> dict[str, object]:
    """Audit a maximal jointly affine five-coordinate shear block."""

    # One-based coordinates (z4,z7,z9,z10,z11).
    shear_indices = (3, 6, 8, 9, 10)
    base_variables = [
        variable
        for index, variable in enumerate(variables)
        if index not in shear_indices
    ]
    quadratics = [
        base_variables[left] * base_variables[right]
        for left, right in itertools.combinations_with_replacement(
            range(len(base_variables)), 2
        )
    ]
    columns: list[list[sp.Expr]] = []
    for index in shear_indices:
        derivatives = [
            sp.expand(sp.diff(component, variables[index]))
            for component in mapping
        ]
        for quadratic in quadratics:
            columns.append(
                [
                    sp.expand(-derivative * quadratic)
                    for derivative in derivatives
                ]
            )

    zero_base = [sp.Integer(0)] * len(mapping)
    high_matrix, _ = high_degree_system(zero_base, columns, variables)
    kernel = high_matrix.nullspace()
    assert high_matrix.shape == (146, 140)
    assert high_matrix.rank() == 126
    assert len(kernel) == 14

    correction_maps = [
        [
            sp.expand(
                sum(
                    vector[column] * columns[column][output]
                    for column in range(len(columns))
                )
            )
            for output in range(len(mapping))
        ]
        for vector in kernel
    ]
    base_cubic = [
        homogeneous_coefficients(component, variables, 3)
        for component in mapping
    ]
    correction_cubics = [
        [
            homogeneous_coefficients(component, variables, 3)
            for component in correction
        ]
        for correction in correction_maps
    ]
    monomials = sorted(
        set().union(
            *(row for row in base_cubic),
            *(
                row
                for correction in correction_cubics
                for row in correction
            ),
        )
    )
    base_matrix = sp.Matrix(
        [
            [row.get(monomial, 0) for monomial in monomials]
            for row in base_cubic
        ]
    )
    correction_matrices = [
        sp.Matrix(
            [
                [row.get(monomial, 0) for monomial in monomials]
                for row in correction
            ]
        )
        for correction in correction_cubics
    ]

    invariant_columns = [
        column
        for column in range(len(monomials))
        if all(
            matrix[row, column] == 0
            for matrix in correction_matrices
            for row in range(6)
        )
    ]
    invariant_block = base_matrix[:6, invariant_columns]
    pivot_column_positions = list(invariant_block.rref()[1])
    pivot_columns = [
        invariant_columns[position] for position in pivot_column_positions
    ]
    pivot_rows = list(
        invariant_block[:, pivot_column_positions].T.rref()[1]
    )
    pivot = base_matrix.extract(pivot_rows, pivot_columns)
    assert pivot.shape == (5, 5)
    assert pivot.det() == 6

    parameters = sp.symbols(f"t0:{len(kernel)}")
    cubic_matrix = base_matrix + sum(
        (
            parameter * correction
            for parameter, correction in zip(
                parameters, correction_matrices, strict=True
            )
        ),
        sp.zeros(*base_matrix.shape),
    )
    pivot_inverse = pivot.inv()
    equations: list[sp.Expr] = []
    for row in range(cubic_matrix.rows):
        if row in pivot_rows:
            continue
        for column in range(cubic_matrix.cols):
            if column in pivot_columns:
                continue
            equation = sp.factor(
                cubic_matrix[row, column]
                - (
                    cubic_matrix.extract([row], pivot_columns)
                    * pivot_inverse
                    * cubic_matrix.extract(pivot_rows, [column])
                )[0]
            )
            if equation:
                polynomial = sp.Poly(equation, *parameters, domain=sp.QQ)
                equations.append(sp.expand(equation / polynomial.LC()))
    unique_equations = list(dict.fromkeys(equations))
    linear_equations = [
        equation
        for equation in unique_equations
        if sp.Poly(equation, *parameters).total_degree() <= 1
    ]
    linear_matrix, linear_target = sp.linear_eq_to_matrix(
        linear_equations, parameters
    )
    linear_rank = linear_matrix.rank()
    augmented_rank = linear_matrix.row_join(linear_target).rank()
    assert len(linear_equations) == 6
    assert (linear_rank, augmented_rank) == (5, 6)

    return {
        "sheared_source_coordinates": [
            str(variables[index]) for index in shear_indices
        ],
        "quadratic_base_variables": [
            str(variable) for variable in base_variables
        ],
        "high_degree_system_shape": list(high_matrix.shape),
        "high_degree_rank": high_matrix.rank(),
        "degree_preserving_kernel_dimension": len(kernel),
        "fixed_cubic_minor": {
            "size": 5,
            "determinant": str(pivot.det()),
        },
        "rank_at_most_five_equations": len(unique_equations),
        "rank_at_most_five_linear_equations": len(linear_equations),
        "rank_at_most_five_linear_rank": linear_rank,
        "rank_at_most_five_augmented_rank": augmented_rank,
        "conclusion": (
            "No member of this 14-parameter degree-preserving coordinated "
            "quadratic source-shear family has cubic-output rank at most 5."
        ),
    }


def nonlinear_z8_pivot_audit(
    variables: tuple[sp.Symbol, ...], mapping: list[sp.Expr]
) -> dict[str, object]:
    """Obstruct quadratic target completion of the non-affine z8 pivot."""

    pivot = 7
    allowed = [
        index
        for index, component in enumerate(mapping)
        if index != pivot and sp.diff(component, variables[pivot]) == 0
    ]
    basis: list[sp.Expr] = []
    for degree in (1, 2):
        for indices in itertools.combinations_with_replacement(
            allowed, degree
        ):
            basis.append(
                sp.expand(sp.prod(mapping[index] for index in indices))
            )
    polynomials = [
        sp.Poly(value, *variables, domain=sp.QQ) for value in basis
    ]
    high_monomials = sorted(
        set().union(
            *(
                {
                    powers
                    for powers, coefficient in polynomial.terms()
                    if coefficient and sum(powers) >= 3
                }
                for polynomial in polynomials
            )
        )
    )
    high_matrix = sp.Matrix(
        [
            [
                polynomial.coeff_monomial(monomial)
                for polynomial in polynomials
            ]
            for monomial in high_monomials
        ]
    )
    kernel = high_matrix.nullspace()
    assert high_matrix.shape == (206, 44)
    assert high_matrix.rank() == 39
    assert len(kernel) == 5

    target_monomial = (1, 1) + (0,) * 10
    degree_two_monomials = sorted(
        {target_monomial}
        | set().union(
            *(
                {
                    powers
                    for powers, coefficient in polynomial.terms()
                    if coefficient and sum(powers) == 2
                }
                for polynomial in polynomials
            )
        )
    )
    degree_two_matrix = sp.Matrix(
        [
            [
                polynomial.coeff_monomial(monomial)
                for polynomial in polynomials
            ]
            for monomial in degree_two_monomials
        ]
    ) * sp.Matrix.hstack(*kernel)
    target = sp.zeros(len(degree_two_monomials), 1)
    target[degree_two_monomials.index(target_monomial), 0] = 1
    rank = degree_two_matrix.rank()
    augmented_rank = degree_two_matrix.row_join(target).rank()
    assert (rank, augmented_rank) == (5, 6)

    # If the completed pivot is g=Y8+P(Y_allowed), put R=P(K_allowed).
    # On g=0, z8=-z1*z2-R and K4=z4+(z1*z2)^2-R^2.  Degree <=3 first
    # forces deg(R)<=2 and then R_2=+/-z1*z2 in the UFD QQ[z].
    return {
        "pivot": "z8",
        "allowed_target_outputs": [index + 1 for index in allowed],
        "target_degree_bound": 2,
        "completion_basis_size": len(basis),
        "degree_at_least_three_system_shape": list(high_matrix.shape),
        "degree_at_least_three_rank": high_matrix.rank(),
        "degree_at_most_two_pullback_dimension": len(kernel),
        "quadratic_part_rank": rank,
        "quadratic_part_augmented_rank_for_z1z2": augmented_rank,
        "identity_on_pivot_slice": (
            "K4=z4+(z1*z2)^2-P(K_allowed)^2"
        ),
        "conclusion": (
            "No target completion of degree at most 2 in the eight "
            "z8-independent output components can make the z8 pivot "
            "restriction degree at most 3."
        ),
    }


def main() -> None:
    g11_variables, g11_mapping, g11_image = traboulsi_g11()
    f12_variables, f12_mapping, f12_image = local_f12()

    assert max(
        sp.Poly(component, *g11_variables).total_degree()
        for component in g11_mapping
    ) == 3
    assert max(
        sp.Poly(component, *f12_variables).total_degree()
        for component in f12_mapping
    ) == 3
    g11_rank = coefficient_rank(g11_mapping, g11_variables, 3)
    f12_rank = coefficient_rank(f12_mapping, f12_variables, 3)
    assert (g11_rank, f12_rank) == (7, 6)

    g11_completions = quadratic_pivot_completion_audit(
        g11_variables, g11_mapping, g11_image
    )
    coordinated_shears = coordinated_source_shear_audit(
        f12_variables, f12_mapping
    )
    nonlinear_z8 = nonlinear_z8_pivot_audit(
        f12_variables, f12_mapping
    )

    artifact = {
        "format": "hvc38-cross-construction-frontier-v1",
        "field": "QQ",
        "public_construction": {
            "source": TRABOULSI_SOURCE,
            "commit": TRABOULSI_COMMIT,
            "degree_three_dimension": 11,
            "cubic_output_rank": g11_rank,
            "compressed_homogeneous_dimension": 19,
            "quartic_HN_dimension": 38,
        },
        "local_independent_construction": {
            "canonical_source": (
                "verified/"
                "TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md"
            ),
            "degree_three_dimension": 12,
            "cubic_output_rank": f12_rank,
            "compressed_homogeneous_dimension": 19,
            "quartic_HN_dimension": 38,
        },
        "next_standard_doubling_target": {
            "required_n_plus_cubic_rank_at_most": 17,
            "required_cubic_homogeneous_dimension_at_most": 18,
            "next_even_quartic_HN_dimension": 36,
            "note": (
                "The ordinary homogeneous cotangent lift doubles dimension, "
                "so this pipeline cannot produce dimension 37."
            ),
        },
        "public_G11_quadratic_pivot_completion_obstructions": g11_completions,
        "local_F12_nonlinear_z8_pivot": nonlinear_z8,
        "local_F12_coordinated_source_shear": coordinated_shears,
        "scope": (
            "Exact bounded obstructions only. They do not exclude nonlinear "
            "coordinate pairs outside the stated bases, higher-degree "
            "automorphisms, non-nested state realizations, or a symmetric "
            "lift that avoids full dimension doubling."
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print("PASS HVC38 cross-audit: public profile (11,7), local profile (12,6)")
    print(
            "PASS public G11: all seven affine quadratic "
            "pivot-completion systems "
            "are inconsistent"
    )
    print(
        "PASS local F12: coordinated shear kernel has dimension 14, "
        "but rank<=5 equations have ranks (5,6)"
    )
    print(
        "PASS local F12: nonlinear z8 pivot has no quadratic target "
        "completion to degree 3"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    sys.exit(main())
