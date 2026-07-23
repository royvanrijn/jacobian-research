#!/usr/bin/env python3
"""Compress the fourth-order obstruction on the residual five-space.

The zero-scale third-order affine equation admits one second correction that
is independent of the residual first direction.  This script constructs that
global base and studies the fixed three-monomial fourth-order residue on the
remaining lower-lift variables.
"""

from __future__ import annotations

import argparse

from sympy.polys.domains import QQ
import sympy as sp
from sympy.polys.matrices.sdm import (
    sdm_irref,
    sdm_nullspace_from_rref,
)

from explore_degree_five_a2_subprincipal import (
    add,
    filtered_monomials,
    pi_power,
    poisson,
    scale,
)
from explore_degree_five_quantum_residue import solve_affine
from explore_rank_two_odd_mixed_quantization import (
    coupling,
    essential_problem,
    extend_sparse_poly,
    linear_combination,
    solve_many_particular,
    split_correction,
)
from explore_rank_two_odd_residual_five_space import (
    RESIDUAL_VECTORS,
    add_sparse,
)
from explore_rank_two_odd_support_three_curves import (
    project_to_third_cokernel,
)
from verify_rank_two_odd_support_three_points import lower_lift_data
from verify_rank_two_odd_mixed_function_field import (
    RESIDUE_WEIGHTS,
    residue,
)


UNIFORM_KERNEL_CORRECTION = {
    18: -QQ(710122059) / QQ(392),
    33: -QQ(922185) / QQ(56),
}

KERNEL_PLANE_RESIDUE_WEIGHTS = {
    (14, 0, 1): QQ(10949113363) / QQ(486999360892800),
    (15, 1, 1): QQ(26951510851) / QQ(19723474116158400),
    (16, 0, 2): QQ(19704634546993) / QQ(966450231691761600),
    (16, 2, 1): QQ(5150073181) / QQ(34516079703277200),
    (16, 4, 0): -QQ(111385294) / QQ(2773613547584775),
    (17, 1, 2): -QQ(8772384231703) / QQ(2899350695075284800),
    (17, 3, 1): QQ(71466318533) / QQ(310644717329494800),
    (17, 5, 0): -QQ(77120461) / QQ(3328336257101730),
    (18, 0, 3): QQ(13486366961461) / QQ(13530303243684662400),
    (18, 2, 2): QQ(892774749523) / QQ(724837673768821200),
    (18, 4, 1): -QQ(1420550843) / QQ(10354823910983160),
    (18, 6, 0): QQ(19983352) / QQ(1664168128550865),
    (19, 1, 3): QQ(1085857489637) / QQ(120806278961470200),
    (19, 3, 2): -QQ(16007235857) / QQ(20709647821966320),
    (19, 5, 1): QQ(31917346) / QQ(554722709516955),
    (19, 7, 0): -QQ(190658) / QQ(47547660815739),
}


def residual_projected_data():
    """Return the five projected coupling rows and affine right sides."""

    projected_couplings, projected_right_sides, _ = (
        project_to_third_cokernel()
    )
    residual_couplings = []
    residual_right_sides = []
    for vector in RESIDUAL_VECTORS:
        columns = []
        for kernel_index in range(42):
            column = {}
            for essential_index, coefficient in vector.items():
                column = add_sparse(
                    column,
                    projected_couplings[essential_index][kernel_index],
                    coefficient,
                )
            columns.append(column)
        residual_couplings.append(columns)

        right_side = {}
        for essential_index, coefficient in vector.items():
            right_side = add_sparse(
                right_side,
                projected_right_sides[essential_index],
                coefficient,
            )
        residual_right_sides.append(right_side)
    return residual_couplings, residual_right_sides


def verify_uniform_correction():
    """Prove that one 42-vector solves all five affine equations."""

    residual_couplings, residual_right_sides = residual_projected_data()
    columns = [
        {
            (residual_index, coordinate): coefficient
            for residual_index in range(5)
            for coordinate, coefficient in residual_couplings[
                residual_index
            ][kernel_index].items()
        }
        for kernel_index in range(42)
    ]
    right_side = {
        (residual_index, coordinate): coefficient
        for residual_index, vector in enumerate(residual_right_sides)
        for coordinate, coefficient in vector.items()
    }
    particular, kernel, rank = solve_affine(columns, right_side, QQ)
    assert particular == UNIFORM_KERNEL_CORRECTION
    assert (rank, len(kernel)) == (15, 27)
    return residual_couplings


def matrix_rank(rows):
    """Return the exact rank of a sparse row dictionary."""

    _, pivots, _ = sdm_irref(
        {row: entries for row, entries in rows.items() if entries}
    )
    return len(pivots)


def residue_pairing_profile(S, T, residual_pairs):
    """Measure the small scalar pairings seen by the fixed residue."""

    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)

    poisson_rows = {}
    poisson_count = 0
    for left, s_monomial in enumerate(s2_monomials):
        row = {}
        for right, t_monomial in enumerate(t2_monomials):
            value = residue(
                poisson(
                    {s_monomial: QQ.one},
                    {t_monomial: QQ.one},
                ),
                QQ,
            )
            if value:
                row[right] = value
                poisson_count += 1
        if row:
            poisson_rows[left] = row

    coupling_rows = {}
    coupling_count = 0
    third_monomials = s3_monomials + t3_monomials
    for residual_index, first_pair in enumerate(residual_pairs):
        row = {}
        for column, monomial in enumerate(third_monomials):
            if column < len(s3_monomials):
                third_pair = ({monomial: QQ.one}, {})
            else:
                third_pair = ({}, {monomial: QQ.one})
            value = residue(coupling(first_pair, third_pair), QQ)
            if value:
                row[column] = value
                coupling_count += 1
        if row:
            coupling_rows[residual_index] = row

    linear_second = {}
    for column, monomial in enumerate(s2_monomials + t2_monomials):
        if column < len(s2_monomials):
            value = residue(
                pi_power({monomial: QQ.one}, T, 3),
                QQ,
            )
        else:
            value = residue(
                pi_power(S, {monomial: QQ.one}, 3),
                QQ,
            )
        if value:
            linear_second[column] = value

    return {
        "poisson_nonzero": poisson_count,
        "poisson_rank": matrix_rank(poisson_rows),
        "coupling_nonzero": coupling_count,
        "coupling_rank": matrix_rank(coupling_rows),
        "linear_second_support": len(linear_second),
    }


def global_zero_scale_base(S, T, residual_pairs):
    """Construct the direction-independent second base and linear thirds."""

    parity_pair, kernel_pairs, third_columns = lower_lift_data(S, T)
    second_base = parity_pair
    for index, coefficient in UNIFORM_KERNEL_CORRECTION.items():
        second_base = (
            add(second_base[0], kernel_pairs[index][0], coefficient),
            add(second_base[1], kernel_pairs[index][1], coefficient),
        )

    third_pairs = []
    for first_pair in residual_pairs:
        rhs = scale(
            add(
                pi_power(first_pair[0], T, 3),
                pi_power(S, first_pair[1], 3),
            ),
            -QQ.one / QQ(24),
        )
        rhs = add(rhs, coupling(first_pair, second_base), -QQ.one)
        vector, _, rank = solve_affine(third_columns, rhs, QQ)
        assert rank == 1034
        s3_monomials = filtered_monomials(23, 2)
        t3_monomials = filtered_monomials(19, 1)
        split = len(s3_monomials)
        third_pairs.append(
            (
                {
                    s3_monomials[index]: value
                    for index, value in vector.items()
                    if index < split
                },
                {
                    t3_monomials[index - split]: value
                    for index, value in vector.items()
                    if index >= split
                },
            )
        )
    return second_base, third_pairs


def rational_residue(poly, weights, field):
    return sum(
        (
            field.convert(weight) * poly.get(monomial, field.zero)
            for monomial, weight in weights.items()
        ),
        field.zero,
    )


def zero_scale_constant_residue(
    S,
    T,
    second_base,
    weights=RESIDUE_WEIGHTS,
):
    """Evaluate the direction-independent fourth defect at u=0."""

    value = poisson(second_base[0], second_base[1])
    cubic = add(
        pi_power(S, second_base[1], 3),
        pi_power(second_base[0], T, 3),
    )
    value = add(value, cubic, QQ.one / QQ(24))
    value = add(
        value,
        pi_power(S, T, 5),
        QQ.one / QQ(1920),
    )
    return rational_residue(value, weights, QQ)


def residual_quadratic_pairs(S, T, residual_pairs):
    """Return particular second corrections for the 15 quadratic terms."""

    s2_monomials = filtered_monomials(25, 3)
    t2_monomials = filtered_monomials(21, 2)
    second_columns = [
        poisson({monomial: QQ.one}, T)
        for monomial in s2_monomials
    ]
    second_columns += [
        poisson(S, {monomial: QQ.one})
        for monomial in t2_monomials
    ]
    indices = []
    right_sides = []
    for left in range(5):
        for right in range(left, 5):
            if left == right:
                bracket = poisson(
                    residual_pairs[left][0],
                    residual_pairs[left][1],
                )
            else:
                bracket = add(
                    poisson(
                        residual_pairs[left][0],
                        residual_pairs[right][1],
                    ),
                    poisson(
                        residual_pairs[right][0],
                        residual_pairs[left][1],
                    ),
                )
            indices.append((left, right))
            right_sides.append(scale(bracket, -QQ.one))
    vectors = solve_many_particular(second_columns, right_sides, QQ)
    return indices, [
        split_correction(vector, s2_monomials, t2_monomials)
        for vector in vectors
    ]


def second_kernel_residue_forms(
    S,
    T,
    second_base,
    kernel_pairs,
    weights=RESIDUE_WEIGHTS,
):
    """Return the scalar linear and quadratic forms on ker(d_2)."""

    linear = []
    diagonal = []
    cross = {}
    for index, kernel_pair in enumerate(kernel_pairs):
        linear_value = add(
            poisson(second_base[0], kernel_pair[1]),
            poisson(kernel_pair[0], second_base[1]),
        )
        cubic = add(
            pi_power(S, kernel_pair[1], 3),
            pi_power(kernel_pair[0], T, 3),
        )
        linear_value = add(
            linear_value,
            cubic,
            QQ.one / QQ(24),
        )
        linear.append(rational_residue(linear_value, weights, QQ))
        diagonal.append(
            rational_residue(
                poisson(kernel_pair[0], kernel_pair[1]),
                weights,
                QQ,
            )
        )
        for right in range(index):
            value = add(
                poisson(kernel_pairs[right][0], kernel_pair[1]),
                poisson(kernel_pair[0], kernel_pairs[right][1]),
            )
            cross[(right, index)] = rational_residue(
                value,
                weights,
                QQ,
            )
    return linear, diagonal, cross


def symbolic_pair(pair, field):
    return (
        {
            monomial: field.to_sympy(coefficient)
            for monomial, coefficient in pair[0].items()
        },
        {
            monomial: field.to_sympy(coefficient)
            for monomial, coefficient in pair[1].items()
        },
    )


def symbolic_pair_combination(pairs, coefficients, field):
    result = ({}, {})
    for pair, coefficient in zip(pairs, coefficients):
        converted = symbolic_pair(pair, field)
        result = (
            add(result[0], converted[0], coefficient),
            add(result[1], converted[1], coefficient),
        )
    return result


def symbolic_residue(poly, weights=RESIDUE_WEIGHTS):
    return sp.expand(
        sum(
            sp.sympify(weight) * poly.get(monomial, sp.S.Zero)
            for monomial, weight in weights.items()
        )
    )


def decompose_parametric_sparse(poly, parameters, field):
    """Split a sparse physical polynomial by parameter monomial."""

    result = {}
    for physical_monomial, expression in poly.items():
        parameter_poly = sp.Poly(
            sp.expand(expression),
            *parameters,
            domain=field,
        )
        for parameter_monomial, coefficient in parameter_poly.terms():
            converted = field.convert(coefficient)
            if converted:
                result.setdefault(parameter_monomial, {})[
                    physical_monomial
                ] = converted
    return result


def solve_parametric_many(columns, right_sides, parameters, field):
    """Apply one fixed linear right inverse to parametric right sides."""

    decompositions = [
        decompose_parametric_sparse(rhs, parameters, field)
        for rhs in right_sides
    ]
    tags = [
        (source, parameter_monomial)
        for source, decomposition in enumerate(decompositions)
        for parameter_monomial in sorted(decomposition)
        if decomposition[parameter_monomial]
    ]
    assert all(
        coefficient
        for column in columns
        for coefficient in column.values()
    )
    assert all(
        coefficient
        for decomposition in decompositions
        for right_side in decomposition.values()
        for coefficient in right_side.values()
    )
    particulars = solve_many_particular(
        columns,
        [
            decompositions[source][parameter_monomial]
            for source, parameter_monomial in tags
        ],
        field,
    )
    results = [{} for _ in right_sides]
    for (source, parameter_monomial), particular in zip(
        tags,
        particulars,
    ):
        parameter_term = sp.prod(
            parameter**exponent
            for parameter, exponent in zip(
                parameters,
                parameter_monomial,
            )
        )
        for index, coefficient in particular.items():
            results[source][index] = (
                results[source].get(index, sp.S.Zero)
                + field.to_sympy(coefficient) * parameter_term
            )
    return results


def symbolic_split(vector, s_monomials, t_monomials):
    split = len(s_monomials)
    return (
        {
            s_monomials[index]: value
            for index, value in vector.items()
            if index < split and value != 0
        },
        {
            t_monomials[index - split]: value
            for index, value in vector.items()
            if index >= split and value != 0
        },
    )


def assert_zero(expression, parameters, field):
    polynomial = sp.Poly(
        sp.expand(expression),
        *parameters,
        domain=field,
    )
    assert polynomial.is_zero, polynomial.as_expr()


def coupling_kernel(coupling_columns, field):
    """Return an exact basis of the kernel of a sparse coupling map."""

    coordinates = sorted(
        set().union(*(set(column) for column in coupling_columns))
    )
    rows = {
        row: {
            column: value
            for column, vector in enumerate(coupling_columns)
            if (value := vector.get(coordinate, field.zero))
        }
        for row, coordinate in enumerate(coordinates)
    }
    reduced, pivots, nonzero = sdm_irref(rows)
    kernel, _ = sdm_nullspace_from_rref(
        reduced,
        field.one,
        len(coupling_columns),
        pivots,
        nonzero,
    )
    return len(pivots), kernel


def evaluate_linear(vector, coefficients, field):
    return sum(
        (
            value * field.convert(coefficients[index])
            for index, value in vector.items()
        ),
        field.zero,
    )


def evaluate_quadratic(vector, diagonal, cross, field):
    value = sum(
        (
            coefficient * coefficient * field.convert(diagonal[index])
            for index, coefficient in vector.items()
        ),
        field.zero,
    )
    indices = sorted(vector)
    for left_position, left in enumerate(indices):
        for right in indices[left_position + 1 :]:
            value += (
                vector[left]
                * vector[right]
                * field.convert(cross[(left, right)])
            )
    return value


def evaluate_bilinear(left, right, diagonal, cross, field):
    value = sum(
        (
            left.get(index, field.zero)
            * right.get(index, field.zero)
            * field.convert(2 * diagonal[index])
            for index in set(left).intersection(right)
        ),
        field.zero,
    )
    for left_index, left_value in left.items():
        for right_index, right_value in right.items():
            if left_index == right_index:
                continue
            pair = tuple(sorted((left_index, right_index)))
            value += (
                left_value
                * right_value
                * field.convert(cross[pair])
            )
    return value


def verify_zero_scale_kernel_residue(
    residual_couplings,
    linear,
    diagonal,
    cross,
):
    """Prove residue invariance on the generic coupling pencil and A-chart."""

    parameter = sp.Symbol("t")
    function_field = QQ.frac_field(parameter)
    parameter_value = function_field.from_sympy(parameter)
    pencil_columns = []
    for kernel_index in range(42):
        column = {}
        for coordinate in set(
            residual_couplings[0][kernel_index]
        ).union(residual_couplings[2][kernel_index]):
            value = (
                function_field.convert(
                    residual_couplings[0][kernel_index].get(
                        coordinate,
                        QQ.zero,
                    )
                )
                * parameter_value
                + function_field.convert(
                    residual_couplings[2][kernel_index].get(
                        coordinate,
                        QQ.zero,
                    )
                )
            )
            if value:
                column[coordinate] = value
        pencil_columns.append(column)

    rank, kernel = coupling_kernel(pencil_columns, function_field)
    assert (rank, len(kernel)) == (11, 31)
    assert all(
        evaluate_linear(vector, linear, function_field)
        == function_field.zero
        for vector in kernel
    )
    assert all(
        evaluate_quadratic(vector, diagonal, cross, function_field)
        == function_field.zero
        for vector in kernel
    )
    assert all(
        evaluate_bilinear(
            kernel[left],
            kernel[right],
            diagonal,
            cross,
            function_field,
        )
        == function_field.zero
        for left in range(len(kernel))
        for right in range(left + 1, len(kernel))
    )

    a_rank, a_kernel = coupling_kernel(residual_couplings[0], QQ)
    assert (a_rank, len(a_kernel)) == (11, 31)
    assert all(
        evaluate_linear(vector, linear, QQ) == QQ.zero
        for vector in a_kernel
    )
    assert all(
        evaluate_quadratic(vector, diagonal, cross, QQ) == QQ.zero
        for vector in a_kernel
    )
    assert all(
        evaluate_bilinear(
            a_kernel[left],
            a_kernel[right],
            diagonal,
            cross,
            QQ,
        )
        == QQ.zero
        for left in range(len(a_kernel))
        for right in range(left + 1, len(a_kernel))
    )
    return len(kernel), len(a_kernel)


def verify_scale_component(
    *,
    name,
    field,
    a_value,
    b_value,
    S,
    T,
    residual_pairs,
    residual_couplings,
    quadratic_indices,
    quadratic_pairs,
    second_base,
    third_bases,
    second_kernel_pairs,
    third_kernel_pairs,
    third_columns,
    project,
):
    """Prove the fourth residue identity on one normalized hyperplane."""

    parameters = sp.symbols(f"{name}_w0:3")
    w0, w1, w2 = parameters
    a_expression = field.to_sympy(field.convert(a_value))
    b_expression = field.to_sympy(field.convert(b_value))
    z = (
        a_expression
        - sp.Rational(4, 3) * w0
        - sp.Rational(2648, 189) * w1
        + sp.Rational(580504, 735) * w2,
        w0,
        b_expression
        - 2 * w1
        + sp.Rational(9838, 105) * w2,
        w1,
        w2,
    )

    field_S = extend_sparse_poly(S, field)
    field_T = extend_sparse_poly(T, field)
    field_residual_pairs = [
        (
            extend_sparse_poly(pair[0], field),
            extend_sparse_poly(pair[1], field),
        )
        for pair in residual_pairs
    ]
    field_quadratic_pairs = [
        (
            extend_sparse_poly(pair[0], field),
            extend_sparse_poly(pair[1], field),
        )
        for pair in quadratic_pairs
    ]
    field_second_base = (
        extend_sparse_poly(second_base[0], field),
        extend_sparse_poly(second_base[1], field),
    )
    field_third_bases = [
        (
            extend_sparse_poly(pair[0], field),
            extend_sparse_poly(pair[1], field),
        )
        for pair in third_bases
    ]
    field_second_kernels = [
        (
            extend_sparse_poly(pair[0], field),
            extend_sparse_poly(pair[1], field),
        )
        for pair in second_kernel_pairs
    ]
    field_third_kernels = [
        (
            extend_sparse_poly(pair[0], field),
            extend_sparse_poly(pair[1], field),
        )
        for pair in third_kernel_pairs
    ]
    field_third_columns = [
        extend_sparse_poly(column, field)
        for column in third_columns
    ]

    first_pair = symbolic_pair_combination(
        field_residual_pairs,
        z,
        field,
    )
    second_quadratic = symbolic_pair_combination(
        field_quadratic_pairs,
        [
            z[left] * z[right]
            for left, right in quadratic_indices
        ],
        field,
    )
    third_base = symbolic_pair_combination(
        field_third_bases,
        z,
        field,
    )

    fixed_couplings = []
    for kernel_index in range(42):
        column = {}
        for coordinate in set(
            residual_couplings[0][kernel_index]
        ).union(residual_couplings[2][kernel_index]):
            value = (
                field.convert(
                    residual_couplings[0][kernel_index].get(
                        coordinate,
                        QQ.zero,
                    )
                )
                * field.convert(a_value)
                + field.convert(
                    residual_couplings[2][kernel_index].get(
                        coordinate,
                        QQ.zero,
                    )
                )
                * field.convert(b_value)
            )
            if value:
                column[coordinate] = value
        fixed_couplings.append(column)
    coupling_rank, zero_u_kernel = coupling_kernel(
        fixed_couplings,
        field,
    )
    assert (coupling_rank, len(zero_u_kernel)) == (11, 31)

    projected_cubic = project(coupling(first_pair, second_quadratic))
    projected_rhs = {
        coordinate: -value
        for coordinate, value in projected_cubic.items()
    }
    kernel_solution = solve_parametric_many(
        fixed_couplings,
        [projected_rhs],
        parameters,
        field,
    )[0]
    second_u = second_quadratic
    for kernel_index, coefficient in kernel_solution.items():
        kernel_pair = symbolic_pair(
            field_second_kernels[kernel_index],
            field,
        )
        second_u = (
            add(second_u[0], kernel_pair[0], coefficient),
            add(second_u[1], kernel_pair[1], coefficient),
        )

    zero_u_second_pairs = []
    for vector in zero_u_kernel:
        zero_u_second_pairs.append(
            symbolic_pair_combination(
                field_second_kernels,
                [
                    field.to_sympy(vector.get(index, field.zero))
                    for index in range(42)
                ],
                field,
            )
        )

    third_right_sides = [
        scale(coupling(first_pair, second_u), -sp.S.One)
    ]
    third_right_sides.extend(
        scale(coupling(first_pair, pair), -sp.S.One)
        for pair in zero_u_second_pairs
    )
    third_vectors = solve_parametric_many(
        field_third_columns,
        third_right_sides,
        parameters,
        field,
    )
    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)
    third_u = symbolic_split(
        third_vectors[0],
        s3_monomials,
        t3_monomials,
    )
    zero_u_third_pairs = [
        symbolic_split(vector, s3_monomials, t3_monomials)
        for vector in third_vectors[1:]
    ]

    diagonal_u = add(
        poisson(second_u[0], second_u[1]),
        coupling(first_pair, third_u),
    )
    assert_zero(
        symbolic_residue(diagonal_u),
        parameters,
        field,
    )

    field_base_symbolic = symbolic_pair(field_second_base, field)
    linear_u = add(
        poisson(field_base_symbolic[0], second_u[1]),
        poisson(second_u[0], field_base_symbolic[1]),
    )
    cubic_u = add(
        pi_power(
            {
                monomial: field.to_sympy(value)
                for monomial, value in field_S.items()
            },
            second_u[1],
            3,
        ),
        pi_power(
            second_u[0],
            {
                monomial: field.to_sympy(value)
                for monomial, value in field_T.items()
            },
            3,
        ),
    )
    linear_u = add(linear_u, cubic_u, sp.Rational(1, 24))
    linear_u = add(linear_u, coupling(first_pair, third_base))
    linear_u = add(
        linear_u,
        pi_power(first_pair[0], first_pair[1], 3),
        sp.Rational(1, 24),
    )
    assert_zero(
        symbolic_residue(linear_u),
        parameters,
        field,
    )

    for second_pair, third_pair in zip(
        zero_u_second_pairs,
        zero_u_third_pairs,
    ):
        cross = add(
            poisson(second_u[0], second_pair[1]),
            poisson(second_pair[0], second_u[1]),
        )
        cross = add(cross, coupling(first_pair, third_pair))
        assert_zero(
            symbolic_residue(cross),
            parameters,
            field,
        )

    for third_kernel_pair in field_third_kernels:
        symbolic_third_kernel = symbolic_pair(
            third_kernel_pair,
            field,
        )
        assert_zero(
            symbolic_residue(
                coupling(first_pair, symbolic_third_kernel)
            ),
            parameters,
            field,
        )
    return {
        "field": str(field),
        "parameters": len(parameters),
        "zero_u_second_kernel": len(zero_u_kernel),
        "third_gauge_kernel": len(third_kernel_pairs),
    }


def verify_kernel_plane(
    *,
    S,
    T,
    residual_pairs,
    quadratic_indices,
    quadratic_pairs,
    second_base,
    third_bases,
    second_kernel_pairs,
    third_kernel_pairs,
    third_columns,
):
    """Prove a second residue identity on the rank-zero kernel plane."""

    weights = KERNEL_PLANE_RESIDUE_WEIGHTS
    parameters = sp.symbols("k_w0:3")
    w0, w1, w2 = parameters
    z = (
        -sp.Rational(4, 3) * w0
        - sp.Rational(2648, 189) * w1
        + sp.Rational(580504, 735) * w2,
        w0,
        -2 * w1 + sp.Rational(9838, 105) * w2,
        w1,
        w2,
    )
    first_pair = symbolic_pair_combination(residual_pairs, z, QQ)
    second_u = symbolic_pair_combination(
        quadratic_pairs,
        [
            z[left] * z[right]
            for left, right in quadratic_indices
        ],
        QQ,
    )
    third_base = symbolic_pair_combination(third_bases, z, QQ)

    assert zero_scale_constant_residue(
        S,
        T,
        second_base,
        weights,
    ) == QQ.one
    linear, diagonal, cross = second_kernel_residue_forms(
        S,
        T,
        second_base,
        second_kernel_pairs,
        weights,
    )
    assert all(value == QQ.zero for value in linear)
    assert all(value == QQ.zero for value in diagonal)
    assert all(value == QQ.zero for value in cross.values())

    s4_monomials = filtered_monomials(21, 1)
    t4_monomials = filtered_monomials(17, 0)
    correction_columns = [
        poisson({monomial: QQ.one}, T)
        for monomial in s4_monomials
    ]
    correction_columns += [
        poisson(S, {monomial: QQ.one})
        for monomial in t4_monomials
    ]
    assert all(
        rational_residue(column, weights, QQ) == QQ.zero
        for column in correction_columns
    )

    symbolic_second_kernels = [
        symbolic_pair(pair, QQ)
        for pair in second_kernel_pairs
    ]
    third_right_sides = [
        scale(coupling(first_pair, second_u), -sp.S.One)
    ]
    third_right_sides.extend(
        scale(coupling(first_pair, pair), -sp.S.One)
        for pair in symbolic_second_kernels
    )
    third_vectors = solve_parametric_many(
        third_columns,
        third_right_sides,
        parameters,
        QQ,
    )
    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)
    third_u = symbolic_split(
        third_vectors[0],
        s3_monomials,
        t3_monomials,
    )
    zero_u_third_pairs = [
        symbolic_split(vector, s3_monomials, t3_monomials)
        for vector in third_vectors[1:]
    ]

    diagonal_u = add(
        poisson(second_u[0], second_u[1]),
        coupling(first_pair, third_u),
    )
    assert_zero(
        symbolic_residue(diagonal_u, weights),
        parameters,
        QQ,
    )

    symbolic_base = symbolic_pair(second_base, QQ)
    linear_u = add(
        poisson(symbolic_base[0], second_u[1]),
        poisson(second_u[0], symbolic_base[1]),
    )
    cubic_u = add(
        pi_power(S, second_u[1], 3),
        pi_power(second_u[0], T, 3),
    )
    linear_u = add(linear_u, cubic_u, sp.Rational(1, 24))
    linear_u = add(linear_u, coupling(first_pair, third_base))
    linear_u = add(
        linear_u,
        pi_power(first_pair[0], first_pair[1], 3),
        sp.Rational(1, 24),
    )
    assert_zero(
        symbolic_residue(linear_u, weights),
        parameters,
        QQ,
    )

    for second_pair, third_pair in zip(
        symbolic_second_kernels,
        zero_u_third_pairs,
    ):
        cross_value = add(
            poisson(second_u[0], second_pair[1]),
            poisson(second_pair[0], second_u[1]),
        )
        cross_value = add(
            cross_value,
            coupling(first_pair, third_pair),
        )
        assert_zero(
            symbolic_residue(cross_value, weights),
            parameters,
            QQ,
        )

    for third_kernel_pair in third_kernel_pairs:
        assert_zero(
            symbolic_residue(
                coupling(
                    first_pair,
                    symbolic_pair(third_kernel_pair, QQ),
                ),
                weights,
            ),
            parameters,
            QQ,
        )
    return {
        "field": "QQ",
        "parameters": len(parameters),
        "zero_u_second_kernel": len(second_kernel_pairs),
        "third_gauge_kernel": len(third_kernel_pairs),
        "residue_support": len(weights),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--component",
        choices=("all", "hyperplane", "quadric", "kernel-plane"),
        default="all",
    )
    arguments = parser.parse_args()
    residual_couplings = verify_uniform_correction()
    _, _, project = project_to_third_cokernel()
    S, T, _, pairs, _, _ = essential_problem()
    residual_pairs = [
        linear_combination(pairs, vector, QQ)
        for vector in RESIDUAL_VECTORS
    ]
    second_base, third_pairs = global_zero_scale_base(
        S,
        T,
        residual_pairs,
    )
    assert len(third_pairs) == 5
    quadratic_indices, quadratic_pairs = residual_quadratic_pairs(
        S,
        T,
        residual_pairs,
    )
    constant_residue = zero_scale_constant_residue(S, T, second_base)
    _, kernel_pairs, third_columns = lower_lift_data(S, T)
    _, third_kernel_vectors, third_rank = solve_affine(
        third_columns,
        {},
        QQ,
    )
    assert (third_rank, len(third_kernel_vectors)) == (1034, 31)
    s3_monomials = filtered_monomials(23, 2)
    t3_monomials = filtered_monomials(19, 1)
    third_kernel_pairs = [
        split_correction(
            vector,
            s3_monomials,
            t3_monomials,
        )
        for vector in third_kernel_vectors
    ]
    linear, diagonal, cross = second_kernel_residue_forms(
        S,
        T,
        second_base,
        kernel_pairs,
    )
    generic_kernel, a_kernel = verify_zero_scale_kernel_residue(
        residual_couplings,
        linear,
        diagonal,
        cross,
    )
    hyperplane_profile = None
    if arguments.component in ("all", "hyperplane"):
        hyperplane_profile = verify_scale_component(
            name="h",
            field=QQ,
            a_value=QQ.one,
            b_value=QQ.zero,
            S=S,
            T=T,
            residual_pairs=residual_pairs,
            residual_couplings=residual_couplings,
            quadratic_indices=quadratic_indices,
            quadratic_pairs=quadratic_pairs,
            second_base=second_base,
            third_bases=third_pairs,
            second_kernel_pairs=kernel_pairs,
            third_kernel_pairs=third_kernel_pairs,
            third_columns=third_columns,
            project=project,
        )
    alpha = sp.sqrt(-2)
    quadratic_field = QQ.algebraic_field(alpha)
    quadric_profile = None
    if arguments.component in ("all", "quadric"):
        quadric_profile = verify_scale_component(
            name="q",
            field=quadratic_field,
            a_value=quadratic_field.from_sympy(
                sp.Rational(2048, 105) + sp.Rational(2, 9) * alpha
            ),
            b_value=quadratic_field.one,
            S=S,
            T=T,
            residual_pairs=residual_pairs,
            residual_couplings=residual_couplings,
            quadratic_indices=quadratic_indices,
            quadratic_pairs=quadratic_pairs,
            second_base=second_base,
            third_bases=third_pairs,
            second_kernel_pairs=kernel_pairs,
            third_kernel_pairs=third_kernel_pairs,
            third_columns=third_columns,
            project=project,
        )
    kernel_plane_profile = None
    if arguments.component in ("all", "kernel-plane"):
        kernel_plane_profile = verify_kernel_plane(
            S=S,
            T=T,
            residual_pairs=residual_pairs,
            quadratic_indices=quadratic_indices,
            quadratic_pairs=quadratic_pairs,
            second_base=second_base,
            third_bases=third_pairs,
            second_kernel_pairs=kernel_pairs,
            third_kernel_pairs=third_kernel_pairs,
            third_columns=third_columns,
        )
    profile = residue_pairing_profile(S, T, residual_pairs)

    print(
        "PASS: one two-coordinate second correction solves the projected "
        "third-order affine equation on the entire residual five-space"
    )
    print(
        "uniform correction="
        f"{sorted(UNIFORM_KERNEL_CORRECTION.items())}"
    )
    print(
        "PASS: the global zero-scale fourth defect has fixed residue "
        f"{constant_residue}"
    )
    print(
        "PASS: the residue kills the linear and quadratic fourth-order "
        "variation on all "
        f"{generic_kernel} generic-pencil and {a_kernel} A-chart "
        "second-kernel directions"
    )
    if arguments.component == "all":
        print(
            "PASS: the unique nonzero-scale lower direction and all its "
            "cross terms have zero residue identically on the hyperplane "
            "and both geometric quadric components"
        )
    else:
        print(
            "PASS: the unique nonzero-scale lower direction and all its "
            f"cross terms have zero residue on the {arguments.component}"
        )
    print(
        f"hyperplane profile={hyperplane_profile}; "
        f"quadric profile={quadric_profile}; "
        f"kernel-plane profile={kernel_plane_profile}"
    )
    print(f"residue pairing profile={profile}")


if __name__ == "__main__":
    main()
