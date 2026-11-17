#!/usr/bin/env python3
"""Lift generic normal-cone directions through the classical quadrics.

At a generic point of L1 or L2 the reduced normal cone has two branches over
Q(sqrt(-5)).  Since the original Maurer--Cartan equations are quadratic, a
formal arc is constructed recursively by one fixed linear solve at every
order.  Higher coefficients are allowed to move along all tangent
directions, including the seven directions tangent to the coordinate
stratum.
"""

from __future__ import annotations

import argparse

import sympy as sp
from sympy.polys.domains import QQ

from explore_degree_five_quantum_residue import column_rank, solve_affine
from explore_rank_two_odd_normal_cones import (
    COORDINATE_SPACES,
    RADICAL_COORDINATES,
    jacobian_rows,
    normal_cone,
    nullspace,
    reduced_quadrics,
)


def extend_vector(vector, field):
    return {
        index: field.convert(value)
        for index, value in vector.items()
        if value
    }


def combine_vectors(vectors, coefficients, field):
    result = {}
    for vector, coefficient in zip(vectors, coefficients):
        if not coefficient:
            continue
        for index, value in vector.items():
            result[index] = (
                result.get(index, field.zero)
                + coefficient * value
            )
            if not result[index]:
                result.pop(index)
    return result


def quadratic_value(quadrics, variables, vector, field):
    coordinates = {
        variable: vector.get(index, field.zero)
        for index, variable in enumerate(variables)
    }
    result = {}
    for row, equation in enumerate(quadrics):
        value = sum(
            (
                coefficient
                * coordinates.get(left, field.zero)
                * coordinates.get(right, field.zero)
                for (left, right), coefficient in equation.items()
            ),
            field.zero,
        )
        if value:
            result[row] = value
    return result


def cross_value(quadrics, variables, left, right, field):
    left_coordinates = {
        variable: left.get(index, field.zero)
        for index, variable in enumerate(variables)
    }
    right_coordinates = {
        variable: right.get(index, field.zero)
        for index, variable in enumerate(variables)
    }
    result = {}
    for row, equation in enumerate(quadrics):
        value = field.zero
        for (first, second), coefficient in equation.items():
            if first == second:
                value += (
                    field(2)
                    * coefficient
                    * left_coordinates.get(first, field.zero)
                    * right_coordinates.get(first, field.zero)
                )
            else:
                value += coefficient * (
                    left_coordinates.get(first, field.zero)
                    * right_coordinates.get(second, field.zero)
                    + right_coordinates.get(first, field.zero)
                    * left_coordinates.get(second, field.zero)
                )
        if value:
            result[row] = value
    return result


def negate(vector, field):
    return {index: -value for index, value in vector.items() if value}


def coefficient(series, order, quadrics, variables, jacobian_columns, field):
    result = {}
    current = series.get(order, {})
    for column, value in current.items():
        for row, coefficient in jacobian_columns[column].items():
            result[row] = (
                result.get(row, field.zero)
                + value * coefficient
            )
    for left in range(1, (order // 2) + 1):
        right = order - left
        if left > right:
            continue
        if left == right:
            contribution = quadratic_value(
                quadrics,
                variables,
                series[left],
                field,
            )
        else:
            contribution = cross_value(
                quadrics,
                variables,
                series[left],
                series[right],
                field,
            )
        for row, value in contribution.items():
            result[row] = result.get(row, field.zero) + value
    return {
        row: value for row, value in result.items() if value
    }


def generic_normal_direction(space, data, field, alpha):
    coefficients = [field(index + 1) for index in range(len(data["tangent_basis"]))]
    root = -field(2) + alpha / field(15)
    if space == "L1":
        coefficients[10] = field.one
        coefficients[8] = root
    else:
        coefficients[11] = field(3)
        coefficients[10] = -field(10)
        coefficients[15] = field.one
        coefficients[13] = root
    return combine_vectors(
        [extend_vector(vector, field) for vector in data["tangent_basis"]],
        coefficients,
        field,
    )


def lift_branch(space, maximum_order):
    alpha_expression = sp.sqrt(-5)
    field = QQ.algebraic_field(alpha_expression)
    alpha = field.from_sympy(alpha_expression)

    rational_data = normal_cone(space, QQ)
    quadrics = [
        {
            monomial: field.convert(coefficient)
            for monomial, coefficient in equation.items()
        }
        for equation in rational_data["quadrics"]
    ]
    active = rational_data["active_variables"]
    coordinate_space = COORDINATE_SPACES[space]
    base_point = {
        variable: field(position + 1)
        for position, variable in enumerate(coordinate_space)
    }
    full_jacobian = jacobian_rows(
        quadrics,
        active,
        base_point,
        field,
    )
    jacobian_rank, tangent_kernel = nullspace(
        full_jacobian,
        len(active),
        field,
    )
    jacobian_columns = [dict() for _ in active]
    for row, entries in full_jacobian.items():
        for column, value in entries.items():
            jacobian_columns[column][row] = value

    normal_direction = generic_normal_direction(
        space,
        rational_data,
        field,
        alpha,
    )
    active_index = {
        variable: index for index, variable in enumerate(active)
    }
    first = {
        active_index[
            rational_data["normal_variable_indices"][normal_index]
        ]: value
        for normal_index, value in normal_direction.items()
    }
    assert not coefficient(
        {1: first},
        1,
        quadrics,
        active,
        jacobian_columns,
        field,
    )

    series = {1: first}
    second_rhs = negate(
        quadratic_value(quadrics, active, first, field),
        field,
    )
    second, _, rank = solve_affine(
        jacobian_columns,
        second_rhs,
        field,
    )
    assert rank == jacobian_rank
    series[2] = second
    assert not coefficient(
        series,
        2,
        quadrics,
        active,
        jacobian_columns,
        field,
    )

    tangent_cross_columns = [
        cross_value(quadrics, active, first, vector, field)
        for vector in tangent_kernel
    ]
    continuation_columns = jacobian_columns + tangent_cross_columns
    continuation_rank = None
    steps = [(1, jacobian_rank), (2, jacobian_rank)]
    for order in range(3, maximum_order + 1):
        residual = coefficient(
            series,
            order,
            quadrics,
            active,
            jacobian_columns,
            field,
        )
        try:
            solution, _, continuation_rank = solve_affine(
                continuation_columns,
                negate(residual, field),
                field,
            )
        except ValueError:
            span_rank = column_rank(continuation_columns)
            augmented_rank = column_rank(
                continuation_columns + [negate(residual, field)]
            )
            return {
                "space": space,
                "field": str(field),
                "jacobian_rank": jacobian_rank,
                "tangent_dimension": len(tangent_kernel),
                "continuation_rank": span_rank,
                "obstruction_order": order,
                "augmented_rank": augmented_rank,
                "steps": steps,
                "support": {
                    coefficient_order: len(vector)
                    for coefficient_order, vector in series.items()
                },
            }
        current = {
            index: value
            for index, value in solution.items()
            if index < len(active)
        }
        correction = {
            index - len(active): value
            for index, value in solution.items()
            if index >= len(active)
        }
        series[order] = current
        series[order - 1] = combine_vectors(
            (series[order - 1], *tangent_kernel),
            (
                field.one,
                *(
                    correction.get(index, field.zero)
                    for index in range(len(tangent_kernel))
                ),
            ),
            field,
        )
        assert not coefficient(
            series,
            order,
            quadrics,
            active,
            jacobian_columns,
            field,
        )
        steps.append((order, continuation_rank))
    return {
        "space": space,
        "field": str(field),
        "jacobian_rank": jacobian_rank,
        "tangent_dimension": len(tangent_kernel),
        "continuation_rank": continuation_rank,
        "steps": steps,
        "support": {
            order: len(vector)
            for order, vector in series.items()
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("space", choices=sorted(COORDINATE_SPACES))
    parser.add_argument("--order", type=int, default=8)
    arguments = parser.parse_args()
    profile = lift_branch(arguments.space, arguments.order)
    if "obstruction_order" in profile:
        print(
            "PASS: the generic reduced normal direction has a genuine "
            "higher Kuranishi obstruction at order "
            f"{profile['obstruction_order']}"
        )
    else:
        print(
            "PASS: generic reduced normal direction formally lifts through "
            f"order {arguments.order}"
        )
    print(profile)


if __name__ == "__main__":
    main()
