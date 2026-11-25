#!/usr/bin/env python3
"""Verify a five-meridian lower bound for T(5,19;2,1).

The verifier uses the literal Artin closure relations and arithmetic over
F_37. The rank-one strand colors are involutions, so the same lower bound
holds for the meridional involution rank of the order-two knot orbifold.
"""

import json

PRIME = 37
SIZE = 5
IDENTITY = tuple(tuple(int(i == j) for j in range(SIZE)) for i in range(SIZE))

AXIS_PRODUCT = (
    (36, 0, 0, 0, 0),
    (0, 17, 29, 0, 0),
    (0, 4, 33, 0, 0),
    (0, 0, 0, 33, 13),
    (0, 0, 0, 12, 7),
)
MERIDIAN_COLUMN = (1, 1, 0, 1, 0)
MERIDIAN_ROW = (5, 13, 26, 17, 34)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def matrix_product(left, right):
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) % PRIME for column in columns)
        for row in left
    )


def matrix_difference(left, right):
    return tuple(
        tuple((a - b) % PRIME for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def eliminate(matrix):
    """Reduce the first five columns and record their rank and determinant."""
    rows = [list(row) for row in matrix]
    next_row = 0
    determinant = 1
    for column in range(SIZE):
        pivot = next((i for i in range(next_row, SIZE) if rows[i][column]), None)
        if pivot is None:
            continue
        if pivot != next_row:
            rows[pivot], rows[next_row] = rows[next_row], rows[pivot]
            determinant = -determinant
        pivot_value = rows[next_row][column]
        determinant = determinant * pivot_value % PRIME
        inverse_pivot = pow(pivot_value, -1, PRIME)
        rows[next_row] = [a * inverse_pivot % PRIME for a in rows[next_row]]
        for i, row in enumerate(rows):
            if i != next_row and row[column]:
                factor = row[column]
                rows[i] = [(a - factor * b) % PRIME for a, b in zip(row, rows[next_row])]
        next_row += 1
    return rows, next_row, determinant % PRIME if next_row == SIZE else 0


def matrix_inverse(matrix):
    rows, rank, _ = eliminate(tuple(row + unit for row, unit in zip(matrix, IDENTITY)))
    require(rank == SIZE, "matrix is singular")
    return tuple(tuple(row[SIZE:]) for row in rows)


def conjugate(conjugator, matrix, inverse_conjugator):
    return matrix_product(matrix_product(conjugator, matrix), inverse_conjugator)


def braid_action(involutions, braid):
    """Apply positive Hurwitz moves (U,V) -> (UVU^-1,U)."""
    bottom = list(involutions)
    for crossing in braid:
        i = crossing - 1
        left, right = bottom[i:i + 2]
        bottom[i], bottom[i + 1] = conjugate(left, right, left), left
    return tuple(bottom)


def verify():
    first = tuple(
        tuple((unit + u * v) % PRIME for unit, v in zip(row, MERIDIAN_ROW))
        for row, u in zip(IDENTITY, MERIDIAN_COLUMN)
    )
    require(matrix_product(first, first) == IDENTITY, "first color is not an involution")
    axis_square = matrix_product(AXIS_PRODUCT, AXIS_PRODUCT)
    axis_fourth_power = matrix_product(axis_square, axis_square)
    inverse_axis_fourth_power = matrix_inverse(axis_fourth_power)
    shift = matrix_product(first, axis_fourth_power)
    inverse_shift = matrix_product(inverse_axis_fourth_power, first)
    second = conjugate(shift, first, inverse_shift)
    third = conjugate(shift, second, inverse_shift)
    fourth = conjugate(shift, third, inverse_shift)
    fifth = conjugate(inverse_axis_fourth_power, first, axis_fourth_power)
    meridians = (first, second, third, fourth, fifth)

    product = IDENTITY
    for meridian in meridians:
        require(matrix_product(meridian, meridian) == IDENTITY,
                "a strand color is not an involution")
        product = matrix_product(product, meridian)
    require(product == AXIS_PRODUCT, "strand meridian product differs from the axis input")

    torus_block = (1, 2, 3, 4) * 19
    full_twist = (1, 1)
    bottom = braid_action(meridians, torus_block)
    require(bottom == meridians, "the torus block changes the coloring")
    require(braid_action(bottom, full_twist) == meridians, "a literal braid closure relation fails")

    adjacent_squares = (1, 2, 3, 4)
    for i in adjacent_squares:
        left, right = meridians[i - 1:i + 1]
        require(matrix_product(left, right) == matrix_product(right, left),
                f"adjacent strand colors {i} and {i + 1} do not commute")
        require(braid_action(meridians, (i, i)) == meridians,
                f"the square of braid generator {i} changes the coloring")

    displacements = tuple(matrix_difference(meridian, IDENTITY) for meridian in meridians)
    ranks = tuple(eliminate(displacement)[1] for displacement in displacements)
    require(ranks == (1,) * SIZE, "a meridian does not have displacement rank one")
    first_rows = tuple(displacement[0] for displacement in displacements)
    _, _, row_determinant = eliminate(first_rows)
    require(row_determinant == 31, "first displacement rows have the wrong determinant")
    columns = tuple(
        next(column for column in zip(*displacement) if any(column))
        for displacement in displacements
    )
    _, span_dimension, _ = eliminate(columns)
    require(span_dimension == SIZE, "the meridian displacement lines do not span")

    return {
        "knot": "T(5,19;2,1)",
        "field": "F_37",
        "braid_crossings": len(torus_block) + len(full_twist),
        "closure_relations_verified": len(meridians),
        "torus_block_fixes_coloring": True,
        "adjacent_braid_generator_squares_fix_coloring": adjacent_squares,
        "meridian_images_are_involutions": True,
        "meridian_displacement_ranks": ranks,
        "first_displacement_rows_determinant": row_determinant,
        "dimension_of_meridian_displacement_span": span_dimension,
        "meridional_rank_lower_bound": SIZE,
        "orbifold_meridional_involution_rank_lower_bound": SIZE,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
