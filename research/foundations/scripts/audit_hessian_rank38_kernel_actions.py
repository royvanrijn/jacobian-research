#!/usr/bin/env python3
"""Audit nonlinear kernel actions of the 42-variable rank-38 HN quartic.

For the stored 21-variable cubic source ``V=x+H``, put

    p(x,y) = y.H(x).

The Hessian of ``p`` has generic corank four.  Only one kernel direction is
constant, but three explicit kernel syzygies are vertical vector fields

    D_b = sum_j b_j(x) d/dy_j

with ``b.H=0``.  The fourth is mixed in ``x,y``.  All four are commuting
locally nilpotent derivations that fix ``p``.  This script verifies the
identities exactly over Q, constructs a rational collision lift where their
action has full rank four, and records why this still is not a global affine
coordinate deletion.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

from audit_fixed_rank_hessian_witness import (
    cotangent_hessian,
    decode_h,
    derivative,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "essential_bcw_21_counterexample.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hessian_rank38_kernel_actions.json"
)

Exponent = tuple[int, ...]
SparsePoly = dict[Exponent, Fraction]


def add(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    answer = dict(left)
    for exponent, coefficient in right.items():
        value = answer.get(exponent, Fraction(0)) + coefficient
        if value:
            answer[exponent] = value
        else:
            answer.pop(exponent, None)
    return answer


def scale(poly: SparsePoly, scalar: Fraction) -> SparsePoly:
    return {
        exponent: scalar * coefficient
        for exponent, coefficient in poly.items()
        if scalar * coefficient
    }


def multiply(left: SparsePoly, right: SparsePoly) -> SparsePoly:
    answer: SparsePoly = {}
    for alpha, coefficient_alpha in left.items():
        for beta, coefficient_beta in right.items():
            exponent = tuple(a + b for a, b in zip(alpha, beta))
            answer = add(
                answer,
                {exponent: coefficient_alpha * coefficient_beta},
            )
    return answer


def monomial(dimension: int, *powers: tuple[int, int]) -> SparsePoly:
    exponent = [0] * dimension
    for variable, power in powers:
        exponent[variable] = power
    return {tuple(exponent): Fraction(1)}


def sum_polynomials(polynomials) -> SparsePoly:
    answer: SparsePoly = {}
    for poly in polynomials:
        answer = add(answer, poly)
    return answer


def matrix_vector_product(
    matrix: list[list[SparsePoly]], vector: list[SparsePoly]
) -> list[SparsePoly]:
    return [
        sum_polynomials(
            multiply(entry, coefficient)
            for entry, coefficient in zip(row, vector)
        )
        for row in matrix
    ]


def lift_x(poly: SparsePoly, dimension: int) -> SparsePoly:
    zero = (0,) * dimension
    return {
        exponent + zero: coefficient for exponent, coefficient in poly.items()
    }


def directional_derivative_of_p(
    h: list[SparsePoly], vector: list[SparsePoly]
) -> SparsePoly:
    dimension = len(h)
    total_dimension = 2 * dimension
    gradient_x: list[SparsePoly] = []
    for variable in range(dimension):
        gradient_x.append(
            sum_polynomials(
                multiply(
                    lift_x(derivative(h[output], variable), dimension),
                    monomial(total_dimension, (dimension + output, 1)),
                )
                for output in range(dimension)
            )
        )
    gradient_y = [lift_x(component, dimension) for component in h]
    return sum_polynomials(
        multiply(coefficient, gradient)
        for coefficient, gradient in zip(vector, gradient_x + gradient_y)
    )


def apply_derivation(poly: SparsePoly, vector: list[SparsePoly]) -> SparsePoly:
    return sum_polynomials(
        multiply(derivative(poly, variable), coefficient)
        for variable, coefficient in enumerate(vector)
    )


def bracket(
    left: list[SparsePoly], right: list[SparsePoly]
) -> list[SparsePoly]:
    return [
        add(
            apply_derivation(right_component, left),
            scale(apply_derivation(left_component, right), Fraction(-1)),
        )
        for left_component, right_component in zip(left, right)
    ]


def set_term(
    vector: list[SparsePoly],
    component: int,
    coefficient: int,
    *powers: tuple[int, int],
) -> None:
    vector[component] = add(
        vector[component],
        scale(monomial(len(vector), *powers), Fraction(coefficient)),
    )


def evaluate(poly: SparsePoly, point: list[Fraction]) -> Fraction:
    answer = Fraction(0)
    for exponent, coefficient in poly.items():
        value = coefficient
        for coordinate, power in zip(point, exponent):
            value *= coordinate**power
        answer += value
    return answer


def rank_over_q(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (
                row
                for row in range(pivot_row, len(work))
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [value / pivot_value for value in work[pivot_row]]
        for row in range(len(work)):
            if row != pivot_row and work[row][column]:
                factor = work[row][column]
                work[row] = [
                    left - factor * right
                    for left, right in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def solve_over_q(
    matrix: list[list[Fraction]], right_hand_side: list[Fraction]
) -> list[Fraction]:
    dimension = len(matrix)
    work = [
        row[:] + [right_hand_side[index]]
        for index, row in enumerate(matrix)
    ]
    for column in range(dimension):
        pivot = next(
            row
            for row in range(column, dimension)
            if work[row][column]
        )
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [value / pivot_value for value in work[column]]
        for row in range(dimension):
            if row != column and work[row][column]:
                factor = work[row][column]
                work[row] = [
                    left - factor * right
                    for left, right in zip(work[row], work[column])
                ]
    return [row[-1] for row in work]


def encode_point(point: list[Fraction]) -> list[str]:
    return [str(value) for value in point]


def encode(poly: SparsePoly) -> list[dict[str, object]]:
    return [
        {
            "coefficient": str(coefficient),
            "monomial": [
                [variable, power]
                for variable, power in enumerate(exponent)
                if power
            ],
        }
        for exponent, coefficient in sorted(poly.items())
    ]


def main() -> None:
    stored = json.loads(SOURCE.read_text())
    dimension = int(stored["dimension"])
    assert dimension == 21
    h = decode_h(stored)
    _jacobian, _upper_left, hessian = cotangent_hessian(h)
    total_dimension = 2 * dimension

    constant = [{} for _ in range(total_dimension)]
    set_term(constant, dimension + 20, 1)

    quadratic = [{} for _ in range(total_dimension)]
    set_term(quadratic, dimension + 7, 1, (0, 2))
    set_term(quadratic, dimension + 8, -1, (0, 1), (3, 1))
    set_term(quadratic, dimension + 9, 1, (1, 1), (3, 1))
    set_term(quadratic, dimension + 10, -1, (0, 1), (1, 1))

    vertical_cubic = [{} for _ in range(total_dimension)]
    set_term(vertical_cubic, dimension + 6, 1, (0, 1), (1, 2))
    set_term(
        vertical_cubic,
        dimension + 8,
        -2,
        (0, 1),
        (1, 1),
        (2, 1),
    )
    set_term(vertical_cubic, dimension + 9, 1, (1, 2), (2, 1))
    set_term(vertical_cubic, dimension + 19, 1, (1, 2), (20, 1))
    set_term(
        vertical_cubic,
        dimension + 17,
        -1,
        (0, 1),
        (2, 1),
        (20, 1),
    )

    mixed_cubic = [{} for _ in range(total_dimension)]
    set_term(mixed_cubic, 16, 504, (1, 1), (6, 1), (22, 1))
    set_term(mixed_cubic, 16, -504, (2, 1), (8, 1), (22, 1))
    set_term(mixed_cubic, 12, 63, (1, 1), (20, 1), (22, 1))
    set_term(mixed_cubic, 13, -18, (2, 1), (20, 1), (22, 1))
    set_term(mixed_cubic, 15, -504, (1, 1), (6, 1), (23, 1))
    set_term(mixed_cubic, 15, 504, (2, 1), (8, 1), (23, 1))
    set_term(mixed_cubic, 10, 42, (1, 1), (20, 1), (23, 1))
    set_term(mixed_cubic, 11, -14, (2, 1), (20, 1), (23, 1))
    set_term(mixed_cubic, 27, 504, (1, 1), (23, 1), (36, 1))
    set_term(mixed_cubic, 29, -504, (2, 1), (23, 1), (36, 1))
    set_term(mixed_cubic, 27, -504, (1, 1), (22, 1), (37, 1))
    set_term(mixed_cubic, 29, 504, (2, 1), (22, 1), (37, 1))

    syzygies = {
        "constant": constant,
        "quadratic": quadratic,
        "mixed_cubic": mixed_cubic,
        "vertical_cubic": vertical_cubic,
    }
    for name, vector in syzygies.items():
        assert all(
            not entry for entry in matrix_vector_product(hessian, vector)
        ), name

    actions = (constant, quadratic, mixed_cubic, vertical_cubic)
    assert all(not directional_derivative_of_p(h, vector) for vector in actions)

    # Every coefficient of each derivation lies in its own kernel.  Thus each
    # derivation is a translation over a fixed coefficient ring and is
    # locally nilpotent.  Pairwise bracket checks prove commutativity.
    assert all(
        not apply_derivation(component, derivation)
        for derivation in actions
        for vector in actions
        for component in vector
    )
    assert all(
        all(not component for component in bracket(left, right))
        for left_index, left in enumerate(actions)
        for right in actions[left_index + 1 :]
    )

    zero_cotangent_collision_ranks = []
    for stored_point in stored["collision_points"]:
        x_point = [Fraction(value) for value in stored_point]
        full_point = x_point + [Fraction(0)] * dimension
        coefficient_matrix = [
            [
                evaluate(vector[output], full_point)
                for vector in actions
            ]
            for output in range(total_dimension)
        ]
        zero_cotangent_collision_ranks.append(rank_over_q(coefficient_matrix))
    assert zero_cotangent_collision_ranks == [1, 3, 3]

    # Lift the second and third source collision points through the gradient
    # map (x,y) -> (x+H(x), y+JH(x)^t y).  Taking the common second output to
    # be e_0 gives small rational y-values and puts both lifted points in the
    # full-rank locus of the four actions.
    jacobian_small = [
        [derivative(h[output], variable) for variable in range(dimension)]
        for output in range(dimension)
    ]
    source_points = [
        [Fraction(value) for value in point]
        for point in stored["collision_points"]
    ]
    eta = [Fraction(0)] * dimension
    eta[0] = Fraction(1)
    lifted_points: list[list[Fraction]] = []
    lifted_action_ranks: list[int] = []
    lifted_gradient_images: list[list[Fraction]] = []
    for x_point in source_points[1:]:
        jacobian_value = [
            [evaluate(entry, x_point) for entry in row]
            for row in jacobian_small
        ]
        second_block = [
            [
                Fraction(row == column) + jacobian_value[column][row]
                for column in range(dimension)
            ]
            for row in range(dimension)
        ]
        y_point = solve_over_q(second_block, eta)
        full_point = x_point + y_point
        lifted_points.append(full_point)
        action_matrix = [
            [evaluate(vector[output], full_point) for vector in actions]
            for output in range(total_dimension)
        ]
        lifted_action_ranks.append(rank_over_q(action_matrix))
        first_image = [
            coordinate + evaluate(component, x_point)
            for coordinate, component in zip(x_point, h)
        ]
        second_image = [
            y_point[row]
            + sum(
                jacobian_value[column][row] * y_point[column]
                for column in range(dimension)
            )
            for row in range(dimension)
        ]
        lifted_gradient_images.append(first_image + second_image)
    assert lifted_action_ranks == [4, 4]
    assert lifted_gradient_images[0] == lifted_gradient_images[1]
    assert lifted_points[0] != lifted_points[1]

    # Rows x_10, y_20, y_7, y_6 give the nonzero determinant
    # 42*x_0^3*x_1^3*x_20*y_2.
    generic_minor = scale(
        monomial(
            total_dimension,
            (0, 3),
            (1, 3),
            (20, 1),
            (dimension + 2, 1),
        ),
        Fraction(42),
    )
    assert all(evaluate(generic_minor, point) for point in lifted_points)

    payload = {
        "format": "hessian-rank38-kernel-actions-v1",
        "source": SOURCE.name,
        "cotangent_dimension": total_dimension,
        "certified_generic_hessian_rank": 38,
        "generic_hessian_corank": 4,
        "syzygy_degrees": [0, 2, 3, 3],
        "p_preserving_kernel_actions": 4,
        "vertical_actions": 3,
        "mixed_actions": 1,
        "action_properties": [
            "commuting",
            "locally nilpotent",
            "fix p=y.H exactly",
            "generic rank four",
        ],
        "generic_rank_minor": encode(generic_minor),
        "zero_cotangent_ranks_at_stored_collision_points": (
            zero_cotangent_collision_ranks
        ),
        "full_rank_collision_lift": {
            "source_collision_indices": [1, 2],
            "common_second_gradient_output": encode_point(eta),
            "points": [encode_point(point) for point in lifted_points],
            "action_ranks": lifted_action_ranks,
            "common_gradient_image": encode_point(lifted_gradient_images[0]),
        },
        "consequence": (
            "The rank-38 quartic has a generically four-dimensional "
            "nonlinear Ga-action preserving p, and a rational collision pair "
            "lies in its full-rank locus. The action is nevertheless not "
            "globally free and its nonlinear quotient is not a linear "
            "orthogonal coordinate deletion."
        ),
        "status": (
            "exact structural computation; no 38- or 39-variable homogeneous "
            "HN counterexample is claimed"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

    print("PASS rank-38 kernel actions: four exact Hessian syzygies")
    print("PASS rank-38 kernel actions: four commuting LNDs fix p")
    print("PASS rank-38 kernel actions: generic action rank is four")
    print("PASS rank-38 kernel actions: zero-cotangent ranks are [1, 3, 3]")
    print("PASS rank-38 kernel actions: exact collision lift has ranks [4, 4]")
    print(f"PASS rank-38 kernel actions: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
