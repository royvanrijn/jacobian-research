#!/usr/bin/env python3
"""Low-degree SL2 invariants missed by the SIC2C4 moment sequence.

The coefficient matrix C represents an element of

    End(Sym^4) = Sym^0 + Sym^2 + Sym^4 + Sym^6 + Sym^8.

Writing A=C^T D, with D_ii=i!(4-i)!, this checker constructs the five
Casimir projectors, the corresponding quadratic trace forms, and the ten
primitive cubic invariants.  It also verifies a full-rank modular Jacobian
certificate for the first 22 moments.
"""

from __future__ import annotations

import itertools
import json
from collections import defaultdict
from math import factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_counterexample_missing_invariant.json"
)
PRIME = 1_000_003
MAX_MOMENT = 22
COMPONENTS = tuple(range(5))
CUBIC_TRIPLES = (
    (1, 1, 2),
    (1, 2, 3),
    (1, 3, 4),
    (2, 2, 2),
    (2, 2, 4),
    (2, 3, 3),
    (2, 3, 4),
    (2, 4, 4),
    (3, 3, 4),
    (4, 4, 4),
)


def coefficient_matrix() -> sp.Matrix:
    return sp.Matrix(
        [
            [-1, 2, 0, 0, 0],
            [-sp.Rational(3, 2), 2, 6, 0, 0],
            [-sp.Rational(1, 2), sp.Rational(3, 2), 6, 6, 0],
            [0, 1, sp.Rational(3, 2), 2, 2],
            [0, 0, -sp.Rational(1, 2), -sp.Rational(3, 2), -1],
        ]
    )


def sl2_matrices() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    raising = sp.zeros(5)
    lowering = sp.zeros(5)
    cartan = sp.zeros(5)
    for index in range(5):
        cartan[index, index] = 2 * index - 4
        if index < 4:
            raising[index + 1, index] = 4 - index
        if index > 0:
            lowering[index - 1, index] = index
    return raising, lowering, cartan


def apolar_matrix() -> sp.Matrix:
    result = sp.zeros(5)
    for index, value in enumerate(
        (1, -sp.Rational(1, 4), sp.Rational(1, 6), -sp.Rational(1, 4), 1)
    ):
        result[index, 4 - index] = value
    raising, lowering, cartan = sl2_matrices()
    for generator in (raising, lowering, cartan):
        assert generator.T * result + result * generator == sp.zeros(5)
    assert result == result.T
    return result


def casimir_projectors() -> dict[int, sp.Matrix]:
    identity = sp.eye(5)

    def adjoint(matrix: sp.Matrix) -> sp.Matrix:
        # Column-major vectorization of XA-AX.
        return sp.kronecker_product(identity, matrix) - sp.kronecker_product(
            matrix.T, identity
        )

    raising, lowering, cartan = sl2_matrices()
    casimir = (
        adjoint(raising) * adjoint(lowering)
        + adjoint(lowering) * adjoint(raising)
        + adjoint(cartan) * adjoint(cartan) / 2
    )
    assert casimir.eigenvals() == {0: 1, 4: 3, 12: 5, 24: 7, 40: 9}

    projectors: dict[int, sp.Matrix] = {}
    for component in COMPONENTS:
        eigenvalue = 2 * component * (component + 1)
        projector = sp.eye(25)
        for other in COMPONENTS:
            if other == component:
                continue
            other_eigenvalue = 2 * other * (other + 1)
            projector *= (casimir - other_eigenvalue * sp.eye(25)) / (
                eigenvalue - other_eigenvalue
            )
        assert projector * projector == projector
        projectors[component] = projector
    assert sum(projectors.values(), sp.zeros(25)) == sp.eye(25)
    return projectors


def quadratic_matrices(
    projectors: dict[int, sp.Matrix],
) -> dict[int, sp.Matrix]:
    # Row-major C and column-major A=C^T D have the same ordering; only
    # the row factorial D_ii rescales the entries.
    coefficient_to_operator = sp.diag(
        *[
            factorial(row) * factorial(4 - row)
            for row in range(5)
            for _column in range(5)
        ]
    )
    transpose_pairing = sp.zeros(25)
    for row in range(5):
        for column in range(5):
            transpose_pairing[row + 5 * column, column + 5 * row] = 1

    result = {}
    for component, projector in projectors.items():
        result[component] = sp.simplify(
            coefficient_to_operator.T
            * projector.T
            * transpose_pairing
            * projector
            * coefficient_to_operator
        )
    return result


def second_moment_matrix() -> sp.Matrix:
    result = sp.zeros(25)
    for dual_left in range(5):
        for coordinate_left in range(5):
            left = 5 * dual_left + coordinate_left
            for dual_right in range(5):
                coordinate_right = (
                    dual_left + dual_right - coordinate_left
                )
                if not 0 <= coordinate_right < 5:
                    continue
                right = 5 * dual_right + coordinate_right
                result[left, right] += (
                    factorial(dual_left + dual_right)
                    * factorial(8 - dual_left - dual_right)
                )
    assert result == result.T
    return result


def vectorized_operator(
    coefficients: sp.Matrix,
) -> tuple[sp.Matrix, sp.Matrix]:
    factorial_diagonal = sp.diag(
        *[factorial(index) * factorial(4 - index) for index in range(5)]
    )
    operator = coefficients.T * factorial_diagonal
    vector = sp.Matrix(
        [
            operator[row, column]
            for column in range(5)
            for row in range(5)
        ]
    )
    return operator, vector


def vectorize(matrix: sp.Matrix) -> sp.Matrix:
    return sp.Matrix(
        [
            matrix[row, column]
            for column in range(5)
            for row in range(5)
        ]
    )


def apolar_adjoint(operator: sp.Matrix) -> sp.Matrix:
    apolar = apolar_matrix()
    return apolar.inv() * operator.T * apolar


def adjoint_coefficients(coefficients: sp.Matrix) -> sp.Matrix:
    factorial_diagonal = sp.diag(
        *[factorial(index) * factorial(4 - index) for index in range(5)]
    )
    operator, _ = vectorized_operator(coefficients)
    return factorial_diagonal.inv() * apolar_adjoint(operator).T


def component_matrices(
    vector: sp.Matrix,
    projectors: dict[int, sp.Matrix],
) -> dict[int, sp.Matrix]:
    result = {}
    for component, projector in projectors.items():
        projected = projector * vector
        result[component] = sp.Matrix(
            5,
            5,
            lambda row, column: projected[row + 5 * column],
        )
    return result


def cubic_invariant(
    components: dict[int, sp.Matrix],
    triple: tuple[int, int, int],
) -> sp.Rational:
    # The odd-sum all-distinct triple is alternating under reversal, so a
    # full symmetrization would cancel it.
    if triple == (2, 3, 4):
        return sp.trace(
            components[2] * components[3] * components[4]
        )
    permutations = sorted(set(itertools.permutations(triple)))
    return sp.factor(
        sum(
            sp.trace(
                components[first]
                * components[second]
                * components[third]
            )
            for first, second, third in permutations
        )
    )


def invariant_hilbert_coefficients(cutoff: int) -> list[int]:
    weights = [
        weight
        for component in COMPONENTS
        for weight in range(-2 * component, 2 * component + 1, 2)
    ]
    coefficients = [defaultdict(int) for _ in range(cutoff + 1)]
    coefficients[0][0] = 1
    for weight in weights:
        updated = [defaultdict(int) for _ in range(cutoff + 1)]
        for degree in range(cutoff + 1):
            for total_weight, multiplicity in coefficients[degree].items():
                for power in range(cutoff - degree + 1):
                    updated[degree + power][
                        total_weight + power * weight
                    ] += multiplicity
        coefficients = updated
    return [
        coefficients[degree][0] - coefficients[degree][2]
        for degree in range(cutoff + 1)
    ]


def refined_invariant_hilbert_coefficients(
    max_component: int,
    cutoff: int,
) -> tuple[list[int], list[int]]:
    coefficients = [defaultdict(int) for _ in range(cutoff + 1)]
    coefficients[0][(0, 0)] = 1
    for component in range(max_component + 1):
        parity = component % 2
        for weight in range(-2 * component, 2 * component + 1, 2):
            updated = [defaultdict(int) for _ in range(cutoff + 1)]
            for degree in range(cutoff + 1):
                for (
                    total_weight,
                    total_parity,
                ), multiplicity in coefficients[degree].items():
                    for power in range(cutoff - degree + 1):
                        updated[degree + power][
                            (
                                total_weight + power * weight,
                                total_parity ^ (parity * (power % 2)),
                            )
                        ] += multiplicity
            coefficients = updated
    even = [
        coefficients[degree][(0, 0)] - coefficients[degree][(2, 0)]
        for degree in range(cutoff + 1)
    ]
    odd = [
        coefficients[degree][(0, 1)] - coefficients[degree][(2, 1)]
        for degree in range(cutoff + 1)
    ]
    return even, odd


def hilbert_numerator(
    hilbert: list[int],
    degrees: tuple[int, ...],
) -> list[int]:
    result = list(hilbert)
    for degree in degrees:
        for index in range(len(result) - 1, degree - 1, -1):
            result[index] -= result[index - degree]
    return result


def weighted_exponents(
    weights: tuple[int, ...],
    total: int,
) -> list[tuple[int, ...]]:
    result: list[tuple[int, ...]] = []

    def recurse(
        index: int,
        remainder: int,
        prefix: tuple[int, ...],
    ) -> None:
        if index == len(weights):
            if remainder == 0:
                result.append(prefix)
            return
        weight = weights[index]
        for exponent in range(remainder // weight + 1):
            recurse(
                index + 1,
                remainder - exponent * weight,
                prefix + (exponent,),
            )

    recurse(0, total, ())
    return result


def d2_moment_values_mod(point: list[list[int]]) -> list[int]:
    factorials = [1]
    for value in range(1, 15):
        factorials.append(factorials[-1] * value % PRIME)
    power = [[1]]
    values = []
    for order in range(1, 8):
        size = 2 * order + 1
        current = [[0] * size for _ in range(size)]
        for left, previous_row in enumerate(power):
            for right, value in enumerate(previous_row):
                if not value:
                    continue
                for row in range(3):
                    for column in range(3):
                        current[left + row][right + column] = (
                            current[left + row][right + column]
                            + value * point[row][column]
                        ) % PRIME
        power = current
        values.append(
            sum(
                factorials[index]
                * factorials[2 * order - index]
                * power[index][index]
                for index in range(size)
            )
            % PRIME
        )
    return values


def moment_jacobian_rows() -> tuple[list[list[int]], list[list[int]]]:
    point = [
        [
            (17 * row + 31 * column + 7 * row * column + 11) % 97 - 48
            for column in range(5)
        ]
        for row in range(5)
    ]
    factorials = [1]
    for value in range(1, 4 * MAX_MOMENT + 1):
        factorials.append(factorials[-1] * value % PRIME)

    # powers[k][i][j] is [x^i y^j] C(x,y)^k modulo PRIME.
    powers = [[[1]]]
    for power in range(1, MAX_MOMENT):
        previous = powers[-1]
        size = 4 * power + 1
        current = [[0] * size for _ in range(size)]
        for left, previous_row in enumerate(previous):
            for right, value in enumerate(previous_row):
                if not value:
                    continue
                for row in range(5):
                    for column in range(5):
                        current[left + row][right + column] = (
                            current[left + row][right + column]
                            + value * point[row][column]
                        ) % PRIME
        powers.append(current)

    rows = []
    for order in range(1, MAX_MOMENT + 1):
        previous = powers[order - 1]
        row_values = []
        for dual_index in range(5):
            for coordinate_index in range(5):
                value = 0
                for total in range(4 * order + 1):
                    left = total - dual_index
                    right = total - coordinate_index
                    if (
                        0 <= left < len(previous)
                        and 0 <= right < len(previous)
                    ):
                        value += (
                            factorials[total]
                            * factorials[4 * order - total]
                            * previous[left][right]
                        )
                row_values.append(order * value % PRIME)
        rows.append(row_values)
    return rows, point


def moment_values_mod(point: list[list[int]]) -> list[int]:
    factorials = [1]
    for value in range(1, 4 * MAX_MOMENT + 1):
        factorials.append(factorials[-1] * value % PRIME)
    power = [[1]]
    values = []
    for order in range(1, MAX_MOMENT + 1):
        size = 4 * order + 1
        current = [[0] * size for _ in range(size)]
        for left, previous_row in enumerate(power):
            for right, value in enumerate(previous_row):
                if not value:
                    continue
                for row in range(5):
                    for column in range(5):
                        current[left + row][right + column] = (
                            current[left + row][right + column]
                            + value * point[row][column]
                        ) % PRIME
        power = current
        values.append(
            sum(
                factorials[index]
                * factorials[4 * order - index]
                * power[index][index]
                for index in range(size)
            )
            % PRIME
        )
    return values


def exact_moment_jacobian(
    coefficients: sp.Matrix,
    cutoff: int,
) -> sp.Matrix:
    powers = [[[sp.Integer(1)]]]
    for power in range(1, cutoff):
        previous = powers[-1]
        size = 4 * power + 1
        current = [[sp.Integer(0)] * size for _ in range(size)]
        for left, previous_row in enumerate(previous):
            for right, value in enumerate(previous_row):
                for row in range(5):
                    for column in range(5):
                        current[left + row][right + column] += (
                            value * coefficients[row, column]
                        )
        powers.append(current)

    rows = []
    for order in range(1, cutoff + 1):
        previous = powers[order - 1]
        row_values = []
        for dual_index in range(5):
            for coordinate_index in range(5):
                value = sp.Integer(0)
                for total in range(4 * order + 1):
                    left = total - dual_index
                    right = total - coordinate_index
                    if (
                        0 <= left < len(previous)
                        and 0 <= right < len(previous)
                    ):
                        value += (
                            factorial(total)
                            * factorial(4 * order - total)
                            * previous[left][right]
                        )
                row_values.append(order * value)
        rows.append(row_values)
    return sp.Matrix(rows)


def rank_and_pivots_mod(
    matrix: list[list[int]],
) -> tuple[int, list[int]]:
    reduced = [row[:] for row in matrix]
    row_count = len(reduced)
    column_count = len(reduced[0])
    rank = 0
    pivots = []
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(rank, row_count)
                if reduced[row][column] % PRIME
            ),
            None,
        )
        if pivot is None:
            continue
        reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
        inverse = pow(reduced[rank][column], PRIME - 2, PRIME)
        reduced[rank] = [
            value * inverse % PRIME for value in reduced[rank]
        ]
        for row in range(row_count):
            if row == rank or not reduced[row][column]:
                continue
            scale = reduced[row][column]
            reduced[row] = [
                (left - scale * right) % PRIME
                for left, right in zip(reduced[row], reduced[rank])
            ]
        pivots.append(column)
        rank += 1
        if rank == row_count:
            break
    return rank, pivots


def determinant_mod(matrix: list[list[int]]) -> int:
    reduced = [row[:] for row in matrix]
    result = 1
    for column in range(len(reduced)):
        pivot = next(
            (
                row
                for row in range(column, len(reduced))
                if reduced[row][column] % PRIME
            ),
            None,
        )
        assert pivot is not None
        if pivot != column:
            reduced[column], reduced[pivot] = (
                reduced[pivot],
                reduced[column],
            )
            result = -result
        pivot_value = reduced[column][column] % PRIME
        result = result * pivot_value % PRIME
        inverse = pow(pivot_value, PRIME - 2, PRIME)
        for row in range(column + 1, len(reduced)):
            scale = reduced[row][column] * inverse % PRIME
            for entry in range(column, len(reduced)):
                reduced[row][entry] = (
                    reduced[row][entry]
                    - scale * reduced[column][entry]
                ) % PRIME
    return result % PRIME


def main() -> None:
    coefficients = coefficient_matrix()
    operator, vector = vectorized_operator(coefficients)
    projectors = casimir_projectors()
    quadratics = quadratic_matrices(projectors)
    components = component_matrices(vector, projectors)

    q_values = {
        component: sp.factor(
            (sp.Matrix(
                [
                    coefficients[row, column]
                    for row in range(5)
                    for column in range(5)
                ]
            ).T
            * quadratic
            * sp.Matrix(
                [
                    coefficients[row, column]
                    for row in range(5)
                    for column in range(5)
                ]
            ))[0]
        )
        for component, quadratic in quadratics.items()
    }
    assert q_values == {0: 0, 1: -864, 2: 2016, 3: 0, 4: 0}
    assert sum(q_values.values()) == sp.trace(operator**2) == 1152

    moment_two = second_moment_matrix()
    moment_two_weights = (126, 84, 36, 9, 1)
    assert moment_two == sum(
        (
            moment_two_weights[component] * quadratics[component]
            for component in COMPONENTS
        ),
        sp.zeros(25),
    )
    assert (
        sum(
            moment_two_weights[component] * q_values[component]
            for component in COMPONENTS
        )
        == 0
    )

    cubic_values = {
        "".join(str(value) for value in triple): cubic_invariant(
            components, triple
        )
        for triple in CUBIC_TRIPLES
    }
    assert cubic_values == {
        "112": -sp.Rational(134784, 5),
        "123": sp.Rational(497664, 5),
        "134": 0,
        "222": -10368,
        "224": 0,
        "233": 0,
        "234": 0,
        "244": 0,
        "334": 0,
        "444": 0,
    }
    generic_coefficients = sp.Matrix(
        5,
        5,
        lambda row, column: (
            (13 * row + 29 * column + 7 * row * column + 5) % 17 - 8
        ),
    )
    _, generic_vector = vectorized_operator(generic_coefficients)
    generic_components = component_matrices(generic_vector, projectors)
    assert all(
        cubic_invariant(generic_components, triple)
        for triple in CUBIC_TRIPLES
    )

    # The apolar adjoint is equivariant and acts by (-1)^r on Sym^(2r).
    generic_operator, _ = vectorized_operator(generic_coefficients)
    generic_adjoint = apolar_adjoint(generic_operator)
    assert apolar_adjoint(generic_adjoint) == generic_operator
    adjoint_components = component_matrices(
        vectorize(generic_adjoint), projectors
    )
    for component in COMPONENTS:
        assert (
            adjoint_components[component]
            == (-1) ** component * generic_components[component]
        )
    assert cubic_invariant(
        adjoint_components, (2, 3, 4)
    ) == -cubic_invariant(generic_components, (2, 3, 4))

    transformed_coefficients = adjoint_coefficients(generic_coefficients)
    expected_transformed_coefficients = sp.Matrix(
        5,
        5,
        lambda row, column: (
            (-1) ** (row + column)
            * generic_coefficients[4 - column, 4 - row]
        ),
    )
    assert transformed_coefficients == expected_transformed_coefficients

    # At F the adjoint is a PGL2 translate, explaining why all odd
    # invariants tested at this special point vanish.
    parity_operator = sp.diag(1, -1, 1, -1, 1)
    assert apolar_adjoint(operator) == (
        parity_operator * operator * parity_operator
    )

    # Compose the adjoint with this PGL2 translate so that the resulting
    # quotient involution fixes F, then split its tangent representation.
    factorial_diagonal = sp.diag(
        *[factorial(index) * factorial(4 - index) for index in range(5)]
    )
    local_involution = sp.zeros(25)
    for basis_index in range(25):
        basis_coefficient = sp.zeros(5)
        basis_coefficient[basis_index // 5, basis_index % 5] = 1
        basis_operator, _ = vectorized_operator(basis_coefficient)
        transformed_operator = (
            parity_operator
            * apolar_adjoint(basis_operator)
            * parity_operator
        )
        transformed_coefficient = (
            factorial_diagonal.inv() * transformed_operator.T
        )
        local_involution[:, basis_index] = sp.Matrix(
            [
                transformed_coefficient[row, column]
                for row in range(5)
                for column in range(5)
            ]
        )
    assert local_involution**2 == sp.eye(25)
    plus_space = sp.Matrix.hstack(
        *(local_involution - sp.eye(25)).nullspace()
    )
    minus_space = sp.Matrix.hstack(
        *(local_involution + sp.eye(25)).nullspace()
    )
    assert (plus_space.cols, minus_space.cols) == (15, 10)

    local_jacobian = exact_moment_jacobian(coefficients, 12)
    assert local_jacobian.rank() == 12
    assert local_jacobian * local_involution == local_jacobian
    assert (
        (local_jacobian * plus_space).rank(),
        (local_jacobian * minus_space).rank(),
    ) == (12, 0)

    orbit_vectors = []
    for generator in sl2_matrices():
        tangent_operator = generator * operator - operator * generator
        tangent_coefficient = factorial_diagonal.inv() * tangent_operator.T
        orbit_vectors.append(
            sp.Matrix(
                [
                    tangent_coefficient[row, column]
                    for row in range(5)
                    for column in range(5)
                ]
            )
        )
    orbit = sp.Matrix.hstack(*orbit_vectors)
    assert orbit.rank() == 3
    orbit_plus_rank = ((sp.eye(25) + local_involution) * orbit).rank()
    orbit_minus_rank = ((sp.eye(25) - local_involution) * orbit).rank()
    assert (orbit_plus_rank, orbit_minus_rank) == (1, 2)
    quotient_tangent_eigenspaces = (
        plus_space.cols - orbit_plus_rank,
        minus_space.cols - orbit_minus_rank,
    )
    quotient_fiber_eigenspaces = (
        plus_space.cols
        - (local_jacobian * plus_space).rank()
        - orbit_plus_rank,
        minus_space.cols - orbit_minus_rank,
    )
    assert quotient_tangent_eigenspaces == (14, 8)
    assert quotient_fiber_eigenspaces == (2, 8)

    hilbert = invariant_hilbert_coefficients(5)
    assert hilbert == [1, 1, 5, 15, 65, 219]
    even_d2, odd_d2 = refined_invariant_hilbert_coefficients(2, 12)
    _, odd_d3 = refined_invariant_hilbert_coefficients(3, 4)
    _, odd_d4 = refined_invariant_hilbert_coefficients(4, 3)
    assert odd_d2[:7] == [0, 0, 0, 0, 0, 0, 1]
    assert odd_d3 == [0, 0, 0, 0, 3]
    assert odd_d4 == [0, 0, 0, 1]
    d2_hilbert = [
        even + odd for even, odd in zip(even_d2, odd_d2)
    ]
    d2_numerator = hilbert_numerator(
        d2_hilbert, tuple(range(1, 7))
    )
    assert [
        (degree, coefficient)
        for degree, coefficient in enumerate(d2_numerator)
        if coefficient
    ] == [
        (0, 1),
        (2, 1),
        (3, 1),
        (4, 1),
        (6, 2),
        (8, 1),
        (9, 1),
        (10, 1),
        (12, 1),
    ]
    assert sum(d2_numerator) == 10
    degree_seven_exponents = weighted_exponents(
        tuple(range(1, 7)), 7
    )
    assert len(degree_seven_exponents) == 14
    d2_degree_seven_evaluations = []
    for sample in range(15):
        point = [
            [
                (
                    17 * row
                    + 31 * column
                    + 7 * row * column
                    + 11
                    + sample * (13 + 5 * row + 3 * column)
                    + sample**2 * (row + 2 * column + 1)
                )
                % 101
                - 50
                for column in range(3)
            ]
            for row in range(3)
        ]
        moments = d2_moment_values_mod(point)
        degree_seven_values = []
        for exponents in degree_seven_exponents:
            value = 1
            for moment, exponent in zip(moments[:6], exponents):
                value = value * pow(moment, exponent, PRIME) % PRIME
            degree_seven_values.append(value)
        degree_seven_values.append(moments[6])
        d2_degree_seven_evaluations.append(degree_seven_values)
    d2_degree_seven_rank, _ = rank_and_pivots_mod(
        d2_degree_seven_evaluations
    )
    assert d2_degree_seven_rank == 15
    d2_degree_seven_determinant = determinant_mod(
        d2_degree_seven_evaluations
    )
    assert d2_degree_seven_determinant != 0

    jacobian, jacobian_point = moment_jacobian_rows()
    rank, pivots = rank_and_pivots_mod(jacobian)
    assert rank == 22
    minor = [[row[column] for column in pivots] for row in jacobian]
    minor_determinant = determinant_mod(minor)
    assert minor_determinant != 0
    adjoint_jacobian_point = [
        [
            (-1) ** (row + column)
            * jacobian_point[4 - column][4 - row]
            for column in range(5)
        ]
        for row in range(5)
    ]
    assert moment_values_mod(jacobian_point) == moment_values_mod(
        adjoint_jacobian_point
    )

    assert coefficients.det() == 48
    assert operator.det() == 3_981_312

    payload = {
        "representation": {
            "decomposition": ["Sym^0", "Sym^2", "Sym^4", "Sym^6", "Sym^8"],
            "hilbert_coefficients_degrees_0_through_5": hilbert,
            "minimal_generator_counts_degrees_1_through_3": [1, 4, 10],
        },
        "quadratic_invariants": {
            "definition": "q_(2r)=tr(P_(2r)(A)^2), A=C^T D",
            "values_at_F": {
                f"q_{2 * component}": str(value)
                for component, value in q_values.items()
            },
            "mu_2_weights_on_q_0_q_2_q_4_q_6_q_8": list(
                moment_two_weights
            ),
            "trace_A_squared_at_F": str(sp.trace(operator**2)),
        },
        "cubic_generators": {
            "component_index_triples": [
                list(triple) for triple in CUBIC_TRIPLES
            ],
            "values_at_F": {
                key: str(value) for key, value in cubic_values.items()
            },
        },
        "determinants": {
            "coefficient_matrix": str(coefficients.det()),
            "operator_matrix": str(operator.det()),
        },
        "moment_jacobian_certificate": {
            "prime": PRIME,
            "orders": [1, MAX_MOMENT],
            "rank": rank,
            "pivot_columns_zero_based": pivots,
            "minor_determinant_mod_prime": minor_determinant,
            "point": jacobian_point,
        },
        "apolar_adjoint": {
            "coefficient_formula": (
                "tau(C)_(ij)=(-1)^(i+j) C_(4-j,4-i)"
            ),
            "component_signs_for_Sym0_Sym2_Sym4_Sym6_Sym8": [
                1,
                -1,
                1,
                -1,
                1,
            ],
            "moments_checked_at_generic_point_through_order": MAX_MOMENT,
            "first_odd_invariant": "c_234 in degree 3",
            "c_234_generic_value": str(
                cubic_invariant(generic_components, (2, 3, 4))
            ),
            "conductor": "zero",
            "at_F": "tau(F) is the diag(1,-1) PGL2 translate of F",
            "local_quotient_tangent_eigenspace_dimensions_at_F": {
                "plus": quotient_tangent_eigenspaces[0],
                "minus": quotient_tangent_eigenspaces[1],
            },
            "all_moment_fiber_quotient_tangent_eigenspaces_at_F": {
                "plus": quotient_fiber_eigenspaces[0],
                "minus": quotient_fiber_eigenspaces[1],
            },
            "first_odd_invariant_dimensions": {
                "End(Sym^2)_degree_6": odd_d2[6],
                "End(Sym^3)_degree_4": odd_d3[4],
                "End(Sym^4)_degree_3": odd_d4[3],
            },
            "degree_of_d2_first_six_moment_parameter_map": sum(
                d2_numerator
            ),
            "d2_mu7_not_in_first_six_parameter_ring_certificate": {
                "prime": PRIME,
                "weighted_degree_seven_monomial_count": len(
                    degree_seven_exponents
                ),
                "evaluation_rank_with_mu7": d2_degree_seven_rank,
                "determinant_mod_prime": d2_degree_seven_determinant,
            },
            "degree_of_d2_full_moment_field": 2,
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

    print("PASS missing invariant: degree-one invariant space is moment-detected")
    print(
        "PASS missing invariant: quadratic values are",
        [q_values[component] for component in COMPONENTS],
    )
    print(
        "PASS missing invariant: mu_2 weights are",
        list(moment_two_weights),
    )
    print("PASS missing invariant: tr(A^2)=1152, before det(C)=48 in degree 5")
    print("PASS missing invariant: low-degree generator counts are 1, 4, 10")
    print(
        "PASS moment algebra: mu_1,...,mu_22 have Jacobian rank 22 modulo",
        PRIME,
    )
    print(
        "PASS moment algebra: apolar adjoint fixes every moment and "
        "has nonzero odd cubic c_234"
    )
    print("PASS moment algebra: invariant fields differ and the conductor is zero")
    print("PASS d=2 moment field: mu_7 generates the degree-five fixed field")


if __name__ == "__main__":
    main()
