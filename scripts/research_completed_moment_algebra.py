#!/usr/bin/env python3
"""Bounded exact tests for completed two-pair moment coordinates.

For V_d = End(Sym^d), d=3,4,5, this exploratory script:

* certifies full modular Jacobian rank for the initial moment prefix;
* constructs every quadratic Casimir q_(2r) and the first convenient
  apolar-odd trace invariant;
* checks the propagated d=4 moment-zero witness in degrees four and five;
* searches over a good finite field for low-weight rational expressions
  of missing quadratic invariants and of the square of the odd invariant.

The relation searches are deliberately one-sided.  A zero intersection
proves that no characteristic-zero relation with the tested support exists.
A positive intersection is only a modular candidate until reconstructed
over QQ and verified symbolically.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from math import comb, factorial
from pathlib import Path

import sympy as sp


DEFAULT_PRIME = 1_000_003
ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "completed_moment_algebra_bounded_tests.json"
)


def matrix_product_mod(
    left: list[list[int]],
    right: list[list[int]],
    prime: int,
) -> list[list[int]]:
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    result = [[0] * columns for _ in range(rows)]
    for row in range(rows):
        for pivot in range(middle):
            coefficient = left[row][pivot] % prime
            if not coefficient:
                continue
            for column in range(columns):
                result[row][column] = (
                    result[row][column]
                    + coefficient * right[pivot][column]
                ) % prime
    return result


def trace_word_mod(matrices: list[list[list[int]]], prime: int) -> int:
    product = matrices[0]
    for matrix in matrices[1:]:
        product = matrix_product_mod(product, matrix, prime)
    return sum(product[index][index] for index in range(len(product))) % prime


def sl2_matrices(d: int) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    size = d + 1
    raising = sp.zeros(size)
    lowering = sp.zeros(size)
    cartan = sp.zeros(size)
    for index in range(size):
        cartan[index, index] = 2 * index - d
        if index < d:
            raising[index + 1, index] = d - index
        if index > 0:
            lowering[index - 1, index] = index
    return raising, lowering, cartan


def casimir_projectors_exact(d: int) -> dict[int, sp.Matrix]:
    size = d + 1
    identity = sp.eye(size)

    def adjoint(matrix: sp.Matrix) -> sp.Matrix:
        return sp.kronecker_product(identity, matrix) - sp.kronecker_product(
            matrix.T, identity
        )

    raising, lowering, cartan = sl2_matrices(d)
    casimir = (
        adjoint(raising) * adjoint(lowering)
        + adjoint(lowering) * adjoint(raising)
        + adjoint(cartan) * adjoint(cartan) / 2
    )
    ambient = (d + 1) ** 2
    projectors: dict[int, sp.Matrix] = {}
    for component in range(d + 1):
        eigenvalue = 2 * component * (component + 1)
        projector = sp.eye(ambient)
        for other in range(d + 1):
            if other == component:
                continue
            other_eigenvalue = 2 * other * (other + 1)
            projector *= (
                casimir - other_eigenvalue * sp.eye(ambient)
            ) / (eigenvalue - other_eigenvalue)
        projectors[component] = projector
    return projectors


def casimir_projectors_mod(
    d: int,
    prime: int,
) -> dict[int, list[list[int]]]:
    return {
        component: [
            [
                int(sp.numer(entry)) % prime
                * pow(int(sp.denom(entry)) % prime, -1, prime)
                % prime
                for entry in row
            ]
            for row in projector.tolist()
        ]
        for component, projector in casimir_projectors_exact(d).items()
    }


def deterministic_point(
    d: int,
    sample: int,
    prime: int,
) -> list[list[int]]:
    # SplitMix64 gives independently varying, reproducible coordinates.
    size = d + 1
    state = (sample + (d << 32)) & ((1 << 64) - 1)
    result = []
    for _row in range(size):
        row_values = []
        for _column in range(size):
            state = (state + 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
            value = state
            value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9
            value &= (1 << 64) - 1
            value = (value ^ (value >> 27)) * 0x94D049BB133111EB
            value &= (1 << 64) - 1
            value ^= value >> 31
            row_values.append(value % prime)
        result.append(row_values)
    return result


def multiply_polynomials_mod(
    previous: list[list[int]],
    point: list[list[int]],
    d: int,
    prime: int,
) -> list[list[int]]:
    size = len(previous) + d
    result = [[0] * size for _ in range(size)]
    for left, previous_row in enumerate(previous):
        for right, value in enumerate(previous_row):
            if not value:
                continue
            for row in range(d + 1):
                for column in range(d + 1):
                    result[left + row][right + column] = (
                        result[left + row][right + column]
                        + value * point[row][column]
                    ) % prime
    return result


def moments_mod(
    point: list[list[int]],
    d: int,
    cutoff: int,
    prime: int,
) -> list[int]:
    factorials = [1]
    for value in range(1, d * cutoff + 1):
        factorials.append(factorials[-1] * value % prime)
    power = [[1]]
    values = []
    for order in range(1, cutoff + 1):
        power = multiply_polynomials_mod(power, point, d, prime)
        values.append(
            sum(
                factorials[index]
                * factorials[d * order - index]
                * power[index][index]
                for index in range(d * order + 1)
            )
            % prime
        )
    return values


def moment_jacobian_mod(
    point: list[list[int]],
    d: int,
    cutoff: int,
    prime: int,
) -> list[list[int]]:
    factorials = [1]
    for value in range(1, d * cutoff + 1):
        factorials.append(factorials[-1] * value % prime)
    previous = [[1]]
    rows = []
    for order in range(1, cutoff + 1):
        derivative = []
        for dual_index in range(d + 1):
            for coordinate_index in range(d + 1):
                value = 0
                for total in range(d * order + 1):
                    left = total - dual_index
                    right = total - coordinate_index
                    if (
                        0 <= left < len(previous)
                        and 0 <= right < len(previous)
                    ):
                        value += (
                            factorials[total]
                            * factorials[d * order - total]
                            * previous[left][right]
                        )
                derivative.append(order * value % prime)
        rows.append(derivative)
        previous = multiply_polynomials_mod(previous, point, d, prime)
    return rows


def component_matrices_mod(
    point: list[list[int]],
    d: int,
    projectors: dict[int, list[list[int]]],
    prime: int,
) -> dict[int, list[list[int]]]:
    size = d + 1
    factorial_diagonal = [
        factorial(index) * factorial(d - index) % prime
        for index in range(size)
    ]
    # Column-major vectorization of A=C^T D.
    vector = [
        point[column][row] * factorial_diagonal[column] % prime
        for column in range(size)
        for row in range(size)
    ]
    components = {}
    for component, projector in projectors.items():
        projected = [
            sum(left * right for left, right in zip(row, vector, strict=True))
            % prime
            for row in projector
        ]
        components[component] = [
            [
                projected[row + size * column]
                for column in range(size)
            ]
            for row in range(size)
        ]
    return components


def invariant_values_mod(
    point: list[list[int]],
    d: int,
    projectors: dict[int, list[list[int]]],
    prime: int,
) -> tuple[list[int], int]:
    components = component_matrices_mod(point, d, projectors, prime)
    quadratics = [
        trace_word_mod([components[index], components[index]], prime)
        for index in range(d + 1)
    ]
    if d == 3:
        odd = trace_word_mod(
            [components[1], components[2], components[3], components[3]],
            prime,
        )
    else:
        odd = trace_word_mod(
            [components[2], components[3], components[4]],
            prime,
        )
    return quadratics, odd


def invariant_values_exact(
    point: list[list[sp.Rational | int]],
    d: int,
) -> tuple[list[sp.Expr], sp.Expr]:
    size = d + 1
    coefficients = sp.Matrix(point)
    factorial_diagonal = sp.diag(
        *[
            factorial(index) * factorial(d - index)
            for index in range(size)
        ]
    )
    operator = coefficients.T * factorial_diagonal
    vector = sp.Matrix(
        [
            operator[row, column]
            for column in range(size)
            for row in range(size)
        ]
    )
    components = {}
    for component, projector in casimir_projectors_exact(d).items():
        projected = projector * vector
        components[component] = sp.Matrix(
            size,
            size,
            lambda row, column: projected[row + size * column],
        )
    quadratics = [
        sp.factor(sp.trace(components[index] ** 2))
        for index in range(d + 1)
    ]
    if d == 3:
        odd = sp.factor(
            sp.trace(
                components[1]
                * components[2]
                * components[3]
                * components[3]
            )
        )
    else:
        odd = sp.factor(
            sp.trace(components[2] * components[3] * components[4])
        )
    return quadratics, odd


def quadratic_jacobian_row_mod(
    point: list[list[int]],
    d: int,
    component: int,
    projectors: dict[int, list[list[int]]],
    prime: int,
) -> list[int]:
    size = d + 1
    base_value = invariant_values_mod(
        point, d, projectors, prime
    )[0][component]
    result = []
    for row in range(size):
        for column in range(size):
            basis = [[0] * size for _ in range(size)]
            basis[row][column] = 1
            basis_value = invariant_values_mod(
                basis, d, projectors, prime
            )[0][component]
            shifted = [current[:] for current in point]
            shifted[row][column] = (
                shifted[row][column] + 1
            ) % prime
            shifted_value = invariant_values_mod(
                shifted, d, projectors, prime
            )[0][component]
            result.append(
                (shifted_value - base_value - basis_value) % prime
            )
    return result


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    if not matrix or not matrix[0]:
        return 0
    reduced = [row[:] for row in matrix]
    row_count = len(reduced)
    column_count = len(reduced[0])
    rank = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(rank, row_count)
                if reduced[row][column] % prime
            ),
            None,
        )
        if pivot is None:
            continue
        reduced[rank], reduced[pivot] = reduced[pivot], reduced[rank]
        inverse = pow(reduced[rank][column] % prime, -1, prime)
        for row in range(rank + 1, row_count):
            if not reduced[row][column] % prime:
                continue
            scale = reduced[row][column] * inverse % prime
            reduced[row] = [
                (left - scale * right) % prime
                for left, right in zip(
                    reduced[row], reduced[rank], strict=True
                )
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


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


def monomial_value(
    values: list[int],
    exponents: tuple[int, ...],
    prime: int,
) -> int:
    result = 1
    for value, exponent in zip(values, exponents, strict=True):
        if exponent:
            result = result * pow(value, exponent, prime) % prime
    return result


def invariant_hilbert_coefficients(d: int, cutoff: int) -> list[int]:
    """Weight-zero-minus-weight-two coefficients for Q[V_d]^SL2."""
    coefficients = [defaultdict(int) for _ in range(cutoff + 1)]
    coefficients[0][0] = 1
    for component in range(d + 1):
        for weight in range(-2 * component, 2 * component + 1, 2):
            # Ascending degree is the unbounded-knapsack update for the
            # symmetric-algebra factor (1-t*z^weight)^(-1).
            for degree in range(1, cutoff + 1):
                for total_weight, multiplicity in list(
                    coefficients[degree - 1].items()
                ):
                    coefficients[degree][total_weight + weight] += (
                        multiplicity
                    )
    return [
        coefficients[degree][0] - coefficients[degree][2]
        for degree in range(cutoff + 1)
    ]


def hilbert_numerator(
    hilbert: list[int],
    degrees: tuple[int, ...],
) -> list[int]:
    result = hilbert[:]
    for degree in degrees:
        for index in range(len(result) - 1, degree - 1, -1):
            result[index] -= result[index - degree]
    return result


def integer_partitions(
    total: int,
    minimum: int = 1,
) -> list[tuple[int, ...]]:
    if total == 0:
        return [()]
    result = []
    for first in range(minimum, total + 1):
        for remainder in integer_partitions(total - first, first):
            result.append((first,) + remainder)
    return result


def hilbert_parameter_tests(
    d: int,
    correction_bound: int = 15,
) -> dict[str, object]:
    quotient_dimension = (d + 1) ** 2 - 3
    natural_top = (
        sum(range(1, quotient_dimension + 1)) - (d + 1) ** 2
    )
    cutoff = natural_top + correction_bound + 5
    hilbert = invariant_hilbert_coefficients(d, cutoff)
    tests = []
    for added_quadratics in range(d):
        last_natural_moment = quotient_dimension - added_quadratics
        fixed_degrees = [1] + [2] * (added_quadratics + 1)
        natural_moments = list(range(3, last_natural_moment + 1))
        first_passing_correction = None
        natural_test = None
        for total_correction in range(correction_bound + 1):
            for partition in integer_partitions(total_correction):
                if len(partition) > len(natural_moments):
                    continue
                if partition:
                    corrected_moments = (
                        natural_moments[: -len(partition)]
                        + [
                            degree + shift
                            for degree, shift in zip(
                                natural_moments[-len(partition) :],
                                partition,
                                strict=True,
                            )
                        ]
                    )
                else:
                    corrected_moments = natural_moments[:]
                degrees = tuple(fixed_degrees + corrected_moments)
                predicted_top = sum(degrees) - (d + 1) ** 2
                numerator = hilbert_numerator(hilbert, degrees)
                first_negative = next(
                    (
                        [degree, coefficient]
                        for degree, coefficient in enumerate(numerator)
                        if coefficient < 0
                    ),
                    None,
                )
                first_tail = next(
                    (
                        [degree, coefficient]
                        for degree, coefficient in enumerate(
                            numerator[predicted_top + 1 :],
                            start=predicted_top + 1,
                        )
                        if coefficient
                    ),
                    None,
                )
                record = {
                    "added_quadratics": added_quadratics,
                    "degrees": list(degrees),
                    "predicted_top_degree": predicted_top,
                    "first_negative_through_cutoff": first_negative,
                    "first_nonzero_after_predicted_top_through_cutoff": (
                        first_tail
                    ),
                    "checked_through_degree": cutoff,
                }
                if total_correction == 0:
                    natural_test = record
                if first_negative is None and first_tail is None:
                    record["numerator_coefficient_sum"] = sum(
                        numerator[: predicted_top + 1]
                    )
                    record["correction_partition"] = list(partition)
                    record["total_degree_correction"] = total_correction
                    first_passing_correction = record
                    break
            if first_passing_correction is not None:
                break
        assert natural_test is not None
        tests.append(
            {
                "natural_sequence": natural_test,
                "first_hilbert_compatible_correction_within_bound": (
                    first_passing_correction
                ),
            }
        )
    return {
        "interpretation": (
            "necessary Hilbert-series tests only; compatibility does not "
            "prove algebraic independence or a nullcone zero fiber"
        ),
        "hilbert_coefficients_degrees_0_through_5": hilbert[:6],
        "tests_by_number_of_added_quadratics": tests,
    }


@dataclass(frozen=True)
class Sample:
    moments: list[int]
    quadratics: list[int]
    odd: int


def relation_intersection(
    samples: list[Sample],
    base_quadratics: tuple[int, ...],
    target: str,
    target_degree: int,
    weight: int,
    prime: int,
) -> tuple[int, int, int, int]:
    base_weights = tuple(range(1, weight + 1)) + (2,) * len(
        base_quadratics
    )
    q_exponents = weighted_exponents(base_weights, weight)
    p_exponents = weighted_exponents(
        base_weights,
        weight - target_degree,
    )
    q_matrix = []
    p_matrix = []
    for sample in samples:
        base_values = sample.moments[:weight] + [
            sample.quadratics[index] for index in base_quadratics
        ]
        if target == "odd_square":
            target_value = sample.odd * sample.odd % prime
        else:
            target_value = sample.quadratics[int(target)]
        q_matrix.append(
            [
                monomial_value(base_values, exponents, prime)
                for exponents in q_exponents
            ]
        )
        p_matrix.append(
            [
                target_value
                * monomial_value(base_values, exponents, prime)
                % prime
                for exponents in p_exponents
            ]
        )
    rank_q = rank_mod(q_matrix, prime)
    rank_p = rank_mod(p_matrix, prime)
    rank_full = rank_mod(
        [
            left + right
            for left, right in zip(q_matrix, p_matrix, strict=True)
        ],
        prime,
    )
    return len(q_exponents), len(p_exponents), rank_full, (
        rank_q + rank_p - rank_full
    )


def propagated_witness(d: int) -> list[list[int]]:
    assert d in (4, 5)
    witness = [
        [-1, 2, 0, 0, 0],
        [-sp.Rational(3, 2), 2, 6, 0, 0],
        [-sp.Rational(1, 2), sp.Rational(3, 2), 6, 6, 0],
        [0, 1, sp.Rational(3, 2), 2, 2],
        [0, 0, -sp.Rational(1, 2), -sp.Rational(3, 2), -1],
    ]
    if d == 4:
        return witness
    result = [[sp.Integer(0)] * 6 for _ in range(6)]
    for row in range(5):
        for column in range(5):
            result[row][column] += witness[row][column]
            result[row + 1][column + 1] += witness[row][column]
    return result


def reduce_point_mod(
    point: list[list[sp.Rational | int]],
    prime: int,
) -> list[list[int]]:
    return [
        [
            int(sp.numer(entry)) % prime
            * pow(int(sp.denom(entry)) % prime, -1, prime)
            % prime
            for entry in row
        ]
        for row in point
    ]


def run_degree(
    d: int,
    max_weight: int,
    extra_samples: int,
    prime: int,
) -> dict[str, object]:
    quotient_dimension = (d + 1) ** 2 - 3
    projectors = casimir_projectors_mod(d, prime)
    jacobian_point = deterministic_point(d, 97, prime)
    jacobian = moment_jacobian_mod(
        jacobian_point, d, quotient_dimension + 15, prime
    )
    jacobian_rank = rank_mod(jacobian[:quotient_dimension], prime)
    print(
        f"d={d} moment_jacobian orders=1..{quotient_dimension} "
        f"rank={jacobian_rank}"
    )
    result: dict[str, object] = {
        "quotient_dimension": quotient_dimension,
        "moment_jacobian_rank": jacobian_rank,
        "quadratic_completion_component_indices": list(range(1, d)),
        "odd_invariant": (
            "tr(A_2 A_4 A_6^2)"
            if d == 3
            else "tr(A_4 A_6 A_8)"
        ),
        "odd_invariant_degree": 4 if d == 3 else 3,
        "mu2_weights_on_q_components": [
            comb(2 * d + 1, d - component)
            for component in range(d + 1)
        ],
        "relation_tests": [],
        "hilbert_parameter_tests": hilbert_parameter_tests(d),
    }
    hilbert_tests = result["hilbert_parameter_tests"]
    assert isinstance(hilbert_tests, dict)
    tests_by_quadratics = hilbert_tests[
        "tests_by_number_of_added_quadratics"
    ]
    assert isinstance(tests_by_quadratics, list)
    parameter_jacobians = []
    quadratic_jacobians = {
        component: quadratic_jacobian_row_mod(
            jacobian_point, d, component, projectors, prime
        )
        for component in range(1, d)
    }
    for added_quadratics, test in enumerate(tests_by_quadratics):
        assert isinstance(test, dict)
        candidate = test[
            "first_hilbert_compatible_correction_within_bound"
        ]
        if not isinstance(candidate, dict):
            continue
        degrees = candidate["degrees"]
        assert isinstance(degrees, list)
        # The fixed prefix is mu_1, then mu_2 and the added quadratic
        # rows.  Every remaining degree names its unique moment.
        moment_orders = [1, 2] + degrees[added_quadratics + 2 :]
        candidate_rows = [
            jacobian[order - 1] for order in moment_orders
        ] + [
            quadratic_jacobians[component]
            for component in range(1, added_quadratics + 1)
        ]
        candidate_rank = rank_mod(candidate_rows, prime)
        parameter_jacobians.append(
            {
                "added_quadratic_components": list(
                    range(1, added_quadratics + 1)
                ),
                "moment_orders": moment_orders,
                "rank": candidate_rank,
            }
        )
    result["hilbert_compatible_parameter_jacobians"] = (
        parameter_jacobians
    )

    odd_degree = 4 if d == 3 else 3
    # q_0 is generated by mu_1^2 and mu_2 supplies one further linear
    # combination.  These d-1 quadratics complete the degree-two space.
    completing_quadratics = tuple(range(1, d))
    targets = [
        ((), str(component), 2, f"q_{2 * component}_over_moments")
        for component in range(1, d + 1)
    ]
    sparse_quadratics = (1,)
    proposed_quadratics = (1, 3)
    for component in range(2, d + 1):
        targets.append(
            (
                sparse_quadratics,
                str(component),
                2,
                f"q_{2 * component}_over_moments_q2",
            )
        )
    targets.extend(
        [
            (
                (),
                "odd_square",
                2 * odd_degree,
                "odd_square_over_moments",
            ),
            (
                sparse_quadratics,
                "odd_square",
                2 * odd_degree,
                "odd_square_over_moments_q2",
            ),
            (
                proposed_quadratics,
                "odd_square",
                2 * odd_degree,
                "odd_square_over_moments_q2_q6",
            ),
            (
                completing_quadratics,
                "odd_square",
                2 * odd_degree,
                "odd_square_over_quadratically_completed_moments",
            ),
        ]
    )

    largest_columns = 0
    for base_quadratics, _target, target_degree, _label in targets:
        if max_weight < target_degree:
            continue
        weights = tuple(range(1, max_weight + 1)) + (2,) * len(
            base_quadratics
        )
        largest_columns = max(
            largest_columns,
            len(weighted_exponents(weights, max_weight))
            + len(
                weighted_exponents(
                    weights, max_weight - target_degree
                )
            ),
        )
    sample_count = largest_columns + extra_samples
    samples = []
    for sample_index in range(1, sample_count + 1):
        point = deterministic_point(d, sample_index, prime)
        quadratics, odd = invariant_values_mod(
            point, d, projectors, prime
        )
        samples.append(
            Sample(
                moments=moments_mod(point, d, max_weight, prime),
                quadratics=quadratics,
                odd=odd,
            )
        )
    assert any(sample.odd for sample in samples)
    mu2_weights = result["mu2_weights_on_q_components"]
    assert isinstance(mu2_weights, list)
    assert all(
        sample.moments[1]
        == sum(
            weight * value
            for weight, value in zip(
                mu2_weights, sample.quadratics, strict=True
            )
        )
        % prime
        for sample in samples
    )

    for base_quadratics, target, target_degree, label in targets:
        if max_weight < target_degree:
            continue
        for weight in range(target_degree, max_weight + 1):
            q_columns, p_columns, rank, intersection = (
                relation_intersection(
                    samples,
                    base_quadratics,
                    target,
                    target_degree,
                    weight,
                    prime,
                )
            )
            print(
                f"d={d} test={label} weight={weight} "
                f"columns={q_columns}+{p_columns} rank={rank} "
                f"relation_intersection={intersection}"
            )
            relation_tests = result["relation_tests"]
            assert isinstance(relation_tests, list)
            relation_tests.append(
                {
                    "label": label,
                    "weight": weight,
                    "base_columns": q_columns,
                    "target_multiple_columns": p_columns,
                    "combined_rank": rank,
                    "relation_intersection": intersection,
                }
            )

    if d in (4, 5):
        rational_witness = propagated_witness(d)
        exact_witness_quadratics, exact_witness_odd = (
            invariant_values_exact(rational_witness, d)
        )
        witness = reduce_point_mod(rational_witness, prime)
        witness_moments = moments_mod(
            witness, d, max(12, max_weight), prime
        )
        witness_quadratics, witness_odd = invariant_values_mod(
            witness, d, projectors, prime
        )
        assert not any(witness_moments)
        assert [
            int(value) % prime for value in exact_witness_quadratics
        ] == witness_quadratics
        assert int(exact_witness_odd) % prime == witness_odd
        print(
            f"d={d} propagated_witness "
            f"q={exact_witness_quadratics} "
            f"odd={exact_witness_odd}"
        )
        result["propagated_moment_zero_witness"] = {
            "moments_checked_through": max(12, max_weight),
            "quadratic_values": [
                str(value) for value in exact_witness_quadratics
            ],
            "odd_invariant_value": str(exact_witness_odd),
        }
    return result


def signed(value: int, prime: int) -> int:
    return value - prime if value > prime // 2 else value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--degrees",
        type=int,
        nargs="+",
        choices=(3, 4, 5),
        default=(3, 4, 5),
    )
    parser.add_argument("--max-weight", type=int, default=10)
    parser.add_argument("--extra-samples", type=int, default=3)
    parser.add_argument("--prime", type=int, default=DEFAULT_PRIME)
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    payload = {
        "status": (
            "bounded modular nonrelation tests; positive intersections "
            "would require characteristic-zero reconstruction"
        ),
        "prime": arguments.prime,
        "max_weight": arguments.max_weight,
        "extra_samples": arguments.extra_samples,
        "degrees": {},
    }
    for d in arguments.degrees:
        result = run_degree(
            d,
            arguments.max_weight,
            arguments.extra_samples,
            arguments.prime,
        )
        payload["degrees"][str(d)] = result
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
