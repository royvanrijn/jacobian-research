#!/usr/bin/env python3
"""Sparse finite-field evaluation and tangent solving for arithmetic DAGs.

The module is intentionally independent of the F2 band geometry.  It accepts
the node protocol used by ``CircuitDAG`` (constant, variable, scale, add, and
multiply), an embedding of the constant keys into a prime field, and a list
of output roots.  At a chosen point it evaluates both the circuit and its
sparse Jacobian, then solves the full linearized system by deterministic
sparse Gaussian elimination.

This is an exploratory good-reduction tool.  A finite-field point is evidence
about one reduction, not a characteristic-zero solution; a failed Newton
step is not an obstruction certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


FieldKey = tuple[tuple[int, int], tuple[int, int]]
ConstantEmbedding = Callable[[FieldKey], int]


@dataclass(frozen=True)
class SparseEvaluation:
    values: list[int]
    gradients: list[dict[int, int]] | None
    variable_nodes: list[int]
    variable_names: list[str]
    gradient_entry_count: int


@dataclass(frozen=True)
class LinearizedSolve:
    rank: int
    inconsistent_rows: list[int]
    correction: dict[str, int]
    pivot_variables: tuple[str, ...]


@dataclass(frozen=True)
class LeftCokernelBasis:
    rank: int
    functionals: tuple[dict[int, int], ...]


def quadratic_key_embedding(prime: int, rho: int) -> ConstantEmbedding:
    """Embed serialized ``b+a*rho`` constants in ``GF(prime)``."""

    def embed(key: FieldKey) -> int:
        (constant_p, constant_q), (rho_p, rho_q) = key
        return (
            constant_p * pow(constant_q, -1, prime)
            + rho_p * pow(rho_q, -1, prime) * rho
        ) % prime

    return embed


def _add_gradients(
    left: dict[int, int],
    right: dict[int, int],
    prime: int,
    left_scale: int = 1,
    right_scale: int = 1,
) -> dict[int, int]:
    result = {
        key: left_scale * value % prime
        for key, value in left.items()
        if left_scale * value % prime
    }
    for key, value in right.items():
        combined = (result.get(key, 0) + right_scale * value) % prime
        if combined:
            result[key] = combined
        else:
            result.pop(key, None)
    return result


def evaluate(
    dag: object,
    variable_values: dict[str, int],
    prime: int,
    constant_embedding: ConstantEmbedding,
    *,
    with_jacobian: bool,
) -> SparseEvaluation:
    """Evaluate a topologically ordered arithmetic DAG over ``GF(prime)``."""

    variable_nodes = [
        index for index, node in enumerate(dag.nodes) if node[0] == "var"
    ]
    variable_names = [str(dag.nodes[node][1]) for node in variable_nodes]
    variable_index = {node: index for index, node in enumerate(variable_nodes)}
    values = [0] * len(dag.nodes)
    gradients: list[dict[int, int]] | None = (
        [dict() for _ in dag.nodes] if with_jacobian else None
    )
    entry_count = 0

    for index, node in enumerate(dag.nodes):
        kind = node[0]
        if kind == "const":
            values[index] = constant_embedding(node[1])
        elif kind == "var":
            values[index] = variable_values.get(str(node[1]), 0) % prime
            if gradients is not None:
                gradients[index] = {variable_index[index]: 1}
        elif kind == "scale":
            scalar = constant_embedding(node[1])
            child = int(node[2])
            values[index] = scalar * values[child] % prime
            if gradients is not None:
                gradients[index] = {
                    key: scalar * value % prime
                    for key, value in gradients[child].items()
                    if scalar * value % prime
                }
        elif kind == "add":
            left, right = int(node[1]), int(node[2])
            values[index] = (values[left] + values[right]) % prime
            if gradients is not None:
                gradients[index] = _add_gradients(
                    gradients[left], gradients[right], prime
                )
        elif kind == "mul":
            left, right = int(node[1]), int(node[2])
            values[index] = values[left] * values[right] % prime
            if gradients is not None and (values[left] or values[right]):
                gradients[index] = _add_gradients(
                    gradients[left],
                    gradients[right],
                    prime,
                    values[right],
                    values[left],
                )
        else:
            raise AssertionError(f"unknown circuit operation {kind!r}")
        if gradients is not None:
            entry_count += len(gradients[index])

    return SparseEvaluation(
        values=values,
        gradients=gradients,
        variable_nodes=variable_nodes,
        variable_names=variable_names,
        gradient_entry_count=entry_count,
    )


def solve_linearization(
    evaluation: SparseEvaluation,
    roots: list[int],
    prime: int,
    *,
    free_values: dict[str, int] | None = None,
    allowed_variables: set[str] | None = None,
    right_hand_side: list[int] | None = None,
    prescribed_values: dict[str, int] | None = None,
) -> LinearizedSolve:
    """Solve ``Jacobian*correction=-value`` over ``GF(prime)``.

    ``allowed_variables`` restricts the Jacobian to a coordinate subspace.
    ``free_values`` may prescribe nonpivot coordinates in the resulting
    affine tangent fiber.  Unspecified free coordinates are zero.  Together
    these options expose sparse tangent charts without constructing a dense
    kernel basis.  By default the system is ``J*correction=-F`` at the
    evaluated point; ``right_hand_side`` replaces ``-F`` for fixed-Jacobian
    deformation lifting.
    ``prescribed_values`` fixes arbitrary coordinates, including would-be
    pivots, before elimination and is useful for formal gauge conditions.
    """

    if evaluation.gradients is None:
        raise ValueError("the evaluation does not contain a Jacobian")
    variable_count = len(evaluation.variable_nodes)
    augmented_column = variable_count
    name_to_index = {
        name: index for index, name in enumerate(evaluation.variable_names)
    }
    allowed_indices = (
        None
        if allowed_variables is None
        else {
            index
            for index, name in enumerate(evaluation.variable_names)
            if name in allowed_variables
        }
    )
    prescribed_by_index: dict[int, int] = {}
    for name, value in (prescribed_values or {}).items():
        index = name_to_index.get(name)
        if index is None:
            raise KeyError(f"unknown circuit variable {name!r}")
        if allowed_indices is not None and index not in allowed_indices:
            raise ValueError(f"variable {name!r} is outside the allowed chart")
        prescribed_by_index[index] = value % prime
    rows: list[dict[int, int]] = []
    if right_hand_side is not None and len(right_hand_side) != len(roots):
        raise ValueError("right-hand-side length does not match the roots")
    for row_index, root in enumerate(roots):
        row = {
            index: value
            for index, value in evaluation.gradients[root].items()
            if (allowed_indices is None or index in allowed_indices)
            and index not in prescribed_by_index
        }
        row_right_hand_side = (
            -evaluation.values[root] % prime
            if right_hand_side is None
            else right_hand_side[row_index] % prime
        )
        row_right_hand_side = (
            row_right_hand_side
            - sum(
                evaluation.gradients[root].get(index, 0) * value
                for index, value in prescribed_by_index.items()
            )
        ) % prime
        if row_right_hand_side:
            row[augmented_column] = row_right_hand_side
        rows.append(row)

    pivots: dict[int, dict[int, int]] = {}
    inconsistent: list[int] = []
    for row_index, row in enumerate(rows):
        while True:
            columns = [
                column
                for column, value in row.items()
                if column < variable_count and value
            ]
            if not columns:
                if row.get(augmented_column, 0):
                    inconsistent.append(row_index)
                break
            column = min(columns)
            pivot = pivots.get(column)
            if pivot is None:
                inverse = pow(row[column], -1, prime)
                pivots[column] = {
                    key: value * inverse % prime
                    for key, value in row.items()
                    if value * inverse % prime
                }
                break
            scalar = row[column]
            for key, value in pivot.items():
                combined = (row.get(key, 0) - scalar * value) % prime
                if combined:
                    row[key] = combined
                else:
                    row.pop(key, None)

    solution = [0] * variable_count
    for index, value in prescribed_by_index.items():
        solution[index] = value
    prescribed = free_values or {}
    pivot_columns = set(pivots)
    for name, value in prescribed.items():
        index = name_to_index.get(name)
        if index is None:
            raise KeyError(f"unknown circuit variable {name!r}")
        if allowed_indices is not None and index not in allowed_indices:
            raise ValueError(f"variable {name!r} is outside the allowed chart")
        if index in prescribed_by_index:
            raise ValueError(f"variable {name!r} is already prescribed")
        if index in pivot_columns:
            raise ValueError(f"cannot prescribe pivot variable {name!r}")
        solution[index] = value % prime
    if not inconsistent:
        for column in sorted(pivots, reverse=True):
            row = pivots[column]
            solution[column] = (
                row.get(augmented_column, 0)
                - sum(
                    value * solution[key]
                    for key, value in row.items()
                    if key < variable_count and key != column
                )
            ) % prime
    correction = {
        evaluation.variable_names[index]: value
        for index, value in enumerate(solution)
        if value
    }
    return LinearizedSolve(
        rank=len(pivots),
        inconsistent_rows=inconsistent,
        correction=correction,
        pivot_variables=tuple(
            evaluation.variable_names[index] for index in pivots
        ),
    )


def left_cokernel_basis(
    evaluation: SparseEvaluation,
    roots: list[int],
    prime: int,
    *,
    allowed_variables: set[str] | None = None,
) -> LeftCokernelBasis:
    """Compute sparse row functionals annihilating the selected Jacobian."""

    if evaluation.gradients is None:
        raise ValueError("the evaluation does not contain a Jacobian")
    allowed_indices = (
        None
        if allowed_variables is None
        else {
            index
            for index, name in enumerate(evaluation.variable_names)
            if name in allowed_variables
        }
    )
    pivots: dict[int, tuple[dict[int, int], dict[int, int]]] = {}
    functionals: list[dict[int, int]] = []
    for row_index, root in enumerate(roots):
        row = {
            index: value
            for index, value in evaluation.gradients[root].items()
            if allowed_indices is None or index in allowed_indices
        }
        combination = {row_index: 1}
        while row:
            column = min(row)
            pivot = pivots.get(column)
            if pivot is None:
                inverse = pow(row[column], -1, prime)
                row = {
                    key: value * inverse % prime
                    for key, value in row.items()
                    if value * inverse % prime
                }
                combination = {
                    key: value * inverse % prime
                    for key, value in combination.items()
                    if value * inverse % prime
                }
                pivots[column] = (row, combination)
                break
            pivot_row, pivot_combination = pivot
            scalar = row[column]
            for key, value in pivot_row.items():
                combined = (row.get(key, 0) - scalar * value) % prime
                if combined:
                    row[key] = combined
                else:
                    row.pop(key, None)
            for key, value in pivot_combination.items():
                combined = (
                    combination.get(key, 0) - scalar * value
                ) % prime
                if combined:
                    combination[key] = combined
                else:
                    combination.pop(key, None)
        else:
            functionals.append(combination)

    if len(functionals) != len(roots) - len(pivots):
        raise AssertionError("left-cokernel dimension mismatch")
    return LeftCokernelBasis(
        rank=len(pivots),
        functionals=tuple(functionals),
    )


def apply_row_functionals(
    functionals: tuple[dict[int, int], ...],
    vector: list[int],
    prime: int,
) -> list[int]:
    """Apply sparse row-coordinate functionals to an equation vector."""

    return [
        sum(coefficient * vector[index] for index, coefficient in row.items())
        % prime
        for row in functionals
    ]


def evaluate_truncated_series(
    dag: object,
    variable_series: dict[str, list[int]],
    prime: int,
    constant_embedding: ConstantEmbedding,
    order: int,
) -> list[tuple[int, ...]]:
    """Evaluate an arithmetic DAG in ``GF(prime)[t]/(t^(order+1))``."""

    if order < 0:
        raise ValueError("the truncation order must be nonnegative")
    zero = (0,) * (order + 1)
    series: list[tuple[int, ...]] = [zero] * len(dag.nodes)
    for index, node in enumerate(dag.nodes):
        kind = node[0]
        if kind == "const":
            series[index] = (constant_embedding(node[1]),) + (0,) * order
        elif kind == "var":
            coefficients = variable_series.get(str(node[1]), [])
            series[index] = tuple(
                coefficients[degree] % prime
                if degree < len(coefficients)
                else 0
                for degree in range(order + 1)
            )
        elif kind == "scale":
            scalar = constant_embedding(node[1])
            child = series[int(node[2])]
            series[index] = tuple(
                scalar * coefficient % prime for coefficient in child
            )
        elif kind == "add":
            left = series[int(node[1])]
            right = series[int(node[2])]
            series[index] = tuple(
                (left[degree] + right[degree]) % prime
                for degree in range(order + 1)
            )
        elif kind == "mul":
            left = series[int(node[1])]
            right = series[int(node[2])]
            series[index] = tuple(
                sum(
                    left[left_degree] * right[degree - left_degree]
                    for left_degree in range(degree + 1)
                )
                % prime
                for degree in range(order + 1)
            )
        else:
            raise AssertionError(f"unknown circuit operation {kind!r}")
    return series


def apply_correction(
    point: dict[str, int], correction: dict[str, int], prime: int
) -> dict[str, int]:
    return apply_scaled_correction(point, correction, 1, prime)


def apply_scaled_correction(
    point: dict[str, int],
    correction: dict[str, int],
    scale: int,
    prime: int,
) -> dict[str, int]:
    """Return ``point + scale*correction`` over ``GF(prime)``."""

    result = dict(point)
    for name, value in correction.items():
        combined = (result.get(name, 0) + scale * value) % prime
        if combined:
            result[name] = combined
        else:
            result.pop(name, None)
    return result


def interpolate_corrections(
    base: dict[str, int],
    unit_point: dict[str, int],
    parameter: int,
    prime: int,
) -> dict[str, int]:
    """Return ``base + parameter*(unit_point-base)`` over ``GF(prime)``."""

    result: dict[str, int] = {}
    for name in base.keys() | unit_point.keys():
        value = (
            base.get(name, 0)
            + parameter * (unit_point.get(name, 0) - base.get(name, 0))
        ) % prime
        if value:
            result[name] = value
    return result


def _trim_polynomial(coefficients: list[int], prime: int) -> list[int]:
    result = [coefficient % prime for coefficient in coefficients]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result or [0]


def _multiply_polynomials(
    left: list[int], right: list[int], prime: int
) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            result[left_degree + right_degree] = (
                result[left_degree + right_degree]
                + left_coefficient * right_coefficient
            ) % prime
    return _trim_polynomial(result, prime)


def evaluate_univariate(
    coefficients: list[int], value: int, prime: int
) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % prime
    return result


def interpolate_consecutive_values(
    values: list[int], prime: int, degree_bound: int
) -> list[int]:
    """Interpolate a bounded-degree polynomial from values at 0,1,... ."""

    if len(values) < degree_bound + 1:
        raise ValueError("not enough samples for the requested degree bound")
    samples = values[: degree_bound + 1]
    result = [0]
    for index, sample in enumerate(samples):
        basis = [1]
        denominator = 1
        for other in range(degree_bound + 1):
            if other == index:
                continue
            basis = _multiply_polynomials(basis, [-other, 1], prime)
            denominator = denominator * (index - other) % prime
        scale = sample * pow(denominator, -1, prime) % prime
        if len(result) < len(basis):
            result.extend([0] * (len(basis) - len(result)))
        for degree, coefficient in enumerate(basis):
            result[degree] = (result[degree] + scale * coefficient) % prime
    result = _trim_polynomial(result, prime)
    if any(
        evaluate_univariate(result, value, prime) != sample % prime
        for value, sample in enumerate(values)
    ):
        raise AssertionError("samples exceed the asserted polynomial degree")
    return result


def polynomial_gcd_modp(
    left: list[int], right: list[int], prime: int
) -> list[int]:
    """Return the monic gcd of two coefficient lists over ``GF(prime)``."""

    def remainder(dividend: list[int], divisor: list[int]) -> list[int]:
        active = _trim_polynomial(dividend, prime)
        divisor = _trim_polynomial(divisor, prime)
        if divisor == [0]:
            raise ZeroDivisionError("polynomial division by zero")
        inverse = pow(divisor[-1], -1, prime)
        while active != [0] and len(active) >= len(divisor):
            shift = len(active) - len(divisor)
            scale = active[-1] * inverse % prime
            for degree, coefficient in enumerate(divisor):
                active[degree + shift] = (
                    active[degree + shift] - scale * coefficient
                ) % prime
            active = _trim_polynomial(active, prime)
        return active

    left = _trim_polynomial(left, prime)
    right = _trim_polynomial(right, prime)
    while right != [0]:
        left, right = right, remainder(left, right)
    if left == [0]:
        return [0]
    inverse = pow(left[-1], -1, prime)
    return [(coefficient * inverse) % prime for coefficient in left]
