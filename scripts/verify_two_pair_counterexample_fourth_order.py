#!/usr/bin/env python3
"""Exact fourth-order continuation at the displayed SIC2C4 point.

This script imports the established local-moduli checker and works only on
the reduced five-plane in its quadratic obstruction cone.  Polynomial
dependence on the five plane coordinates is retained as sparse homogeneous
coefficient maps; no expanded SymPy contraction is used.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from pathlib import Path

import sympy as sp
from flint import fmpq_mat


ROOT = Path(__file__).resolve().parents[1]
LOCAL_SCRIPT = ROOT / "scripts" / "verify_two_pair_counterexample_local_moduli.py"
LOCAL_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_counterexample_local_moduli.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_counterexample_fourth_order.json"
)

spec = importlib.util.spec_from_file_location("sic2_local", LOCAL_SCRIPT)
assert spec and spec.loader
local = importlib.util.module_from_spec(spec)
spec.loader.exec_module(local)

seed = local.seed
BASIS = local.BASIS
Exponent5 = tuple[int, int, int, int, int]
ScalarMap = dict[Exponent5, sp.Rational]
VectorMap = dict[Exponent5, sp.Matrix]
MAX_ORDER = 57
TAIL_ORDER_COUNT = MAX_ORDER - 11
MASTER_REPLAY_ORDER = 71


def add_exponents(left: Exponent5, right: Exponent5) -> Exponent5:
    return tuple(left[index] + right[index] for index in range(5))


def compositions(total: int, slots: int = 5) -> list[Exponent5]:
    result: list[Exponent5] = []

    def recurse(prefix: tuple[int, ...], remaining: int) -> None:
        if len(prefix) == slots - 1:
            result.append(prefix + (remaining,))
            return
        for value in range(remaining + 1):
            recurse(prefix + (value,), remaining - value)

    recurse((), total)
    return result


def add_scalar_maps(*terms: tuple[sp.Rational, ScalarMap]) -> ScalarMap:
    result: ScalarMap = {}
    for scale, polynomial in terms:
        for exponent, coefficient in polynomial.items():
            result[exponent] = (
                result.get(exponent, sp.Rational(0))
                + scale * coefficient
            )
    return {
        exponent: sp.factor(coefficient)
        for exponent, coefficient in result.items()
        if coefficient
    }


def bilinear_map(
    left: VectorMap,
    right: VectorMap,
    matrix: sp.Matrix,
) -> ScalarMap:
    result: ScalarMap = {}
    for left_exponent, left_vector in left.items():
        for right_exponent, right_vector in right.items():
            exponent = add_exponents(left_exponent, right_exponent)
            coefficient = (left_vector.T * matrix * right_vector)[0]
            result[exponent] = (
                result.get(exponent, sp.Rational(0)) + coefficient
            )
    return {
        exponent: sp.factor(coefficient)
        for exponent, coefficient in result.items()
        if coefficient
    }


def vector_to_polynomial(vector: sp.Matrix) -> seed.Polynomial:
    result: seed.Polynomial = {}
    for index, coefficient in enumerate(vector):
        if not coefficient:
            continue
        value = sp.Rational(coefficient)
        result[BASIS[index]] = Fraction(
            int(value.p),
            int(value.q),
        )
    return result


def add_numerator(
    target: dict[Exponent5, seed.Polynomial],
    exponent: Exponent5,
    numerator: seed.Polynomial,
) -> None:
    target[exponent] = seed.add(target.get(exponent, {}), numerator)


def contraction_weights(
    f_power: seed.Polynomial,
    numerator_degree: int,
) -> dict[seed.Exponent, sp.Rational]:
    """Contraction of each balanced numerator monomial against one F power."""
    by_difference: dict[int, list[tuple[seed.Exponent, Fraction]]] = {}
    for exponent, coefficient in f_power.items():
        by_difference.setdefault(exponent[0] - exponent[2], []).append(
            (exponent, coefficient)
        )
    result: dict[seed.Exponent, sp.Rational] = {}
    for xi1 in range(numerator_degree + 1):
        for z1 in range(numerator_degree + 1):
            numerator_exponent = (
                xi1,
                numerator_degree - xi1,
                z1,
                numerator_degree - z1,
            )
            value = Fraction(0)
            for f_exponent, coefficient in by_difference.get(
                z1 - xi1,
                [],
            ):
                total_xi1 = xi1 + f_exponent[0]
                total_xi2 = numerator_degree - xi1 + f_exponent[1]
                value += (
                    coefficient
                    * local.factorial(total_xi1)
                    * local.factorial(total_xi2)
                )
            if value:
                result[numerator_exponent] = local.rational(value)
    return result


@lru_cache(maxsize=None)
def even_height_coefficients(
    numerator_degree: int,
    xi1: int,
    z1: int,
) -> tuple[tuple[int, int], ...]:
    """Even coefficients of (1+t)^z1 (1-t)^(d-xi1)."""
    coefficients = [0] * (numerator_degree - xi1 + z1 + 1)
    for left in range(z1 + 1):
        for right in range(numerator_degree - xi1 + 1):
            coefficients[left + right] += (
                comb(z1, left)
                * comb(numerator_degree - xi1, right)
                * (-1) ** right
            )
    return tuple(
        (degree, coefficient)
        for degree, coefficient in enumerate(coefficients)
        if degree % 2 == 0 and coefficient
    )


@lru_cache(maxsize=None)
def generating_weights(
    numerator_degree: int,
    f_order: int,
) -> dict[seed.Exponent, sp.Rational]:
    """All-order Hopf/beta coefficient identity for E(A F^n)."""
    result: dict[seed.Exponent, sp.Rational] = {}
    radial_factor = factorial(4 * f_order + numerator_degree + 1)
    for xi1 in range(numerator_degree + 1):
        for z1 in range(numerator_degree + 1):
            target = f_order - xi1 + z1
            value = Fraction(0)
            for index in range(f_order + 1):
                if not 0 <= target <= f_order + 2 * index:
                    continue
                height_integral = sum(
                    Fraction(coefficient, 2 * index + degree + 1)
                    for degree, coefficient in even_height_coefficients(
                        numerator_degree,
                        xi1,
                        z1,
                    )
                )
                value += (
                    (-1) ** index
                    * comb(f_order, index)
                    * comb(f_order + 2 * index, target)
                    * height_integral
                )
            value *= Fraction(
                radial_factor,
                2 ** (
                    f_order
                    + z1
                    + numerator_degree
                    - xi1
                ),
            )
            if value:
                result[
                    (
                        xi1,
                        numerator_degree - xi1,
                        z1,
                        numerator_degree - z1,
                    )
                ] = local.rational(value)
    return result


def contraction_map(
    numerators: dict[Exponent5, seed.Polynomial],
    weights: dict[seed.Exponent, sp.Rational],
) -> ScalarMap:
    result: ScalarMap = {}
    for exponent, numerator in numerators.items():
        value = sum(
            local.rational(coefficient) * weights.get(monomial, 0)
            for monomial, coefficient in numerator.items()
        )
        if value:
            result[exponent] = sp.factor(value)
    return result


def coefficient_vector(
    maps: list[ScalarMap],
    exponent: Exponent5,
) -> sp.Matrix:
    return sp.Matrix(
        [polynomial.get(exponent, sp.Rational(0)) for polynomial in maps]
    )


def main() -> None:
    stored = json.loads(LOCAL_ARTIFACT.read_text())
    witness = seed.witness()[0]

    # Check the generating formula against direct contraction before using it
    # for the higher numerator degrees.
    witness_power = seed.monomial(seed.ZERO)
    for order in range(4):
        assert generating_weights(4, order) == contraction_weights(
            witness_power,
            4,
        )
        witness_power = seed.multiply(witness_power, witness)

    # The complete fourth coefficient is a single shifted row sequence on
    # numerator degrees 4, 8, 12, and 16.  Beta expansion of
    # generating_weights shows symbolically that its normalized entries have
    # a common denominator and numerator degree at most 41.  The exact replay
    # below checks the resulting 42-row basis and eighteen subsequent rows.
    def master_row(order: int) -> list[str]:
        blocks: list[sp.Rational] = []
        for degree, power_order, scale in (
            (4, order, sp.Rational(1)),
            (8, order - 1, sp.Rational(order)),
            (12, order - 2, sp.Rational(order * (order - 1), 2)),
            (
                16,
                order - 3,
                sp.Rational(order * (order - 1) * (order - 2), 24),
            ),
        ):
            weights = generating_weights(degree, power_order)
            blocks.extend(
                scale
                * weights.get(
                    (xi1, degree - xi1, z1, degree - z1),
                    sp.Rational(0),
                )
                for xi1 in range(degree + 1)
                for z1 in range(degree + 1)
            )
        return [str(value) for value in blocks]

    master_rows = [
        master_row(order)
        for order in range(12, MASTER_REPLAY_ORDER + 1)
    ]
    master_basis_rank = fmpq_mat(master_rows[:42]).rank()
    master_replay_rank = fmpq_mat(master_rows).rank()
    assert master_basis_rank == master_replay_rank == 42
    tangent_basis = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in stored["tangent_basis_columns"]
        ]
    )
    tangent_change = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in stored["tangent_coordinate_change_orbit_radial_first"]
        ]
    )
    residual_parametrization = sp.zeros(9, 5)
    for index in range(5):
        residual_parametrization[4 + index, index] = 1
    residual_parametrization[1, :] = sp.Matrix(
        [[sp.Rational(8, 35), 0, 2, -3, 3]]
    )
    residual_parametrization[2, :] = sp.Matrix(
        [[sp.Rational(43, 28), -6, 12, -12, 12]]
    )
    residual_parametrization[3, :] = sp.Matrix(
        [[sp.Rational(251, 105), -8, 12, -12, 13]]
    )
    residual_vectors = (
        tangent_basis
        * tangent_change[:, 4:]
        * residual_parametrization
    )

    unit_exponents = [
        tuple(1 if index == position else 0 for index in range(5))
        for position in range(5)
    ]
    h_map: VectorMap = {
        unit_exponents[index]: residual_vectors[:, index]
        for index in range(5)
    }
    h_polynomials = [
        vector_to_polynomial(residual_vectors[:, index])
        for index in range(5)
    ]

    rows = [
        [
            generating_weights(4, order).get(exponent, sp.Rational(0))
            for exponent in BASIS
        ]
        for order in range(MAX_ORDER + 1)
    ]
    first_rows = sp.Matrix(rows[:12])
    pivot_columns = list(first_rows.rref()[1])
    pivot_inverse = first_rows[:, pivot_columns].inv()
    row_coordinates = [
        (
            sp.Matrix([[row[index] for index in pivot_columns]])
            * pivot_inverse
        )
        for row in rows
    ]
    bilinear_matrices = []
    for order in range(MAX_ORDER):
        weights = generating_weights(8, order)
        matrix = sp.zeros(25)
        for left, left_exponent in enumerate(BASIS):
            for right in range(left, 25):
                right_exponent = BASIS[right]
                exponent = tuple(
                    left_exponent[index] + right_exponent[index]
                    for index in range(4)
                )
                matrix[left, right] = weights.get(
                    exponent,
                    sp.Rational(0),
                )
                matrix[right, left] = matrix[left, right]
        bilinear_matrices.append(matrix)

    # The particular second-order correction K_0(h), before its tangent
    # correction, is obtained from the first twelve moment equations.
    q_maps: list[ScalarMap] = [{}]
    for order in range(1, 12):
        q_maps.append(
            add_scalar_maps(
                (
                    sp.Rational(order, 2),
                    bilinear_map(
                        h_map,
                        h_map,
                        bilinear_matrices[order - 1],
                    ),
                )
            )
        )
    degree_two = compositions(2)
    second_map: VectorMap = {}
    for exponent in degree_two:
        pivot_values = -pivot_inverse * coefficient_vector(q_maps, exponent)
        vector = sp.zeros(25, 1)
        for index, pivot in enumerate(pivot_columns):
            vector[pivot] = pivot_values[index]
        second_map[exponent] = vector

    # Add the polynomial tangent correction certified by the local checker.
    h_symbols = sp.symbols("h0:5")
    tangent_correction_polynomials = [
        sp.Poly(sp.sympify(value), *h_symbols)
        for value in stored["polynomial_tangent_correction_for_third_order"]
    ]
    for exponent in degree_two:
        monomial = sp.prod(
            h_symbols[index] ** exponent[index] for index in range(5)
        )
        tangent_coordinates = sp.Matrix(
            [
                polynomial.coeff_monomial(monomial)
                for polynomial in tangent_correction_polynomials
            ]
        )
        second_map[exponent] += tangent_basis * tangent_coordinates

    # Build H^2 and H^3 numerator coefficient maps once.
    h_squared_numerators: dict[Exponent5, seed.Polynomial] = {}
    for left_index, left_exponent in enumerate(unit_exponents):
        for right_index, right_exponent in enumerate(unit_exponents):
            add_numerator(
                h_squared_numerators,
                add_exponents(left_exponent, right_exponent),
                seed.multiply(
                    h_polynomials[left_index],
                    h_polynomials[right_index],
                ),
            )
    h_cubed_numerators: dict[Exponent5, seed.Polynomial] = {}
    for squared_exponent, squared_numerator in h_squared_numerators.items():
        for index, unit_exponent in enumerate(unit_exponents):
            add_numerator(
                h_cubed_numerators,
                add_exponents(squared_exponent, unit_exponent),
                seed.multiply(squared_numerator, h_polynomials[index]),
            )
    degree_twelve_weights = [
        generating_weights(12, order)
        for order in range(MAX_ORDER - 1)
    ]
    h_cubed_contractions = [
        contraction_map(
            h_cubed_numerators,
            degree_twelve_weights[order],
        )
        for order in range(MAX_ORDER - 1)
    ]

    # Solve the first twelve cubic equations for a particular L(h).
    third_maps: list[ScalarMap] = [{}]
    for order in range(1, 12):
        terms = [
            (
                sp.Rational(order),
                bilinear_map(
                    h_map,
                    second_map,
                    bilinear_matrices[order - 1],
                ),
            )
        ]
        if order >= 2:
            terms.append(
                (
                    sp.Rational(order * (order - 1), 6),
                    h_cubed_contractions[order - 2],
                )
            )
        third_maps.append(add_scalar_maps(*terms))
    degree_three = compositions(3)
    third_map: VectorMap = {}
    for exponent in degree_three:
        pivot_values = -pivot_inverse * coefficient_vector(
            third_maps,
            exponent,
        )
        vector = sp.zeros(25, 1)
        for index, pivot in enumerate(pivot_columns):
            vector[pivot] = pivot_values[index]
        third_map[exponent] = vector

    # Precompute H^2 K and H^4 numerator coefficient maps.
    second_polynomials = {
        exponent: vector_to_polynomial(vector)
        for exponent, vector in second_map.items()
    }
    h2_second_numerators: dict[Exponent5, seed.Polynomial] = {}
    for h2_exponent, h2_numerator in h_squared_numerators.items():
        for second_exponent, second_polynomial in second_polynomials.items():
            add_numerator(
                h2_second_numerators,
                add_exponents(h2_exponent, second_exponent),
                seed.multiply(h2_numerator, second_polynomial),
            )
    h_fourth_numerators: dict[Exponent5, seed.Polynomial] = {}
    for left_exponent, left_numerator in h_squared_numerators.items():
        for right_exponent, right_numerator in h_squared_numerators.items():
            add_numerator(
                h_fourth_numerators,
                add_exponents(left_exponent, right_exponent),
                seed.multiply(left_numerator, right_numerator),
            )
    h2_second_contractions = [
        contraction_map(
            h2_second_numerators,
            degree_twelve_weights[order],
        )
        for order in range(MAX_ORDER - 1)
    ]
    degree_sixteen_weights = [
        generating_weights(16, order)
        for order in range(MAX_ORDER - 2)
    ]
    h_fourth_contractions = [
        contraction_map(
            h_fourth_numerators,
            degree_sixteen_weights[order],
        )
        for order in range(MAX_ORDER - 2)
    ]

    # Assemble the known fourth-order coefficient maps.
    fourth_maps: list[ScalarMap] = [{}]
    for order in range(1, MAX_ORDER + 1):
        terms = [
            (
                sp.Rational(order),
                bilinear_map(
                    h_map,
                    third_map,
                    bilinear_matrices[order - 1],
                ),
            ),
            (
                sp.Rational(order, 2),
                bilinear_map(
                    second_map,
                    second_map,
                    bilinear_matrices[order - 1],
                ),
            ),
        ]
        if order >= 2:
            terms.append(
                (
                    sp.Rational(order * (order - 1), 2),
                    h2_second_contractions[order - 2],
                )
            )
        if order >= 3:
            terms.append(
                (
                    sp.Rational(
                        order * (order - 1) * (order - 2),
                        24,
                    ),
                    h_fourth_contractions[order - 3],
                )
            )
        fourth_maps.append(add_scalar_maps(*terms))

    fourth_obstructions: list[ScalarMap] = []
    for order in range(12, MAX_ORDER + 1):
        terms = [(sp.Rational(1), fourth_maps[order])]
        terms.extend(
            (
                -row_coordinates[order][index],
                fourth_maps[index],
            )
            for index in range(12)
        )
        fourth_obstructions.append(add_scalar_maps(*terms))

    # The coefficient of the new tangent correction is the same linear
    # matrix M(h) that appeared at cubic order.
    tangent_linear: list[list[list[sp.Rational]]] = []
    for order in range(12, MAX_ORDER + 1):
        order_rows: list[list[sp.Rational]] = []
        for tangent_index in range(13):
            coefficients = []
            tangent_vector = tangent_basis[:, tangent_index]
            for h_index in range(5):
                value = (
                    order
                    * (
                        residual_vectors[:, h_index].T
                        * bilinear_matrices[order - 1]
                        * tangent_vector
                    )[0]
                )
                for prefix in range(12):
                    if prefix == 0:
                        continue
                    value -= (
                        row_coordinates[order][prefix]
                        * prefix
                        * (
                            residual_vectors[:, h_index].T
                            * bilinear_matrices[prefix - 1]
                            * tangent_vector
                        )[0]
                    )
                coefficients.append(sp.factor(value))
            order_rows.append(coefficients)
        tangent_linear.append(order_rows)

    degree_four = compositions(4)
    column_pairs = [
        (tangent_index, exponent)
        for tangent_index in range(13)
        for exponent in degree_three
    ]
    system_rows: list[list[str]] = []
    for obstruction_index, obstruction in enumerate(fourth_obstructions):
        for target_exponent in degree_four:
            row: list[sp.Rational] = []
            for tangent_index, cubic_exponent in column_pairs:
                coefficient = sp.Rational(0)
                for h_index, unit_exponent in enumerate(unit_exponents):
                    if add_exponents(cubic_exponent, unit_exponent) == target_exponent:
                        coefficient += tangent_linear[
                            obstruction_index
                        ][tangent_index][h_index]
                row.append(coefficient)
            right = -obstruction.get(target_exponent, sp.Rational(0))
            system_rows.append([str(value) for value in row + [right]])

    assert len(system_rows) == TAIL_ORDER_COUNT * 70
    assert len(system_rows[0]) == 13 * 35 + 1
    augmented = fmpq_mat(system_rows)
    reduced, rank = augmented.rref()
    unknown_count = len(column_pairs)
    inconsistent_rows = []
    pivots = []
    for row in range(rank):
        pivot = next(
            (
                column
                for column in range(unknown_count)
                if reduced[row, column]
            ),
            None,
        )
        if pivot is None:
            if reduced[row, unknown_count]:
                inconsistent_rows.append(row)
            continue
        pivots.append(pivot)

    sample_values = {
        "generic_1_2_3_4_5": (1, 2, 3, 4, 5),
        "axis_h0": (1, 0, 0, 0, 0),
        "axis_h1": (0, 1, 0, 0, 0),
        "axis_h2": (0, 0, 1, 0, 0),
        "axis_h3": (0, 0, 0, 1, 0),
        "axis_h4": (0, 0, 0, 0, 1),
    }
    sample_reports = {}
    for label, values in sample_values.items():
        sample_matrix = sp.Matrix(
            [
                [
                    sum(
                        tangent_linear[order][tangent_index][h_index]
                        * values[h_index]
                        for h_index in range(5)
                    )
                    for tangent_index in range(13)
                ]
                for order in range(TAIL_ORDER_COUNT)
            ]
        )
        sample_constant = sp.Matrix(
            [
                sum(
                    coefficient
                    * sp.prod(
                        values[index] ** exponent[index]
                        for index in range(5)
                    )
                    for exponent, coefficient in obstruction.items()
                )
                for obstruction in fourth_obstructions
            ]
        )
        sample_reports[label] = {
            "correction_rank": sample_matrix.rank(),
            "augmented_rank": sample_matrix.row_join(
                -sample_constant
            ).rank(),
        }

    residual_symbols = sp.symbols("h0:5")
    symbolic_tangent_matrix = sp.Matrix(
        [
            [
                sum(
                    tangent_linear[order][tangent_index][h_index]
                    * residual_symbols[h_index]
                    for h_index in range(5)
                )
                for tangent_index in range(13)
            ]
            for order in range(TAIL_ORDER_COUNT)
        ]
    )
    symbolic_fourth_constant = sp.Matrix(
        [
            sum(
                coefficient
                * sp.prod(
                    residual_symbols[index] ** exponent[index]
                    for index in range(5)
                )
                for exponent, coefficient in obstruction.items()
            )
            for obstruction in fourth_obstructions
        ]
    )
    generic_substitution = dict(
        zip(residual_symbols, sample_values["generic_1_2_3_4_5"])
    )
    generic_matrix = symbolic_tangent_matrix.subs(generic_substitution)
    generic_columns = list(generic_matrix.rref()[1])
    assert len(generic_columns) == 2
    generic_rows = list(generic_matrix[:, generic_columns].T.rref()[1])
    assert len(generic_rows) == 2
    obstruction_row = next(
        row
        for row in range(TAIL_ORDER_COUNT)
        if row not in generic_rows
        and generic_matrix[
            generic_rows + [row],
            generic_columns,
        ].row_join(
            -symbolic_fourth_constant[
                generic_rows + [row],
                :,
            ].subs(generic_substitution)
        ).rank()
        == 3
    )
    selected_rows = generic_rows + [obstruction_row]
    selected_matrix = symbolic_tangent_matrix[
        selected_rows,
        generic_columns,
    ].row_join(
        -symbolic_fourth_constant[selected_rows, :]
    )
    chart_obstruction = sp.factor(selected_matrix.det())
    chart_pivot = sp.factor(
        symbolic_tangent_matrix[
            generic_rows,
            generic_columns,
        ].det()
    )
    assert chart_obstruction.subs(generic_substitution) != 0

    # Restore the full eleven-dimensional freedom in the cubic lift at one
    # exact generic residual direction.  The fourth-order compatibility
    # equations become twelve quadrics in those eleven free coordinates.
    generic_values = sample_values["generic_1_2_3_4_5"]

    def evaluate_vector_map(
        polynomial: VectorMap,
        values: tuple[int, int, int, int, int],
    ) -> sp.Matrix:
        result = sp.zeros(25, 1)
        for exponent, vector in polynomial.items():
            result += vector * sp.prod(
                values[index] ** exponent[index] for index in range(5)
            )
        return result

    generic_h = residual_vectors * sp.Matrix(generic_values)
    generic_second = evaluate_vector_map(second_map, generic_values)
    generic_tangent_matrix = symbolic_tangent_matrix.subs(
        generic_substitution
    )
    generic_tangent_kernel = sp.Matrix.hstack(
        *generic_tangent_matrix.nullspace()
    )
    assert generic_tangent_kernel.shape == (13, 11)
    second_variation = tangent_basis * generic_tangent_kernel

    third_variation_base = sp.zeros(12, 11)
    for order in range(1, 12):
        third_variation_base[order, :] = (
            order
            * generic_h.T
            * bilinear_matrices[order - 1]
            * second_variation
        )
    third_variation = sp.zeros(25, 11)
    third_pivot_variation = -pivot_inverse * third_variation_base
    for index, pivot in enumerate(pivot_columns):
        third_variation[pivot, :] = third_pivot_variation[index, :]

    generic_h_squared = seed.monomial(seed.ZERO)
    generic_h_polynomial = vector_to_polynomial(generic_h)
    generic_h_squared = seed.multiply(
        generic_h_polynomial,
        generic_h_polynomial,
    )
    h2_functionals = []
    for order in range(MAX_ORDER - 1):
        weights = degree_twelve_weights[order]
        functional = []
        for basis_exponent in BASIS:
            numerator = seed.multiply(
                generic_h_squared,
                seed.monomial(basis_exponent),
            )
            functional.append(
                sum(
                    local.rational(coefficient)
                    * weights.get(exponent, 0)
                    for exponent, coefficient in numerator.items()
                )
            )
        h2_functionals.append(sp.Matrix(functional))

    raw_linear = [sp.zeros(1, 11)]
    raw_quadratic = [sp.zeros(11)]
    for order in range(1, MAX_ORDER + 1):
        bilinear = bilinear_matrices[order - 1]
        linear = (
            order * generic_h.T * bilinear * third_variation
            + order
            * generic_second.T
            * bilinear
            * second_variation
        )
        if order >= 2:
            linear += (
                sp.Rational(order * (order - 1), 2)
                * h2_functionals[order - 2].T
                * second_variation
            )
        raw_linear.append(linear)
        raw_quadratic.append(
            sp.Rational(order, 2)
            * second_variation.T
            * bilinear
            * second_variation
        )

    tail_linear = []
    tail_quadratic = []
    for order in range(12, MAX_ORDER + 1):
        linear = raw_linear[order].copy()
        quadratic = raw_quadratic[order].copy()
        for prefix in range(12):
            linear -= row_coordinates[order][prefix] * raw_linear[prefix]
            quadratic -= (
                row_coordinates[order][prefix] * raw_quadratic[prefix]
            )
        tail_linear.append(linear)
        tail_quadratic.append(quadratic)

    left_kernel = sp.Matrix.hstack(
        *generic_tangent_matrix.T.nullspace()
    ).T
    assert left_kernel.shape == (
        TAIL_ORDER_COUNT - 2,
        TAIL_ORDER_COUNT,
    )
    generic_fixed_constant = symbolic_fourth_constant.subs(
        generic_substitution
    )
    free_variables = sp.symbols("u0:11")
    free_vector = sp.Matrix(free_variables)
    fiber_equations = []
    for row in range(left_kernel.rows):
        equation = sp.Integer(0)
        for tail_index in range(TAIL_ORDER_COUNT):
            equation += left_kernel[row, tail_index] * (
                generic_fixed_constant[tail_index]
                + (tail_linear[tail_index] * free_vector)[0]
                + (
                    free_vector.T
                    * tail_quadratic[tail_index]
                    * free_vector
                )[0]
            )
        fiber_equations.append(sp.factor(equation))
    singular_equations = [
        str(equation).replace("**", "^")
        for equation in fiber_equations
        if equation
    ]
    fiber_script = "\n".join(
        [
            "ring r=0,(u0,u1,u2,u3,u4,u5,u6,u7,u8,u9,u10),dp;",
            f"ideal I={','.join(singular_equations)};",
            "ideal G=std(I);",
            'print("FIBER_DIM"); print(dim(G));',
            'print("FIBER_DEG"); print(deg(G));',
            'print("FIBER_SIZE"); print(size(G));',
            'print("FIBER_GB"); print(G); print("FIBER_END");',
        ]
    )
    singular = shutil.which("Singular")
    assert singular, "Singular is required for the generic fourth-order fiber"
    fiber_result = subprocess.run(
        [singular, "-q"],
        input=fiber_script,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    fiber_lines = [
        line.strip() for line in fiber_result.stdout.splitlines()
    ]

    def value_after(label: str) -> int:
        index = fiber_lines.index(label)
        return int(fiber_lines[index + 1])

    generic_fiber_report = {
        "direction": list(generic_values),
        "free_cubic_lift_parameters": 11,
        "compatibility_equation_count": len(singular_equations),
        "dimension": value_after("FIBER_DIM"),
        "degree": value_after("FIBER_DEG"),
        "groebner_basis_size": value_after("FIBER_SIZE"),
        "groebner_basis": fiber_lines[
            fiber_lines.index("FIBER_GB") + 1:
            fiber_lines.index("FIBER_END")
        ],
    }
    assert generic_fiber_report["dimension"] >= 0
    free_symbols = sp.symbols("u0:11")
    symbol_table = {str(symbol): symbol for symbol in free_symbols}
    groebner_polynomials = [
        sp.sympify(
            polynomial.rstrip(",").replace("^", "**"),
            locals=symbol_table,
        )
        for polynomial in generic_fiber_report["groebner_basis"]
    ]
    quadric = next(
        polynomial
        for polynomial in groebner_polynomials
        if sp.Poly(polynomial, *free_symbols).total_degree() == 2
    )
    quadric_variables = list(free_symbols[1:])
    quadratic_rank = (sp.hessian(quadric, quadric_variables) / 2).rank()
    assert quadratic_rank == 1
    quadratic_matrix = sp.hessian(quadric, quadric_variables) / 2
    linear_part = sp.Matrix(
        [
            sp.diff(quadric, variable).subs(
                {entry: 0 for entry in quadric_variables}
            )
            for variable in quadric_variables
        ]
    )
    assert all(
        (linear_part.T * null_vector)[0] == 0
        for null_vector in quadratic_matrix.nullspace()
    )
    zero_specialization = {
        variable: 0
        for variable in quadric_variables
        if variable != quadric_variables[0]
    }
    univariate_quadric = sp.Poly(
        quadric.subs(zero_specialization),
        quadric_variables[0],
    )
    discriminant = int(sp.discriminant(univariate_quadric.as_expr()))
    discriminant_factorization = sp.factorint(discriminant)
    square_class = 1
    for prime, exponent in discriminant_factorization.items():
        if exponent % 2:
            square_class *= prime
    assert square_class == 41
    generic_fiber_report["quadric_rank"] = quadratic_rank
    generic_fiber_report["discriminant"] = str(discriminant)
    generic_fiber_report["discriminant_square_class"] = square_class
    generic_fiber_report["geometric_components"] = 2
    generic_fiber_report["splitting_field"] = "Q(sqrt(41))"
    generic_fiber_report["has_Q_point"] = False

    # Choose one explicit point on the positive sqrt(41) component, then
    # reconstruct the complete H,K,L,M jet and verify its moments directly.
    univariate_coefficients = univariate_quadric.all_coeffs()
    root_u1 = sp.factor(
        (
            -univariate_coefficients[1]
            - sp.sqrt(discriminant)
        )
        / (2 * univariate_coefficients[0])
    )
    explicit_u = sp.zeros(11, 1)
    explicit_u[1] = root_u1
    linear_equation = next(
        polynomial
        for polynomial in groebner_polynomials
        if sp.Poly(polynomial, *free_symbols).total_degree() == 1
    )
    explicit_u[0] = sp.factor(
        sp.solve(
            linear_equation.subs(
                {
                    free_symbols[index]: explicit_u[index]
                    for index in range(1, 11)
                }
            ),
            free_symbols[0],
        )[0]
    )
    explicit_substitution = {
        free_symbols[index]: explicit_u[index]
        for index in range(11)
    }
    assert all(
        sp.factor(equation.subs(explicit_substitution)) == 0
        for equation in fiber_equations
    )

    explicit_constant = sp.Matrix(
        [
            generic_fixed_constant[index]
            + (tail_linear[index] * explicit_u)[0]
            + (
                explicit_u.T
                * tail_quadratic[index]
                * explicit_u
            )[0]
            for index in range(TAIL_ORDER_COUNT)
        ]
    )
    cubic_tangent_coordinates = sp.zeros(13, 1)
    cubic_pivot_solution = -generic_tangent_matrix[
        generic_rows,
        generic_columns,
    ].inv() * explicit_constant[generic_rows, :]
    for index, column in enumerate(generic_columns):
        cubic_tangent_coordinates[column] = sp.factor(
            cubic_pivot_solution[index]
        )

    generic_third = evaluate_vector_map(third_map, generic_values)
    explicit_second = sp.simplify(
        generic_second + second_variation * explicit_u
    )
    explicit_third = sp.simplify(
        generic_third
        + third_variation * explicit_u
        + tangent_basis * cubic_tangent_coordinates
    )

    h_cubed_polynomial = seed.multiply(
        generic_h_squared,
        generic_h_polynomial,
    )
    h_fourth_polynomial = seed.multiply(
        generic_h_squared,
        generic_h_squared,
    )

    def polynomial_contraction(
        polynomial: seed.Polynomial,
        degree: int,
        f_order: int,
    ) -> sp.Expr:
        weights = generating_weights(degree, f_order)
        return sp.factor(
            sum(
                sp.Rational(coefficient)
                * weights.get(exponent, sp.Rational(0))
                for exponent, coefficient in polynomial.items()
            )
        )

    fourth_known_prefix = []
    for order in range(12):
        if order == 0:
            fourth_known_prefix.append(sp.Rational(0))
            continue
        value = (
            order
            * (generic_h.T * bilinear_matrices[order - 1] * explicit_third)[0]
            + sp.Rational(order, 2)
            * (
                explicit_second.T
                * bilinear_matrices[order - 1]
                * explicit_second
            )[0]
        )
        if order >= 2:
            value += (
                sp.Rational(order * (order - 1), 2)
                * (h2_functionals[order - 2].T * explicit_second)[0]
            )
        if order >= 3:
            value += (
                sp.Rational(order * (order - 1) * (order - 2), 24)
                * polynomial_contraction(
                    h_fourth_polynomial,
                    16,
                    order - 3,
                )
            )
        fourth_known_prefix.append(sp.factor(value))
    fourth_pivot_values = -pivot_inverse * sp.Matrix(
        fourth_known_prefix
    )
    explicit_fourth = sp.zeros(25, 1)
    for index, pivot in enumerate(pivot_columns):
        explicit_fourth[pivot] = sp.factor(fourth_pivot_values[index])

    for order in range(MAX_ORDER + 1):
        assert sp.factor((sp.Matrix([rows[order]]) * generic_h)[0]) == 0
        second_value = (sp.Matrix([rows[order]]) * explicit_second)[0]
        if order >= 1:
            second_value += (
                sp.Rational(order, 2)
                * (
                    generic_h.T
                    * bilinear_matrices[order - 1]
                    * generic_h
                )[0]
            )
        assert sp.factor(second_value) == 0

        third_value = (sp.Matrix([rows[order]]) * explicit_third)[0]
        if order >= 1:
            third_value += (
                order
                * (
                    generic_h.T
                    * bilinear_matrices[order - 1]
                    * explicit_second
                )[0]
            )
        if order >= 2:
            third_value += (
                sp.Rational(order * (order - 1), 6)
                * polynomial_contraction(
                    h_cubed_polynomial,
                    12,
                    order - 2,
                )
            )
        assert sp.factor(third_value) == 0

        fourth_value = (sp.Matrix([rows[order]]) * explicit_fourth)[0]
        if order >= 1:
            fourth_value += (
                order
                * (
                    generic_h.T
                    * bilinear_matrices[order - 1]
                    * explicit_third
                )[0]
                + sp.Rational(order, 2)
                * (
                    explicit_second.T
                    * bilinear_matrices[order - 1]
                    * explicit_second
                )[0]
            )
        if order >= 2:
            fourth_value += (
                sp.Rational(order * (order - 1), 2)
                * (h2_functionals[order - 2].T * explicit_second)[0]
            )
        if order >= 3:
            fourth_value += (
                sp.Rational(order * (order - 1) * (order - 2), 24)
                * polynomial_contraction(
                    h_fourth_polynomial,
                    16,
                    order - 3,
                )
            )
        assert sp.factor(fourth_value) == 0

    def sparse_vector(vector: sp.Matrix) -> dict[str, str]:
        return {
            "^".join(map(str, BASIS[index])): str(sp.factor(value))
            for index, value in enumerate(vector)
            if value
        }

    generic_fiber_report["explicit_component_point"] = {
        "free_second_lift_coordinates": [
            str(sp.factor(value)) for value in explicit_u
        ],
        "cubic_tangent_coordinates": [
            str(sp.factor(value)) for value in cubic_tangent_coordinates
        ],
        "jet": {
            "H": sparse_vector(generic_h),
            "K": sparse_vector(explicit_second),
            "L": sparse_vector(explicit_third),
            "M": sparse_vector(explicit_fourth),
        },
        "verified_moment_orders": list(range(MAX_ORDER + 1)),
        "all_order_completion": "combined 42-row beta-tail",
    }

    result = {
        "format": "two-pair-counterexample-fourth-order-v1",
        "residual_parameter_count": 5,
        "fourth_order_equation_count": len(system_rows),
        "cubic_correction_coefficient_count": unknown_count,
        "coefficient_rank": len(pivots),
        "augmented_rank": rank,
        "consistent": not inconsistent_rows,
        "inconsistent_rref_row_count": len(inconsistent_rows),
        "fixed_polynomial_section_sample_fibers": sample_reports,
        "fixed_section_generic_chart": {
            "tail_row_indices": selected_rows,
            "tangent_column_indices": generic_columns,
            "pivot": str(chart_pivot),
            "fourth_obstruction": str(chart_obstruction),
        },
        "generic_direction_full_cubic_freedom": generic_fiber_report,
        "moment_tail_orders": list(range(12, MAX_ORDER + 1)),
        "all_order_master_tail": {
            "formula": (
                "E(A4 F^r)+r E(A8 F^(r-1))"
                "+r(r-1)/2 E(A12 F^(r-2))"
                "+r(r-1)(r-2)/24 E(A16 F^(r-3))"
            ),
            "symbolic_beta_numerator_degree_bound": 41,
            "basis_orders": list(range(12, 54)),
            "basis_rank": master_basis_rank,
            "replay_orders": list(range(12, MASTER_REPLAY_ORDER + 1)),
            "replay_rank": master_replay_rank,
        },
        "source_local_artifact": str(LOCAL_ARTIFACT.relative_to(ROOT)),
    }
    if not inconsistent_rows:
        solution = [sp.Rational(0) for _ in range(unknown_count)]
        for row, pivot in enumerate(pivots):
            solution[pivot] = sp.Rational(str(reduced[row, unknown_count]))
        correction = []
        for tangent_index in range(13):
            terms = {}
            for exponent_index, exponent in enumerate(degree_three):
                value = solution[
                    tangent_index * len(degree_three) + exponent_index
                ]
                if value:
                    terms["^".join(map(str, exponent))] = str(value)
            correction.append(terms)
        result["solution_nonzero_coefficient_count"] = sum(
            len(polynomial) for polynomial in correction
        )
        result["universal_cubic_tangent_correction"] = correction
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")

    print(
        "PASS SIC2C4 fourth: assembled exact "
        f"{len(system_rows)}x{unknown_count} universal correction system"
    )
    if inconsistent_rows:
        print(
            "PASS SIC2C4 fourth: universal polynomial lift is obstructed "
            f"by {len(inconsistent_rows)} reduced rows"
        )
    else:
        print(
            "PASS SIC2C4 fourth: every reduced quadratic direction "
            "lifts through fourth order"
        )
    print(f"PASS SIC2C4 fourth: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
