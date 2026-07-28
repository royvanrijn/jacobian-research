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
from fractions import Fraction
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


def contraction_map(
    numerators: dict[Exponent5, seed.Polynomial],
    f_power: seed.Polynomial,
) -> ScalarMap:
    result: ScalarMap = {}
    for exponent, numerator in numerators.items():
        value = local.scalar_contraction(seed.multiply(numerator, f_power))
        if value:
            result[exponent] = value
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

    f, _q, _generators = seed.witness()
    powers = [seed.monomial(seed.ZERO)]
    for _ in range(25):
        powers.append(seed.multiply(powers[-1], f))
    rows = [local.linear_row(power) for power in powers]
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
    bilinear_matrices = [
        local.bilinear_matrix(powers[order], order)
        for order in range(25)
    ]

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
    h_cubed_contractions = [
        contraction_map(h_cubed_numerators, powers[order])
        for order in range(24)
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
        contraction_map(h2_second_numerators, powers[order])
        for order in range(24)
    ]
    h_fourth_contractions = [
        contraction_map(h_fourth_numerators, powers[order])
        for order in range(23)
    ]

    # Assemble the known fourth-order coefficient maps.
    fourth_maps: list[ScalarMap] = [{}]
    for order in range(1, 26):
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
    for order in range(12, 26):
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
    for order in range(12, 26):
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

    assert len(system_rows) == 14 * 70
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

    result = {
        "format": "two-pair-counterexample-fourth-order-v1",
        "residual_parameter_count": 5,
        "fourth_order_equation_count": len(system_rows),
        "cubic_correction_coefficient_count": unknown_count,
        "augmented_rank": rank,
        "consistent": not inconsistent_rows,
        "inconsistent_rref_row_count": len(inconsistent_rows),
        "moment_tail_orders": list(range(12, 26)),
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
