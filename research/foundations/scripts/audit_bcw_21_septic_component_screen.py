#!/usr/bin/env python3
"""Screen base septics by the five vertical-ideal components.

Let E be the coefficient of X_20 in the stored Lie derivation restricted to
B=Q[X_0,...,X_13].  For a base septic P_0, E(P_0) modulo the vertical
coefficient ideal is the first necessary obstruction to lifting P_0 to
P_0+X_20*R.

The five substitutions below are the minimal components certified by
audit_bcw_21_vertical_ideal.sing.  Nonvanishing on any component proves that
the obstruction is nonzero over Q.  Modulo the good prime 1000003, this audit
checks all 77520 degree-seven base monomials and excludes 71588 of them.
It then stacks all five restrictions in each of 29 torus sectors. Unique-row
peeling plus exact modular core ranks give total rank 61060. Thus over Q the
subspace of arbitrary base septics surviving the radical-level screen has
dimension at most 16460. Embedded torsion and higher lifts remain open.
"""

from __future__ import annotations

from collections import defaultdict, deque
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "essential_bcw_21_counterexample.json"
)
PRIME = 1_000_003
DIMENSION = 14
WEIGHT_1 = (-1, 1, 2, 1, 1, 0, 1, 2, 0, -2, 0, 1, 1, 2)
WEIGHT_2 = (-4, 2, 5, 2, 2, -1, 2, 5, -1, -7, -1, 2, 2, 5)

Exponent = tuple[int, ...]
Polynomial = dict[Exponent, int]
StackedPolynomial = dict[tuple[int, ...], int]


def add_term(poly: Polynomial, exponent: Exponent, coefficient: int) -> None:
    updated = (poly.get(exponent, 0) + coefficient) % PRIME
    if updated:
        poly[exponent] = updated
    else:
        poly.pop(exponent, None)


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    answer: Polynomial = {}
    for alpha, coefficient_alpha in left.items():
        for beta, coefficient_beta in right.items():
            exponent = tuple(a + b for a, b in zip(alpha, beta))
            add_term(answer, exponent, coefficient_alpha * coefficient_beta)
    return answer


def variable(index: int, coefficient: int = 1) -> Polynomial:
    exponent = [0] * DIMENSION
    exponent[index] = 1
    return {tuple(exponent): coefficient % PRIME}


def decode_first_base_derivation(
    source: dict[str, object],
) -> list[Polynomial]:
    components = []
    for component in source["H"][:DIMENSION]:
        decoded: Polynomial = {}
        for term in component:
            s_power = next(
                (
                    power
                    for index, power in term["monomial"]
                    if index == 20
                ),
                0,
            )
            if s_power != 1:
                continue
            exponent = [0] * DIMENSION
            for index, power in term["monomial"]:
                if index < DIMENSION:
                    exponent[index] = power
            coefficient = Fraction(term["coefficient"])
            modular = (
                coefficient.numerator
                * pow(coefficient.denominator, -1, PRIME)
                % PRIME
            )
            add_term(decoded, tuple(exponent), modular)
        components.append(decoded)
    return components


def component_substitutions() -> list[list[Polynomial]]:
    zero: Polynomial = {}

    def coordinate_plane(zero_indices: set[int]) -> list[Polynomial]:
        return [
            zero if index in zero_indices else variable(index)
            for index in range(DIMENSION)
        ]

    substitutions = [
        coordinate_plane({1, 2, 6}),
        coordinate_plane({0, 5, 8}),
    ]

    oblique = coordinate_plane({0})
    oblique[5] = variable(8, -3)
    oblique[1] = {}
    inverse_seven = pow(7, -1, PRIME)
    for index, coefficient in ((3, 6), (4, 2), (6, 3)):
        for exponent, value in variable(
            index, coefficient * inverse_seven
        ).items():
            add_term(oblique[1], exponent, value)
    substitutions.append(oblique)

    substitutions.extend(
        [
            coordinate_plane({0, 3, 8}),
            coordinate_plane({0, 1}),
        ]
    )
    return substitutions


def substitute_monomial(
    exponent: Exponent, substitutions: list[Polynomial]
) -> Polynomial:
    answer = {(0,) * DIMENSION: 1}
    for index, power in enumerate(exponent):
        for _ in range(power):
            answer = multiply(answer, substitutions[index])
            if not answer:
                return {}
    return answer


def substitute(
    polynomial: Polynomial,
    substitutions: list[Polynomial],
    cache: dict[Exponent, Polynomial],
) -> Polynomial:
    answer: Polynomial = {}
    for exponent, coefficient in polynomial.items():
        if exponent not in cache:
            cache[exponent] = substitute_monomial(exponent, substitutions)
        for image_exponent, image_coefficient in cache[exponent].items():
            add_term(
                answer,
                image_exponent,
                coefficient * image_coefficient,
            )
    return answer


def derivative_column(
    exponent: Exponent, components: list[Polynomial]
) -> Polynomial:
    answer: Polynomial = {}
    for index, power in enumerate(exponent):
        if not power:
            continue
        derivative = list(exponent)
        derivative[index] -= 1
        for component_exponent, coefficient in components[index].items():
            target = tuple(
                left + right
                for left, right in zip(derivative, component_exponent)
            )
            add_term(answer, target, power * coefficient)
    return answer


def reduce_column(
    column: StackedPolynomial,
    pivots: dict[tuple[int, ...], StackedPolynomial],
) -> StackedPolynomial:
    while column:
        pivot = min(column)
        coefficient = column[pivot]
        if pivot not in pivots:
            return column
        for row, value in pivots[pivot].items():
            updated = (column.get(row, 0) - coefficient * value) % PRIME
            if updated:
                column[row] = updated
            else:
                column.pop(row, None)
    return {}


def column_rank(columns: list[StackedPolynomial]) -> int:
    pivots: dict[tuple[int, ...], StackedPolynomial] = {}
    for original in columns:
        residual = reduce_column(dict(original), pivots)
        if not residual:
            continue
        pivot = min(residual)
        inverse = pow(residual[pivot], -1, PRIME)
        pivots[pivot] = {
            row: coefficient * inverse % PRIME
            for row, coefficient in residual.items()
        }
    return len(pivots)


def peel_unique_rows(
    columns: list[StackedPolynomial],
) -> tuple[int, list[StackedPolynomial]]:
    row_columns: defaultdict[tuple[int, ...], set[int]] = defaultdict(set)
    for column_index, column in enumerate(columns):
        for row in column:
            row_columns[row].add(column_index)

    active = [True] * len(columns)
    queue = deque(row for row, support in row_columns.items() if len(support) == 1)
    peeled = 0
    while queue:
        row = queue.popleft()
        support = row_columns[row]
        if len(support) != 1:
            continue
        column_index = next(iter(support))
        if not active[column_index]:
            continue
        active[column_index] = False
        peeled += 1
        for incident_row in columns[column_index]:
            incident_support = row_columns[incident_row]
            incident_support.discard(column_index)
            if len(incident_support) == 1:
                queue.append(incident_row)

    core_rows = {row for row, support in row_columns.items() if support}
    core = [
        {
            row: coefficient
            for row, coefficient in column.items()
            if row in core_rows
        }
        for column, is_active in zip(columns, active)
        if is_active
    ]
    return peeled, core


def stacked_component_rank(
    components: list[Polynomial],
    substitutions: list[list[Polynomial]],
) -> tuple[int, int, int, int, int]:
    buckets: defaultdict[tuple[int, int], list[Exponent]] = defaultdict(list)
    for indices in combinations_with_replacement(range(DIMENSION), 7):
        exponent = [0] * DIMENSION
        for index in indices:
            exponent[index] += 1
        sector = (
            sum(WEIGHT_1[index] for index in indices),
            sum(WEIGHT_2[index] for index in indices),
        )
        buckets[sector].append(tuple(exponent))

    caches: list[dict[Exponent, Polynomial]] = [
        {} for _ in substitutions
    ]
    total_columns = 0
    total_peeled = 0
    total_core_rank = 0
    residual_sectors = 0
    for exponents in buckets.values():
        columns: list[StackedPolynomial] = []
        for exponent in exponents:
            raw = derivative_column(exponent, components)
            stacked: StackedPolynomial = {}
            for component_index, component in enumerate(substitutions):
                image = substitute(
                    raw, component, caches[component_index]
                )
                for image_exponent, coefficient in image.items():
                    stacked[(component_index,) + image_exponent] = coefficient
            columns.append(stacked)
        peeled, core = peel_unique_rows(columns)
        total_columns += len(columns)
        total_peeled += peeled
        total_core_rank += column_rank(core)
        residual_sectors += bool(core)
    return (
        len(buckets),
        total_columns,
        total_peeled,
        residual_sectors,
        total_core_rank,
    )


def main() -> None:
    source = json.loads(SOURCE.read_text())
    components = decode_first_base_derivation(source)
    substitutions = component_substitutions()
    caches: list[dict[Exponent, Polynomial]] = [
        {} for _ in substitutions
    ]

    component_counts = [0] * len(substitutions)
    detected = 0
    undetected = 0
    zero_columns = []
    total = 0

    for indices in combinations_with_replacement(range(DIMENSION), 7):
        exponent = [0] * DIMENSION
        for index in indices:
            exponent[index] += 1
        exponent_tuple = tuple(exponent)
        column = derivative_column(exponent_tuple, components)
        if not column:
            zero_columns.append(indices)

        flags = []
        for component_index, component in enumerate(substitutions):
            nonzero = bool(
                substitute(column, component, caches[component_index])
            )
            flags.append(nonzero)
            component_counts[component_index] += nonzero
        detected += any(flags)
        undetected += not any(flags)
        total += 1

    expected_zero = [
        tuple(sorted((3,) * (7 - power) + (5,) * power))
        for power in range(8)
    ]
    assert total == 77_520
    assert component_counts == [32_032, 25_740, 65_520, 25_740, 45_136]
    assert detected == 71_588
    assert undetected == 5_932
    assert zero_columns == expected_zero

    stacked = stacked_component_rank(components, substitutions)
    assert stacked == (29, 77_520, 28_764, 23, 32_296)
    total_rank = stacked[2] + stacked[4]
    assert total_rank == 61_060
    assert total - total_rank == 16_460

    print("PASS septic component screen: 77520 base monomials checked")
    print(f"PASS component detections: {component_counts}")
    print("THEOREM: 71588 support-one base septics have no s-adic lift")
    print("PASS stacked component map: rank 61060/77520 mod 1000003")
    print("THEOREM: radical-level surviving subspace has dimension at most 16460")
    print("OPEN: embedded torsion and higher s-adic lifts")


if __name__ == "__main__":
    main()
