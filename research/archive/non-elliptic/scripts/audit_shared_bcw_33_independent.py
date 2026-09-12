#!/usr/bin/env python3
"""Dependency-free replay of the shared-factor 33-variable BCW artifact.

The generator uses SymPy to verify a factor trace frozen from a deterministic
search.  This audit uses only Fraction-valued sparse dictionaries.  It replays
every stored source factor exposure and elementary target shear, reconstructs
the homogeneous cubic map, and checks the transported three-point collision.
"""

from __future__ import annotations

from fractions import Fraction as Q
from itertools import permutations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "artifacts" / "generated-results" / "shared_bcw_33_counterexample.json"
Poly = dict[tuple[int, ...], Q]


def clean(poly: Poly) -> Poly:
    return {exponents: coefficient for exponents, coefficient in poly.items() if coefficient}


def add(*polys: Poly) -> Poly:
    result: Poly = {}
    for poly in polys:
        for exponents, coefficient in poly.items():
            result[exponents] = result.get(exponents, Q(0)) + coefficient
    return clean(result)


def scale(poly: Poly, coefficient: Q | int) -> Poly:
    return clean({exponents: Q(coefficient) * value for exponents, value in poly.items()})


def multiply(left: Poly, right: Poly) -> Poly:
    result: Poly = {}
    for first, first_coefficient in left.items():
        for second, second_coefficient in right.items():
            exponents = tuple(a + b for a, b in zip(first, second))
            result[exponents] = result.get(exponents, Q(0)) + first_coefficient * second_coefficient
    return clean(result)


def power(poly: Poly, exponent: int) -> Poly:
    dimension = len(next(iter(poly)))
    result: Poly = {(0,) * dimension: Q(1)}
    base = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        remaining //= 2
    return result


def variable(dimension: int, index: int) -> Poly:
    exponents = [0] * dimension
    exponents[index] = 1
    return {tuple(exponents): Q(1)}


def constant(dimension: int, value: Q | int) -> Poly:
    return {(0,) * dimension: Q(value)} if value else {}


def pad(poly: Poly, amount: int) -> Poly:
    return {exponents + (0,) * amount: coefficient for exponents, coefficient in poly.items()}


def monomial(dimension: int, exponents: tuple[int, ...]) -> Poly:
    return {exponents + (0,) * (dimension - len(exponents)): Q(1)}


def derivative(poly: Poly, index: int) -> Poly:
    result: Poly = {}
    for exponents, coefficient in poly.items():
        if exponents[index]:
            reduced = list(exponents)
            reduced[index] -= 1
            result[tuple(reduced)] = coefficient * exponents[index]
    return result


def determinant3(matrix: list[list[Poly]]) -> Poly:
    dimension = len(next(iter(matrix[0][0])))
    result: Poly = {}
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        term = constant(dimension, -1 if inversions % 2 else 1)
        for row in range(3):
            term = multiply(term, matrix[row][permutation[row]])
        result = add(result, term)
    return result


def evaluate(poly: Poly, point: list[Q]) -> Q:
    total = Q(0)
    for exponents, coefficient in poly.items():
        value = coefficient
        for coordinate, exponent in zip(point, exponents):
            value *= coordinate**exponent
        total += value
    return total


def dense(support: list[list[int]], dimension: int) -> tuple[int, ...]:
    exponents = [0] * dimension
    for index, exponent in support:
        exponents[index] = exponent
    return tuple(exponents)


def main() -> None:
    stored = json.loads(ARTIFACT.read_text())
    assert stored["format"] == "shared-bcw-sparse-cubic-homogeneous-map-v1"
    assert stored["degree_reduced_dimension"] == 16
    assert stored["dimension"] == 33
    assert len(stored["degree_lowering_steps"]) == 17

    dimension = 3
    x, y, z = (variable(dimension, index) for index in range(3))
    v = add(constant(dimension, 1), scale(multiply(x, y), 2))
    long_map = [
        add(
            multiply(power(v, 3), z),
            scale(
                multiply(multiply(power(y, 2), v), add(constant(dimension, 2), scale(multiply(x, y), 3))),
                4,
            ),
        ),
        add(
            y,
            scale(multiply(multiply(x, power(v, 2)), z), 3),
            scale(
                multiply(multiply(x, power(y, 2)), add(constant(dimension, 2), scale(multiply(x, y), 3))),
                12,
            ),
        ),
        add(scale(x, -1), scale(multiply(power(x, 2), y), 3), multiply(power(x, 3), z)),
    ]
    jacobian = [[derivative(long_map[row], column) for column in range(3)] for row in range(3)]
    assert determinant3(jacobian) == {(0, 0, 0): Q(1)}
    polynomials = [scale(long_map[2], -1), long_map[1], long_map[0]]

    points = [
        [Q(0), Q(0), -Q(1, 8)],
        [Q(1), -Q(3, 4), Q(13, 4)],
        [-Q(1), Q(3, 4), Q(13, 4)],
    ]
    for point in points:
        assert [evaluate(poly, point) for poly in long_map] == [-Q(1, 8), Q(0), Q(0)]

    registry: dict[tuple[int, ...], int] = {}
    introduced = 0
    selected_degrees: list[int] = []
    for step in stored["degree_lowering_steps"]:
        maximum_degree = max(sum(exponents) for poly in polynomials for exponents in poly)
        component = step["component"]
        removed = dense(step["removed_monomial"], dimension)
        coefficient = Q(step["removed_coefficient"])
        first = dense(step["first_factor"], dimension)
        second = dense(step["second_factor"], dimension)
        assert step["selected_degree"] == maximum_degree == sum(removed)
        assert tuple(a + b for a, b in zip(first, second)) == removed
        assert polynomials[component].get(removed) == coefficient
        assert registry.get(first) != component and registry.get(second) != component

        missing: list[tuple[int, ...]] = []
        for factor in (first, second):
            if factor not in registry and factor not in missing:
                missing.append(factor)
        records = step["new_factors"]
        assert len(records) == len(missing)
        old_dimension = dimension
        additions_by_point = [[] for _ in points]
        for offset, (factor, record) in enumerate(zip(missing, records)):
            assert dense(record["factor"], old_dimension) == factor
            assert record["new_variable"] == dimension
            assert record["output"] == len(polynomials) + offset
            for additions, point in zip(additions_by_point, points):
                additions.append(-evaluate(monomial(old_dimension, factor), point))
            dimension += 1

        added = len(missing)
        introduced += added
        polynomials = [pad(poly, added) for poly in polynomials]
        registry = {factor + (0,) * added: output for factor, output in registry.items()}
        padded_missing = [factor + (0,) * added for factor in missing]
        for padded_factor, record in zip(padded_missing, records):
            new_variable = variable(dimension, record["new_variable"])
            factor_poly = monomial(dimension, padded_factor)
            registry[padded_factor] = len(polynomials)
            polynomials.append(add(new_variable, factor_poly))
        for point, additions in zip(points, additions_by_point):
            point.extend(additions)

        first += (0,) * added
        second += (0,) * added
        outputs = [registry[first], registry[second]]
        assert outputs == step["factor_outputs"]
        assert component not in outputs
        polynomials[component] = add(
            polynomials[component],
            scale(multiply(polynomials[outputs[0]], polynomials[outputs[1]]), -coefficient),
        )
        registry = {factor: output for factor, output in registry.items() if output != component}
        selected_degrees.append(maximum_degree)

    assert selected_degrees == [7, 6, 6] + [5] * 4 + [4] * 10
    assert introduced == 13
    assert dimension == 16
    assert max(sum(exponents) for poly in polynomials for exponents in poly) == 3
    for index, poly in enumerate(polynomials):
        linear = {exponents: coefficient for exponents, coefficient in poly.items() if sum(exponents) == 1}
        assert linear == {tuple(1 if j == index else 0 for j in range(dimension)): Q(1)}
        assert poly.get((0,) * dimension, Q(0)) == 0

    expected_target = [Q(0), Q(0), -Q(1, 8)] + [Q(0)] * 13
    for point in points:
        assert [evaluate(poly, point) for poly in polynomials] == expected_target

    quadratic = [{e: c for e, c in poly.items() if sum(e) == 2} for poly in polynomials]
    cubic = [{e: c for e, c in poly.items() if sum(e) == 3} for poly in polynomials]
    final_dimension = 2 * dimension + 1
    t_index = final_dimension - 1
    h_components: list[Poly] = [{} for _ in range(final_dimension)]
    for index, (quadratic_part, cubic_part) in enumerate(zip(quadratic, cubic)):
        exponents = [0] * final_dimension
        exponents[dimension + index], exponents[t_index] = 1, 2
        h_components[index][tuple(exponents)] = Q(1)
        for old, coefficient in quadratic_part.items():
            exponents = list(old) + [0] * (dimension + 1)
            exponents[t_index] = 1
            h_components[index][tuple(exponents)] = coefficient
        for old, coefficient in cubic_part.items():
            h_components[dimension + index][old + (0,) * (dimension + 1)] = -coefficient

    stored_h = [
        {dense(term["monomial"], final_dimension): Q(term["coefficient"]) for term in component}
        for component in stored["H"]
    ]
    assert h_components == stored_h
    assert all(sum(exponents) == 3 for component in h_components for exponents in component)

    final_points: list[list[Q]] = []
    for point in points:
        cubic_value = [evaluate(poly, point) for poly in cubic]
        final_points.append(point + cubic_value + [Q(1)])
    assert final_points == [[Q(value) for value in point] for point in stored["collision_points"]]

    common_image = [Q(value) for value in stored["common_image"]]
    for point in final_points:
        image = [coordinate + evaluate(component, point) for coordinate, component in zip(point, h_components)]
        assert image == common_image
    assert len({tuple(point) for point in final_points}) == 3

    print("PASS (stdlib) shared BCW: replayed 17 cancellations with 13 variables")
    print("PASS (stdlib) shared BCW: reconstructed the 33-variable cubic map")
    print("PASS (stdlib) shared BCW: verified three stored points and their common image")


if __name__ == "__main__":
    main()
