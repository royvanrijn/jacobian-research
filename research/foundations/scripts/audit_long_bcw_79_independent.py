#!/usr/bin/env python3
"""Independent standard-library replay of Long's 79-variable BCW artifact.

The generator uses SymPy.  This audit deliberately does not: sparse
polynomials use dictionaries and every coefficient is a Fraction.  It starts
from Long's displayed three-variable map, independently selects all balanced
degree-lowering steps, reconstructs the cubic H, and checks the stored
79-variable collision.
"""

from fractions import Fraction as Q
from itertools import permutations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "artifacts" / "generated-results" / "long_bcw_79_counterexample.json"
Poly = dict[tuple[int, ...], Q]


def clean(poly: Poly) -> Poly:
    return {exponents: coefficient for exponents, coefficient in poly.items() if coefficient}


def add(*polys: Poly) -> Poly:
    out: Poly = {}
    for poly in polys:
        for exponents, coefficient in poly.items():
            out[exponents] = out.get(exponents, Q(0)) + coefficient
    return clean(out)


def scale(poly: Poly, coefficient: Q) -> Poly:
    return clean({exponents: coefficient * value for exponents, value in poly.items()})


def multiply(left: Poly, right: Poly) -> Poly:
    out: Poly = {}
    for a, ac in left.items():
        for b, bc in right.items():
            exponents = tuple(x + y for x, y in zip(a, b))
            out[exponents] = out.get(exponents, Q(0)) + ac * bc
    return clean(out)


def power(poly: Poly, exponent: int) -> Poly:
    dimension = len(next(iter(poly)))
    out = {tuple([0] * dimension): Q(1)}
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = multiply(out, base)
        base = multiply(base, base)
        n //= 2
    return out


def variable(dimension: int, index: int) -> Poly:
    exponents = [0] * dimension
    exponents[index] = 1
    return {tuple(exponents): Q(1)}


def constant(dimension: int, value: Q | int) -> Poly:
    return {tuple([0] * dimension): Q(value)} if value else {}


def pad(poly: Poly, amount: int = 2) -> Poly:
    return {exponents + (0,) * amount: coefficient for exponents, coefficient in poly.items()}


def monomial(dimension: int, exponents: tuple[int, ...], coefficient: Q | int = 1) -> Poly:
    return {exponents + (0,) * (dimension - len(exponents)): Q(coefficient)}


def derivative(poly: Poly, index: int) -> Poly:
    out: Poly = {}
    for exponents, coefficient in poly.items():
        if exponents[index]:
            reduced = list(exponents)
            reduced[index] -= 1
            out[tuple(reduced)] = coefficient * exponents[index]
    return out


def determinant3(matrix: list[list[Poly]]) -> Poly:
    dimension = len(next(iter(matrix[0][0])))
    out: Poly = {}
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        term = constant(dimension, -1 if inversions % 2 else 1)
        for row in range(3):
            term = multiply(term, matrix[row][permutation[row]])
        out = add(out, term)
    return out


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


def balanced_factor(exponents: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    first_degree = sum(exponents) // 2
    first = [0] * len(exponents)
    second = list(exponents)
    remaining = first_degree
    for index in range(len(exponents)):
        take = min(second[index], remaining)
        first[index] = take
        second[index] -= take
        remaining -= take
        if remaining == 0:
            break
    assert sum(first) == first_degree
    return tuple(first), tuple(second)


def encoded(poly: Poly) -> list[dict[str, object]]:
    return [
        {
            "coefficient": str(coefficient),
            "monomial": [[index, exponent] for index, exponent in enumerate(exponents) if exponent],
        }
        for exponents, coefficient in sorted(poly.items(), reverse=True)
        if coefficient
    ]


def main() -> None:
    stored = json.loads(ARTIFACT.read_text())
    assert stored["format"] == "long-bcw-sparse-cubic-homogeneous-map-v1"
    assert stored["dimension"] == 79
    assert stored["degree_reduced_dimension"] == 39
    assert len(stored["degree_lowering_steps"]) == 18

    dimension = 3
    x, y, z = (variable(dimension, index) for index in range(3))
    v = add(constant(dimension, 1), scale(multiply(x, y), Q(2)))
    long_map = [
        add(
            multiply(power(v, 3), z),
            scale(multiply(multiply(power(y, 2), v), add(constant(dimension, 2), scale(multiply(x, y), 3))), 4),
        ),
        add(
            y,
            scale(multiply(multiply(x, power(v, 2)), z), 3),
            scale(multiply(multiply(x, power(y, 2)), add(constant(dimension, 2), scale(multiply(x, y), 3))), 12),
        ),
        add(scale(x, -1), scale(multiply(power(x, 2), y), 3), multiply(power(x, 3), z)),
    ]
    jacobian = [[derivative(long_map[row], column) for column in range(3)] for row in range(3)]
    assert determinant3(jacobian) == {(0, 0, 0): Q(1)}
    polynomials = [scale(long_map[2], -1), long_map[1], long_map[0]]

    high_support = {
        str(degree): sum(
            1 for poly in polynomials for exponents in poly if sum(exponents) == degree
        )
        for degree in range(4, 8)
    }
    assert high_support == stored["high_degree_support"] == {"4": 3, "5": 2, "6": 2, "7": 1}

    points = [
        [Q(0), Q(0), -Q(1, 8)],
        [Q(1), -Q(3, 4), Q(13, 4)],
        [-Q(1), Q(3, 4), Q(13, 4)],
    ]
    for point in points:
        assert [evaluate(poly, point) for poly in long_map] == [-Q(1, 8), Q(0), Q(0)]

    selected_degrees: list[int] = []
    for stored_step in stored["degree_lowering_steps"]:
        maximum_degree = max(sum(exponents) for poly in polynomials for exponents in poly)
        selected = None
        for component, poly in enumerate(polynomials):
            candidates = [
                (exponents, coefficient)
                for exponents, coefficient in poly.items()
                if sum(exponents) == maximum_degree
            ]
            if candidates:
                exponents, coefficient = max(candidates, key=lambda item: item[0])
                selected = component, exponents, coefficient
                break
        assert selected is not None
        component, exponents, coefficient = selected
        a_exponents, b_exponents = balanced_factor(exponents)

        expected_step = {
            "selected_degree": maximum_degree,
            "component": component,
            "removed_coefficient": str(coefficient),
            "removed_monomial": [[i, e] for i, e in enumerate(exponents) if e],
            "a_exponents": [[i, e] for i, e in enumerate(a_exponents) if e],
            "b_exponents": [[i, e] for i, e in enumerate(b_exponents) if e],
            "new_variables": [dimension, dimension + 1],
        }
        assert expected_step == stored_step

        for point in points:
            a_value = coefficient
            b_value = Q(1)
            for coordinate, exponent in zip(point, a_exponents):
                a_value *= coordinate**exponent
            for coordinate, exponent in zip(point, b_exponents):
                b_value *= coordinate**exponent
            point.extend([-a_value, -b_value])

        dimension += 2
        polynomials = [pad(poly) for poly in polynomials]
        a = monomial(dimension, a_exponents, coefficient)
        b = monomial(dimension, b_exponents)
        new_y = variable(dimension, dimension - 2)
        new_z = variable(dimension, dimension - 1)
        polynomials[component] = add(
            polynomials[component], scale(multiply(add(new_y, a), add(new_z, b)), -1)
        )
        polynomials.extend([add(new_y, a), add(new_z, b)])
        selected_degrees.append(maximum_degree)

    assert selected_degrees == [7, 6, 6, 5, 5, 5] + [4] * 12
    assert dimension == 39
    assert max(sum(exponents) for poly in polynomials for exponents in poly) == 3

    expected_target = [Q(0), Q(0), -Q(1, 8)] + [Q(0)] * 36
    for point in points:
        assert [evaluate(poly, point) for poly in polynomials] == expected_target

    quadratic = [{e: c for e, c in poly.items() if sum(e) == 2} for poly in polynomials]
    cubic = [{e: c for e, c in poly.items() if sum(e) == 3} for poly in polynomials]
    final_dimension = 2 * dimension + 1
    t_index = final_dimension - 1
    h_components: list[Poly] = [{} for _ in range(final_dimension)]
    for index, (q_part, c_part) in enumerate(zip(quadratic, cubic)):
        exponent = [0] * final_dimension
        exponent[dimension + index], exponent[t_index] = 1, 2
        h_components[index][tuple(exponent)] = Q(1)
        for old, coefficient in q_part.items():
            exponent = list(old) + [0] * (dimension + 1)
            exponent[t_index] = 1
            h_components[index][tuple(exponent)] = coefficient
        for old, coefficient in c_part.items():
            h_components[dimension + index][old + (0,) * (dimension + 1)] = -coefficient

    stored_h = [
        {dense(term["monomial"], final_dimension): Q(term["coefficient"]) for term in component}
        for component in stored["H"]
    ]
    assert h_components == stored_h
    assert all(sum(exponents) == 3 for component in h_components for exponents in component)

    final_points: list[list[Q]] = []
    for point in points:
        c_value = [evaluate(poly, point) for poly in cubic]
        final_points.append(point + c_value + [Q(1)])
    assert final_points == [[Q(value) for value in point] for point in stored["collision_points"]]

    common_image = [Q(value) for value in stored["common_image"]]
    for point in final_points:
        image = [coordinate + evaluate(component, point) for coordinate, component in zip(point, h_components)]
        assert image == common_image
    assert len({tuple(point) for point in final_points}) == 3

    print("PASS (stdlib) Long BCW: regenerated the determinant-one starting map")
    print("PASS (stdlib) Long BCW: independently replayed all 18 balanced steps")
    print("PASS (stdlib) Long BCW: reconstructed the exact 79-variable sparse H")
    print("PASS (stdlib) Long BCW: verified three stored points and their common image")


if __name__ == "__main__":
    main()
