#!/usr/bin/env python3
"""Dependency-free sparse replay of the rank-compressed 24-variable artifact."""

from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "shared_bcw_33_counterexample.json"
ARTIFACT = ROOT / "artifacts" / "generated-results" / "rank_compressed_bcw_24_counterexample.json"
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


def variable(dimension: int, index: int) -> Poly:
    exponents = [0] * dimension
    exponents[index] = 1
    return {tuple(exponents): Q(1)}


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


def pad(poly: Poly, amount: int) -> Poly:
    return {exponents + (0,) * amount: coefficient for exponents, coefficient in poly.items()}


def coefficient_rows(polynomials: list[Poly]) -> tuple[list[tuple[int, ...]], list[list[Q]]]:
    monomials = sorted({exponents for poly in polynomials for exponents in poly})
    return monomials, [[poly.get(exponents, Q(0)) for exponents in monomials] for poly in polynomials]


def rref(matrix: list[list[Q]]) -> tuple[list[list[Q]], list[int]]:
    work = [row[:] for row in matrix]
    if not work:
        return work, []
    rows, columns = len(work), len(work[0])
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(columns):
        selected = next((row for row in range(pivot_row, rows) if work[row][column]), None)
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        divisor = work[pivot_row][column]
        work[pivot_row] = [value / divisor for value in work[pivot_row]]
        for row in range(rows):
            if row != pivot_row and work[row][column]:
                factor = work[row][column]
                work[row] = [a - factor * b for a, b in zip(work[row], work[pivot_row])]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break
    return work, pivot_columns


def independent_rows(matrix: list[list[Q]]) -> list[int]:
    transpose = [list(column) for column in zip(*matrix)] if matrix and matrix[0] else []
    return rref(transpose)[1]


def solve_square(matrix: list[list[Q]], target: list[Q]) -> list[Q]:
    augmented = [row[:] + [value] for row, value in zip(matrix, target)]
    reduced, pivots = rref(augmented)
    size = len(matrix)
    assert pivots[:size] == list(range(size))
    return [reduced[index][-1] for index in range(size)]


def extract_source(source: dict[str, object]) -> tuple[list[Poly], list[Poly]]:
    n = source["degree_reduced_dimension"]
    old_dimension = source["dimension"]
    assert old_dimension == 2 * n + 1
    t_index = old_dimension - 1
    quadratic: list[Poly] = []
    cubic: list[Poly] = []
    for index in range(n):
        q_part: Poly = {}
        for term in source["H"][index]:
            exponents = dense(term["monomial"], old_dimension)
            coefficient = Q(term["coefficient"])
            if exponents[n + index] == 1 and exponents[t_index] == 2:
                assert coefficient == 1 and sum(exponents) == 3
                continue
            assert exponents[t_index] == 1 and not any(exponents[n:2 * n])
            old = exponents[:n]
            assert sum(old) == 2
            q_part[old] = q_part.get(old, Q(0)) + coefficient
        c_part: Poly = {}
        for term in source["H"][n + index]:
            exponents = dense(term["monomial"], old_dimension)
            assert not any(exponents[n:]) and sum(exponents[:n]) == 3
            old = exponents[:n]
            c_part[old] = c_part.get(old, Q(0)) - Q(term["coefficient"])
        quadratic.append(clean(q_part))
        cubic.append(clean(c_part))
    return quadratic, cubic


def main() -> None:
    source = json.loads(SOURCE.read_text())
    stored = json.loads(ARTIFACT.read_text())
    assert stored["format"] == "rank-compressed-bcw-sparse-cubic-homogeneous-map-v1"
    assert stored["source_artifact"] == str(SOURCE.relative_to(ROOT))
    quadratic, cubic = extract_source(source)
    n = len(cubic)

    monomials, rows = coefficient_rows(cubic)
    basis_components = independent_rows(rows)
    assert basis_components == [0, 1, 2, 3, 4, 6, 8]
    k = len(basis_components)
    assert k == stored["cubic_output_rank"] == 7

    # Express every cubic row in the independent component rows.  Pivot
    # monomial columns give a nonsingular k-by-k coordinate minor.
    basis_rows = [rows[index] for index in basis_components]
    pivot_monomials = rref(basis_rows)[1]
    minor = [[basis_rows[row][column] for row in range(k)] for column in pivot_monomials]
    B: list[list[Q]] = []
    for row in rows:
        target = [row[column] for column in pivot_monomials]
        coefficients = solve_square(minor, target)
        assert all(
            sum(coefficients[j] * basis_rows[j][column] for j in range(k)) == row[column]
            for column in range(len(monomials))
        )
        B.append(coefficients)
    assert B == [[Q(value) for value in row] for row in stored["rank_factorization"]["B"]]
    assert basis_components == stored["rank_factorization"]["basis_components"]

    final_dimension = n + k + 1
    t_index = final_dimension - 1
    t = variable(final_dimension, t_index)
    x = [variable(final_dimension, index) for index in range(n)]
    y = [variable(final_dimension, n + index) for index in range(k)]
    q = [pad(poly, k + 1) for poly in quadratic]
    all_c = [pad(poly, k + 1) for poly in cubic]
    c = [all_c[index] for index in basis_components]

    # Sparse replay of P_t o (E_t x id) o A_t = id+tN.
    a_second = [add(y[j], scale(multiply(t, c[j]), -1)) for j in range(k)]
    e_first = [
        add(x[i], multiply(t, q[i]), multiply(multiply(t, t), all_c[i]))
        for i in range(n)
    ]
    right_first = [
        add(
            e_first[i],
            multiply(t, add(*(scale(a_second[j], B[i][j]) for j in range(k)))),
        )
        for i in range(n)
    ]
    left_first = [
        add(
            x[i],
            multiply(
                t,
                add(q[i], *(scale(y[j], B[i][j]) for j in range(k))),
            ),
        )
        for i in range(n)
    ]
    assert right_first == left_first
    assert a_second == [add(y[j], scale(multiply(t, c[j]), -1)) for j in range(k)]

    h_components: list[Poly] = [{} for _ in range(final_dimension)]
    for i in range(n):
        h_components[i] = add(
            multiply(t, q[i]),
            multiply(
                multiply(t, t),
                add(*(scale(y[j], B[i][j]) for j in range(k))),
            ),
        )
    for j in range(k):
        h_components[n + j] = scale(c[j], -1)
    stored_h = [
        {dense(term["monomial"], final_dimension): Q(term["coefficient"]) for term in component}
        for component in stored["H"]
    ]
    assert h_components == stored_h
    assert all(sum(exponents) == 3 for component in h_components for exponents in component)

    old_points = [[Q(value) for value in point[:n]] for point in source["collision_points"]]
    old_target = [Q(value) for value in source["common_image"][:n]]
    final_points = [point + [evaluate(cubic[index], point) for index in basis_components] + [Q(1)] for point in old_points]
    assert final_points == [[Q(value) for value in point] for point in stored["collision_points"]]
    common_image = [Q(value) for value in stored["common_image"]]
    assert common_image == old_target + [Q(0)] * k + [Q(1)]
    for point in final_points:
        image = [coordinate + evaluate(component, point) for coordinate, component in zip(point, h_components)]
        assert image == common_image
    assert len({tuple(point) for point in final_points}) == 3

    print("PASS (stdlib) rank-compressed BCW: independently recovered exact QQ-rank 7")
    print("PASS (stdlib) rank-compressed BCW: replayed id+tN factorization sparsely")
    print("PASS (stdlib) rank-compressed BCW: reconstructed the 24-variable cubic map")
    print("PASS (stdlib) rank-compressed BCW: verified three stored points and their common image")


if __name__ == "__main__":
    main()
