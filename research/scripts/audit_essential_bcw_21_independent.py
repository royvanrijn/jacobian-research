#!/usr/bin/env python3
"""Dependency-free replay of the frozen 3 -> 17 -> 24 -> 21 BCW route."""

from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path

from audit_constant_kernel_bcw_22_independent import (
    inverse,
    matmul,
    matvec_polynomial,
    nullspace,
    rref,
    sparse_rows,
    substitute_linear,
    transpose,
)
from audit_rank_compressed_bcw_24_independent import (
    coefficient_rows,
    independent_rows,
    solve_square,
)
from audit_shared_bcw_33_independent import (
    Poly,
    add,
    clean,
    constant,
    dense,
    derivative,
    determinant3,
    evaluate,
    monomial,
    multiply,
    pad,
    power,
    scale,
    variable,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "artifacts" / "generated-results" / "essential_bcw_21_counterexample.json"


def decode_h(stored: dict[str, object]) -> list[Poly]:
    dimension = stored["dimension"]
    return [
        {
            dense(term["monomial"], dimension): Q(term["coefficient"])
            for term in component
        }
        for component in stored["H"]
    ]


def replay_trace(stored: dict[str, object]) -> tuple[list[Poly], list[list[Q]]]:
    dimension = 3
    x, y, z = (variable(dimension, index) for index in range(3))
    v = add(constant(dimension, 1), scale(multiply(x, y), 2))
    long_map = [
        add(
            multiply(power(v, 3), z),
            scale(
                multiply(
                    multiply(power(y, 2), v),
                    add(constant(dimension, 2), scale(multiply(x, y), 3)),
                ),
                4,
            ),
        ),
        add(
            y,
            scale(multiply(multiply(x, power(v, 2)), z), 3),
            scale(
                multiply(
                    multiply(x, power(y, 2)),
                    add(constant(dimension, 2), scale(multiply(x, y), 3)),
                ),
                12,
            ),
        ),
        add(scale(x, -1), scale(multiply(power(x, 2), y), 3), multiply(power(x, 3), z)),
    ]
    jacobian = [[derivative(long_map[i], j) for j in range(3)] for i in range(3)]
    assert determinant3(jacobian) == {(0, 0, 0): Q(1)}
    polynomials = [scale(long_map[2], -1), long_map[1], long_map[0]]
    points = [
        [Q(0), Q(0), -Q(1, 8)],
        [Q(1), -Q(3, 4), Q(13, 4)],
        [-Q(1), Q(3, 4), Q(13, 4)],
    ]
    registry: dict[tuple[int, ...], int] = {}
    for step in stored["degree_lowering_steps"]:
        component = step["component"]
        first = dense(step["first_factor"], dimension)
        second = dense(step["second_factor"], dimension)
        removed = dense(step["removed_monomial"], dimension)
        coefficient = Q(step["removed_coefficient"])
        assert tuple(a + b for a, b in zip(first, second)) == removed
        assert polynomials[component].get(removed) == coefficient
        missing = []
        for factor in (first, second):
            if factor not in registry and factor not in missing:
                missing.append(factor)
        records = step["new_factors"]
        assert len(records) == len(missing)
        old_dimension = dimension
        additions = [[] for _ in points]
        for offset, (factor, record) in enumerate(zip(missing, records)):
            assert dense(record["factor"], old_dimension) == factor
            assert record["new_variable"] == dimension
            assert record["output"] == len(polynomials) + offset
            for extra, point in zip(additions, points):
                extra.append(-evaluate(monomial(old_dimension, factor), point))
            dimension += 1
        added = len(missing)
        polynomials = [pad(poly, added) for poly in polynomials]
        registry = {factor + (0,) * added: output for factor, output in registry.items()}
        for factor, record in zip((factor + (0,) * added for factor in missing), records):
            registry[factor] = len(polynomials)
            polynomials.append(add(variable(dimension, record["new_variable"]), monomial(dimension, factor)))
        for point, extra in zip(points, additions):
            point.extend(extra)
        first += (0,) * added
        second += (0,) * added
        outputs = [registry[first], registry[second]]
        assert outputs == step["factor_outputs"] and component not in outputs
        polynomials[component] = add(
            polynomials[component],
            scale(multiply(polynomials[outputs[0]], polynomials[outputs[1]]), -coefficient),
        )
        registry = {factor: output for factor, output in registry.items() if output != component}

    assert dimension == stored["degree_reduced_dimension"] == 17
    assert max(sum(exponents) for poly in polynomials for exponents in poly) == 3
    for index, poly in enumerate(polynomials):
        linear = {e: c for e, c in poly.items() if sum(e) == 1}
        assert linear == {tuple(1 if j == index else 0 for j in range(dimension)): Q(1)}
    target = [Q(0), Q(0), -Q(1, 8)] + [Q(0)] * 14
    assert all([evaluate(poly, point) for poly in polynomials] == target for point in points)
    return polynomials, points


def main() -> None:
    stored = json.loads(ARTIFACT.read_text())
    assert stored["format"] == "essential-constant-kernel-bcw-sparse-cubic-homogeneous-map-v1"
    polynomials, points = replay_trace(stored)
    n = len(polynomials)
    quadratic = [{e: c for e, c in poly.items() if sum(e) == 2} for poly in polynomials]
    cubic = [{e: c for e, c in poly.items() if sum(e) == 3} for poly in polynomials]

    monomials, rows = coefficient_rows(cubic)
    basis = independent_rows(rows)
    assert basis == stored["rank_factorization"]["basis_components"]
    k = len(basis)
    assert k == 6
    basis_rows = [rows[index] for index in basis]
    pivot_monomials = rref(basis_rows)[1]
    minor = [[basis_rows[row][column] for row in range(k)] for column in pivot_monomials]
    Bc = []
    for row in rows:
        coefficients = solve_square(minor, [row[column] for column in pivot_monomials])
        assert all(
            sum(coefficients[j] * basis_rows[j][column] for j in range(k)) == row[column]
            for column in range(len(monomials))
        )
        Bc.append(coefficients)
    assert Bc == [[Q(value) for value in row] for row in stored["rank_factorization"]["B"]]

    ambient_dimension = n + k + 1
    assert ambient_dimension == stored["ambient_homogeneous_dimension"] == 24
    t_index = ambient_dimension - 1
    ambient_h: list[Poly] = [{} for _ in range(ambient_dimension)]
    for i in range(n):
        for old, coefficient in quadratic[i].items():
            exponents = list(old) + [0] * (k + 1)
            exponents[t_index] = 1
            ambient_h[i][tuple(exponents)] = coefficient
        for j in range(k):
            if Bc[i][j]:
                exponents = [0] * ambient_dimension
                exponents[n + j], exponents[t_index] = 1, 2
                ambient_h[i][tuple(exponents)] = Bc[i][j]
    for j, component in enumerate(basis):
        ambient_h[n + j] = {e + (0,) * (k + 1): -c for e, c in cubic[component].items()}

    jacobian = [[derivative(ambient_h[i], j) for j in range(ambient_dimension)] for i in range(ambient_dimension)]
    derivative_monomials = sorted({e for row in jacobian for poly in row for e in poly})
    coefficient_matrix = [
        [jacobian[i][j].get(e, Q(0)) for j in range(ambient_dimension)]
        for i in range(ambient_dimension) for e in derivative_monomials
    ]
    computed_kernel = nullspace(coefficient_matrix)
    assert len(computed_kernel[0]) == 3
    factor = stored["quotient_factorization"]
    K = transpose(sparse_rows(factor["kernel_columns"], ambient_dimension))
    assert rref(transpose(computed_kernel))[0] == rref(transpose(K))[0]
    B = sparse_rows(factor["B_rows"], ambient_dimension)
    C = sparse_rows(factor["C_rows"], 21)
    assert matmul(B, K) == [[Q(0)] * 3 for _ in range(21)]
    assert matmul(B, C) == [[Q(i == j) for j in range(21)] for i in range(21)]

    ambient_variables = [variable(ambient_dimension, i) for i in range(ambient_dimension)]
    cb_forms = matvec_polynomial(matmul(C, B), ambient_variables)
    assert [substitute_linear(poly, cb_forms) for poly in ambient_h] == ambient_h
    q_variables = [variable(21, i) for i in range(21)]
    section_forms = matvec_polynomial(C, q_variables)
    quotient_h = matvec_polynomial(
        B, [substitute_linear(poly, section_forms) for poly in ambient_h]
    )
    assert quotient_h == decode_h(stored)
    assert all(sum(e) == 3 for poly in quotient_h for e in poly)

    ambient_points = []
    for point in points:
        ambient_points.append(point + [evaluate(cubic[index], point) for index in basis] + [Q(1)])
    quotient_points = [
        [sum(B[i][j] * point[j] for j in range(ambient_dimension)) for i in range(21)]
        for point in ambient_points
    ]
    assert quotient_points == [[Q(value) for value in point] for point in stored["collision_points"]]
    assert len({tuple(point) for point in quotient_points}) == 3
    common = [Q(value) for value in stored["common_image"]]
    for point in quotient_points:
        assert [point[i] + evaluate(quotient_h[i], point) for i in range(21)] == common

    S = [C[row] + K[row] for row in range(ambient_dimension)]
    Sinv = inverse(S)
    assert Sinv[:21] == B
    conjugated = matvec_polynomial(
        Sinv,
        [substitute_linear(poly, matvec_polynomial(S, ambient_variables)) for poly in ambient_h],
    )
    assert all(not any(any(e[21:]) for e in poly) for poly in conjugated)

    print("PASS (stdlib) essential BCW 21: replayed 17 exact stable cancellations")
    print("PASS (stdlib) essential BCW 21: independently recovered cubic-output rank 6")
    print("PASS (stdlib) essential BCW 21: reconstructed the 24D homogeneous map")
    print("PASS (stdlib) essential BCW 21: recomputed the three-dimensional constant kernel")
    print("PASS (stdlib) essential BCW 21: rebuilt the 21D quotient and collision")
    print("PASS (stdlib) essential BCW 21: verified triangular determinant factorization")


if __name__ == "__main__":
    main()
