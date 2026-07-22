#!/usr/bin/env python3
"""Dependency-free replay of the constant-kernel 24D-to-22D quotient."""

from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "rank_compressed_bcw_24_counterexample.json"
ARTIFACT = ROOT / "artifacts" / "generated-results" / "constant_kernel_bcw_22_counterexample.json"
Poly = dict[tuple[int, ...], Q]
Matrix = list[list[Q]]


def clean(poly: Poly) -> Poly:
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def add(*polys: Poly) -> Poly:
    result: Poly = {}
    for poly in polys:
        for monomial, coefficient in poly.items():
            result[monomial] = result.get(monomial, Q(0)) + coefficient
    return clean(result)


def scale(poly: Poly, coefficient: Q) -> Poly:
    return clean({monomial: coefficient * value for monomial, value in poly.items()})


def multiply(left: Poly, right: Poly) -> Poly:
    result: Poly = {}
    for a, ca in left.items():
        for b, cb in right.items():
            monomial = tuple(x + y for x, y in zip(a, b))
            result[monomial] = result.get(monomial, Q(0)) + ca * cb
    return clean(result)


def power(poly: Poly, exponent: int) -> Poly:
    result = {(0,) * len(next(iter(poly), ())): Q(1)}
    for _ in range(exponent):
        result = multiply(result, poly)
    return result


def variable(dimension: int, index: int) -> Poly:
    monomial = [0] * dimension
    monomial[index] = 1
    return {tuple(monomial): Q(1)}


def decode_h(stored: dict[str, object]) -> list[Poly]:
    dimension = stored["dimension"]
    result = []
    for component in stored["H"]:
        poly: Poly = {}
        for term in component:
            monomial = [0] * dimension
            for index, exponent in term["monomial"]:
                monomial[index] = exponent
            key = tuple(monomial)
            poly[key] = poly.get(key, Q(0)) + Q(term["coefficient"])
        result.append(clean(poly))
    return result


def derivative(poly: Poly, index: int) -> Poly:
    result: Poly = {}
    for monomial, coefficient in poly.items():
        if monomial[index]:
            lowered = list(monomial)
            lowered[index] -= 1
            result[tuple(lowered)] = coefficient * monomial[index]
    return clean(result)


def substitute_linear(poly: Poly, linear_forms: list[Poly]) -> Poly:
    dimension = len(next(iter(linear_forms[0])))
    result: Poly = {}
    for monomial, coefficient in poly.items():
        term = {(0,) * dimension: coefficient}
        for index, exponent in enumerate(monomial):
            if exponent:
                term = multiply(term, power(linear_forms[index], exponent))
        result = add(result, term)
    return result


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return [
        [sum((left[i][k] * right[k][j] for k in range(len(right))), Q(0))
         for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def matvec_polynomial(matrix: Matrix, vector: list[Poly]) -> list[Poly]:
    return [add(*(scale(poly, coefficient) for coefficient, poly in zip(row, vector))) for row in matrix]


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def rref(matrix: Matrix) -> tuple[Matrix, list[int]]:
    work = [row[:] for row in matrix]
    if not work:
        return work, []
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(len(work[0])):
        selected = next((row for row in range(pivot_row, len(work)) if work[row][column]), None)
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        pivot = work[pivot_row][column]
        work[pivot_row] = [value / pivot for value in work[pivot_row]]
        for row in range(len(work)):
            if row != pivot_row and work[row][column]:
                factor = work[row][column]
                work[row] = [a - factor * b for a, b in zip(work[row], work[pivot_row])]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    return work, pivot_columns


def nullspace(matrix: Matrix) -> Matrix:
    reduced, pivots = rref(matrix)
    columns = len(matrix[0])
    free = [column for column in range(columns) if column not in pivots]
    basis: Matrix = []
    for free_column in free:
        vector = [Q(0)] * columns
        vector[free_column] = Q(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][free_column]
        basis.append(vector)
    return transpose(basis)


def inverse(matrix: Matrix) -> Matrix:
    size = len(matrix)
    augmented = [row[:] + [Q(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    reduced, pivots = rref(augmented)
    assert pivots[:size] == list(range(size))
    return [row[size:] for row in reduced]


def sparse_rows(stored: list[list[list[object]]], columns: int) -> Matrix:
    rows: Matrix = []
    for sparse in stored:
        row = [Q(0)] * columns
        for index, value in sparse:
            row[index] = Q(value)
        rows.append(row)
    return rows


def evaluate(poly: Poly, point: list[Q]) -> Q:
    return sum(
        (coefficient * __import__("math").prod(x**e for x, e in zip(point, monomial))
         for monomial, coefficient in poly.items()),
        Q(0),
    )


def main() -> None:
    source = json.loads(SOURCE.read_text())
    stored = json.loads(ARTIFACT.read_text())
    assert source["dimension"] == 24 and stored["dimension"] == 22
    h24 = decode_h(source)

    # Recompute the simultaneous constant kernel of all coefficient matrices
    # of JH, without trusting the quotient artifact's declared basis.
    jacobian = [[derivative(h24[i], j) for j in range(24)] for i in range(24)]
    monomials = sorted({m for row in jacobian for poly in row for m in poly})
    coefficient_matrix = [
        [jacobian[i][j].get(monomial, Q(0)) for j in range(24)]
        for i in range(24) for monomial in monomials
    ]
    K = nullspace(coefficient_matrix)
    assert len(K) == 24 and len(K[0]) == 2
    expected = [[Q(0), Q(0)] for _ in range(24)]
    expected[6][0], expected[10][0] = Q(-3), Q(1)
    expected[12][1], expected[13][1], expected[15][1] = Q(-2, 11), Q(-3, 11), Q(1)
    assert rref(transpose(K))[0] == rref(transpose(expected))[0]
    K = expected

    factorization = stored["quotient_factorization"]
    B = sparse_rows(factorization["B_rows"], 24)
    C = sparse_rows(factorization["C_rows"], 22)
    assert matmul(B, K) == [[Q(0)] * 2 for _ in range(22)]
    assert matmul(B, C) == [[Q(i == j) for j in range(22)] for i in range(22)]

    x22 = [variable(22, i) for i in range(22)]
    cb = matmul(C, B)
    cb_forms = matvec_polynomial(cb, [variable(24, i) for i in range(24)])
    assert [substitute_linear(poly, cb_forms) for poly in h24] == h24

    section_forms = matvec_polynomial(C, x22)
    quotient_h = matvec_polynomial(B, [substitute_linear(poly, section_forms) for poly in h24])
    assert all(sum(monomial) == 3 for poly in quotient_h for monomial in poly)
    assert quotient_h == decode_h(stored)

    points = [[Q(value) for value in point] for point in source["collision_points"]]
    quotient_points = [[sum((B[i][j] * point[j] for j in range(24)), Q(0)) for i in range(22)] for point in points]
    assert quotient_points == [[Q(value) for value in point] for point in stored["collision_points"]]
    assert len({tuple(point) for point in quotient_points}) == 3
    common = [Q(value) for value in stored["common_image"]]
    for point in quotient_points:
        assert [point[i] + evaluate(quotient_h[i], point) for i in range(22)] == common

    # In S=(C|K) coordinates, the nonlinear part has no fiber variables.
    # The conjugated full map therefore has Jacobian [[Df,0],[Dg,I_2]], so
    # det D(S^-1 F S)=det Df exactly.
    S = [C[row] + K[row] for row in range(24)]
    Sinv = inverse(S)
    assert Sinv[:22] == B
    qr = [variable(24, i) for i in range(24)]
    s_forms = matvec_polynomial(S, qr)
    conjugated_h = matvec_polynomial(Sinv, [substitute_linear(poly, s_forms) for poly in h24])
    assert all(not any(any(monomial[22:]) for monomial in poly) for poly in conjugated_h)
    assert [
        {monomial[:22]: coefficient for monomial, coefficient in poly.items()}
        for poly in conjugated_h[:22]
    ] == quotient_h

    print("PASS (stdlib) constant-kernel BCW: recomputed dim ker(JH_24)=2")
    print("PASS (stdlib) constant-kernel BCW: verified BK=0, BC=I, and H=HCB")
    print("PASS (stdlib) constant-kernel BCW: rebuilt the cubic-homogeneous 22D quotient")
    print("PASS (stdlib) constant-kernel BCW: verified the descended three-point collision")
    print("PASS (stdlib) constant-kernel BCW: verified block-triangular determinant factorization")


if __name__ == "__main__":
    main()
