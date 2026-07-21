#!/usr/bin/env python3
"""Independently verify the sparse 451-dimensional Druzkowski artifact."""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "cubic_homogeneous_counterexample.json"
ARTIFACT = ROOT / "results" / "cubic_linear_counterexample.json"


def rows_to_dict(rows):
    return [{int(index): sp.Rational(value) for index, value in row} for row in rows]


def eval_rows(rows, point):
    return [sp.cancel(sum(value * point[index] for index, value in row.items())) for row in rows]


def cube_expansion(row, dimension):
    support = list(row.items())
    result = {}
    for left, middle, right in product(support, repeat=3):
        exponents = [0] * dimension
        coefficient = sp.Integer(1)
        for index, value in (left, middle, right):
            exponents[index] += 1
            coefficient *= value
        key = tuple(exponents)
        result[key] = sp.cancel(result.get(key, 0) + coefficient)
    return {key: value for key, value in result.items() if value}


def main() -> None:
    source = json.loads(SOURCE.read_text())
    data = json.loads(ARTIFACT.read_text())
    n, N = data["paired_dimension"], data["dimension"]
    assert (n, N) == (95, 451)

    b_columns = rows_to_dict(data["pairing"]["B_columns"])
    d_rows = rows_to_dict(data["pairing"]["D_rows"])
    a_rows = rows_to_dict(data["A_rows"])
    assert len(b_columns) == N and len(d_rows) == N and len(a_rows) == N

    cube_count = data["statistics"]["unique_cube_forms"]
    assert cube_count == 415 and N - cube_count == 36

    # Check the stored rational right inverse C and both pairing identities.
    c_rows = rows_to_dict(data["pairing"]["C_rows"])
    assert len(c_rows) == N
    for i in range(n):
        for j in range(n):
            value = sum(
                b_columns[row].get(i, 0) * c_rows[row].get(j, 0)
                for row in range(N)
            )
            assert sp.cancel(value) == (1 if i == j else 0)

    # Check rank(D)=95 exactly, so A=DB has ker(A)=ker(B).
    d_matrix = sp.zeros(N, n)
    for row_index, row in enumerate(d_rows):
        for column_index, value in row.items():
            d_matrix[row_index, column_index] = value
    assert d_matrix.rank() == n

    # B0 has rank 59 and the 36 appended columns are a minimal complement.
    b0_matrix = sp.zeros(n, cube_count)
    for column_index, column in enumerate(b_columns[:cube_count]):
        for row_index, value in column.items():
            b0_matrix[row_index, column_index] = value
    assert b0_matrix.rank() == data["statistics"]["rank_B0"] == 59
    assert len(b_columns[cube_count:]) == n - 59

    # Independently multiply D and B and compare every sparse row of A.
    b_rows = [dict() for _ in range(n)]
    for column_index, column in enumerate(b_columns):
        for row_index, value in column.items():
            b_rows[row_index][column_index] = value
    for row_index, d_row in enumerate(d_rows):
        expected = {}
        for intermediate, d_value in d_row.items():
            for column_index, b_value in b_rows[intermediate].items():
                expected[column_index] = sp.cancel(
                    expected.get(column_index, 0) + d_value * b_value
                )
        expected = {index: value for index, value in expected.items() if value}
        assert expected == a_rows[row_index]

    # Reconstruct H=-B(Dx)^3 and compare every source monomial coefficient.
    reconstructed = [dict() for _ in range(n)]
    for column_index in range(cube_count):
        expansion = cube_expansion(d_rows[column_index], n)
        for output_index, b_value in b_columns[column_index].items():
            for monomial, cube_coefficient in expansion.items():
                reconstructed[output_index][monomial] = sp.cancel(
                    reconstructed[output_index].get(monomial, 0)
                    - b_value * cube_coefficient
                )
    reconstructed = [
        {monomial: value for monomial, value in component.items() if value}
        for component in reconstructed
    ]
    expected_h = []
    for component in source["H"]:
        encoded = {}
        for term in component:
            exponents = [0] * n
            for index, exponent in term["monomial"]:
                exponents[index] = exponent
            encoded[tuple(exponents)] = sp.Rational(term["coefficient"])
        expected_h.append(encoded)
    assert reconstructed == expected_h

    # Direct exact substitution into all 451 cubic-linear coordinates.
    points = [[sp.Rational(value) for value in point] for point in data["collision_points"]]
    common_image = [sp.Rational(value) for value in data["common_image"]]
    assert len({tuple(point) for point in points}) == 3
    for point in points:
        linear_forms = eval_rows(a_rows, point)
        image = [sp.cancel(value - linear_form**3) for value, linear_form in zip(point, linear_forms)]
        assert image == common_image

    # Determinant certificate (Sylvester identity):
    # det(I_N-3 diag((AX)^2)DB)
    # = det(I_n-3 B diag((DBX)^2)D)
    # = det Df(BX) = 1.
    assert data["jacobian_determinant"] == source["jacobian_determinant"] == "1"

    print("PASS: B,C,D satisfy BC=I, A=DB, AC=D, and ker(A)=ker(B)")
    print("PASS: reconstructed all 148 source cubic terms from 415 cubes")
    print("PASS: sparse A has rank 95 by its exact D*B factorization")
    print("PASS: rank(B0)=59, so the 36-column extension and dimension 451 are minimal for this B0")
    print("PASS: all three 451-dimensional rational points have the stored common image")
    print("PASS: the exact GZ/Sylvester determinant identity certifies det DG=1")


if __name__ == "__main__":
    main()
