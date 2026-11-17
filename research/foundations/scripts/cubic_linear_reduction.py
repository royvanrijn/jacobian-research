#!/usr/bin/env python3
"""Pair the explicit Yagzhev counterexample with a cubic-linear map.

Starting from the sparse 95-dimensional map f=x+H, write

    H(x) = -B(Dx)^{*3}

using the standard polarization identities for cubic monomials.  After adding
a smallest coordinate complement to B, set A=DB.  Then ker(A)=ker(B), and

    G(X) = X - (AX)^{*3}

is Gorni--Zampieri paired with f.  The resulting Druzkowski map has dimension
451.  The script also transports the three-point collision explicitly.
"""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "cubic_homogeneous_counterexample.json"
OUTPUT = ROOT / "artifacts" / "generated-results" / "cubic_linear_counterexample.json"


def qtext(value: sp.Expr) -> str:
    value = sp.cancel(value)
    assert value.is_Rational
    return str(value)


def canonical_form(form: list[int], coefficient: sp.Rational) -> tuple[tuple[int, ...], sp.Rational]:
    first = next(value for value in form if value)
    if first < 0:
        return tuple(-value for value in form), -coefficient
    return tuple(form), coefficient


def add_cube(
    cubes: dict[tuple[int, ...], dict[int, sp.Rational]],
    form: list[int],
    output_index: int,
    coefficient: sp.Rational,
) -> None:
    key, coefficient = canonical_form(form, coefficient)
    column = cubes.setdefault(key, {})
    column[output_index] = sp.cancel(column.get(output_index, 0) + coefficient)
    if column[output_index] == 0:
        del column[output_index]


def add_monomial_decomposition(cubes, dimension, output_index, coefficient, support):
    """Add B-coefficients whose cube sum equals -coefficient*monomial."""
    exponents = {index: exponent for index, exponent in support}
    indices = sorted(exponents)
    target = -sp.Rational(coefficient)
    if len(indices) == 1:
        form = [0] * dimension
        form[indices[0]] = 1
        add_cube(cubes, form, output_index, target)
        return

    if len(indices) == 2:
        single = next(index for index in indices if exponents[index] == 1)
        squared = next(index for index in indices if exponents[index] == 2)
        # a*b^2=((a+b)^3+(a-b)^3-2a^3)/6.
        for a_sign, b_sign, factor in (
            (1, 1, sp.Rational(1, 6)),
            (1, -1, sp.Rational(1, 6)),
            (1, 0, -sp.Rational(1, 3)),
        ):
            form = [0] * dimension
            form[single] = a_sign
            form[squared] = b_sign
            add_cube(cubes, form, output_index, target * factor)
        return

    assert len(indices) == 3 and all(exponents[index] == 1 for index in indices)
    a, b, c = indices
    # abc=((a+b+c)^3+(a-b-c)^3-(a+b-c)^3-(a-b+c)^3)/24.
    for signs, factor in (
        ((1, 1, 1), sp.Rational(1, 24)),
        ((1, -1, -1), sp.Rational(1, 24)),
        ((1, 1, -1), -sp.Rational(1, 24)),
        ((1, -1, 1), -sp.Rational(1, 24)),
    ):
        form = [0] * dimension
        for index, sign in zip((a, b, c), signs):
            form[index] = sign
        add_cube(cubes, form, output_index, target * factor)


def sparse_vector(entries: dict[int, sp.Expr]) -> list[list[object]]:
    return [[index, qtext(value)] for index, value in sorted(entries.items()) if value != 0]


def eval_rows(rows, point):
    return [
        sp.cancel(sum(sp.Rational(value) * point[index] for index, value in row))
        for row in rows
    ]


def cube_expansion(form: tuple[int, ...]) -> dict[tuple[int, ...], sp.Rational]:
    support = [(index, value) for index, value in enumerate(form) if value]
    result: dict[tuple[int, ...], sp.Rational] = {}
    for left, middle, right in product(support, repeat=3):
        exponents = [0] * len(form)
        coefficient = sp.Integer(1)
        for index, value in (left, middle, right):
            exponents[index] += 1
            coefficient *= value
        key = tuple(exponents)
        result[key] = sp.cancel(result.get(key, 0) + coefficient)
    return {key: value for key, value in result.items() if value}


def main() -> None:
    source = json.loads(SOURCE.read_text())
    n = source["dimension"]
    assert n == 95

    cubes: dict[tuple[int, ...], dict[int, sp.Rational]] = {}
    for output_index, component in enumerate(source["H"]):
        for term in component:
            add_monomial_decomposition(
                cubes, n, output_index, term["coefficient"], term["monomial"]
            )
    cubes = {form: column for form, column in cubes.items() if column}
    forms = sorted(cubes)
    m = len(forms)
    assert m == 415

    # D0 already has full column rank.  B0 has rank 59, so append only a
    # 36-column coordinate complement, rather than the redundant I_n used by
    # the first 510-dimensional construction.
    b_columns: list[dict[int, sp.Rational]] = [cubes[form] for form in forms]
    b0 = sp.zeros(n, m)
    for column_index, column in enumerate(b_columns):
        for row_index, value in column.items():
            b0[row_index, column_index] = value
    b0_rank = b0.rank()
    assert b0_rank == 59
    current = b0
    complement_indices = []
    current_rank = b0_rank
    for index in range(n):
        candidate = current.row_join(sp.eye(n)[:, index])
        candidate_rank = candidate.rank()
        if candidate_rank > current_rank:
            complement_indices.append(index)
            current = candidate
            current_rank = candidate_rank
    assert current_rank == n and len(complement_indices) == n - b0_rank == 36
    b_columns.extend([{index: sp.Integer(1)} for index in complement_indices])
    d_rows: list[dict[int, sp.Rational]] = [
        {index: sp.Integer(value) for index, value in enumerate(form) if value}
        for form in forms
    ]
    d_rows.extend([{} for _ in complement_indices])
    N = m + len(complement_indices)
    assert N == 451

    d_matrix = sp.zeros(N, n)
    for row_index, row in enumerate(d_rows):
        for column_index, value in row.items():
            d_matrix[row_index, column_index] = value
    assert d_matrix.rank() == n

    # Build row access for B, then A=DB without dense 451x451 multiplication.
    b_rows: list[dict[int, sp.Rational]] = [dict() for _ in range(n)]
    for column_index, column in enumerate(b_columns):
        for row_index, value in column.items():
            b_rows[row_index][column_index] = value
    a_rows: list[dict[int, sp.Rational]] = []
    for d_row in d_rows:
        a_row: dict[int, sp.Rational] = {}
        for intermediate, d_value in d_row.items():
            for column_index, b_value in b_rows[intermediate].items():
                a_row[column_index] = sp.cancel(
                    a_row.get(column_index, 0) + d_value * b_value
                )
        a_rows.append({index: value for index, value in a_row.items() if value})

    # Reconstruct H=-B(Dx)^3 directly as a sparse polynomial certificate.
    reconstructed: list[dict[tuple[int, ...], sp.Rational]] = [dict() for _ in range(n)]
    for column_index, form in enumerate(forms):
        expansion = cube_expansion(form)
        for output_index, b_value in b_columns[column_index].items():
            for monomial, cube_coefficient in expansion.items():
                reconstructed[output_index][monomial] = sp.cancel(
                    reconstructed[output_index].get(monomial, 0)
                    - b_value * cube_coefficient
                )
    reconstructed = [
        {monomial: coefficient for monomial, coefficient in component.items() if coefficient}
        for component in reconstructed
    ]
    expected: list[dict[tuple[int, ...], sp.Rational]] = []
    for component in source["H"]:
        encoded = {}
        for term in component:
            exponents = [0] * n
            for index, exponent in term["monomial"]:
                exponents[index] = exponent
            encoded[tuple(exponents)] = sp.Rational(term["coefficient"])
        expected.append(encoded)
    assert reconstructed == expected

    b_json = [sparse_vector(column) for column in b_columns]
    d_json = [sparse_vector(row) for row in d_rows]
    a_json = [sparse_vector(row) for row in a_rows]

    # Construct a sparse rational right inverse C from pivot columns of B.  If
    # f(p_i) is constant, then
    # G(Cp_i) differ by kernel vectors.  Translate each Cp_i by that difference
    # to obtain literal collisions for G.
    b_matrix = sp.zeros(n, N)
    for column_index, column in enumerate(b_columns):
        for row_index, value in column.items():
            b_matrix[row_index, column_index] = value
    _, pivot_columns = b_matrix.rref()
    assert len(pivot_columns) == n
    pivot_inverse = b_matrix[:, list(pivot_columns)].inv()
    c_matrix = sp.zeros(N, n)
    for local_row, global_row in enumerate(pivot_columns):
        c_matrix[global_row, :] = pivot_inverse[local_row, :]
    assert b_matrix * c_matrix == sp.eye(n)
    c_rows = [
        {column: c_matrix[row, column] for column in range(n) if c_matrix[row, column]}
        for row in range(N)
    ]

    source_points = [[sp.Rational(value) for value in point] for point in source["collision_points"]]
    base_points = [
        [sp.cancel(sum(value * point[index] for index, value in row.items())) for row in c_rows]
        for point in source_points
    ]

    def G(point):
        linear_forms = eval_rows(a_json, point)
        return [sp.cancel(value - linear_form**3) for value, linear_form in zip(point, linear_forms)]

    base_images = [G(point) for point in base_points]
    common_image = base_images[0]
    collision_points = []
    for point, image in zip(base_points, base_images):
        delta = [sp.cancel(want - got) for want, got in zip(common_image, image)]
        assert all(value == 0 for value in eval_rows(a_json, delta))
        lifted = [sp.cancel(value + shift) for value, shift in zip(point, delta)]
        assert G(lifted) == common_image
        collision_points.append(lifted)
    assert len({tuple(point) for point in collision_points}) == 3

    artifact = {
        "format": "sparse-cubic-linear-map-v1",
        "construction": "Gorni--Zampieri pairing / Druzkowski reduction",
        "source_cubic_homogeneous_artifact": str(SOURCE.relative_to(ROOT)),
        "paired_dimension": n,
        "dimension": N,
        "map": "G(X)=X-(A X)^{*3}",
        "jacobian_determinant": "1",
        "rank_A": n,
        "pairing": {
            "B_shape": [n, N],
            "B_columns": b_json,
            "D_shape": [N, n],
            "D_rows": d_json,
            "C_shape": [N, n],
            "C_rows": [sparse_vector(row) for row in c_rows],
            "identities": ["BC=I_95", "A=DB", "AC=D", "ker(A)=ker(B)"],
        },
        "A_rows": a_json,
        "collision_points": [[qtext(value) for value in point] for point in collision_points],
        "common_image": [qtext(value) for value in common_image],
        "statistics": {
            "source_cubic_terms": sum(len(component) for component in source["H"]),
            "unique_cube_forms": m,
            "rank_B0": b0_rank,
            "complement_columns": len(complement_indices),
            "complement_coordinate_indices": complement_indices,
            "dimension": N,
            "nonzero_rows_A": sum(bool(row) for row in a_rows),
            "nonzero_entries_A": sum(len(row) for row in a_rows),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(f"PASS: decomposed 148 cubic terms into {m} unique cubes of linear forms")
    print(f"PASS: rank(B0)={b0_rank}; appended {len(complement_indices)} complement columns")
    print(f"PASS: constructed cubic-linear G in dimension {N} with rank(A)={n}")
    print("PASS: exact GZ pairing identities certify det DG=1")
    print("PASS: transported three distinct rational points to one common image")
    print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
