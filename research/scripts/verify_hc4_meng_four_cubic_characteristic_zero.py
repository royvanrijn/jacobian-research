#!/usr/bin/env python3
"""Promote the four-cubic Meng descent obstruction to characteristic zero.

The three-cubic QQ certificate supplies the 234 exact quartic principal
parts and the fixed Hessian tables.  This checker exhausts all 4,845
quadruples of cubic monomials over every quartic.  It computes the combined
determinant-degree-seven and degree-one rank over QQ, delegates coordinate
boundaries to the three-cubic theorem, and proves that every genuine null
line or higher-dimensional nullspace has unit full-determinant equations
over QQ.
"""

from __future__ import annotations

import contextlib
import io
from itertools import combinations
from pathlib import Path
import runpy

import sympy as sp


THREE_CUBIC = Path(__file__).with_name(
    "verify_hc4_meng_three_cubic_characteristic_zero.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    prior = runpy.run_path(str(THREE_CUBIC))

quartics = prior["quartics"]
quartic_exponents = prior["quartic_exponents"]
cubic_exponents = prior["cubic_exponents"]
base_hessian = prior["base_hessian"]
verification_points = prior["verification_points"]
verification_quartic_hessians = prior[
    "verification_quartic_hessians"
]
verification_cubic_hessians = prior["verification_cubic_hessians"]
degree_cubic_hessians = prior["degree_cubic_hessians"]
degree_one_signatures = prior["degree_one_signatures"]


def rank_and_null_basis(
    columns: tuple[tuple[sp.Rational, ...], ...],
) -> tuple[int, tuple[tuple[sp.Rational, ...], ...]]:
    """Exact RREF of a ten-row, four-column homogeneous system."""

    column_count = len(columns)
    rows = [list(row) for row in zip(*columns, strict=True)]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, len(rows))
                if rows[row][column] != 0
            ),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [
            sp.cancel(entry / scale) for entry in rows[pivot_row]
        ]
        for row in range(len(rows)):
            if row == pivot_row or rows[row][column] == 0:
                continue
            scale = rows[row][column]
            rows[row] = [
                sp.cancel(
                    rows[row][index]
                    - scale * rows[pivot_row][index]
                )
                for index in range(column_count)
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == column_count:
            break

    free_columns = [
        column for column in range(column_count) if column not in pivots
    ]
    basis = []
    for free in free_columns:
        vector = [sp.Rational(0)] * column_count
        vector[free] = sp.Rational(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = sp.factor(-rows[row][free])
        basis.append(tuple(vector))
    return len(pivots), tuple(basis)


def direction_hessians(
    cubic_quadruple: tuple[int, ...],
    directions: tuple[tuple[sp.Rational, ...], ...],
) -> tuple[tuple[sp.Matrix, ...], ...]:
    result = []
    for direction in directions:
        point_hessians = []
        for point_index in range(len(verification_points)):
            hessian = sp.zeros(4)
            for coefficient, cubic_index in zip(
                direction, cubic_quadruple, strict=True
            ):
                hessian += (
                    coefficient
                    * verification_cubic_hessians[point_index][cubic_index]
                )
            point_hessians.append(hessian)
        result.append(tuple(point_hessians))
    return tuple(result)


def unit_ideal_prefix(
    base_matrices: tuple[sp.Matrix, ...],
    directed_hessians: tuple[tuple[sp.Matrix, ...], ...],
    parameters: tuple[sp.Symbol, ...],
) -> int | None:
    equations = []
    for point_index in range(len(verification_points)):
        candidate = base_matrices[point_index].copy()
        for parameter, hessians in zip(
            parameters, directed_hessians, strict=True
        ):
            candidate += parameter * hessians[point_index]
        equation = sp.Poly(
            sp.expand(candidate.det(method="berkowitz") - 64),
            *parameters,
            domain=sp.QQ,
        )
        if not equation.is_zero:
            equations.append(equation.as_expr())
        if point_index < 2:
            continue
        basis = sp.groebner(
            equations,
            *parameters,
            domain=sp.QQ,
            order="grevlex",
        )
        if (
            len(basis.polys) == 1
            and sp.expand(basis.polys[0].as_expr()) == 1
        ):
            return point_index + 1
    return None


parameters = sp.symbols("lambda mu nu xi")
rank_counts = {rank: 0 for rank in range(5)}
boundary_counts = {1: 0, 2: 0, 3: 0}
genuine_counts = {1: 0, 2: 0, 3: 0}
maximum_points = {0: 0, 1: 0, 2: 0, 3: 0}

for support, coefficients in quartics:
    quartic_hessians = []
    base_matrices = []
    for point_index in range(len(verification_points)):
        quartic_hessian = sp.zeros(4)
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            quartic_hessian += (
                coefficient
                * verification_quartic_hessians[point_index][monomial_index]
            )
        quartic_hessians.append(quartic_hessian)
        base_matrices.append(base_hessian + quartic_hessian)
    base_matrices_tuple = tuple(base_matrices)

    signatures = []
    for cubic_index in range(len(cubic_exponents)):
        degree_seven = tuple(
            sp.factor(
                sp.trace(
                    quartic_hessians[point_index + 4].adjugate()
                    * degree_cubic_hessians[point_index][cubic_index]
                )
            )
            for point_index in range(5)
        )
        signatures.append(
            degree_seven + degree_one_signatures[cubic_index]
        )

    for cubic_quadruple in combinations(range(len(cubic_exponents)), 4):
        rank, null_basis = rank_and_null_basis(
            tuple(signatures[index] for index in cubic_quadruple)
        )
        rank_counts[rank] += 1
        if rank == 4:
            continue

        nullity = 4 - rank
        if rank in (1, 2, 3) and any(
            all(vector[index] == 0 for vector in null_basis)
            for index in range(4)
        ):
            boundary_counts[rank] += 1
            continue
        if rank in (1, 2, 3):
            genuine_counts[rank] += 1

        directed = direction_hessians(cubic_quadruple, null_basis)
        if nullity == 1:
            common_gcd: sp.Poly | None = None
            for point_index in range(len(verification_points)):
                equation = sp.Poly(
                    sp.expand(
                        (
                            base_matrices[point_index]
                            + parameters[0] * directed[0][point_index]
                        ).det(method="berkowitz")
                        - 64
                    ),
                    parameters[0],
                    domain=sp.QQ,
                )
                if equation.is_zero:
                    continue
                common_gcd = (
                    equation
                    if common_gcd is None
                    else sp.gcd(common_gcd, equation)
                )
                if common_gcd.degree() == 0:
                    maximum_points[rank] = max(
                        maximum_points[rank], point_index + 1
                    )
                    break
            assert common_gcd is not None and common_gcd.degree() == 0
            continue

        prefix = unit_ideal_prefix(
            base_matrices_tuple,
            directed,
            parameters[:nullity],
        )
        assert prefix is not None
        maximum_points[rank] = max(maximum_points[rank], prefix)


assert rank_counts == {
    0: 5_430,
    1: 79_396,
    2: 353_740,
    3: 504_818,
    4: 190_346,
}
assert boundary_counts == {
    1: 71_440,
    2: 347_658,
    3: 504_352,
}
assert genuine_counts == {
    1: 7_956,
    2: 6_082,
    3: 466,
}

print(
    "PASS: the exact QQ odd-layer ranks are "
    "r0=5430, r1=79396, r2=353740, r3=504818, r4=190346"
)
print(
    "PASS: all 466 genuine rank-three lines have unit determinant gcd in QQ"
)
print(
    "PASS: all 6082 genuine rank-two planes have unit determinant ideal "
    "in QQ[lambda,mu]"
)
print(
    "PASS: all 7956 genuine rank-one three-spaces have unit determinant "
    "ideal in QQ[lambda,mu,nu]"
)
print(
    "PASS: all 5430 rank-zero four-spaces have unit determinant ideal "
    "in QQ[lambda,mu,nu,xi]"
)
print(
    "DETAIL: maximum evaluation-prefix lengths for ranks 3,2,1,0 are "
    f"{maximum_points[3]},{maximum_points[2]},"
    f"{maximum_points[1]},{maximum_points[0]}"
)
print(
    "SCOPE: cubic support <=4 is excluded in characteristic zero; full "
    "cubic-kernel and dense-quartic corrections remain open"
)
