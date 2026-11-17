#!/usr/bin/env python3
"""Promote the three-cubic Meng descent obstruction to characteristic zero.

The parent sparse-quartic certificate leaves 234 quartic principal parts.
This checker reconstructs every one over QQ and exhausts all 1,140 triples
of cubic monomials.  Determinant degrees seven and one give ten homogeneous
linear equations in the three cubic coefficients.  Their ranks over QQ are
computed without modular inference.

Rank three is inconsistent.  Rank-two null lines are either boundaries of
the already-certified two-cubic calculation or have unit exact univariate
determinant gcd.  Rank-one null planes and rank-zero three-spaces are tested
by exact Groebner bases in QQ[lambda,mu] and QQ[lambda,mu,nu].  Evaluation
equations are added only until the unit ideal is reached.
"""

from __future__ import annotations

import contextlib
import io
from itertools import combinations
from pathlib import Path
import runpy

import sympy as sp


PARENT = Path(__file__).with_name(
    "verify_hc4_meng_sparse_quartic_obstruction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    parent = runpy.run_path(str(PARENT))

VARIABLE_COUNT: int = parent["VARIABLE_COUNT"]
quartic_exponents = parent["quartic_exponents"]
cubic_exponents = parent["cubic_exponents"]
base_hessian: sp.Matrix = parent["base_hessian"]
collision_point = parent["collision_point"]
scaled_collision_target = parent["scaled_collision_target"]
collision_monomials_exact = parent["collision_monomials_exact"]
principal_unique_survivors = parent["principal_unique_survivors"]
principal_family_survivors = parent["principal_family_survivors"]
monomial_hessian_exact = parent["monomial_hessian_exact"]
tau = parent["tau"]


def reconstruct_quartics() -> tuple[
    tuple[tuple[int, ...], tuple[sp.Rational, ...]], ...
]:
    quartics = []
    target = sp.Matrix(scaled_collision_target)
    for support, _ in principal_unique_survivors:
        exponent_matrix = sp.Matrix.hstack(
            *(sp.Matrix(quartic_exponents[index]) for index in support)
        )
        solution = next(iter(sp.linsolve((exponent_matrix, target))))
        assert not set().union(
            *(entry.free_symbols for entry in solution)
        )
        coefficients = tuple(
            sp.factor(
                solution[position] / collision_monomials_exact[index]
            )
            for position, index in enumerate(support)
        )
        quartics.append((support, coefficients))

    for support, coefficient_family, principal_gcd in (
        principal_family_survivors
    ):
        roots = sp.solve(principal_gcd.as_expr(), tau)
        assert len(roots) == 1
        coefficients = tuple(
            sp.factor(coefficient.subs(tau, roots[0]))
            for coefficient in coefficient_family
        )
        quartics.append((support, coefficients))
    assert len(quartics) == 234
    return tuple(quartics)


def exact_rank_data(
    columns: tuple[tuple[sp.Rational, ...], ...],
) -> tuple[
    int,
    tuple[sp.Rational, ...] | None,
    tuple[
        tuple[sp.Rational, ...], tuple[sp.Rational, ...]
    ] | None,
]:
    """Return rank, a corank-one vector, or a rank-one null-plane basis."""

    rows = tuple(zip(*columns, strict=True))
    first = next((tuple(row) for row in rows if any(row)), None)
    if first is None:
        return 0, None, None

    second = next(
        (
            tuple(row)
            for row in rows
            if any(
                first[left] * row[right] - first[right] * row[left]
                for left in range(3)
                for right in range(left)
            )
        ),
        None,
    )
    if second is None:
        pivot = next(index for index, entry in enumerate(first) if entry)
        basis = []
        for free in range(3):
            if free == pivot:
                continue
            vector = [sp.Rational(0)] * 3
            vector[free] = sp.Rational(1)
            vector[pivot] = -first[free] / first[pivot]
            basis.append(tuple(vector))
        assert len(basis) == 2
        return 1, first, (basis[0], basis[1])

    cross = (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )
    if all(
        sum(row[index] * cross[index] for index in range(3)) == 0
        for row in rows
    ):
        return 2, tuple(sp.factor(entry) for entry in cross), None
    return 3, None, None


def is_unit_groebner(basis: sp.GroebnerBasis) -> bool:
    return (
        len(basis.polys) == 1
        and sp.expand(basis.polys[0].as_expr()) == 1
    )


quartics = reconstruct_quartics()
verification_points = tuple(
    dict.fromkeys(
        tuple(
            sp.Rational(coordinate)
            for coordinate in point
        )
        for point in (
            tuple(parent["sample_points"])
            + tuple(parent["mixed_sample_points"])
        )
    )
)
degree_points = verification_points[4:9]
verification_quartic_hessians = tuple(
    tuple(
        monomial_hessian_exact(exponents, point)
        for exponents in quartic_exponents
    )
    for point in verification_points
)
verification_cubic_hessians = tuple(
    tuple(
        monomial_hessian_exact(exponents, point)
        for exponents in cubic_exponents
    )
    for point in verification_points
)
degree_quartic_hessians = verification_quartic_hessians[4:9]
degree_cubic_hessians = verification_cubic_hessians[4:9]
base_adjugate = base_hessian.adjugate()
degree_one_signatures = tuple(
    tuple(
        sp.trace(
            base_adjugate
            * degree_cubic_hessians[point_index][cubic_index]
        )
        for point_index in range(len(degree_points))
    )
    for cubic_index in range(len(cubic_exponents))
)

rank_counts = {0: 0, 1: 0, 2: 0, 3: 0}
rank_two_boundary_lines = 0
rank_two_genuine_lines = 0
rank_one_boundary_planes = 0
rank_one_genuine_planes = 0
maximum_rank_two_points = 0
maximum_rank_one_points = 0
maximum_rank_zero_points = 0
lambda_parameter, mu_parameter, nu_parameter = sp.symbols(
    "lambda mu nu"
)

for support, coefficients in quartics:
    quartic_hessians = []
    base_matrices = []
    for point_index in range(len(verification_points)):
        quartic_hessian = sp.zeros(VARIABLE_COUNT)
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            quartic_hessian += (
                coefficient
                * verification_quartic_hessians[point_index][monomial_index]
            )
        quartic_hessians.append(quartic_hessian)
        base_matrices.append(base_hessian + quartic_hessian)

    odd_layer_signatures = []
    for cubic_index in range(len(cubic_exponents)):
        degree_seven = tuple(
            sp.factor(
                sp.trace(
                    quartic_hessians[point_index + 4].adjugate()
                    * degree_cubic_hessians[point_index][cubic_index]
                )
            )
            for point_index in range(len(degree_points))
        )
        odd_layer_signatures.append(
            degree_seven + degree_one_signatures[cubic_index]
        )

    for cubic_triple in combinations(range(len(cubic_exponents)), 3):
        rank, line_direction, plane_basis = exact_rank_data(
            tuple(
                odd_layer_signatures[cubic_index]
                for cubic_index in cubic_triple
            )
        )
        rank_counts[rank] += 1
        if rank == 3:
            continue

        if rank == 2:
            assert line_direction is not None
            if 0 in line_direction:
                rank_two_boundary_lines += 1
                continue
            rank_two_genuine_lines += 1
            common_gcd: sp.Poly | None = None
            for point_index in range(len(verification_points)):
                directed_hessian = sp.zeros(VARIABLE_COUNT)
                for coefficient, cubic_index in zip(
                    line_direction, cubic_triple, strict=True
                ):
                    directed_hessian += (
                        coefficient
                        * verification_cubic_hessians[point_index][
                            cubic_index
                        ]
                    )
                equation = sp.Poly(
                    sp.expand(
                        (
                            base_matrices[point_index]
                            + lambda_parameter * directed_hessian
                        ).det(method="berkowitz")
                        - 64
                    ),
                    lambda_parameter,
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
                    maximum_rank_two_points = max(
                        maximum_rank_two_points, point_index + 1
                    )
                    break
            assert common_gcd is not None and common_gcd.degree() == 0
            continue

        if rank == 1:
            assert line_direction is not None and plane_basis is not None
            if sum(entry != 0 for entry in line_direction) == 1:
                rank_one_boundary_planes += 1
                continue
            rank_one_genuine_planes += 1
            direction_hessians = []
            for direction in plane_basis:
                point_hessians = []
                for point_index in range(len(verification_points)):
                    directed_hessian = sp.zeros(VARIABLE_COUNT)
                    for coefficient, cubic_index in zip(
                        direction, cubic_triple, strict=True
                    ):
                        directed_hessian += (
                            coefficient
                            * verification_cubic_hessians[point_index][
                                cubic_index
                            ]
                        )
                    point_hessians.append(directed_hessian)
                direction_hessians.append(point_hessians)

            equations = []
            unit = False
            for point_index in range(len(verification_points)):
                equation = sp.Poly(
                    sp.expand(
                        (
                            base_matrices[point_index]
                            + lambda_parameter
                            * direction_hessians[0][point_index]
                            + mu_parameter
                            * direction_hessians[1][point_index]
                        ).det(method="berkowitz")
                        - 64
                    ),
                    lambda_parameter,
                    mu_parameter,
                    domain=sp.QQ,
                )
                if not equation.is_zero:
                    equations.append(equation.as_expr())
                if point_index < 2:
                    continue
                unit = is_unit_groebner(
                    sp.groebner(
                        equations,
                        lambda_parameter,
                        mu_parameter,
                        domain=sp.QQ,
                        order="grevlex",
                    )
                )
                if unit:
                    maximum_rank_one_points = max(
                        maximum_rank_one_points, point_index + 1
                    )
                    break
            assert unit
            continue

        assert rank == 0
        equations = []
        unit = False
        for point_index in range(len(verification_points)):
            candidate = base_matrices[point_index].copy()
            for parameter, cubic_index in zip(
                (lambda_parameter, mu_parameter, nu_parameter),
                cubic_triple,
                strict=True,
            ):
                candidate += (
                    parameter
                    * verification_cubic_hessians[point_index][cubic_index]
                )
            equation = sp.Poly(
                sp.expand(candidate.det(method="berkowitz") - 64),
                lambda_parameter,
                mu_parameter,
                nu_parameter,
                domain=sp.QQ,
            )
            if not equation.is_zero:
                equations.append(equation.as_expr())
            if point_index < 2:
                continue
            unit = is_unit_groebner(
                sp.groebner(
                    equations,
                    lambda_parameter,
                    mu_parameter,
                    nu_parameter,
                    domain=sp.QQ,
                    order="grevlex",
                )
            )
            if unit:
                maximum_rank_zero_points = max(
                    maximum_rank_zero_points, point_index + 1
                )
                break
        assert unit


assert rank_counts == {
    0: 5_480,
    1: 53_364,
    2: 130_508,
    3: 77_408,
}
assert rank_two_boundary_lines == 129_588
assert rank_two_genuine_lines == 920
assert rank_one_boundary_planes == 50_412
assert rank_one_genuine_planes == 2_952

print(
    "PASS: the exact QQ odd-layer ranks are "
    "r0=5480, r1=53364, r2=130508, r3=77408"
)
print(
    "PASS: all 920 genuine rank-two lines have unit determinant gcd in QQ"
)
print(
    "PASS: all 2952 genuine rank-one planes have unit determinant ideal "
    "in QQ[lambda,mu]"
)
print(
    "PASS: all 5480 rank-zero spaces have unit determinant ideal "
    "in QQ[lambda,mu,nu]"
)
print(
    "DETAIL: maximum evaluation-prefix lengths for ranks 2,1,0 are "
    f"{maximum_rank_two_points},{maximum_rank_one_points},"
    f"{maximum_rank_zero_points}"
)
print(
    "SCOPE: cubic support <=3 is excluded in characteristic zero; cubic "
    "support >=4, dense quartics, and non-coordinate reductions remain open"
)
