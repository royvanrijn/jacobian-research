#!/usr/bin/env python3
"""Close the rank-zero four-cubic coefficient spaces over F_1000003.

For each of the 234 quartic principal parts, identify the cubic monomials
whose determinant-degree-seven and degree-one signatures vanish.  Their
four-subsets are exactly the 5,430 rank-zero cases left by
``verify_hc4_meng_four_cubic_rank_gate.py``.  This checker forms the full
four-parameter determinant evaluation ideal for each case and enlarges the
evaluation prefix until its Groebner basis is {1}.
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

PRIME: int = parent["PRIME"]
quartics = parent["all_principal_quartics_mod"]
cubic_exponents = parent["cubic_exponents"]
sample_points = parent["sample_points"]
points = parent["two_cubic_points"]
quartic_monomial_hessians = parent["two_cubic_quartic_hessians_mod"]
cubic_monomial_hessians = parent["two_cubic_cubic_hessians_mod"]
base_hessian = parent["base_hessian_mod"]
determinant_mod = parent["determinant_mod"]
interpolate_degree_four_mod = parent["interpolate_degree_four_mod"]
scaling_coefficient_mod = parent["scaling_coefficient_mod"]


def add_scaled(target, scale, contribution) -> None:
    for row in range(4):
        for column in range(4):
            target[row][column] = (
                target[row][column]
                + scale * contribution[row][column]
            ) % PRIME


parameters = sp.symbols("lambda mu nu xi")
degree_point_indices = tuple(range(4, len(sample_points)))
zero_set_sizes = {6: 0, 12: 0}
rank_zero_count = 0
maximum_points = 0
survivors = []

for quartic_index, (support, coefficients) in enumerate(quartics):
    quartic_hessians = []
    base_matrices = []
    for point_index in range(len(points)):
        quartic_hessian = [[0] * 4 for _ in range(4)]
        for coefficient, monomial_index in zip(
            coefficients, support, strict=True
        ):
            add_scaled(
                quartic_hessian,
                coefficient,
                quartic_monomial_hessians[point_index][monomial_index],
            )
        quartic_hessians.append(quartic_hessian)
        base_matrices.append(
            [
                [
                    (
                        base_hessian[row][column]
                        + quartic_hessian[row][column]
                    )
                    % PRIME
                    for column in range(4)
                ]
                for row in range(4)
            ]
        )

    zero_cubics = []
    for cubic_index in range(len(cubic_exponents)):
        signature = []
        for point_index in degree_point_indices:
            values = []
            for coefficient in range(5):
                values.append(
                    determinant_mod(
                        [
                            [
                                (
                                    quartic_hessians[point_index][row][column]
                                    + coefficient
                                    * cubic_monomial_hessians[point_index][
                                        cubic_index
                                    ][row][column]
                                )
                                % PRIME
                                for column in range(4)
                            ]
                            for row in range(4)
                        ]
                    )
                )
            polynomial = interpolate_degree_four_mod(values)
            signature.append(
                polynomial[1] if len(polynomial) > 1 else 0
            )
        for point_index in degree_point_indices:
            signature.append(
                scaling_coefficient_mod(
                    quartic_hessians[point_index],
                    cubic_monomial_hessians[point_index][cubic_index],
                    1,
                )
            )
        if not any(signature):
            zero_cubics.append(cubic_index)

    assert len(zero_cubics) in zero_set_sizes
    zero_set_sizes[len(zero_cubics)] += 1

    for cubic_quadruple in combinations(zero_cubics, 4):
        rank_zero_count += 1
        equations = []
        unit = False
        for point_index in range(len(points)):
            candidate = sp.Matrix(base_matrices[point_index])
            for parameter, cubic_index in zip(
                parameters, cubic_quadruple, strict=True
            ):
                candidate += parameter * sp.Matrix(
                    cubic_monomial_hessians[point_index][cubic_index]
                )
            equation = sp.Poly(
                sp.expand(candidate.det(method="berkowitz") - 64),
                *parameters,
                modulus=PRIME,
            )
            if not equation.is_zero:
                equations.append(equation.as_expr())
            if point_index < 2:
                continue
            basis = sp.groebner(
                equations,
                *parameters,
                modulus=PRIME,
                order="grevlex",
            )
            unit = (
                len(basis.polys) == 1
                and sp.expand(basis.polys[0].as_expr()) == 1
            )
            if unit:
                maximum_points = max(maximum_points, point_index + 1)
                break
        if not unit:
            survivors.append(
                (quartic_index, cubic_quadruple)
            )


assert zero_set_sizes == {6: 230, 12: 4}
assert rank_zero_count == 5_430
assert not survivors

print(
    "PASS: the zero-signature cubic sets have size 6 for 230 quartics "
    "and size 12 for four quartics"
)
print(
    "PASS: all 5430 rank-zero four-parameter determinant ideals are units "
    f"modulo 1000003 using at most {maximum_points} evaluation points"
)
print(
    "SCOPE: together with the rank-gate checker this closes cubic support "
    "<=4 over the certificate field; characteristic-zero promotion remains"
)
