#!/usr/bin/env python3
"""Degree-four jet obstruction for unrestricted quadratic-cubic graph shears.

For a selected Lagrangian chart, take the fully general potential

    V(u) = V_2(u) + V_3(u),

with ten quadratic and twenty homogeneous-cubic coefficients.  Expand the
Jacobian determinant of B=m_0+grad V(q_0) through total graph degree four.
Adjoin z*c_0-1, where c_0 is the constant determinant term, to saturate by
the requirement c_0 != 0.  Singular then tests whether the jet coefficient
ideal is the unit ideal.

Seven charts are certified over Q: 0001, 0100, 0101, 1010, and 1111 at
degree four; 1011 and 1110 at degree five.  Other charts can be explored with
the same CLI, including modular runs.
"""

from __future__ import annotations

import argparse
from itertools import permutations, product
import subprocess
import runpy

import sympy as sp


graph_namespace = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
h_variables = graph_namespace["h_variables"]
Q = list(graph_namespace["position_coordinates_h"])
M = list(graph_namespace["momentum_coordinates_h"])

u = sp.symbols("u0:4")
k = sp.symbols("k0:10")
c = sp.symbols("c0:20")
z_saturation = sp.symbols("z_saturation")
parameters = k + c + (z_saturation,)
t = sp.symbols("t")

quadratic_entries = [
    (row, column) for row in range(4) for column in range(row, 4)
]
K = sp.zeros(4)
for coefficient, (row, column) in zip(k, quadratic_entries, strict=True):
    K[row, column] = coefficient
    K[column, row] = coefficient

cubic_monomials = [
    sp.prod(u[index] ** exponents[index] for index in range(4))
    for exponents in product(range(4), repeat=4)
    if sum(exponents) == 3
]
cubic_potential = sum(
    coefficient * monomial
    for coefficient, monomial in zip(c, cubic_monomials, strict=True)
)
cubic_gradient = sp.Matrix(
    [sp.diff(cubic_potential, variable) for variable in u]
)
cubic_hessian = sp.hessian(
    cubic_potential,
    u,
)


def truncated_multiply(
    left: list[sp.Expr], right: list[sp.Expr], order: int
) -> list[sp.Expr]:
    output = [sp.Integer(0)] * order
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            if left_degree + right_degree < order:
                output[left_degree + right_degree] += (
                    left_coefficient * right_coefficient
                )
    return output


def determinant_jet(mask: tuple[int, ...], order: int = 5) -> list[sp.Expr]:
    q0 = sp.Matrix([M[index] if mask[index] else Q[index] for index in range(4)])
    m0 = sp.Matrix([-Q[index] if mask[index] else M[index] for index in range(4)])
    jacobian = m0.jacobian(h_variables) + (
        K
        + cubic_hessian.subs(dict(zip(u, q0, strict=True)), simultaneous=True)
    ) * q0.jacobian(h_variables)
    scaling = {variable: t * variable for variable in h_variables}

    entry_jets = []
    for row in range(4):
        jet_row = []
        for column in range(4):
            scaled_entry = sp.expand(jacobian[row, column].subs(scaling))
            jet_row.append(
                [scaled_entry.coeff(t, degree) for degree in range(order)]
            )
        entry_jets.append(jet_row)

    determinant = [sp.Integer(0)] * order
    for permutation in permutations(range(4)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(4)
            for right in range(left + 1, 4)
        )
        term = [sp.Integer(1)] + [sp.Integer(0)] * (order - 1)
        for row, column in enumerate(permutation):
            term = truncated_multiply(term, entry_jets[row][column], order)
        sign = -1 if inversions % 2 else 1
        for degree in range(order):
            determinant[degree] += sign * term[degree]
    return determinant


def jet_equations(mask: tuple[int, ...], order: int = 5) -> list[sp.Expr]:
    determinant = determinant_jet(mask, order)
    equations = []
    for degree in range(1, order):
        homogeneous_part = sp.Poly(
            sp.expand(determinant[degree]), *h_variables
        )
        equations.extend(
            coefficient for _, coefficient in homogeneous_part.terms()
        )
    equations.append(z_saturation * sp.expand(determinant[0]) - 1)
    return equations


def collision_equations(
    mask: tuple[int, ...], collision: str
) -> list[sp.Expr]:
    if collision == "none":
        return []

    q_values = graph_namespace["position_collision_values"]
    m_values = graph_namespace["momentum_collision_values"]
    projected_values = []
    for point_index in range(3):
        q_value = sp.Matrix(
            [
                m_values[point_index][index]
                if mask[index]
                else q_values[point_index][index]
                for index in range(4)
            ]
        )
        m_value = sp.Matrix(
            [
                -q_values[point_index][index]
                if mask[index]
                else m_values[point_index][index]
                for index in range(4)
            ]
        )
        projected_values.append(
            m_value
            + K * q_value
            + cubic_gradient.subs(
                dict(zip(u, q_value, strict=True)), simultaneous=True
            )
        )

    pairs = {
        "01": ((0, 1),),
        "02": ((0, 2),),
        "12": ((1, 2),),
        "all": ((0, 1), (0, 2)),
    }[collision]
    return [
        sp.expand(equation)
        for left, right in pairs
        for equation in projected_values[right] - projected_values[left]
    ]


def collision_parameterization(
    equations: list[sp.Expr],
) -> tuple[dict[sp.Symbol, sp.Expr] | None, tuple[sp.Symbol, ...]]:
    if not equations:
        return {}, k + c

    coefficient_matrix, constant_vector = sp.linear_eq_to_matrix(
        equations, k + c
    )
    reduced, pivots = coefficient_matrix.row_join(constant_vector).rref()
    parameter_count = len(k + c)
    if parameter_count in pivots:
        return None, ()

    pivot_columns = set(pivots)
    free_indices = [
        index for index in range(parameter_count) if index not in pivot_columns
    ]
    free_parameters = tuple((k + c)[index] for index in free_indices)
    substitution = {}
    for row, pivot in enumerate(pivots):
        substitution[(k + c)[pivot]] = sp.expand(
            reduced[row, -1]
            - sum(
                reduced[row, free_index] * (k + c)[free_index]
                for free_index in free_indices
            )
        )
    return substitution, free_parameters


def singular_expression(
    expression: sp.Expr,
    characteristic: int,
    ring_parameters: tuple[sp.Symbol, ...],
) -> str:
    cleared = sp.Poly(
        expression, *ring_parameters
    ).clear_denoms()[1].as_expr()
    if characteristic:
        cleared = sp.Poly(
            cleared, *ring_parameters, modulus=characteristic
        ).as_expr()
    return str(cleared).replace("**", "^")


def singular_remainder(
    equations: list[sp.Expr],
    ring_parameters: tuple[sp.Symbol, ...],
    characteristic: int,
    timeout: int,
    algorithm: str,
) -> str:
    script = (
        f"ring r={characteristic},({','.join(map(str, ring_parameters))}),dp;\n"
        f"ideal I={','.join(singular_expression(eq, characteristic, ring_parameters) for eq in equations)};\n"
        "option(redSB);\n"
        f"ideal G={algorithm}(I);\n"
        "reduce(1,G);\n"
    )
    result = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    return result.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chart", default="1111")
    parser.add_argument("--characteristic", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="std")
    parser.add_argument("--max-degree", type=int, default=4)
    parser.add_argument(
        "--collision",
        choices=("none", "01", "02", "12", "all"),
        default="none",
        help=(
            "retain a selected pair, or all three points, of the certified "
            "PC(2) collision"
        ),
    )
    parser.add_argument(
        "--collision-reduction",
        choices=("substitute", "equations"),
        default="substitute",
        help=(
            "either eliminate the linear collision equations over Q before "
            "the jet calculation, or give them directly to Singular"
        ),
    )
    args = parser.parse_args()

    if len(args.chart) != 4 or any(bit not in "01" for bit in args.chart):
        raise ValueError("chart must be a four-bit mask")
    mask = tuple(int(bit) for bit in args.chart)
    collision_constraints = collision_equations(mask, args.collision)
    substitution, free_parameters = collision_parameterization(
        collision_constraints
    )
    if substitution is None:
        print(
            f"PASS: chart {args.chart} cannot retain collision "
            f"{args.collision} in the quadratic-cubic shear class"
        )
        return
    equations = jet_equations(mask, args.max_degree + 1)
    if collision_constraints and args.collision_reduction == "substitute":
        equations = [
            substituted
            for equation in equations
            if (
                substituted := sp.expand(
                    equation.subs(substitution, simultaneous=True)
                )
            )
            != 0
        ]
        ring_parameters = free_parameters + (z_saturation,)
    else:
        equations = collision_constraints + equations
        ring_parameters = k + c + (z_saturation,)
    print(
        f"chart {args.chart}: {len(equations)} saturated degree-{args.max_degree} "
        f"jet equations over characteristic {args.characteristic}; "
        f"collision={args.collision}; {len(free_parameters)} shear parameters"
    )
    remainder = singular_remainder(
        equations,
        ring_parameters,
        args.characteristic,
        args.timeout,
        args.algorithm,
    )
    print(f"remainder of 1: {remainder}")
    certified_degrees = {
        "0001": 4,
        "0100": 4,
        "0101": 4,
        "1010": 4,
        "1011": 5,
        "1110": 5,
        "1111": 4,
    }
    if (
        args.characteristic == 0
        and args.chart in certified_degrees
        and args.max_degree >= certified_degrees[args.chart]
    ):
        assert remainder == "0"
        print(
            f"PASS: unrestricted quadratic-cubic potentials in chart "
            f"{args.chart} cannot have constant nonzero determinant"
        )


if __name__ == "__main__":
    main()
