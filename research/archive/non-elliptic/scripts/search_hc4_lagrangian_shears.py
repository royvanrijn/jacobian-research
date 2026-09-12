#!/usr/bin/env python3
"""Search collision-preserving Lagrangian linear projections of the PC(2) graph.

Write the ambient graph coordinates as four Darboux pairs (Q,M), with

    Q=(x,q,R,T),  M=(-p,-z,D,S).

Every symmetric matrix K defines another Lagrangian momentum projection

    B = M + K Q.

The two nontrivial points in the certified collision differ in Q by a
multiple of (3,2,0,0) and have equal M.  Thus that collision survives exactly
when K(3,2,0,0)^t=0.  The symmetric integral matrices satisfying this are
parameterized by six integers:

    [[ 4a,-6a, 2b, 2c],
     [-6a, 9a,-3b,-3c],
     [ 2b,-3b,  d,  e],
     [ 2c,-3c,  e,  f]].

This script exhausts a user-selected coefficient box.  A modular evaluation
filter rejects determinants that vary with the graph parameter, and every
survivor is then verified symbolically over Q.
"""

from __future__ import annotations

import argparse
from itertools import permutations, product
import runpy

import sympy as sp


def rational_mod(value: sp.Expr, prime: int) -> int:
    value = sp.Rational(value)
    return (int(value.p) % prime) * pow(int(value.q) % prime, -1, prime) % prime


def determinant_mod(matrix: list[list[int]], prime: int) -> int:
    total = 0
    for permutation in permutations(range(4)):
        inversions = sum(
            permutation[i] > permutation[j] for i in range(4) for j in range(i + 1, 4)
        )
        term = 1
        for row, column in enumerate(permutation):
            term = term * matrix[row][column] % prime
        total = total - term if inversions % 2 else total + term
    return total % prime


def k_matrix(parameters: tuple[int, ...]) -> sp.Matrix:
    a, b, c, d, e, f = parameters
    return sp.Matrix(
        [
            [4 * a, -6 * a, 2 * b, 2 * c],
            [-6 * a, 9 * a, -3 * b, -3 * c],
            [2 * b, -3 * b, d, e],
            [2 * c, -3 * c, e, f],
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=2)
    parser.add_argument("--prime", type=int, default=1_000_003)
    args = parser.parse_args()

    namespace = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
    h_variables = namespace["h_variables"]
    q_coordinates = sp.Matrix(namespace["position_coordinates_h"])
    m_coordinates = sp.Matrix(namespace["momentum_coordinates_h"])
    jq = q_coordinates.jacobian(h_variables)
    jm = m_coordinates.jacobian(h_variables)

    sample_points = (
        (0, 0, 0, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 1),
        (1, 1, 1, 1),
        (-1, 2, -2, 1),
        (2, -1, 1, -2),
        (3, 1, -1, 2),
    )
    evaluated = []
    for point in sample_points:
        substitution = dict(zip(h_variables, point, strict=True))
        jq_value = [
            [rational_mod(jq[i, j].subs(substitution), args.prime) for j in range(4)]
            for i in range(4)
        ]
        jm_value = [
            [rational_mod(jm[i, j].subs(substitution), args.prime) for j in range(4)]
            for i in range(4)
        ]
        evaluated.append((jq_value, jm_value))

    modular_survivors: list[tuple[int, ...]] = []
    parameter_range = range(-args.bound, args.bound + 1)
    for parameters in product(parameter_range, repeat=6):
        k_value = [
            [int(entry) % args.prime for entry in row]
            for row in k_matrix(parameters).tolist()
        ]
        determinants = []
        for jq_value, jm_value in evaluated:
            matrix = [
                [
                    (
                        jm_value[i][j]
                        + sum(k_value[i][r] * jq_value[r][j] for r in range(4))
                    )
                    % args.prime
                    for j in range(4)
                ]
                for i in range(4)
            ]
            determinants.append(determinant_mod(matrix, args.prime))
        if determinants[0] != 0 and len(set(determinants)) == 1:
            modular_survivors.append(parameters)

    symbolic_survivors = []
    for parameters in modular_survivors:
        k_value = k_matrix(parameters)
        determinant = sp.factor((jm + k_value * jq).det(method="berkowitz"))
        if determinant != 0 and not determinant.free_symbols:
            symbolic_survivors.append((parameters, determinant))

    total = (2 * args.bound + 1) ** 6
    print(f"coefficient box: [-{args.bound},{args.bound}]^6 ({total} matrices)")
    print(f"modular survivors: {len(modular_survivors)}")
    print(f"exact constant-nonzero survivors: {len(symbolic_survivors)}")
    for parameters, determinant in symbolic_survivors:
        print(f"parameters={parameters}, determinant={determinant}")


if __name__ == "__main__":
    main()
