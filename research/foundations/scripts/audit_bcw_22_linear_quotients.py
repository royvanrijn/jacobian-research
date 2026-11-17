#!/usr/bin/env python3
"""Rule out further collision-preserving linear quotients of the 22D BCW map.

For a quotient with row space W to carry an induced map, W must be invariant
under every coefficient matrix of JH.  This script proves that the only
common invariant row spaces are 0, the homogenizing covector, and the full
22-dimensional space.  The one-dimensional option is constant on all stored
collision points, so no proper linear quotient preserves the collision.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts" / "generated-results" / "constant_kernel_bcw_22_counterexample.json"
PRIME = 1_000_003


def decode_h(stored: dict[str, object]) -> tuple[list[sp.Symbol], list[sp.Expr]]:
    variables = list(sp.symbols(f"x0:{stored['dimension']}"))
    components = []
    for component in stored["H"]:
        expression = sp.Integer(0)
        for term in component:
            expression += sp.Rational(term["coefficient"]) * sp.prod(
                variables[index] ** exponent for index, exponent in term["monomial"]
            )
        components.append(sp.expand(expression))
    return variables, components


def jacobian_coefficient_matrices(
    components: list[sp.Expr], variables: list[sp.Symbol]
) -> list[sp.Matrix]:
    jacobian = sp.Matrix(components).jacobian(variables)
    monomials = sorted(
        {
            exponents
            for entry in jacobian
            for exponents, coefficient in sp.Poly(entry, *variables, domain=sp.QQ).terms()
            if coefficient
        }
    )
    result = []
    for exponents in monomials:
        matrix = sp.Matrix(
            [
                [
                    sp.Poly(jacobian[i, j], *variables, domain=sp.QQ).coeff_monomial(
                        exponents
                    )
                    for j in range(len(variables))
                ]
                for i in range(len(components))
            ]
        )
        if matrix != sp.zeros(len(variables)):
            result.append(matrix)
    return result


def modular_flat(matrix: sp.Matrix) -> tuple[int, ...]:
    return tuple(
        int(value.p) * pow(int(value.q), -1, PRIME) % PRIME for value in matrix
    )


def modular_product(left: tuple[int, ...], right: tuple[int, ...], n: int) -> tuple[int, ...]:
    result = [0] * (n * n)
    for i in range(n):
        for k in range(n):
            coefficient = left[i * n + k]
            if coefficient:
                for j in range(n):
                    result[i * n + j] = (
                        result[i * n + j] + coefficient * right[k * n + j]
                    ) % PRIME
    return tuple(result)


def generated_algebra_dimension(generators: list[tuple[int, ...]], n: int) -> int:
    pivots: dict[int, tuple[int, ...]] = {}
    queue: list[tuple[int, ...]] = []

    def add(vector: tuple[int, ...]) -> bool:
        reduced = list(vector)
        for pivot in sorted(pivots):
            if reduced[pivot]:
                factor = reduced[pivot]
                reduced = [
                    (left - factor * right) % PRIME
                    for left, right in zip(reduced, pivots[pivot])
                ]
        try:
            pivot = next(index for index, value in enumerate(reduced) if value)
        except StopIteration:
            return False
        inverse = pow(reduced[pivot], -1, PRIME)
        normalized = tuple(value * inverse % PRIME for value in reduced)
        pivots[pivot] = normalized
        queue.append(normalized)
        return True

    identity = tuple(1 if i == j else 0 for i in range(n) for j in range(n))
    add(identity)
    for generator in generators:
        add(generator)
    head = 0
    while head < len(queue) and len(pivots) < n * n:
        element = queue[head]
        head += 1
        for generator in generators:
            add(modular_product(element, generator, n))
            if len(pivots) == n * n:
                break
    return len(pivots)


def main() -> None:
    stored = json.loads(SOURCE.read_text())
    variables, h = decode_h(stored)
    assert len(variables) == 22
    matrices = jacobian_coefficient_matrices(h, variables)
    assert len(matrices) == 66

    # L=<e_21^T> is invariant because the last component of H is zero.
    last_row = sp.zeros(1, 22)
    last_row[0, 21] = 1
    assert all(last_row * matrix == sp.zeros(1, 22) for matrix in matrices)
    points = [sp.Matrix([sp.Rational(value) for value in point]) for point in stored["collision_points"]]
    assert len({point[21] for point in points}) == 1

    # On V/L the generated coefficient algebra is all M_21.  Full rank modulo
    # one good prime certifies full rank over Q, hence V/L is a simple module.
    quotient_generators = [modular_flat(matrix[:21, :21]) for matrix in matrices]
    algebra_dimension = generated_algebra_dimension(quotient_generators, 21)
    assert algebra_dimension == 21 * 21

    # A 21D invariant row space would be an invariant hyperplane.  Its
    # annihilator would be a common invariant column line.  Of the 66
    # coefficient matrices, 65 are nilpotent; on such a line their eigenvalue
    # must be zero.  Their common kernel is zero, excluding that line.
    nilpotent = [matrix for matrix in matrices if matrix**22 == sp.zeros(22)]
    exceptional = [matrix for matrix in matrices if matrix**22 != sp.zeros(22)]
    assert len(nilpotent) == 65 and len(exceptional) == 1
    assert sp.Matrix.vstack(*nilpotent).nullspace() == []
    assert sp.factor(exceptional[0].charpoly().as_expr()) == sp.Symbol("lambda") ** 20 * (
        sp.Symbol("lambda") - 6
    ) * (sp.Symbol("lambda") + 6)

    print("PASS 22D quotient audit: 66 exact Jacobian coefficient matrices")
    print("PASS 22D quotient audit: quotient algebra mod <e_21^T> is full M_21(Q)")
    print("PASS 22D quotient audit: no common invariant column line")
    print("PASS 22D quotient audit: only proper row module cannot separate the collision")
    print("PASS 22D quotient audit: no further collision-preserving linear quotient")


if __name__ == "__main__":
    main()
