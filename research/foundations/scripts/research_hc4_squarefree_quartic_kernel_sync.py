#!/usr/bin/env python3
"""Explore constant-polar synchronization on four-line quintic packets.

This is an experiment, not a theorem checker.  For a fixed arrangement of
four lines and a small projective atlas of constant directions v_i, it solves

    D_{v_i}(h_5) in (L_i^2)  if v_i is tangent to L_i,
    D_{v_i}(h_5) in (L_i^3)  if v_i is transverse to L_i.

The first condition gives a constant kernel of Hess(h_5) along L_i.  The
extra transverse order is equivalent, on the generic-corank-one locus, to
the second determinant factor.  Thus the synchronized system implies P^2
divides det(Hess(h_5)).
"""

from __future__ import annotations

import argparse
import itertools

import sympy as sp


x, y, z, n = sp.symbols("x y z n")
variables = (x, y, z)
monomial_exponents = [
    (a, b, 5 - a - b)
    for a in range(6)
    for b in range(6 - a)
]
monomials = [x**a * y**b * z**c for a, b, c in monomial_exponents]
coefficients = sp.symbols(f"c0:{len(monomials)}")
h_five = sum(c * m for c, m in zip(coefficients, monomials))


ARRANGEMENTS = {
    "general": (x, y, z, x + y + z),
    "triple": (x, y, x + y, z),
    "pencil": (x, y, x + y, x - y),
}


def projective_directions(prime: int | None = None) -> list[tuple[int, int, int]]:
    """Return a small characteristic-zero atlas or all points over F_p."""
    if prime is not None:
        return (
            [(1, second, third) for second in range(prime) for third in range(prime)]
            + [(0, 1, third) for third in range(prime)]
            + [(0, 0, 1)]
        )
    directions: list[tuple[int, int, int]] = []
    for vector in itertools.product((-1, 0, 1), repeat=3):
        if vector == (0, 0, 0):
            continue
        first = next(entry for entry in vector if entry)
        normalized = tuple(first * entry for entry in vector)
        if normalized not in directions:
            directions.append(normalized)
    return directions


def normal_substitution(line: sp.Expr) -> tuple[dict[sp.Symbol, sp.Expr], tuple[sp.Symbol, ...]]:
    """Choose a linear chart in which line is the normal coordinate n."""
    if line == x:
        return {x: n}, (n, y, z)
    if line == y:
        return {y: n}, (x, n, z)
    if line == z:
        return {z: n}, (x, y, n)
    if line == x + y + z:
        return {z: n - x - y}, (x, y, n)
    if line == x + y:
        return {y: n - x}, (x, n, z)
    if line == x - y:
        return {y: x - n}, (x, n, z)
    raise ValueError(f"unsupported line: {line}")


def constraint_block(
    line: sp.Expr,
    direction: tuple[int, int, int],
    prime: int | None = None,
) -> sp.Matrix:
    derivative = sum(
        entry * sp.diff(h_five, variable)
        for entry, variable in zip(direction, variables)
    )
    substitution, chart_variables = normal_substitution(line)
    transformed = sp.Poly(sp.expand(derivative.subs(substitution)), *chart_variables)
    line_on_direction = int(line.subs(dict(zip(variables, direction))))
    tangent = line_on_direction == 0 if prime is None else line_on_direction % prime == 0
    required_order = 2 if tangent else 3
    rows: list[list[int]] = []
    for term, expression in transformed.terms():
        normal_exponent = term[chart_variables.index(n)]
        if normal_exponent >= required_order:
            continue
        rows.append([int(sp.expand(expression).coeff(c)) for c in coefficients])
    return sp.Matrix(rows)


def extend_echelon(
    basis: tuple[tuple[int, list[int]], ...],
    rows: tuple[tuple[int, ...], ...],
    prime: int = 1000003,
) -> tuple[tuple[int, list[int]], ...]:
    """Increment a modular row-echelon basis without rebuilding its matrix."""
    extended = [(pivot, row.copy()) for pivot, row in basis]
    for input_row in rows:
        row = [int(entry) % prime for entry in input_row]
        for pivot, basis_row in extended:
            multiplier = row[pivot]
            if multiplier:
                row = [
                    (left - multiplier * right) % prime
                    for left, right in zip(row, basis_row)
                ]
        pivot = next((column for column, entry in enumerate(row) if entry), None)
        if pivot is None:
            continue
        inverse = pow(row[pivot], prime - 2, prime)
        row = [(entry * inverse) % prime for entry in row]
        extended.append((pivot, row))
        extended.sort(key=lambda item: item[0])
        if len(extended) == len(monomials):
            break
    return tuple(extended)


def modular_hessian_witness(
    basis: tuple[tuple[int, list[int]], ...],
    prime: int = 1000003,
) -> list[int] | None:
    """Find a small modular null vector whose Hessian is nonzero somewhere."""
    pivots = {pivot for pivot, _ in basis}
    free_columns = [column for column in range(len(monomials)) if column not in pivots]
    points = ((2, 3, 5), (1, 4, 7), (3, 8, 2))
    assignments = (
        tuple(index + 1 for index in range(len(free_columns))),
        tuple(1 for _ in free_columns),
        tuple((index + 1) ** 2 for index in range(len(free_columns))),
    )
    for assignment in assignments:
        vector = [0] * len(monomials)
        for column, value in zip(free_columns, assignment):
            vector[column] = value % prime
        for pivot, row in reversed(basis):
            vector[pivot] = -sum(
                row[column] * vector[column]
                for column in range(pivot + 1, len(monomials))
            ) % prime
        for point in points:
            numeric = [[0] * 3 for _ in range(3)]
            for coefficient, exponents in zip(vector, monomial_exponents):
                for first in range(3):
                    for second in range(first, 3):
                        multiplier = exponents[first] * (
                            exponents[second] - int(first == second)
                        )
                        if not multiplier:
                            continue
                        value = coefficient * multiplier
                        for index, (coordinate, exponent) in enumerate(zip(point, exponents)):
                            value *= coordinate ** (
                                exponent
                                - int(index == first)
                                - int(index == second)
                            )
                        numeric[first][second] = (
                            numeric[first][second] + value
                        ) % prime
                        numeric[second][first] = numeric[first][second]
            a, b, c = numeric[0]
            _, d, e = numeric[1]
            _, _, f = numeric[2]
            determinant = a * d * f + 2 * b * c * e - a * e * e - d * c * c - f * b * b
            if determinant % prime:
                return vector
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arrangement", choices=ARRANGEMENTS, default="general")
    parser.add_argument("--show", type=int, default=20)
    parser.add_argument(
        "--finite-field-prime",
        type=int,
        help="exhaust all projective kernel directions over this prime field",
    )
    args = parser.parse_args()

    lines = ARRANGEMENTS[args.arrangement]
    field_prime = args.finite_field_prime
    directions = projective_directions(field_prime)
    calculation_prime = field_prime or 1000003
    exact_blocks = {
        (index, direction): constraint_block(line, direction, field_prime)
        for index, line in enumerate(lines)
        for direction in directions
    }
    blocks = {
        key: tuple(tuple(int(entry) for entry in row) for row in block.tolist())
        for key, block in exact_blocks.items()
    }

    survivors: list[
        tuple[
            int,
            tuple[tuple[int, int, int], ...],
            tuple[tuple[int, list[int]], ...],
        ]
    ] = []
    determinant_survivors: list[
        tuple[
            int,
            tuple[tuple[int, int, int], ...],
            tuple[tuple[int, list[int]], ...],
        ]
    ] = []
    rank_histogram: dict[int, int] = {}

    def visit(
        index: int,
        direction_tuple: tuple[tuple[int, int, int], ...],
        basis: tuple[tuple[int, list[int]], ...],
    ) -> None:
        if len(basis) == len(monomials):
            completions = len(directions) ** (4 - index)
            rank_histogram[len(monomials)] = (
                rank_histogram.get(len(monomials), 0) + completions
            )
            return
        if index == 4:
            rank = len(basis)
            rank_histogram[rank] = rank_histogram.get(rank, 0) + 1
            packet = (rank, direction_tuple, basis)
            survivors.append(packet)
            if modular_hessian_witness(basis, calculation_prime) is not None:
                determinant_survivors.append(packet)
            return
        for direction in directions:
            visit(
                index + 1,
                direction_tuple + (direction,),
                extend_echelon(basis, blocks[index, direction], calculation_prime),
            )

    visit(0, (), ())

    print(f"ARRANGEMENT: {args.arrangement}")
    print(f"FIELD: {'Q-atlas' if field_prime is None else f'F_{field_prime}'}")
    print(f"DIRECTIONS: {len(directions)}")
    print(f"PACKETS TESTED: {len(directions) ** 4}")
    print(f"MODULAR RANK HISTOGRAM: {sorted(rank_histogram.items())}")
    print(f"SYNCHRONIZED QUINTIC SPACES: {len(survivors)}")
    print(f"MODULAR NONZERO-HESSIAN WITNESSES: {len(determinant_survivors)}")
    print("WITNESS TEST: three nullspace assignments at three source points per space")
    if field_prime is not None:
        for rank, direction_tuple, _ in determinant_survivors[: args.show]:
            print("FINITE-FIELD SURVIVOR", f"rank={rank}", f"directions={direction_tuple}")
        return
    denominator = sp.prod(lines)
    nonzero_determinant_count = 0
    displayed_packets = determinant_survivors or survivors
    for rank, direction_tuple, _ in sorted(displayed_packets)[: args.show]:
        stacked = sp.Matrix.vstack(
            *(exact_blocks[index, direction] for index, direction in enumerate(direction_tuple))
        )
        nullspace = stacked.nullspace()
        basis_forms = [
            sp.factor(sum(entry * monomial for entry, monomial in zip(vector, monomials)))
            for vector in nullspace
        ]
        test_form = sum(
            (index + 1) * form for index, form in enumerate(basis_forms)
        )
        determinant = sp.factor(sp.hessian(test_form, variables).det())
        quotient = sp.cancel(determinant / denominator**2) if determinant else 0
        if determinant:
            assert sp.denom(quotient) == 1
            nonzero_determinant_count += 1
        print(
            "SURVIVOR",
            f"rank={rank}",
            f"directions={direction_tuple}",
            f"basis={basis_forms}",
            f"test_quotient={sp.factor(quotient)}",
        )
    print(
        "NONZERO TEST DETERMINANTS AMONG DISPLAYED SPACES:",
        nonzero_determinant_count,
    )


if __name__ == "__main__":
    main()
