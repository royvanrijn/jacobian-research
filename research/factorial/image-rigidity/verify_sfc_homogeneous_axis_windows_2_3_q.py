#!/usr/bin/env python3
"""Check the fractional mark and, on one side, every ordinary degree."""

import argparse
from math import factorial, prod

from sympy import Poly, QQ, ZZ, ground_roots, symbols
from sympy.polys.matrices import DomainMatrix

D = symbols("D")
R = ZZ.poly_ring(D)
degree_parameter = R.gens[0]


def monomials(degree):
    for i in range(degree + 1):
        for j in range(degree - i + 1):
            yield (i, j, degree - i - j)


def cleared_macaulay(offsets, opposite_third, terminal):
    maximum = max(offsets)
    critical = terminal + 3
    target = list(monomials(critical))
    row_for = {alpha: row for row, alpha in enumerate(target)}
    columns = []
    for order in (2, 3, terminal):
        terms = []
        for alpha in monomials(order):
            if opposite_third and 0 < alpha[2] < order:
                continue
            count = sum(a * b for a, b in zip(alpha, offsets))
            multinomial = factorial(order) // prod(factorial(a) for a in alpha)
            # Clear the common row denominator (order*D)_(order*maximum).
            coefficient = R(multinomial * factorial(count))
            for j in range(count, order * maximum):
                coefficient *= order * degree_parameter - j
            terms.append((alpha, coefficient))
        for multiplier in monomials(critical - order):
            column = [R.zero] * len(target)
            for alpha, coefficient in terms:
                row = row_for[tuple(a + b for a, b in zip(alpha, multiplier))]
                column[row] = coefficient
            columns.append(column)
    rows = [list(row) for row in zip(*columns)]
    return DomainMatrix(rows, (len(target), len(columns)), R)


def evaluate(matrix, degree):
    evaluated = [
        [coefficient.evaluate(degree_parameter, R.domain.convert(degree))
         for coefficient in row]
        for row in matrix.to_list()
    ]
    return DomainMatrix(evaluated, matrix.shape, QQ)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("offsets", type=int, nargs=3)
    parser.add_argument("--terminal", type=int, default=4)
    parser.add_argument("--opposite-third", action="store_true")
    args = parser.parse_args()
    offsets = tuple(args.offsets)
    on_one_side = offsets[:2] if args.opposite_third else offsets
    if min(offsets) < 0 or len(set(on_one_side)) != len(on_one_side):
        parser.error("Offsets must be nonnegative and distinct on each side.")
    if args.terminal < 4:
        parser.error("The terminal order must be at least four.")

    matrix = cleared_macaulay(offsets, args.opposite_third, args.terminal)
    mark = QQ(args.terminal * max(offsets) - 1, args.terminal)
    # Evaluate via QQ, since the mark need not be an integer.
    field_ring = QQ.poly_ring(D)
    fractional = [
        [field_ring.convert(coefficient, R).evaluate(field_ring.gens[0], mark)
         for coefficient in row]
        for row in matrix.to_list()
    ]
    _, pivots = DomainMatrix(fractional, matrix.shape, QQ).rref()
    if len(pivots) != matrix.shape[0]:
        raise ArithmeticError("Fractional-degree specialization lost rank.")
    print(f"Fractional mark {mark} has full rank {len(pivots)}.", flush=True)
    if args.opposite_third:
        return

    minor = matrix.extract(list(range(matrix.shape[0])), list(pivots)).det()
    determinant = Poly(R.to_sympy(minor), D)
    failures = [
        int(root) for root in ground_roots(determinant)
        if root.q == 1 and root >= max(offsets)
        and evaluate(matrix, int(root)).rank() != matrix.shape[0]
    ]
    print(f"Failing ordinary degrees: {sorted(failures)}.")


if __name__ == "__main__":
    main()
