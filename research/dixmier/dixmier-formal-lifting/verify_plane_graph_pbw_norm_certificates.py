#!/usr/bin/env python3
"""Replay the exact reduced-norm filtration procedure on prime-field graphs."""

from sympy import GF
from sympy.polys.rings import ring


def diagonal(size, value, zero):
    return [
        [value if row == col else zero for col in range(size)]
        for row in range(size)
    ]


def multiply(left, right, zero):
    size = len(left)
    return [
        [
            sum((left[i][j] * right[j][k] for j in range(size)), zero)
            for k in range(size)
        ]
        for i in range(size)
    ]


def powers(matrix, size, zero, one):
    output = [diagonal(size, one, zero)]
    for _ in range(size - 1):
        output.append(multiply(output[-1], matrix, zero))
    return output


def determinant(matrix, zero, one):
    """Bareiss elimination in the polynomial ring over the prime field."""
    work = [row[:] for row in matrix]
    sign, denominator = one, one
    size = len(work)
    for i in range(size - 1):
        pivot_row = next((j for j in range(i, size) if work[j][i]), None)
        if pivot_row is None:
            return zero
        if pivot_row != i:
            work[i], work[pivot_row] = work[pivot_row], work[i]
            sign = -sign
        pivot = work[i][i]
        for row in range(i + 1, size):
            for col in range(i + 1, size):
                numerator = pivot * work[row][col] - work[row][i] * work[i][col]
                work[row][col] = numerator.exquo(denominator)
            work[row][i] = zero
        denominator = pivot
    return sign * work[-1][-1]


def symbol(matrix, prime, zero, one):
    """Coordinates are source p-roots; Frobenius is the identity on GF(p)."""
    norm = determinant(matrix, zero, one)
    if not norm:
        # The graph algebra is a domain, so zero norm must be the zero row.
        assert not any(entry for row in matrix for entry in row)
        return {}
    assert all(i % prime == j % prime == 0 for i, j in norm)
    degree = max(i + j for i, j in norm)
    return {
        (i // prime, j // prime): c
        for (i, j), c in norm.items()
        if i + j == degree
    }


def filtration(prime, cutoff, get_map):
    center, x0, y0 = ring("x0, y0", GF(prime))
    zero, one = center.zero, center.one
    f1, f2 = get_map(x0, y0)
    jacobian = f1.diff(x0) * f2.diff(y0) - f1.diff(y0) * f2.diff(x0)
    assert jacobian == one
    d1, d2 = (max(i + j for i, j in image) for image in (f1, f2))
    assert max(d1, d2) < prime

    u, v = diagonal(prime, f1, zero), diagonal(prime, f2, zero)
    for j in range(1, prime):
        u[j][j - 1] = one
        v[j - 1][j] = center.ground_new(j)
    uv, vu = multiply(u, v, zero), multiply(v, u, zero)
    assert [
        [vu[i][j] - uv[i][j] for j in range(prime)] for i in range(prime)
    ] == diagonal(prime, one, zero)
    assert symbol(u, prime, zero, one) == {
        monomial: coefficient
        for monomial, coefficient in f1.items()
        if sum(monomial) == d1
    }
    assert symbol(v, prime, zero, one) == {
        monomial: coefficient
        for monomial, coefficient in f2.items()
        if sum(monomial) == d2
    }
    u_powers, v_powers = powers(u, prime, zero, one), powers(v, prime, zero, one)
    assert multiply(u_powers[-1], u, zero) == diagonal(prime, f1**prime, zero)
    assert multiply(v_powers[-1], v, zero) == diagonal(prime, f2**prime, zero)
    complement = prime - 1
    pivots = {}
    for i in range(prime):
        for j in range(prime):
            allowance = cutoff + (complement - i) * d1 + (complement - j) * d2 - 2 * complement
            center_degree = allowance // prime
            if center_degree < 0:
                continue
            pbw = multiply(u_powers[i], v_powers[j], zero)
            for r in range(center_degree + 1):
                for s in range(center_degree - r + 1):
                    monomial = x0 ** (prime * r) * y0 ** (prime * s)
                    row = [[monomial * entry for entry in line] for line in pbw]
                    eliminate(row, pivots, prime, zero, one)
    increments = [0] * (cutoff + 1)
    for degree, _ in pivots:
        if degree <= cutoff:
            increments[degree] += 1
    return increments


def eliminate(row, pivots, prime, zero, one):
    previous = None
    while head := symbol(row, prime, zero, one):
        monomial = max(head)
        pivot = sum(monomial), monomial
        if previous is not None:
            assert pivot < previous
        previous = pivot
        coefficient = head[monomial]
        if pivot not in pivots:
            inverse = coefficient ** (prime - 2)
            pivots[pivot] = [[inverse * entry for entry in line] for line in row]
            return
        stored = pivots[pivot]
        row = [
            [entry - coefficient * other for entry, other in zip(line, stored_line)]
            for line, stored_line in zip(row, stored)
        ]


def main():
    tests = (
        (2, 2, lambda x, y: (x, y)),
        (3, 3, lambda x, y: (x, y + x**2)),
        (5, 3, lambda x, y: (x + y**2, y + (x + y**2)**2)),
    )
    for prime, cutoff, get_map in tests:
        answer = filtration(prime, cutoff, get_map)
        assert answer == list(range(1, cutoff + 2))
        print(f"p={prime}, through degree {cutoff}: increments {answer}")


if __name__ == "__main__":
    main()
