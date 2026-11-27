"""Verify the finite weighted slice recursion with exact rational arithmetic."""

from collections import defaultdict
from fractions import Fraction


def bracket(left, right):
    coefficients = defaultdict(Fraction)
    for (i, j), a in left.items():
        for (u, v), b in right.items():
            coefficients[i + u - 1, j + v - 1] += a * b * (i*v - j*u)
    return {monomial: value for monomial, value in coefficients.items() if value}


def solve_order(symbol, solution, order):
    """Return the kernel-free solution and the unsolvable Laurent terms."""
    rhs = defaultdict(Fraction)
    for i in range(1, order+1):
        for monomial, coefficient in bracket(
            symbol.get(i, {}), solution[order-i],
        ).items():
            rhs[monomial] -= coefficient
    next_coefficient = {}
    resonant = {}
    for (i, j), coefficient in sorted(rhs.items()):
        if coefficient:
            u, v = i-5, j-2
            factor = 6*v-3*u
            if factor:
                next_coefficient[u, v] = coefficient / factor
            else:
                resonant[i, j] = coefficient
    return next_coefficient, resonant


def check():
    symbol = {
        0: {(6, 3): 1},
        1: {(3, 2): 1, (10, 4): 1},
        4: {(15, 5): 1},
        9: {(21, 6): 1},
    }
    solution = {0: {(-5, -2): Fraction(1, 3)}}
    assert bracket(symbol[0], solution[0]) == {(0, 0): 1}
    solution[1], resonance = solve_order(symbol, solution, 1)
    assert not resonance
    assert solution[1] == {(-8, -3): Fraction(-2, 9)}
    _, resonance = solve_order(symbol, solution, 2)
    assert resonance == {(1, 0): Fraction(4, 9)}
    symbol[2] = {(7, 3): 1}
    _, resonance = solve_order(symbol, solution, 2)
    assert resonance == {(1, 0): Fraction(1, 9)}
    symbol[2] = {(7, 3): Fraction(4, 3)}
    for order in range(2, 8):
        solution[order], resonance = solve_order(symbol, solution, order)
        assert not resonance, (order, resonance)
    _, resonance = solve_order(symbol, solution, 8)
    assert resonance == {(-3, -2): Fraction(8030, 59049)}
    print("Exact weighted slice residues verified at orders two and eight.")


if __name__ == "__main__":
    check()
