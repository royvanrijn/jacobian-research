#!/usr/bin/env python3
"""Replay Laurent-mode witnesses for the reducible square conic pencil."""

from sympy import Matrix, Rational, expand, symbols

from verify_quintic_rank_two_tangent_double_conic import determinant_piece, r, x


def laurent_piece(pieces, drop):
    return expand(x*determinant_piece(pieces, drop))


leading = Matrix([r**2, r, r**2])
U = Matrix([r, 3, 0])
W = Matrix([Rational(5, 6), 0, 1])
zero = Matrix.zeros(3, 1)
a, s, t, u = symbols("a s t u")

assert laurent_piece([leading, a*x**2*U], 2) == 24*r**2*a**2*x**4
assert laurent_piece([leading, a*x*U], 2) == 12*r**2*a**2*x**2

second = Matrix([5*s/6 - a**2, 0, s])*x**2
pieces = [leading, a*x*U, second]
assert laurent_piece(pieces, 2) == 0
assert expand(laurent_piece(pieces, 3) + 15*a*r*x**3*(s + 2*a**2)) == 0
assert expand(laurent_piece(pieces, 4) - 5*s*x**4*(s - 12*a**2)/3) == 0

assert laurent_piece([leading, zero, s*x**2*W], 4) == 5*s**2*x**4/3
assert expand(laurent_piece([leading, zero, x*(u*U + t*W)], 4)
              - x**2*(72*r**2*u**2 - 60*r*t*u + 5*t**2)/6) == 0
assert laurent_piece([leading, zero, zero, s*x*W], 6) == 5*s**2*x**2/6

print("Reducible square conic: exact Laurent-mode witnesses verified.")
