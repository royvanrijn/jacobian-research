#!/usr/bin/env python3
"""Replay the Laurent-mode coefficients of the zero-parameter mixed pencil."""

from sympy import Matrix, Poly, Rational, expand, symbols

from verify_quintic_rank_two_tangent_double_conic import determinant_piece, r, x


def laurent_piece(pieces, drop):
    return expand(x*determinant_piece(pieces, drop))


s, t, a = symbols("s t a")
leading = Matrix([r**2, r, r])
direction = Matrix([2*r, 1, 1])
second_direction = Matrix([r, Rational(-1, 3), 0])
extra_direction = Matrix([0, 5*r/3, r])
constant_direction = direction - 2*second_direction
zero = Matrix.zeros(3, 1)
euler = Matrix([5*r**2, 5*r, 4*r])
assert euler.cross(leading.diff(r)) == Matrix([r, 3*r**2, -5*r**2])

quadratic = 3*s**2 + 3*s*t + 2*t**2
cubic = s*(2*s + t)*(3*s - t)
kernel_mode = s*direction + t*second_direction

assert laurent_piece([leading, s*x**3*constant_direction], 2) == -10*r*s**2*x**6
assert expand(laurent_piece([leading, x**2*kernel_mode], 2)
              + 4*r*x**4*quadratic/3) == 0
assert expand(laurent_piece([leading, x**2*kernel_mode], 3)
              + 2*x**6*cubic/3) == 0

first_mode = x*(a*extra_direction + t*second_direction)
first_quadratic = Poly(laurent_piece([leading, first_mode], 2), x, r)
assert first_quadratic.coeff_monomial(x**2*r**3) == -10*a**2/3
assert laurent_piece([leading, t*x*second_direction], 2) == -4*r*t**2*x**2/3
second_mode = x**2*Matrix([2*t**2/3, 5*s/3, s])
assert laurent_piece([leading, t*x*second_direction, second_mode], 2) == 0
assert expand(laurent_piece([leading, t*x*second_direction, second_mode], 3)
              - 5*t*x**3*(9*r*s - 2*t**2)/9) == 0

assert laurent_piece([leading, zero, s*x**2*constant_direction], 4) == -20*r*s**2*x**4/3
assert expand(laurent_piece([leading, zero, x*kernel_mode], 4)
              + 2*r*x**2*quadratic/3) == 0
assert expand(laurent_piece([leading, zero, x*kernel_mode], 6)
              + x**3*cubic/3) == 0
assert laurent_piece([leading, zero, zero, s*x*constant_direction], 6) == -10*r*s**2*x**2/3

print("Zero mixed reduced-line pencil: exact Laurent-mode coefficients verified.")
