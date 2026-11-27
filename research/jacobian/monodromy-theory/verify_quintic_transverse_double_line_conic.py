#!/usr/bin/env python3
"""Replay Laurent-mode coefficients for the transverse double-line conic restriction."""

from sympy import Matrix, Poly, Rational, expand, symbols

from verify_quintic_rank_two_tangent_double_conic import determinant_piece, r, x


def laurent_piece(pieces, drop):
    return expand(x*determinant_piece(pieces, drop))


b, rho, delta, nu, s, t = symbols("b rho delta nu s t")
zero = Matrix.zeros(3, 1)

pure = Matrix([1, b*r + rho*r**2, r])
pure_direction = Matrix([0, b + 2*rho*r, 1])
assert Matrix([5, 5*(b*r + rho*r**2), 4*r]).cross(pure.diff(r)).applyfunc(expand) == Matrix([
    b*r - 3*rho*r**2, -5, 5*(b + 2*rho*r),
])
assert laurent_piece([pure, x**2*pure_direction], 2) == 20*rho*x**4
assert laurent_piece([pure, zero, x*pure_direction], 4) == 10*rho*x**2

mixed = Matrix([r**2 + delta, r + nu, r])
mixed_direction = Matrix([2*r, 1, 1])
cofactor = Matrix([5*(r**2 + delta), 5*(r + nu), 4*r]).cross(mixed.diff(r)).applyfunc(expand)
assert cofactor == Matrix([r + 5*nu, 3*r**2 - 5*delta, -5*r**2 - 10*nu*r + 5*delta])
assert expand(laurent_piece([mixed, x**2*mixed_direction], 2)
              + 4*x**4*(r + 5*nu)) == 0
assert expand(laurent_piece([mixed, zero, x*mixed_direction], 4)
              + 2*x**2*(r + 5*nu)) == 0

exceptional = mixed.subs(delta, 15*nu**2)
second_direction = Matrix([r - 5*nu, Rational(-1, 3), 0])
constant_direction = mixed_direction - 2*second_direction
assert expand(laurent_piece([exceptional, s*x**3*constant_direction], 2)
              + 10*s**2*(r - 3*nu)*x**6) == 0
assert cofactor.subs(delta, 15*nu**2).dot(second_direction).expand() == 0
quadratic_form = 15*nu*s*(s + t) + r*(3*s**2 + 3*s*t + 2*t**2)
direction = s*mixed_direction + t*second_direction
assert expand(laurent_piece([exceptional, x**2*direction], 2)
              + 4*x**4*quadratic_form/3) == 0
assert expand(laurent_piece([exceptional, x*t*second_direction], 2)
              + 4*r*t**2*x**2/3) == 0
assert expand(laurent_piece([exceptional, zero, x**2*constant_direction], 4)
              + 20*x**4*(r - 3*nu)/3) == 0
assert expand(laurent_piece([exceptional, zero, x*direction], 4)
              + 2*x**2*quadratic_form/3) == 0
assert expand(laurent_piece([exceptional, zero, zero, x*constant_direction], 6)
              + 10*x**2*(r - 3*nu)/3) == 0


amplitude = symbols("amplitude")
pure_extra = Matrix([Rational(10, 3), 5*b*r/3, r])
pure_second = amplitude**2*x**2*Matrix([Rational(35, 9), -5*b**2/(18*rho), -5*b/(18*rho)])
pure_modes = [pure, amplitude*x*pure_extra, pure_second]
assert laurent_piece(pure_modes, 2) == 0
assert expand(laurent_piece(pure_modes, 3)
              - 5*amplitude**3*x**3*(-5*b + 14*r*rho)*(b + 2*r*rho)/(9*rho)) == 0

mixed_extra = Matrix([10*delta/3, 5*(r + 2*nu)/3, r])
assert expand(cofactor.dot(mixed_extra)) == 0
generic_quadratic = Poly(laurent_piece([mixed, amplitude*x*mixed_extra], 2), x, r)
assert generic_quadratic.coeff_monomial(x**2*r**3) == -10*amplitude**2/3

exceptional_extra = amplitude*mixed_extra.subs(delta, 15*nu**2) + t*second_direction
exceptional_quadratic = Poly(laurent_piece([exceptional, x*exceptional_extra], 2), x, r)
assert exceptional_quadratic.coeff_monomial(x**2*r**3) == -10*amplitude**2/3

print("Transverse double line and conic: exact Laurent-mode witnesses verified.")
