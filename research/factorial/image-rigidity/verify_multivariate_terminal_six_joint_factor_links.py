#!/usr/bin/env python3
"""Verify ordinary terminal-six Gauss-factor links by rational differentiation."""

import sympy as sym

from verify_multivariate_terminal_six_private_angles import require


x, y, z = sym.symbols("x y z")
Q = sym.Rational

alpha_derivative = (
    2 * x * (18 * x**2 + 20 * x + 5)
    / (5 * (x + 1) * (2 * x + 1) * (3 * x + 1) * (4 * x + 1))
)
beta_derivative = (
    y * (5 * y + 2) * (24 * y**2 + 25 * y + 5)
    / (3 * (y + 1) * (2 * y + 1) * (3 * y + 1) * (4 * y + 1) * (5 * y + 1))
)
v_derivative = 3 * x**2 / ((1 + x) * (1 + 3 * x) * (1 + 4 * x))
y_derivative = 6 * y**2 / ((1 + y) * (1 + 4 * y) * (1 + 5 * y))
z_derivative = 3 * y**2 / ((1 + 2 * y) * (1 + 3 * y) * (1 + 5 * y))


def check_nonnegative_on_interval(expression):
    mapped = sym.cancel(expression.subs(x, 2 * z / (1 + z)))
    numerator, denominator = map(lambda p: sym.Poly(p, z), sym.fraction(mapped))
    if denominator.LC() < 0:
        numerator, denominator = -numerator, -denominator
    require(denominator.nth(0) > 0 and min(denominator.all_coeffs()) >= 0,
            ("denominator sign", denominator))
    require(numerator.LC() > 0 and min(numerator.all_coeffs()) >= 0,
            ("numerator sign", numerator))


def verify():
    alpha = (2 * sum(sym.log(1 + k * x) for k in (1, 2, 3)) - 3 * sym.log(1 + 4 * x)) / 10
    beta = (sum(sym.log(1 + k * y) for k in (1, 2, 3, 4)) - 2 * sym.log(1 + 5 * y)) / 6
    defect_five = (sym.log(1 + x) - sym.log(1 + 3 * x)) / 2 + sym.log(1 + 4 * x) / 4
    defect_six_first = (sym.log(1 + y) - sym.log(1 + 4 * y)) / 2 + Q(3, 10) * sym.log(1 + 5 * y)
    defect_six_second = (sym.log(1 + 2 * y) - sym.log(1 + 3 * y)) / 2 + sym.log(1 + 5 * y) / 10
    for function, variable, expected in (
        (alpha, x, alpha_derivative), (beta, y, beta_derivative),
        (defect_five, x, v_derivative), (defect_six_first, y, y_derivative),
        (defect_six_second, y, z_derivative),
    ):
        require(sym.cancel(sym.diff(function, variable) - expected) == 0,
                ("derivative formula", function))
    items = (
        ("outer mean", 2, Q(25, 144), 2 / (1 + 4 * x), Q(5, 2) / (1 + 5 * y), Q(25, 24)),
        ("outer skew", 3, Q(125, 648), alpha_derivative, beta_derivative, Q(125, 108)),
        ("first shape defect", 4, Q(125, 432), v_derivative, y_derivative, Q(125, 108)),
        ("second shape defect", 4, Q(125, 864), v_derivative, z_derivative, Q(125, 216)),
    )
    for name, power, cap, lower, upper, leading in items:
        for gamma, sign in ((0, 1), (1, -1)):
            new_x = 5 * x / (6 - gamma * x)
            derivative = upper.subs(y, new_x) * sym.diff(new_x, x) - leading * lower
            check_nonnegative_on_interval(power * cap * x ** (power - 1) + sign * derivative)
        print(f"{name}: both exact derivative signs verified on [0,2].")


if __name__ == "__main__":
    verify()
