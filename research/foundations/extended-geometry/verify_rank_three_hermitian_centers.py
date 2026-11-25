"""Exact interlacing-polytope proof of the Hermitian-center moment formulas."""

from functools import cache
from math import comb
import sympy as sp

x, y, z, lower, upper, entry = sp.symbols("x y z lower upper entry")
c1, c2, c3, g, h = sp.symbols("c1 c2 c3 g h")
v, w, tau = sp.symbols("v w tau")
variables = (x, y, z, lower, upper, entry)
total = c1 + c2 + c3
traceless_square = sum((c - total / 3)**2 for c in (c1, c2, c3))
traceless_cube = sum((c - total / 3)**3 for c in (c1, c2, c3))
diagonal = sp.Poly(
    (c1 - c2) * entry + (c2 - c3) * (lower + upper)
    + c3 * (x + y + z) - total / 2,
    *variables,
)
first_norm = sp.Poly(x * y * z, *variables)
second_norm = sp.Poly((1 - x) * (1 - y) * (1 - z), *variables)


def simplex(i, j, k):
    return sp.Rational(1, (i + 1) * (i + j + 2) * (i + j + k + 3))


def rectangle(i, j, k, a, b):
    return (
        simplex(i, j + a + 1, k + b + 1)
        - simplex(i + a + 1, j, k + b + 1)
        - simplex(i, j + a + b + 2, k)
        + simplex(i + a + 1, j + b + 1, k)
    ) / ((a + 1) * (b + 1))


@cache
def interlacing(i, j, k, a, b, d):
    return 360 * (
        rectangle(i, j, k, a, b + d + 1)
        - rectangle(i, j, k, a + d + 1, b)
    ) / (d + 1)


@cache
def mixed(center_power, first_power, second_power):
    polynomial = diagonal**center_power * first_norm**first_power * second_norm**second_power
    return sp.expand(sum(coefficient * interlacing(*exponents)
                         for exponents, coefficient in polynomial.terms()))


def moment(order):
    total = 0
    for first in range(order // 2 + 1):
        for second in range(order // 2 - first + 1):
            factor = (
                (-(g + h))**first * (-(g - h))**second
                * comb(order, 2 * first) * comb(order - 2 * first, 2 * second)
                * sp.Rational(comb(2 * first, first) * comb(2 * second, second),
                              4**(first + second))
            )
            total += factor * mixed(order - 2 * first - 2 * second, first, second)
    return sp.expand(total)


def verify():
    assert interlacing(0, 0, 0, 0, 0, 0) == 1
    assert moment(1) == 0
    substitution = {tau: total, v: traceless_square, w: traceless_cube}
    second = (-6*g + tau**2 + 3*v) / 84
    third = -3*h*tau / 140
    fourth = (
        567*g**2 - 396*g*tau**2 - 540*g*v + 243*h**2 + 19*tau**4
        + 90*tau**2*v - 72*tau*w + 162*v**2
    ) / 45360
    for order, expected in ((2, second), (3, third), (4, fourth)):
        assert sp.expand(moment(order) - expected.subs(substitution)) == 0

    value_g = (tau**2 + 3*v) / 6
    sixth_balanced = -(
        2746*tau**6 - 4869*tau**4*v - 10836*tau**3*w - 51192*tau**2*v**2
        - 2268*tau*v*w + 2187*v**3 - 18144*w**2
    ) / 150885504
    sixth_reduced = moment(6).subs(h, 0).subs(g, value_g.subs(substitution))
    assert sp.expand(sixth_reduced - sixth_balanced.subs(substitution)) == 0

    assert sp.expand(fourth.subs(g, value_g).subs(tau, 0) - (36*h**2 + 5*v**2)/6720) == 0
    value_w = (135*v**2 - 414*v - 125)/288
    assert sp.expand(fourth.subs(g, value_g).subs({tau: 1, h: 0, w: value_w})) == 0
    quartic = 14175*v**4 - 90936*v**3 + 295542*v**2 + 38928*v - 14333
    assert sp.expand(sixth_balanced.subs({tau: 1, w: value_w}) - quartic/536481792) == 0
    print("Interlacing moments verify the six-moment Hermitian-center exclusion.")


if __name__ == "__main__":
    verify()
