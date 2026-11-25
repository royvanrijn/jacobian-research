"""Exact finite moment exclusion in the rank-three signed spectral model."""

from functools import cache
from math import comb

import sympy as sp


x, y, z = sp.symbols("x y z")
a, b, c, g, h = sp.symbols("a b c g h")
r, inverse = sp.symbols("r inverse")

spectral_sum = sp.Poly(x + y + z - sp.Rational(3, 2), x, y, z)
first_norm = sp.Poly(x * y * z, x, y, z)
second_norm = sp.Poly((1 - x) * (1 - y) * (1 - z), x, y, z)
density = sp.Poly(180 * (y - x) * (z - x) * (z - y), x, y, z)


@cache
def spectral_moment(sum_power, first_power, second_power):
    """Integrate a monomial of the three spectral variables."""
    polynomial = (
        density
        * spectral_sum**sum_power
        * first_norm**first_power
        * second_norm**second_power
    )
    return sum(
        coefficient / ((i + 1) * (i + j + 2) * (i + j + k + 3))
        for (i, j, k), coefficient in polynomial.terms()
    )


@cache
def centered_moment(sum_power, difference_power, total_power, first_power, second_power):
    if difference_power:
        arguments = (sum_power, difference_power - 1, total_power)
        return centered_moment(*arguments, first_power + 1, second_power) - centered_moment(
            *arguments, first_power, second_power + 1
        )
    if total_power:
        arguments = (sum_power, 0, total_power - 1)
        return (
            centered_moment(*arguments, first_power + 1, second_power)
            + centered_moment(*arguments, first_power, second_power + 1)
            - centered_moment(*arguments, first_power, second_power) / 7
        )
    return spectral_moment(sum_power, first_power, second_power)


def circle_moment(half_power):
    return sp.Rational(comb(2 * half_power, half_power), 4**half_power)


def phase_moment(order):
    answer = 0
    for first in range(order // 2 + 1):
        for second in range(order // 2 - first + 1):
            left = order - 2 * first - 2 * second
            angular_weight = (
                (-(g + h)) ** first
                * (-(g - h)) ** second
                * comb(order, 2 * first)
                * comb(order - 2 * first, 2 * second)
                * circle_moment(first)
                * circle_moment(second)
            )
            for i in range(left + 1):
                for j in range(left - i + 1):
                    k = left - i - j
                    answer += (
                        angular_weight
                        * comb(left, i)
                        * comb(left - i, j)
                        * a**i
                        * b**j
                        * c**k
                        * centered_moment(i, j, k, first, second)
                    )
    return sp.expand(answer)


def reduced_numerator(moment, order, value_g, value_h):
    numerator = sp.fraction(sp.cancel(moment.subs(g, value_g).subs(h, value_h)))[0]
    if order % 2:
        numerator = sp.cancel(numerator / c)
    reduced = sp.expand(numerator).subs(c**2, r)
    assert not reduced.has(c)
    return reduced


def verify():
    assert spectral_moment(0, 0, 0) == 1
    assert spectral_moment(0, 1, 0) == spectral_moment(0, 0, 1) == sp.Rational(1, 14)
    moments = {order: phase_moment(order) for order in range(1, 7)}
    assert moments[1] == 0
    value_g = 3 * a**2 / 2 + 6 * a * b / 5 + 4 * b**2 / 15 + 4 * c**2 / 35
    assert sp.expand(moments[2].subs(g, value_g)) == 0
    assert sp.expand(moments[2].coeff(g) + sp.Rational(1, 14)) == 0
    m3 = moments[3].subs(g, value_g)
    value_h = (
        2 * c * (8085 * a**2 + 11319 * a * b + 3087 * b**2 - 257 * c**2)
        / (8085 * (9 * a + 4 * b))
    )
    assert sp.cancel(m3.subs(h, value_h)) == 0
    assert sp.expand(m3.coeff(h) + sp.Rational(1, 140) * (9 * a + 4 * b)) == 0

    singular = {b: -sp.Rational(9, 4) * a}
    assert sp.expand(m3.subs(singular) + c * (4112 * c**2 + 28077 * a**2) / 9055200) == 0
    fourth_singular = moments[4].subs(g, value_g).subs(singular).subs(c, 0)
    assert sp.expand(fourth_singular - (2402400 * h**2 + 123091 * a**4) / 448448000) == 0

    fourth, fifth, sixth = (
        reduced_numerator(moments[order], order, value_g, value_h)
        for order in (4, 5, 6)
    )
    for scale, variables, denominator_constraint in (
        ({a: 1}, (inverse, r, b), (inverse * (9 + 4 * b) - 1,)),
        ({a: 0, b: 1}, (r,), ()),
    ):
        equations = [eq.subs(scale) for eq in (fourth, fifth, sixth)]
        assert list(sp.groebner(equations + list(denominator_constraint), *variables)) == [1]
        equations = [eq.subs(scale).subs(r, 0) for eq in (fourth, sixth)]
        assert list(sp.groebner(equations + list(denominator_constraint), *variables)) == [1]
    print("First six moments exclude the radial-affine determinant phases.")


if __name__ == "__main__":
    verify()
