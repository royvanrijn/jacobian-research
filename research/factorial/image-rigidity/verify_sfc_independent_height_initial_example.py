#!/usr/bin/env python3
"""Verify the initial window for all common translates of {1, X, Y²}."""

from math import comb, prod

from sympy import Poly, ZZ, resultant, symbols


def stripped_moment(order, coefficients, height_x, height_y):
    anchor, linear, quadratic = coefficients
    terms = []
    for i in range(order + 1):
        for j in range(order - i + 1):
            terms.append(
                comb(order, i)
                * comb(order - i, j)
                * anchor ** (order - i - j)
                * linear**i
                * quadratic**j
                * prod(order * height_x + h for h in range(1, i + 1))
                * prod(order * height_y + h for h in range(1, 2 * j + 1))
            )
    return sum(terms)


def main():
    u, v, z = symbols("u v z")
    anchor = -(u + 1) - (v + 1) * (v + 2) * z
    coefficients = (anchor, 1, z)
    quadratic = Poly(stripped_moment(2, coefficients, u, v), z)
    cubic = Poly(stripped_moment(3, coefficients, u, v), z)
    assert quadratic.degree() == 2 and cubic.degree() == 3

    at_infinity = Poly(quadratic.LC(), u, v, domain=ZZ)
    assert at_infinity.eval({u: 0, v: 0}) == 20
    assert all(coefficient > 0 for coefficient in at_infinity.coeffs())

    eliminant = Poly(resultant(quadratic, cubic, z), u, v, domain=ZZ)
    assert eliminant.total_degree() == 18
    assert len(eliminant.terms()) == 91
    assert eliminant.eval({u: 0, v: 0}) == 382464
    assert all(coefficient > 0 for coefficient in eliminant.coeffs())
    print("Initial common-translation example verified by 91 positive resultant terms.")


if __name__ == "__main__":
    main()
