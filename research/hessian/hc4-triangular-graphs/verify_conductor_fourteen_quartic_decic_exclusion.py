"""Check exact Laurent and collision coefficients for sequence (10,4,7)."""

from functools import cache
from itertools import product

from sympy import Poly, Rational, binomial, diff, expand, factorial, resultant, symbols, together


def coefficient_factory(a, b, c):
    @cache
    def series_coefficient(power, degree):
        if degree < 0:
            return 0
        value = 0
        for exponents in product(*(range(degree//weight + 1) for weight in (2, 3, 4))):
            if sum(weight*index for weight, index in zip((2, 3, 4), exponents)) != degree:
                continue
            total = sum(exponents)
            term = binomial(power, total)*factorial(total)
            for parameter, index in zip((a, b, c), exponents):
                term *= parameter**index/factorial(index)
            value += term
        return value
    return series_coefficient


def check_laurent_and_collision_coefficients():
    z, a, b, delta, s, product_of_endpoints = symbols("z a b Delta s pi")
    c = (a**2-delta)/4
    quartic = z**4+a*z**2+b*z+c
    alpha = -5*(3*delta**2-8*a*b**2)/64
    beta = -5*(delta**3+8*a*b**2*delta-4*b**4)/256
    coefficient = coefficient_factory(a, b, c)

    def square_root_coefficient(k):
        return (coefficient(Rational(5, 2), k)
                + alpha*coefficient(Rational(1, 2), k-8)/2
                + beta*coefficient(Rational(-1, 2), k-12)/2)

    assert expand(square_root_coefficient(11)) == 0
    assert expand(square_root_coefficient(12)) == 0
    assert expand(square_root_coefficient(13)-5*b**3*delta/128) == 0
    decic = Poly(sum(square_root_coefficient(k)*z**(10-k) for k in range(11)), z)
    root = Poly(decic.as_expr()**2-quartic**5-alpha*quartic**3-beta*quartic**2, z)
    assert root.degree() == 7
    assert expand(root.nth(7)+5*b**3*delta/64) == 0

    differences = [0, 1]
    for _ in range(2, 11):
        differences.append(expand(s*differences[-1]-product_of_endpoints*differences[-2]))
    quartic_difference = sum(Poly(quartic, z).nth(i)*differences[i] for i in range(1, 5))
    decic_difference = sum(decic.nth(i)*differences[i] for i in range(1, 11))
    product_value = (s**3+a*s+b)/(2*s)
    assert expand(quartic_difference.subs(product_of_endpoints, product_value)) == 0
    sum_polynomial = (2*s**7 + 8*a*s**5 - 2*b*s**4 + (5*delta+8*a**2)*s**3
                      - 4*a*b*s**2 + (10*delta*a-2*b**2)*s - 5*delta*b)
    difference = 32*decic_difference.subs(product_of_endpoints, product_value) - s**2*sum_polynomial
    assert expand(together(difference).as_numer_denom()[0]) == 0

    shape = {a: Rational(3, 8), b: 1, delta: 1}
    quartic_shape = quartic.subs(shape)
    decic_shape = decic.as_expr().subs(shape)
    assert resultant(diff(quartic_shape, z), diff(decic_shape, z), z) == 208000


if __name__ == "__main__":
    check_laurent_and_collision_coefficients()
    print("Verified the quartic–decic Laurent coefficients, collision equation and leading-shape resultant.")
