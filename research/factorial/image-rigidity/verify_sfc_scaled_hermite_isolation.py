"""Replay the finite Frobenius and valuation inputs in the Hermite isolation proof."""

from fractions import Fraction
from math import factorial
from warnings import catch_warnings, simplefilter

from sympy import Poly, symbols
from sympy.utilities.exceptions import SymPyDeprecationWarning

z = symbols("z")

# Each entry is (k, parity, prime, degrees of squarefree factors).
FROBENIUS = (
    (5, 0, 37, (2, 3)),
    (5, 1, 17, (2, 3)),
    (6, 0, 13, (1, 5)),
    (6, 0, 31, (1, 2, 3)),
    (6, 1, 31, (1, 5)),
    (6, 1, 19, (1, 2, 3)),
    (7, 0, 19, (3, 4)),
    (7, 1, 29, (2, 2, 3)),
    (8, 0, 19, (1, 2, 5)),
    (8, 1, 29, (3, 5)),
    (9, 0, 41, (1, 3, 5)),
    (9, 1, 41, (1, 3, 5)),
    (10, 0, 29, (3, 7)),
    (10, 1, 61, (3, 7)),
    (11, 0, 59, (4, 7)),
    (11, 1, 37, (1, 1, 2, 7)),
    (12, 0, 97, (1, 1, 3, 7)),
    (12, 1, 43, (1, 2, 2, 7)),
)

VALUATIONS = (
    (7, 7, (Fraction(-1, 3),)),
    (8, 7, (Fraction(1, 8),)),
    (9, 7, (Fraction(0), Fraction(1, 6))),
    (5, 5, (Fraction(-1, 4),)),
    (6, 5, (Fraction(1, 6),)),
    (9, 5, (Fraction(0), Fraction(1, 4))),
    (4, 3, (Fraction(1, 4),)),
    (9, 3, (Fraction(-2, 3), Fraction(-1, 2))),
)


def square_polynomial(order):
    """Return He_order(x) / x**(order % 2) as a polynomial in z=x**2."""
    half = order // 2
    terms = [
        (-1) ** h
        * (factorial(order) // (2**h * factorial(h) * factorial(order - 2 * h)))
        * z ** (half - h)
        for h in range(half + 1)
    ]
    return Poly(sum(terms), z)


def integer_valuation(number, prime):
    assert number != 0
    result = 0
    while number % prime == 0:
        number //= prime
        result += 1
    return result


def slope(left, right):
    return Fraction(right[1] - left[1], right[0] - left[0])


def node_valuations(order, prime):
    # Replacing z in the unscaled square polynomial by order*z
    # adds i*v_p(order) to its i-th coefficient valuation.
    polygon = []
    polynomial = square_polynomial(order)
    scale_valuation = integer_valuation(order, prime)
    for degree in range(polynomial.degree() + 1):
        coefficient = int(polynomial.nth(degree))
        if coefficient == 0:
            continue
        point = (
            degree,
            integer_valuation(coefficient, prime) + degree * scale_valuation,
        )
        while len(polygon) >= 2 and slope(polygon[-2], polygon[-1]) >= slope(
            polygon[-1], point
        ):
            polygon.pop()
        polygon.append(point)
    return tuple(sorted({-slope(a, b) / 2 for a, b in zip(polygon, polygon[1:])}))


def main():
    with catch_warnings():
        simplefilter("ignore", SymPyDeprecationWarning)
        for half, parity, prime, expected in FROBENIUS:
            reduced = square_polynomial(2 * half + parity).set_modulus(prime)
            factors = reduced.factor_list()[1]
            assert all(multiplicity == 1 for _, multiplicity in factors)
            assert tuple(sorted(factor.degree() for factor, _ in factors)) == expected

    for order, prime, expected in VALUATIONS:
        assert node_valuations(order, prime) == expected, (order, prime)

    # Nonzero squared roots here are simple and nonzero modulo prime.
    # Their signs lift unramified as well, since prime is odd.
    for prime, last_order in ((7, 6), (5, 4), (3, 3)):
        for order in range(2, last_order + 1):
            if prime == 3 and order == 3:
                # This normalized signed polynomial is x**2 - 1.
                assert square_polynomial(3).as_expr() == z - 3
                continue
            assert order % prime
            reduced = square_polynomial(order).set_modulus(prime)
            assert reduced.TC() != 0
            assert reduced.gcd(reduced.diff()).degree() == 0

    ninth = square_polynomial(9)
    for prime, zero_multiplicity in ((7, 3), (5, 2)):
        remainder = ninth.set_modulus(prime).exquo(
            Poly(z**zero_multiplicity, z, modulus=prime)
        )
        assert remainder.TC() != 0
        assert remainder.gcd(remainder.diff()).degree() == 0

    print("18 Frobenius certificates, 8 Newton valuations, and all unit-lift checks passed.")


if __name__ == "__main__":
    main()
