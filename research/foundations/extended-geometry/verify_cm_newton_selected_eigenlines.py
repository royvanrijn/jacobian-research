"""Verify the explicit nonsplit transition and the supporting-face residues."""

from collections import defaultdict
from fractions import Fraction
from math import comb, gcd


ROOTS = ((0, 1), (1, 1), (2, 2), (3, 2), (4, 2))


def root_polynomial(roots, modulus=None):
    coefficients = [1]
    for root, multiplicity in roots:
        for _ in range(multiplicity):
            updated = [0] * (len(coefficients) + 1)
            for index, value in enumerate(coefficients):
                updated[index] -= root * value
                updated[index + 1] += value
            coefficients = updated if modulus is None else [v % modulus for v in updated]
    return coefficients


def coherent_dimension(order, character, branches):
    residues = sum(character * exponent % order for exponent in branches)
    assert residues % order == 0
    return residues // order - 1


def frobenius_matrix(order, prime, roots, infinity, character):
    target = prime * character % order
    branches = tuple(exponent for _, exponent in roots) + (infinity,)
    source_dimension = coherent_dimension(order, character, branches)
    target_dimension = coherent_dimension(order, target, branches)
    factors = ((root, prime * (character * exponent % order) // order)
               for root, exponent in roots)
    coefficient = root_polynomial(factors, prime)

    def entry(source_index, target_index):
        index = prime * source_index - target_index
        return coefficient[index] if 0 <= index < len(coefficient) else 0

    return tuple(tuple(entry(a, b) for a in range(1, source_dimension + 1))
                 for b in range(1, target_dimension + 1))


def check_order_seven():
    assert pow(5, 6, 7) == 1
    assert frobenius_matrix(7, 5, ROOTS, 6, 1) == ((0,), (1,))
    assert frobenius_matrix(7, 5, ROOTS, 6, 5) == ((1, 0),)
    assert frobenius_matrix(7, 5, ROOTS, 6, 3) == ((3, 2, 3),)
    print("Order seven: return zero, selected transition 3 modulo 5.")


def laurent_lift(order):
    coefficients = defaultdict(Fraction)
    for index, value in enumerate(root_polynomial(ROOTS[1:])):
        angle = index - 6
        if angle >= 0:
            coefficients[angle, order - 6 * angle] += value
        else:
            pole = -angle
            for k in range(pole + 1):
                coefficients[angle, order + 6 * k] += (
                    value * (-1) ** (pole + k) * comb(pole, k)
                )
    return coefficients


def second_moment(coefficients):
    return sum(first * second / (j + ell + 1)
               for (i, j), first in coefficients.items()
               for (k, ell), second in coefficients.items()
               if i + k == 0 and (j + ell) % 2 == 0)


def residue(value, prime):
    assert value.denominator % prime != 0
    return value.numerator * pow(value.denominator, -1, prime) % prime


def check_family():
    selected = root_polynomial((root, 2 * exponent) for root, exponent in ROOTS[1:])[12]
    assert selected % 5 == 3
    for iterate in (1, 2, 3, 4):
        order = (3 * 5 ** iterate - 1) // 2
        branches = (order - 6, 1, 2, 2, 2, order - 1)
        assert sum(branches) == 2 * order and gcd(order, *branches) == 1
        different_dimension = 3 * 5 ** (iterate - 1)
        assert different_dimension < order and 5 * different_dimension % order == 1
        assert coherent_dimension(order, different_dimension, branches) == 3
        assert coherent_dimension(order, 1, branches) == 1
        phase = laurent_lift(order)
        phase[0, 1] += Fraction(7, 2)
        assert all((j + 6 * i - order) % 3 == 0 for i, j in phase)
        moment = second_moment(phase)
        assert residue(3 * 5 ** iterate * moment, 5) == selected % 5
        print(f"Connected order {order}, degree {order + 30}, normalized second moment 3 modulo 5.")

    # This lower monomial violates the support congruence. Its squared
    # interval term contributes at the smaller denominator 5.
    phase = laurent_lift(7)
    phase[0, 2] += Fraction(2)
    assert residue(15 * second_moment(phase), 5) == 0
    print("Negative control: a lower term outside the congruence cancels the residue.")


def check_zero_return_progression():
    for iterate in (4, 10, 16):
        order = (3 * 5 ** iterate - 1) // 2
        character = (17 * order + 1) // 27
        assert character == pow(5, 3 * iterate, order)
        assert Fraction(3, 5) < Fraction(character, order) < Fraction(19, 30)
        roots = ((0, order - 6), *ROOTS[1:])
        assert frobenius_matrix(order, 5, roots, order - 1, character) == ((0,),)
        print(f"Connected order {order}: zero map in the character-one return orbit.")


def main():
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    check_order_seven()
    check_family()
    check_zero_return_progression()


if __name__ == "__main__":
    main()
