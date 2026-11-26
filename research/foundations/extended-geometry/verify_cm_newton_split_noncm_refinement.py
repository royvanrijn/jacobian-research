#!/usr/bin/env python3
"""Replay the selected and failing split lines at c = -15."""

from fractions import Fraction
from math import comb

from verify_spherical_critical_face_p_power_moments import power


def middle_polynomial_coefficient(exponent, height, c, prime):
    if not 0 <= height <= 2*exponent:
        return 0
    coefficient = sum(comb(exponent, j)*comb(exponent, height-j)*c**(exponent-height+j)
                      for j in range(max(0, height-exponent), min(exponent, height)+1))
    return ((-1)**(2*exponent-height) * coefficient) % prime


def refinement_dimensions(i):
    return tuple(((8*t-5*i) % 40 + 10*i + (-8*t-5*i) % 40)//40 - 1
                 for t in range(5))


def original_matrix(i, c, prime):
    dimension = 1 + 2*((5*i)//8)
    z_exponent = prime * (3*i % 8)//8
    exponent = prime * i//8

    def coefficient(height):
        index, remainder = divmod(height-z_exponent, 5)
        return 0 if remainder else middle_polynomial_coefficient(exponent, index, c, prime)

    return tuple(tuple(coefficient(prime*source-target) for source in range(1, dimension+1))
                 for target in range(1, dimension+1))


def check_split_prime(prime, expected_quotient, selected, failing_middle):
    c = -15
    order = (prime-1)//8
    assert prime % 8 == 1 and prime % 5 in (2, 3)
    quotient = tuple(middle_polynomial_coefficient((5*i % 8)*order, (5*i % 8)*order,
                                                   c, prime) for i in (1, 3, 5, 7))
    assert quotient == expected_quotient
    assert original_matrix(1, c, prime) == ((selected,),)
    assert original_matrix(3, c, prime) == ((0, 0, 0), (0, failing_middle, 0), (0, 0, 0))
    assert quotient[2] == selected and quotient[3] == failing_middle
    print(f"p={prime}: quotient {quotient}; selected {selected}; character 3 has rank 1 of 3")


def check_sphere_moment():
    c, prime = -15, 17
    # Literal coefficient of y^5 t^3 after multiplying by -32 c.
    phase = {(-5, 3+2*j): -c*(-1)**j*comb(5, j) for j in range(6)}
    phase.update({(0, 8): -(1+c), (5, 3): 1})
    moment = sum(Fraction(value, height+1)
                 for (angle, height), value in power(phase, 2).items()
                 if angle == 0 and height % 2 == 0)
    normalized = prime * moment
    residue = normalized.numerator * pow(normalized.denominator, -1, prime) % prime
    assert residue == 13
    print("p=17: normalized order-two sphere moment is 13")


if __name__ == "__main__":
    assert (1-2**4) % 9 == 3
    assert refinement_dimensions(1) == (1, 0, 0, 0, 0)
    assert refinement_dimensions(3) == (1, 1, 0, 0, 1)
    check_split_prime(17, (13, 13, 13, 13), 13, 13)
    check_split_prime(73, (57, 20, 20, 57), 20, 57)
    check_sphere_moment()
