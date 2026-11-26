#!/usr/bin/env python3
"""Replay the second pole and the residue-orbit valuation."""

from fractions import Fraction
from math import comb, gcd

from verify_spherical_critical_face_p_power_moments import multiply, power


def factorial_valuation(number, prime):
    valuation = 0
    while number:
        number //= prime
        valuation += number
    return valuation


def literal_phase(d, b, a):
    def factor(shift, constant):
        return {(d, 0): 1, (0, 1): -1, (1, 0): -shift, (0, 0): constant}

    phase = {(-b, 0): 1}
    for shift, constant, multiplicity in ((1, -1, a), (1, 1, b-a),
                                           (2, -1, b-a), (2, 1, a)):
        phase = multiply(phase, power(factor(shift, constant), multiplicity))
    return phase


def check_literal_moment(B):
    d, a, prime, r = 2, 1, 5, 2
    b, n = d * B, (2*d-1) * B
    q = prime**r
    assert (q-1) % (d*n) == 0
    order = (q-1) // n
    phase = literal_phase(d, b, a)
    peak = max(j + Fraction(i, d) for i, j in phase)
    assert peak == n
    assert all(j >= 0 and n-2*b <= j + Fraction(i, d) for i, j in phase)
    moment = sum(Fraction(coefficient, height+1)
                 for (angle, height), coefficient in power(phase, order).items()
                 if angle == 0 and height % 2 == 0)
    beta = Fraction((-1)**(a*order) * 2**(b*order),
                    (b*order+1)*comb(b*order, a*order))
    assert moment == beta
    normalized = prime * moment
    assert normalized.numerator * pow(normalized.denominator, -1, prime) % prime == 2
    assert (q * moment).numerator % prime == 0
    print(f"B={B}, p={prime}, r={r}: order {order}, second residue 2, top residue 0")


def residue_orbit(multiplier, modulus):
    orbit = [1]
    current = multiplier % modulus
    while current != 1:
        orbit.append(current)
        current = current * multiplier % modulus
    return orbit


def check_density_and_depth():
    for height in (1, 2, 3, 4):
        n = 3**height
        b = 2*n//3
        units = [unit for unit in range(1, n) if gcd(unit, n) == 1]
        qualifying = 0
        for unit in units:
            orbit = residue_orbit(unit, n)
            if n-1 in orbit:
                assert unit % 3 == 2
                carries = sum((entry % n) + ((b-1)*entry % n) >= n for entry in orbit)
                assert carries * 2 == len(orbit)
                qualifying += 1
        assert qualifying * 2 == len(units)

    n, b, a, prime = 9, 6, 1, 5
    r = len(residue_orbit(prime, 2*n))
    order = (prime**r - 1) // n
    valuation = (factorial_valuation(b*order, prime)
                 - factorial_valuation(a*order, prime)
                 - factorial_valuation((b-a)*order, prime))
    assert valuation == r//2 == 3
    print("Density classes checked at 3, 9, 27, 81; B=3, p=5: depth 3 at r=6")


if __name__ == "__main__":
    check_literal_moment(1)
    check_literal_moment(2)
    check_density_and_depth()
