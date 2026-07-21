"""Exact checks for the universal factorization-slice note.

The finite-field enumerations use homogeneous binary forms: two coefficient
vectors have a common root at infinity when both top coefficients vanish.
"""

from itertools import product

import sympy as sp


def trim(poly):
    values = list(poly)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return values


def remainder(dividend, divisor, prime):
    dividend = trim(dividend)
    divisor = trim(divisor)
    inverse = pow(divisor[-1], -1, prime)
    while any(dividend) and len(dividend) >= len(divisor):
        shift = len(dividend) - len(divisor)
        coefficient = dividend[-1] * inverse % prime
        for index, value in enumerate(divisor):
            dividend[index + shift] = (
                dividend[index + shift] - coefficient * value
            ) % prime
        dividend = trim(dividend)
    return dividend


def homogeneous_coprime(first, second, prime):
    # The affine gcd misses the common root at infinity.
    if first[-1] == 0 and second[-1] == 0:
        return False
    left = list(first)
    right = list(second)
    while any(right):
        left, right = right, remainder(left, right, prime)
    return len(trim(left)) == 1


def projective_forms(degree, prime):
    for coefficients in product(range(prime), repeat=degree + 1):
        if not any(coefficients):
            continue
        first_nonzero = next(index for index, value in enumerate(coefficients) if value)
        if coefficients[first_nonzero] == 1:
            yield coefficients


def product_coefficients(first, second, prime):
    output = [0] * (len(first) + len(second) - 1)
    for i, left in enumerate(first):
        for j, right in enumerate(second):
            output[i + j] = (output[i + j] + left * right) % prime
    return output


def factorization_count(a, b, functional, prime):
    count = 0
    for first in projective_forms(a, prime):
        for second in projective_forms(b, prime):
            if not homogeneous_coprime(first, second, prime):
                continue
            coefficients = product_coefficients(first, second, prime)
            value = sum(
                coefficient * entry
                for coefficient, entry in zip(functional, coefficients)
            ) % prime
            if value:
                count += 1
    return count


# Restriction of a coefficient functional to the small diagonal.
alpha, beta = sp.symbols("alpha beta")
cube = [alpha**3, 3 * alpha**2 * beta, 3 * alpha * beta**2, beta**3]
assert sp.factor(cube[1]) == 3 * alpha**2 * beta
assert sp.factor(cube[1] - cube[2]) == 3 * alpha * beta * (alpha - beta)
assert cube[0] == alpha**3
print("PASS: the cubic hyperplanes have contact types (2,1), (1,1,1), and (3)")


# Resultant scaling has weight b-a under (Q,R) -> (lambda Q,lambda^-1 R).
lam = sp.symbols("lambda", nonzero=True)
q0, q1, r0, r1, r2 = sp.symbols("q0 q1 r0 r1 r2")
linear = q0 * sp.Symbol("W") + q1
quadratic = r0 * sp.Symbol("W") ** 2 + r1 * sp.Symbol("W") + r2
resultant = sp.resultant(linear, quadratic, sp.Symbol("W"))
scaled = sp.resultant(lam * linear, quadratic / lam, sp.Symbol("W"))
assert sp.simplify(scaled - lam * resultant) == 0
print("PASS: the consecutive-degree resultant supplies a unique torus gauge")


# Cubic contact-type counts. Coefficients are in increasing affine degree;
# reversing U,V only interchanges a contact partition with itself.
for prime in (5, 7, 11):
    transverse = factorization_count(1, 2, (0, 1, -1, 0), prime)
    tangent = factorization_count(1, 2, (0, 1, 0, 0), prime)
    osculating = factorization_count(1, 2, (1, 0, 0, 0), prime)
    assert transverse == prime**3 - prime
    assert tangent == prime**3
    assert osculating == prime**3 - prime**2
print("PASS: cubic contact slices have counts q^3-q, q^3, and q^3-q^2")


# The (2,3) middle-coefficient slice. Both middle coefficients give the same
# count after interchanging U and V.
for prime in (2, 3, 5, 7, 11):
    functional = (0, 0, 1, 0, 0, 0)
    count = factorization_count(2, 3, functional, prime)
    assert count == prime**5 - prime**3 + prime**2
print("PASS: the (2,3) candidate has count q^5-q^3+q^2, not q^5")


# Symbolic audit of the three Moebius/rank contributions.
q = sp.symbols("q")
degree_zero = (q**2 + q + 1) * q**3
degree_one = -(q**3 + q * (q + 1) * q**2)
degree_two = q**2
assert sp.expand(degree_zero + degree_one + degree_two) == q**5 - q**3 + q**2
print("PASS: homogeneous Moebius inversion gives the announced count polynomial")


# Weighted targets are the pullback of the universal polynomial along sigma_H.
A, B, C, c, W = sp.symbols("A B C c W")
H = W**5 - W**2
E = H - B * C * W + c * A * C**2
assert sp.diff(E, A) == c * C**2
assert sp.diff(E, B) == -C * W
assert sp.diff(E, C) == -B * W + 2 * c * A * C
print("PASS: weighted inverse polynomials form the stated three-parameter pullback")
