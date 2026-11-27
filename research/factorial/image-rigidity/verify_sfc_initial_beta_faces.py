#!/usr/bin/env python3
"""Exact controls for initial fractional faces and the Beta Gram obstruction.

The general face theorem is a written square-integral argument. These finite
controls check its normalization and resultant specialization. The fractional
gap witnesses are certified using rational intervals, without floating point.
"""

from fractions import Fraction
from math import comb, factorial, gcd, prod

import sympy as sp


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def rising(start, length):
    return prod(start + j for j in range(length))


def face_coefficients(a, w, h, other_columns):
    """Coefficients of U_1 and the symmetric coefficient matrix of U_2."""
    baseline = Fraction(a) - Fraction(h, 3)
    first = [
        rising(baseline + 1, k * w)
        * prod(factorial(row[k]) for row in other_columns)
        for k in (0, 1)
    ]
    second = [
        rising(2 * baseline + 1, k * w)
        * prod(factorial((2 - k) * row[0] + k * row[1]) for row in other_columns)
        for k in (0, 1, 2)
    ]
    return first, second


def verify_faces():
    checked = 0
    profiles = [(), ((0, 1),), ((2, 0),), ((0, 2), (3, 0))]
    for a in range(1, 4):
        for w in range(4):
            for h in range(1, a + 1):
                if gcd(h, 3) != 1:
                    continue
                for rows in profiles:
                    if w == 0 and all(row[0] == row[1] for row in rows):
                        continue
                    (p, r), (u, v, z) = face_coefficients(a, w, h, rows)
                    result = p * p * z - 2 * p * r * v + r * r * u
                    require(result > 0, "The initial face resultant is not positive.")
                    integral_result = result * 3 ** (4 * w)
                    require(integral_result.denominator == 1, "Face clearing failed.")
                    require(u * z - v * v > 0, "The face Gram determinant failed.")
                    checked += 1

    # A repeated column must fail strictness; distinctness is indispensable.
    (p, r), (u, v, z) = face_coefficients(1, 0, 1, ((2, 2),))
    require(p * p * z - 2 * p * r * v + r * r * u == 0,
            "The repeated-column negative control failed.")
    return checked


def full_form(order, a, w, rows, baseline, variables):
    result = 0
    for i in range(order + 1):
        for j in range(order - i + 1):
            k = order - i - j
            weight = factorial(order) // (factorial(i) * factorial(j) * factorial(k))
            offset = a * j + (a + w) * k
            weight *= rising(order * baseline + 1, offset)
            weight *= prod(factorial(i * row[0] + j * row[1] + k * row[2])
                           for row in rows)
            result += weight * variables[0] ** i * variables[1] ** j * variables[2] ** k
    return sp.expand(result)


def verify_specializations():
    x, y, z = sp.symbols("X Y Z")
    count = 0
    for a, w, h, rows in (
        (1, 0, 1, ((2, 0, 1),)),
        (1, 1, 1, ((2, 1, 0),)),
        (2, 1, 2, ((1, 3, 0),)),
        (2, 2, 1, ((0, 1, 2), (1, 0, 1))),
    ):
        baseline = -sp.Rational(h, 3)
        f1, f2, f3 = [full_form(q, a, w, rows, baseline, (x, y, z))
                      for q in (1, 2, 3)]
        anchor = prod(factorial(3 * row[0]) for row in rows)
        require(f3 == anchor * x ** 3, "The terminal mark did not isolate the anchor.")
        coefficient = f1.coeff(x)
        substitution = -(f1 - coefficient * x) / coefficient
        q2 = sp.expand(f2.subs(x, substitution).subs(z, 1))
        q3 = sp.expand(f3.subs(x, substitution).subs(z, 1))
        full_result = coefficient ** 6 * sp.resultant(q2, q3, y)

        face_rows = tuple((row[1], row[2]) for row in rows)
        (p, r), (u, v, zz) = face_coefficients(a, w, h, face_rows)
        face_result = p * p * zz - 2 * p * r * v + r * r * u
        rh = 3 ** (4 * w) * face_result
        b1 = prod(3 * nu - h for nu in range(1, a + 1))
        b2 = prod(3 * nu - 2 * h for nu in range(1, 2 * a + 1))
        integer = anchor ** 2 * rh ** 3 * b1 ** 6 * b2 ** 3
        require(sp.Rational(integer) == 3 ** (12 * (a + w)) * full_result,
                "The marked full resultant and divisor integer disagree.")
        count += 1
    return count


def interval_add(left, right):
    return left[0] + right[0], left[1] + right[1]


def interval_multiply(left, right):
    values = [a * b for a in left for b in right]
    return min(values), max(values)


def arctan_interval(inverse, terms=18):
    partial = sum((Fraction((-1) ** k, (2 * k + 1) * inverse ** (2 * k + 1))
                   for k in range(terms)), Fraction(0))
    correction = Fraction((-1) ** terms,
                          (2 * terms + 1) * inverse ** (2 * terms + 1))
    return min(partial, partial + correction), max(partial, partial + correction)


def pi_interval():
    # Machin's identity: pi = 16 atan(1/5) - 4 atan(1/239).
    return interval_add(interval_multiply((16, 16), arctan_interval(5)),
                        interval_multiply((-4, -4), arctan_interval(239)))


def polynomial_interval(expression, variable, bounds):
    result = (Fraction(0), Fraction(0))
    for coefficient in sp.Poly(expression, variable).all_coeffs():
        coefficient = Fraction(coefficient)
        result = interval_add(interval_multiply(result, bounds),
                              (coefficient, coefficient))
    return result


def rational_interval(expression):
    numerator, denominator = sp.fraction(sp.cancel(expression))
    bounds = pi_interval()
    nn = polynomial_interval(numerator, sp.pi, bounds)
    dd = polynomial_interval(denominator, sp.pi, bounds)
    require(dd[0] > 0, "The rational-function denominator is not certified positive.")
    return interval_multiply(nn, (1 / dd[1], 1 / dd[0]))


def angular_coefficients(left, right, middle, span):
    offsets = (0, middle, span)

    def beta(a, b):
        return sp.simplify(sp.gamma(a) * sp.gamma(b) / sp.gamma(a + b))

    norms = [beta(left + e + 1, right + span - e + 1) for e in offsets]

    def moment(i, j, k):
        order = i + j + k
        offset = j * middle + k * span
        return sp.simplify(beta(order * left + offset + 1,
                                order * (right + span) - offset + 1)
                           / (norms[0] ** i * norms[1] ** j * norms[2] ** k))

    def contrast(order, y_degree):
        return sp.simplify(sum(
            comb(order, y_degree) * comb(order - y_degree, i) * comb(y_degree, j)
            * (-1) ** (order - i - j) * moment(i, j, order - i - j)
            for i in range(order - y_degree + 1) for j in range(y_degree + 1)))

    aa, bb, cc = [contrast(2, k) for k in range(3)]
    dd, ff, gg = [contrast(3, k) for k in (0, 2, 3)]
    gram = sp.Matrix([[moment(2, 0, 0), moment(1, 1, 0), moment(1, 0, 1)],
                      [moment(1, 1, 0), moment(0, 2, 0), moment(0, 1, 1)],
                      [moment(1, 0, 1), moment(0, 1, 1), moment(0, 0, 2)]])
    return sp.factor(cc * cc * dd + aa * bb * gg - aa * cc * ff), gram


def verify_fractional_obstruction():
    small, small_gram = angular_coefficients(1, 0, sp.Rational(1, 2), 1)
    expected = 9 * (143613 * sp.pi ** 4 - 996096 * sp.pi ** 2 - 4194304) / (
        448000 * sp.pi ** 4)
    require(sp.cancel(small - expected) == 0, "The small fractional witness changed.")
    require(rational_interval(small)[1] < 0, "The small signed remainder is not negative.")
    require(rational_interval(small_gram.det())[0] > 0,
            "The small Beta Gram determinant is not positive.")

    value, gram = angular_coefficients(19, 0, sp.Rational(3, 2), 2)
    low, high = rational_interval(value)
    require(Fraction(-4, 3) < low and high < Fraction(-13, 10),
            "The span-two fractional witness escaped its rational enclosure.")
    require(rational_interval(gram.det())[0] > 0,
            "The span-two Beta Gram determinant is not positive.")
    return 2


def main():
    faces = verify_faces()
    specializations = verify_specializations()
    obstructions = verify_fractional_obstruction()
    print(f"Passed {faces} exact fractional-face controls and a repeated-column control.")
    print(f"Passed {specializations} full resultant specializations, including varying unmarked rows.")
    print(f"Certified {obstructions} negative fractional-gap remainders with positive Beta Gram determinants.")
    print("At (L,H,a,d)=(19,0,3/2,2), -4/3 < S_angular < -13/10.")


if __name__ == "__main__":
    main()
