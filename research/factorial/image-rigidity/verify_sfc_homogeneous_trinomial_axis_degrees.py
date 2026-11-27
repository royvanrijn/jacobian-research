#!/usr/bin/env python3
"""Decide the initial factorial window for one homogeneous axis profile."""

import argparse
from fractions import Fraction
from math import comb, factorial, prod

import sympy as sp
from sympy.polys.matrices import DomainMatrix

D = sp.symbols("D")
SYLVESTER_ENTRY_BOUND = 500
RESULTANT_ERROR_BOUND = 600 * SYLVESTER_ENTRY_BOUND**5


def moment_rows(degree, offsets, opposite_third, *, decoupled):
    rows = []
    for order in (1, 2, 3):
        terms = {}
        for a in range(order + 1):
            for b in range(order - a + 1):
                c = order - a - b
                powers = (a, b, c)
                mixed = opposite_third and 0 < c < order
                if mixed and decoupled:
                    continue
                if mixed:
                    x_degree = a * offsets[0] + b * offsets[1] + c * (degree - offsets[2])
                    ratio = Fraction(1, comb(order * degree, x_degree))
                else:
                    low_degree = sum(k * exponent for k, exponent in zip(powers, offsets))
                    denominator = sp.prod(order * degree - i for i in range(low_degree))
                    ratio = sp.S(factorial(low_degree)) / denominator
                multinomial = factorial(order) // prod(factorial(k) for k in powers)
                terms[powers] = multinomial * ratio
        rows.append(terms)
    return rows


def restrict_to_first_zero(rows):
    first, second, third = rows
    x, y, z = (first[power] for power in ((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    forms = []
    for order, row in ((2, second), (3, third)):
        coefficients = [0] * (order + 1)
        for (a, b, c), coefficient in row.items():
            # X=-y*U-z*V, Y=x*U, Z=x*V.
            for i in range(a + 1):
                coefficients[b + i] += (
                    coefficient * (-1)**a * comb(a, i) * y**i * z**(a - i) * x**(b + c)
                )
        forms.append(coefficients)
    return forms


def binary_resultant(forms):
    quadratic, cubic = (list(reversed(form)) for form in forms)
    rows = [[0] * i + quadratic + [0] * (2 - i) for i in range(3)]
    rows.extend([0] * i + cubic + [0] * (1 - i) for i in range(2))
    matrix = DomainMatrix.from_Matrix(sp.Matrix(rows)).to_field()
    return matrix.domain.to_sympy(matrix.det())


def degree_cutoff(offsets, opposite_third):
    decoupled = moment_rows(D, offsets, opposite_third, decoupled=True)
    resultant = sp.cancel(binary_resultant(restrict_to_first_zero(decoupled)))
    if resultant == 0:
        raise ArithmeticError("The decoupled resultant vanished identically.")
    numerator, denominator = (sp.Poly(p, D) for p in sp.fraction(resultant))
    leading = abs(numerator.LC())
    tail = sum(abs(coef) for coef in numerator.all_coeffs()[1:])
    denominator_bound = sum(abs(coef) for coef in denominator.all_coeffs())
    lower_bound = Fraction(leading / (2 * denominator_bound))
    power = max(0, denominator.degree() - numerator.degree())
    degree = max(6 * max(offsets), 1 + int(2 * tail / leading), 4 * power)
    if opposite_third:
        # Beyond 4*power the error decreases relative to the lower bound
        # separately on the even and odd degrees.
        while any(
            RESULTANT_ERROR_BOUND * n**power >= lower_bound * 2**(n // 2)
            for n in (degree, degree + 1)
        ):
            degree += 1
    return degree


def failing_degrees(offsets, opposite_third, cutoff):
    failures = []
    for degree in range(max(offsets), cutoff):
        first_exponents = (
            (offsets[0], offsets[1], degree - offsets[2]) if opposite_third else offsets
        )
        if len(set(first_exponents)) != 3:
            continue
        rows = moment_rows(degree, offsets, opposite_third, decoupled=False)
        if binary_resultant(restrict_to_first_zero(rows)) == 0:
            failures.append(degree)
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("offsets", type=int, nargs=3)
    parser.add_argument("--opposite-third", action="store_true")
    args = parser.parse_args()
    offsets = tuple(args.offsets)
    opposite_third = args.opposite_third
    if min(offsets) < 0:
        parser.error("Offsets must be nonnegative.")
    same_side = offsets[:2] if opposite_third else offsets
    if len(same_side) != len(set(same_side)):
        parser.error("Offsets on the same side must be distinct.")

    cutoff = degree_cutoff(offsets, opposite_third)
    failures = failing_degrees(offsets, opposite_third, cutoff)
    print(f"Degrees D >= {cutoff} are excluded.")
    print(f"Failing degrees below the cutoff: {failures}")


if __name__ == "__main__":
    main()
