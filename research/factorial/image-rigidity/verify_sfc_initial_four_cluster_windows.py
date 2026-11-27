#!/usr/bin/env python3
"""Check the initial bivariate Gaussian windows through four columns."""

import argparse
from fractions import Fraction
from math import factorial
from subprocess import run

from sympy import Matrix, Rational


def counts(total, size):
    if size == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in counts(total - first, size - 1):
            yield (first,) + tail


def gaussian_terms(order, powers):
    for alpha in counts(order, len(powers)):
        degree = sum(k * n for k, n in zip(powers, alpha))
        if degree % 2:
            continue
        half = degree // 2
        denominator = (2 * order) ** half * factorial(half)
        multinomial = factorial(order)
        for n in alpha:
            denominator *= factorial(n)
        yield alpha, Fraction(multinomial * factorial(degree), denominator)


def triple_margin(orders):
    degree = sum(order - 1 for order in orders) + 1
    monomials = list(counts(degree, 3))
    index = {alpha: i for i, alpha in enumerate(monomials)}
    columns = []
    for order in orders:
        terms = tuple(gaussian_terms(order, (0, 1, 2)))
        for multiplier in counts(degree - order, 3):
            column = [0] * len(monomials)
            for alpha, value in terms:
                total = tuple(a + b for a, b in zip(alpha, multiplier))
                column[index[total]] = Rational(value.numerator, value.denominator)
            columns.append(column)
    matrix = Matrix(columns).T
    pure_powers = Matrix([
        [int(alpha[j] == degree) for j in range(3)] for alpha in monomials
    ])
    solution, parameters = matrix.gauss_jordan_solve(pure_powers)
    solution = solution.subs({parameter: 0 for parameter in parameters})
    if matrix * solution != pure_powers:
        raise RuntimeError(f"Invalid homogeneous ideal identities for {orders}.")
    coefficient_bound = max(sum(abs(x) for x in solution[:, j]) for j in range(3))
    return Fraction(1 / coefficient_bound)


def normalized_profile_forms(profile, prime):
    degrees, powers = zip(*profile)
    for order in (1, 2, 3, 4):
        monomials = []
        for alpha, value in gaussian_terms(order, powers):
            exponent = sum(d * n for d, n in zip(degrees, alpha))
            residue = value.numerator * pow(value.denominator, -1, prime)
            residue = residue * pow(order, exponent, prime) % prime
            if residue:
                monomial = '*'.join(f'x{i + 1}^{n}' for i, n in enumerate(alpha) if n)
                monomials.append(f'{residue}*{monomial}')
        yield '+'.join(monomials) or '0'


def check_profiles(profiles, prime, singular):
    program = [f'ring R={prime},(x1,x2,x3,x4),dp;']
    for profile in profiles:
        forms = ','.join(normalized_profile_forms(profile, prime))
        program.extend([f'ideal I={forms};', 'dim(std(I))==0;', 'kill I;'])
    program.append('quit;')
    result = run([singular, '-q'], input='\n'.join(program), text=True,
                 capture_output=True, check=True)
    results = result.stdout.splitlines()
    if results != ['1'] * len(profiles):
        failures = [profile for profile, status in zip(profiles, results) if status != "1"]
        raise RuntimeError(f"Got {len(results)} responses with failures {failures}. {result.stderr}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--singular", default="Singular")
    args = parser.parse_args()

    margins = (triple_margin((1, 2, 3)), triple_margin((2, 3, 4)))
    lower_bound = min(Fraction(1), *margins)
    ratio = 40000 // lower_bound + 1
    cutoff = 1
    while Fraction(4, 3) ** cutoff <= ratio * ratio:
        cutoff += 1
    print(f'Triple margins {margins}, R={ratio}, gap cutoff={cutoff}.', flush=True)
    profiles = [tuple((0, power) for power in range(4))]
    for gap in range(1, cutoff):
        triple = ((0, 0), (0, 1), (0, 2))
        profiles.append(triple + ((gap, 0),))
        profiles.append(((0, 0),) + tuple((gap, power) for power in range(3)))
    check_profiles(profiles, 65537, args.singular)
    print(f'Checked {len(profiles)} bounded profiles at prime 65537.')


if __name__ == '__main__':
    main()
