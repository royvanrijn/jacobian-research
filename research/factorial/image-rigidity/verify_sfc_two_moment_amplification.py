#!/usr/bin/env python3
"""Exact controls for two-moment multiplication and modular rank transfer.

The actual two-term example is a modular obstruction to an unconditional
local transfer. It is not a characteristic-zero SFC counterexample.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, factorial, prod
from pathlib import Path
import platform

import flint
from flint import fmpz_poly, nmod_poly


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def coefficients(alpha, beta, order, marker=None):
    marker = marker or (0,)*len(alpha)
    return [comb(order, k)*prod(factorial((order-k)*a+k*b+r)
                              for a, b, r in zip(alpha, beta, marker))
            for k in range(order+1)]


def direct_moment(alpha, beta, order, ratio, marker=None):
    """Expand powers in the original spatial polynomial ring."""
    marker = marker or (0,)*len(alpha)
    power = {(0,)*len(alpha): 1}
    for _ in range(order):
        following = {}
        for exponents, value in power.items():
            for column, scalar in ((alpha, 1), (beta, ratio)):
                key = tuple(a+b for a, b in zip(exponents, column))
                following[key] = following.get(key, 0)+scalar*value
        power = following
    return sum(value*prod(factorial(a+r) for a, r in zip(exponents, marker))
               for exponents, value in power.items())


def rank_mod(matrix, prime):
    matrix = [[int(x) % prime for x in row] for row in matrix]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [x*inverse % prime for x in matrix[rank]]
        for i in range(rank+1, len(matrix)):
            multiplier = matrix[i][column]
            matrix[i] = [(x-multiplier*y) % prime
                         for x, y in zip(matrix[i], matrix[rank])]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def multiplication_rows(left, right, degree):
    rows = []
    for polynomial in (left, right):
        for shift in range(degree-len(polynomial)+2):
            row = [0]*(degree+1)
            row[shift:shift+len(polynomial)] = polynomial
            rows.append(row)
    return rows


def multiplication_controls():
    prime = 127
    t = nmod_poly([0, 1], prime)
    cases = [((t+1)**2, t**3+2),
             ((t+1)*(t+2), (t+1)*(t**2+3)),
             ((t+1)**2*(t+2), (t+1)**2*(t**2+3))]
    records = []
    for left, right in cases:
        m, n = left.degree(), right.degree()
        common_degree = left.gcd(right).degree()
        for degree in range(m+n-1, m+n+3):
            rows = multiplication_rows([int(x) for x in left], [int(x) for x in right], degree)
            rank = rank_mod(rows, prime)
            require(rank == degree+1-common_degree, "Macaulay corank differs from gcd degree.")
            records.append({"degrees": [m, n], "coefficient_degree": degree,
                            "gcd_degree": common_degree, "rank": rank,
                            "corank": degree+1-rank})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name(
        'sfc_two_moment_amplification_verification.json'))
    args = parser.parse_args()
    alpha, beta, n, prime, ratio = (1, 2), (3, 1), 14, 127, 17
    p, q = coefficients(alpha, beta, n), coefficients(alpha, beta, n+1)
    require(flint.fmpz(prime).is_prime(), 'The residue characteristic must be prime.')
    pm, qm = nmod_poly(p, prime), nmod_poly(q, prime)
    require(pm.degree() == n and qm.degree() == n+1 and p[0] % prime and q[0] % prime,
            'An endpoint coefficient vanished.')
    require(pm.gcd(qm) == nmod_poly([-ratio, 1], prime), 'The stated modular gcd changed.')
    values = []
    for order, marker in ((n, (0, 0)), (n+1, (0, 0)), (n, (1, 0)), (n, (2, 0))):
        direct = direct_moment(alpha, beta, order, ratio, marker)
        polynomial = coefficients(alpha, beta, order, marker)
        require(direct == sum(c*ratio**k for k, c in enumerate(polynomial)),
                'Spatial expansion and coefficient moment disagree.')
        values.append(direct % prime)
    require(values[:3] == [0, 0, 97], 'The modular marked-moment obstruction changed.')
    marked_matrix = [[values[0], values[2]], [values[2], values[3]]]
    determinant = (values[0]*values[3]-values[2]**2) % prime
    require(determinant == 116 and rank_mod(marked_matrix, prime) == 2,
            'The marked Hankel matrix is not full rank.')
    maximum_argument = max((n+1)*max(alpha+beta), n*max(alpha+beta)+2)
    require(maximum_argument < prime, 'Factorial reduction is no longer good.')
    require(fmpz_poly(p).gcd(fmpz_poly(q)).degree() == 0,
            'This local example unexpectedly has a characteristic-zero common factor.')

    ranks = []
    for degree in (28, 29, 32, 42):
        rows = multiplication_rows(p, q, degree)
        rank = rank_mod(rows, prime)
        require(degree+1-rank == 1, 'Multiplication closure gained an extra corank.')
        ranks.append({'coefficient_degree': degree, 'rows': len(rows),
                      'columns': degree+1, 'rank': rank, 'corank': 1})
    curvature = sum(Fraction((a-b)**2, max(a, b)) for a, b in zip(alpha, beta))
    require(curvature == Fraction(11, 6) and all(alpha) and all(beta),
            'The profile left the fully shared residual shape class.')

    result = {'status': 'PASS', 'python': platform.python_version(), 'python_flint': flint.__version__,
        'scope': 'local modular non-closure and finite controls of the written multiplication-rank theorem',
        'actual_profile': {'alpha': alpha, 'beta': beta, 'terminal': n+1,
            'fully_shared': True, 'shared_curvature': str(curvature),
            'prime': prime, 'coefficient_ratio': ratio, 'modular_gcd': 't-17',
            'M14_M15_XM14_X2M14_mod_127': values,
            'marked_hankel_matrix': marked_matrix, 'marked_hankel_determinant': determinant,
            'maximum_factorial_argument': maximum_argument,
            'characteristic_zero_gcd_degree': 0,
            'lower_moment_coefficients_mod_127': [int(x) for x in pm],
            'upper_moment_coefficients_mod_127': [int(x) for x in qm]},
        'actual_multiplication_ranks': ranks,
        'synthetic_multiplication_controls': multiplication_controls(),
        'negative_control': 'the proposed vanishing of the first spatial mark is rejected by residue 97',
        'not_established': ['a characteristic-zero SFC counterexample',
            'a transfer conditional on a global common factor and eligible primes',
            'Matsusaka seed propagation or optimized SFC determinant bounds'],
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'assurance': 'internal exact computation; no independent replay, formal verification or external review'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: actual moments 14,15 vanish modulo 127; first mark=97; Hankel determinant=116.')
    print('PASS: multiplication closure has corank one at all four tested coefficient degrees.')
    print('PASS: independent spatial expansion and 12 synthetic Macaulay controls.')


if __name__ == '__main__':
    main()
