#!/usr/bin/env python3
"""Exact controls for the elliptic class-specific prime eligibility sieve."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import gcd, isqrt
from pathlib import Path
import json
import platform


Q0 = 2709173
CLASSES = (1, 2, 5*Q0, 10*Q0)
D = 51730
A2, A4 = -10656, 27091730
P = (F(373248, 49), F(59459616, 343))
S = (F(92179201, 6400), F(553221611649, 512000))
T = (F(0), F(0))


def is_prime(integer):
    return integer >= 2 and all(integer % divisor
                               for divisor in range(2, isqrt(integer)+1))


def legendre(integer, prime):
    residue = pow(integer % prime, (prime-1)//2, prime)
    return 0 if residue == 0 else (1 if residue == 1 else -1)


def rational_square(value):
    return value >= 0 and isqrt(value.numerator)**2 == value.numerator and \
        isqrt(value.denominator)**2 == value.denominator


def valuation(value, prime):
    value = F(value)
    assert value
    numerator, denominator = value.numerator, value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def on_curve(point):
    if point is None:
        return True
    x, y = point
    return y*y == x*(x*x+A2*x+A4)


def add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if y+v == 0:
            return None
        slope = (3*x*x+2*A2*x+A4)/(2*y)
    else:
        slope = (v-y)/(u-x)
    new_x = slope*slope-A2-x-u
    return new_x, -y+slope*(x-new_x)


def multiply(point, integer):
    if integer < 0:
        return multiply((point[0], -point[1]), -integer)
    result = None
    while integer:
        if integer & 1:
            result = add(result, point)
        point = add(point, point)
        integer >>= 1
    return result


def run():
    checks = []

    def check(label, assertion):
        if not assertion:
            raise AssertionError(label)
        checks.append(label)

    x0 = F(25865, 128)
    check('exact_branch_Q_square', x0*x0+A2*x0+A4 == F(639735,128)**2)
    check('branch_x_squareclass', x0 == D*F(1,16)**2)
    check('branch_constant_factorization', 639735 == 3*5*42649)
    check('squarefree_branch_class_factorization', D == 2*5*7*739)
    check('constant_primality', all(is_prime(p) for p in (2,3,5,7,739,Q0,42649)))
    check('branch_class_prime_congruences', 7 % 3 == 1 and 739 % 3 == 1)
    check('ramified_class_prime_is_inert_and_five_mod_eight', Q0 % 3 == 2 and Q0 % 8 == 5)

    # Prime-coordinate order is (-1,2,3,5,7,739,q0).
    vectors = (0b0000101, 0b0000010, 0b1001000, 0b0111010)
    sums = set()
    for coefficients in product((0,1), repeat=4):
        value = 0
        for coefficient, vector in zip(coefficients, vectors):
            if coefficient:
                value ^= vector
        sums.add(value)
    check('four_squareclasses_are_independent', len(sums) == 16)
    patterns = list(product((-1,1), repeat=4))
    density = {}
    for a in CLASSES:
        denominator = numerator = 0
        for chi_minus3, chi_2, chi_5q, chi_d in patterns:
            chi_a = {1: 1, 2: chi_2, 5*Q0: chi_5q, 10*Q0: chi_2*chi_5q}[a]
            if chi_minus3 == -1 and chi_a == 1:
                denominator += 1
                numerator += chi_d == 1
        density[a] = (F(denominator,16), F(numerator,16))
    check('prime_eligibility_density_arithmetic', density == {
        1: (F(1,2),F(1,4)),
        2: (F(1,4),F(1,8)),
        5*Q0: (F(1,4),F(1,8)),
        10*Q0: (F(1,4),F(1,8)),
    })

    check('eleven_eligibility_and_absent_numerator_branch',
          [legendre(a,11) for a in CLASSES] == [1,-1,1,-1] and legendre(D,11) == -1)
    check('seventeen_eligibility_and_present_numerator_branch',
          [legendre(a,17) for a in CLASSES] == [1,1,-1,-1] and legendre(D,17) == 1)
    check('ramified_q0_and_class_two_cases',
          [legendre(a,Q0) for a in CLASSES] == [1,-1,0,0])

    primes = [p for p in range(11,252) if p % 3 == 2 and is_prime(p)]
    finite_field_rows = []
    for p in primes:
        branches = []
        for x in range(p):
            if (128*x-25865) % p == 0:
                rhs = x*(x*x+A2*x+A4) % p
                ys = [y for y in range(p) if y*y % p == rhs]
                branches.extend((x,y) for y in ys)
                if ys:
                    assert legendre(x,p) == 1 and legendre(D,p) == 1
        expected = 2 if legendre(D,p) == 1 else 0
        assert len(branches) == expected
        finite_field_rows.append({'p': p, 'branch_points': [list(point) for point in branches]})
    check('complete_finite_field_branch_eligibility_scans', bool(finite_field_rows))

    check('exact_rational_generator_points', all(on_curve(point) for point in (P,S,T)))
    rational_points = set()
    for m,n,e in product(range(-2,3), range(-2,3), range(2)):
        point = add(add(multiply(P,m),multiply(S,n)),multiply(T,e))
        assert on_curve(point)
        if point is not None and point[0] and 128*point[0]-25865:
            rational_points.add(point)
    seen_classes, cases = set(), 0
    for x,y in rational_points:
        classes = [a for a in CLASSES if rational_square(x/a)]
        assert len(classes) == 1
        a = classes[0]
        seen_classes.add(a)
        m, c = x.numerator, isqrt(x.denominator)
        assert c*c == x.denominator and gcd(m,c) == 1
        n = 128*m-25865*c*c
        f = 128*x-25865
        for p in primes+[Q0]:
            k = valuation(f,p)
            if legendre(a,p) != 1:
                assert k == 0 and c % p and n % p
            if legendre(D,p) == -1:
                assert n % p
            if k < 0:
                assert k == -2*valuation(c,p) and legendre(a,p) == 1
            elif k > 0:
                assert k == valuation(n,p) and legendre(a,p) == legendre(D,p) == 1
            cases += 1
    check('bounded_rational_subgroup_contains_all_four_classes', seen_classes == set(CLASSES))
    check('bounded_rational_integer_sieve_controls', cases > 0)

    base = Path(__file__).resolve().parent
    proof = base/'ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md'
    source_names = [
        'ELLIPTIC_EXCEPTION_ARITHMETIC.md',
        'ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md',
        'ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md',
        'ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md',
        'sources/milne-ant300.pdf',
    ]
    check('retained_proof_inputs_exist', all((base/name).is_file() for name in source_names))
    return {
        'status': 'pass',
        'checks_passed': len(checks),
        'exact_controls': checks,
        'proof_path': proof.name,
        'proof_sha256': sha256(proof.read_bytes()).hexdigest(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'command': '.venv/bin/python research/arithmetic/verify_elliptic_exception_prime_eligibility.py --write-receipt',
        'python_version': platform.python_version(),
        'source_hashes': {name: sha256((base/name).read_bytes()).hexdigest() for name in source_names},
        'finite_field_controls': finite_field_rows,
        'rational_subgroup_domain': {'m': [-2,2], 'n': [-2,2], 'e': [0,1],
                                     'generators': ['P','S','T'],
                                     'nonzero_distinct_points': len(rational_points),
                                     'prime_point_cases': cases,
                                     'additional_prime': Q0},
        'eligibility_density_arithmetic': {str(a): [str(v) for v in pair]
                                          for a,pair in density.items()},
        'external_imports': ['Milne ANT300 Proposition7.31/Theorem7.32 Hensel',
                             'Milne ANT300 Theorem8.31/Corollary8.32 Chebotarev'],
        'unrestricted_proof_not_finite_experiment': True,
        'formal_verification': False,
        'external_review': False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    receipt = run()
    if args.write_receipt:
        path = Path(__file__).with_name('elliptic_exception_prime_eligibility_checks.json')
        path.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(f"PASS: {receipt['checks_passed']} exact controls; "
          f"{receipt['rational_subgroup_domain']['prime_point_cases']} bounded prime/point cases")
