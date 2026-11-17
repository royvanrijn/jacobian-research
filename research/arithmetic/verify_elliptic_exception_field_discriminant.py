#!/usr/bin/env python3
"""Verify the EEH cubic maximal-order discriminant using exact standard-library arithmetic."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
import math
from pathlib import Path
import platform


def multiply(a, b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def evaluate(coeff, x):
    value = 0
    for c in reversed(coeff):
        value = value*x+c
    return value


def cubic_discriminant(coeff):
    d, c, b, a = coeff
    return b*b*c*c-4*a*c**3-4*b**3*d-27*a*a*d*d+18*a*b*c*d


def determinant3(m):
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])
            - m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, math.isqrt(n)+1, 2))


def factor_product(factors):
    return math.prod(int(p)**int(e) for p, e in factors)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path,
                        default=Path(__file__).with_name('elliptic_exception_field_discriminant_certificate.json'))
    parser.add_argument('--write-receipt', action='store_true')
    args = parser.parse_args()
    cert = json.loads(args.certificate.read_text(encoding='utf-8'))
    here = Path(__file__).resolve().parent
    old_path = here/cert['original_hasse_certificate']
    old = json.loads(old_path.read_text(encoding='utf-8'))
    checks = []

    def require(name, value, **details):
        if not value:
            raise AssertionError(name)
        checks.append({'name': name, 'passed': True, **details})

    R = list(map(F, cert['reduced_cubic_low_to_high']))
    H = list(map(F, old['integral_cubic_coefficients_low_to_high']))
    old_h = list(map(F, old['cubic_coefficients_low_to_high']))
    old_scale = F(old['integer_scale'])
    require('same_original_integral_model', H == [a*old_scale**(3-i) for i, a in enumerate(old_h)])
    require('monic_cubic_models', len(R) == len(H) == 4 and R[-1] == H[-1] == 1)
    p = int(cert['irreducibility_control_prime'])
    require('irreducibility_control_prime', is_prime(p) and p == int(old['irreducible_mod_prime']))
    require('reduced_cubic_irreducible_mod_prime', all(evaluate(R, r) % p for r in range(p)), prime=p)
    require('original_cubic_irreducible_mod_prime', all(evaluate(H, r) % p for r in range(p)), prime=p)

    def reduce(poly):
        value = list(map(F, poly))
        while len(value) > 3:
            degree = len(value)-1
            top = value.pop()
            for j in range(3):
                value[degree-3+j] -= top*R[j]
        return value+[F(0)]*(3-len(value))

    def product(a, b):
        return reduce(multiply(a, b))

    denominator = int(cert['original_integral_root_denominator'])
    alpha = [F(int(a), denominator) for a in cert['original_integral_root_numerators_low_to_high']]
    image = [F(0)]*3
    for coefficient in reversed(H):
        image = product(image, alpha)
        image[0] += coefficient
    require('exact_original_root_bridge_mod_R', image == [0, 0, 0])

    I = int(cert['basis_denominator'])
    a = int(cert['basis_linear_subtrahend'])
    b = int(cert['basis_constant_subtrahend'])
    basis = [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(-b, I), F(-a, I), F(1, I)]]
    require('rank_three_basis', I > 0)

    def coordinates(poly):
        v0, v1, v2 = reduce(poly)
        return [v0+b*v2, v1+a*v2, I*v2]

    table = {tuple(row['pair']): list(map(int, row['coordinates'])) for row in cert['multiplication_table']}
    require('complete_nontrivial_basis_products', set(table) == {(1, 1), (1, 2), (2, 2)})
    for pair, expected in table.items():
        actual = coordinates(product(basis[pair[0]], basis[pair[1]]))
        require('integral_order_product_'+str(pair), actual == expected and all(v.denominator == 1 for v in actual))

    matrices = []
    for beta in basis:
        columns = [coordinates(product(beta, gamma)) for gamma in basis]
        matrices.append([[columns[col][row] for col in range(3)] for row in range(3)])
    trace_pairing = [[sum(matrices[i][r][k]*matrices[j][k][r] for r in range(3) for k in range(3))
                      for j in range(3)] for i in range(3)]
    order_disc = determinant3(trace_pairing)
    D = int(cert['field_discriminant'])
    require('trace_pairing_discriminant', order_disc == D, discriminant=str(D))
    require('power_basis_discriminant_comparison', cubic_discriminant(R) == D*I*I)

    factors = cert['field_discriminant_factorization']
    primes = [int(p) for p, _ in factors]
    require('complete_order_discriminant_factorization', cert['field_discriminant_sign'] == -1
            and -factor_product(factors) == D and primes == sorted(set(primes)))
    require('all_discriminant_primes_are_prime', all(is_prime(p) for p in primes), primes=primes)
    require('positive_factor_exponents', all(int(e) > 0 for _, e in factors))
    squarefree = -math.prod(int(p) for p, e in factors if int(e) % 2)
    dK = squarefree if squarefree % 4 == 1 else 4*squarefree
    require('quadratic_fundamental_discriminant', dK == int(cert['quadratic_fundamental_discriminant'])
            and squarefree == int(old['quadratic_radicand']))

    f = int(cert['conductor_factor'])
    ff = cert['conductor_factorization']
    require('field_discriminant_squarepart_identity', D == dK*f*f)
    require('squarepart_complete_factorization', factor_product(ff) == f
            and all(is_prime(int(p)) and int(e) == 1 for p, e in ff))

    eis_primes = set()
    for row in cert['eisenstein_shifts']:
        prime, r = int(row['prime']), int(row['shift'])
        unit = int(row['constant_over_prime_mod_prime'])
        constant = int(evaluate(R, r))
        linear = int(R[1]+2*R[2]*r+3*r*r)
        quadratic = int(R[2]+3*r)
        require('Eisenstein_shift_'+str(prime), is_prime(prime) and 0 <= r < prime
                and constant % prime == 0 and constant//prime % prime == unit
                and unit % prime != 0 and linear % prime == 0 and quadratic % prime == 0,
                prime=prime, shift=str(r), constant_unit=str(unit))
        eis_primes.add(prime)
    exponents = {int(p): int(e) for p, e in factors}
    automatic_odd = {p for p, e in exponents.items() if p % 2 and e == 1}
    require('maximality_prime_coverage', eis_primes | automatic_odd | {2} == set(primes)
            and all(exponents[p] == 2 for p in eis_primes),
            eisenstein_primes=sorted(eis_primes), valuation_one_odd_primes=sorted(automatic_odd))
    require('two_adic_Stickelberger_index_obstruction', exponents[2] == 3 and (D//4) % 4 == 2)

    original_index = factor_product(cert['original_order_index_factorization'])
    require('original_order_index', cubic_discriminant(H) == D*original_index**2,
            index_factorization=cert['original_order_index_factorization'])
    require('Hasse_discriminant_predicate', all(int(p) % 3 == 1 for p, _ in ff) and f % 3 != 0)
    source = here/'ELLIPTIC_EXCEPTION_FIELD_DISCRIMINANT.md'
    milne = here/'sources/milne-ant300.pdf'
    report = {'status': 'passed', 'python': platform.python_version(),
              'command': 'python3 research/arithmetic/verify_elliptic_exception_field_discriminant.py --write-receipt',
              'scope': ('Exact field-isomorphism, integral multiplication table, trace discriminant, '
                        'prime factorization and local Eisenstein certificates. Maximality uses the '
                        'written local-ring and discriminant theorems. No PARI computation, formal '
                        'verification, external review or elliptic Hasse-infinitude claim is assumed.'),
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'source_sha256': sha256(source.read_bytes()).hexdigest(),
              'certificate_sha256': sha256(args.certificate.read_bytes()).hexdigest(),
              'original_certificate_sha256': sha256(old_path.read_bytes()).hexdigest(),
              'primary_Milne_sha256': sha256(milne.read_bytes()).hexdigest(),
              'check_count': len(checks), 'checks': checks,
              'field_discriminant': str(D), 'quadratic_discriminant': str(dK), 'squarepart': str(f),
              'squarepart_factors': ff, 'formal_verification': False,
              'independent_whole_proof_replay': False, 'external_review': False}
    if args.write_receipt:
        (here/'elliptic_exception_field_discriminant_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': report['status'], 'check_count': len(checks),
                      'field_discriminant': str(D), 'squarepart': str(f)}))


if __name__ == '__main__':
    main()
