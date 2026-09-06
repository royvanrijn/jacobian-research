#!/usr/bin/env python3
"""Exact Brumer--Kramer class-2-rank target for the MW16 fibre3/17.

This applies an existing theorem; it does not compute a class-group upper
bound or claim a new general descent algorithm. The ramification proof uses
only exact integer valuations, independent of the profiled PARI maximal order.
"""
import argparse
import json
from math import gcd
from pathlib import Path

import certify_small_conductor_curve as original

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
PROOF = ART / 'small_conductor_rank22_proof_v1.json'
OUTPUT = ART / 'small_conductor_class_rank_target_v1.json'


def valuation(value, prime):
    if not value:
        raise ArithmeticError('zero coefficient needs separate valuation logic')
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def expected():
    proof = original.cert.read(PROOF)
    original.verify(proof)
    if proof['integral_model'][:3] != ['1', '0', '0'] or proof['rank_lower_bound'] != 22:
        raise ArithmeticError('fixed model and rank certificate required')
    a, b = map(int, proof['integral_model'][3:])
    coefficients = [64*b, 16*a, 1, 1]
    c, d = coefficients[0], coefficients[1]
    discriminant = d*d - 4*d**3 - 4*c - 27*c*c + 18*d*c
    if discriminant != 256*int(proof['discriminant']):
        raise ArithmeticError('2-division cubic discriminant identity failed')
    shift = 11
    translated = [c+d*shift+shift**2+shift**3,
                  d+2*shift+3*shift**2, 1+3*shift, 1]
    vals = [valuation(v, 17) for v in translated]
    if vals != [2, 2, 1, 0]:
        raise ArithmeticError('fixed Newton polygon differs')
    # All middle points lie strictly above the segment (0,2)--(3,0),
    # whose relatively prime rise and run force root valuation2/3.
    if not (3*vals[1] > 4 and 3*vals[2] > 2 and gcd(vals[0], 3) == 1):
        raise ArithmeticError('one irreducible slope required')
    multiplicative_even = []
    additive = []
    for row in proof['local_proof']:
        prime = int(row['prime'])
        if row['conductor_valuation'] == 1:
            if row['discriminant_valuation'] % 2 == 0:
                multiplicative_even.append(prime)
        elif row['conductor_valuation'] == 2:
            additive.append(prime)
        else:
            raise ArithmeticError('unhandled local reduction')
    if additive != [17] or multiplicative_even != [2, 5, 13, 19, 71]:
        raise ArithmeticError('complete bad-place partition differs')
    u = 2 if int(proof['discriminant']) > 0 else 1
    n = len(multiplicative_even) + (1 - 1)
    offset = u + n
    lower = proof['rank_lower_bound'] - offset
    return {
        'schema': 'elliptic-curves.small-conductor-brumer-kramer-target.v1',
        'status': 'PASS', 'family': proof['family'], 'parameter': proof['parameter'],
        'sources': {str(p.relative_to(ROOT)): original.cert.hashed(p)
                    for p in [PROOF, Path(__file__).resolve()]},
        'integral_model': proof['integral_model'],
        'two_division_cubic_ascending': list(map(str, coefficients)),
        'cubic_discriminant': str(discriminant), 'curve_discriminant_ratio': 256,
        'generator_transport': 'If alpha is the x-coordinate of a nonzero2-torsion point on the displayed model, theta=4*alpha satisfies F(theta)=0.',
        'prime17_certificate': {
            'translation': shift, 'translated_coefficients_ascending': list(map(str, translated)),
            'coefficient_valuations': vals, 'root_valuation': '2/3',
            'ramification_index': 3, 'residue_degree': 1, 'number_of_primes_above17': 1,
            'argument': 'The sole lower Newton-polygon edge joins(0,2)and(3,0). Every root has valuation2/3. A rational17-adic root is impossible, so the cubic is irreducible overQ17 and hence Q. The root valuation forces the degree3 extension to have ramification index3, leaving one prime above17.'},
        'brumer_kramer': {
            'theorem_reference': 'Brumer--Kramer, The rank of elliptic curves, Duke Math.J.44(1977), Proposition7.1; restated as Proposition3.1 in Klagsbrun--Sherman--Weigandt.',
            'reference_url': 'https://arxiv.org/html/1606.07178#S3.SS1',
            'multiplicative_even_discriminant_primes': multiplicative_even,
            'additive_primes': additive, 'u': u, 'n': n, 'offset': offset,
            'inequality': '22 <= rank(E/Q) <= dim_F2 Sel_2(E/Q) <= dim_F2 Cl(K)[2] + 7'},
        'unconditional_class_group_two_rank_lower_bound': lower,
        'sufficient_class_group_two_rank_upper_bound_for_exact_rank22': lower,
        'necessary_class_group_two_rank_lower_bound_for_rank23': 23-offset,
        'class_group_two_rank_upper_bound': None, 'rank_upper_bound': None,
        'claim_boundary': 'Unconditional class2-rank lower bound15 and a sufficient exact-rank criterion using an existing theorem. No class-group upper bound, full2-Selmer computation, exact rank22, nonexistence of a23rd point or new general descent algorithm is established.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = expected()
    if args.check:
        if json.loads(json.dumps(result)) != original.cert.read(OUTPUT):
            raise ArithmeticError('class-rank target certificate differs')
    else:
        if OUTPUT.exists():
            raise FileExistsError('preserve class-rank target proof')
        original.cert.write(OUTPUT, result)
    print('PROVED: Cl(K)[2] dimension >=15; upper bound15 would prove exact rank22')
