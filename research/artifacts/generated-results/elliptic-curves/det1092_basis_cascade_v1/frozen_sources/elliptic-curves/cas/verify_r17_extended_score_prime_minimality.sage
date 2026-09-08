#!/usr/bin/env sage-python
"""Independent polynomial resultants and complete Sage prime-range valuations."""
import argparse
import json
import sys
from pathlib import Path
from sage.all import ZZ, PolynomialRing, prime_range

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'elliptic-curves/cas'))
import extend_r17_score_prime_minimality as audit
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint

OUT = audit.ART / 'r17_extended_score_prime_minimality_sage_replay_v1.json'


def expected():
    data = cert.read(audit.OUT)
    if data['status'] != 'PASS' or data['sources'] != audit.sources():
        raise ArithmeticError('source-bound range proof required')
    protocol = cert.read(audit.D / 'protocol.json')
    if data['protocol_sha256'] != cert.hashed(audit.D / 'protocol.json') or protocol['sources'] != audit.sources():
        raise ArithmeticError('protocol binding differs')
    primes = list(prime_range(998, 131072))
    if list(map(int, primes)) != data['additional_primes']:
        raise ArithmeticError('independent complete prime interval differs')
    families = cert.read(audit.support.INPUT)['families']
    old = cert.read(audit.support.OUT)
    if len(families) != 6 or len(data['rows']) != 6:
        raise ArithmeticError('all six families required')
    R = PolynomialRing(ZZ, 't')
    rows = []
    for f, previous, row in zip(families, old['rows'], data['rows']):
        A = R(list(map(ZZ, f['A_coefficients_low_to_high'])))
        B = R(list(map(ZZ, f['B_coefficients_low_to_high'])))
        result = A.resultant(B)
        if A.degree() != 8 or B.degree() != 12 or not result or str(result) != row['homogeneous_resultant']:
            raise ArithmeticError('independent homogeneous resultant differs')
        if f['family'] != row['family'] or previous['family'] != row['family']:
            raise ArithmeticError('family correspondence differs')
        remaining = abs(result)
        for p in prime_range(2, 998):
            remaining //= p ** remaining.valuation(p)
        if str(remaining) != previous['remaining_cofactor']:
            raise ArithmeticError('earlier prime stripping differs')
        factors = []
        for p in primes:
            exponent = result.valuation(p)
            if exponent >= 4:
                raise ArithmeticError('unexcluded scale prime')
            if exponent:
                factors.append([int(p), int(exponent)])
                remaining //= p ** exponent
        if factors != row['additional_trial_factors'] or str(remaining) != row['remaining_cofactor']:
            raise ArithmeticError('new prime valuations or remainder differ')
        rows.append({'family': f['family'], 'additional_trial_factors': factors,
                     'independent_polynomial_resultant': True, 'prime_pairs_checked': len(primes)})
    if data['additional_family_prime_pairs'] != 6 * len(primes) or data['total_non13_family_prime_pairs'] != 6 * len(list(prime_range(5, 131072))) - 6:
        raise ArithmeticError('prime-pair accounting differs')
    return {'schema': 'elliptic-curves.r17-score-prime-range-sage-replay.v1', 'status': 'PASS',
            'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in [Path(__file__).resolve(), audit.OUT, audit.support.INPUT, audit.support.OUT]},
            'rows': rows, 'claim_boundary': 'Independent Sage polynomial-resultant and complete prime-valuation replay for the six new intervals1009through131071. Earlier residue-tree exclusions through997 and the exceptional13 classification are preserved inputs, not rerun here. No new parameter or point search, exact conductor or new rank.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--check', action='store_true'); args = p.parse_args()
    data = expected()
    if args.check:
        if cert.read(OUT) != json.loads(json.dumps(data)):
            raise ArithmeticError('independent prime-range replay differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve independent range replay')
        checkpoint(OUT, data)
    print('INDEPENDENT SIX R17 FULL-SCORE PRIME INTERVALS PASS', flush=True)
