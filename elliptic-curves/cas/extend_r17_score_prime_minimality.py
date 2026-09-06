#!/usr/bin/env python3
"""Extend the proved R17 scaling-prime exclusion across the full score range."""
import argparse
import json
from pathlib import Path

import audit_r17_scaling_prime_support as support
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint

ROOT = support.ROOT
ART = support.ART
OLD = ART / 'r17_small_prime_minimality_v1.json'
OUT = ART / 'r17_extended_score_prime_minimality_v1.json'
D = ROOT / 'artifacts/local/elliptic-curves/r17-extended-score-prime-minimality-v1'


def sources():
    paths = [Path(__file__).resolve(), Path(support.__file__).resolve(),
             support.INPUT, support.OUT, OLD,
             ROOT / 'elliptic-curves/cas/mod2_reduction_independence.py',
             ROOT / 'elliptic-curves/cas/verify_r17_extended_score_prime_minimality.sage']
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}


def prepare():
    if (D / 'protocol.json').exists():
        raise FileExistsError('preserve fixed prime-range extension')
    checkpoint(D / 'protocol.json', {
        'schema': 'elliptic-curves.r17-score-prime-range-protocol.v1',
        'sources': sources(), 'old_maximum': 997, 'new_maximum': 131071,
        'families': 6, 'wall_seconds_per_stage': 120, 'rss_bytes': 1073741824,
        'scope': 'Extend exact resultant trial valuations from997 through131071 on all six existing compact R17 models. A removable scale at p>=5 requires valuation of the homogeneous coefficient resultant at least4. Combine only with the preserved proved exclusions through997; the separate13 classification is unchanged. Recompute all six determinants and independently verify polynomial resultants and prime valuations in Sage. No parameter scan, point search, new cache, factorization escalation, score change or new-rank claim. Primes2,3 and above131071 remain outside scope; unfactored residual cofactors remain UNKNOWN.'})


def expected():
    protocol = cert.read(D / 'protocol.json')
    if protocol['sources'] != sources():
        raise ArithmeticError('fixed extension sources changed')
    original = support.expected()
    if json.loads(json.dumps(original)) != cert.read(support.OUT):
        raise ArithmeticError('original exact resultants differ')
    old = cert.read(OLD)
    if old['status'] != 'PASS' or old['prime_family_pairs'] != 990:
        raise ArithmeticError('proved earlier prime interval required')
    for data in [old, original]:
        if any(cert.hashed(ROOT / n) != h for n, h in data['sources'].items()):
            raise ArithmeticError('prior proof binding changed')
    primes = [p for p in _primes_up_to(protocol['new_maximum']) if p > protocol['old_maximum']]
    rows = []
    for row in original['rows']:
        remaining = int(row['remaining_cofactor'])
        factors = []
        for prime in primes:
            exponent = 0
            while remaining % prime == 0:
                remaining //= prime
                exponent += 1
            if exponent:
                factors.append([prime, exponent])
        if any(e >= 4 for p, e in factors):
            raise ArithmeticError('new scaling candidate requires a separate residue proof')
        rows.append({'family': row['family'], 'homogeneous_resultant': row['homogeneous_resultant'],
                     'additional_trial_factors': factors, 'remaining_cofactor': str(remaining),
                     'all_new_prime_valuations_below4': True})
    if len(rows) != 6 or [r['family'] for r in rows] != [r['family'] for r in old['rows']]:
        raise ArithmeticError('complete earlier six-family correspondence required')
    return {'schema': 'elliptic-curves.r17-extended-score-prime-minimality.v1',
            'status': 'PASS', 'sources': sources(), 'protocol_sha256': cert.hashed(D / 'protocol.json'),
            'additional_primes': primes, 'additional_family_prime_pairs': 6 * len(primes),
            'total_non13_family_prime_pairs': 990 + 6 * len(primes), 'rows': rows,
            'claim_boundary': protocol['scope']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'build', 'check'])
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare()
    else:
        data = expected()
        if args.stage == 'build':
            if OUT.exists():
                raise FileExistsError('preserve exact range-extension certificate')
            checkpoint(OUT, data)
        elif cert.read(OUT) != data:
            raise ArithmeticError('extended prime proof differs')
        print('R17 SCORE-PRIME EXCLUSION PASS', data['total_non13_family_prime_pairs'], flush=True)
