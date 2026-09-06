#!/usr/bin/env python3
"""Exact MW16 homogeneous resultants and bounded scaling-prime support."""
import argparse
from pathlib import Path
import audit_r17_scaling_prime_support as determinant_engine
import compact_mw16_specialization as spec
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
D = ROOT / 'artifacts/local/elliptic-curves/mw16-scaling-prime-support-v1'
OUT = ART / 'mw16_scaling_prime_support_v1.json'
INPUT = spec.ATLAS


def sources():
    paths = [Path(__file__).resolve(), Path(determinant_engine.__file__), INPUT,
             ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',
             ROOT/'elliptic-curves/cas/verify_mw16_scaling_prime_support.sage']
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}


def prepare():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve MW16 support audit')
    checkpoint(D/'protocol.json', {
        'schema': 'elliptic-curves.mw16-scaling-prime-support.v1',
        'sources': sources(), 'families': 5, 'maximum_trial_prime': 131071,
        'wall_seconds': 120, 'rss_bytes': 1073741824,
        'gate': 'The current target-free higher-parameter scanners omit singular displayed reductions. '
                'For primitive homogeneous parameters, a prime dividing both short coefficients must '
                'divide their homogeneous resultant. Bound possible removable-scaling support across '
                'all five families before interpreting displayed singularity as intrinsic bad reduction.',
        'scope': 'Exactly five degree8/12 homogeneous Sylvester determinants, bounded trial division '
                 'through131071 covering all active selection and withheld prime bands, and a separate '
                 'Sage polynomial-resultant replay. No integer factorization beyond trial division, '
                 'parameter enumeration, point search, new score or change to a running selection. '
                 'A remaining cofactor is preserved as UNKNOWN. Resultant divisibility is necessary, '
                 'not sufficient, for a common projective root or removable scale. Primes2,3 require '
                 'separate local-model treatment. No exact rank or new curve follows.'})


def expected():
    p = cert.read(D/'protocol.json')
    if p['sources'] != sources():
        raise ArithmeticError('frozen support inputs changed')
    primes = _primes_up_to(p['maximum_trial_prime'])
    rows = []
    for f in cert.read(INPUT)['families']:
        a = list(map(int, f['A_coefficients_low_to_high']))
        b = list(map(int, f['B_coefficients_low_to_high']))
        if len(a) != 9 or len(b) != 13 or not a[-1] or not b[-1]:
            raise ArithmeticError('full degree8/12 atlas required')
        result = determinant_engine.determinant(determinant_engine.sylvester(a, b))
        if not result:
            raise ArithmeticError('zero resultant leaves unrestricted prime support')
        remaining = abs(result)
        factors = []
        for q in primes:
            e = 0
            while remaining % q == 0:
                remaining //= q
                e += 1
            if e:
                factors.append([q, e])
        product = remaining
        for q, e in factors:
            product *= q**e
        if product != abs(result):
            raise ArithmeticError('factor reconstruction differs')
        rows.append({'family': f['fibration_id'], 'homogeneous_resultant': str(result),
                     'trial_prime_factors': factors, 'remaining_cofactor': str(remaining),
                     'complete_prime_support': remaining == 1,
                     'score_range_support': [q for q, e in factors if q >= 5]})
    if len(rows) != 5:
        raise ArithmeticError('complete five-family audit required')
    return {'schema': 'elliptic-curves.mw16-scaling-prime-support-result.v1',
            'status': 'PASS_EXACT_BOUNDED_SUPPORT', 'sources': sources(),
            'protocol_sha256': cert.hashed(D/'protocol.json'), 'rows': rows,
            'claim_boundary': p['scope']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'build', 'check'])
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare()
    else:
        data = expected()
        if args.stage == 'check':
            if cert.read(OUT) != data:
                raise ArithmeticError('MW16 support replay differs')
        else:
            if OUT.exists():
                raise FileExistsError('preserve MW16 support certificate')
            checkpoint(OUT, data)
        print(data['status'], [(r['family'], r['score_range_support'],
                              r['complete_prime_support']) for r in data['rows']], flush=True)
