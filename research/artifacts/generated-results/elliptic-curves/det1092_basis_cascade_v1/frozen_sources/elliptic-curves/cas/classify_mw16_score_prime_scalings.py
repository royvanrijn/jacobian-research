#!/usr/bin/env python3
"""Classify the 33 possible MW16 scaling pairs in the active score-prime range."""
import argparse
from pathlib import Path
import audit_mw16_scaling_prime_support as support
import classify_r17_other_small_prime_scalings as engine
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint

ROOT = support.ROOT
D = ROOT/'artifacts/local/elliptic-curves/mw16-score-prime-scalings-v1'
OUT = support.ART/'mw16_score_prime_scalings_v1.json'


def sources():
    paths = [Path(__file__).resolve(), Path(engine.__file__), Path(support.__file__),
             support.INPUT, support.OUT, support.D/'ledger.json',
             ROOT/'elliptic-curves/cas/verify_mw16_score_prime_scalings.py']
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}


def roster():
    return [{'family': r['family'], 'prime': q, 'resultant_valuation': e}
            for r in cert.read(support.OUT)['rows']
            for q, e in r['trial_prime_factors'] if q >= 5 and e >= 4]


def prepare():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve fixed scaling classification')
    if cert.read(support.D/'ledger.json')['status'] != 'PASS' or len(roster()) != 33:
        raise ArithmeticError('independent support proof and all33 pairs required')
    checkpoint(D/'protocol.json', {
        'schema': 'elliptic-curves.mw16-score-prime-scalings.v1',
        'sources': sources(), 'pairs': roster(), 'maximum_depth': 6,
        'maximum_live_residues_per_chart': 4096, 'maximum_candidates_per_level': 200000,
        'wall_seconds': 300, 'rss_bytes': 2147483648,
        'gate': 'All five exact homogeneous resultants and an independent polynomial replay pass. '
                'Only33 family/prime pairs in5..131071 have resultant valuation at least4, a necessary '
                'condition for a removable short-model scale. All lie at primes at most101. '
                'Classify those necessary coefficient congruences before assessing omitted good primes.',
        'scope': 'Exactly33 finite coefficient residue trees, both affine and infinity charts, '
                 'depth at most6 with fixed node caps. A branch is excluded only by a necessary '
                 'coefficient congruence; an entire ball is admitted only by exact polynomial '
                 'divisibility by p^4 and p^6. Independent replay checks every branch partition and '
                 'admitted polynomial by Horner composition. Capped branches remain UNKNOWN. '
                 'No parameter search, point search, score change or adaptive extension. '
                 'Primes2,3 and primes beyond131071 are outside this classification. '
                 'An admitted removable scale does not by itself establish good reduction after scaling.'})


def protocol():
    p = cert.read(D/'protocol.json')
    if p['sources'] != sources() or p['pairs'] != roster():
        raise ArithmeticError('frozen residue classification changed')
    return p


def expected():
    p = protocol()
    families = {f['fibration_id']: f for f in cert.read(support.INPUT)['families']}
    rows = []
    for pair in p['pairs']:
        f = families[pair['family']]
        a = list(map(int, f['A_coefficients_low_to_high']))
        b = list(map(int, f['B_coefficients_low_to_high']))
        charts = [{'chart': name, **engine.chart(a, b, pair['prime'], infinity, p)}
                  for name, infinity in [('affine', False), ('infinity', True)]]
        complete = all(c['status'] == 'COMPLETE_RESIDUE_CLASSIFICATION' for c in charts)
        count = sum(len(c['scale_balls']) for c in charts)
        status = ('CLASSIFIED_SCALE_BALLS' if count else 'NO_REMOVABLE_SCALE') if complete else 'UNKNOWN_INCOMPLETE_CLASSIFICATION'
        rows.append({**pair, 'status': status, 'scale_balls': count, 'charts': charts})
    return {'schema': 'elliptic-curves.mw16-score-prime-scalings-result.v1',
            'status': 'PASS_BOUNDED_EXACT_CLASSIFICATION', 'sources': sources(),
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
                raise ArithmeticError('exact coefficient tree replay differs')
        else:
            if OUT.exists():
                raise FileExistsError('preserve scaling classification')
            checkpoint(OUT, data)
        print(data['status'], [(r['family'], r['prime'], r['status'], r['scale_balls']) for r in data['rows']], flush=True)
