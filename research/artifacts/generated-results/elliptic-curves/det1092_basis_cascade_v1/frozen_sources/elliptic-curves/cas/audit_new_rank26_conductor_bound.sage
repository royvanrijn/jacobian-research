#!/usr/bin/env sage-python
"""One bounded local-conductor audit; no discriminant factorization campaign."""
from pathlib import Path
import sys
from sage.all import EllipticCurve, QQ, ZZ, prime_range
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint


def audit():
    source = ROOT/'artifacts/generated-results/elliptic-curves/new_compact_rank26_proof_v1.json'
    database = ROOT/'artifacts/local/elliptic-curves/next12-current-catalogue-v1/database.json'
    out = ROOT/'artifacts/local/elliptic-curves/new-rank26-conductor-bound-v1/result.json'
    if out.exists():
        raise FileExistsError('preserve conductor bound')
    data = cert.read(source); E = EllipticCurve([QQ(a) for a in data['minimal_curve']])
    if E.c4().gcd(E.c6()) != 27:
        raise ArithmeticError('minimality input differs')
    remaining = abs(ZZ(E.discriminant())); known = ZZ(1); local = []
    for p in prime_range(2, 10001):
        e = int(remaining.valuation(p))
        if not e:
            continue
        remaining //= p**e
        # Away from3, c4 is a unit at every bad prime, hence multiplicative.
        if p != 3:
            if E.c4() % p == 0:
                raise ArithmeticError('unexpected additive prime')
            f, kind = 1, 'multiplicative'
        else:
            ld = E.local_data(p, proof=True)
            f, kind = int(ld.conductor_valuation()), str(ld.kodaira_symbol())
        known *= p**f
        local.append({'prime': str(p), 'discriminant_valuation': e,
                      'conductor_valuation': f, 'reduction': kind})
    upper = known*remaining
    curves = cert.read(database)['curves']
    listed = sorted((ZZ(r['conductor']), r['id']) for r in curves
                    if r['rank_lower_bound'] >= 26 and r['conductor'])
    record = {'schema': 'elliptic-curves.new-rank26-conductor-bound.v1',
        'status': 'COMPLETE_BOUNDED_LOCAL_AUDIT', 'curve': data['minimal_curve'],
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (Path(__file__).resolve(), source, database)},
        'discriminant': str(E.discriminant()), 'local_data': local,
        'unprocessed_cofactor': str(remaining), 'conductor_upper_bound': str(upper),
        'exact_conductor': remaining == 1,
        'listed_rank26_minimum': {'id': listed[0][1], 'conductor': str(listed[0][0])},
        'upper_bound_beats_listed_minimum': bool(upper < listed[0][0]),
        'missing_rank26_conductors': [r['id'] for r in curves if r['rank_lower_bound'] >= 26 and not r['conductor']],
        'claim_boundary': 'The exact conductor divides this upper bound: all residual primes have multiplicative reduction. The cofactor is not factored or asserted prime. This bound does not establish a record; failure to beat the listed minimum does not exclude a smaller actual conductor.'}
    checkpoint(out, record)
    print('EXACT LOCAL CONDUCTOR UPPER BOUND', upper, 'unprocessed digits', len(str(remaining)),
          'beats listed minimum', record['upper_bound_beats_listed_minimum'], flush=True)


if __name__ == '__main__':
    audit()
