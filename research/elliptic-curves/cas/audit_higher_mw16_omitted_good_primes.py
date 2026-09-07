#!/usr/bin/env python3
"""Measure omitted good-prime terms on a frozen scalar roster, without reselection."""
import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import certify_compact_r17_candidates as cert
import classify_mw16_score_prime_scalings as scaling
from research_runtime.store import checkpoint

ROOT = scaling.ROOT
ART = scaling.support.ART
D = ROOT/'artifacts/local/elliptic-curves/joint-mw16-higher-scores-v1'
OUT = ART/'higher_mw16_omitted_good_primes_v1.json'


def trace(a, b, prime):
    a, b = a % prime, b % prime
    if (4*a**3+27*b*b) % prime == 0:
        return None
    character = 0
    for x in range(prime):
        z = (x**3+a*x+b) % prime
        s = pow(z, (prime-1)//2, prime)
        character += 0 if s == 0 else 1 if s == 1 else -1
    ap = -character
    count = 1 + sum((y*y-x**3-a*x-b) % prime == 0
                    for x in range(prime) for y in range(prime))
    if ap != prime+1-count or ap*ap > 4*prime:
        raise ArithmeticError('independent full finite-point enumeration differs')
    return ap


def expected():
    paths = [Path(__file__).resolve(), Path(cert.__file__), scaling.OUT,
             scaling.D/'ledger.json', scaling.support.OUT, scaling.support.D/'ledger.json',
             D/'protocol.json', D/'result.json', D/'controller/ledger.json',
             ART/'joint_mw16_higher_selection_v1.json']
    if any(cert.read(p)['status'] != 'PASS' for p in
           [scaling.D/'ledger.json', scaling.support.D/'ledger.json', D/'controller/ledger.json']):
        raise ArithmeticError('complete support, classification and scalar gates required')
    p = cert.read(D/'protocol.json')
    scores = cert.read(D/'result.json')
    selected = cert.read(paths[-1])
    if any(cert.hashed(ROOT/n) != h for n, h in p['sources'].items()):
        raise ArithmeticError('frozen scalar inputs changed')
    if (scores['status'] != 'COMPLETE_FROZEN10240' or len(p['rows']) != 10240
            or len(scores['rows']) != 10240 or len(selected['selected']) != 60
            or scores['protocol_sha256'] != cert.hashed(D/'protocol.json')
            or selected['scores_sha256'] != cert.hashed(D/'result.json')):
        raise ArithmeticError('complete frozen scalar correspondence required')
    candidates = defaultdict(list)
    for r in cert.read(scaling.OUT)['rows']:
        if r['status'] == 'CLASSIFIED_SCALE_BALLS':
            candidates[r['family']].append(r['prime'])
        elif r['status'] != 'NO_REMOVABLE_SCALE':
            raise ArithmeticError('incomplete prime classification')
    if sorted({q for qs in candidates.values() for q in qs}) != [5, 13]:
        raise ArithmeticError('fixed small-prime support changed')
    rows, groups = [], defaultdict(list)
    traces = {}
    for row, score in zip(p['rows'], scores['rows']):
        if row['id'] != score['id']:
            raise ArithmeticError('scalar ordering differs')
        original_a, original_b = map(int, row['model'][3:])
        if not 4*original_a**3+27*original_b**2:
            raise ArithmeticError('singular curve')
        restored, scales = [], []
        for prime in candidates[row['family']]:
            a, b, exponent = original_a, original_b, 0
            while a % prime**4 == 0 and b % prime**6 == 0:
                a //= prime**4
                b //= prime**6
                exponent += 1
            if not exponent:
                continue
            scales.append([prime, exponent])
            key = (a % prime, b % prime, prime)
            if key not in traces:
                traces[key] = trace(*key)
            ap = traces[key]
            if ap is not None:
                # Preserve the existing short-band floating operation order and quantization.
                units = round((2-ap)/(prime+1-ap)*math.log(prime)*10**12)
                restored.append({'prime': prime, 'scale_exponent': exponent,
                                 'reduced_A': a % prime, 'reduced_B': b % prime,
                                 'trace': ap, 'score_units': units})
        correction = sum(r['score_units'] for r in restored)
        result = {'id': row['id'], 'family': row['family'], 'band': row['band'],
                  'parameter': row['parameter'], 'removable_scales': scales,
                  'restored_terms': restored, 'correction_units': correction}
        rows.append(result)
        groups[(row['band'], row['family'])].append({
            **row, **score, 'audit_correction_units': correction,
            'audit_restored_good': len(restored)})
    comparisons = []
    for (band, family), pool in sorted(groups.items()):
        if len(pool) != 1024:
            raise ArithmeticError('fixed per-group scalar allocation differs')
        actual = [r['id'] for r in selected['selected'] if (r['band'], r['family']) == (band, family)]
        original = sorted(pool, key=lambda r: (-r['combined_late_units'], -r['combined_late_good'],
                                              r['denominator'], r['numerator']))[:6]
        if actual != [r['id'] for r in original]:
            raise ArithmeticError('frozen six-fibre ordering differs')
        diagnostic = sorted(pool, key=lambda r: (
            -r['combined_late_units']-r['audit_correction_units'],
            -r['combined_late_good']-r['audit_restored_good'], r['denominator'], r['numerator']))[:6]
        alternative = [r['id'] for r in diagnostic]
        comparisons.append({'band': band, 'family': family, 'scalar_candidates': len(pool),
                            'candidates_with_restored_good_terms': sum(r['audit_restored_good'] > 0 for r in pool),
                            'frozen_finalists': actual, 'diagnostic_six_with_restored_terms': alternative,
                            'shared_finalists': len(set(actual) & set(alternative)),
                            'correction_minimum_units': min(r['audit_correction_units'] for r in pool),
                            'correction_maximum_units': max(r['audit_correction_units'] for r in pool)})
    return {'schema': 'elliptic-curves.higher-mw16-omitted-good-primes.v1',
            'status': 'PASS_EXACT_LOCAL_REDUCTION_DIAGNOSTIC',
            'sources': {str(path.relative_to(ROOT)): cert.hashed(path) for path in paths},
            'models_checked': len(rows),
            'candidates_with_restored_good_terms': sum(bool(r['restored_terms']) for r in rows),
            'restored_terms_by_prime': dict(sorted(Counter(str(t['prime']) for r in rows for t in r['restored_terms']).items())),
            'distinct_reduced_models_checked': len(traces), 'groups': comparisons, 'rows': rows,
            'claim_boundary': 'Exact local p-scaling and restored good-prime traces on the frozen10240 '
                'narrow higher-MW16 scalar candidates. All five family resultants and complete33-pair '
                'classifications restrict removable scales in5..131071 to the listed5/13 cases. '
                'Traces use both character sums and exhaustive finite-point enumeration. '
                'Displayed-reduction scores omit these terms; their saved values remain correct for '
                'the stated displayed-model score. The diagnostic ordering reads no point outcomes '
                'or catalogue and does not alter frozen finalists, validation, point budgets or either '
                'live campaign. It does not recover initial-retention omissions, certify a better rank '
                'predictor or authorize a corrected-score parameter search. Primes2,3 are excluded.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = expected()
    if args.check:
        if cert.read(OUT) != json.loads(json.dumps(data)):
            raise ArithmeticError('local-reduction diagnostic replay differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve local-reduction diagnostic')
        checkpoint(OUT, data)
    print(data['status'], data['candidates_with_restored_good_terms'],
          data['restored_terms_by_prime'], flush=True)
