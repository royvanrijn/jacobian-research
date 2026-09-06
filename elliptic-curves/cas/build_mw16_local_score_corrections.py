#!/usr/bin/env python3
"""Exact local correction tables on P1(Z/125) and P1(Z/169)."""
import argparse
import math
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
import audit_mw16_postscale_reduction as reduction
import audit_higher_mw16_omitted_good_primes as diagnostic
from research_runtime.store import checkpoint

ROOT = reduction.ROOT
OUT = reduction.first.support.ART/'mw16_local_score_corrections_v1.json'


def index(n, d, prime, modulus):
    if d % prime:
        return n*pow(d, -1, modulus) % modulus
    if n % prime == 0:
        raise ArithmeticError('nonprimitive local pair')
    s = d*pow(n, -1, modulus) % modulus
    if s % prime:
        raise ArithmeticError('infinity coordinate outside pZ')
    return modulus+s//prime


def coefficients(f, n, d):
    return [sum(int(c)*n**i*d**(weight-i) for i, c in enumerate(f[key]))
            for key, weight in [('A_coefficients_low_to_high', 8), ('B_coefficients_low_to_high', 12)]]


def expected():
    first = cert.read(reduction.first.OUT)
    after = cert.read(reduction.OUT)
    replay = reduction.first.support.ART/'mw16_postscale_reduction_replay_v1.json'
    if (cert.read(replay)['status'] != 'PASS'
            or cert.read(reduction.D/'ledger.json')['status'] != 'PASS'
            or after['status'] != 'PASS_SINGLE_SCALE_CLASSIFICATION'
            or cert.read(reduction.first.D/'ledger.json')['status'] != 'PASS'):
        raise ArithmeticError('complete independent family-wide reduction proofs required')
    for artifact in (first, after):
        if any(cert.hashed(ROOT/n) != h for n, h in artifact['sources'].items()):
            raise ArithmeticError('family-wide proof source changed')
    tables = []
    traces = {}
    for f in cert.read(spec.ATLAS)['families']:
        for prime, modulus in [(5, 125), (13, 169)]:
            entries = []
            for i in range(modulus+modulus//prime):
                chart, coordinate = ('affine', i) if i < modulus else ('infinity', prime*(i-modulus))
                n, d = (coordinate, 1) if chart == 'affine' else (1, coordinate)
                if index(n, d, prime, modulus) != i:
                    raise ArithmeticError('complete projective-ring frame differs')
                a, b = coefficients(f, n, d)
                exponent = 0
                while a % prime**4 == 0 and b % prime**6 == 0:
                    a //= prime**4
                    b //= prime**6
                    exponent += 1
                if exponent > 1:
                    raise ArithmeticError('representative contradicts second-scale exclusion')
                # Match the complete family-wide polynomial proof, not only this representative.
                matches = [r for r in after['rows']
                           if (r['family'], r['prime'], r['chart']) == (f['fibration_id'], prime, chart)
                           and (coordinate-r['first_residue']) % r['first_modulus'] == 0]
                if len(matches) != exponent:
                    raise ArithmeticError('representative scale differs from universal ball proof')
                ap = None
                if matches:
                    r = matches[0]
                    z = ((coordinate-r['first_residue'])//r['first_modulus']) % prime
                    cell = r['cells'][z]
                    if (a % prime, b % prime) != (cell['reduced_A'], cell['reduced_B']):
                        raise ArithmeticError('local coefficient table differs from exact divided polynomial')
                    key = (a % prime, b % prime, prime)
                    if key not in traces:
                        traces[key] = diagnostic.trace(*key)
                    ap = traces[key]
                    if (ap is not None) != cell['good_reduction_after_one_scale']:
                        raise ArithmeticError('local good-reduction classification differs')
                units = 0 if ap is None else round((2-ap)/(prime+1-ap)*math.log(prime)*10**12)
                entries.append({'index': i, 'restored_good': ap is not None,
                                'restored_trace': ap, 'correction_units': units})
            tables.append({'family': f['fibration_id'], 'prime': prime, 'modulus': modulus,
                           'entries': entries})
    # New full-height parameters check the lookup's unit-coordinate normalization.
    p = cert.read(diagnostic.D/'protocol.json')
    measured = cert.read(diagnostic.OUT)
    if measured['status'] != 'PASS_EXACT_LOCAL_REDUCTION_DIAGNOSTIC' or len(measured['rows']) != 10240:
        raise ArithmeticError('complete scalar-model diagnostic required')
    if any(cert.hashed(ROOT/n) != h for n, h in measured['sources'].items()):
        raise ArithmeticError('frozen diagnostic source changed')
    table_index = {(t['family'], t['prime']): t for t in tables}
    restored_terms = restored_candidates = 0
    for row, proof in zip(p['rows'], measured['rows']):
        if row['id'] != proof['id']:
            raise ArithmeticError('scalar-model correspondence differs')
        actual = []
        for prime in (5, 13):
            t = table_index[(row['family'], prime)]
            entry = t['entries'][index(row['numerator'], row['denominator'], prime, t['modulus'])]
            if entry['restored_good']:
                actual.append((prime, entry['restored_trace'], entry['correction_units']))
        expected_terms = [(r['prime'], r['trace'], r['score_units']) for r in proof['restored_terms']]
        if actual != expected_terms or sum(r[2] for r in actual) != proof['correction_units']:
            raise ArithmeticError('local projective lookup differs from exact full-height scaling')
        restored_terms += len(actual)
        restored_candidates += bool(actual)
    if len(p['rows']) != 10240:
        raise ArithmeticError('all10240 scalar models required')
    paths = [Path(__file__).resolve(), Path(diagnostic.__file__), spec.ATLAS,
             reduction.first.OUT, reduction.first.D/'ledger.json', reduction.OUT,
             reduction.D/'ledger.json', replay, diagnostic.OUT, diagnostic.D/'protocol.json']
    return {'schema': 'elliptic-curves.mw16-local-score-corrections.v1', 'status': 'PASS',
            'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths},
            'tables': tables, 'projective_ring_entries': sum(len(t['entries']) for t in tables),
            'full_height_models_checked': len(p['rows']), 'restored_candidates': restored_candidates,
            'restored_terms': restored_terms,
            'claim_boundary': 'Exact reusable additive good-prime corrections on the complete local '
                'projective rings modulo125 and169 for all five MW16 families. Universal first- and '
                'second-scale proofs justify constancy on each local cell, including infinity; '
                'all10240 frozen higher scalar models independently reproduce the lookup corrections. '
                'This does not change the running scans or prove selection optimality or a new rank. '
                'It supplies corrected scoring data for a separately declared future campaign.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = expected()
    if args.check:
        if cert.read(OUT) != data:
            raise ArithmeticError('local correction table replay differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve local score tables')
        checkpoint(OUT, data)
    print('MW16 LOCAL CORRECTION TABLES PASS', data['projective_ring_entries'],
          data['full_height_models_checked'], data['restored_terms'], flush=True)
