#!/usr/bin/env python3
"""Finite second-scale exclusion and reduction cells on the ten admitted balls."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import classify_mw16_score_prime_scalings as first
import classify_r17_other_small_prime_scalings as tree
from research_runtime.store import checkpoint

ROOT = first.ROOT
D = ROOT/'artifacts/local/elliptic-curves/mw16-postscale-reduction-v1'
OUT = first.support.ART/'mw16_postscale_reduction_v1.json'


def sources():
    paths = [Path(__file__).resolve(), Path(tree.__file__), first.OUT, first.D/'ledger.json']
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}


def balls():
    return [{'family': r['family'], 'prime': r['prime'], 'chart': c['chart'], **b}
            for r in cert.read(first.OUT)['rows'] for c in r['charts'] for b in c['scale_balls']]


def prepare():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve finite postscale audit')
    if cert.read(first.D/'ledger.json')['status'] != 'PASS' or len(balls()) != 10:
        raise ArithmeticError('all ten independently verified first-scale balls required')
    checkpoint(D/'protocol.json', {'sources': sources(), 'balls': balls(),
        'maximum_depth': 6, 'maximum_live_residues_per_chart': 4096,
        'maximum_candidates_per_level': 200000, 'wall_seconds': 300,
        'rss_bytes': 2147483648,
        'scope': 'Exactly ten already certified p-adic first-scale balls, at primes5 or13. '
                 'Test all next residues for good displayed reduction and run one bounded six-level '
                 'coefficient tree on each divided polynomial to decide whether a second scale is '
                 'possible anywhere in that ball. No recursive third-scale analysis, parameter search, '
                 'point search or live score alteration. Any second-scale ball or resource cap remains '
                 'unresolved for the one-scale minimality claim.'})


def expected():
    p = cert.read(D/'protocol.json')
    if p['sources'] != sources() or p['balls'] != balls():
        raise ArithmeticError('frozen postscale scope changed')
    rows = []
    for ball in p['balls']:
        q = ball['prime']
        a = list(map(int, ball['A_divided_coefficients']))
        b = list(map(int, ball['B_divided_coefficients']))
        classification = tree.chart(a, b, q, False, p)
        no_second = (classification['status'] == 'COMPLETE_RESIDUE_CLASSIFICATION'
                     and not classification['scale_balls'])
        cells = []
        for z in range(q):
            aa = sum(c*pow(z, i, q) for i, c in enumerate(a)) % q
            bb = sum(c*pow(z, i, q) for i, c in enumerate(b)) % q
            good = (4*aa**3+27*bb**2) % q != 0
            cells.append({'next_digit': z, 'reduced_A': aa, 'reduced_B': bb,
                          'good_reduction_after_one_scale': good,
                          'parameter_residue': ball['residue']+ball['modulus']*z,
                          'parameter_modulus': ball['modulus']*q})
        rows.append({'family': ball['family'], 'prime': q, 'chart': ball['chart'],
                     'first_residue': ball['residue'], 'first_modulus': ball['modulus'],
                     'no_second_scale': no_second, 'second_scale_tree': classification,
                     'cells': cells})
    return {'schema': 'elliptic-curves.mw16-postscale-reduction.v1',
            'status': 'PASS_SINGLE_SCALE_CLASSIFICATION' if all(r['no_second_scale'] for r in rows) else 'UNKNOWN_SECOND_SCALE',
            'sources': sources(), 'protocol_sha256': cert.hashed(D/'protocol.json'),
            'rows': rows, 'claim_boundary': p['scope']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'build', 'check'])
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare()
    else:
        data = expected()
        if args.stage == 'check':
            if cert.read(OUT) != json.loads(json.dumps(data)):
                raise ArithmeticError('postscale reduction replay differs')
        else:
            if OUT.exists():
                raise FileExistsError('preserve postscale proof')
            checkpoint(OUT, data)
        print(data['status'], [(r['family'], r['prime'], r['no_second_scale']) for r in data['rows']], flush=True)
