#!/usr/bin/env python3
"""Exact native M25 and productive generic anchors for curve200; no point search."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from research_runtime.store import checkpoint
from memory_rank_certificate import checked_rank
from future_point_admission import FinitePointAdmission
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity
import compact_mw16_specialization as spec

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
ART = ROOT/'artifacts/generated-results/elliptic-curves'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(folder):
    from sage.all import ZZ, matrix, pari
    folder.mkdir(exist_ok=False)
    paths = [ART/'strata60_mw16_results_v1.json', spec.ATLAS,
        LOCAL/'strata60-mw16-pari-v1/strata-a1-fibration-02-b2-s-1-top/proof-input.json',
        LOCAL/'prospective-mw16-h4096-v1/a1-fibration-02/generic-census.json']
    inputs = [json.loads(p.read_text()) for p in paths]
    proofs, atlas, initial, census = inputs
    proof = next(r for r in proofs['curves'] if r['id'] == 'strata-a1-fibration-02-b2-s-1-top')
    family = next(r for r in atlas['families'] if r['fibration_id'] == 'a1-fibration-02')
    protocol = {'family': 'a1-fibration-02', 'parameter': '-32999/14074', 'initial_rank': 25,
        'generic_rank': 16, 'prime_bound': 1000, 'cvp_node_limit': 2000000,
        'point_searches': 0,
        'rule': 'Freshly certify every gain in the retained attempt. Retain all productive '
                'generic parities. Solve only their generic CVPs, with independent rational '
                'replay; no full-shell completeness claim is needed for this subset.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'sources': {str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in (
            'prepare_curve200_productive.py', 'compact_mw16_specialization.py',
            'future_point_admission.py', 'memory_rank_certificate.py',
            'visibility_lattice_fast.py', 'visibility_lattice_v2.py')}}
    checkpoint(folder/'protocol.json', protocol)
    model = tuple(map(F, proof['curve']))
    basis = tuple(tuple(map(F, p)) for p in proof['points'])
    original, generic = spec.specialize(family, '-32999/14074')
    scale = F(initial['family_to_curve_scale_u'])
    if (len(basis) != 25 or proof['parameter'] != '-32999/14074'
            or tuple(v/scale**e for v, e in zip(original, (1,2,3,4,6))) != model
            or tuple((x/scale**2, y/scale**3) for x, y in generic) != basis[:16]
            or initial['generic_points'] != proof['points'][:16]):
        raise ArithmeticError('exact native model/generic basis identity failed')
    old = proof['rank_certificate']
    fresh = checked_rank(model, basis, [s['prime'] for s in old['signatures']],
                         old['no_rational_2_torsion_prime'])
    if json.loads(json.dumps(fresh)) != old:
        raise ArithmeticError('existing M25 certificate differs')
    checkpoint(folder/'seed-M25.json', {'curve': list(map(str, model)),
        'points': [list(map(str, p)) for p in basis], 'generic_rank': 16,
        'rank_lower_bound': 25, 'proof': fresh})
    gains = []
    for name, data, start in (('initial', initial, basis[:16]),):
        if tuple(map(F, data['curve'])) != model:
            raise ArithmeticError('historical attempt uses another curve')
        admission = FinitePointAdmission(model, start)
        for i, chart in enumerate(data['charts']):
            if chart['index'] != i:
                raise ArithmeticError('historical chart order differs')
            before = len(admission.points)
            for p in chart['search']['finite_curve_points']:
                admission.consider((F(p['x']), F(p['y'])))
            if len(admission.points) > before:
                certificate = checked_rank(model, admission.points, admission.primes,
                    old['no_rational_2_torsion_prime'])
                mask = sum((int(w) % 2) << j for j, w in enumerate(chart['centre']['representative'][:16]))
                gains.append({'arm': name, 'chart_index': i, 'mask': mask, 'before': before,
                    'after': len(admission.points), 'proof': certificate,
                    'points': [list(map(str, p)) for p in admission.points]})
        if len(admission.points) != 25:
            raise ArithmeticError('retained whole-cloud bound differs')
    checkpoint(folder/'gains.json', gains)
    gram = [[F(x) for x in row] for row in family['generic_height_gram']]
    if gram != [[F(x) for x in row] for row in census['gram']]:
        raise ArithmeticError('generic Gram differs from census coordinate frame')
    if any((2*x).denominator != 1 for row in gram for x in row):
        raise ArithmeticError('unexpected generic height denominator')
    g = matrix(ZZ, [[int(2*x) for x in row] for row in gram])
    u = matrix(ZZ, pari(g).qflllgram()).transpose()
    if abs(u.det()) != 1:
        raise ArithmeticError('nonunimodular generic LLL')
    inv = u.inverse()
    fast, reference = IntegerExactParity((u*g*u.transpose()).rows()), ExactParity((u*g*u.transpose()).rows())
    masks = sorted({row['mask'] for row in gains})
    seeds = {row['mask']: row for row in census['records']}
    rows, checks = [], []
    for mask in masks:
        w = matrix(ZZ, 1, 16, seeds[mask]['representative'])
        reduced_seed = tuple(map(int, (w*inv).row(0)))
        residue = tuple(x % 2 for x in reduced_seed)
        cert = fast.solve(residue, reduced_seed, protocol['cvp_node_limit'])
        if cert != reference.solve(residue, reduced_seed, protocol['cvp_node_limit']):
            raise ArithmeticError('independent generic CVP replay differs')
        words = []
        for v in cert['minima']:
            z = tuple(map(int, (matrix(ZZ, 1, 16, v)*u).row(0)))
            if next(x for x in z if x) < 0:
                z = tuple(-x for x in z)
            if sum((x % 2) << j for j, x in enumerate(z)) != mask:
                raise ArithmeticError('generic parity transport differs')
            words.append(z)
        rows.append({'mask': mask, 'word': list(min(words)), 'norm': cert['norm']})
        checks.append({'mask': mask, 'reduced_seed': list(reduced_seed), 'proof': cert})
    checkpoint(folder/'generic-cvp-proofs.json', {'gram': [list(map(int, r)) for r in g.rows()],
        'LLL': [list(map(int, r)) for r in u.rows()], 'checks': checks})
    bank = {'status': 'COMPLETE_FROZEN_PRODUCTIVE_SUBSET', 'dimension': 16,
        'generic_gram_scale': 2, 'rows': rows, 'shells': sorted({r['norm'] for r in rows}),
        'protocol_sha256': sha(folder/'protocol.json'), 'gains_sha256': sha(folder/'gains.json'),
        'generic_cvp_sha256': sha(folder/'generic-cvp-proofs.json'),
        'claim_boundary': 'Exact historically productive subset, not the complete generic shell bank.'}
    checkpoint(folder/'anchor-bank.json', bank)
    for category in ('inputs', 'sources'):
        if any(sha(ROOT/p) != h for p, h in protocol[category].items()):
            raise ArithmeticError('preparation binding changed')
    checkpoint(folder/'prepared.json', {'status': 'PASS_NATIVE_M25_AND_PRODUCTIVE_ANCHORS',
        'point_searches': 0, 'anchor_count': len(rows),
        'files': {p.name: sha(p) for p in sorted(folder.glob('*.json'))}})
    print('CURVE200_PREPARED', len(rows), masks, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder', type=Path, required=True)
    run(p.parse_args().folder)
