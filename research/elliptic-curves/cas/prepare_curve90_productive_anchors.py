#!/usr/bin/env python3
"""Freeze every historically productive curve90 generic anchor, with fresh proofs.

Retrospective discovery design, not a blinded selection. Replays retained point
clouds without rerunning searches and exports no unobserved-rank claim.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from future_point_admission import FinitePointAdmission
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
ART = ROOT/'artifacts/generated-results/elliptic-curves'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(folder):
    folder.mkdir(exist_ok=False)
    visibility_path = ART/'mw16_rank27_visibility_input_v1.json'
    inputs = [visibility_path,
        LOCAL/'extended20-mw16-pari-v1/a1-fibration-01-052/result.json',
        LOCAL/'mw16-new26-a1-fibration-01-052-adaptive-v1/result.json',
        LOCAL/'curve90-upper-shell-audit-v1/parent-bank.json']
    visibility, initial, adaptive, bank = [json.loads(p.read_text()) for p in inputs]
    for p in inputs[1:3]:
        if sha(p) != visibility['original_sources'][str(p.relative_to(ROOT))]:
            raise ArithmeticError('historical discovery input differs from retained visibility proof')
    protocol = {'point_searches': 0, 'prime_bound': 1000,
        'rule': 'Replay both complete historical attempts with their actual ordered initial '
                'subgroups. Retain every generic parity from a chart with a freshly certified '
                'finite-rank gain. No per-chart tuning, refill or new point search.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        'sources': {str(p.relative_to(ROOT)): sha(p) for p in (
            Path(__file__), CAS/'future_point_admission.py', CAS/'memory_rank_certificate.py')}}
    checkpoint(folder/'protocol.json', protocol)
    seed = visibility['point_proof']
    model = tuple(map(F, seed['discovery_curve']))
    no_two = seed['rank_certificate']['no_rational_2_torsion_prime']
    gains = []
    for name, data, raw_basis in (('initial', initial, visibility['generic_points']),
                                 ('adaptive', adaptive, visibility['initial_points'])):
        basis = tuple(tuple(map(F, p)) for p in raw_basis)
        admission = FinitePointAdmission(model, basis, prime_bound=1000)
        for index, chart in enumerate(data['charts']):
            if chart['index'] != index:
                raise ArithmeticError('historical chart ordering differs')
            before = len(admission.points)
            for p in chart['search']['finite_curve_points']:
                admission.consider((F(p['x']), F(p['y'])))
            if len(admission.points) > before:
                proof = checked_rank(model, admission.points, admission.primes, no_two)
                word = chart['centre']['representative']
                mask = sum((int(w) % 2) << j for j, w in enumerate(word[:16]))
                gains.append({'arm': name, 'chart_index': index, 'generic_mask': mask,
                    'before': before, 'after': len(admission.points),
                    'points': [list(map(str, p)) for p in admission.points], 'proof': proof})
                checkpoint(folder/'gains.json', gains)
        expected = 26 if name == 'initial' else 27
        if len(admission.points) != expected:
            raise ArithmeticError('retained complete-cloud lower bound differs')
    masks = sorted({r['generic_mask'] for r in gains})
    by_mask = {r['mask']: r for r in bank['rows']}
    if any(m not in by_mask for m in masks):
        raise ArithmeticError('productive anchor missing from exact upper-shell bank')
    result = {'status': 'COMPLETE_FROZEN_PRODUCTIVE_SUBSET', 'dimension': 16,
        'generic_gram_scale': 2, 'rows': [by_mask[m] for m in masks],
        'shells': sorted({by_mask[m]['norm'] for m in masks}), 'productive_charts': len(gains),
        'point_searches': 0, 'protocol_sha256': sha(folder/'protocol.json'),
        'gains_sha256': sha(folder/'gains.json'),
        'claim_boundary': 'All freshly rechecked historical gains define a finite anchor subset; '
                          'not the complete generic shell bank, a new rank or a success predictor.'}
    for category in ('inputs', 'sources'):
        if any(sha(ROOT/p) != h for p, h in protocol[category].items()):
            raise ArithmeticError('productive-bank input changed')
    checkpoint(folder/'anchor-bank.json', result)
    print('PRODUCTIVE_ANCHORS', masks, len(gains), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True, type=Path)
    run(parser.parse_args().folder)
