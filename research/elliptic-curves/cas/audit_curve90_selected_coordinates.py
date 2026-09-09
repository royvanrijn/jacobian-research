#!/usr/bin/env python3
"""Post-seal pointwise comparison of two mapping policies; no point search.

This is retrospective calibration. Neither point coordinates nor map-visibility
labels feed back into the frozen centre selector. Reference CVP status is kept
separate from this exact rational map/point-visibility audit.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

from v3_warm_engine import load, certified_state
from pointed_quartic_search import PointedQuarticSearch
from search_observability import point_visibility
from research_runtime.store import checkpoint
import pari_pointed_backend as backend
from run_future_visibility_upper_audit import guard

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
PREP = ROOT/'artifacts/local/elliptic-curves/curve90-v3-preparation-v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def run(folder, output):
    if output.exists():
        raise FileExistsError('preserve coordinate audit')
    seal, protocol = read(folder/'sealed.json'), read(folder/'protocol.json')
    guard(protocol)
    if (seal['status'] != 'SEALED_MASKED_M26_SELECTION'
            or seal['selection_sha256'] != sha(folder/'landscape/selection.json')
            or seal['protocol_sha256'] != sha(folder/'protocol.json')):
        raise ArithmeticError('selection not sealed')
    seed, point_record = read(PREP/'seed-M26.json'), read(PREP/'withheld-point.json')
    prepared = read(PREP/'prepared.json')
    if any(sha(PREP/n) != prepared['files'][n] for n in ('seed-M26.json', 'withheld-point.json')):
        raise ArithmeticError('seed or evaluation point changed')
    selection = read(folder/'landscape/selection.json')
    if selection['basis'] != seed['points']:
        raise ArithmeticError('selected seed differs')
    model = tuple(map(F, seed['curve']))
    basis = tuple(tuple(map(F, p)) for p in seed['points'])
    point = tuple(map(F, point_record['point']))
    state = certified_state(model, basis, seed['proof'])
    policies = {'factor_free': 'factor_free_pari_mapping.sage',
                'quartic_minimized': 'prepare_extended20_mw16_pari_batch.sage'}
    modules = {key: load('coordinate_audit_'+key, CAS/name) for key, name in policies.items()}
    for module in modules.values():
        module.pari.allocatemem(256000000, silent=True)
    result = {'status': 'RUNNING', 'point_searches': 0,
        'selection_sha256': seal['selection_sha256'],
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in (
            folder/'sealed.json', folder/'protocol.json', folder/'landscape/selection.json',
            PREP/'seed-M26.json', PREP/'withheld-point.json')},
        'sources': {**backend.sources(), **{str(p.relative_to(ROOT)): sha(p) for p in (
            Path(__file__), *(CAS/n for n in policies.values()))}},
        'height': 125000, 'mapping_seconds': {key: 0.0 for key in policies}, 'rows': []}
    checkpoint(output, result)
    for index, centre in enumerate(selection['centres']):
        for key, module in modules.items():
            start = time.monotonic()
            mapping = module.mapping(model, basis, centre)
            result['mapping_seconds'][key] += time.monotonic()-start
            search = PointedQuarticSearch(state=state,
                centre={'coefficients': centre['representative']},
                coordinate_policy=mapping['coordinate_policy'])
            backend.validate_map(search, mapping)
            for sign in (1, -1):
                observed = point_visibility(search.chart_record(), (point[0], sign*point[1]))
                height = observed.get('minimum_affine_height')
                result['rows'].append({'chart_index': index, 'policy': key, 'sign': sign,
                    'height': height, 'coordinate': observed.get('coordinate'),
                    'at_parameter_infinity': observed.get('at_parameter_infinity', False),
                    'visible': ((height is not None and height <= 125000)
                                or observed.get('at_parameter_infinity', False))})
        if (index+1) % 100 == 0:
            checkpoint(output, result)
            print('COORDINATE_AUDIT', index+1, '/', len(selection['centres']), flush=True)
    guard(protocol)
    for category in ('inputs', 'sources'):
        if any(sha(ROOT/p) != h for p, h in result[category].items()):
            raise ArithmeticError('coordinate audit binding changed')
    result['status'] = 'PASS_EXACT_SELECTED_COORDINATE_AUDIT'
    result['visible_counts'] = {key: sum(r['visible'] for r in result['rows'] if r['policy'] == key)
                                for key in policies}
    result['claim_boundary'] = ('Exact pointwise geometry after centre selection. '
        'Reference CVP replay of the expanded selection remains separate. No new point/rank result.')
    checkpoint(output, result)
    print('COORDINATE_RESULT', result['visible_counts'], result['mapping_seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.folder, args.output)
