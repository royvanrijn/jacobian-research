#!/usr/bin/env python3
"""Describe two frozen higher-parameter cohorts without changing selection."""
import argparse
import ast
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / 'elliptic-curves/cas'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'higher_mw16_population_comparison_v1.json'


def summary(values):
    values = sorted(values)
    if not values:
        raise ArithmeticError('empty comparison group')
    n = len(values)
    return {'count': n, 'minimum': values[0], 'maximum': values[-1],
            'mean': str(Fraction(sum(values), n)),
            'median': str(Fraction(values[(n-1)//2] + values[n//2], 2))}


def functions(path):
    return {n.name: ast.dump(n, include_attributes=False)
            for n in ast.parse(path.read_text()).body
            if isinstance(n, ast.FunctionDef)}


def load_arm(prefix):
    score = LOCAL / (prefix + '-mw16-higher-scores-v1')
    scan = LOCAL / (prefix + '-mw16-higher-annuli-v1')
    paths = [score / 'protocol.json', score / 'result.json',
             score / 'controller/ledger.json', score / 'fresh-validation.json',
             scan / 'protocol.json', scan / 'replay.json',
             ART / (prefix + '_mw16_higher_selection_v1.json')]
    p, scores, ledger, validation, population, replay, selection = map(cert.read, paths)
    if (ledger['status'] != 'PASS' or replay['status'] != 'PASS'
            or scores['status'] != 'COMPLETE_FROZEN10240'
            or selection['status'] != 'PASS_FROZEN60_SELECTION'
            or validation['status'] != 'PASS'):
        raise ArithmeticError('both complete scalar and disjoint validation gates required')
    if any(cert.hashed(ROOT/n) != h for n, h in p['sources'].items()):
        raise ArithmeticError('frozen scalar sources changed')
    if (scores['protocol_sha256'] != cert.hashed(paths[0])
            or selection['protocol_sha256'] != cert.hashed(paths[0])
            or selection['scores_sha256'] != cert.hashed(paths[1])
            or validation['selection_sha256'] != cert.hashed(paths[-1])
            or validation['protocol_sha256'] != cert.hashed(paths[0])):
        raise ArithmeticError('score/selection/validation binding differs')
    if (len(p['rows']) != 10240 or len(scores['rows']) != 10240
            or len(selection['selected']) != 60 or len(validation['rows']) != 60):
        raise ArithmeticError('fixed later allocations differ')
    initial = defaultdict(list)
    final = defaultdict(list)
    full = {}
    for row, score_row in zip(p['rows'], scores['rows']):
        if row['id'] != score_row['id'] or row['id'] in full:
            raise ArithmeticError('scalar row ordering or uniqueness differs')
        full[row['id']] = {**row, **score_row}
        initial[(row['band'], row['family'])].append(row)
    seen = set()
    for row, check in zip(selection['selected'], validation['rows']):
        if (row['id'] != check['id'] or row['id'] in seen
                or row != full[row['id']] or check['status'] != 'PASS'):
            raise ArithmeticError('frozen finalist/validation correspondence differs')
        seen.add(row['id'])
        final[(row['band'], row['family'])].append({**row, **check})
    records = {}
    for group, rows in sorted(initial.items()):
        finalists = final[group]
        if len(rows) != 1024 or len(finalists) != 6:
            raise ArithmeticError('fixed per-group allocation differs')
        records[group] = {
            'initial_selected_score_units': summary([r['combined_selection_units'] for r in rows]),
            'final_selection_score_units': summary([r['combined_late_units'] for r in finalists]),
            'disjoint_validation_score_units': summary([r['validation_units'] for r in finalists]),
            'finalist_parameter_heights': summary([max(abs(r['numerator']), r['denominator']) for r in finalists]),
            'finalists': [{'id': r['id'], 'parameter': r['parameter']} for r in finalists]}
    if len(records) != 10:
        raise ArithmeticError('ten band/family groups required')
    return p, population, records, paths


def expected():
    narrow_source = CAS / 'score_joint_mw16_higher.py'
    broad_source = CAS / 'score_broad_mw16_higher.py'
    narrow_functions, broad_functions = map(functions, (narrow_source, broad_source))
    same = ['trace', 'selected_rows', 'fresh_program']
    if any(narrow_functions[n] != broad_functions[n] for n in same):
        raise ArithmeticError('scalar score or finalist selector implementation changed')
    arms = [load_arm(prefix) for prefix in ('joint', 'broad')]
    np, ns, nr, npaths = arms[0]
    bp, bs, br, bpaths = arms[1]
    fields = ['ordering', 'fresh_validation_primes', 'selected_curves', 'seconds_per_curve',
              'score_wall_seconds', 'validation_wall_seconds', 'maximum_workers',
              'checkpoint_block', 'cost_gate_rows', 'cost_gate_maximum_projected_serial_seconds']
    if any(np[k] != bp[k] for k in fields):
        raise ArithmeticError('selection ordering, prime band or later budget differs')
    def address_frame(row):
        return tuple(row[k] for k in ('band', 'family', 'sign', 'shards', 'shard'))
    narrow_frames = {address_frame(r) for r in ns['rows']}
    broad_frames = {address_frame(r) for r in bs['rows']}
    if (len(narrow_frames) != 20 or len(broad_frames) != 320
            or narrow_frames & broad_frames):
        raise ArithmeticError('disjoint fixed parameter slices required')
    groups = []
    for band, family in sorted(nr):
        n, b = nr[(band, family)], br[(band, family)]
        deltas = {}
        for metric in ('initial_selected_score_units', 'final_selection_score_units',
                       'disjoint_validation_score_units'):
            deltas[metric] = {k: str(Fraction(b[metric][k]) - Fraction(n[metric][k]))
                              for k in ('minimum', 'mean', 'median', 'maximum')}
        groups.append({'band': band, 'family': family, 'narrow': n, 'broad': b,
                       'broad_minus_narrow': deltas})
    paths = [Path(__file__).resolve(), Path(cert.__file__), narrow_source, broad_source,
             *npaths, *bpaths]
    return {'schema': 'elliptic-curves.higher-mw16-population-comparison.v1',
            'status': 'PASS_DESCRIPTIVE_FIXED_COHORT_COMPARISON',
            'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths},
            'identical_function_asts': same, 'identical_protocol_fields': fields,
            'primitive_addresses': {'narrow': ns['total_primitive_population'],
                                    'broad': bs['total_primitive_population']},
            'groups': groups,
            'claim_boundary': 'Descriptive comparison of two disjoint finite parameter populations. '
                'The broad trial changes slice coverage while holding later selection rules and budgets fixed. '
                'These are not nested populations or replicated randomized experiments. Score improvement '
                'is not a new rank, proof of optimal selection, or calibrated prediction of discovery yield. '
                'Wholly disjoint prime scores were recorded only after finalists were frozen and cannot '
                'change selection or point budgets. No point outcome, record label or catalogue is read.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = expected()
    if args.check:
        if cert.read(OUT) != data:
            raise ArithmeticError('population comparison replay differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve population comparison')
        checkpoint(OUT, data)
    print(data['status'], flush=True)
