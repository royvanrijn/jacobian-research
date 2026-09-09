#!/usr/bin/env python3
"""Zero-search upper-shell MW16 calibration; initial norm8/10 audit is preserved."""
import argparse
import csv
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from v3_warm_engine import load, certified_state
from research_runtime.store import checkpoint
from visibility_future_upper import landscape

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
PREP = LOCAL/'curve90-v3-preparation-v1'
CONTROL = LOCAL/'curve302-seed-universality-panel-v1/residual-strict-03'
UPPER = LOCAL/'curve90-upper-shell-audit-v1'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    return json.loads(p.read_text())


def sources(old):
    names = ('visibility_future_upper.py', 'visibility_lattice_fast.py', 'v3_warm_engine.py',
             'v3_warm_support.py', 'run_future_visibility_upper_audit.py')
    return {**old.sources(), **{str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in names}}


def guard(protocol):
    for category in ('inputs', 'sources'):
        if any(sha(ROOT/p) != h for p, h in protocol[category].items()):
            raise ArithmeticError('frozen audit binding changed')


def run(mode, folder):
    folder.mkdir(exist_ok=False)
    old = load('audit_reference_v3', CAS/'adaptive_visibility_cascade_v3.sage')
    base_policy = read(CONTROL/'protocol.json')
    policy = {k: base_policy[k] for k in ('anchors_per_shell', 'canonical_per_shell',
                                         'exact_cvp_node_limit', 'height')}
    if mode == 'calibrate302':
        seed, proof = read(CONTROL/'seed-input.json'), read(CONTROL/'seed-proof.json')
        rows = []
        with old.ORBITS.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                q = int(row['minimum_norm'])
                if q in (8, 10):
                    rows.append({'mask': int(row['orbit_mask']), 'norm': 2*q,
                                 'word': list(map(int, row['parent_MW17_w'].split()))})
        bank = {'dimension': 17, 'generic_gram_scale': 2, 'shells': [16, 20],
                'status': 'COMPLETE_EXACT_GENERIC_SHELL_BANK', 'rows': rows}
        # This is an adapter view of the retained 302 bank, not a new enumeration.
        inputs = [CONTROL/'seed-input.json', CONTROL/'seed-proof.json', old.ORBITS,
                  CONTROL/'replay-M17/epoch-00/selection.json', CONTROL/'protocol.json']
        policy['generic_rank'] = 17
    else:
        prepared = read(PREP/'prepared.json')
        if prepared['status'] != 'PASS_EXACT_MW16_BANK_AND_SEEDS':
            raise ArithmeticError('exact MW16 preparation required')
        if any(sha(PREP/n) != h for n, h in prepared['files'].items()):
            raise ArithmeticError('prepared input changed')
        seed = read(PREP/'seed-M26.json')
        proof = seed['proof']
        guard(read(UPPER/'protocol.json'))
        bank = read(UPPER/'parent-bank.json')
        if not bank['all_nonzero_parities_represented']:
            raise ArithmeticError('full upper-shell parity census required')
        # No read of the withheld point, M27 packet or historical chart data.
        inputs = [PREP/'seed-M26.json', UPPER/'parent-bank.json', UPPER/'protocol.json', PREP/'prepared.json',
                  CONTROL/'protocol.json']
        policy['generic_rank'] = 16
    policy['scaled_shells'] = bank['shells']
    protocol = {'mode': mode, 'policy': policy, 'point_searches': 0,
        'sources': sources(old), 'inputs': {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        'oracle_boundary': 'Selection uses only its frozen certified seed and generic bank; '
                           'withheld-point evaluation is a separate post-seal command.'}
    checkpoint(folder/'protocol.json', protocol)
    model = tuple(map(F, seed['curve']))
    basis = tuple(tuple(map(F, p)) for p in seed['points'])
    certified_state(model, basis, proof)
    result = landscape(model, basis, set(), folder/'landscape', policy, bank)
    if mode == 'calibrate302':
        comparison = dict(result)
        comparison.pop('generic_rank')
        for anchor in comparison['anchors']:
            anchor.pop('maps_sha256')
        expected = read(CONTROL/'replay-M17/epoch-00/selection.json')
        if json.loads(json.dumps(comparison)) != expected:
            raise ArithmeticError('future selector differs from retained 302 M18 selection')
    guard(protocol)
    checkpoint(folder/'sealed.json', {'status': 'PASS_RETAINED_302_SELECTION_EQUIVALENCE'
        if mode == 'calibrate302' else 'SEALED_MASKED_M26_SELECTION',
        'protocol_sha256': sha(folder/'protocol.json'),
        'selection_sha256': sha(folder/'landscape/selection.json'),
        'initial_rank': len(basis), 'point_searches': 0, 'centres': len(result['centres'])})
    print('SEALED', mode, len(result['centres']), flush=True)


def evaluate(folder):
    import numpy as np
    from sage.all import ZZ, matrix
    from visibility_lattice_v2 import ExactParity
    from pointed_quartic_search import PointedQuarticSearch
    from search_observability import point_visibility
    import pari_pointed_backend as backend
    if (folder/'evaluation.json').exists():
        raise FileExistsError('preserve masked evaluation')
    seal, protocol = read(folder/'sealed.json'), read(folder/'protocol.json')
    if (seal['status'] != 'SEALED_MASKED_M26_SELECTION'
            or sha(folder/'protocol.json') != seal['protocol_sha256']
            or sha(folder/'landscape/selection.json') != seal['selection_sha256']):
        raise ArithmeticError('masked selection must be sealed first')
    guard(protocol)
    data = read(folder/'landscape/selection.json')
    seed = read(PREP/'seed-M26.json')
    model = tuple(map(F, seed['curve']))
    basis = tuple(tuple(map(F, p)) for p in seed['points'])
    if data['basis'] != seed['points']:
        raise ArithmeticError('masked basis differs')
    state = certified_state(model, basis, seed['proof'])
    g, u = (matrix(ZZ, data[k]) for k in ('rounded_gram', 'LLL'))
    exact = ExactParity((u*g*u.transpose()).rows())
    for ai, anchor in enumerate(data['anchors']):
        path = folder/'landscape'/f'anchor-{ai:02d}-full.npz'
        if sha(path) != anchor['full_scores_sha256']:
            raise ArithmeticError('retained scores changed')
        with np.load(path, allow_pickle=False) as arrays:
            for row in anchor['refined']:
                j = row['extension']
                proof = exact.solve(arrays['reduced_residues'][j], arrays['babai_words'][j],
                                    protocol['policy']['exact_cvp_node_limit'])
                if json.loads(json.dumps(proof)) != row['cvp']:
                    raise ArithmeticError('independent rational CVP replay differs')
        checkpoint(folder/'evaluation-progress.json', {'phase': 'REFERENCE_CVP_REPLAY', 'anchors': ai+1})
        print('REFERENCE_CVP_REPLAY', ai+1, flush=True)
    # Only now open the withheld point, after selection and reference CVP checks.
    withheld_path = PREP/'withheld-point.json'
    prepared = read(PREP/'prepared.json')
    if sha(withheld_path) != prepared['files']['withheld-point.json']:
        raise ArithmeticError('withheld point changed')
    point = tuple(map(F, read(withheld_path)['point']))
    mapper = load('evaluation_factor_free_mapper', CAS/'factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    rows = []
    for i, centre in enumerate(data['centres']):
        mapping = mapper.mapping(model, basis, centre)
        search = PointedQuarticSearch(state=state, centre={'coefficients': centre['representative']},
                                     coordinate_policy=mapping['coordinate_policy'])
        backend.validate_map(search, mapping)
        for sign in (1, -1):
            located = point_visibility(search.chart_record(), (point[0], sign*point[1]))
            height = located.get('minimum_affine_height')
            rows.append({'chart_index': i, 'sign': sign, 'height': height,
                         'coordinate': located.get('coordinate'),
                         'visible': height is not None and height <= protocol['policy']['height']})
    guard(protocol)
    checkpoint(folder/'evaluation.json', {'status': 'PASS_REFERENCE_CVP_AND_MASKED_VISIBILITY_AUDIT',
        'selection_sha256': seal['selection_sha256'], 'withheld_point_sha256': sha(withheld_path),
        'point_searches': 0, 'rows': rows,
        'visible_signed_representatives': sum(r['visible'] for r in rows),
        'claim_boundary': 'Pointwise retrospective visibility only; no point constructor or new rank gain.'})


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('calibrate302', 'masked90', 'evaluate'))
    p.add_argument('--folder', type=Path, required=True)
    args = p.parse_args()
    evaluate(args.folder) if args.mode == 'evaluate' else run(args.mode, args.folder)
