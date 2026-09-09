#!/usr/bin/env python3
"""Retrospective preconditioned-map visibility at a fixed certified winning centre."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from v3_warm_engine import load, certified_state
from pointed_quartic_search import PointedQuarticSearch
from search_observability import point_visibility
from research_runtime.store import checkpoint
import pari_pointed_backend as backend

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]


def run(epoch, output):
    read = lambda p: json.loads(p.read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    if output.exists():
        raise FileExistsError('preserve audit')
    gain, seed = read(epoch/'gain.json'), read(epoch/'seed.json')
    receipt = epoch/f"chart-{gain['chart_index']:04d}.json"
    chart = read(receipt)
    if sha(receipt) != gain['chart_sha256']:
        raise ArithmeticError('gain receipt changed')
    verified = read(epoch/'reference-verified.json')
    if (verified['status'] != 'PASS_FULL_PRODUCTIVE_EPOCH_REPLAY'
            or verified['stage_sha256'] != sha(epoch/'stage.json')):
        raise ArithmeticError('independent epoch replay required')
    model = tuple(map(F, seed['curve']))
    basis = tuple(tuple(map(F, p)) for p in seed['points'])
    state = certified_state(model, basis, seed['proof'])
    certified_state(model, tuple(tuple(map(F, p)) for p in gain['points']), gain['proof'])
    point = tuple(map(F, gain['points'][-1]))
    rows = []
    names = ('preconditioned_pari_mapping.sage',)
    for name in names:
        module = load('gain_audit_'+name, CAS/name)
        module.pari.allocatemem(256000000, silent=True)
        mapping = module.mapping(model, basis, chart['mapping']['centre'])
        search = PointedQuarticSearch(state=state,
            centre={'coefficients': chart['mapping']['centre']['representative']},
            coordinate_policy=mapping['coordinate_policy'])
        backend.validate_map(search, mapping)
        for sign in (1, -1):
            observed = point_visibility(search.chart_record(), (point[0], sign*point[1]))
            rows.append({'mapper': name, 'sign': sign, 'observation': observed})
    paths = [epoch/n for n in ('gain.json', 'seed.json', 'stage.json', 'reference-verified.json')]
    paths += [receipt, Path(__file__), CAS/'factor_free_pari_mapping.sage', *(CAS/n for n in names)]
    result = {'status': 'PASS_EXACT_WINNING_CENTRE_COORDINATES', 'point_searches': 0,
        'bindings': {str(p.relative_to(ROOT)): sha(p) for p in paths}, 'rows': rows,
        'claim_boundary': 'Post-discovery exact pointwise visibility only. No prospective '
            'selector input, absence result, or generic prediction.'}
    checkpoint(output, result)
    print(json.dumps(rows), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epoch', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.epoch.resolve(), a.output.resolve())
