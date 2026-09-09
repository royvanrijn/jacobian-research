#!/usr/bin/env python3
"""Compare future integer CVP against frozen V3 on eight retained cosets.

Uses the first anchor at M18/M26/M30 and first/middle/last refined entries.
No point search, selector change, or write to a historical campaign.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from research_runtime.store import checkpoint
from visibility_lattice_fast import IntegerExactParity
from visibility_lattice_v2 import ExactParity

ROOT = Path(__file__).resolve().parents[2]
CAS = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    if output.exists():
        raise FileExistsError('preserve earlier benchmark')
    folder = ROOT / 'artifacts/local/elliptic-curves/curve302-seed-universality-panel-v1/residual-strict-03/replay-M17'
    result = {'status': 'RUNNING', 'point_searches': 0, 'inputs': {}, 'rows': [],
              'sources': {str(p.relative_to(ROOT)): sha(p) for p in (
                  Path(__file__), CAS / 'visibility_lattice_fast.py', CAS / 'visibility_lattice_v2.py')}}
    checkpoint(output, result)
    for epoch in (0, 8, 12):
        wd = folder / f'epoch-{epoch:02d}'
        selection = wd / 'selection.json'
        full = wd / 'anchor-00-full.npz'
        for p in (selection, full):
            result['inputs'][str(p.relative_to(ROOT))] = sha(p)
        data = json.loads(selection.read_text())
        anchor = data['anchors'][0]
        if sha(full) != anchor['full_scores_sha256']:
            raise ArithmeticError('retained anchor array binding differs')
        g, u = (np.array(data[k], dtype=object) for k in ('rounded_gram', 'LLL'))
        reduced = (u @ g @ u.T).tolist()
        old, new = ExactParity(reduced), IntegerExactParity(reduced)
        rows = anchor['refined']
        with np.load(full, allow_pickle=False) as arrays:
            for position in sorted({0, len(rows)//2, len(rows)-1}):
                row = rows[position]
                index = row['extension']
                residue, seed = arrays['reduced_residues'][index], arrays['babai_words'][index]
                t = time.perf_counter()
                expected = old.solve(residue, seed)
                old_seconds = time.perf_counter() - t
                t = time.perf_counter()
                actual = new.solve(residue, seed)
                new_seconds = time.perf_counter() - t
                if actual != expected or json.loads(json.dumps(actual)) != row['cvp']:
                    raise ArithmeticError('norm, minimizers, node count or retained proof differs')
                result['rows'].append({'rank': data['rank'], 'extension': index,
                                       'nodes': actual['nodes'], 'reference_seconds': old_seconds,
                                       'integer_seconds': new_seconds,
                                       'speed_ratio': old_seconds / new_seconds,
                                       'exact_reference_and_retained_match': True})
                checkpoint(output, result)
                print(result['rows'][-1], flush=True)
    for key in ('inputs', 'sources'):
        if any(sha(ROOT / p) != h for p, h in result[key].items()):
            raise ArithmeticError('benchmark binding changed')
    if len(result['rows']) != 8:
        raise ArithmeticError('fixed eight-coset roster differs')
    result['status'] = 'PASS_EIGHT_RETAINED_CVP_EQUIVALENCE_BENCHMARK'
    result['aggregate_speed_ratio'] = (sum(r['reference_seconds'] for r in result['rows']) /
                                       sum(r['integer_seconds'] for r in result['rows']))
    result['claim_boundary'] = 'Microbenchmark on eight fixed cosets; not an end-to-end speedup or new rank result.'
    checkpoint(output, result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output)
