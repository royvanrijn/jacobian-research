#!/usr/bin/env python3
"""Exact bounded audit of omitted MW16 half-integral upper shells; zero searches."""
import argparse
import hashlib
import json
from pathlib import Path
import time
from visibility_generic_bank import enumerate_bank
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
CAS = Path(__file__).resolve().parent
PREP = ROOT/'artifacts/local/elliptic-curves/curve90-v3-preparation-v1'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(folder):
    from sage.all import ZZ, matrix
    folder.mkdir(exist_ok=False)
    paths = [PREP/'generic-metric.json', PREP/'prepared.json']
    protocol = {'scaled_bound': 23, 'generic_gram_scale': 2,
        'scaled_upper_shells': [16, 19, 20, 23], 'node_limit': 20000000,
        'point_searches': 0,
        'reason': 'Retrospective parent coverage audit after the norm8/10 masked null. '
                  'Retain all generic minima at least8 through the historical maximum11.5; '
                  'no singled-out winning mask or known point enters enumeration.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'sources': {str(p.relative_to(ROOT)): sha(p) for p in (
            Path(__file__), CAS/'visibility_generic_bank.py', CAS/'visibility_lattice_fast.py',
            CAS/'visibility_lattice_v2.py')}}
    checkpoint(folder/'protocol.json', protocol)
    metric, prepared = [json.loads(p.read_text()) for p in paths]
    if sha(paths[0]) != prepared['files']['generic-metric.json']:
        raise ArithmeticError('prepared generic metric changed')
    g, u = (matrix(ZZ, metric[k]) for k in ('gram', 'LLL'))
    if abs(u.det()) != 1:
        raise ArithmeticError('invalid parent lattice transport')
    t = time.monotonic()
    bank = enumerate_bank([list(map(int, r)) for r in (u*g*u.transpose()).rows()],
        [list(map(int, r)) for r in u.rows()], bound=23, shells=(16, 19, 20, 23),
        node_limit=protocol['node_limit'])
    bank['wall_seconds'] = time.monotonic() - t
    bank['generic_gram_scale'] = 2
    bank['all_nonzero_parities_represented'] = len(bank['represented_parity_minima']) == 2**16-1
    bank['upper_shell_counts'] = {str(q): sum(r['norm'] == q for r in bank['rows'])
                                  for q in (16, 19, 20, 23)}
    for category in ('inputs', 'sources'):
        if any(sha(ROOT/p) != h for p, h in protocol[category].items()):
            raise ArithmeticError('upper-shell audit binding changed')
    checkpoint(folder/'parent-bank.json', bank)
    print('UPPER_SHELLS', bank['upper_shell_counts'], bank['wall_seconds'],
          bank['all_nonzero_parities_represented'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--folder', type=Path, required=True)
    run(p.parse_args().folder)
