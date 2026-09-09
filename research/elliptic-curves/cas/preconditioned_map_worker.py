#!/usr/bin/env python3
"""Construct one exact V3 coordinate map in an externally bounded process.

The supervisor, not this worker, enforces wall/RSS limits. A timeout supplies
no map and no point-search coverage. This worker makes no rank claim.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from v3_warm_engine import load
from research_runtime.store import checkpoint

CAS = Path(__file__).resolve().parent
MAPPERS = {'preconditioned_full': 'preconditioned_full_pari_mapping.sage',
           'factor_free': 'factor_free_pari_mapping.sage'}


def run(input_path, output):
    if output.exists():
        raise FileExistsError('preserve prior map result')
    raw = input_path.read_bytes()
    data = json.loads(raw)
    module_path = CAS/MAPPERS[data['policy']]
    module = load('bounded_one_map', module_path)
    module.pari.allocatemem(256000000, silent=True)
    mapping = module.mapping(tuple(map(F, data['curve'])),
                             tuple(tuple(map(F, p)) for p in data['points']), data['centre'])
    checkpoint(output, dict(status='EXACT_MAP_CONSTRUCTED', mapping=mapping,
                           input_sha256=hashlib.sha256(raw).hexdigest(),
                           mapper_sha256=hashlib.sha256(module_path.read_bytes()).hexdigest(),
                           worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print('EXACT_MAP_CONSTRUCTED', data['policy'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.input.resolve(), args.output.resolve())
