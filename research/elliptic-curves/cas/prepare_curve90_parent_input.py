#!/usr/bin/env python3
"""Consolidate the sealed native M27 and productive bank for fresh parent derivation."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT/'artifacts/local/elliptic-curves'


def run(output):
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    read = lambda p: json.loads(p.read_text())
    prep = LOCAL/'curve90-v3-preparation-v1'
    source = LOCAL/'curve90-productive-anchor-bank-v1'
    sealed = read(prep/'prepared.json')
    assert all(sha(prep/n) == h for n, h in sealed['files'].items())
    bank = read(source/'anchor-bank.json')
    assert bank['status'] == 'COMPLETE_FROZEN_PRODUCTIVE_SUBSET'
    assert bank['protocol_sha256'] == sha(source/'protocol.json')
    assert bank['gains_sha256'] == sha(source/'gains.json')
    old = read(source/'protocol.json')
    assert all(sha(ROOT/n) == h for n, h in old['inputs'].items())
    seed = read(prep/'seed-M27.json')
    assert seed['rank_lower_bound'] == 27 and seed['generic_rank'] == bank['dimension'] == 16
    proof = seed['proof']
    fresh = checked_rank(tuple(map(Fraction, seed['curve'])),
                         tuple(tuple(map(Fraction, p)) for p in seed['points']),
                         [r['prime'] for r in proof['signatures']],
                         proof['no_rational_2_torsion_prime'])
    assert json.loads(json.dumps(fresh)) == proof
    geometry = read(prep/'generic-metric.json')
    assert geometry['gram_scale'] == bank['generic_gram_scale'] == 2
    g = geometry['gram']
    for row in bank['rows']:
        w = row['word']
        assert sum((x % 2) << j for j, x in enumerate(w)) == row['mask']
        assert sum(w[i]*g[i][j]*w[j] for i in range(16) for j in range(16)) == row['norm']
    native = read(prep/'protocol.json')
    paths = [prep/n for n in ('prepared.json','protocol.json','seed-M27.json','generic-metric.json')]
    paths += [source/n for n in ('protocol.json','gains.json','anchor-bank.json')]
    protocol = dict(family=native['family'], parameter=native['parameter'], generic_rank=16,
                    initial_rank=27, point_searches=0,
                    rule='Verbatim certified M27 and historical productive words; no masked oracle input.',
                    inputs={str(p.relative_to(ROOT)):sha(p) for p in paths},
                    source_sha256=sha(Path(__file__)))
    output.mkdir(exist_ok=False)
    checkpoint(output/'protocol.json', protocol)
    checkpoint(output/'seed-M27.json', seed)
    checkpoint(output/'gains.json', read(source/'gains.json'))
    checkpoint(output/'generic-cvp-proofs.json', geometry)
    bank['protocol_sha256'] = sha(output/'protocol.json')
    bank['gains_sha256'] = sha(output/'gains.json')
    checkpoint(output/'anchor-bank.json', bank)
    checkpoint(output/'prepared.json', dict(status='PASS_NATIVE_M27_PRODUCTIVE_INPUT_CONSOLIDATION',
        point_searches=0, files={p.name:sha(p) for p in sorted(output.glob('*.json'))}))
    print('CURVE90_PARENT_INPUT_PREPARED', len(bank['rows']), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
