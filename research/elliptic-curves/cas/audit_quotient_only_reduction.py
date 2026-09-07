#!/usr/bin/env python3
"""Compare a bounded quotient-only cache with an independently replayed point cloud."""
from dataclasses import asdict
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT/'artifacts/generated-results/elliptic-curves/recorded_rank25_mod2_recertification_0_v1.json'

def run(output):
    if output.exists():
        raise FileExistsError('preserve quotient-only validation')
    data = cert.read(INPUT)
    cache = QuotientOnlyReductionCache(MemoryFactStore(), point_cache_limit=512)
    points = [tuple(map(cert.F, p)) for p in data['points']]
    count = 0
    for old in data['signatures']:
        actual = asdict(cache.signature(data['curve'], points, old['prime']))
        if json.dumps(actual, sort_keys=True) != json.dumps(old, sort_keys=True):
            raise ArithmeticError('signature differs from the independent batch proof')
        count += len(points)
    facts = cache.store.snapshot()
    if any(f['record']['key']['namespace'] != 'finite-field/mod2-quotient' for f in facts['facts']):
        raise ArithmeticError('per-point fact unexpectedly persisted')
    if len(cache._points) > 512:
        raise ArithmeticError('point cache limit exceeded')
    paths = [Path(__file__).resolve(), INPUT,
             ROOT/'elliptic-curves/cas/research_runtime/quotient_only_reduction.py',
             ROOT/'elliptic-curves/cas/research_runtime/finite_reduction.py',
             ROOT/'elliptic-curves/tests/test_quotient_only_reduction.py']
    checkpoint(output, {'schema':'elliptic-curves.quotient-only-cache-validation.v1',
        'status':'PASS', 'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'rational_points':len(points), 'point_prime_masks_checked':count,
        'retained_quotient_facts':len(facts['facts']), 'retained_point_facts':0,
        'in_memory_point_cache_entries':len(cache._points), 'point_cache_limit':512,
        'claim_boundary':'Arithmetic equivalence on this fixed finite cloud; no performance theorem or new curve.'})
    print('QUOTIENT-ONLY CACHE PASS', count, 'masks;', len(facts['facts']), 'quotient facts', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args().output)
