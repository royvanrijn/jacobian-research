#!/usr/bin/env python3
"""Compare exact old/new state records on a retained long observation history."""
import argparse
from pathlib import Path
import time
import certify_compact_r17_candidates as cert
from research_runtime.cached_observation_state import CachedObservationMWState, _contains
from research_runtime.mw_state import MWState
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.store import digest

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT/'artifacts/local/elliptic-curves/retained-admission-profile-v1/snapshot.json'

def build(output):
    if output.exists():
        raise FileExistsError('preserve cache validation')
    data = cert.read(INPUT)
    old_cache = QuotientOnlyReductionCache(MemoryFactStore())
    new_cache = QuotientOnlyReductionCache(MemoryFactStore())
    for cache in (old_cache, new_cache):
        cache.store.import_snapshot(data['arithmetic_facts'])
    old = MWState.from_record(data['final_state'], cache=old_cache)
    new = CachedObservationMWState.from_record(data['final_state'], cache=new_cache)
    if old.record() != new.record():
        raise ArithmeticError('initial records differ')
    rows = []
    for i in range(3):
        P = old.basis[i]
        t = time.monotonic(); old = old.adjoin(P, cache=old_cache); old_time = time.monotonic()-t
        t = time.monotonic(); new = new.adjoin(P, cache=new_cache); new_time = time.monotonic()-t
        if old.record() != new.record():
            raise ArithmeticError('portable state records or hash chain differ')
        rows.append({'point': list(P), 'old_seconds': old_time, 'cached_seconds': new_time,
                     'identical_record_digest': digest(old.record())})
    sources = (INPUT, Path(__file__).resolve(), ROOT/'elliptic-curves/cas/research_runtime/cached_observation_state.py',
               ROOT/'elliptic-curves/cas/research_runtime/mw_state.py',
               ROOT/'elliptic-curves/cas/research_runtime/quotient_only_reduction.py')
    cert.write(output, {'schema': 'elliptic-curves.cached-observation-state-validation.v1',
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in sources},
        'original_observations': len(data['final_state']['state']['observations']),
        'status': 'PASS_IDENTICAL_EXACT_RECORDS', 'rows': rows, 'cache_info': _contains.cache_info()._asdict(),
        'claim_boundary': 'Exact record equality on three retained long-history known-point admissions. Timings are a small diagnostic, not a general speed guarantee. Separate tests cover ambiguity, new independent points and invalid observations.'})
    print('CACHE VALIDATED', len(rows), 'identical records;', rows, flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    build(p.parse_args().output.resolve())
