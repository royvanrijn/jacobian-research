#!/usr/bin/env python3
"""Profile one known-point admission on a hash-pinned retained MWState.

This is a performance diagnostic, not a new point search or rank claim.
"""
import argparse
import cProfile
from pathlib import Path
import pstats
import time
import certify_compact_r17_candidates as cert
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.mw_state import MWState
from research_runtime.store import checkpoint

def profile(path, expected_hash, output):
    if cert.hashed(path) != expected_hash or output.exists():
        raise ArithmeticError('immutable profile input/output gate differs')
    data = cert.read(path)
    cache = ReductionCache(MemoryFactStore())
    cache.store.import_snapshot(data['arithmetic_facts'])
    state = MWState.from_record(data['final_state'], cache=cache)
    profiler = cProfile.Profile()
    start = time.monotonic()
    result = profiler.runcall(state.adjoin, state.basis[0], cache=cache)
    elapsed = time.monotonic()-start
    if result.basis != state.basis or len(result.observations) != len(state.observations)+1:
        raise ArithmeticError('known-point admission semantics differ')
    raw = output.with_suffix('.pstats')
    profiler.dump_stats(raw)
    rows = []
    for (file, line, function), (primitive, total, own, cumulative, callers) in pstats.Stats(profiler).stats.items():
        rows.append({'file': file, 'line': line, 'function': function, 'primitive_calls': primitive,
                     'total_calls': total, 'self_seconds': own, 'cumulative_seconds': cumulative})
    rows.sort(key=lambda r: -r['cumulative_seconds'])
    checkpoint(output, {'schema': 'elliptic-curves.retained-admission-profile.v1',
        'input_sha256': expected_hash, 'worker_sha256': cert.hashed(Path(__file__).resolve()),
        'observations_before': len(state.observations), 'rank_lower_bound': state.rank,
        'known_point': list(state.basis[0]), 'wall_seconds_with_profiler': elapsed,
        'functions': rows, 'profile_sha256': cert.hashed(raw),
        'claim_boundary': 'One exact known-point admission, with unchanged basis. Profiled timings are diagnostic and include profiler overhead.'})
    print('PROFILED ONE KNOWN POINT', elapsed, 'seconds;', len(state.observations), 'retained observations', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--input-sha256', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    profile(a.input, a.input_sha256, a.output)
