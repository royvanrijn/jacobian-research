#!/usr/bin/env python3
"""Replay compact-atlas chart witnesses and finite-rank admissions without resieving."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import digest
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.mw_state import MWState


def replay(path):
    data=cert.read(path);protocol=cert.read(path.parents[2]/'protocol.json')
    if data['protocol_hash']!=digest(protocol):raise ArithmeticError('protocol mismatch')
    f=next(r for r in cert.read(spec.ATLAS)['families'] if r['family']==data['family'])
    model=tuple(map(cert.F,data['curve']));points=tuple(tuple(map(cert.F,p)) for p in data['generic_points'])
    cert.checked_rank(model,points);spec.family_check(f,data['parameter'],model,points)
    if sorted(r['mask'] for r in data['centres'])!=protocol['generic_masks'][data['family']]:raise ArithmeticError('generic centre set changed')
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(data['arithmetic_facts'])
    state=MWState.from_record(data['initial_state'],cache=cache)
    if state.basis!=tuple(tuple(map(str,p)) for p in points):raise ArithmeticError('initial basis mismatch')
    for i,row in enumerate(data['charts']):
        if row['centre']!=data['centres'][i] or row['admission_prime_bound']!=protocol['admission_prime_bound']:raise ArithmeticError('chart plan changed')
        rep=row['centre']['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        result=search.verify_record(row['search'])
        primes=tuple(p for p in range(3,row['admission_prime_bound']+1) if _is_prime(p))
        for point in result.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        if state.rank!=row['rank_lower_bound']:raise ArithmeticError('rank admission changed')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound']:raise ArithmeticError('final state changed')
    print('REPLAYED ATLAS SEARCH',data['family'],data['parameter'],'charts',len(data['charts']),'rank >=',state.rank,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('result',type=Path)
    replay(p.parse_args().result.resolve())
