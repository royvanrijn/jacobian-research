#!/usr/bin/env python3
"""Replay retained compact-atlas adaptive centres, charts and rank admissions."""
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
    data=cert.read(path);protocol=cert.read(path.parents[1]/'protocol.json')
    if data['protocol_hash']!=digest(protocol):raise ArithmeticError('protocol mismatch')
    initial=data['initial_state']['state'];model=tuple(map(cert.F,data['curve']))
    points=tuple(tuple(map(cert.F,p)) for p in initial['reductions']['points']);dimension=len(points)
    if dimension!=data['initial_dimension'] or dimension!=protocol['initial_rank_lower_bound']:raise ArithmeticError('initial dimension changed')
    cert.checked_rank(model,points,initial['reductions']['primes'],initial['no_two_torsion_prime'])
    family=next(r for r in cert.read(spec.ATLAS)['families'] if r['family']==data['family'])
    spec.family_check(family,data['parameter'],model,points)
    words=sorted(range(1,1<<(dimension-17)),key=lambda w:(w.bit_count(),w))
    expected={(g['mask'],words[i%len(words)],g['mask']|(words[i%len(words)]<<17)) for i,g in enumerate(protocol['generic_pool'])}
    actual={(c['generic_mask'],c['quotient_word'],c['parity']) for c in data['centres']}
    if len(data['centres'])!=301 or len(actual)!=301 or actual!=expected:raise ArithmeticError('adaptive class pool changed')
    for c in data['centres']:
        rep=c['representative']
        if len(rep)!=dimension or any((rep[j]-(c['parity']>>j))%2 for j in range(dimension)):raise ArithmeticError('adaptive centre dimension or parity changed')
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(data['arithmetic_facts'])
    state=MWState.from_record(data['initial_state'],cache=cache)
    primes=tuple(p for p in range(3,protocol['admission_prime_bound']+1) if _is_prime(p))
    for i,row in enumerate(data['charts']):
        if row['centre']!=data['centres'][i] or row['admission_prime_bound']!=protocol['admission_prime_bound']:raise ArithmeticError('chart order or allowance changed')
        rep=row['centre']['representative']+[0]*(state.rank-dimension)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        result=search.verify_record(row['search'])
        for point in result.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        if state.rank!=row['rank_lower_bound']:raise ArithmeticError('rank admission changed')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound']:raise ArithmeticError('final state mismatch')
    if data['status']=='COMPLETE_DECLARED_PILOT' and len(data['charts'])!=301:raise ArithmeticError('incomplete declared plan')
    print('REPLAYED RETAINED ATLAS FOLLOWUP',data['family'],data['parameter'],'charts',len(data['charts']),'rank >=',state.rank,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('result',type=Path)
    replay(p.parse_args().result.resolve())
