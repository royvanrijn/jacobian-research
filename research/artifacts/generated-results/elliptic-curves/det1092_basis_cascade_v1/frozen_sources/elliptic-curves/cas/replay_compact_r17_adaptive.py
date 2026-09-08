#!/usr/bin/env python3
"""Replay retained adaptive charts and rank admissions without another sieve."""
import argparse
from pathlib import Path
import sys
from fractions import Fraction as F

ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT/'elliptic-curves/cas'),str(ROOT/'elliptic-curves')]
import certify_compact_r17_candidates as certificate
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.memory_store import MemoryFactStore
from research_runtime.finite_reduction import ReductionCache
from research_runtime.mw_state import MWState
from research_runtime.store import digest


def replay(path):
    result=certificate.read(path)
    protocol=certificate.read(path.parent.parent/'protocol.json')
    if digest(protocol)!=result['protocol_hash']:raise ArithmeticError('protocol binding changed')
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(result['arithmetic_facts'])
    state=MWState.from_record(result['initial_state'],cache=cache)
    model=tuple(map(F,result['curve']));points=[tuple(map(F,p)) for p in state.basis]
    certificate.checked_rank(model,points);certificate.family_check(result['parameter'],model,points)
    primes=tuple(p for p in range(3,protocol['admission_prime_bound']+1) if _is_prime(p))
    completed=[0]*len(result['rounds'])
    for row in result['charts']:
        plan=result['rounds'][row['round']];i=row['chart_index']
        if i!=completed[row['round']] or row['centre']!=plan['centres'][i]:raise ArithmeticError('chart order changed')
        if i==0 and (state.key!=plan['state_key'] or state.rank!=plan['rank_before']):raise ArithmeticError('adaptive round basis changed')
        rep=row['centre']['representative']+[0]*(state.rank-plan['rank_before'])
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.verify_record(row['search'])
        for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        if state.rank!=row['rank_lower_bound']:raise ArithmeticError('rank admission changed')
        completed[row['round']]+=1
    if completed!=[r['completed'] for r in result['rounds']]:raise ArithmeticError('completion count changed')
    if state.record()!=result['final_state'] or state.rank!=result['rank_lower_bound']:raise ArithmeticError('final MWState changed')
    print('REPLAYED ADAPTIVE',result['parameter'],'charts',len(result['charts']),'rank >=',state.rank,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('result',type=Path);a=p.parse_args();replay(a.result)
