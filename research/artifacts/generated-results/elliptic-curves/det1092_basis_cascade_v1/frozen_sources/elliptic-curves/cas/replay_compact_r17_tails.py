#!/usr/bin/env python3
"""Replay exact coordinate matching, contiguous tails and rank admissions."""
import argparse
from pathlib import Path
from fractions import Fraction as F
import sys

ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT/'elliptic-curves/cas'),str(ROOT/'elliptic-curves')]
import certify_compact_r17_candidates as certificate
from research_runtime.store import digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.finite_reduction import ReductionCache
from research_runtime.mw_state import MWState
from pointed_quartic_search import PointedQuarticSearch
from mod2_reduction_independence import _is_prime


def replay(path):
    d=certificate.read(path);protocol=certificate.read(path.parent.parent/'protocol.json')
    if digest(protocol)!=d['protocol_hash']:raise ArithmeticError('tail protocol changed')
    row=next(r for r in protocol['rows'] if r['parameter']==d['parameter'])
    source=ROOT/row['input_path']
    if certificate.hashed(source)!=row['input_sha256']:raise ArithmeticError('retained prefix changed')
    old=certificate.read(source);cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(d['arithmetic_facts'])
    state=MWState.from_record(d['initial_state'],cache=cache)
    certificate.family_check(d['parameter'],tuple(map(F,d['curve'])),[tuple(map(F,p)) for p in state.basis])
    primes=tuple(p for p in range(3,252) if _is_prime(p))
    if len(d['charts'])!=len(old['charts']):raise ArithmeticError('tail chart count mismatch')
    for i,r in enumerate(d['charts']):
        prefix=old['charts'][i]['search']
        if prefix['denominator_start']!=1 or prefix['denominator_end']!=protocol['height'] or prefix['height_bound']!=protocol['height']:
            raise ArithmeticError('retained prefix has a different box')
        if r['chart_index']!=i or r['centre']!=old['charts'][i]['centre'] or r['retained_prefix_end']!=prefix['completed_denominator']:
            raise ArithmeticError('tail addressing mismatch')
        rep=r['centre']['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        for key in ('base_point','coefficients','horizontal_matrix','ordinate_scale','pointed_chart','short_model','short_model_x_shift'):
            if search.chart_record()[key]!=prefix[key]:raise ArithmeticError('tail coordinate map changed')
        if 'search' in r:
            tail=r['search']
            if tail['denominator_start']!=r['retained_prefix_end']+1 or tail['denominator_end']!=prefix['denominator_end'] or tail['height_bound']!=prefix['height_bound']:
                raise ArithmeticError('tail is not contiguous in the original box')
            outcome=search.verify_record(tail)
            for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
            if r['union_completed_denominator']!=tail['completed_denominator']:raise ArithmeticError('tail endpoint changed')
        elif r['union_completed_denominator']!=r['retained_prefix_end'] or r['retained_prefix_end']!=prefix['denominator_end']:
            raise ArithmeticError('unfinished prefix has no tail witness')
        if state.rank!=r['rank_lower_bound']:raise ArithmeticError('tail rank changed')
    if state.record()!=d['final_state'] or state.rank!=d['rank_lower_bound']:raise ArithmeticError('tail final state changed')
    complete=sum(r['union_completed_denominator']==protocol['height'] for r in d['charts'])
    if complete!=d['full_boxes_completed']:raise ArithmeticError('full-box count changed')
    print('REPLAYED TAILS',d['parameter'],'full_boxes',complete,'rank >=',state.rank,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('result',type=Path);a=p.parse_args();replay(a.result)
