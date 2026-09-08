#!/usr/bin/env python3
"""Sage-free replay of prospective MW16 transports, chart witnesses and rank admissions."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import digest
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.mw_state import MWState

def replay(path):
    data=cert.read(path);protocol=cert.read(path.parents[2]/'point-protocol.json')
    if data['protocol_hash']!=digest(protocol):raise ArithmeticError('protocol mismatch')
    f=next(r for r in cert.read(spec.ATLAS)['families'] if r['fibration_id']==data['family'])
    original,generic=spec.specialize(f,data['parameter']);u=cert.F(data['family_to_curve_scale_u'])
    model=tuple(map(cert.F,data['curve']));points=tuple(tuple(map(cert.F,p)) for p in data['generic_points'])
    if not u or model!=(cert.F(0),cert.F(0),cert.F(0),original[3]/u**4,original[4]/u**6):raise ArithmeticError('model transport failed')
    if tuple((x/u**2,y/u**3) for x,y in generic)!=points:raise ArithmeticError('point transport failed')
    initial=data['initial_state']['state'];initial_points=tuple(tuple(map(cert.F,p)) for p in initial['reductions']['points'])
    cert.checked_rank(model,initial_points,initial['reductions']['primes'],initial['no_two_torsion_prime'])
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(data['arithmetic_facts']);state=MWState.from_record(data['initial_state'],cache=cache)
    if data['status']=='INCOMPLETE_GENERIC_MOD2_CERTIFICATE':
        if state.rank>=16 or data['charts'] or data['centres']:raise ArithmeticError('censored input semantics changed')
    elif data['centres']:
        if initial_points!=points or state.rank!=16:raise ArithmeticError('initial generic basis mismatch')
        if sorted(c['mask'] for c in data['centres'])!=sorted(protocol['generic_masks'][data['family']]):raise ArithmeticError('generic centre set changed')
        for c in data['centres']:
            rep=c['representative']
            if len(rep)!=16 or any((rep[j]-(c['mask']>>j))%2 for j in range(16)):raise ArithmeticError('centre parity failed')
    primes=tuple(p for p in range(3,protocol['admission_prime_bound']+1) if _is_prime(p))
    for i,row in enumerate(data['charts']):
        if row['centre']!=data['centres'][i] or row['admission_prime_bound']!=protocol['admission_prime_bound']:raise ArithmeticError('chart plan changed')
        rep=row['centre']['representative']+[0]*(state.rank-16)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.verify_record(row['search'])
        for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        if state.rank!=row['rank_lower_bound']:raise ArithmeticError('rank admission changed')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound']:raise ArithmeticError('final state changed')
    if data['status']=='COMPLETE_DECLARED_PILOT' and len(data['charts'])!=43:raise ArithmeticError('incomplete declared plan')
    print('REPLAYED RETAINED PROSPECTIVE MW16',data['family'],data['parameter'],'charts',len(data['charts']),'rank >=',state.rank,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('result',type=Path);replay(p.parse_args().result.resolve())
