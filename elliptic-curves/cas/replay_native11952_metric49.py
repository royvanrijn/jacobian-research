#!/usr/bin/env python3
"""Exact replay of the native11952 generic-point-only49-chart control."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import digest
from research_runtime.pointed_orbit_compression import compress
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.cached_observation_state import CachedObservationMWState as MWState

def replay(path):
    data=cert.read(path);protocol=cert.read(path.parents[1]/'protocol.json')
    if data['protocol_hash']!=digest(protocol):raise ArithmeticError('protocol mismatch')
    source_path=path.parents[1]/'redacted-input.json';source=cert.read(source_path)
    if cert.hashed(source_path)!=protocol['input_sha256']:raise ArithmeticError('redacted control input differs')
    model=tuple(map(cert.F,data['curve']));points=tuple(tuple(map(cert.F,p)) for p in data['generic_points'])
    original=tuple(map(cert.F,source['boundary_short_model']));u=cert.F(source['boundary_to_curve_scale_u'])
    if not u or model!=(cert.F(0),cert.F(0),cert.F(0),original[3]/u**4,original[4]/u**6):raise ArithmeticError('control model transport differs')
    if points!=tuple((cert.F(x)/u**2,cert.F(y)/u**3) for x,y in source['boundary_generic_points']):raise ArithmeticError('generic-only point transport differs')
    initial=data['initial_state']['state'];initial_points=tuple(tuple(map(cert.F,p)) for p in initial['reductions']['points'])
    cert.checked_rank(model,initial_points,initial['reductions']['primes'],initial['no_two_torsion_prime'])
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(data['arithmetic_facts']);state=MWState.from_record(data['initial_state'],cache=cache)
    if state.rank!=17 or initial_points!=points:raise ArithmeticError('initial generic basis differs')
    for phase,key in [('initial43','initial43_masks'),('additional6','additional6_masks')]:
        block=[c for c in data['centres'] if c['phase']==phase]
        if sorted(c['mask'] for c in block)!=sorted(protocol[key]) or block!=sorted(block,key=lambda c:(-c['metric_norm'],c['mask'])):raise ArithmeticError('fixed control phase differs')
    if [c['phase'] for c in data['centres']]!=['initial43']*43+['additional6']*6:raise ArithmeticError('control phase order differs')
    G=[[round(cert.F(x)*1000000) for x in row] for row in data['metric_gram']]
    for c in data['centres']:
        rep=c['representative']
        if len(rep)!=17 or any(type(v) is not int for v in rep) or any((rep[j]-(c['mask']>>j))%2 for j in range(17)):raise ArithmeticError('centre parity differs')
        if sum(rep[i]*G[i][j]*rep[j] for i in range(17) for j in range(17))!=c['metric_norm']:raise ArithmeticError('metric norm differs')
    primes=tuple(p for p in range(3,protocol['admission_prime_bound']+1) if _is_prime(p))
    for i,row in enumerate(data['charts']):
        if row['centre']!=data['centres'][i] or row['admission_prime_bound']!=protocol['admission_prime_bound']:raise ArithmeticError('chart plan changed')
        rep=row['centre']['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.verify_record(row['search'])
        compression=compress(model,state.basis,rep,outcome.curve_points)
        if compression!=row['admission_compression']:raise ArithmeticError('orbit witness differs')
        for j in compression['kept_indices']:
            state=state.adjoin(outcome.curve_points[j],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        if state.rank!=row['rank_lower_bound']:raise ArithmeticError('rank admission changed')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound']:raise ArithmeticError('final state changed')
    if data['status']=='COMPLETE_DECLARED_PILOT' and len(data['charts'])!=49:raise ArithmeticError('incomplete declared plan')
    print('REPLAYED NATIVE11952 METRIC49 CONTROL',data['family'],data['parameter'],'charts',len(data['charts']),'rank >=',state.rank,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('result',type=Path);replay(p.parse_args().result.resolve())
