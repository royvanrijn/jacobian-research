#!/usr/bin/env python3
"""Sage-free exact model-map and point-history replay for the PARI-map control."""
import argparse
from pathlib import Path
from math import isqrt
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch
from half_lattice_pointed_sieve import linear_combination
from search_observability import transform,multiply,prepare_chart
from research_runtime.store import digest
from research_runtime.pointed_orbit_compression import compress
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from search_native11952_pari49 import ROOT,PARENT,MAPS,sources

def replay(path):
    data=cert.read(path);protocol=cert.read(path.parents[1]/'protocol.json');parent=cert.read(PARENT/'candidate-00/result.json');maps=cert.read(MAPS)
    if data['protocol_hash']!=digest(protocol) or protocol['sources']!=sources() or cert.hashed(MAPS)!=protocol['maps_sha256'] or cert.hashed(PARENT/'candidate-00/result.json')!=protocol['parent_sha256']:raise ArithmeticError('control binding changed')
    if cert.hashed(PARENT/'redacted-input.json')!=protocol['generic_input_sha256']:raise ArithmeticError('generic-only input differs')
    source=cert.read(PARENT/'redacted-input.json');model=tuple(map(cert.F,data['curve']));points=tuple(tuple(map(cert.F,p)) for p in data['generic_points'])
    original=tuple(map(cert.F,source['boundary_short_model']));u=cert.F(source['boundary_to_curve_scale_u'])
    if model!=(0,0,0,original[3]/u**4,original[4]/u**6) or points!=tuple((cert.F(x)/u**2,cert.F(y)/u**3) for x,y in source['boundary_generic_points']):raise ArithmeticError('generic transport differs')
    if data['initial_state']!=parent['initial_state'] or data['centres']!=parent['centres']:raise ArithmeticError('initial generic state or fixed centres changed')
    initial=data['initial_state']['state'];cert.checked_rank(model,points,initial['reductions']['primes'],initial['no_two_torsion_prime'])
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(data['arithmetic_facts']);state=MWState.from_record(data['initial_state'],cache=cache)
    if state.rank!=17 or tuple(tuple(map(cert.F,p)) for p in state.basis)!=points:raise ArithmeticError('generic basis changed')
    primes=tuple(p for p in range(3,998) if _is_prime(p))
    for i,row in enumerate(data['charts']):
        saved=maps['rows'][i];c=data['centres'][i]
        if row['centre']!=c or saved['centre']!=c or row['coordinate_policy']!=saved['coordinate_policy']:raise ArithmeticError('fixed chart plan differs')
        x,y=linear_combination(model,points,c['representative']);raw=(-3*x*x-4*model[3],-8*y,-6*x,cert.F(0),cert.F(1))
        if list(map(str,raw))!=saved['raw_coefficients']:raise ArithmeticError('raw quartic differs')
        M=tuple(map(cert.F,saved['matrix']));first=tuple(map(cert.F,saved['first_matrix']));second=tuple(map(cert.F,saved['second_matrix']))
        if M!=multiply(first,second) or saved['coordinate_policy']!={'kind':'raw','matrix':saved['matrix']}:raise ArithmeticError('horizontal composition differs')
        P,Q=[list(map(cert.F,saved[k])) for k in ('reduced_P','reduced_Q')]
        disc=tuple(4*P[j]+sum(Q[k]*Q[j-k] for k in range(3) if 0<=j-k<3) for j in range(5));ratio=cert.F(saved['square_ratio'])
        if list(map(str,disc))!=saved['discriminant_quartic'] or ratio<=0 or isqrt(ratio.numerator)**2!=ratio.numerator or isqrt(ratio.denominator)**2!=ratio.denominator or transform(raw,M)!=tuple(ratio*f for f in disc):raise ArithmeticError('exact reduced quartic identity differs')
        rep=c['representative']+[0]*(state.rank-17);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=saved['coordinate_policy']);prepare_chart(row['search']);outcome=search.verify_record(row['search']);compression=compress(model,state.basis,rep,outcome.curve_points)
        if compression!=row['admission_compression']:raise ArithmeticError('orbit witness differs')
        for j in compression['kept_indices']:
            state=state.adjoin(outcome.curve_points[j],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        if state.rank!=row['rank_lower_bound']:raise ArithmeticError('chart rank differs')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound']:raise ArithmeticError('final state differs')
    if data['status']=='COMPLETE_DECLARED_PILOT' and len(data['charts'])!=49:raise ArithmeticError('incomplete49 plan')
    print('REPLAYED NATIVE11952 PARI49 MAPS AND POINTS',len(data['charts']),'rank >=',state.rank,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('result',type=Path);replay(p.parse_args().result.resolve())
