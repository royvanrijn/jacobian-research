#!/usr/bin/env python3
"""Sage/PARI-free replay of the fast301 backend maps, squares and admissions."""
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
import pari_pointed_backend as backend
from search_small_conductor_fast301 import ROOT,SOURCE,MAPS,sources

def replay(path):
    data=cert.read(path);protocol=cert.read(path.parents[1]/'protocol.json');parent=cert.read(SOURCE);maps=cert.read(MAPS)
    if data['protocol_hash']!=digest(protocol) or protocol['sources']!=sources() or cert.hashed(SOURCE)!=protocol['source_sha256'] or cert.hashed(MAPS)!=protocol['maps_sha256']:raise ArithmeticError('frozen follow-up binding differs')
    from audit_prospective_mw16_ambiguities import family_check
    model=tuple(map(cert.F,data['curve']));points=tuple(tuple(map(cert.F,p)) for p in data['generic_points'])
    family_check(data,model,points)
    if data['initial_state']!=parent['initial_state'] or data['centres']!=parent['centres'] or data['initial_dimension']!=22:raise ArithmeticError('initial subgroup or centres differ')
    initial=data['initial_state']['state'];initial_points=tuple(tuple(map(cert.F,p)) for p in initial['reductions']['points']);cert.checked_rank(model,initial_points,initial['reductions']['primes'],initial['no_two_torsion_prime'])
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(data['arithmetic_facts']);state=MWState.from_record(data['initial_state'],cache=cache)
    if state.rank!=22 or initial_points[:16]!=points:raise ArithmeticError('generic prefix differs')
    primes=tuple(p for p in range(3,998) if _is_prime(p))
    for i,row in enumerate(data['charts']):
        mapping=maps['rows'][i];c=data['centres'][i]
        if row['centre']!=c or mapping['centre']!=c or row['admission_prime_bound']!=997 or row['search']['height_bound']!=100000 or row['search']['timeout_seconds']!=5 or row['search']['gp_binary_sha256']!=protocol['gp_binary_sha256']:raise ArithmeticError('fixed search plan differs')
        rep=c['representative']+[0]*(state.rank-22);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=mapping['coordinate_policy']);found=backend.replay(search,mapping,row['search']);compression=compress(model,state.basis,rep,found)
        if compression!=row['admission_compression']:raise ArithmeticError('orbit witnesses differ')
        for j in compression['kept_indices']:
            state=state.adjoin(found[j],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        if state.rank!=row['rank_lower_bound']:raise ArithmeticError('chart rank differs')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound']:raise ArithmeticError('final state differs')
    if data['status']=='COMPLETE_DECLARED_PILOT' and len(data['charts'])!=301:raise ArithmeticError('incomplete planned roster')
    print('REPLAYED SMALL CONDUCTOR FAST301',len(data['charts']),'rank >=',state.rank,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('result',type=Path);replay(p.parse_args().result.resolve())
