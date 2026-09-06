#!/usr/bin/env python3
"""Current GMP point search on49 fixed generic-only PARI coordinate maps."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.pointed_orbit_compression import compress
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';LOCAL=ROOT/'artifacts/local/elliptic-curves'
D=LOCAL/'native11952-pari49-control-v1';PARENT=LOCAL/'native11952-metric49-control-v1';MAPS=LOCAL/'native11952-pari-maps-v1/maps.json'

def sources():
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),CAS/'research_runtime/pointed_orbit_compression.py',CAS/'research_runtime/quotient_only_reduction.py',CAS/'research_runtime/cached_observation_state.py')}}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve coordinate comparison')
    maps=cert.read(MAPS);source=cert.read(PARENT/'candidate-00/result.json')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or len(maps['rows'])!=49 or [r['centre'] for r in maps['rows']]!=source['centres']:raise ArithmeticError('fixed generic centres changed')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native11952-pari49.v1','sources':sources(),'parent_sha256':cert.hashed(PARENT/'candidate-00/result.json'),'maps_sha256':cert.hashed(MAPS),'generic_input_sha256':cert.hashed(PARENT/'redacted-input.json'),'height':100000,'seconds_per_chart':8,'worker_wall_seconds':500,'worker_rss_bytes':1610612736,'maximum_workers':1,'admission_prime_bound':997,'gate':'Completed metric49 control certified25 while the older blind generic43 points certify a larger subgroup. Test coordinate reduction with identical generic17 centres, all49 fixed before searching, and the current GMP sieve. PARI performs model transforms only. No exceptional points enter map selection or search.','claim_boundary':'Known-control calibration, not prospective discovery. Finite coverage and exact point independence only. All49 charts attempted regardless of gain.'})

def run():
    protocol=cert.read(D/'protocol.json')
    if sources()!=protocol['sources'] or cert.hashed(MAPS)!=protocol['maps_sha256'] or cert.hashed(PARENT/'candidate-00/result.json')!=protocol['parent_sha256']:raise ArithmeticError('fixed comparison changed')
    source=cert.read(PARENT/'candidate-00/result.json');maps=cert.read(MAPS);cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(source['arithmetic_facts']);state=MWState.from_record(source['initial_state'],cache=cache)
    if state.rank!=17:raise ArithmeticError('must start with only generic17')
    output=D/'candidate-00/result.json'
    if output.exists():raise FileExistsError('preserve comparison attempt')
    data={k:source[k] for k in ('family','parameter','curve','generic_points','initial_state','initial_dimension','metric_gram','centres')};data.update(schema='elliptic-curves.native11952-pari49-result.v1',protocol_hash=digest(protocol),charts=[],status='RUNNING',rank_lower_bound=17,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(output,data)
    primes=tuple(p for p in range(3,998) if _is_prime(p));model=tuple(map(cert.F,data['curve']))
    for row in maps['rows']:
        c=row['centre'];rep=c['representative']+[0]*(state.rank-17);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=row['coordinate_policy'])
        outcome=search.search(100000,8,checkpoint_dir=output.parent/'charts'/state.key);compression=compress(model,state.basis,rep,outcome.curve_points)
        for index in compression['kept_indices']:
            state=state.adjoin(outcome.curve_points[index],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        data['charts'].append({'centre':c,'coordinate_policy':row['coordinate_policy'],'search':outcome.record,'admission_compression':compression,'rank_lower_bound':state.rank,'admission_prime_bound':997});data.update(rank_lower_bound=state.rank,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(output,data);print('NATIVE11952 PARI MAP GMP',len(data['charts']),'rank',state.rank,'through',outcome.record['completed_denominator'],flush=True)
    data['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,data)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run']);a=p.parse_args();prepare() if a.stage=='prepare' else run()
