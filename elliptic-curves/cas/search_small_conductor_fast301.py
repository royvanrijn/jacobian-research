#!/usr/bin/env python3
"""Fixed301 adaptive-centre follow-up with exact PARI maps and a calibrated backend."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.pointed_orbit_compression import compress
import pari_pointed_backend as backend
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';LOCAL=ROOT/'artifacts/local/elliptic-curves';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=LOCAL/'small-conductor-fast301-v1';SOURCE=LOCAL/'prospective-mw16-small-conductor-followup-v2/candidate-00/result.json';MAPS=LOCAL/'small-conductor-adaptive-pari-maps-v1/maps.json'

def sources():
    return {**backend.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),CAS/'research_runtime/pointed_orbit_compression.py',CAS/'research_runtime/quotient_only_reduction.py',CAS/'research_runtime/cached_observation_state.py')}}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fast301 protocol')
    source,maps=cert.read(SOURCE),cert.read(MAPS);gate=cert.read(ART/'new_rank26_engine_comparison_v2.json')
    if source['rank_lower_bound']!=22 or len(source['centres'])!=301 or maps['status']!='COMPLETE_DECLARED_MAPS' or [r['centre'] for r in maps['rows']]!=source['centres']:raise ArithmeticError('fixed301 geometry differs')
    if len(gate['rows'])!=43 or any(r['status']!='PASS' for r in gate['rows']):raise ArithmeticError('backend equality gate failed')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.small-conductor-fast301.v1','sources':sources(),'source_sha256':cert.hashed(SOURCE),'maps_sha256':cert.hashed(MAPS),'backend_gate_sha256':cert.hashed(ART/'new_rank26_engine_comparison_v2.json'),'gp_binary_sha256':cert.hashed(Path('/usr/bin/gp')),'gp_version':'2.15.4','height':100000,'charts':301,'seconds_per_chart':5,'worker_wall_seconds':900,'worker_rss_bytes':1610612736,'maximum_workers':1,'admission_prime_bound':997,'gate':'The new curve at3/17 has certified22 independent points and an exact76-digit conductor; one additional independent point would improve its listed rank23 conductor comparison. Keep its certified22-point subgroup and all301 previously fixed adaptive centres. The earlier attempt reached127 charts. Separately calibrated PARI coordinates and backend improve control recovery and preserve92 complete affine hit sets. No known-record points or new parameter selection enter the search.','claim_boundary':'No rank upper bound. Each timeout has no invented denominator prefix; raw output and exact infinity witnesses remain. Stop only at rank>=23 pending independent replay, otherwise attempt all301 within the outer cap.'})

def run():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or cert.hashed(SOURCE)!=p['source_sha256'] or cert.hashed(MAPS)!=p['maps_sha256']:raise ArithmeticError('frozen worker binding changed')
    source,maps=cert.read(SOURCE),cert.read(MAPS);cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(source['arithmetic_facts']);state=MWState.from_record(source['initial_state'],cache=cache)
    if state.rank!=22:raise ArithmeticError('initial subgroup differs')
    output=D/'candidate-00/result.json'
    if output.exists():raise FileExistsError('preserve fast301 attempt')
    data={k:source[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u','initial_state','initial_dimension','centres','metric_gram')};data.update(schema='elliptic-curves.small-conductor-fast301-result.v1',protocol_hash=digest(p),charts=[],status='RUNNING',rank_lower_bound=22,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(output,data)
    model=tuple(map(cert.F,data['curve']));primes=tuple(q for q in range(3,998) if _is_prime(q))
    for mapping in maps['rows']:
        c=mapping['centre'];rep=c['representative']+[0]*(state.rank-22);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=mapping['coordinate_policy']);record,points=backend.execute(search,mapping,100000,5,p['gp_binary_sha256']);compression=compress(model,state.basis,rep,points)
        for i in compression['kept_indices']:
            state=state.adjoin(points[i],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        data['charts'].append({'centre':c,'search':record,'admission_compression':compression,'rank_lower_bound':state.rank,'admission_prime_bound':997});data.update(rank_lower_bound=state.rank,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(output,data);print('SMALL CONDUCTOR FAST301',len(data['charts']),record['status'],'rank',state.rank,flush=True)
        if state.rank>=23:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(output,data);return
    data['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,data)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run']);a=p.parse_args();prepare() if a.stage=='prepare' else run()
