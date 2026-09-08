#!/usr/bin/env python3
"""Exact replay of fixed rank26 tails and their union with the original intervals."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import digest,checkpoint
from research_runtime.pointed_orbit_compression import compress
from search_new_rank26_tails import GEOMETRY_KEYS
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.cached_observation_state import CachedObservationMWState as MWState

def replay(path,output):
    if output.exists():raise FileExistsError('preserve tail coverage proof')
    data=cert.read(path);protocol=cert.read(path.parents[1]/'protocol.json')
    if data['protocol_hash']!=digest(protocol):raise ArithmeticError('protocol mismatch')
    f=next(r for r in cert.read(spec.ATLAS)['families'] if r['family']==data['family'])
    original,generic=spec.specialize(f,data['parameter']);u=cert.F(data['family_to_curve_scale_u'])
    model=tuple(map(cert.F,data['curve']));points=tuple(tuple(map(cert.F,p)) for p in data['generic_points'])
    if not u or model!=(cert.F(0),cert.F(0),cert.F(0),original[3]/u**4,original[4]/u**6):raise ArithmeticError('model transport failed')
    if tuple((x/u**2,y/u**3) for x,y in generic)!=points:raise ArithmeticError('point transport failed')
    initial=data['initial_state']['state'];initial_points=tuple(tuple(map(cert.F,p)) for p in initial['reductions']['points'])
    cert.checked_rank(model,initial_points,initial['reductions']['primes'],initial['no_two_torsion_prime'])
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(data['arithmetic_facts']);state=MWState.from_record(data['initial_state'],cache=cache)
    root=Path(__file__).resolve().parents[2];source_path=root/protocol['input_path']
    if cert.hashed(source_path)!=protocol['input_sha256']:raise ArithmeticError('parent source differs')
    source=cert.read(source_path)
    if data['initial_state']!=source['final_state'] or data['centres']!=source['centres']:raise ArithmeticError('initial state or centres differ')
    if state.rank!=data['initial_dimension'] or initial_points[:17]!=points:raise ArithmeticError('generic prefix differs')
    coverage=[]
    primes=tuple(p for p in range(3,protocol['admission_prime_bound']+1) if _is_prime(p))
    for i,row in enumerate(data['charts']):
        tail=protocol['tails'][i];j=tail['parent_chart'];parent=source['charts'][j]
        if row['parent_chart']!=j or row['centre']!=parent['centre'] or row['admission_prime_bound']!=protocol['admission_prime_bound']:raise ArithmeticError('tail plan differs')
        r=row['search'];old=parent['search']
        if r['denominator_start']!=tail['denominator_start'] or r['denominator_start']!=old['completed_denominator']+1 or r['denominator_end']!=tail['denominator_end'] or r['denominator_end']!=old['denominator_end'] or r['height_bound']!=old['height_bound'] or old['denominator_start']!=1:raise ArithmeticError('interval union differs')
        rep=row['centre']['representative']+[0]*(state.rank-source['initial_dimension'])
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        current=search.chart_record()
        if any(current[k]!=old[k] for k in GEOMETRY_KEYS):raise ArithmeticError('original chart map differs')
        outcome=search.verify_record(row['search'])
        compression=compress(model,state.basis,rep,outcome.curve_points)
        if compression!=row['admission_compression']:raise ArithmeticError('orbit witness differs')
        for kept_index in compression['kept_indices']:
            state=state.adjoin(outcome.curve_points[kept_index],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        if state.rank!=row['rank_lower_bound']:raise ArithmeticError('rank admission changed')
        count=old['integer_pairs_covered']+r['integer_pairs_covered']
        if count!=r['completed_denominator']*(2*r['height_bound']+1):raise ArithmeticError('integer pair coverage differs')
        coverage.append({'parent_chart':j,'completed_denominator':r['completed_denominator'],'height':r['height_bound'],'integer_pairs_covered':count,'combined_box_complete':r['completed_denominator']==r['denominator_end']})
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound']:raise ArithmeticError('final state changed')
    if data['status']=='COMPLETE_DECLARED_PILOT' and len(data['charts'])!=len(protocol['tails']):raise ArithmeticError('incomplete declared plan')
    checkpoint(output,{'schema':'elliptic-curves.new-rank26-combined-box-coverage.v2','status':'PASS',
        'input_path':str(path.relative_to(root)),'input_sha256':cert.hashed(path),'parent_input_sha256':cert.hashed(source_path),
        'sources':{str(p.relative_to(root)):cert.hashed(p) for p in (Path(__file__).resolve(),root/'elliptic-curves/cas/search_new_rank26_tails.py')},
        'rank_lower_bound':state.rank,'charts_replayed':len(coverage),'planned_tails':len(protocol['tails']),
        'complete_combined_boxes':sum(r['combined_box_complete'] for r in coverage),'coverage':coverage,
        'claim_boundary':'Exact retained maps, point admissions and contiguous interval union. Completeness relies on the pinned sieve executions, not a second sieve enumeration. Even an exhausted finite box supplies no upper rank bound.'})
    print('REPLAYED FIXED NEW RANK26 TAILS',data['family'],data['parameter'],'charts',len(data['charts']),'rank >=',state.rank,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();replay(a.input.resolve(),a.output.resolve())
