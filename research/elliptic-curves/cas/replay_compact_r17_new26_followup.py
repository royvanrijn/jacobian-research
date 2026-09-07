#!/usr/bin/env python3
"""Replay the retained new-rank26 R17 follow-up with exact point and admission checks."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import digest
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.cached_observation_state import CachedObservationMWState as MWState


def replay(path):
    data=cert.read(path);protocol=cert.read(path.parents[1]/'protocol.json')
    if data['protocol_hash']!=digest(protocol):raise ArithmeticError('protocol mismatch')
    initial=data['initial_state']['state'];model=tuple(map(cert.F,data['curve']))
    points=tuple(tuple(map(cert.F,p)) for p in initial['reductions']['points']);dimension=len(points)
    if dimension!=data['initial_dimension'] or dimension!=protocol['initial_rank_lower_bound']:raise ArithmeticError('initial dimension changed')
    cert.checked_rank(model,points,initial['reductions']['primes'],initial['no_two_torsion_prime'])
    family=next(r for r in cert.read(spec.ATLAS)['families'] if r['family']==data['family'])
    original,generic=spec.specialize(family,data['parameter']);u=cert.F(data['family_to_curve_scale_u'])
    if not u or model!=(cert.F(0),cert.F(0),cert.F(0),original[3]/u**4,original[4]/u**6):raise ArithmeticError('family model transport differs')
    transported=tuple((x/u**2,y/u**3) for x,y in generic)
    if transported!=points[:17] or transported!=tuple(tuple(map(cert.F,p)) for p in data['generic_points']):raise ArithmeticError('generic basis prefix differs')
    root=Path(__file__).resolve().parents[2]
    for name,h in protocol['sources'].items():
        if cert.hashed(root/name)!=h:raise ArithmeticError('frozen follow-up source differs')
    census_path=root/protocol['generic_census_path']
    if cert.hashed(census_path)!=protocol['generic_census_sha256']:raise ArithmeticError('generic census changed')
    census=cert.read(census_path)
    pool=[{'mask':r['mask'],'generic_norm':r['norm']} for r in sorted(census['records'][1:],key=lambda r:(-cert.F(r['norm']),r['mask']))[:301]]
    if pool!=protocol['generic_pool']:raise ArithmeticError('generic class pool differs')
    words=sorted(range(1,1<<(dimension-17)),key=lambda w:(w.bit_count(),w))
    expected={(g['mask'],words[i%len(words)],g['mask']|(words[i%len(words)]<<17)) for i,g in enumerate(protocol['generic_pool'])}
    actual={(c['generic_mask'],c['quotient_word'],c['parity']) for c in data['centres']}
    if len(data['centres'])!=301 or len(actual)!=301 or actual!=expected:raise ArithmeticError('adaptive class pool changed')
    G=[[round(cert.F(q)*1000000) for q in row] for row in data['metric_gram']]
    if len(G)!=dimension or any(len(r)!=dimension for r in G):raise ArithmeticError('metric dimension differs')
    if sorted(data['centres'],key=lambda c:(-c['metric_norm'],c['generic_mask'],c['quotient_word']))!=data['centres']:raise ArithmeticError('centre order differs')
    for c in data['centres']:
        rep=c['representative']
        if any(type(x) is not int for x in rep):raise ArithmeticError('centre coefficients are not integers')
        if sum(rep[i]*G[i][j]*rep[j] for i in range(dimension) for j in range(dimension))!=c['metric_norm']:raise ArithmeticError('rounded-metric norm differs')
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
    print('REPLAYED RETAINED NEW RANK26 R17 FOLLOWUP',data['family'],data['parameter'],'charts',len(data['charts']),'rank >=',state.rank,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('result',type=Path)
    replay(p.parse_args().result.resolve())
