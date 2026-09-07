#!/usr/bin/env python3
"""Continue fixed charts with a preloaded prime bank and archived observations."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.pointed_orbit_compression import compress
from research_runtime.preloaded_prime_state import preload
from research_runtime.rotated_observation_state import rotate
import pari_pointed_backend as backend
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';LOCAL=ROOT/'artifacts/local/elliptic-curves'
CASES={'rank26':('new-rank26-fast301-v1','new-rank26-adaptive-pari-maps-v1',26,17,32),'small-conductor':('small-conductor-fast301-v1','small-conductor-adaptive-pari-maps-v1',22,16,23)}

def paths(case):
    parent,maps,rank,generic,target=CASES[case]
    return LOCAL/(case+'-primebank-continuation-v1'),LOCAL/parent,LOCAL/maps/'maps.json',rank,generic,target

def sources():
    return {**backend.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),CAS/'research_runtime/preloaded_prime_state.py',CAS/'research_runtime/rotated_observation_state.py',CAS/'research_runtime/pointed_orbit_compression.py',CAS/'research_runtime/cached_observation_state.py',CAS/'research_runtime/quotient_only_reduction.py')}}

def prepare(case):
    d,parent,maps,rank,generic,target=paths(case);seed_path=parent/'candidate-00/result.json'
    if (d/'protocol.json').exists():raise FileExistsError('preserve continuation protocol')
    terminal,seed=cert.read(parent/'terminal.json'),cert.read(seed_path)
    if terminal['result_sha256']!=cert.hashed(seed_path) or terminal['supervision']['outcome']=='running':raise ArithmeticError('seed not terminal')
    start=len(seed['charts']);mapping=cert.read(maps)
    if not 0<start<301 or seed['rank_lower_bound']!=rank or mapping['status']!='COMPLETE_DECLARED_MAPS' or [r['centre'] for r in mapping['rows']]!=seed['centres']:raise ArithmeticError('fixed continuation roster differs')
    bench=LOCAL/'preloaded-prime-state-benchmark-v1/result.json'
    if cert.read(bench)['status']!='PASS':raise ArithmeticError('prime-bank calibration failed')
    checkpoint(d/'protocol.json',{'schema':'elliptic-curves.fixed-pari-primebank-continuation.v1','case':case,'sources':sources(),'seed_path':str(seed_path.relative_to(ROOT)),'seed_sha256':cert.hashed(seed_path),'seed_terminal_sha256':cert.hashed(parent/'terminal.json'),'maps_path':str(maps.relative_to(ROOT)),'maps_sha256':cert.hashed(maps),'primebank_benchmark_sha256':cert.hashed(bench),'initial_rank':rank,'generic_dimension':generic,'start_chart':start,'end_chart':301,'height':100000,'seconds_per_chart':5,'admission_prime_bound':997,'gp_binary_sha256':cert.hashed(Path('/usr/bin/gp')),'target_rank':target,'worker_wall_seconds':600,'worker_rss_bytes':1610612736,'maximum_workers':1,'gate':'An exact81-point comparison preserves every basis/rank/status/relation admission while preloading finite primes once. Observation metadata do not enter the reviewed admission decision. Archive the complete prior state before each chart, preserve the source prefix, and search only the remaining fixed maps. Tests include a newly independent point, dependence, signs, archive restoration and infinity.','claim_boundary':'Continuation with changed state keys and explicit archived observation history. No completed source chart is rerun. Point-search timeouts have no denominator-prefix claim. All states, raw outputs, points and admissions require replay; no upper rank bound.'})

def initial(protocol,cache):
    seed=cert.read(ROOT/protocol['seed_path']);cache.store.import_snapshot(seed['arithmetic_facts']);state=MWState.from_record(seed['final_state'],cache=cache);points=tuple(tuple(map(cert.F,p)) for p in state.basis)
    cert.checked_rank(tuple(map(cert.F,seed['curve'])),points,state.reductions.primes,state.no_two_torsion_prime)
    if state.rank!=protocol['initial_rank']:raise ArithmeticError('seed lower bound differs')
    state,bank=preload(state,cache,protocol['admission_prime_bound'])
    return seed,state,bank

def validate(protocol):
    if protocol['sources']!=sources() or cert.hashed(ROOT/protocol['seed_path'])!=protocol['seed_sha256'] or cert.hashed(ROOT/protocol['maps_path'])!=protocol['maps_sha256']:raise ArithmeticError('frozen continuation binding differs')

def run(case):
    d,*_=paths(case);protocol=cert.read(d/'protocol.json');validate(protocol);cache=ReductionCache(MemoryFactStore());seed,state,bank=initial(protocol,cache);maps=cert.read(ROOT/protocol['maps_path']);output=d/'candidate-00/result.json'
    if output.exists():raise FileExistsError('preserve continuation result')
    data={k:seed[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u','centres')};data.update(schema='elliptic-curves.fixed-pari-primebank-result.v1',protocol_hash=digest(protocol),initial_dimension=protocol['initial_rank'],initial_state=state.record(),prime_bank=bank,charts=[],status='RUNNING',rank_lower_bound=state.rank,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(output,data);model=tuple(map(cert.F,data['curve']))
    for index in range(protocol['start_chart'],protocol['end_chart']):
        mapping=maps['rows'][index];centre=mapping['centre'];state,archive=rotate(state);archive_path=d/'states'/f'{index:03}.json';checkpoint(archive_path,archive)
        rep=centre['representative']+[0]*(state.rank-protocol['initial_rank']);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=mapping['coordinate_policy']);record,points=backend.execute(search,mapping,protocol['height'],protocol['seconds_per_chart'],protocol['gp_binary_sha256']);compression=compress(model,state.basis,rep,points)
        for i in compression['kept_indices']:state=state.adjoin(points[i],cache=cache)
        final=state.record();data['charts'].append({'parent_chart':index,'centre':centre,'archive_path':str(archive_path.relative_to(ROOT)),'archive_sha256':cert.hashed(archive_path),'search':record,'admission_compression':compression,'admission_observations':final['state']['observations'],'rank_lower_bound':state.rank,'state_key':state.key});data.update(rank_lower_bound=state.rank,final_state=final,arithmetic_facts=cache.store.snapshot());checkpoint(output,data);print('PRIMEBANK CONTINUATION',case,index+1,record['status'],'rank',state.rank,flush=True)
        if state.rank>=protocol['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(output,data);return
    data['status']='COMPLETE_DECLARED_CONTINUATION';checkpoint(output,data)

def replay(case):
    d,*_=paths(case);protocol=cert.read(d/'protocol.json');validate(protocol);data=cert.read(d/'candidate-00/result.json');cache=ReductionCache(MemoryFactStore());seed,state,bank=initial(protocol,cache);maps=cert.read(ROOT/protocol['maps_path'])
    if any(data[k]!=seed[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u')):raise ArithmeticError('seed curve metadata differs')
    if data['protocol_hash']!=digest(protocol) or data['initial_state']!=state.record() or data['prime_bank']!=bank or data['centres']!=seed['centres']:raise ArithmeticError('initial bank/state differs')
    model=tuple(map(cert.F,data['curve']));generic=tuple(tuple(map(cert.F,p)) for p in data['generic_points'])
    if protocol['generic_dimension']==16:
        from audit_prospective_mw16_ambiguities import family_check
        family_check(data,model,generic)
    else:
        import compact_atlas_specialization as spec
        family=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']==data['family']);original,points=spec.specialize(family,data['parameter']);u=cert.F(data['family_to_curve_scale_u'])
        if not u or model!=(0,0,0,original[3]/u**4,original[4]/u**6) or generic!=tuple((x/u**2,y/u**3) for x,y in points):raise ArithmeticError('family transport differs')
    if tuple(tuple(map(cert.F,p)) for p in state.basis[:protocol['generic_dimension']])!=generic:raise ArithmeticError('generic prefix differs')
    for offset,row in enumerate(data['charts']):
        index=protocol['start_chart']+offset;mapping=maps['rows'][index]
        if row['parent_chart']!=index or row['centre']!=mapping['centre']:raise ArithmeticError('fixed remaining chart order differs')
        state,archive=rotate(state);path=ROOT/row['archive_path']
        if cert.hashed(path)!=row['archive_sha256'] or cert.read(path)!=archive:raise ArithmeticError('archived state history differs')
        rep=row['centre']['representative']+[0]*(state.rank-protocol['initial_rank']);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=mapping['coordinate_policy']);r=row['search']
        if r['height_bound']!=protocol['height'] or r['timeout_seconds']!=protocol['seconds_per_chart'] or r['gp_binary_sha256']!=protocol['gp_binary_sha256']:raise ArithmeticError('search budget/backend changed')
        points=backend.replay(search,mapping,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('orbit compression differs')
        for i in compression['kept_indices']:state=state.adjoin(points[i],cache=cache)
        if state.rank!=row['rank_lower_bound'] or state.key!=row['state_key'] or state.record()['state']['observations']!=row['admission_observations']:raise ArithmeticError('admission decisions/history differ')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound']:raise ArithmeticError('final state differs')
    if data['status']=='COMPLETE_DECLARED_CONTINUATION' and len(data['charts'])!=protocol['end_chart']-protocol['start_chart']:raise ArithmeticError('continuation roster incomplete')
    print('REPLAYED PRIMEBANK CONTINUATION',case,len(data['charts']),'rank >=',state.rank,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);p.add_argument('--case',choices=CASES,required=True);a=p.parse_args();globals()[a.stage](a.case)
