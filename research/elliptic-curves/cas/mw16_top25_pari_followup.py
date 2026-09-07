#!/usr/bin/env python3
"""Same-height coordinate/backend follow-up on both previously new MW16 rank25 curves."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
import pari_pointed_backend as backend
from memory_rank_certificate import checked_rank
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.preloaded_prime_state import preload
from research_runtime.rotated_observation_state import rotate
from research_runtime.pointed_orbit_compression import compress
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'mw16-top25-pari-followup-v1';INDEX=ART/'new_high_rank_curve_index_v7.json'
ROSTER=(('mw16-04','prospective-mw16-h4096-v1/a1-fibration-04/candidate-01/result.json','prospective_mw16_wide_results_v1.json','-1647/91'),('mw16-05','prospective-mw16-h1024-v1/a1-fibration-05/candidate-03/result.json','prospective_mw16_results_v1.json','307/206'))
def sources():
    names=['mw16_top25_pari_followup.py','prepare_mw16_top25_pari.sage','prepare_fresh_r17_pari_batch.sage','memory_rank_certificate.py','compact_mw16_specialization.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/cached_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py']
    return {**backend.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [*(CAS/n for n in names),spec.ATLAS,INDEX]}}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve MW16 coordinate followup')
    rows=[];inventory=cert.read(INDEX);selected=[r for r in inventory['curves'] if r['family'].startswith('a1-') and r['rank_lower_bound']>=25]
    if {(r['family'],r['parameter'],r['rank_lower_bound']) for r in selected}!={('a1-fibration-04','-1647/91',25),('a1-fibration-05','307/206',25)}:raise ArithmeticError('fixed top MW16 roster differs')
    for identifier,name,proofname,t in ROSTER:
        path=LOCAL/name;seed=cert.read(path);proof=next(r for r in cert.read(ART/proofname)['curves'] if r['parameter']==t and r['family']==seed['family']);points=seed['final_state']['state']['reductions']['points']
        if seed['rank_lower_bound']!=25 or seed['curve']!=proof['curve'] or points!=proof['points'] or proof['icarm_matches'] or proof['previous_matches'] or len(seed['centres'])!=43 or len(seed['generic_points'])!=16:raise ArithmeticError('independently new25 input differs')
        q=proof['rank_certificate'];checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in points],[s['prime'] for s in q['signatures']],q['no_rational_2_torsion_prime'])
        rows.append({'id':identifier,'family':seed['family'],'parameter':t,'seed_path':str(path.relative_to(ROOT)),'seed_sha256':cert.hashed(path),'proof_path':str((ART/proofname).relative_to(ROOT)),'proof_sha256':cert.hashed(ART/proofname),'old_completed_boxes':sum(r['search']['status']=='bounded_search_complete' for r in seed['charts'])})
    gate=ART/'native29_million_chart_benchmark_v1.json'
    if cert.read(gate)['status']!='PASS':raise ArithmeticError('completed backend gate required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.mw16-top25-pari-followup.v1','sources':sources(),'rows':rows,'height':100000,'seconds_per_chart':10,'target_rank':28,'generic_rank':16,'initial_rank':25,'map_wall_seconds':120,'worker_wall_seconds':600,'replay_wall_seconds':300,'rss_bytes':1610612736,'maximum_workers':1,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'gate_sha256':cert.hashed(gate),'gate':'The two highest previously new compactMW16 curves have25 certified points, but their old generic43 GMP boxes at100000 were not completed under the old metric coordinate policy. PARI coordinates/backend now have exact rational-map and point-proof controls on R17, including a completed29th direction. Test the same point centres and same100000 height on both MW16 examples, freezing PARI maps before any search.','policy':'Exactly two pre-existing curves, all43 original generic16 centres on each, unchanged point models and centre order. Start from each own25-point subgroup. Change only the horizontal coordinate reduction/backend and per-chart time cap. One worker alongside the independent million-height pilot. Stop separately at28 pending replay. No public oracle, new candidate sampler, refills or automatic higher-height continuation.','boundaries':'At most86 completed or censored chart attempts. The comparison changes coordinate/backend/time cap together and makes no single-factor causal claim. No exact rank, upper bound, saturation or point absence.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen MW16 sources differ')
    for r in p['rows']:
        for k in ('seed','proof'):
            if cert.hashed(ROOT/r[k+'_path'])!=r[k+'_sha256']:raise ArithmeticError('frozen MW16 input differs')
    return p
def initial(row,cache):
    seed=cert.read(ROOT/row['seed_path']);model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['final_state']['state']['reductions']['points']);state=raw_state(model,points,cache=cache,prime_bound=1000);state=MWState.from_record(state.record(),cache=cache);state,bank=preload(state,cache,997)
    if state.rank!=25 or [list(map(str,P)) for P in points[:16]]!=seed['generic_points']:raise ArithmeticError('initial25 or generic16 prefix differs')
    f=next(f for f in cert.read(spec.ATLAS)['families'] if f['fibration_id']==row['family']);original,generic=spec.specialize(f,row['parameter']);u=cert.F(seed['family_to_curve_scale_u'])
    if model!=(0,0,0,original[3]/u**4,original[4]/u**6) or points[:16]!=tuple((x/u**2,y/u**3) for x,y in generic):raise ArithmeticError('exact generic transport differs')
    return seed,model,state,bank
def worker(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];maps=cert.read(folder/'maps.json');out=folder/'result.json';cache=ReductionCache(MemoryFactStore());seed,model,state,bank=initial(row,cache)
    if out.exists():raise FileExistsError('preserve MW16 attempt')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or maps['protocol_hash']!=digest(p) or maps['centres']!=seed['centres'] or len(maps['rows'])!=43:raise ArithmeticError('fixed43 maps incomplete')
    data={k:seed[k] for k in ('curve','family','parameter','generic_points','family_to_curve_scale_u')};data.update(protocol_hash=digest(p),maps_sha256=cert.hashed(folder/'maps.json'),initial_state=state.record(),prime_preload=bank,centres=maps['centres'],charts=[],status='RUNNING',rank_lower_bound=25,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=folder/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-16);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record();data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':final['state']['observations'],'rank_lower_bound':state.rank,'state_key':state.key});data.update(final_state=final,rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('MW16 TOP25 PARI',row['id'],i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=p['target_rank']:data['status']='TARGET_REACHED_PENDING_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_POINT_ATTEMPT';checkpoint(out,data)
def replay(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];data=cert.read(folder/'result.json');maps=cert.read(folder/'maps.json');cache=ReductionCache(MemoryFactStore());seed,model,state,bank=initial(row,cache)
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(folder/'maps.json') or maps['protocol_hash']!=digest(p) or data['initial_state']!=state.record() or data['prime_preload']!=bank or data['centres']!=seed['centres'] or maps['centres']!=seed['centres']:raise ArithmeticError('MW16 initial binding differs')
    for i,r in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/r['archive_path']
        if r['index']!=i or r['centre']!=m['centre'] or m['centre']!=seed['centres'][i] or cert.hashed(ap)!=r['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('MW16 chart/archive differs')
        rep=m['centre']['representative']+[0]*(state.rank-16);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);s=r['search']
        if s['height_bound']!=p['height'] or s['timeout_seconds']!=p['seconds_per_chart'] or s['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('MW16 budget differs')
        points=backend.replay(search,m,s);compression=compress(model,state.basis,rep,points)
        if compression!=r['admission_compression']:raise ArithmeticError('MW16 compression differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=r['state_key'] or state.rank!=r['rank_lower_bound'] or state.record()['state']['observations']!=r['admission_observations']:raise ArithmeticError('MW16 exact admission differs')
    if data['final_state']!=state.record() or data['rank_lower_bound']!=state.rank or (data['status']=='COMPLETE_DECLARED_POINT_ATTEMPT' and len(data['charts'])!=43):raise ArithmeticError('MW16 terminal state differs')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED MW16 TOP25 PARI',row['id'],len(data['charts']),state.rank,flush=True)
def launch():
    p=protocol();out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve MW16 ledger')
    ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger)
    for index,row in enumerate(p['rows']):
        folder=D/row['id'];stages=[]
        for name,args,seconds in [('maps',['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'prepare_mw16_top25_pari.sage'),'--index',str(index)],120),('worker',[sys.executable,str(Path(__file__).resolve()),'worker','--index',str(index)],600),('replay',[sys.executable,str(Path(__file__).resolve()),'replay','--index',str(index)],300)]:
            r=run(args,limits=Limits(seconds,p['rss_bytes']),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=ROOT);stages.append({'stage':name,'supervision':r});print('MW16 STAGE',row['id'],name,r['outcome'],r['returncode'],flush=True)
            if r['outcome']!='completed' or r['returncode']!=0:break
        ok=len(stages)==3 and all(r['supervision']['outcome']=='completed' and r['supervision']['returncode']==0 for r in stages);ledger['rows'].append({'id':row['id'],'status':'PASS' if ok else 'FAILED_OR_CENSORED','stages':stages});checkpoint(out,ledger)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch','worker','replay']);p.add_argument('--index',type=int);a=p.parse_args();globals()[a.stage](a.index) if a.stage in ('worker','replay') else globals()[a.stage]()
