#!/usr/bin/env python3
"""One fixed49-map million-height pilot on new11952 at4286/1881."""
import argparse,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import certify_compact_r17_candidates as cert
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
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'full11952-million-height-v1'
def sources():
    names=['full11952_million_height.py','memory_rank_certificate.py','research_runtime/search_state.py','research_runtime/cached_observation_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve one full11952 million-height protocol')
    gates=[ART/'native29_million_chart_benchmark_v1.json',ART/'full11952_high_rank_models_portable_replay_v1.json',ART/'full11952_specialized_followup_v1.json',ART/'higher26_million_higher_11952_coverage_v1.json']
    g=[cert.read(q) for q in gates]
    if any(r['status']!='PASS' for r in g) or not any(r['arm']=='pari' and r['status']=='bounded_search_complete' and r['rank_lower_bound']==29 and r['wall_seconds']<=60 for r in g[0]['arms']):raise ArithmeticError('exact native29 cost/visibility and completed own27 proof gates required')
    if g[2]['completed_point_boxes']!=49 or g[2]['rows'][0]['rank_lower_bound']!=27:raise ArithmeticError('completed own27 initial specialized boxes required')
    initial=LOCAL/'full11952-64-r17-pari-v1/11952-0962587/result.json';prior=LOCAL/'full11952-specialized-followup-v1/new-20260906-186';path=prior/'result.json';mp=prior/'maps.json';proofpath=ART/'full11952_64_r17_results_v1.json'
    seed,maps=cert.read(path),cert.read(mp);proof=next(r for r in cert.read(proofpath)['curves'] if r['id']=='11952-0962587')
    if seed['parameter']!='4286/1881' or seed['rank_lower_bound']!=27 or seed['curve']!=proof['curve'] or seed['final_state']['state']['reductions']['points']!=proof['points'] or proof['icarm_matches'] or proof['previous_matches'] or maps['status']!='COMPLETE_DECLARED_MAPS' or maps['centres']!=seed['centres'] or len(maps['rows'])!=49 or any(len(m['centre']['representative'])!=27 for m in maps['rows']):raise ArithmeticError('fixed own27 seed and49 specialized maps required')
    old=proof['rank_certificate'];checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in proof['points']],[r['prime'] for r in old['signatures']],old['no_rational_2_torsion_prime'])
    row={'id':'new-20260906-186','seed_path':str(path.relative_to(ROOT)),'seed_sha256':cert.hashed(path),'maps_path':str(mp.relative_to(ROOT)),'maps_sha256':cert.hashed(mp),'proof_path':str(proofpath.relative_to(ROOT)),'proof_sha256':cert.hashed(proofpath),'charts':49}
    prior_inputs=[initial,path,prior/'replay.supervisor.json',LOCAL/'full11952-specialized-followup-controller-v1/ledger.json']
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.full11952-million-height.v1','sources':sources(),'rows':[row],'height':1000000,'seconds_per_chart':60,'worker_wall_seconds':3300,'replay_wall_seconds':600,'target_rank':28,'rss_bytes':1610612736,'maximum_workers':1,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'height_control_bindings':{str(q.relative_to(ROOT)):cert.hashed(q) for q in gates},'selection_proof_hashes':{str(q.relative_to(ROOT)):cert.hashed(q) for q in [proofpath,*prior_inputs]},'prior_point_histories':[str(q.relative_to(ROOT)) for q in [initial,path]],'selection':{'chosen':'new-20260906-186','parameter':'4286/1881','reason':'The new globally minimal27-point curve completes its original49 generic and separate49 specialized boxes at125000. Its specialized maps have not been tested at one million. Reuse exactly those49 frozen own27 maps to test this explicit remaining exposure limit.'},'gate':'The retrospectively chosen native control chart completed at one million within60seconds and recovered29 from its blind28-point subgroup. This is a finite cost/visibility check, not a prospective sensitivity guarantee. Earlier million-height trials on an older27-point curve and the higher-parameter26-point curve were null; neither guarantees success here. This one newly discovered curve has an independently verified minimal model and its27-point seed, all64 initial cohort proofs and the fixed49 specialized follow-up are complete.','policy':'One worker, exactly49 unchanged specialized maps in unchanged order, own27-point seed, uniform height1000000 and60seconds per chart. Retain every completion/timeout, stop at a provisional28 pending independent point proof. No public points, new centre choices, parameter scan, refill or automatic extra wave.','boundaries':'At most49 larger boxes on this one curve. Original49 generic and49 specialized histories remain for the complete witness union. A timeout or miss supplies no rank upper bound, saturation, exact rank or point absence.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen source changed')
    for name,h in {**p['selection_proof_hashes'],**p['height_control_bindings']}.items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('selection pool source changed')
    for r in p['rows']:
        for k in ('seed','maps','proof'):
            if cert.hashed(ROOT/r[k+'_path'])!=r[k+'_sha256']:raise ArithmeticError('frozen follow-up input changed')
    return p
def initial(cache,row):
    seed=cert.read(ROOT/row['seed_path']);model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['final_state']['state']['reductions']['points']);state=raw_state(model,points,cache=cache,prime_bound=1000);state=MWState.from_record(state.record(),cache=cache);state,info=preload(state,cache,997)
    if state.rank!=27 or list(map(list,map(lambda P:list(map(str,P)),points[:17])))!=seed['generic_points']:raise ArithmeticError('certified27 generic prefix differs')
    return seed,model,state,info
def worker(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];out=folder/'result.json';maps=cert.read(ROOT/row['maps_path']);cache=ReductionCache(MemoryFactStore());seed,model,state,info=initial(cache,row)
    if out.exists():raise FileExistsError('preserve follow-up result')
    data={k:seed[k] for k in ('family','parameter','curve','generic_points')};data.update(schema='elliptic-curves.full11952-million-height-result.v1',height=p['height'],protocol_hash=digest(p),maps_sha256=row['maps_sha256'],initial_state=state.record(),initial_dimension=27,centres=maps['centres'],prime_preload=info,charts=[],status='RUNNING',rank_lower_bound=27,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=folder/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-27);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':state.record()['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=state.record(),rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('FULL11952 MILLION',row['id'],i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=p['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_HEIGHT_FOLLOWUP';checkpoint(out,data)
def replay(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];data=cert.read(folder/'result.json');maps=cert.read(ROOT/row['maps_path']);cache=ReductionCache(MemoryFactStore());seed,model,state,info=initial(cache,row)
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=row['maps_sha256'] or data['initial_state']!=state.record() or data['prime_preload']!=info or data['centres']!=maps['centres'] or any(data[k]!=seed[k] for k in ('curve','generic_points','family','parameter')):raise ArithmeticError('follow-up initial binding differs')
    for i,r in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/r['archive_path']
        if r['index']!=i or r['centre']!=m['centre'] or cert.hashed(ap)!=r['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('follow-up archive differs')
        rep=m['centre']['representative']+[0]*(state.rank-27);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);s=r['search']
        if s['height_bound']!=p['height'] or s['timeout_seconds']!=p['seconds_per_chart'] or s['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('follow-up budget differs')
        points=backend.replay(search,m,s);compression=compress(model,state.basis,rep,points)
        if compression!=r['admission_compression']:raise ArithmeticError('follow-up compression differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=r['state_key'] or state.rank!=r['rank_lower_bound'] or state.record()['state']['observations']!=r['admission_observations']:raise ArithmeticError('follow-up admission differs')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound'] or (data['status']=='COMPLETE_DECLARED_HEIGHT_FOLLOWUP' and len(data['charts'])!=row['charts']):raise ArithmeticError('follow-up final state differs')
    if data['status'] not in ('COMPLETE_DECLARED_HEIGHT_FOLLOWUP','TARGET_REACHED_PENDING_INDEPENDENT_REPLAY') or (data['status']=='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY' and state.rank<p['target_rank']):raise ArithmeticError('unsupported terminal million-height status')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED FULL11952 MILLION',row['id'],len(data['charts']),state.rank,flush=True)
def launch():
    p=protocol();ledger=D/'ledger.json'
    if ledger.exists():raise FileExistsError('preserve follow-up ledger')
    result={'status':'RUNNING','rows':[]};checkpoint(ledger,result)
    def one(index):
        row=p['rows'][index];folder=D/row['id'];stages=[]
        for name,seconds in [('worker',p['worker_wall_seconds']),('replay',p['replay_wall_seconds'])]:
            s=run([sys.executable,str(Path(__file__).resolve()),name,'--index',str(index)],limits=Limits(seconds,p['rss_bytes']),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=ROOT);stages.append({'name':name,'supervision':s})
            if s['outcome']!='completed' or s['returncode']!=0:return {'id':row['id'],'status':'FAILED_OR_CENSORED','stages':stages}
        return {'id':row['id'],'status':'PASS','stages':stages,'rank_lower_bound':cert.read(folder/'result.json')['rank_lower_bound'],'result_sha256':cert.hashed(folder/'result.json')}
    with ThreadPoolExecutor(max_workers=1) as pool:
        pending={pool.submit(one,i):i for i in range(len(p['rows']))};got={}
        for f in as_completed(pending):got[pending[f]]=f.result();result['rows']=[got[i] for i in sorted(got)];checkpoint(ledger,result);print('FULL11952 MILLION TERMINAL',f.result()['id'],f.result()['status'],f.result().get('rank_lower_bound'),flush=True)
    result['status']='PASS' if all(r['status']=='PASS' for r in result['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(ledger,result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch','worker','replay']);p.add_argument('--index',type=int);a=p.parse_args();globals()[a.stage](a.index) if a.stage in ('worker','replay') else globals()[a.stage]()
