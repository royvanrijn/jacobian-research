#!/usr/bin/env python3
"""One fixed49-map million-height pilot on the smallest-coefficient new27 curve."""
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
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'new27-million-height-pilot-v1'
ROSTER=(('retention-11952','retention24-r17-pari-v1/11952-113','retention_high_rank_minimal_proofs_v1.json','11952-113'),)
PROOFS=('paired_high_rank_minimal_proofs_v2.json','next24_high_rank_minimal_proofs_v4.json','retention_high_rank_minimal_proofs_v1.json')
def sources():
    names=['new27_million_height_pilot.py','memory_rank_certificate.py','research_runtime/search_state.py','research_runtime/cached_observation_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve new27 height protocol')
    gates=[ART/'native29_million_chart_benchmark_v1.json',ART/'retention_rank27_followup_portable_replay_v1.json'];g=[cert.read(p) for p in gates]
    if any(r['status']!='PASS' for r in g) or not any(r['arm']=='pari' and r['status']=='bounded_search_complete' and r['rank_lower_bound']==29 and r['wall_seconds']<=60 for r in g[0]['arms']):raise ArithmeticError('completed native29 million-height gate required')
    candidates=[r for name in PROOFS for r in cert.read(ART/name)['curves'] if r['rank_lower_bound']==27 and not r['icarm_matches'] and not r['previous_matches']]
    candidates.sort(key=lambda r:(max(abs(int(v)).bit_length() for v in r['minimal_curve']),r['family'],r['parameter']))
    if len(candidates)!=5 or candidates[0]['id']!='11952-113':raise ArithmeticError('smallest-coefficient new27 choice differs')
    rows=[]
    for name,dirname,proofname,identifier in ROSTER:
        path=LOCAL/dirname/'result.json';mp=LOCAL/dirname/'maps.json';seed,maps=cert.read(path),cert.read(mp);proof=next(r for r in cert.read(ART/proofname)['curves'] if r['id']==identifier)
        if seed['rank_lower_bound']!=27 or seed['curve']!=proof['discovery_curve'] or seed['final_state']['state']['reductions']['points']!=proof['discovery_points'] or proof['icarm_matches'] or proof['previous_matches'] or maps['status']!='COMPLETE_DECLARED_MAPS' or maps['centres']!=seed['centres']:raise ArithmeticError('certified new27 seed/map gate failed')
        old=proof['rank_certificate'];checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in proof['discovery_points']],[s['prime'] for s in old['signatures']],old['no_rational_2_torsion_prime'])
        rows.append({'id':name,'seed_path':str(path.relative_to(ROOT)),'seed_sha256':cert.hashed(path),'maps_path':str(mp.relative_to(ROOT)),'maps_sha256':cert.hashed(mp),'proof_path':str((ART/proofname).relative_to(ROOT)),'proof_sha256':cert.hashed(ART/proofname),'charts':len(maps['rows'])})
    if sum(r['charts'] for r in rows)!=49:raise ArithmeticError('fixed49 maps differ')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.new27-million-height-pilot.v1','sources':sources(),'rows':rows,'height':1000000,'seconds_per_chart':60,'worker_wall_seconds':3300,'replay_wall_seconds':600,'target_rank':28,'rss_bytes':1610612736,'maximum_workers':1,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'height_control_bindings':{str(p.relative_to(ROOT)):cert.hashed(p) for p in gates},'selection_proof_hashes':{str((ART/n).relative_to(ROOT)):cert.hashed(ART/n) for n in PROOFS},'selection':{'pool':[{'id':r['id'],'family':r['family'],'parameter':r['parameter'],'largest_minimal_coefficient_bits':max(abs(int(v)).bit_length() for v in r['minimal_curve'])} for r in candidates],'chosen':'11952-113','reason':'Smallest largest minimal coefficient among the five exactly certified new27 curves, used as a cost-first pilot without claiming coefficient size predicts rank.'},'gate':'PARI completed the original native29 control chart12 at one million in38.013 seconds and independently recovered29 from its blind28 subgroup; the same GMP arm remained censored after60seconds. The chart choice was retrospective, so this is a finite cost/visibility gate, not a blind general recovery guarantee. The selected prospective curve has completed generic and301 adaptive125000 boxes without a28th certified direction.','policy':'Exactly the49 original generic maps of new11952 at2012/211, unchanged order and coordinates. Start from its independently certified27-point subgroup. Uniform height1000000 and60seconds per chart, checkpoint every completed or censored chart. Stop immediately when the worker certifies28 pending independent replay. No oracle points, known-control chart choices, new geometry, refill or automatic continuation enter the worker.','boundaries':'At most49 point attempts on one previously new curve. Incomplete outputs remain censored, and no miss supplies an upper rank bound, saturation or point absence.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen source changed')
    for name,h in p['selection_proof_hashes'].items():
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
    data={k:seed[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u')};data.update(schema='elliptic-curves.new27-million-height-result.v1',height=p['height'],protocol_hash=digest(p),maps_sha256=row['maps_sha256'],initial_state=state.record(),initial_dimension=27,centres=maps['centres'],prime_preload=info,charts=[],status='RUNNING',rank_lower_bound=27,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=folder/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-17);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':state.record()['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=state.record(),rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('NEW27 MILLION',row['id'],i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=p['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_HEIGHT_FOLLOWUP';checkpoint(out,data)
def replay(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];data=cert.read(folder/'result.json');maps=cert.read(ROOT/row['maps_path']);cache=ReductionCache(MemoryFactStore());seed,model,state,info=initial(cache,row)
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=row['maps_sha256'] or data['initial_state']!=state.record() or data['prime_preload']!=info or data['centres']!=maps['centres'] or any(data[k]!=seed[k] for k in ('curve','generic_points','family_to_curve_scale_u','family','parameter')):raise ArithmeticError('follow-up initial binding differs')
    for i,r in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/r['archive_path']
        if r['index']!=i or r['centre']!=m['centre'] or cert.hashed(ap)!=r['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('follow-up archive differs')
        rep=m['centre']['representative']+[0]*(state.rank-17);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);s=r['search']
        if s['height_bound']!=p['height'] or s['timeout_seconds']!=p['seconds_per_chart'] or s['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('follow-up budget differs')
        points=backend.replay(search,m,s);compression=compress(model,state.basis,rep,points)
        if compression!=r['admission_compression']:raise ArithmeticError('follow-up compression differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=r['state_key'] or state.rank!=r['rank_lower_bound'] or state.record()['state']['observations']!=r['admission_observations']:raise ArithmeticError('follow-up admission differs')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound'] or (data['status']=='COMPLETE_DECLARED_HEIGHT_FOLLOWUP' and len(data['charts'])!=row['charts']):raise ArithmeticError('follow-up final state differs')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED NEW27 MILLION',row['id'],len(data['charts']),state.rank,flush=True)
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
        for f in as_completed(pending):got[pending[f]]=f.result();result['rows']=[got[i] for i in sorted(got)];checkpoint(ledger,result);print('NEW27 MILLION TERMINAL',f.result()['id'],f.result()['status'],f.result().get('rank_lower_bound'),flush=True)
    result['status']='PASS' if all(r['status']=='PASS' for r in result['rows']) else 'COMPLETE_WITH_FAILURES';checkpoint(ledger,result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch','worker','replay']);p.add_argument('--index',type=int);a=p.parse_args();globals()[a.stage](a.index) if a.stage in ('worker','replay') else globals()[a.stage]()
