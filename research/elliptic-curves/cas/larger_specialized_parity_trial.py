#!/usr/bin/env python3
"""One fixed specialized-parity follow-up on the new full11952 rank27 curves."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
import pari_pointed_backend as backend
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run as supervise,Limits
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.rotated_observation_state import rotate
from research_runtime.pointed_orbit_compression import compress
from research_runtime.search_state import raw_state
from research_runtime.preloaded_prime_state import preload
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
BATCH=LOCAL/'larger-specialized-parity-trial-v1';INDEX=ART/'new_high_rank_curve_index_v14.json';GATE=ART/'native28_specialized_parity_adaptive_coverage_v1.json'

def configure(index):
    global D,SEED,ROW
    p=cert.read(BATCH/'protocol.json')
    if not 0<=index<len(p['rows']):raise ValueError('fixed eligible roster only')
    ROW=p['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json'

def sources():
    names=['larger_specialized_parity_trial.py','prepare_larger_specialized_parity_trial.sage','audit_larger_specialized_parity_trial.py','memory_rank_certificate.py','prepare_fresh_r17_pari_batch.sage','fresh_r17_pari_batch.py','prospective_half_lattice_v2.sage']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def prepare():
    if (BATCH/'protocol.json').exists():raise FileExistsError('preserve larger-parity point trial')
    import audit_larger_specialized_parity as geometry
    audit=geometry.expected()
    if cert.read(geometry.OUT)!=audit or audit['status']!='PASS_FIXED_GEOMETRY_AUDIT':raise ArithmeticError('complete65536 geometry audit required')
    rows=[];seeds={};inputs={}
    paths=[geometry.OUT,geometry.D/'protocol.json',geometry.D/'ledger.json',geometry.D/'check.supervisor.json']
    if cert.read(geometry.D/'ledger.json')['status']!='PASS':raise ArithmeticError('complete geometry replay required')
    for case in geometry.CASES:
        identifier=case['id'];folder=BATCH/identifier;seed_path=folder/'seed.json'
        original=LOCAL/'native11952-height125-control-v1/125000/result.json' if identifier=='native28' else LOCAL/'full11952-specialized-followup-v1/new-20260906-186/seed.json'
        data=cert.read(original);points=data['final_state']['state']['reductions']['points'] if identifier=='native28' else data['points']
        if len(points)!=case['rank']:raise ArithmeticError('original recovered seed only')
        payload={k:data[k] for k in ('family','parameter','curve','generic_points')};payload['points']=points
        if seed_path.exists():raise FileExistsError('preserve seed')
        checkpoint(seed_path,payload);seeds[str(seed_path.relative_to(ROOT))]=cert.hashed(seed_path)
        sample=geometry.D/identifier/'sample.json';paths.extend([original,ROOT/case['maps'],sample])
        rows.append({'id':identifier,'initial_rank':case['rank'],'target_rank':case['rank']+1,'sample_path':str(sample.relative_to(ROOT)),'old_maps_path':case['maps']})
    checkpoint(BATCH/'protocol.json',{'schema':'elliptic-curves.larger-specialized-parity-trial.v1','sources':sources(),
        'inputs':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'seed_hashes':seeds,'rows':rows,
        'sample_size':65536,'charts':49,'height':125000,'seconds_per_chart':10,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),
        'geometry_wall_seconds':120,'worker_wall_seconds':900,'replay_wall_seconds':600,'audit_wall_seconds':900,'rss_bytes':2147483648,
        'maximum_workers':1,'maximum_point_boxes':98,
        'gate':'The fixed65536-mask geometry audit gives46 new top49 classes on the known28 control and47 on own27 ID186, with all exact parities and rounded norms replayed. Test whether this changed exposure still recovers the known29th direction before applying it to the fixed own curve. Both map rosters precede all points. The unrecovered public point and oracle words are never input.',
        'boundaries':'Run the known28 control first, stopping at provisional29; require exact history, geometry and cloud proofs modulo2,3,5 before the own27 search. If it does not recover29, the own curve remains unsearched under this protocol. Otherwise run its49 charts and stop at provisional28. Maximum98 boxes, unchanged125000 height and10 seconds per chart, one point worker, no new parameter, score change, further sample expansion, automatic retry or inference of exact rank/absence/covering optimality. The separate live late64 cohort is unchanged.'})

def protocol():
    p=cert.read(BATCH/'protocol.json')
    if p['sources']!=sources() or any(cert.hashed(ROOT/n)!=h for n,h in {**p['inputs'],**p['seed_hashes']}.items()):raise ArithmeticError('frozen larger parity inputs differ')
    return p

def masks(p):return [r['parity'] for r in cert.read(ROOT/ROW['sample_path'])['sample']]

def initial(cache):
    seed=cert.read(SEED);model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points']);raw=raw_state(model,points,cache=cache,prime_bound=1000);state=MWState.from_record(raw.record(),cache=cache);state,info=preload(state,cache,997)
    if state.rank!=ROW['initial_rank']:raise ArithmeticError('certified own27 initial subgroup differs')
    return seed,state

def worker():
    p=protocol();maps=cert.read(D/'maps.json');out=D/'result.json';cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if out.exists():raise FileExistsError('preserve adaptive point attempt')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or maps['protocol_hash']!=digest(p) or len(maps['rows'])!=49:raise ArithmeticError('fixed49 maps incomplete')
    data={k:seed[k] for k in ('family','parameter','curve','generic_points')};data.update(protocol_hash=digest(p),maps_sha256=cert.hashed(D/'maps.json'),initial_state=state.record(),initial_dimension=ROW['initial_rank'],centres=maps['centres'],metric_gram=maps['metric_gram'],charts=[],status='RUNNING',rank_lower_bound=ROW['initial_rank'],final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=D/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-ROW['initial_rank']);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record();data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':final['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=final,rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('ADAPTIVE LARGER PARITY',i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=ROW['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT';checkpoint(out,data)

def replay():
    p=protocol();data=cert.read(D/'result.json');maps=cert.read(D/'maps.json');cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(D/'maps.json') or data['initial_state']!=state.record() or any(data[k]!=seed[k] for k in ('family','parameter','curve','generic_points')) or data['centres']!=maps['centres']:raise ArithmeticError('adaptive initial binding differs')
    candidates=maps['sample'];expected=sorted(candidates,key=lambda c:(-c['metric_norm'],c['parity']))[:49]
    if [c['parity'] for c in candidates]!=masks(p) or maps['centres']!=expected or [m['centre'] for m in maps['rows']]!=expected:raise ArithmeticError('fixed sample and selected roster differ')
    for c in candidates:
        if len(c['representative'])!=ROW['initial_rank'] or any((c['representative'][j]-(c['parity']>>j))%2 for j in range(ROW['initial_rank'])):raise ArithmeticError('specialized parity differs')
    for i,row in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/row['archive_path']
        if row['index']!=i or row['centre']!=m['centre'] or cert.hashed(ap)!=row['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('adaptive archive differs')
        rep=m['centre']['representative']+[0]*(state.rank-ROW['initial_rank']);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r=row['search']
        if r['height_bound']!=p['height'] or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('adaptive budget differs')
        points=backend.replay(search,m,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('adaptive orbit differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=row['state_key'] or state.rank!=row['rank_lower_bound'] or state.record()['state']['observations']!=row['admission_observations']:raise ArithmeticError('adaptive admissions differ')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound'] or (data['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(data['charts'])!=49):raise ArithmeticError('adaptive final state differs')
    if data['status'] not in ('COMPLETE_DECLARED_ADAPTIVE_ATTEMPT','TARGET_REACHED_PENDING_INDEPENDENT_REPLAY') or (data['status']=='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY' and state.rank<ROW['target_rank']):raise ArithmeticError('unsupported target stop or nonterminal point status')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED ADAPTIVE LARGER PARITY',len(data['charts']),'rank >=',state.rank,flush=True)

def launch():
    p=protocol();out=BATCH/'ledger.json'
    if out.exists():raise FileExistsError('preserve fixed follow-up ledger')
    ledger={'status':'RUNNING_GEOMETRY','maps':[],'rows':[]};checkpoint(out,ledger)
    # Every geometry choice precedes every new point result.
    for index in range(len(p['rows'])):
        configure(index);name='maps';path=D/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve prospective geometry supervision')
        s=supervise(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'prepare_larger_specialized_parity_trial.sage'),'--index',str(index)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0
        ledger['maps'].append({'id':ROW['id'],'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(out,ledger)
        print(ROW['id'],name,s['outcome'],s['returncode'],flush=True)
        if not ok:raise ArithmeticError('complete frozen maps required')
    ledger['status']='RUNNING_POINTS';checkpoint(out,ledger)
    for index in range(len(p['rows'])):
        configure(index);entry={'id':ROW['id'],'status':'RUNNING','stages':[]};ledger['rows'].append(entry);checkpoint(out,ledger)
        for name,seconds in [('worker',p['worker_wall_seconds']),('replay',p['replay_wall_seconds'])]:
            path=D/(name+'.supervisor.json')
            if path.exists():raise FileExistsError('preserve prospective search supervision')
            s=supervise([sys.executable,str(Path(__file__).resolve()),name,'--index',str(index)],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT)
            ok=s['outcome']=='completed' and s['returncode']==0
            entry['stages'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(out,ledger)
            print(ROW['id'],name,s['outcome'],s['returncode'],flush=True)
            if not ok:raise ArithmeticError('prospective worker/replay failed or censored')
        path=D/'audit.supervisor.json'
        if path.exists():raise FileExistsError('preserve geometry/cloud audit')
        q=supervise([sys.executable,str(CAS/'audit_larger_specialized_parity_trial.py'),'--index',str(index)],limits=Limits(p['audit_wall_seconds'],p['rss_bytes']),log_path=D/'audit.log',checkpoint_path=path,cwd=ROOT)
        if q['outcome']!='completed' or q['returncode']!=0:raise ArithmeticError('independent geometry/cloud audit failed')
        result=cert.read(D/'result.json');entry.update(status='PASS',rank_lower_bound=result['rank_lower_bound'],result_sha256=cert.hashed(D/'result.json'));checkpoint(out,ledger)
        if index==0:
            coverage=cert.read(ART/'larger_parity_native28_coverage_v1.json')
            if coverage['mod2_lower_bound']<29 or min(coverage['odd_modulus_lower_bounds'].values())<29:
                ledger['status']='CONTROL_NOT_RECOVERED_OWN_CURVE_UNSEARCHED';checkpoint(out,ledger);return
    ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','worker','replay','launch']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.stage in ('worker','replay'):configure(a.index)
    try:globals()[a.stage]()
    except Exception as error:
        if a.stage=='launch' and (BATCH/'ledger.json').exists():
            d=cert.read(BATCH/'ledger.json');d.update(status='FAILED_OR_CENSORED',reason=str(error));checkpoint(BATCH/'ledger.json',d)
        raise
