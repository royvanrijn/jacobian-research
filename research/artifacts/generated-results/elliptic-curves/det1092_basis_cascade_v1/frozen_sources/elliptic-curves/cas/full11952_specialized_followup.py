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
BATCH=LOCAL/'full11952-specialized-followup-v1';INDEX=ART/'new_high_rank_curve_index_v14.json';GATE=ART/'native28_specialized_parity_adaptive_coverage_v1.json'

def configure(index):
    global D,SEED,ROW
    p=cert.read(BATCH/'protocol.json')
    if not 0<=index<len(p['rows']):raise ValueError('fixed eligible roster only')
    ROW=p['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json'

def sources():
    names=['audit_full11952_specialized_followup.py','replay_full11952_specialized_geometry.sage','export_new_high_rank_curve_index_v14.py','memory_rank_certificate.py','full11952_specialized_followup.py','prepare_full11952_specialized_followup.sage','prepare_fresh_r17_pari_batch.sage','fresh_r17_pari_batch.py','prospective_half_lattice_v2.sage','research_runtime/search_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/cached_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def prepare():
    if (BATCH/'protocol.json').exists():raise FileExistsError('preserve one full11952 follow-up protocol')
    import export_new_high_rank_curve_index_v14 as inventory
    inventory.promotion_gate()
    replay=ART/'new_high_rank_curve_index_v14_memory_replay_v1.json'
    verified=cert.read(replay)
    if verified['status']!='PASS' or any(cert.hashed(ROOT/n)!=h for n,h in verified['sources'].items()):raise ArithmeticError('bound exact V14 replay required')
    gate=cert.read(GATE);rows=[];seeds={};oldpath=ART/'new_high_rank_curve_index_v13.json';paths=[INDEX,GATE,replay,oldpath]
    if gate['mod2_lower_bound']!=29 or gate['odd_modulus_lower_bounds']!={'3':29,'5':29} or gate['adaptive_completed_boxes']!=18:raise ArithmeticError('exact known29 recovery gate required')
    gd=LOCAL/'native28-specialized-parity-control-v1'
    for label in ('replay','cloud-audit','geometry'):
        path=gd/(label+'.supervisor.json');r=cert.read(path);paths.append(path)
        if r['outcome']!='completed' or r['returncode']!=0:raise ArithmeticError('complete control replays required')
    previous_ids={r['id'] for r in cert.read(oldpath)['curves']}
    index=cert.read(INDEX)
    eligible=sorted((r for r in index['curves'] if r['id'] not in previous_ids and r['rank_lower_bound']==27 and r['source_certificate']=='full11952_64_r17_results_v1.json'),key=lambda r:r['id'])
    if len(eligible)!=1 or eligible[0]['id']!='new-20260906-186' or eligible[0]['parameter']!='4286/1881':raise ArithmeticError('exactly the new full11952 rank27 curve4286/1881 required')
    for r in eligible:
        source=ART/r['source_certificate'];paths.append(source);q=cert.read(source)['curves'][r['source_curve_index']];g=q['generic_points'];rank=r['rank_lower_bound']
        if q['icarm_matches'] or q['previous_matches'] or r['points'][:len(g)]!=g or len(g)!=17 or rank!=len(r['points']):raise ArithmeticError('new catalogue-unmatched own27 seed and generic17 prefix required')
        proof=r['rank_certificate'];actual=checked_rank(tuple(map(cert.F,r['curve'])),[tuple(map(cert.F,P)) for P in r['points']],[a['prime'] for a in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if cert.read(source)['curves'][r['source_curve_index']]['rank_certificate']!=proof or digest(actual)!=digest(proof):raise ArithmeticError('exact seed proof differs')
        seed={k:q[k] for k in ('family','parameter','curve','generic_points','points','rank_certificate')};path=BATCH/r['id']/'seed.json'
        if path.exists():raise FileExistsError('preserve prospective seed')
        checkpoint(path,seed);seeds[str(path.relative_to(ROOT))]=cert.hashed(path)
        rows.append({**{k:r[k] for k in ('id','family','parameter','source_certificate','source_curve_index')},'initial_rank':rank})
    checkpoint(BATCH/'protocol.json',{'schema':'elliptic-curves.full11952-specialized-followup.v1','sources':sources(),'inputs':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'seed_hashes':seeds,'rows':rows,'sample_size':2048,'sample_domain':'full11952-specialized-followup-v1','charts':49,'height':125000,'seconds_per_chart':10,'target_rank':28,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'geometry_wall_seconds':300,'worker_wall_seconds':900,'replay_wall_seconds':600,'rss_bytes':2147483648,'maximum_workers':1,'maximum_point_boxes':49*len(rows),'gate':'The known native control recovered its29th direction from28 known points using the same2048-mask specialized-parity policy at125000 on chart18. The fixed full11952 initial searches use only17 generic-point centres; every new27 seed supplies additional exact independent directions not used to define those initial centres. Test this changed point exposure on the entire eligible new full11952 roster. The control proves finite recovery, not prospective sensitivity, rank density or solubility of other covers.','centre_policy':'Use only each curve own certified27-point seed. Generate2048 distinct SHA256 masks of that dimension with nonzero quotient above the17-point generic prefix, using the fixed sample domain and increasing index. Compute384-bit canonical heights, round at10^6, use a unimodular LLL change and numerical CVP, exactly check parities and rounded norms, then select49 largest computed norms with mask ties. Freeze every cohort map file before any point search. No covering or optimality theorem is claimed.','selection':'All V14 entries newly added beyond V13 from full11952_64_r17_results_v1.json with certified lower bound27, fixed to4286/1881 and stable local ID186. Require the whole eligible roster to contain exactly one curve; otherwise stop before writing any seed. No public points, oracle words, new parameter scans, score changes, selection among previous adaptive outcomes or result-dependent refill.','boundaries':'Exactly one finite follow-up of49 charts per eligible curve, at most49 boxes,125000 height and ten seconds each, one point worker. Stop each curve at its first28-point lower bound pending independent replay; retain all timeouts and partial records. No second wave, larger heights, automatic retry, exact-rank, upper-bound, saturation, record or universal-novelty claim.'})

def protocol():
    p=cert.read(BATCH/'protocol.json')
    if p['sources']!=sources() or any(cert.hashed(ROOT/n)!=h for n,h in {**p['inputs'],**p['seed_hashes']}.items()):raise ArithmeticError('frozen prospective cohort inputs differ')
    return p

def masks(p):
    result=[];i=0;generic=len(cert.read(SEED)['generic_points'])
    while len(result)<p['sample_size']:
        m=int(digest([p['sample_domain'],i]),16)%(1<<ROW['initial_rank']);i+=1
        if m>>generic and m not in result:result.append(m)
    return result

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
        final=state.record();data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':final['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=final,rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('ADAPTIVE FULL11952 SPECIALIZED',i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=p['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(out,data);return
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
    if data['status'] not in ('COMPLETE_DECLARED_ADAPTIVE_ATTEMPT','TARGET_REACHED_PENDING_INDEPENDENT_REPLAY') or (data['status']=='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY' and state.rank<p['target_rank']):raise ArithmeticError('unsupported target stop or nonterminal point status')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED ADAPTIVE FULL11952 SPECIALIZED',len(data['charts']),'rank >=',state.rank,flush=True)

def launch():
    p=protocol();out=BATCH/'ledger.json'
    if out.exists():raise FileExistsError('preserve fixed follow-up ledger')
    ledger={'status':'RUNNING_GEOMETRY','maps':[],'rows':[]};checkpoint(out,ledger)
    # Every geometry choice precedes every new point result.
    for index in range(len(p['rows'])):
        configure(index);name='maps';path=D/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve prospective geometry supervision')
        s=supervise(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'prepare_full11952_specialized_followup.sage'),'--index',str(index)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT)
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
        result=cert.read(D/'result.json');entry.update(status='PASS',rank_lower_bound=result['rank_lower_bound'],result_sha256=cert.hashed(D/'result.json'));checkpoint(out,ledger)
    ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','worker','replay','launch']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.stage in ('worker','replay'):configure(a.index)
    try:globals()[a.stage]()
    except Exception as error:
        if a.stage=='launch' and (BATCH/'ledger.json').exists():
            d=cert.read(BATCH/'ledger.json');d.update(status='FAILED_OR_CENSORED',reason=str(error));checkpoint(BATCH/'ledger.json',d)
        raise
