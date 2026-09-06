#!/usr/bin/env python3
"""A bounded point experiment on all21 nonsingular compact endpoints."""
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
BATCH=LOCAL/'endpoint-specialized-parity-v1';INDEX=ART/'compact_atlas_endpoints_v2.json';GATE=ART/'compact_endpoint_summary_v1.json'

def configure(index):
    global D,SEED,ROW
    p=cert.read(BATCH/'protocol.json')
    if not 0<=index<len(p['rows']):raise ValueError('fixed endpoint roster only')
    ROW=p['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json'

def dimension():
    return ROW['initial_rank']

def sources():
    names=['audit_endpoint_specialized_trial.py','replay_endpoint_specialized_geometry.sage','replay_retention24_geometry.py','audit_recorded_point_mod2_rank_v3.py','audit_retained_cloud_modl.py','memory_rank_certificate.py','endpoint_specialized_parity_trial.py','prepare_endpoint_specialized_parity.sage','prepare_fresh_r17_pari_batch.sage','fresh_r17_pari_batch.py','prospective_half_lattice_v2.sage','research_runtime/search_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/cached_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def prepare():
    if (BATCH/'protocol.json').exists():raise FileExistsError('preserve endpoint protocol')
    rows=[];seeds={};paths=[INDEX,GATE,ART/'skew_endpoint_portable_replay_v1.json']
    replay=cert.read(paths[-1])
    if replay['status']!='PASS' or replay['logical_stages']!=54:raise ArithmeticError('all54 isolated endpoint/skew checks required')
    index=cert.read(INDEX)
    for r in index['rows']:
        if r['status']!='CERTIFIED_SPECIALIZED_SUBGROUP':continue
        if r['catalogue_matches'] or r['previous_matches']:raise ArithmeticError('exact unmatched endpoint required')
        proof=r['rank_certificate'];rank=r['rank_lower_bound']
        if len(r['points'])!=rank or not 11<=rank<=17:raise ArithmeticError('certified independent endpoint subset required')
        checked_rank(tuple(map(cert.F,r['curve'])),[tuple(map(cert.F,P)) for P in r['points']],[a['prime'] for a in proof['signatures']],proof['no_rational_2_torsion_prime'])
        ident=r['family']+'-'+r['endpoint'];seed={k:r[k] for k in ('family','curve','generic_points','points','rank_certificate')};seed['parameter']=r['endpoint']
        path=BATCH/ident/'seed.json';checkpoint(path,seed);seeds[str(path.relative_to(ROOT))]=cert.hashed(path)
        rows.append({'id':ident,'family':r['family'],'endpoint':r['endpoint'],'initial_rank':rank})
    if len(rows)!=21:raise ArithmeticError('all21 nonsingular endpoints required')
    checkpoint(BATCH/'protocol.json',{'schema':'elliptic-curves.endpoint-specialized-parity.v1','sources':sources(),'inputs':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'seed_hashes':seeds,'rows':rows,'sample_size':256,'sample_domain':'compact-endpoint-specialized-parity-v1','charts':12,'height':125000,'seconds_per_chart':10,'target_rank':22,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'geometry_wall_seconds':180,'worker_wall_seconds':300,'replay_wall_seconds':180,'rss_bytes':2147483648,'maximum_workers':1,'gate':'The exact all22 endpoint audit exposes21 nonsingular curves omitted from both compact parameter scans. Each has its own11to17-point independently certified subgroup, with full generic-section cloud bounds agreeing modulo2,3,5. Test point visibility on all21 before any larger or score-selected endpoint campaign. The existing specialized-parity control demonstrates finite recovery without increasing height but does not guarantee sensitivity here.','centre_policy':'Use each endpoint own certified independent subset only. SHA256 of [fixed domain,index], modulo2^r, supplies the first256 distinct nonzero r-bit masks. Compute384-bit numerical canonical heights, round at10^6, and transport a unimodular LLL basis change and numerical CVP representatives exactly. Select12 largest computed norms with mask ties. Freeze every endpoint map file before any point search. This finite sample has no covering or CVP optimality claim.','selection':'All21 nonsingular zero/infinity fibres from the six R17 and five MW16 compact atlases, in the exact audit order. No record target, public point, score fit, catalogue rank or outcome-dependent refill enters geometry or execution. Only the independently certified subset seeds each state; do not presume independence of all specialized generic sections.','boundaries':'At most252 point boxes, height125000 and ten seconds each, one worker. Stop each curve at22 certified points pending independent replay; preserve every partial or censored record. Geometry and histories have explicit supervisor limits. No automatic refill, high-height escalation, exact-rank, rank-upper-bound, saturation, point-absence or universal-novelty claim.'})

def protocol():
    p=cert.read(BATCH/'protocol.json')
    if p['sources']!=sources() or any(cert.hashed(ROOT/n)!=h for n,h in {**p['inputs'],**p['seed_hashes']}.items()):raise ArithmeticError('frozen prospective cohort inputs differ')
    return p

def masks(p):
    result=[];i=0
    while len(result)<p['sample_size']:
        m=int(digest([p['sample_domain'],i]),16)%(1<<dimension());i+=1
        if m and m not in result:result.append(m)
    return result

def initial(cache):
    seed=cert.read(SEED);model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points']);raw=raw_state(model,points,cache=cache,prime_bound=1000);state=MWState.from_record(raw.record(),cache=cache);state,info=preload(state,cache,997)
    if state.rank!=dimension():raise ArithmeticError('certified endpoint initial subgroup differs')
    return seed,state

def worker():
    p=protocol();maps=cert.read(D/'maps.json');out=D/'result.json';cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if out.exists():raise FileExistsError('preserve adaptive point attempt')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or maps['protocol_hash']!=digest(p) or len(maps['rows'])!=p['charts']:raise ArithmeticError('fixed endpoint maps incomplete')
    data={k:seed[k] for k in ('family','parameter','curve','generic_points')};data.update(protocol_hash=digest(p),maps_sha256=cert.hashed(D/'maps.json'),initial_state=state.record(),initial_dimension=dimension(),centres=maps['centres'],metric_gram=maps['metric_gram'],charts=[],status='RUNNING',rank_lower_bound=dimension(),final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=D/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-dimension());search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record();data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':final['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=final,rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('ENDPOINT SPECIALIZED',i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=p['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT';checkpoint(out,data)

def replay():
    p=protocol();data=cert.read(D/'result.json');maps=cert.read(D/'maps.json');cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(D/'maps.json') or data['initial_state']!=state.record() or any(data[k]!=seed[k] for k in ('family','parameter','curve','generic_points')) or data['centres']!=maps['centres']:raise ArithmeticError('adaptive initial binding differs')
    candidates=maps['sample'];expected=sorted(candidates,key=lambda c:(-c['metric_norm'],c['parity']))[:p['charts']]
    if [c['parity'] for c in candidates]!=masks(p) or maps['centres']!=expected or [m['centre'] for m in maps['rows']]!=expected:raise ArithmeticError('fixed sample and selected roster differ')
    for c in candidates:
        if len(c['representative'])!=dimension() or any((c['representative'][j]-(c['parity']>>j))%2 for j in range(dimension())):raise ArithmeticError('specialized parity differs')
    for i,row in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/row['archive_path']
        if row['index']!=i or row['centre']!=m['centre'] or cert.hashed(ap)!=row['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('adaptive archive differs')
        rep=m['centre']['representative']+[0]*(state.rank-dimension());search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r=row['search']
        if r['height_bound']!=p['height'] or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('adaptive budget differs')
        points=backend.replay(search,m,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('adaptive orbit differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=row['state_key'] or state.rank!=row['rank_lower_bound'] or state.record()['state']['observations']!=row['admission_observations']:raise ArithmeticError('adaptive admissions differ')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound'] or (data['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(data['charts'])!=p['charts']):raise ArithmeticError('adaptive final state differs')
    if data['status'] not in ('COMPLETE_DECLARED_ADAPTIVE_ATTEMPT','TARGET_REACHED_PENDING_INDEPENDENT_REPLAY') or (data['status']=='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY' and state.rank<p['target_rank']):raise ArithmeticError('terminal status differs')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED ENDPOINT SPECIALIZED',len(data['charts']),'rank >=',state.rank,flush=True)

def launch():
    p=protocol()
    # Complete every geometry choice before observing any prospective point result.
    for index in range(len(p['rows'])):
        configure(index);name='maps';path=D/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve prospective geometry supervision')
        s=supervise(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'prepare_endpoint_specialized_parity.sage'),'--index',str(index)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT);print(ROW['id'],name,s['outcome'],s['returncode'],flush=True)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('complete frozen maps required')
    for index in range(len(p['rows'])):
        configure(index)
        for name,seconds in [('worker',p['worker_wall_seconds']),('replay',p['replay_wall_seconds'])]:
            path=D/(name+'.supervisor.json')
            if path.exists():raise FileExistsError('preserve prospective search supervision')
            s=supervise([sys.executable,str(Path(__file__).resolve()),name,'--index',str(index)],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT);print(ROW['id'],name,s['outcome'],s['returncode'],flush=True)
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('prospective worker/replay failed or censored')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','worker','replay','launch']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.stage in ('worker','replay'):configure(a.index)
    globals()[a.stage]()
