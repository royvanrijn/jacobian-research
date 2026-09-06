#!/usr/bin/env python3
"""A frozen specialized-parity experiment on all six new27 inventory curves."""
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
BATCH=LOCAL/'new27-specialized-parity-six-v1';INDEX=ART/'new_high_rank_curve_index_v12.json';GATE=ART/'native28_specialized_parity_adaptive_coverage_v1.json'

def configure(index):
    global D,SEED,ROW
    p=cert.read(BATCH/'protocol.json')
    if not 0<=index<6:raise ValueError('fixed six curves only')
    ROW=p['rows'][index];D=BATCH/ROW['id'];SEED=D/'seed.json'

def sources():
    names=['memory_rank_certificate.py','new27_specialized_parity_six.py','prepare_new27_specialized_parity.sage','prepare_fresh_r17_pari_batch.sage','fresh_r17_pari_batch.py','prospective_half_lattice_v2.sage','research_runtime/search_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/cached_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def prepare():
    if (BATCH/'protocol.json').exists():raise FileExistsError('preserve prospective six-curve protocol')
    gate=cert.read(GATE);rows=[];seeds={};paths=[INDEX,GATE]
    if gate['mod2_lower_bound']!=29 or gate['odd_modulus_lower_bounds']!={'3':29,'5':29} or gate['adaptive_completed_boxes']!=18:raise ArithmeticError('exact known29 recovery gate required')
    gd=LOCAL/'native28-specialized-parity-control-v1'
    for label in ('replay','cloud-audit','geometry'):
        path=gd/(label+'.supervisor.json');r=cert.read(path);paths.append(path)
        if r['outcome']!='completed' or r['returncode']!=0:raise ArithmeticError('complete control replays required')
    index=cert.read(INDEX)
    for r in sorted((r for r in index['curves'] if r['rank_lower_bound']==27),key=lambda r:r['id']):
        source=ART/r['source_certificate'];paths.append(source);q=cert.read(source)['curves'][r['source_curve_index']];g=q['generic_points']
        if q['icarm_matches'] or q['previous_matches'] or r['points'][:len(g)]!=g or len(g) not in (16,17):raise ArithmeticError('catalogue-unmatched own27 and generic prefix required')
        proof=r['rank_certificate'];checked_rank(tuple(map(cert.F,r['curve'])),[tuple(map(cert.F,P)) for P in r['points']],[a['prime'] for a in proof['signatures']],proof['no_rational_2_torsion_prime'])
        seed={k:q[k] for k in ('family','parameter','curve','generic_points','points','rank_certificate')};path=BATCH/r['id']/'seed.json';checkpoint(path,seed);seeds[str(path.relative_to(ROOT))]=cert.hashed(path)
        rows.append({k:r[k] for k in ('id','family','parameter','source_certificate','source_curve_index')})
    if len(rows)!=6:raise ArithmeticError('all and only six inventory27 curves required')
    checkpoint(BATCH/'protocol.json',{'schema':'elliptic-curves.new27-specialized-parity-six.v1','sources':sources(),'inputs':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'seed_hashes':seeds,'rows':rows,'initial_rank':27,'sample_size':2048,'sample_domain':'new27-specialized-parity-v1','charts':49,'height':125000,'seconds_per_chart':10,'target_rank':28,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'geometry_wall_seconds':300,'worker_wall_seconds':900,'replay_wall_seconds':300,'rss_bytes':2147483648,'maximum_workers':1,'gate':'The separately frozen specialized-parity control recovers the known29th direction from28 previously recovered points at125000 on its18th chart; exact admissions, full-cloud modulo2,3,5 and2048 parity/norm transports replay. This demonstrates a finite new exposure policy without increasing height. Test it on all six existing catalogue-unmatched inventory curves with27 certified points, without selecting among their previous follow-up outcomes.','centre_policy':'For each own27-point seed, generate2048 distinct27-bit SHA256-derived masks with nonzero quotient above its exact generic prefix of16 or17 points. Digest [fixed sample domain,index] modulo2^27, choose the first2048 eligible distinct masks. Compute384-bit numerical canonical heights, round at10^6, and use a unimodular LLL basis change for numerical CVP. Transport representatives exactly back to the original seed and verify parity and rounded norms. Select49 largest computed norms, mask ties, and freeze all six map files before any point search. This is a finite sample, without a claim of coset or covering optimality.','selection':'All and only the six rank27 entries in the frozen V12 inventory, sorted by stable ID. Only their own certified27 points and generic prefixes enter geometry. No public curve, record point, oracle word, score fit, parameter scan, trace extension or result-dependent refill. The known recovery control calibrates this policy but is not a new curve or an independent population test.','boundaries':'At most294 point boxes,125000 height and ten seconds each, with one worker. Stop each curve at its first28-point lower bound pending independent replay; preserve every partial or censored record. No automatic additional samples, larger heights, exact-rank, upper-bound, saturation or universal-novelty claim.'})

def protocol():
    p=cert.read(BATCH/'protocol.json')
    if p['sources']!=sources() or any(cert.hashed(ROOT/n)!=h for n,h in {**p['inputs'],**p['seed_hashes']}.items()):raise ArithmeticError('frozen prospective cohort inputs differ')
    return p

def masks(p):
    result=[];i=0;generic=len(cert.read(SEED)['generic_points'])
    while len(result)<p['sample_size']:
        m=int(digest([p['sample_domain'],i]),16)%(1<<27);i+=1
        if m>>generic and m not in result:result.append(m)
    return result

def initial(cache):
    seed=cert.read(SEED);model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points']);raw=raw_state(model,points,cache=cache,prime_bound=1000);state=MWState.from_record(raw.record(),cache=cache);state,info=preload(state,cache,997)
    if state.rank!=27:raise ArithmeticError('certified own27 initial subgroup differs')
    return seed,state

def worker():
    p=protocol();maps=cert.read(D/'maps.json');out=D/'result.json';cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if out.exists():raise FileExistsError('preserve adaptive point attempt')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or maps['protocol_hash']!=digest(p) or len(maps['rows'])!=49:raise ArithmeticError('fixed49 maps incomplete')
    data={k:seed[k] for k in ('family','parameter','curve','generic_points')};data.update(protocol_hash=digest(p),maps_sha256=cert.hashed(D/'maps.json'),initial_state=state.record(),initial_dimension=27,centres=maps['centres'],metric_gram=maps['metric_gram'],charts=[],status='RUNNING',rank_lower_bound=27,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=D/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-27);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record();data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':final['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=final,rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('ADAPTIVE NEW27 SPECIALIZED',i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=p['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT';checkpoint(out,data)

def replay():
    p=protocol();data=cert.read(D/'result.json');maps=cert.read(D/'maps.json');cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(D/'maps.json') or data['initial_state']!=state.record() or any(data[k]!=seed[k] for k in ('family','parameter','curve','generic_points')) or data['centres']!=maps['centres']:raise ArithmeticError('adaptive initial binding differs')
    candidates=maps['sample'];expected=sorted(candidates,key=lambda c:(-c['metric_norm'],c['parity']))[:49]
    if [c['parity'] for c in candidates]!=masks(p) or maps['centres']!=expected or [m['centre'] for m in maps['rows']]!=expected:raise ArithmeticError('fixed sample and selected roster differ')
    for c in candidates:
        if len(c['representative'])!=27 or any((c['representative'][j]-(c['parity']>>j))%2 for j in range(27)):raise ArithmeticError('specialized parity differs')
    for i,row in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/row['archive_path']
        if row['index']!=i or row['centre']!=m['centre'] or cert.hashed(ap)!=row['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('adaptive archive differs')
        rep=m['centre']['representative']+[0]*(state.rank-27);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r=row['search']
        if r['height_bound']!=p['height'] or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('adaptive budget differs')
        points=backend.replay(search,m,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('adaptive orbit differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=row['state_key'] or state.rank!=row['rank_lower_bound'] or state.record()['state']['observations']!=row['admission_observations']:raise ArithmeticError('adaptive admissions differ')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound'] or (data['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(data['charts'])!=49):raise ArithmeticError('adaptive final state differs')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED ADAPTIVE NEW27 SPECIALIZED',len(data['charts']),'rank >=',state.rank,flush=True)

def launch():
    p=protocol()
    # Complete every geometry choice before observing any prospective point result.
    for index in range(6):
        configure(index);name='maps';path=D/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve prospective geometry supervision')
        s=supervise(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'prepare_new27_specialized_parity.sage'),'--index',str(index)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT);print(ROW['id'],name,s['outcome'],s['returncode'],flush=True)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('complete frozen maps required')
    for index in range(6):
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
