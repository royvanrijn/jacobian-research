#!/usr/bin/env python3
"""A fixed49-chart specialized-parity control from28 already recovered points."""
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
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'native28-specialized-parity-control-v1';SEED=LOCAL/'native11952-height125-control-v1/125000/result.json';VERIFY=SEED.parent/'verification.json'

def sources():
    names=['memory_rank_certificate.py','native28_specialized_parity_control.py','prepare_native28_specialized_parity.sage','prepare_fresh_r17_pari_batch.sage','fresh_r17_pari_batch.py','prospective_half_lattice_v2.sage','research_runtime/rotated_observation_state.py','research_runtime/cached_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve specialized-parity control')
    seed=cert.read(SEED);v=cert.read(VERIFY)
    if seed['status']!='COMPLETE_DECLARED_HEIGHT_ARM' or seed['rank_lower_bound']!=28 or v['status']!='PASS' or v['completed_boxes']!=49 or len(seed['final_state']['state']['reductions']['points'])!=28:raise ArithmeticError('completed independent known28 gate required')
    gate=ART/'local_feature_portable_replay_v1.json'
    if cert.read(gate)['status']!='PASS':raise ArithmeticError('complete local-feature audit required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native28-specialized-parity-control.v1','sources':sources(),'seed_path':str(SEED.relative_to(ROOT)),'seed_sha256':cert.hashed(SEED),'verification_sha256':cert.hashed(VERIFY),'gate_sha256':cert.hashed(gate),'initial_rank':28,'generic_rank':17,'sample_size':2048,'sample_domain':'native28-specialized-parity-v1','charts':49,'height':125000,'seconds_per_chart':10,'target_rank':29,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'geometry_wall_seconds':300,'worker_wall_seconds':900,'rss_bytes':2147483648,'maximum_workers':1,'gate':'The current301-label policies pair a selected generic label with one quotient word, leaving most specialized parity classes unmeasured. The known native11952 control previously recovered28 from17 generic seeds at125000 but required a retrospectively chosen million-height chart for its29th direction. Test a different finite centre policy at125000 using only its28 previously recovered points. The unrecovered public point and all oracle words/visibility minima stay outside preparation and execution.','centre_policy':'Take the first2048 distinct28-bit SHA256-derived masks with nonzero quotient above bit16, digest of [fixed sample domain,index] modulo2^28. Compute384-bit numerical canonical heights and round at10^6. Use a unimodular LLL basis change only for numerical CVP, then transport every representative exactly back to the original28-point subgroup and check its original parity and rounded norm. Select49 largest computed norms, mask ties. Freeze all rational maps before point search. Neither numerical coset optimality nor complete specialized covering is claimed.','boundaries':'This is recovery on a previously known29 curve, not a new curve, independent population validation or a new rank record. Its choice is informed by prior control history. A recovery measures this fixed policy only. Preserve partial/censored attempts and stop upon29 pending exact replay; no automatic larger sample, height increase or prospective escalation.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or cert.hashed(SEED)!=p['seed_sha256'] or cert.hashed(VERIFY)!=p['verification_sha256']:raise ArithmeticError('frozen specialized control inputs differ')
    return p

def masks(p):
    result=[];i=0
    while len(result)<p['sample_size']:
        m=int(digest([p['sample_domain'],i]),16)%(1<<28);i+=1
        if m>>17 and m not in result:result.append(m)
    return result

def initial(cache):
    seed=cert.read(SEED);cache.store.import_snapshot(seed['arithmetic_facts']);state=MWState.from_record(seed['final_state'],cache=cache)
    if state.rank!=28:raise ArithmeticError('adaptive initial subgroup differs')
    return seed,state

def worker():
    p=protocol();maps=cert.read(D/'maps.json');out=D/'result.json';cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if out.exists():raise FileExistsError('preserve adaptive point attempt')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or maps['protocol_hash']!=digest(p) or len(maps['rows'])!=49:raise ArithmeticError('fixed49 maps incomplete')
    data={k:seed[k] for k in ('family','parameter','curve','generic_points')};data.update(protocol_hash=digest(p),maps_sha256=cert.hashed(D/'maps.json'),initial_state=state.record(),initial_dimension=28,centres=maps['centres'],metric_gram=maps['metric_gram'],charts=[],status='RUNNING',rank_lower_bound=28,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=D/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-28);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record();data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':final['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=final,rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('ADAPTIVE KNOWN28 SPECIALIZED',i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=p['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT';checkpoint(out,data)

def replay():
    p=protocol();data=cert.read(D/'result.json');maps=cert.read(D/'maps.json');cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(D/'maps.json') or data['initial_state']!=state.record() or any(data[k]!=seed[k] for k in ('family','parameter','curve','generic_points')) or data['centres']!=maps['centres']:raise ArithmeticError('adaptive initial binding differs')
    candidates=maps['sample'];expected=sorted(candidates,key=lambda c:(-c['metric_norm'],c['parity']))[:49]
    if [c['parity'] for c in candidates]!=masks(p) or maps['centres']!=expected or [m['centre'] for m in maps['rows']]!=expected:raise ArithmeticError('fixed sample and selected roster differ')
    for c in candidates:
        if len(c['representative'])!=28 or any((c['representative'][j]-(c['parity']>>j))%2 for j in range(28)):raise ArithmeticError('specialized parity differs')
    for i,row in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/row['archive_path']
        if row['index']!=i or row['centre']!=m['centre'] or cert.hashed(ap)!=row['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('adaptive archive differs')
        rep=m['centre']['representative']+[0]*(state.rank-28);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r=row['search']
        if r['height_bound']!=p['height'] or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('adaptive budget differs')
        points=backend.replay(search,m,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('adaptive orbit differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=row['state_key'] or state.rank!=row['rank_lower_bound'] or state.record()['state']['observations']!=row['admission_observations']:raise ArithmeticError('adaptive admissions differ')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound'] or (data['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(data['charts'])!=49):raise ArithmeticError('adaptive final state differs')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED ADAPTIVE KNOWN28 SPECIALIZED',len(data['charts']),'rank >=',state.rank,flush=True)

def launch():
    p=protocol()
    for name,cmd,seconds in [('maps',['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'prepare_native28_specialized_parity.sage')],p['geometry_wall_seconds']),('worker',[sys.executable,str(Path(__file__).resolve()),'worker'],p['worker_wall_seconds'])]:
        path=D/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve adaptive supervisor')
        s=supervise(cmd,limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT);print('KNOWN28 SPECIALIZED ADAPTIVE STAGE',name,s['outcome'],s['returncode'],flush=True)
        if s['outcome']!='completed' or s['returncode']!=0:return
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','worker','replay','launch']);a=p.parse_args();globals()[a.stage]()
