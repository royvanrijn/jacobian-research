#!/usr/bin/env python3
"""One fixed301-chart follow-up using the new27th MW16 direction oddly in every centre."""
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
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';PROOF=ART/'new_mw16_rank27_minimal_proof_v1.json'

IDS=('a1-fibration-01-052',)
def configure(identifier):
    global ID,D,SEED,CENSUS
    if identifier not in IDS:raise ValueError('only the independently certified new27 MW16 curve is in scope')
    ID=identifier;D=LOCAL/('mw16-new27-'+ID+'-newdirection-v1');SEED=LOCAL/('mw16-new26-'+ID+'-adaptive-v1')/'result.json';CENSUS=LOCAL/'prospective-mw16-h4096-v1'/ID.rsplit('-',1)[0]/'generic-census.json'

def sources():
    names=['memory_rank_certificate.py','followup_new_mw16_rank27_direction.py','prepare_new_mw16_rank27_direction.sage','prepare_fresh_r17_pari_batch.sage','fresh_r17_pari_batch.py','prospective_half_lattice_v2.sage','research_runtime/rotated_observation_state.py','research_runtime/cached_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve adaptive rank27 protocol')
    gate=ART/'million_height_mw16_portable_replay_v1.json'
    if cert.read(gate)['status']!='PASS':raise ArithmeticError('isolated MW16 proof gate required')
    seed,census=cert.read(SEED),cert.read(CENSUS);proof=next(r for r in cert.read(PROOF)['curves'] if r['id']==ID)
    if seed['status']!='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' or seed['rank_lower_bound']!=27 or proof['icarm_matches'] or proof['previous_matches'] or seed['curve']!=proof['discovery_curve']:raise ArithmeticError('certified new27 gate failed')
    points=seed['final_state']['state']['reductions']['points']
    if points!=proof['discovery_points']:raise ArithmeticError('certified ordered subgroup differs')
    p=proof['rank_certificate'];checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in points],[r['prime'] for r in p['signatures']],p['no_rational_2_torsion_prime'])
    pool=sorted(census['records'][1:],key=lambda r:(-cert.F(r['norm']),r['mask']))[:301]
    if len(pool)!=301 or len({r['mask'] for r in pool})!=301:raise ArithmeticError('generic pool differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.new-mw16-rank27-adaptive.v1','sources':sources(),'seed_path':str(SEED.relative_to(ROOT)),'seed_sha256':cert.hashed(SEED),'proof_sha256':cert.hashed(PROOF),'completed_mw16_gate_sha256':cert.hashed(gate),'generic_census_path':str(CENSUS.relative_to(ROOT)),'generic_census_sha256':cert.hashed(CENSUS),'generic_pool':[{'mask':r['mask'],'computed_generic_norm':r['norm']} for r in pool],'initial_rank':27,'generic_rank':16,'charts':301,'height':100000,'seconds_per_chart':10,'target_rank':28,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'geometry_wall_seconds':300,'worker_wall_seconds':3600,'rss_bytes':1610612736,'maximum_workers':1,'gate':'The first fixed301-centre MW16 follow-up produced a27th independently certified point on a1-fibration-01 at-1867/270. Its original43 and first301 centres all have even coefficient of this last new basis point. Test exactly301 parities with odd coefficient of that new direction, at the unchanged100000 height, from the exactly certified27-point subgroup. The earlier wave is terminal and its full cloud and minimal-model proof replay. One curve and one worker only.','centre_policy':'Retain the same301 generic-pool labels by computed norm and mask. Pair them with words1024+w, where w runs through ten-bit words including0 ordered by Hamming weight then integer. Every centre has bit26=1; all original and first-wave centres had bit26=0 when embedded in this ordered27-point basis. Compute384-bit numerical heights and CVP representatives, then freeze exact rational maps before searching. These301 pairings are not every quotient class or a specialized minimum theorem.','boundaries':'Only points discovered on this prospective curve enter centres. Freeze all rational maps before point searching. Preserve all histories and timeout outputs; no automatic continuation or rank upper bound. Stop at the first worker lower bound28 pending exact replay.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or cert.hashed(SEED)!=p['seed_sha256'] or cert.hashed(PROOF)!=p['proof_sha256'] or cert.hashed(CENSUS)!=p['generic_census_sha256']:raise ArithmeticError('frozen adaptive inputs differ')
    return p

def initial(cache):
    seed=cert.read(SEED);cache.store.import_snapshot(seed['arithmetic_facts']);state=MWState.from_record(seed['final_state'],cache=cache)
    if state.rank!=27:raise ArithmeticError('adaptive initial subgroup differs')
    return seed,state

def worker():
    p=protocol();maps=cert.read(D/'maps.json');out=D/'result.json';cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if out.exists():raise FileExistsError('preserve adaptive point attempt')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or maps['protocol_hash']!=digest(p) or len(maps['rows'])!=301:raise ArithmeticError('fixed301 maps incomplete')
    data={k:seed[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u')};data.update(protocol_hash=digest(p),maps_sha256=cert.hashed(D/'maps.json'),initial_state=state.record(),initial_dimension=27,centres=maps['centres'],metric_gram=maps['metric_gram'],charts=[],status='RUNNING',rank_lower_bound=27,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=D/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-27);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record();data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':final['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=final,rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('ADAPTIVE NEW27',i+1,r['status'],'rank',state.rank,flush=True)
        if state.rank>=p['target_rank']:data['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT';checkpoint(out,data)

def replay():
    p=protocol();data=cert.read(D/'result.json');maps=cert.read(D/'maps.json');cache=ReductionCache(MemoryFactStore());seed,state=initial(cache);model=tuple(map(cert.F,seed['curve']))
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(D/'maps.json') or data['initial_state']!=state.record() or any(data[k]!=seed[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u')) or data['centres']!=maps['centres']:raise ArithmeticError('adaptive initial binding differs')
    words=[1024+w for w in sorted(range(1024),key=lambda w:(w.bit_count(),w))];expected={(r['mask'],words[i%len(words)]) for i,r in enumerate(p['generic_pool'])}
    if len(maps['centres'])!=301 or {(c['generic_mask'],c['quotient_word']) for c in maps['centres']}!=expected or [m['centre'] for m in maps['rows']]!=maps['centres']:raise ArithmeticError('adaptive centre roster differs')
    for c in maps['centres']:
        if c['parity']!=c['generic_mask']|(c['quotient_word']<<16) or len(c['representative'])!=27 or any((c['representative'][j]-(c['parity']>>j))%2 for j in range(27)):raise ArithmeticError('adaptive parity differs')
    for i,row in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/row['archive_path']
        if row['index']!=i or row['centre']!=m['centre'] or cert.hashed(ap)!=row['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('adaptive archive differs')
        rep=m['centre']['representative']+[0]*(state.rank-27);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r=row['search']
        if r['height_bound']!=p['height'] or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('adaptive budget differs')
        points=backend.replay(search,m,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('adaptive orbit differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=row['state_key'] or state.rank!=row['rank_lower_bound'] or state.record()['state']['observations']!=row['admission_observations']:raise ArithmeticError('adaptive admissions differ')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound'] or (data['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT' and len(data['charts'])!=301):raise ArithmeticError('adaptive final state differs')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED ADAPTIVE NEW27',len(data['charts']),'rank >=',state.rank,flush=True)

def launch():
    p=protocol()
    for name,cmd,seconds in [('maps',['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'prepare_new_mw16_rank27_direction.sage'),'--id',ID],p['geometry_wall_seconds']),('worker',[sys.executable,str(Path(__file__).resolve()),'worker','--id',ID],p['worker_wall_seconds'])]:
        path=D/(name+'.supervisor.json')
        if path.exists():raise FileExistsError('preserve adaptive supervisor')
        s=supervise(cmd,limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=path,cwd=ROOT);print('NEW27 ADAPTIVE STAGE',name,s['outcome'],s['returncode'],flush=True)
        if s['outcome']!='completed' or s['returncode']!=0:return
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','worker','replay','launch']);p.add_argument('--id',required=True);a=p.parse_args();configure(a.id);globals()[a.stage]()
