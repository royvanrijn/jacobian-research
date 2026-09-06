#!/usr/bin/env python3
"""Fixed top four per family from the wider retention extended-prime policy."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import select_higher13_scaled_candidates as extension
import pari_pointed_backend as backend
from memory_rank_certificate import checked_rank as memory_checked_rank
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run as supervise,Limits
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.preloaded_prime_state import preload
from research_runtime.rotated_observation_state import rotate
from research_runtime.pointed_orbit_compression import compress
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/scaled13_24-r17-pari-v1';PARITY=ART/'r17_exact_maximum_parity_classes_v1.json';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def sources():
    names=['memory_rank_certificate.py','select_higher13_scaled_candidates.py','verify_scaled13_24_r17_pari_batch.py','scaled13_24_r17_pari_batch.py','prepare_scaled13_24_r17_pari_batch.sage','prospective_half_lattice_v2.sage','compact_atlas_specialization.py','research_runtime/search_state.py','research_runtime/cached_observation_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [*(CAS/n for n in names),spec.ATLAS,PARITY]}}


def expected_roster():
    selected=cert.read(extension.D/'result.json');rows=[]
    if selected['status']!='COMPLETE_FROZEN_CONDITIONAL_SELECTION':raise ArithmeticError('terminal score extension required')
    for f in sorted(selected['selection']):
        available=[r for r in selected['rows'] if r['family']==f]
        available.sort(key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
        if [r['retained_index'] for r in available[:4]]!=selected['selection'][f]:raise ArithmeticError('fixed selection differs')
        for r in available[:4]:rows.append({'id':f+f"-{r['retained_index']:03}",'family':f,'parameter':r['parameter'],'retained_index':r['retained_index'],'arms':['exact13_scaled_s1_top_four'],'combined_selection_units':r['combined_selection_units']})
    if len(rows)!=24 or len({r['id'] for r in rows})!=24:raise ArithmeticError('fixed24 roster differs')
    return rows

def freeze():
    if (D/'protocol.json').exists():raise FileExistsError('preserve scaled13_24 protocol')
    trace=cert.read(extension.D/'replay.supervisor.json')
    control=cert.read(ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json')
    if trace['outcome']!='completed' or trace['returncode']!=0 or control['status']!='PASS' or control['completed_boxes']!=49 or control['rank_lower_bound']!=28:raise ArithmeticError('complete score replay and visibility control required')
    rows=expected_roster();masks={f['family']:[r['mask'] for r in f['classes']] for f in cert.read(PARITY)['families']}
    if sum(len(masks[r['family']]) for r in rows)!=1080:raise ArithmeticError('fixed generic chart count differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.scaled13_24-r17-pari.v1','sources':sources(),'selection_sha256':cert.hashed(extension.D/'result.json'),'selection_protocol_sha256':cert.hashed(extension.D/'protocol.json'),'selection_replay_sha256':cert.hashed(extension.D/'replay.supervisor.json'),'control_sha256':cert.hashed(ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json'),'rows':rows,'generic_masks':masks,'height':125000,'seconds_per_chart':10,'admission_prime_bound':997,'stop_rank':28,'map_wall_seconds':120,'worker_wall_seconds':600,'rss_bytes':1610612736,'maximum_workers':2,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'gate':'The original6144 S1 finalists have complete retained-score and trace replay. The exact local-model audit and all55 product-first experiment checks now pass. A separate frozen selector restricts to exact13-scaled models and uniformly excludes addresses in all three completed higher-height point protocols, leaving454 eligible rows. Four per family retain the original S1 score order. This tests a measured model-size difference, without asserting that it causes higher rank or better point visibility. The prior identical49-map native control reached28 at125000 versus27 at100000; this is a finite visibility and cost gate, not a predictor for these candidates.','selection':'Exactly four per family from the454 eligible13-scaled addresses using the unchanged S1 score through32749, good count, denominator, numerator. Validation32771..65521 is excluded even from ties. Public equations, points, ranks and catalogue labels are absent from selection and execution. No refill of repeats or failures.','geometry':'All43 or49 exact generic maximum parity classes. Numerical384-bit heights choose representatives only; rational centres and maps are checked exactly.','boundaries':'Exactly24 attempts and at most1080 initial boxes, stop each curve on a provisional28-point independent subgroup pending replay. Each starts from17 generic sections only. Retain all checkpoints, failures and raw points. Catalogue and previous-equation comparison follows terminal batch replay. No complete32768 coverage, rank upper bound, exact rank, record or universal novelty.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['selection_sha256']!=cert.hashed(extension.D/'result.json') or p['selection_protocol_sha256']!=cert.hashed(extension.D/'protocol.json') or p['selection_replay_sha256']!=cert.hashed(extension.D/'replay.supervisor.json') or p['rows']!=expected_roster():raise ArithmeticError('frozen scaled13_24 inputs or roster differ')
    return p

def initial(data,cache):
    model=tuple(map(cert.F,data['curve']));points=tuple(tuple(map(cert.F,P)) for P in data['generic_points']);s=raw_state(model,points,cache=cache,prime_bound=1000);s=MWState.from_record(s.record(),cache=cache)
    if s.rank!=17:raise ArithmeticError('generic17 certificate incomplete')
    return preload(s,cache,997)

def worker(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];maps=cert.read(folder/'maps.json');out=folder/'result.json'
    if out.exists():raise FileExistsError('preserve point attempt')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or maps['protocol_hash']!=digest(p):raise ArithmeticError('maps not complete/bound')
    cache=ReductionCache(MemoryFactStore());state,bank=initial(maps,cache);model=tuple(map(cert.F,maps['curve']));data={k:maps[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u','centres','metric_gram')};data.update(protocol_hash=digest(p),maps_sha256=cert.hashed(folder/'maps.json'),initial_state=state.record(),prime_bank=bank,charts=[],status='RUNNING',rank_lower_bound=17,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    for i,mapping in enumerate(maps['rows']):
        state,archive=rotate(state);archive_path=folder/'states'/f'{i:03}.json';checkpoint(archive_path,archive);rep=mapping['centre']['representative']+[0]*(state.rank-17);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=mapping['coordinate_policy']);record,points=backend.execute(search,mapping,p['height'],p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record();data['charts'].append({'index':i,'centre':mapping['centre'],'archive_path':str(archive_path.relative_to(ROOT)),'archive_sha256':cert.hashed(archive_path),'search':record,'admission_compression':compression,'admission_observations':final['state']['observations'],'rank_lower_bound':state.rank,'state_key':state.key});data.update(rank_lower_bound=state.rank,final_state=final,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('SCALED13_24 R17 PARI',row['id'],i+1,'rank',state.rank,flush=True)
        if state.rank>=p['stop_rank']:data['status']='TARGET_REACHED_PENDING_REPLAY';checkpoint(out,data);return
    data['status']='COMPLETE_DECLARED_POINT_ATTEMPT';checkpoint(out,data)

def replay(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];maps=cert.read(folder/'maps.json');data=cert.read(folder/'result.json');cache=ReductionCache(MemoryFactStore());state,bank=initial(maps,cache)
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(folder/'maps.json') or data['initial_state']!=state.record() or data['prime_bank']!=bank:raise ArithmeticError('initial point binding differs')
    if any(data[k]!=maps[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u','centres','metric_gram')):raise ArithmeticError('maps metadata differs')
    f=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']==data['family']);original,points=spec.specialize(f,data['parameter']);u=cert.F(data['family_to_curve_scale_u']);model=tuple(map(cert.F,data['curve']));generic=tuple(tuple(map(cert.F,P)) for P in data['generic_points'])
    if model!=(0,0,0,original[3]/u**4,original[4]/u**6) or generic!=tuple((x/u**2,y/u**3) for x,y in points):raise ArithmeticError('specialization transport differs')
    if sorted(c['mask'] for c in maps['centres'])!=sorted(p['generic_masks'][data['family']]) or [m['centre'] for m in maps['rows']]!=maps['centres']:raise ArithmeticError('complete mask roster differs')
    for c in maps['centres']:
        if len(c['representative'])!=17 or any((c['representative'][j]-(c['mask']>>j))%2 for j in range(17)):raise ArithmeticError('centre parity differs')
    for i,row in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);path=ROOT/row['archive_path']
        if row['index']!=i or row['centre']!=m['centre'] or cert.hashed(path)!=row['archive_sha256'] or cert.read(path)!=archive:raise ArithmeticError('archive or chart roster differs')
        rep=m['centre']['representative']+[0]*(state.rank-17);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r=row['search']
        if r['height_bound']!=p['height'] or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('chart budget differs')
        points=backend.replay(search,m,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('orbit compression differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if row['state_key']!=state.key or row['rank_lower_bound']!=state.rank or row['admission_observations']!=state.record()['state']['observations']:raise ArithmeticError('admission replay differs')
    if data['final_state']!=state.record() or data['rank_lower_bound']!=state.rank or (data['status']=='COMPLETE_DECLARED_POINT_ATTEMPT' and len(data['charts'])!=len(maps['rows'])):raise ArithmeticError('final point state differs')
    memory_checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED SCALED13_24 R17',p['rows'][index]['id'],len(data['charts']),'rank >=',state.rank,flush=True)

def batch():
    p=protocol();path=D/'ledger.json'
    if path.exists():raise FileExistsError('preserve paired batch')
    ledger={'status':'RUNNING','rows':[{**r,'status':'PENDING'} for r in p['rows']]};checkpoint(path,ledger)
    def one(index):
        row=p['rows'][index];folder=D/row['id'];mapping=supervise([SAGE,str(CAS/'prepare_scaled13_24_r17_pari_batch.sage'),'--index',str(index)],limits=Limits(p['map_wall_seconds'],p['rss_bytes']),log_path=folder/'maps.log',checkpoint_path=folder/'maps.supervisor.json',cwd=ROOT)
        r={**row,'status':'MAP_FAILED_OR_CENSORED','map_supervision':mapping}
        if mapping['outcome']!='completed' or mapping['returncode']!=0:return r
        search=supervise(['/usr/bin/python3',str(Path(__file__).resolve()),'worker','--index',str(index)],limits=Limits(p['worker_wall_seconds'],p['rss_bytes']),log_path=folder/'worker.log',checkpoint_path=folder/'worker.supervisor.json',cwd=ROOT);r.update(status='POINT_FAILED_OR_CENSORED',point_supervision=search)
        if (folder/'result.json').exists():
            d=cert.read(folder/'result.json');r.update(rank_lower_bound=d['rank_lower_bound'],charts=len(d['charts']),result_sha256=cert.hashed(folder/'result.json'))
            if search['outcome']=='completed' and search['returncode']==0:r['status']=d['status']
        return r
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(one,i):i for i in range(len(p['rows']))}
        for future in as_completed(futures):
            i=futures[future];r=future.result();ledger['rows'][i]=r;checkpoint(path,ledger);print('SCALED13_24 R17 ATTEMPT',r['id'],r['status'],r.get('rank_lower_bound'),flush=True)
    ledger['status']='COMPLETE_FIXED_BATCH_ATTEMPTS';checkpoint(path,ledger)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['freeze','batch','worker','replay']);a.add_argument('--index',type=int);v=a.parse_args();globals()[v.stage](v.index) if v.stage in ('worker','replay') else globals()[v.stage]()
