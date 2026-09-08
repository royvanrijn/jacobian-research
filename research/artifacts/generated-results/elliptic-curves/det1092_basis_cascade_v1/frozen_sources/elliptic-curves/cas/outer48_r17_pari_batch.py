#!/usr/bin/env python3
"""Fixed eight per family outside the old H32768 parameter envelope."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import extend_outer131072_r17 as extension
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
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/outer48-r17-pari-v1';PARITY=ART/'r17_exact_maximum_parity_classes_v1.json';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def sources():
    names=['memory_rank_certificate.py','extend_outer131072_r17.py','verify_outer48_r17_pari_batch.py','outer48_r17_pari_batch.py','prepare_outer48_r17_pari_batch.sage','replay_outer48_geometry.py','audit_outer48_clouds_modl.py','stream_outer48_verification.py','prospective_half_lattice_v2.sage','compact_atlas_specialization.py','research_runtime/search_state.py','research_runtime/cached_observation_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [*(CAS/n for n in names),spec.ATLAS,PARITY,ART/'public_compact_parameter_heights_v1.json']}}


def expected_roster():
    selected=cert.read(extension.D/'result.json');rows=[]
    if selected['status']!='COMPLETE_FROZEN_TRACE_EXTENSION':raise ArithmeticError('terminal outer score extension required')
    expected=extension.selection(selected['rows'])
    if selected['selection']!=expected:raise ArithmeticError('fixed outer selection differs')
    indexed={(r['family'],r['retained_index']):r for r in selected['rows']}
    if len(indexed)!=len(selected['rows']):raise ArithmeticError('duplicate retained address')
    for family in sorted(expected):
        for index in expected[family]:
            r=indexed[family,index]
            if not 32768<max(abs(r['numerator']),r['denominator'])<=131072:raise ArithmeticError('outer height envelope differs')
            rows.append({'id':family+f"-{index:03}",'family':family,'parameter':r['parameter'],'retained_index':index,'arms':['outer131072_top8'],'combined_selection_units':r['combined_selection_units']})
    if len(rows)!=48 or len({r['id'] for r in rows})!=48:raise ArithmeticError('fixed48 roster differs')
    return rows

def freeze():
    if (D/'protocol.json').exists():raise FileExistsError('preserve outer48 protocol')
    trace=cert.read(extension.D/'replay.supervisor.json')
    control=cert.read(ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json')
    if trace['outcome']!='completed' or trace['returncode']!=0 or control['status']!='PASS' or control['completed_boxes']!=49 or control['rank_lower_bound']!=28:raise ArithmeticError('complete score replay and visibility control required')
    rows=expected_roster();masks={f['family']:[r['mask'] for r in f['classes']] for f in cert.read(PARITY)['families']}
    if sum(len(masks[r['family']]) for r in rows)!=2160:raise ArithmeticError('fixed generic chart count differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.outer48-r17-pari.v1','sources':sources(),'selection_sha256':cert.hashed(extension.D/'result.json'),'selection_protocol_sha256':cert.hashed(extension.D/'protocol.json'),'selection_replay_sha256':cert.hashed(extension.D/'replay.supervisor.json'),'control_sha256':cert.hashed(ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json'),'rows':rows,'generic_masks':masks,'height':125000,'seconds_per_chart':10,'admission_prime_bound':997,'stop_rank':28,'map_wall_seconds':120,'worker_wall_seconds':600,'map_roster_wall_seconds':1200,'point_batch_wall_seconds':15000,'geometry_before_search':True,'rss_bytes':1610612736,'maximum_workers':2,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'gate':'The exact public native11952 known29 control has compact parameter89074/31895 outside32768 and within131072; its generic17 initial detector reaches28 at point height125000. This motivates testing candidate incidence in a larger parameter envelope, without claiming a density law. The frozen stratified scan contains122433806 primitive addresses and retains6144 short-score rows;5798 outer rows receive the unchanged extended scores. Six benchmark curves passed direct trace and timing gates. The fixed48 selector excludes all741 previously measured equations and earlier selected equations by exact rational isomorphism.', 'selection':'Exactly eight per family after the frozen outer score replay, unchanged score through32749, good count, denominator and signed numerator. Validation32771through65521 never enters ordering. Equation-only exclusions precede all maps and point execution; no public ranks or points enter the selector or maps. No refill of failures or catalogue matches.','geometry':'All43 or49 exact generic maximum parity classes. Numerical384-bit heights choose representatives only; rational centres and maps are checked exactly.','boundaries':'Exactly48 attempts and at most2160 initial boxes, stop each curve on a provisional28-point independent subgroup pending replay. Each starts from17 generic sections only. Retain all checkpoints, failures and raw points. Catalogue and previous-equation comparison follows terminal batch replay. No full131072 parameter coverage, rank upper bound, exact rank, record or universal novelty.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['control_sha256']!=cert.hashed(ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json') or p['gp_sha256']!=cert.hashed(Path('/usr/bin/gp')) or p['sources']!=sources() or p['selection_sha256']!=cert.hashed(extension.D/'result.json') or p['selection_protocol_sha256']!=cert.hashed(extension.D/'protocol.json') or p['selection_replay_sha256']!=cert.hashed(extension.D/'replay.supervisor.json') or p['rows']!=expected_roster():raise ArithmeticError('frozen outer48 inputs or roster differ')
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
        final=state.record();data['charts'].append({'index':i,'centre':mapping['centre'],'archive_path':str(archive_path.relative_to(ROOT)),'archive_sha256':cert.hashed(archive_path),'search':record,'admission_compression':compression,'admission_observations':final['state']['observations'],'rank_lower_bound':state.rank,'state_key':state.key});data.update(rank_lower_bound=state.rank,final_state=final,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('OUTER48 R17 PARI',row['id'],i+1,'rank',state.rank,flush=True)
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
    if data['status'] not in ('COMPLETE_DECLARED_POINT_ATTEMPT','TARGET_REACHED_PENDING_REPLAY') or (data['status']=='TARGET_REACHED_PENDING_REPLAY' and state.rank<p['stop_rank']):raise ArithmeticError('terminal point status differs')
    memory_checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED OUTER48 R17',p['rows'][index]['id'],len(data['charts']),'rank >=',state.rank,flush=True)

def maps():
    p=protocol();path=D/'maps-ledger.json'
    if path.exists():raise FileExistsError('preserve full48 geometry roster')
    ledger={'status':'RUNNING','rows':[{**r,'status':'PENDING'} for r in p['rows']]};checkpoint(path,ledger)
    def one(index):
        row=p['rows'][index];folder=D/row['id']
        r=supervise([SAGE,str(CAS/'prepare_outer48_r17_pari_batch.sage'),'--index',str(index)],limits=Limits(p['map_wall_seconds'],p['rss_bytes']),log_path=folder/'maps.log',checkpoint_path=folder/'maps.supervisor.json',cwd=ROOT)
        return {**row,'status':'PASS' if r['outcome']=='completed' and r['returncode']==0 else 'FAILED_OR_CENSORED','supervision':r,'maps_sha256':cert.hashed(folder/'maps.json') if (folder/'maps.json').exists() else None}
    with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
        pending={pool.submit(one,i):i for i in range(len(p['rows']))}
        for f in as_completed(pending):
            i=pending[f];ledger['rows'][i]=f.result();checkpoint(path,ledger);print('OUTER48 FROZEN MAP',p['rows'][i]['id'],ledger['rows'][i]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES_OR_CENSORING';checkpoint(path,ledger)

def batch():
    p=protocol();path=D/'ledger.json';geometry=cert.read(D/'maps-ledger.json')
    if path.exists():raise FileExistsError('preserve outer48 point batch')
    if geometry['status']!='PASS' or [r['id'] for r in geometry['rows']]!=[r['id'] for r in p['rows']]:raise ArithmeticError('all48 geometry choices must precede searches')
    for row in geometry['rows']:
        if cert.hashed(D/row['id']/'maps.json')!=row['maps_sha256']:raise ArithmeticError('frozen map changed')
    ledger={'status':'RUNNING','rows':[{**r,'status':'PENDING'} for r in p['rows']]};checkpoint(path,ledger)
    def one(index):
        row=p['rows'][index];folder=D/row['id'];r={**row,'status':'POINT_FAILED_OR_CENSORED','map_supervision':geometry['rows'][index]['supervision']}
        search=supervise(['/usr/bin/python3',str(Path(__file__).resolve()),'worker','--index',str(index)],limits=Limits(p['worker_wall_seconds'],p['rss_bytes']),log_path=folder/'worker.log',checkpoint_path=folder/'worker.supervisor.json',cwd=ROOT);r['point_supervision']=search
        if (folder/'result.json').exists():
            data=cert.read(folder/'result.json');r.update(rank_lower_bound=data['rank_lower_bound'],charts=len(data['charts']),result_sha256=cert.hashed(folder/'result.json'))
            if search['outcome']=='completed' and search['returncode']==0:r['status']=data['status']
        return r
    with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
        pending={pool.submit(one,i):i for i in range(len(p['rows']))}
        for f in as_completed(pending):
            i=pending[f];r=f.result();ledger['rows'][i]=r;checkpoint(path,ledger);print('OUTER48 POINT ATTEMPT',r['id'],r['status'],r.get('rank_lower_bound'),flush=True)
    ledger['status']='COMPLETE_FIXED_BATCH_ATTEMPTS';checkpoint(path,ledger)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['freeze','maps','batch','worker','replay']);a.add_argument('--index',type=int);v=a.parse_args();globals()[v.stage](v.index) if v.stage in ('worker','replay') else globals()[v.stage]()
