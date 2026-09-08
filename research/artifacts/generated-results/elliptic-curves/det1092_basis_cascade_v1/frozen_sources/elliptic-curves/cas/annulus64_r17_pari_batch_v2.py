#!/usr/bin/env python3
"""Fixed64 full11952 point attempts after complete retained-score verification."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import score_11952_new_annulus_v2 as extension
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
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/annulus64-r17-pari-v2';PARITY=ART/'r17_exact_maximum_parity_classes_v1.json';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def sources():
    names=['memory_rank_certificate.py','score_11952_new_annulus_v2.py','verify_annulus64_r17_pari_batch_v2.py','annulus64_r17_pari_batch_v2.py','prepare_annulus64_r17_pari_batch_v2.sage','replay_annulus64_geometry_v2.py','audit_annulus64_clouds_modl_v2.py','prospective_half_lattice_v2.sage','compact_atlas_specialization.py','research_runtime/search_state.py','research_runtime/cached_observation_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [*(CAS/n for n in names),spec.ATLAS,PARITY]}}


CONTROL=extension.D/'controller'

def expected_roster():
    selected=cert.read(extension.OUT);rows=[]
    if selected['status']!='PASS_FROZEN64_SELECTION' or len(selected['selected'])!=64:raise ArithmeticError('terminal fixed64 selection required')
    for r in selected['selected']:
        if r['family']!='11952' or r['id']!=f"11952-{r['retained_index']:07}" or not 131072<max(abs(r['numerator']),r['denominator'])<=524288 or cert.F(r['numerator'],r['denominator'])!=cert.F(r['parameter']):raise ArithmeticError('full11952 retained address differs')
        rows.append({k:r[k] for k in ('id','family','parameter','retained_index','combined_selection_units','combined_late_units')})
        rows[-1]['arms']=['new11952_annulus_fresh4096_top64']
    if len({r['id'] for r in rows})!=64:raise ArithmeticError('fixed64 unique source roster required')
    return rows

def freeze():
    if (D/'protocol.json').exists():raise FileExistsError('preserve full11952 point protocol')
    ledger=cert.read(CONTROL/'ledger.json');trace=cert.read(CONTROL/'validation-check.supervisor.json');validation=cert.read(extension.D/'fresh-validation.json');selection_check=cert.read(CONTROL/'selection-check.supervisor.json')
    control=cert.read(ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json')
    if ledger['status']!='PASS' or len(ledger['rows'])!=10 or trace['outcome']!='completed' or trace['returncode']!=0 or selection_check['outcome']!='completed' or selection_check['returncode']!=0 or validation['status']!='PASS' or len(validation['rows'])!=64 or control['status']!='PASS' or control['completed_boxes']!=49 or control['rank_lower_bound']!=28:raise ArithmeticError('complete4096 score pipeline, fresh validation and point-visibility gate required')
    rows=expected_roster();masks={f['family']:[r['mask'] for r in f['classes']] for f in cert.read(PARITY)['families']}
    if sum(len(masks[r['family']]) for r in rows)!=3136:raise ArithmeticError('fixed3136 chart roster differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.annulus64-r17-pari.v1','sources':sources(),'selection_sha256':cert.hashed(extension.OUT),'selection_protocol_sha256':cert.hashed(extension.D/'protocol.json'),'selection_replay_sha256':cert.hashed(CONTROL/'selection-check.supervisor.json'),'validation_sha256':cert.hashed(extension.D/'fresh-validation.json'),'validation_replay_sha256':cert.hashed(CONTROL/'validation-check.supervisor.json'),'score_controller_sha256':cert.hashed(CONTROL/'ledger.json'),'control_sha256':cert.hashed(ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json'),'rows':rows,'generic_masks':masks,'height':125000,'seconds_per_chart':10,'admission_prime_bound':997,'stop_rank':28,'map_wall_seconds':120,'worker_wall_seconds':600,'map_roster_wall_seconds':1800,'point_batch_wall_seconds':20000,'geometry_before_search':True,'rss_bytes':1610612736,'maximum_workers':2,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'gate':'The separately frozen4096 new equations all pass scalar4099..65521 traces, prior cached-score agreement and8192 direct character sums. The64 finalists are selected through65521 and then pass wholly disjoint65537..131071 validation. Exact exclusions remove the frozen977 earlier cohort and593 pinned catalogue equations before the4096 roster; selected equations are distinct. The native generic17-only control recovers28 at125000. Two completed11952 cohorts supplied new27-point fibres inside the former complete131072 square. This finite trial tests fresh131072<H<=524288 addresses, without assuming a rank-density or sensitivity law.','selection':'The64 frozen Q-isomorphism-distinct nonsingular finalists from the fresh4096 second-band selector. Combined S1 through65521, total good count, denominator and signed numerator determine ordering. Fresh validation and public ranks/points never enter ordering or maps. Catalogue equations are exclusions only. The4096 roster is an explicit truncation of the32768 retained new-annulus addresses; no refill or further parameter scan.','geometry':'Every one of the49 generic maximum-parity classes on each curve. Numerical384-bit heights and rounded CVP choose representatives; exact parities, rational centres and quartic transports are checked. Every one of the64 map files precedes all point searches.','boundaries':'Exactly64 point attempts and at most3136 initial boxes, height125000 and ten seconds each, generic17 seeds and a per-curve provisional28 stop pending exact proof. No adaptive wave, automatic budget extension, exact rank, upper bound, saturation, point absence, score optimality or universal novelty.'})

def protocol():
    p=cert.read(D/'protocol.json')
    bindings={'selection_sha256':extension.OUT,'selection_protocol_sha256':extension.D/'protocol.json','selection_replay_sha256':CONTROL/'selection-check.supervisor.json','validation_sha256':extension.D/'fresh-validation.json','validation_replay_sha256':CONTROL/'validation-check.supervisor.json','score_controller_sha256':CONTROL/'ledger.json','control_sha256':ROOT/'artifacts/local/elliptic-curves/native11952-height125-control-v1/125000/verification.json','gp_sha256':Path('/usr/bin/gp')}
    if p['sources']!=sources() or any(p[k]!=cert.hashed(v) for k,v in bindings.items()) or p['rows']!=expected_roster():raise ArithmeticError('frozen full11952 point inputs or roster differ')
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
        final=state.record();data['charts'].append({'index':i,'centre':mapping['centre'],'archive_path':str(archive_path.relative_to(ROOT)),'archive_sha256':cert.hashed(archive_path),'search':record,'admission_compression':compression,'admission_observations':final['state']['observations'],'rank_lower_bound':state.rank,'state_key':state.key});data.update(rank_lower_bound=state.rank,final_state=final,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('ANNULUS64 R17 PARI',row['id'],i+1,'rank',state.rank,flush=True)
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
    memory_checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED ANNULUS64 R17',p['rows'][index]['id'],len(data['charts']),'rank >=',state.rank,flush=True)

def maps():
    p=protocol();path=D/'maps-ledger.json'
    if path.exists():raise FileExistsError('preserve full64 geometry roster')
    ledger={'status':'RUNNING','rows':[{**r,'status':'PENDING'} for r in p['rows']]};checkpoint(path,ledger)
    def one(index):
        row=p['rows'][index];folder=D/row['id']
        r=supervise([SAGE,str(CAS/'prepare_annulus64_r17_pari_batch_v2.sage'),'--index',str(index)],limits=Limits(p['map_wall_seconds'],p['rss_bytes']),log_path=folder/'maps.log',checkpoint_path=folder/'maps.supervisor.json',cwd=ROOT)
        return {**row,'status':'PASS' if r['outcome']=='completed' and r['returncode']==0 else 'FAILED_OR_CENSORED','supervision':r,'maps_sha256':cert.hashed(folder/'maps.json') if (folder/'maps.json').exists() else None}
    with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
        pending={pool.submit(one,i):i for i in range(len(p['rows']))}
        for f in as_completed(pending):
            i=pending[f];ledger['rows'][i]=f.result();checkpoint(path,ledger);print('ANNULUS64 FROZEN MAP',p['rows'][i]['id'],ledger['rows'][i]['status'],flush=True)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in ledger['rows']) else 'COMPLETE_WITH_FAILURES_OR_CENSORING';checkpoint(path,ledger)

def batch():
    p=protocol();path=D/'ledger.json';geometry=cert.read(D/'maps-ledger.json')
    if path.exists():raise FileExistsError('preserve annulus64 point batch')
    if geometry['status']!='PASS' or [r['id'] for r in geometry['rows']]!=[r['id'] for r in p['rows']]:raise ArithmeticError('all64 geometry choices must precede searches')
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
            i=pending[f];r=f.result();ledger['rows'][i]=r;checkpoint(path,ledger);print('ANNULUS64 POINT ATTEMPT',r['id'],r['status'],r.get('rank_lower_bound'),flush=True)
    ledger['status']='COMPLETE_FIXED_BATCH_ATTEMPTS';checkpoint(path,ledger)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['freeze','maps','batch','worker','replay']);a.add_argument('--index',type=int);v=a.parse_args();globals()[v.stage](v.index) if v.stage in ('worker','replay') else globals()[v.stage]()
