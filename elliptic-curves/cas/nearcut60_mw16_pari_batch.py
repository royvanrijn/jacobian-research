#!/usr/bin/env python3
"""Fixed exposure on highest unsearched retained near-finalists."""
import argparse,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
import select_retained_mw16_nearcutoff as extension
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
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/nearcut60-mw16-pari-v1';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def sources():
    names=['memory_rank_certificate.py','nearcut60_mw16_pari_batch.py','prepare_nearcut60_mw16_pari_batch.sage','prospective_half_lattice_v2.sage','compact_mw16_specialization.py','research_runtime/search_state.py','research_runtime/cached_observation_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [*(CAS/n for n in names),spec.ATLAS,Path(extension.__file__),*census_paths()]}}

PARENT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h4096-v1'
def census_paths():
    return [PARENT/f'a1-fibration-{i:02}'/'generic-census.json' for i in range(1,6)]
def masks():
    result={}
    for path in census_paths():
        d=cert.read(path)
        if d['status']!='COMPLETE_DECLARED_CENSUS' or len(d['records'])!=65536 or len(d['selected'])!=43:raise ArithmeticError('fixed recorded MW16 census differs')
        result[d['family']]=[r['mask'] for r in d['selected']]
    return result

def expected_roster():
    d=cert.read(extension.OUT)
    if d['status']!='PASS_FROZEN60_RETAINED_SELECTION' or len(d['selected'])!=60:
        raise ArithmeticError('complete frozen retained60 roster required')
    return d['selected']

def freeze():
    if (D/'protocol.json').exists():raise FileExistsError('preserve retained point protocol')
    extension.completion_gate();extension.protocol()
    gate=cert.read(extension.D/'controller/ledger.json')
    if gate['status']!='PASS_FROZEN60_RETAINED_SELECTION':raise ArithmeticError('retained selection replay required')
    rows=expected_roster();roster=masks()
    if sum(len(roster[r['family']]) for r in rows)!=2580:raise ArithmeticError('fixed2580 maps differ')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.nearcut60-mw16-pari.v1',
        'sources':sources(),'selection_sha256':cert.hashed(extension.OUT),
        'selection_protocol_sha256':cert.hashed(extension.D/'protocol.json'),
        'selection_ledger_sha256':cert.hashed(extension.D/'controller/ledger.json'),
        'rows':rows,'generic_masks':roster,'height':125000,'seconds_per_chart':10,
        'admission_prime_bound':997,'stop_rank':None,'map_wall_seconds':120,
        'baseline_wall_seconds':120,'worker_wall_seconds':600,'rss_bytes':1610612736,
        'maximum_workers':2,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),
        'scope':'Existing corrected late scores only; six highest unsearched near-finalists per family/band, from ranks7..32. All60 maps and exactly independent generic16 baselines before points. Exactly43 generic parity charts per curve; identical fixed exposure, no rank stop, adaptive wave, retry or refill. Retain terminal partial prefixes and independent exposure checkpoints. Certify rank gains separately from completion and measured computation. No validation-prime selection or new parameter sweep.'})

def protocol():
    p=cert.read(D/'protocol.json')
    bindings=[('selection_sha256',extension.OUT),('selection_protocol_sha256',extension.D/'protocol.json'),('selection_ledger_sha256',extension.D/'controller/ledger.json')]
    if p['sources']!=sources() or any(p[k]!=cert.hashed(v) for k,v in bindings) or p['rows']!=expected_roster() or p['generic_masks']!=masks():raise ArithmeticError('frozen retained inputs differ')
    return p

def initial(data,cache):
    model=tuple(map(cert.F,data['curve']));points=tuple(tuple(map(cert.F,P)) for P in data['generic_points']);s=raw_state(model,points,cache=cache,prime_bound=1000);s=MWState.from_record(s.record(),cache=cache)
    if s.rank!=16:raise ArithmeticError('generic16 certificate incomplete')
    return preload(s,cache,997)

def worker(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];maps=cert.read(folder/'maps.json');out=folder/'result.json'
    if out.exists():raise FileExistsError('preserve point attempt')
    if maps['status']!='COMPLETE_DECLARED_MAPS' or maps['protocol_hash']!=digest(p):raise ArithmeticError('maps not complete/bound')
    cache=ReductionCache(MemoryFactStore());state,bank=initial(maps,cache);model=tuple(map(cert.F,maps['curve']));data={k:maps[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u','centres','metric_gram')};data.update(protocol_hash=digest(p),maps_sha256=cert.hashed(folder/'maps.json'),initial_state=state.record(),prime_bank=bank,charts=[],status='RUNNING',rank_lower_bound=16,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(out,data)
    exposure={'protocol_hash':digest(p),'maps_sha256':cert.hashed(folder/'maps.json'),'charts':[]}
    checkpoint(folder/'exposure.json',exposure)
    for i,mapping in enumerate(maps['rows']):
        state,archive=rotate(state);archive_path=folder/'states'/f'{i:03}.json';checkpoint(archive_path,archive);rep=mapping['centre']['representative']+[0]*(state.rank-16);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=mapping['coordinate_policy']);exposure['charts'].append({'index':i,'centre':mapping['centre'],'status':'STARTED','started_monotonic':time.monotonic()});checkpoint(folder/'exposure.json',exposure)
        record,points=backend.execute(search,mapping,p['height'],p['seconds_per_chart'],p['gp_sha256'])
        exposure['charts'][-1].update(status='RETURNED',search=record);checkpoint(folder/'exposure.json',exposure)
        compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        final=state.record();data['charts'].append({'index':i,'centre':mapping['centre'],'archive_path':str(archive_path.relative_to(ROOT)),'archive_sha256':cert.hashed(archive_path),'search':record,'admission_compression':compression,'admission_observations':final['state']['observations'],'rank_lower_bound':state.rank,'state_key':state.key});data.update(rank_lower_bound=state.rank,final_state=final,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('NEARCUT60 MW16 PARI',row['id'],i+1,'rank',state.rank,flush=True)
    data['status']='COMPLETE_DECLARED_POINT_ATTEMPT';checkpoint(out,data)

def replay(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];maps=cert.read(folder/'maps.json');data=cert.read(folder/'result.json');cache=ReductionCache(MemoryFactStore());state,bank=initial(maps,cache)
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(folder/'maps.json') or data['initial_state']!=state.record() or data['prime_bank']!=bank:raise ArithmeticError('initial point binding differs')
    if any(data[k]!=maps[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u','centres','metric_gram')):raise ArithmeticError('maps metadata differs')
    f=next(f for f in cert.read(spec.ATLAS)['families'] if f['fibration_id']==data['family']);original,points=spec.specialize(f,data['parameter']);u=cert.F(data['family_to_curve_scale_u']);model=tuple(map(cert.F,data['curve']));generic=tuple(tuple(map(cert.F,P)) for P in data['generic_points'])
    if model!=(0,0,0,original[3]/u**4,original[4]/u**6) or generic!=tuple((x/u**2,y/u**3) for x,y in points):raise ArithmeticError('specialization transport differs')
    if sorted(c['mask'] for c in maps['centres'])!=sorted(p['generic_masks'][data['family']]) or [m['centre'] for m in maps['rows']]!=maps['centres']:raise ArithmeticError('complete mask roster differs')
    for c in maps['centres']:
        if len(c['representative'])!=16 or any((c['representative'][j]-(c['mask']>>j))%2 for j in range(16)):raise ArithmeticError('centre parity differs')
    for i,row in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);path=ROOT/row['archive_path']
        if row['index']!=i or row['centre']!=m['centre'] or cert.hashed(path)!=row['archive_sha256'] or cert.read(path)!=archive:raise ArithmeticError('archive or chart roster differs')
        rep=m['centre']['representative']+[0]*(state.rank-16);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r=row['search']
        if r['height_bound']!=p['height'] or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('chart budget differs')
        points=backend.replay(search,m,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('orbit compression differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if row['state_key']!=state.key or row['rank_lower_bound']!=state.rank or row['admission_observations']!=state.record()['state']['observations']:raise ArithmeticError('admission replay differs')
    if data['status'] not in ('RUNNING','COMPLETE_DECLARED_POINT_ATTEMPT') or len(data['charts'])>43:raise ArithmeticError('unexpected partial/terminal state')
    if data['final_state']!=state.record() or data['rank_lower_bound']!=state.rank or (data['status']=='COMPLETE_DECLARED_POINT_ATTEMPT' and len(data['charts'])!=len(maps['rows'])):raise ArithmeticError('final point state differs')
    memory_checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED NEARCUT60 MW16',p['rows'][index]['id'],len(data['charts']),'rank >=',state.rank,flush=True)

def maps():
    p=protocol();path=D/'map-ledger.json'
    if path.exists():raise FileExistsError('preserve all-map stage')
    ledger={'status':'RUNNING','rows':[]};checkpoint(path,ledger);got={}
    def one(i):
        row=p['rows'][i];folder=D/row['id']
        supervision=supervise([SAGE,str(CAS/'prepare_nearcut60_mw16_pari_batch.sage'),'--index',str(i)],limits=Limits(p['map_wall_seconds'],p['rss_bytes']),log_path=folder/'maps.log',checkpoint_path=folder/'maps.supervisor.json',cwd=ROOT)
        ok=supervision['outcome']=='completed' and supervision['returncode']==0
        return {'id':row['id'],'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':supervision,'maps_sha256':cert.hashed(folder/'maps.json') if ok else None}
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(one,i):i for i in range(len(p['rows']))}
        for f in as_completed(pending):
            i=pending[f];got[i]=f.result();ledger['rows']=[got[j] for j in sorted(got)];checkpoint(path,ledger);print('NEARCUT60 MAP',got[i]['id'],got[i]['status'],flush=True)
    ledger['status']='PASS' if len(got)==60 and all(r['status']=='PASS' for r in got.values()) else 'FAILED_OR_CENSORED';checkpoint(path,ledger)
    if ledger['status']!='PASS':raise ArithmeticError('some maps failed; no points launched')

def baseline(index):
    p=protocol();row=p['rows'][index];folder=D/row['id'];out=folder/'baseline.json'
    if out.exists():raise FileExistsError('preserve baseline proof')
    maps=cert.read(folder/'maps.json');cache=ReductionCache(MemoryFactStore());state,bank=initial(maps,cache)
    f=next(f for f in cert.read(spec.ATLAS)['families'] if f['fibration_id']==row['family'])
    original,points=spec.specialize(f,row['parameter']);u=cert.F(maps['family_to_curve_scale_u'])
    model=tuple(map(cert.F,maps['curve']));generic=tuple(tuple(map(cert.F,P)) for P in maps['generic_points'])
    if model!=(0,0,0,original[3]/u**4,original[4]/u**6) or generic!=tuple((x/u**2,y/u**3) for x,y in points):raise ArithmeticError('baseline transport differs')
    memory_checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime)
    checkpoint(out,{'status':'PASS','protocol_hash':digest(p),'maps_sha256':cert.hashed(folder/'maps.json'),'rank':16,'state':state.record()})

def baselines():
    p=protocol();out=D/'baseline-ledger.json'
    if out.exists():raise FileExistsError('preserve baseline stage')
    if cert.read(D/'map-ledger.json')['status']!='PASS':raise ArithmeticError('all maps first')
    ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger);got={}
    def one(i):
        row=p['rows'][i];folder=D/row['id']
        s=supervise(['/usr/bin/python3',str(Path(__file__).resolve()),'baseline','--index',str(i)],limits=Limits(p['baseline_wall_seconds'],p['rss_bytes']),log_path=folder/'baseline.log',checkpoint_path=folder/'baseline.supervisor.json',cwd=ROOT)
        return {'id':row['id'],'status':'PASS' if s['outcome']=='completed' and s['returncode']==0 else 'FAILED_OR_CENSORED','supervision':s,'baseline_sha256':cert.hashed(folder/'baseline.json') if (folder/'baseline.json').exists() else None}
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(one,i):i for i in range(60)}
        for f in as_completed(pending):
            got[pending[f]]=f.result();ledger['rows']=[got[i] for i in sorted(got)];checkpoint(out,ledger)
    ledger['status']='PASS' if all(r['status']=='PASS' for r in got.values()) else 'FAILED_OR_CENSORED';checkpoint(out,ledger)
    if ledger['status']!='PASS':raise ArithmeticError('baseline gate failed; no points')

def batch():
    p=protocol();path=D/'ledger.json';mapledger=cert.read(D/'map-ledger.json')
    if path.exists():raise FileExistsError('preserve point batch')
    baseline_ledger=cert.read(D/'baseline-ledger.json')
    if baseline_ledger['status']!='PASS' or len(baseline_ledger['rows'])!=60:raise ArithmeticError('all60 baseline proofs first')
    for row,b in zip(p['rows'],baseline_ledger['rows']):
        if b['id']!=row['id'] or b['status']!='PASS' or b['baseline_sha256']!=cert.hashed(D/row['id']/'baseline.json'):raise ArithmeticError('certified baseline changed')
    if mapledger['status']!='PASS' or len(mapledger['rows'])!=60:raise ArithmeticError('all60 maps must precede points')
    for row,binding in zip(p['rows'],mapledger['rows']):
        if binding['id']!=row['id'] or binding['status']!='PASS' or binding['maps_sha256']!=cert.hashed(D/row['id']/'maps.json'):raise ArithmeticError('premade map differs')
    ledger={'status':'RUNNING','map_ledger_sha256':cert.hashed(D/'map-ledger.json'),'rows':[{**r,'status':'PENDING'} for r in p['rows']]};checkpoint(path,ledger)
    def one(index):
        row=p['rows'][index];folder=D/row['id'];r={**row,'status':'POINT_FAILED_OR_CENSORED','map_supervision':mapledger['rows'][index]['supervision']}
        search=supervise(['/usr/bin/python3',str(Path(__file__).resolve()),'worker','--index',str(index)],limits=Limits(p['worker_wall_seconds'],p['rss_bytes']),log_path=folder/'worker.log',checkpoint_path=folder/'worker.supervisor.json',cwd=ROOT);r['point_supervision']=search
        if (folder/'result.json').exists():
            d=cert.read(folder/'result.json');r.update(rank_lower_bound=d['rank_lower_bound'],charts=len(d['charts']),result_sha256=cert.hashed(folder/'result.json'))
            if search['outcome']=='completed' and search['returncode']==0:r['status']=d['status']
        return r
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(one,i):i for i in range(len(p['rows']))}
        for future in as_completed(futures):
            i=futures[future];r=future.result();ledger['rows'][i]=r;checkpoint(path,ledger);print('NEARCUT60 MW16 ATTEMPT',r['id'],r['status'],r.get('rank_lower_bound'),flush=True)
    ledger['status']='COMPLETE_FIXED_BATCH_ATTEMPTS';checkpoint(path,ledger)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['freeze','maps','batch','worker','replay','baseline','baselines']);a.add_argument('--index',type=int);v=a.parse_args();globals()[v.stage](v.index) if v.stage in ('worker','replay','baseline') else globals()[v.stage]()
