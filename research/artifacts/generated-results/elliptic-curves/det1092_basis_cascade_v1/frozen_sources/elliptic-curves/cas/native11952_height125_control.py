#!/usr/bin/env python3
"""A bounded generic17-only125000 control following the exact113933 visibility witness."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
from memory_rank_certificate import checked_rank
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.preloaded_prime_state import preload
from research_runtime.rotated_observation_state import rotate
from research_runtime.pointed_orbit_compression import compress
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'native11952-height125-control-v1';MAPS=LOCAL/'native11952-pari-maps-v1/maps.json';HEIGHTS=(125000,)
def sources():
    names=['native11952_height125_control.py','memory_rank_certificate.py','research_runtime/search_state.py','research_runtime/cached_observation_state.py','research_runtime/preloaded_prime_state.py','research_runtime/rotated_observation_state.py','research_runtime/pointed_orbit_compression.py','research_runtime/quotient_only_reduction.py','research_runtime/supervisor.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve height pair')
    gate=ART/'native11952_translated_visibility_replay_v1.json';g=cert.read(gate);maps=cert.read(MAPS)
    if g['status']!='PASS' or [r['visibility']['minimum_affine_height'] for r in g['best_per_direction']]!=[113933,918522] or maps['status']!='COMPLETE_DECLARED_MAPS' or len(maps['rows'])!=49 or len(maps['generic_points'])!=17:raise ArithmeticError('exact retrospective height gate failed')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native11952-height125-control.v1','sources':sources(),'maps_sha256':cert.hashed(MAPS),'height_bounds':list(HEIGHTS),'seconds_per_chart':10,'worker_wall_seconds':300,'replay_wall_seconds':300,'rss_bytes':1610612736,'maximum_workers':1,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'gate_source_sha256':cert.hashed(gate),'gate':'A separate196-word retrospective exact audit locates representatives of both missing known-control directions at113933 and918522. The first direction has a witness just above100000. Test height125000, the next25000 increment beyond113933, in a generic17-only control. The larger million-height arm has timed out on its first ten charts at10seconds each; preserve its fixed full attempt. This separate bounded control changes no existing protocol. This is calibration of known points, never new-curve discovery.','policy':'One arm starts from only the17 generic points in the frozen49-map file, with the completed100000 arm retained as its baseline. At most three concurrent subprocess campaigns overall: this control, the already running million-height control, and the saved-shard trace extension. Use identical map order, backend, admission primes and10second per-chart limits. Attempt all49 charts regardless of gain; no oracle coordinates, words, published point list or rank29 stopping rule enter either worker.','boundaries':'No exact rank, global upper bound or general prospective recovery guarantee. Preserve all partial outputs and failures. A higher-height prospective policy requires completed exact control replays first.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['maps_sha256']!=cert.hashed(MAPS) or p['height_bounds']!=list(HEIGHTS):raise ArithmeticError('height-pair bindings differ')
    return p
def initial(cache,maps):
    model=tuple(map(cert.F,maps['curve']));points=tuple(tuple(map(cert.F,P)) for P in maps['generic_points']);state=raw_state(model,points,cache=cache,prime_bound=1000);state=MWState.from_record(state.record(),cache=cache);state,info=preload(state,cache,997)
    if state.rank!=17:raise ArithmeticError('generic17 input required')
    return model,state,info
def worker(height):
    p=protocol()
    if height not in HEIGHTS:raise ValueError('undeclared arm')
    maps=cert.read(MAPS);folder=D/str(height);out=folder/'result.json'
    if out.exists():raise FileExistsError('preserve height attempt')
    cache=ReductionCache(MemoryFactStore());model,state,info=initial(cache,maps)
    data={'schema':'elliptic-curves.native11952-height125-result.v1','family':'11952-control','parameter':'fixed-native-rank29-control','curve':maps['curve'],'generic_points':maps['generic_points'],'height':height,'protocol_hash':digest(p),'maps_sha256':cert.hashed(MAPS),'initial_state':state.record(),'initial_dimension':17,'centres':[m['centre'] for m in maps['rows']],'prime_preload':info,'charts':[],'status':'RUNNING','rank_lower_bound':17,'final_state':state.record(),'arithmetic_facts':cache.store.snapshot()};checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        state,archive=rotate(state);ap=folder/'states'/f'{i:03}.json';checkpoint(ap,archive);rep=m['centre']['representative']+[0]*(state.rank-17);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r,points=backend.execute(search,m,height,p['seconds_per_chart'],p['gp_sha256']);compression=compress(model,state.basis,rep,points)
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        data['charts'].append({'index':i,'centre':m['centre'],'archive_path':str(ap.relative_to(ROOT)),'archive_sha256':cert.hashed(ap),'search':r,'admission_compression':compression,'admission_observations':state.record()['state']['observations'],'state_key':state.key,'rank_lower_bound':state.rank});data.update(final_state=state.record(),rank_lower_bound=state.rank,arithmetic_facts=cache.store.snapshot());checkpoint(out,data);print('HEIGHT PAIR',height,i+1,r['status'],'rank',state.rank,flush=True)
    data['status']='COMPLETE_DECLARED_HEIGHT_ARM';checkpoint(out,data)
def replay(height):
    p=protocol();maps=cert.read(MAPS);folder=D/str(height);data=cert.read(folder/'result.json');cache=ReductionCache(MemoryFactStore());model,state,info=initial(cache,maps)
    if height not in HEIGHTS or data['height']!=height or data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(MAPS) or data['initial_state']!=state.record() or data['prime_preload']!=info or data['generic_points']!=maps['generic_points'] or data['curve']!=maps['curve'] or data['centres']!=[m['centre'] for m in maps['rows']]:raise ArithmeticError('initial arm differs')
    for i,row in enumerate(data['charts']):
        m=maps['rows'][i];state,archive=rotate(state);ap=ROOT/row['archive_path']
        if row['index']!=i or row['centre']!=m['centre'] or cert.hashed(ap)!=row['archive_sha256'] or cert.read(ap)!=archive:raise ArithmeticError('archived state differs')
        rep=m['centre']['representative']+[0]*(state.rank-17);search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=m['coordinate_policy']);r=row['search']
        if r['height_bound']!=height or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('arm budget differs')
        points=backend.replay(search,m,r);compression=compress(model,state.basis,rep,points)
        if compression!=row['admission_compression']:raise ArithmeticError('compression differs')
        for j in compression['kept_indices']:state=state.adjoin(points[j],cache=cache)
        if state.key!=row['state_key'] or state.rank!=row['rank_lower_bound'] or state.record()['state']['observations']!=row['admission_observations']:raise ArithmeticError('admissions differ')
    if state.record()!=data['final_state'] or state.rank!=data['rank_lower_bound'] or (data['status']=='COMPLETE_DECLARED_HEIGHT_ARM' and len(data['charts'])!=49):raise ArithmeticError('final state differs')
    checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime);print('REPLAYED HEIGHT PAIR',height,len(data['charts']),state.rank,flush=True)
def launch():
    p=protocol()
    for height in HEIGHTS:
        folder=D/str(height)
        for label,seconds in [('worker',p['worker_wall_seconds']),('replay',p['replay_wall_seconds'])]:
            target=folder/(label+'.supervisor.json')
            if target.exists():raise FileExistsError('preserve control supervision')
            s=run([sys.executable,str(Path(__file__).resolve()),label,'--height',str(height)],limits=Limits(seconds,p['rss_bytes']),log_path=folder/(label+'.log'),checkpoint_path=target,cwd=ROOT);print('HEIGHT STAGE',height,label,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
            if s['outcome']!='completed' or s['returncode']!=0:return
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch','worker','replay']);p.add_argument('--height',type=int);a=p.parse_args();globals()[a.stage](a.height) if a.stage in ('worker','replay') else globals()[a.stage]()
