#!/usr/bin/env python3
"""Search the fixed unvisited denominator intervals on the new native11952 curve."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _is_prime
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.pointed_orbit_compression import compress
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
PARENT=ROOT/'artifacts/local/elliptic-curves/native11952-metric49-control-v1'
D=ROOT/'artifacts/local/elliptic-curves/native11952-metric49-tails-v1'
GEOMETRY_KEYS=('base_point','short_model','short_model_x_shift','pointed_chart','horizontal_matrix','ordinate_scale','coefficients')


def sources():
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (
        Path(__file__).resolve(),Path(cert.__file__).resolve(),CAS/'research_runtime/memory_store.py',
        CAS/'research_runtime/quotient_only_reduction.py',CAS/'research_runtime/cached_observation_state.py',
        CAS/'research_runtime/pointed_orbit_compression.py',CAS/'alternate_quartic_covers.py')}}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fixed-tail protocol')
    source_path=PARENT/'candidate-00/result.json';source=cert.read(source_path);terminal=cert.read(PARENT/'terminal.json')
    verification=cert.read(PARENT/'verification/result.json')
    if terminal['result_sha256']!=cert.hashed(source_path) or source['status']!='COMPLETE_DECLARED_PILOT' or len(source['charts'])!=49 or source['rank_lower_bound']!=25:raise ArithmeticError('native11952 parent differs')
    if len(verification['rows'])!=3 or any(r['status']!='PASS' for r in verification['rows']):raise ArithmeticError('parent replay incomplete')
    omitted=ROOT/'artifacts/local/elliptic-curves/r17-omitted-generic-classes-v1/verification/ledger.json'
    v=cert.read(omitted)
    if len(v['rows'])!=8 or any(r['status']!='PASS' for r in v['rows']):raise ArithmeticError('first orbit-compressed batch not replayed')
    tails=[]
    for i,row in enumerate(source['charts']):
        r=row['search']
        if r['denominator_start']!=1 or r['denominator_end']!=100000 or r['height_bound']!=100000:raise ArithmeticError('parent box changed')
        first=r['completed_denominator']+1
        if first<=100000:tails.append({'parent_chart':i,'denominator_start':first,'denominator_end':100000})
    if len(tails)!=49:raise ArithmeticError('fixed49-tail roster changed')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native11952-fixed-tails.v1','sources':sources(),
        'input_path':str(source_path.relative_to(ROOT)),'input_sha256':cert.hashed(source_path),
        'parent_verification_sha256':cert.hashed(PARENT/'verification/result.json'),'orbit_batch_verification_sha256':cert.hashed(omitted),
        'tails':tails,'height':100000,'seconds_per_chart':2,'admission_prime_bound':997,
        'coordinate_policy':{'kind':'metric','weight':'16'},'worker_wall_seconds':200,'worker_rss_bytes':1610612736,'maximum_workers':1,
        'mathematical_gate':'The generic-point-only native11952 control recovered eight directions from49 partially completed boxes. Complete exactly the missing intervals, holding rational maps and height fixed, to remove incomplete coverage as a confound when comparing the older PARI implementation. This known-curve calibration is separate from prospective selection; exceptional points remain unavailable.',
        'admission':'Preserve all raw points and exact pointed-involution skip witnesses. Independent complete-cloud audit through997 after termination.',
        'claim_boundary':'A further finite attempt, not an upper rank bound. Parent plus tail exhaust a particular recorded box only when the tail reaches its declared endpoint with the exact original geometry. No automatic completion or gain is assumed. Attempt all49 tails regardless of gains.'})
    print('FROZEN NATIVE11952 CONTROL TAILS',len(tails),flush=True)


def run():
    protocol=cert.read(D/'protocol.json');path=ROOT/protocol['input_path']
    if protocol['sources']!=sources() or cert.hashed(path)!=protocol['input_sha256']:raise ArithmeticError('fixed-tail input changed')
    source=cert.read(path);output=D/'candidate-00/result.json';folder=output.parent
    if output.exists():raise FileExistsError('preserve tail attempt')
    cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(source['arithmetic_facts'])
    state=MWState.from_record(source['final_state'],cache=cache);model=tuple(map(cert.F,source['curve']))
    if state.rank!=25:raise ArithmeticError('initial rank differs')
    data={'schema':'elliptic-curves.native11952-fixed-tail-result.v1','protocol_hash':digest(protocol),
        'family':source['family'],'parameter':source['parameter'],'curve':source['curve'],'generic_points':source['generic_points'],
        'initial_dimension':state.rank,'initial_state':state.record(),
        'centres':source['centres'],'charts':[],'status':'RUNNING','rank_lower_bound':state.rank,'final_state':state.record(),'arithmetic_facts':cache.store.snapshot()}
    checkpoint(output,data);primes=tuple(p for p in range(3,998) if _is_prime(p))
    for tail in protocol['tails']:
        i=tail['parent_chart'];parent=source['charts'][i];centre=parent['centre']
        rep=centre['representative']+[0]*(state.rank-source['initial_dimension'])
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        chart=search.chart_record()
        if any(chart[k]!=parent['search'][k] for k in GEOMETRY_KEYS):raise ArithmeticError('tail changed its original chart geometry')
        outcome=search.search(100000,2,denominator_start=tail['denominator_start'],denominator_end=tail['denominator_end'],checkpoint_dir=folder/'charts'/state.key)
        compression=compress(model,state.basis,rep,outcome.curve_points)
        for j in compression['kept_indices']:
            state=state.adjoin(outcome.curve_points[j],cache=cache,extra_primes=primes)
            if not isinstance(state,MWState):state=MWState.from_record(state.record(),cache=cache)
        data['charts'].append({'parent_chart':i,'centre':centre,'search':outcome.record,'admission_compression':compression,'rank_lower_bound':state.rank,'admission_prime_bound':997})
        data.update(rank_lower_bound=state.rank,final_state=state.record(),arithmetic_facts=cache.store.snapshot());checkpoint(output,data)
        print('NATIVE11952 CONTROL TAIL',i+1,'through',outcome.record['completed_denominator'],'rank',state.rank,flush=True)
    data['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,data)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);a=p.parse_args();prepare() if a.stage=='prepare' else run()
