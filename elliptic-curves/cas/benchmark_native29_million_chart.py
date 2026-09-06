#!/usr/bin/env python3
"""Two frozen-engine cost attempts on the same retrospective control chart."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import pari_pointed_backend as pari_backend
from pointed_quartic_search import PointedQuarticSearch
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'native29-million-chart-benchmark-v1';MAPS=LOCAL/'native11952-pari-maps-v1/maps.json';SEED=LOCAL/'native11952-height125-control-v1/125000/result.json';ORACLE=ART/'native11952_rank28_coset_visibility_replay_v1.json';OUT=ART/'native29_million_chart_benchmark_v1.json'
def sources():
    return {**pari_backend.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),MAPS,SEED,ROOT/'elliptic-curves/cas/memory_rank_certificate.py',ROOT/'elliptic-curves/cas/research_runtime/cached_observation_state.py',ROOT/'elliptic-curves/cas/research_runtime/quotient_only_reduction.py')}}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve million-chart benchmark')
    gate=cert.read(ORACLE);seed=cert.read(SEED)
    if gate['status']!='PASS' or gate['best_per_direction'][0]['chart_index']!=12 or gate['best_per_direction'][0]['visibility']['minimum_affine_height']!=918522 or seed['rank_lower_bound']!=28:raise ArithmeticError('exact retrospective control gate differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native29-million-chart-benchmark.protocol.v1','sources':sources(),'oracle_gate_sha256':cert.hashed(ORACLE),'chart_index':12,'height':1000000,'seconds_per_chart':60,'arms':['pari','gmp'],'rss_bytes':1610612736,'outer_seconds_per_arm':90,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'gate':'The exact retrospective audit exhibits a representative of the last known29th direction at height918522 in original generic chart12. Existing10-second million-height PARI calls were censored. Measure two unchanged engines on this same chart and full million-height box, with a fixed60seconds each; retain every partial transcript.','boundaries':'The chart choice and budget are informed by a published-point diagnostic, so recovery is retrospective calibration, not blind validation or a new curve. Workers consume only the existing28-point seed and old generic map, never an oracle point or word. No additional charts, retries, point absence, exact rank or automatic prospective escalation.'})
def initial():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('benchmark sources differ')
    seed=cert.read(SEED);cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(seed['arithmetic_facts']);state=MWState.from_record(seed['final_state'],cache=cache)
    if state.rank!=28:raise ArithmeticError('control28 basis differs')
    mapping=cert.read(MAPS)['rows'][p['chart_index']];rep=mapping['centre']['representative']+[0]*11;search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=mapping['coordinate_policy'])
    return p,seed,cache,state,mapping,search
def run(arm):
    p,seed,cache,state,mapping,search=initial();out=D/(arm+'.json')
    if arm not in p['arms'] or out.exists():raise FileExistsError('preserve declared benchmark arm')
    if arm=='pari':record,points=pari_backend.execute(search,mapping,p['height'],p['seconds_per_chart'],p['gp_sha256'])
    else:
        found=search.search(p['height'],p['seconds_per_chart'],checkpoint_dir=D/'gmp-checkpoints');record,points=found.record,found.curve_points
    checkpoint(out,{'status':'RAW_ATTEMPT_RETAINED','protocol_hash':digest(p),'arm':arm,'initial_state':state.record(),'search':record})
    for P in points:state=state.adjoin(P,cache=cache)
    model=tuple(map(cert.F,seed['curve']));proof=checked_rank(model,state.basis,state.reductions.primes,state.no_two_torsion_prime)
    checkpoint(out,{'status':'TERMINAL_BENCHMARK_ARM','protocol_hash':digest(p),'arm':arm,'initial_state':search.state.record(),'search':record,'rank_lower_bound':state.rank,'independent_points':[list(map(str,P)) for P in state.basis],'rank_certificate':proof,'final_state':state.record()});print('MILLION CONTROL',arm,record['status'],record.get('completed_denominator'),record['wall_seconds'],'rank',state.rank,flush=True)
def replay():
    if OUT.exists():raise FileExistsError('preserve benchmark replay')
    rows=[]
    for arm in ('pari','gmp'):
        p,seed,cache,state,mapping,search=initial();r=cert.read(D/(arm+'.json'));s=r['search']
        if r['protocol_hash']!=digest(p) or r['initial_state']!=state.record() or r['status']!='TERMINAL_BENCHMARK_ARM' or s['height_bound']!=p['height'] or s['timeout_seconds']!=p['seconds_per_chart']:raise ArithmeticError('fixed benchmark record differs')
        points=pari_backend.replay(search,mapping,s) if arm=='pari' else search.verify_record(s).curve_points
        for P in points:state=state.adjoin(P,cache=cache)
        proof=checked_rank(tuple(map(cert.F,seed['curve'])),state.basis,state.reductions.primes,state.no_two_torsion_prime)
        if state.record()!=r['final_state'] or r['rank_lower_bound']!=state.rank or [list(map(str,P)) for P in state.basis]!=r['independent_points'] or digest(proof)!=digest(r['rank_certificate']):raise ArithmeticError('control point proof differs')
        rows.append({'arm':arm,'status':s['status'],'wall_seconds':s['wall_seconds'],'completed_denominator':s.get('completed_denominator'),'points_returned':len(points),'rank_lower_bound':state.rank,'input_sha256':cert.hashed(D/(arm+'.json'))})
    checkpoint(OUT,{'schema':'elliptic-curves.native29-million-chart-benchmark.v1','status':'PASS','sources':sources(),'protocol':p,'arms':rows,'claim_boundary':p['boundaries']});print('REPLAYED TWO MILLION-HEIGHT CONTROL ARMS',rows,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);p.add_argument('--arm',choices=['pari','gmp']);a=p.parse_args();run(a.arm) if a.stage=='run' else globals()[a.stage]()
