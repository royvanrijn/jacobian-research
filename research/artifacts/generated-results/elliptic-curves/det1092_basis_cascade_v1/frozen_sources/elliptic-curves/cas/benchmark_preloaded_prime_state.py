#!/usr/bin/env python3
"""Exact fixed-point admission comparison; no point search or new rank claim."""
from pathlib import Path
from time import monotonic
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import _is_prime
from research_runtime.memory_store import MemoryFactStore
from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as ReductionCache
from research_runtime.cached_observation_state import CachedObservationMWState as MWState
from research_runtime.preloaded_prime_state import preload
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/preloaded-prime-state-benchmark-v1'

def main():
    protocol=cert.read(D/'protocol.json');path=D/'input.json'
    if cert.hashed(path)!=protocol['input_sha256']:raise ArithmeticError('snapshot changed')
    data=cert.read(path);r=data['charts'][0];points=[(cert.F(p['x']),cert.F(p['y'])) for p in r['search']['finite_curve_points']];points=[points[i] for i in r['admission_compression']['kept_indices']];primes=tuple(p for p in range(3,998) if _is_prime(p));outputs=[]
    for mode in ('existing','preloaded'):
        cache=ReductionCache(MemoryFactStore());cache.store.import_snapshot(data['arithmetic_facts']);state=MWState.from_record(data['initial_state'],cache=cache);n=len(state.observations);t=monotonic();bank=None
        if mode=='preloaded':state,bank=preload(state,cache)
        prep=monotonic()-t;t=monotonic();history=[]
        for p in points:
            state=state.adjoin(p,cache=cache,extra_primes=primes if mode=='existing' else ())
            history.append({'rank':state.rank,'basis':state.basis,'observation':{'point':state.observations[-1].point,'status':state.observations[-1].status,'relation':state.observations[-1].finite_relation_mask}})
        elapsed=monotonic()-t
        result={'mode':mode,'initial_observations':n,'point_admissions':len(points),'preparation_seconds':prep,'admission_seconds':elapsed,'bank':bank,'history':history,'final_state':state.record()};outputs.append(result);checkpoint(D/(mode+'.json'),result);print('ADMISSION BENCHMARK',mode,len(points),prep,elapsed,flush=True)
    if outputs[0]['history']!=outputs[1]['history']:raise ArithmeticError('finite admission outcomes differ')
    state=outputs[1]['final_state']['state'];model=tuple(map(cert.F,data['curve']));basis=[tuple(map(cert.F,p)) for p in state['reductions']['points']];rank=cert.checked_rank(model,basis,state['reductions']['primes'],state['no_two_torsion_prime'])
    checkpoint(D/'result.json',{'status':'PASS','point_admissions':len(points),'identical_basis_rank_observation_relation_history':True,'rank_certificate':rank,'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ROOT/'elliptic-curves/cas/research_runtime/preloaded_prime_state.py')},'rows':[{'mode':r['mode'],'preparation_seconds':r['preparation_seconds'],'admission_seconds':r['admission_seconds']} for r in outputs],'claim_boundary':'One fixed retained chart. Same ordered points, exact finite-admission outcomes and certified subgroup; full prime banks change state keys and finite certificate presentation. No search or new independent point.'})
if __name__=='__main__':main()
