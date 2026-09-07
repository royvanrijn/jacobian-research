#!/usr/bin/env python3
"""Select independent columns from the retained control cloud; no search rerun."""
import argparse
from pathlib import Path
from dataclasses import asdict
import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
from half_lattice_pointed_sieve import linear_combination_python
from pointed_quartic_search import PointedQuarticSearch,point_record
from memory_rank_certificate import checked_rank
from audit_recorded_point_mod2_rank_v3 import signature,insert,_primes_up_to
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/factor-free-known28-control-v1'
OLD=ROOT/'artifacts/local/elliptic-curves/inventory188-own27-geometry-control-v1'
OUT=ART/'factor_free_known28_control_v1.json'

def expected():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen control inputs differ')
    seed=cert.read(OLD/'seed.json');m=cert.read(D/'map.json');r=cert.read(D/'search.json')
    assert m['centre']==cert.read(OLD/'maps.json')['centres'][4]
    model=tuple(map(cert.F,seed['curve']));points=[tuple(map(cert.F,P)) for P in seed['points']]
    C=linear_combination_python(model,points,m['centre']['representative'])
    search=PointedQuarticSearch(curve=model,subgroup=[],centre={'point':point_record(C)},coordinate_policy=m['coordinate_policy'])
    assert r['height_bound']==p['height'] and r['timeout_seconds']==p['seconds'] and r['gp_binary_sha256']==p['gp_sha256']
    found=backend.replay(search,m,r);seen={(x,abs(y)) for x,y in points}
    for P in found:
        if (P[0],abs(P[1])) not in seen:seen.add((P[0],abs(P[1])));points.append(P)
    cache=ReductionCache(MemoryFactStore());pivots={};sigs=[]
    for prime in _primes_up_to(997):
        if prime==2:continue
        try:s=signature(cache,model,points,prime)
        except ValueError:continue
        before=len(pivots)
        for row in s.rows:insert(pivots,row)
        if len(pivots)>before:sigs.append(asdict(s))
    chosen=sorted(pivots);basis=[points[i] for i in chosen]
    rank=checked_rank(model,basis,[s['prime'] for s in sigs],seed['rank_certificate']['no_rational_2_torsion_prime'])
    paths=[Path(__file__).resolve(),D/'protocol.json',D/'map.json',D/'search.json',D/'worker.log',D/'worker.supervisor.json']
    return {'status':'PASS' if r['status']=='bounded_search_complete' and len(chosen)>=28 else 'CONTROL_NOT_RECOVERED',
        'sources':{**p['sources'],**{str(x.relative_to(ROOT)):cert.hashed(x) for x in paths}},
        'curve':seed['curve'],'points':[list(map(str,P)) for P in points],'signatures':sigs,
        'independent_column_indices':chosen,'independent_points':[list(map(str,P)) for P in basis],
        'rank_certificate':rank,'rank_lower_bound':len(chosen),'completed_boxes':int(r['status']=='bounded_search_complete'),
        'point_seconds':r['wall_seconds'],'search_cpu_ms':r['search_cpu_ms'],
        'previous_failure':'The control search completed; its wrapper wrongly required every retained cloud point to be independent. This separate audit selects independent columns and replays the retained transcript without repeating point search.',
        'boundary':p['scope']}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=expected()
    if a.check:assert digest(r)==digest(cert.read(OUT))
    else:
        if OUT.exists():raise FileExistsError('preserve control audit')
        checkpoint(OUT,r)
    print(r['status'],r['rank_lower_bound'],len(r['points']))
