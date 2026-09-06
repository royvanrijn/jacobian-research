#!/usr/bin/env python3
"""Independent rational group replay of98 retrospective oracle translations."""
from pathlib import Path
import argparse
from collections import Counter
from functools import lru_cache
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from search_nagao_u42_skew_height import short_add,short_multiply
from search_observability import point_visibility
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';INPUT=ART/'native11952_rank28_coset_visibility_v1.json';OUT=ART/'native11952_rank28_coset_visibility_replay_v1.json';CONTROL=ROOT/'artifacts/local/elliptic-curves/native11952-pari49-control-v1/candidate-00/result.json'
def main(output):
    if output.exists():raise FileExistsError('preserve oracle replay')
    data=cert.read(INPUT);source=cert.read(CONTROL)
    for p,h in data['sources'].items():
        if cert.hashed(ROOT/p)!=h:raise ArithmeticError('oracle source changed')
    model=tuple(map(cert.F,data['curve']));basis=[tuple(map(cert.F,p)) for p in data['basis']];oracles=[tuple(map(cert.F,p)) for p in data['oracle_points']];proof=data['rank29_witness']['rank_certificate']
    if len(basis)!=28 or len(oracles)!=1 or source['curve']!=data['curve']:raise ArithmeticError('oracle/baseline dimensions differ')
    checked_rank(model,basis+oracles,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    @lru_cache(None)
    def multiple(j,k):return short_multiply(model[3],basis[j],k)
    counts=Counter();roster=[]
    for r in data['proposals']:
        i,k,sign=r['chart_index'],r['oracle_index'],r['sign'];roster.append((i,k,sign))
        if not 0<=i<49 or not 0<=k<1 or sign not in (1,-1) or len(r['word'])!=28 or any(type(x)!=int for x in r['word']):raise ArithmeticError('invalid fixed group word')
        point=oracles[k][0],sign*oracles[k][1]
        for j,n in enumerate(r['word']):
            if n:point=short_add(model[3],point,multiple(j,n))
        if point!=tuple(map(cert.F,r['point'])) or not cert.is_on_weierstrass_curve(model,point):raise ArithmeticError('exact translated point differs')
        v=point_visibility(source['charts'][i]['search'],point)
        if v!=r['visibility'] or r['published_index']!=[27][k]:raise ArithmeticError('exact visibility differs')
        counts[v['status']]+=1
    if roster!=[(i,k,sign) for i in range(49) for k in range(1) for sign in (1,-1)] or dict(counts)!=data['status_counts']:raise ArithmeticError('fixed98 roster differs')
    best=[min((r for r in data['proposals'] if r['oracle_index']==k),key=lambda r:r['visibility'].get('minimum_affine_height') or 10**1000) for k in range(1)]
    if data['status']!='COMPLETE_DECLARED_ORACLE_AUDIT' or best!=data['best_per_direction']:raise ArithmeticError('terminal best proposals differ')
    checkpoint(output,{'schema':'elliptic-curves.native11952-rank28-coset-visibility-replay.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),INPUT,CONTROL,ROOT/'elliptic-curves/cas/search_nagao_u42_skew_height.py',ROOT/'elliptic-curves/cas/search_observability.py',ROOT/'elliptic-curves/cas/memory_rank_certificate.py')},'exact_group_words_checked':98,'status_counts':dict(counts),'best_per_direction':best,'claim_boundary':'Exact rational group words and chart-square/visibility identities. Numerical CVP optimality is not asserted. No prospective oracle leakage, new curve or upper bound.'})
    print('REPLAYED98 NATIVE29 TRANSLATIONS',[(r['published_index'],r['visibility']['minimum_affine_height']) for r in best],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUT);a=p.parse_args();main(a.output.resolve())
