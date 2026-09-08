#!/usr/bin/env python3
"""Exact finite neighbour audit around the98 retained control-only CVP proposals."""
import argparse
from pathlib import Path
from collections import Counter
from functools import lru_cache
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from search_nagao_u42_skew_height import short_add,short_multiply
from search_observability import point_visibility
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'native29-cvp-neighbours-v1';INPUT=ART/'native11952_rank28_coset_visibility_v1.json';CONTROL=LOCAL/'native11952-pari49-control-v1/candidate-00/result.json';HIGH=LOCAL/'native11952-height125-control-v1/125000/result.json';OUT=ART/'native29_cvp_neighbours_v1.json'
def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),INPUT,CONTROL,HIGH,ART/'native11952_rank28_coset_visibility_replay_v1.json',ROOT/'elliptic-curves/cas/search_nagao_u42_skew_height.py',ROOT/'elliptic-curves/cas/search_observability.py',ROOT/'elliptic-curves/cas/memory_rank_certificate.py')}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve neighbour protocol')
    d=cert.read(INPUT);gate=cert.read(ART/'native11952_rank28_coset_visibility_replay_v1.json')
    if gate['status']!='PASS' or len(d['proposals'])!=98 or len(d['basis'])!=28 or len(d['lll_matrix'])!=28:raise ArithmeticError('exact98-proposal baseline required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.native29-cvp-neighbours.protocol.v1','sources':sources(),'proposals':5488,'baseline_proposals':98,'neighbours_per_baseline':56,'wall_seconds':300,'rss_bytes':1610612736,'gate':'The known remaining29th direction has one exact translated representative at918522, while the original49 charts completed at125000 and certify28. A single canonical-height CVP proposal does not minimize rational chart-coordinate height. Check every plus/minus row of the previously recorded unimodular28-dimensional LLL matrix around all98 retained proposals.','boundaries':'Exactly5488 group translates, no point enumeration, no new numerical metric or adaptive neighbourhood enlargement. Exact words and visibility replay. Retrospective known-control diagnostics only; no oracle enters prospective centres or selection. A smaller observed coordinate is not a global minimum or a new curve.'})
def context():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen neighbour sources differ')
    d=cert.read(INPUT);control=cert.read(CONTROL);high=cert.read(HIGH);model=tuple(map(cert.F,d['curve']));basis=[tuple(map(cert.F,P)) for P in d['basis']];proof=d['rank29_witness']['rank_certificate'];oracle=[tuple(map(cert.F,P)) for P in d['oracle_points']]
    checked_rank(model,basis+oracle,[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime'])
    if len(high['charts'])!=49 or any(r['search']['status']!='bounded_search_complete' or r['search']['height_bound']!=125000 for r in high['charts']):raise ArithmeticError('completed125000 coverage missing')
    high_records=[]
    for i,row in enumerate(high['charts']):
        if row['search']['base_point']!=control['charts'][i]['search']['base_point'] or row['search']['coefficients']!=control['charts'][i]['search']['coefficients'] or row['search']['horizontal_matrix']!=control['charts'][i]['search']['horizontal_matrix']:raise ArithmeticError('old/new fixed chart mismatch')
        # A completed whole-box PARI call covers every denominator through H.
        # This adapter is never applied to a timeout or failure.
        high_records.append({**row['search'],'completed_denominator':125000})
    return p,d,control,high_records,model,basis,oracle
def run():
    if OUT.exists():raise FileExistsError('preserve neighbour output')
    p,d,control,high,model,basis,oracle=context();directions=[]
    for word in d['lll_matrix']:
        P=None
        for q,n in zip(basis,word):
            if n:P=short_add(model[3],P,short_multiply(model[3],q,n))
        if P is None:raise ArithmeticError('zero LLL direction')
        directions.append(P)
    result={'schema':'elliptic-curves.native29-cvp-neighbours.v1','status':'RUNNING','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':[]};checkpoint(OUT,result);counts=Counter()
    for baseline_index,r in enumerate(d['proposals']):
        base=tuple(map(cert.F,r['point']))
        for j,Q in enumerate(directions):
            for sign in (1,-1):
                point=short_add(model[3],base,(Q[0],sign*Q[1]));word=[x+sign*y for x,y in zip(r['word'],d['lll_matrix'][j])]
                v=point_visibility(high[r['chart_index']],point)
                if v['status'] in ('VISIBLE_NOT_RECORDED','UNSEARCHED_INTERVAL'):raise ArithmeticError('completed125000 control discrepancy')
                result['rows'].append({'baseline_index':baseline_index,'chart_index':r['chart_index'],'lll_row':j,'lll_sign':sign,'oracle_sign':r['sign'],'word':word,'point':list(map(str,point)),'visibility':v});counts[v['status']]+=1
        checkpoint(OUT,result)
        if (baseline_index+1)%14==0:print('NATIVE29 NEIGHBOURS',len(result['rows']),dict(counts),flush=True)
    best=min(result['rows'],key=lambda r:r['visibility'].get('minimum_affine_height') or 10**1000)
    result.update(status='COMPLETE_DECLARED_NEIGHBOURS',best=best,status_counts=dict(counts),claim_boundary=p['boundaries']);checkpoint(OUT,result);print('BEST NEIGHBOUR',best['chart_index'],best['visibility']['minimum_affine_height'],flush=True)
def replay(output):
    if output.exists():raise FileExistsError('preserve neighbour replay')
    p,d,control,high,model,basis,oracle=context();data=cert.read(OUT);roster=[];counts=Counter()
    if data['sources']!=sources() or data['protocol_sha256']!=cert.hashed(D/'protocol.json'):raise ArithmeticError('neighbour result binding differs')
    @lru_cache(None)
    def multiple(j,n):return short_multiply(model[3],basis[j],n)
    @lru_cache(None)
    def translated(sign,word):
        P=oracle[0][0],sign*oracle[0][1]
        for j,n in enumerate(word):
            if n:P=short_add(model[3],P,multiple(j,n))
        return P
    for r in data['rows']:
        i,j,sign=r['baseline_index'],r['lll_row'],r['lll_sign'];roster.append((i,j,sign));baseline=d['proposals'][i];word=[x+sign*y for x,y in zip(baseline['word'],d['lll_matrix'][j])]
        if r['word']!=word or r['oracle_sign']!=baseline['sign'] or r['chart_index']!=baseline['chart_index']:raise ArithmeticError('neighbour word/roster differs')
        P=translated(r['oracle_sign'],tuple(word))
        if list(map(str,P))!=r['point'] or point_visibility(high[r['chart_index']],P)!=r['visibility']:raise ArithmeticError('independent rational replay differs')
        counts[r['visibility']['status']]+=1
    best=min(data['rows'],key=lambda r:r['visibility'].get('minimum_affine_height') or 10**1000)
    if data['status']!='COMPLETE_DECLARED_NEIGHBOURS' or roster!=[(i,j,s) for i in range(98) for j in range(28) for s in (1,-1)] or best!=data['best'] or dict(counts)!=data['status_counts']:raise ArithmeticError('terminal5488 neighbour result differs')
    checkpoint(output,{'schema':'elliptic-curves.native29-cvp-neighbours-replay.v1','status':'PASS','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in (Path(__file__).resolve(),OUT,D/'protocol.json')},'exact_words':5488,'status_counts':dict(counts),'best':best,'claim_boundary':p['boundaries']});print('REPLAYED5488 EXACT CONTROL TRANSLATES',best['visibility']['minimum_affine_height'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);p.add_argument('--output',type=Path,default=ART/'native29_cvp_neighbours_replay_v1.json');a=p.parse_args();replay(a.output.resolve()) if a.stage=='replay' else globals()[a.stage]()
