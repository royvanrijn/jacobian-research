#!/usr/bin/env python3
"""Exact finite mod2 lower bounds on an explicitly hash-pinned retained snapshot.

A retained snapshot can come from a censored search; its points still prove lower
bounds. This audit does not claim the source search completed its chart plan.
The v1 checker and its already emitted certificates remain untouched.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
from mod2_reduction_independence import Mod2ReductionSignature,_primes_up_to
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2]

def sources():
    paths=[Path(__file__).resolve(),Path(cert.__file__).resolve()]
    paths += [ROOT/'elliptic-curves/cas'/p for p in ('mod2_reduction_independence.py','elliptic_candidate_record.py',
              'research_runtime/finite_reduction.py','research_runtime/memory_store.py','research_runtime/arithmetic.py','research_runtime/store.py')]
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def signature(cache,model,points,prime):
    record,table,_=cache.quotient(model,prime);masks=[]
    for point in points:
        reduced=None if any(x.denominator%prime==0 for x in point) else tuple(x.numerator*pow(x.denominator,-1,prime)%prime for x in point)
        masks.append(table[reduced])
    return Mod2ReductionSignature(prime,record['group_order'],record['doubled_order'],record['dimension'],
        tuple(tuple((mask>>j)&1 for mask in masks) for j in range(record['dimension'])))

def insert(pivots,row):
    value=sum(int(bit)<<i for i,bit in enumerate(row))
    while value:
        j=(value&-value).bit_length()-1
        if j not in pivots:pivots[j]=value;return True
        value^=pivots[j]
    return False

def build(path,output,bound,expected_hash):
    if output.exists():raise FileExistsError('preserve recorded-point audit')
    if cert.hashed(path)!=expected_hash:raise ArithmeticError('explicit input snapshot hash changed')
    data=cert.read(path)
    if not data.get('charts') or not data.get('final_state'):raise ArithmeticError('retained point transcript required')
    model=tuple(map(cert.F,data['curve']))
    if any(model[:3]):raise ArithmeticError('short model required for direct projective reductions')
    basis=[tuple(map(cert.F,p)) for p in data['final_state']['state']['reductions']['points']];points=list(basis);seen={(x,abs(y)) for x,y in points}
    for chart in data['charts']:
        for raw in chart['search']['finite_curve_points']:
            p=cert.F(raw['x']),cert.F(raw['y']);key=p[0],abs(p[1])
            if key not in seen:seen.add(key);points.append(p)
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('retained point membership failed')
    result={'schema':'elliptic-curves.recorded-point-mod2-rank.v2','sources':sources(),'input_path':str(path.resolve().relative_to(ROOT)),
        'input_sha256':expected_hash,'input_status':data['status'],'family':data.get('family','published-R17'),'parameter':data['parameter'],'curve':data['curve'],
        'points':[list(map(str,p)) for p in points],'original_rank_lower_bound':len(basis),'prime_bound':bound,
        'rank_lower_bound':0,'signatures':[],'prime_progress':[],'status':'RUNNING',
        'claim_boundary':'Exact finite quotient lower bounds on the retained point cloud only. No new point search, public oracle point, exact rank or rank upper bound.'}
    checkpoint(output,result);cache=ReductionCache(MemoryFactStore());pivots={};dimension=0
    for prime in _primes_up_to(bound):
        if prime==2:continue
        try:sig=signature(cache,model,points,prime)
        except ValueError:continue
        before=len(pivots);dimension+=sig.quotient_dimension
        for row in sig.rows:insert(pivots,row)
        result['prime_progress'].append({'prime':prime,'available_quotient_dimension':dimension,'finite_rank':len(pivots)})
        if len(pivots)>before:
            result['signatures'].append(asdict(sig));result['rank_lower_bound']=len(pivots);checkpoint(output,result)
            print('RECORDED MOD2',result['family'],result['parameter'],'prime',prime,'rank',len(pivots),'available rows',dimension,flush=True)
        if len(pivots)>=32:break
    indices=sorted(pivots);selected=[points[i] for i in indices];primes=[r['prime'] for r in result['signatures']]
    proof=cert.checked_rank(model,selected,primes,data['final_state']['state']['no_two_torsion_prime'])
    result.update(status='COMPLETE_DECLARED_FINITE_AUDIT',independent_column_indices=indices,
        independent_points=[list(map(str,p)) for p in selected],rank_certificate=proof,rank_lower_bound=len(indices))
    if cert.hashed(path)!=expected_hash:raise ArithmeticError('input changed during finite audit')
    checkpoint(output,result);print('RECERTIFIED RECORDED POINTS',result['family'],result['parameter'],'rank >=',len(indices),'previous',len(basis),flush=True)

def check(path):
    data=cert.read(path)
    if data['sources']!=sources():raise ArithmeticError('recorded-point checker sources changed')
    model=tuple(map(cert.F,data['curve']));points=[tuple(map(cert.F,p)) for p in data['points']]
    if any(not cert.is_on_weierstrass_curve(model,p) for p in points):raise ArithmeticError('point membership failed')
    cache=ReductionCache(MemoryFactStore());pivots={}
    for old in data['signatures']:
        sig=signature(cache,model,points,old['prime'])
        if json.dumps(asdict(sig),sort_keys=True)!=json.dumps(old,sort_keys=True):raise ArithmeticError('finite point cloud signature changed')
        for row in sig.rows:insert(pivots,row)
    indices=sorted(pivots)
    if indices!=data['independent_column_indices'] or [data['points'][i] for i in indices]!=data['independent_points']:raise ArithmeticError('independent subset changed')
    old=data['rank_certificate'];actual=cert.checked_rank(model,[points[i] for i in indices],[r['prime'] for r in old['signatures']],old['no_rational_2_torsion_prime'])
    if json.dumps(actual,sort_keys=True)!=json.dumps(old,sort_keys=True) or len(indices)!=data['rank_lower_bound']:raise ArithmeticError('standalone rank proof changed')
    print('REPLAYED RECORDED-POINT RANK',data['family'],data['parameter'],'rank >=',len(indices),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-sha256');p.add_argument('--input',type=Path);p.add_argument('--output',type=Path);p.add_argument('--prime-bound',type=int,default=1000);p.add_argument('--check',type=Path);a=p.parse_args()
    check(a.check) if a.check else build(a.input,a.output,a.prime_bound,a.input_sha256)
