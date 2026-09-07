#!/usr/bin/env python3
"""Retrospective character reweighting of a frozen complete Frobenius ledger."""
import argparse
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path
from sage.all import GF, PolynomialRing
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'PAIRED_CHARACTER_MOMENTS_PROTOCOL.json'
CASES=r.OUT/'rank_jump_soluble_quartet_compression_inputs_v1.json'
TRACES=r.OUT/'rank_jump_native_twist_frobenius_v1.json'
VERIFIED=r.OUT/'rank_jump_native_twist_frobenius_verification_v2.json'
PARITY=r.OUT/'rank_jump_native_twist_moment_parity_v1.json'
OUTPUT=r.OUT/'rank_jump_paired_character_moments_v1.json'


def char(p,x):
    x%=p
    if not x:return 0
    ans=pow(x,(p-1)//2,p);assert ans in(1,p-1)
    return 1 if ans==1 else -1


def qchar(form,n,a,b):
    p=131;c0,c1,c2=[x%p for x in form]
    if a==-1:u,v=c2,0
    else:u,v=(c0+c1*a+c2*(a*a+2*b*b))%p,(c1*b+2*c2*a*b)%p
    return char(p,u) if n==1 else char(p,u*u-2*v*v)


def compute():
    cases=r.read(CASES);base=r.read(TRACES);verify=r.read(VERIFIED);parity=r.read(PARITY)
    for data in (cases,base,verify,parity):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    roster={c['label']:c['form'] for case in cases['cases'] for c in case['covers']};assert len(roster)==12
    targets=[]
    for case in cases['cases']:
        labels=[c['label'] for c in case['covers']]
        for mask in range(1,16):
            targets.append({'system':case['id'],'mask':mask,'labels':[label for i,label in enumerate(labels) if mask>>i&1]})
    targets.append({'system':'cross_group_FD','mask':3,'labels':['orbit-1795d','orbit-11278']});assert len(targets)==46
    R=PolynomialRing(GF(131),'t');A,B,_=[R(x) for x in base['geometry']['modular_coefficients']]
    delta=-16*(4*A**3+27*B**2)
    factors=[R(row['factor']) for row in parity['multiplicative_places']]
    assert prod(factors)==delta.monic() and all(f.is_irreducible() for f in factors)
    chars={label:[qchar(form,n,a,b) for n,a,b,mult,ap,unused in base['fibre_trace_ledger']] for label,form in roster.items()}
    norms={label:[int(f.resultant(R(form))) for f in factors] for label,form in roster.items()}
    rows=[]
    for target in targets:
        labels=target['labels'];k=len(labels);q=prod(R(roster[label]) for label in labels)
        gate={'expected_degree':2*k,'actual_degree':int(q.degree()),'squarefree':bool(q.is_squarefree()),
              'gcd_with_original_discriminant_degree':int(q.gcd(delta).degree())}
        good=gate['actual_degree']==2*k and gate['squarefree'] and gate['gcd_with_original_discriminant_degree']==0
        row={**target,'number_of_quadratics':k,'good_reduction_gate':gate,
             'arithmetic_generic_rank_lower_bound':1 if k==1 else 0}
        if not good:
            rows.append(row|{'status':'UNKNOWN_BAD_REDUCTION','arithmetic_generic_rank_upper_bound':'UNKNOWN'});continue
        traces=[0,0]
        for j,(n,a,b,mult,ap,unused) in enumerate(base['fibre_trace_ledger']):
            traces[n-1]-=mult*ap*prod(chars[label][j] for label in labels)
        N=20+4*k;s1=F(traces[0],131);s2=F(traces[1],131**2);M2=(N+s2)/2
        assert -N<=s1<N and s1*s1/N<=M2<=N
        c=(s1-M2)/(N-s1);bound=(M2-2*c*s1+N*c*c)/(1-c)**2
        upper=bound.numerator//bound.denominator
        local=[]
        for i,old in enumerate(parity['multiplicative_places']):
            nq=prod(norms[label][i] for label in labels)%131;assert nq
            local.append({'factor_index':i,'norm_of_product':nq,
                          'root_number':old['local_root_numbers'][0]*char(131,nq)})
        W=prod(x['root_number'] for x in local)
        refined=max(m for m in range(upper+1) if (-1)**m==W)
        assert refined>=row['arithmetic_generic_rank_lower_bound']
        rows.append(row|{'status':'PASS','cohomology_dimension':N,'Frobenius_traces':traces,
                         'quadratic_center':str(c),'exact_moment_bound':str(bound),'moment_integer_bound':upper,
                         'multiplicative_root_numbers':local,'additive_root_number_product':1,'global_root_number':W,
                         'arithmetic_generic_rank_upper_bound':refined})
    native=next(row for row in rows if row['labels']==['orbit-1795d'])
    assert native['Frobenius_traces']==[122,33710] and native['arithmetic_generic_rank_upper_bound']==7
    comparisons=[]
    for name,labels,solubility,auxrank in [
        ('positive_FG',['orbit-0911e','orbit-1795d'],'YES',2),
        ('cross_group_FD',['orbit-1795d','orbit-11278'],'YES',3),
        ('obstructed_AD',['orbit-030cb','orbit-11278'],'NO',2)]:
        selected=[next(row for row in rows if row['labels']==[label]) for label in labels]
        selected.append(next(row for row in rows if set(row['labels'])==set(labels)))
        upper=17+sum(row['arithmetic_generic_rank_upper_bound'] for row in selected) if all(row['status']=='PASS' for row in selected) else 'UNKNOWN'
        comparisons.append({'id':name,'labels':labels,'generic_pullback_rank_lower_bound':19,
                            'generic_pullback_rank_upper_bound':upper,'carrier_global_solubility':solubility,
                            'carrier_Jacobian_exact_rank':auxrank,
                            'twist_bounds':[{'labels':row['labels'],'lower':row['arithmetic_generic_rank_lower_bound'],'upper':row['arithmetic_generic_rank_upper_bound']} for row in selected]})
    singletons=[row for row in rows if row['number_of_quadratics']==1]
    return {'schema':'rank-jump.paired-character-moments.v1','status':'PASS','layer':'incidence',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,CASES,TRACES,VERIFIED,PARITY,Path(__file__),HERE/'retrospective.py')},
            'rows':rows,'paired_carrier_comparisons':comparisons,
            'good_rows':sum(row['status']=='PASS' for row in rows),'bad_reduction_rows':sum(row['status']!='PASS' for row in rows),
            'singleton_bounds':[{k:row[k] for k in('system','labels','status','arithmetic_generic_rank_upper_bound')} for row in singletons],
            'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);mode=p.parse_args().mode
    result=compute()
    if mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print('PASS;',result['good_rows'],'good rows;',result['bad_reduction_rows'],'bad-reduction UNKNOWN')
    for row in result['singleton_bounds']:print(row)
    for row in result['paired_carrier_comparisons']:print(row)
