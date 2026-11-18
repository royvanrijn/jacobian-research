#!/usr/bin/env python3
"""Exact toric Hensel charts at2,3 and a real Sturm gate for the Q26 curve."""
from __future__ import annotations
import argparse
from functools import reduce
import hashlib
import json
import math
from pathlib import Path
import platform
import flint

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHARTS=[{'prime':2,'m':0,'a':0,'b':1,'root':1},
        {'prime':3,'m':18,'a':-2,'b':2,'root':2}]

def require(test,message):
    if not test:raise AssertionError(message)

def vp(n,p):
    require(n!=0,'Valuation of zero')
    k=0
    while n%p==0:n//=p;k+=1
    return k

def sign(x):return 1 if x>0 else -1 if x<0 else 0

def variations(signs):
    ss=[s for s in signs if s]
    return sum(a!=b for a,b in zip(ss,ss[1:]))

def evaluate(coeff,x):
    out=0
    for c in reversed(coeff):out=out*x+c
    return out

def chart_certificate(terms,chart):
    p,m,a,b,root=[chart[k] for k in ('prime','m','a','b','root')]
    transformed=[];valuations=[]
    for i,k,c in terms:
        power=m+a*i+b*k
        valuations.append(vp(c,p)+power)
        if power>=0:d=c*p**power
        else:
            divisor=p**(-power)
            require(c%divisor==0,'Nonintegral toric coefficient')
            d=c//divisor
        transformed.append((i,k,d))
    require(min(valuations)==0,'Toric primitive valuation')
    reduction=sorted((i,k,c%p) for i,k,c in transformed if c%p)
    require(reduction==[(11,0,1),(12,0,1)],'Full global chart reduction')
    value=(pow(root,11,p)+pow(root,12,p))%p
    derivative=(11*pow(root,10,p)+12*pow(root,11,p))%p
    require(value==0 and derivative!=0,'Simple residue root')
    # Concrete finite precision controls atV1 supplement the uniform Hensel proof.
    specialized=[sum(c for j,k,c in transformed if j==i) for i in range(13)]
    dcoeff=[i*specialized[i] for i in range(1,13)]
    lifts=[];x=root
    for precision in range(1,9):
        modulus=p**precision
        require(evaluate(specialized,x)%modulus==0,'Finite Hensel approximation')
        require(evaluate(dcoeff,x)%p!=0,'Unit derivative lost')
        lifts.append({'precision':precision,'X_mod_p_power':x})
        candidates=[x+digit*modulus for digit in range(p)
                    if evaluate(specialized,x+digit*modulus)%(modulus*p)==0]
        require(len(candidates)==1,'Nonunique simple Hensel digit')
        x=candidates[0]
    text=json.dumps(transformed,separators=(',',':')).encode()
    return {**chart,'G_terms':len(terms),'transformed_coefficients_all_integral':True,
            'minimum_transformed_coefficient_valuation':min(valuations),
            'full_reduction_terms_i_k_coefficient':reduction,
            'partial_X_at_residue_root':derivative,
            'transformed_coefficients_sha256':hashlib.sha256(text).hexdigest(),
            'unit_parameter_valuations':{'t':-b,'s':a-b},
            'V1_unique_Hensel_controls':lifts,
            'uniformity':'For every V in Zp, one unique X congruent root modp; finite tests supplement written Hensel proof.'}

def sturm_certificate(coeff):
    content=math.gcd(*coeff)
    aa=[c//content for c in coeff]
    if aa[-1]<0:aa=[-c for c in aa]
    A=flint.fmpq_poly(aa)
    require(A.degree()==12 and A.gcd(A.derivative()).degree()==0,'Simple infinity roots')
    chain=[A,A.derivative()]
    while chain[-1].degree()>0:
        rem=-(chain[-2]%chain[-1])
        require(not rem.is_zero(),'Nontrivial Sturm gcd')
        chain.append(rem)
    plus=[sign(f[f.degree()]) for f in chain]
    minus=[sign(f[f.degree()])*(-1)**f.degree() for f in chain]
    count=variations(minus)-variations(plus)
    require(count==2,'Real infinity root count')
    return {'primitive_leading_polynomial_ascending':aa,'coefficient_content_removed':content,
            'degree':12,'gcd_with_derivative_degree':0,
            'Sturm_chain_degrees':[f.degree() for f in chain],
            'signs_at_negative_infinity':minus,'signs_at_positive_infinity':plus,
            'variations_negative':variations(minus),'variations_positive':variations(plus),
            'simple_real_roots':count}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-receipt',type=Path,
                        default=HERE/'q26_small_prime_local_branches_verification.json')
    args=parser.parse_args()
    data=json.loads((HERE/'q26_high_genus_affine_trace_curve_inputs.json').read_text())['curve']
    coeff=data['coefficients_ascending_s']
    require([len(a)-1 for a in coeff]==list(range(144,131,-1)),'Infinity degree binding')
    terms=[(i,144-i-j,c) for i,a in enumerate(coeff) for j,c in enumerate(a) if c]
    require(all(k>=0 for i,k,c in terms),'Nonpolynomial infinity transform')
    leading=[a[-1] for a in coeff]
    require([c%2 for c in leading]==[0]*11+[1,1],'Q2 leading residue binding')
    result={'schema_version':1,'result':'PASS','python':platform.python_version(),
            'python_flint':flint.__version__,
            'toric_charts':[chart_certificate(terms,c) for c in CHARTS],
            'real_infinity':sturm_certificate(leading),
            'naive_projective_point':{'S_T_U':[1,0,0],
                'exact_root':True,'minimum_local_T_U_degree':144-12,
                'normalization_point_not_inferred_from_this_alone':True},
            'conclusion':'Infinitely many smooth good normalized affine-trace models over Q2,Q3,R; no Q-rational point or all-critical-real claim.',
            'formal_verification':False,'independent_external_replay':False,
            'external_human_review':False}
    paths=[Path(__file__).resolve(),HERE/'Q26_SMALL_PRIME_LOCAL_BRANCHES.md',
           HERE/'q26_high_genus_affine_trace_curve_inputs.json',
           HERE/'Q26_HIGH_GENUS_AFFINE_TRACE_CURVE.md',
           HERE/'ODD_AFFINE_CRITICAL_TRACE_CURVES.md',HERE/'sources/milne-ant300.pdf']
    result['source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    args.write_receipt.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: exact Q2/Q3 toric Hensel charts and two simple real infinity branches.')

if __name__=='__main__':main()
