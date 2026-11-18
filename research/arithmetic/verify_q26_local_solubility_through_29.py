#!/usr/bin/env python3
"""Exact direct Q5 certificate and small-prime affine Hensel witnesses."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import flint

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def require(test,message):
    if not test:raise AssertionError(message)

def evaluate(a,x):
    out=0
    for c in reversed(a):out=out*x+c
    return out

def vp(n,p):
    require(n!=0,'Valuation of zero')
    k=0
    while n%p==0:n//=p;k+=1
    return k

def direct_p5(coeff,data):
    f=flint.fmpz_poly([evaluate(a,data['t_base']) for a in coeff])
    transformed=f(flint.fmpz_poly([data['s_base'],data['s_step']]))
    raw=[int(c) for c in transformed]
    values=[vp(c,5) for c in raw]
    require(values==data['raw_coefficient_valuations'],'Exact p5 coefficient valuations')
    exponent=min(values)
    require(exponent==data['exponent']==40,'Normalization exponent')
    require(all(c%5**exponent==0 for c in raw),'Nonintegral p5 coefficient')
    normalized=[c//5**exponent for c in raw]
    require(normalized==data['normalized_coefficients_ascending'],'Direct identity binding')
    residue=[c%5 for c in normalized]
    require(residue==data['reduction_mod5_ascending']==[0,3,4]+[0]*10,
            'Entire p5 polynomial reduction')
    derivative=[j*normalized[j] for j in range(1,len(normalized))]
    require(data['t_family_step_exponent']==exponent+1==41,'Uniform slope step')
    controls=[]
    for root,der in zip(data['simple_roots'],data['derivatives']):
        require(evaluate(normalized,root)%5==0,'Residue root')
        require(evaluate(derivative,root)%5==der!=0,'Unit derivative')
        x=root;digits=[]
        for precision in range(1,9):
            modulus=5**precision
            require(evaluate(normalized,x)%modulus==0,'Hensel digit')
            digits.append({'precision':precision,'X_mod_5_power':x})
            candidates=[x+d*modulus for d in range(5)
                        if evaluate(normalized,x+d*modulus)%(modulus*5)==0]
            require(len(candidates)==1,'Unique Hensel lift')
            x=candidates[0]
        controls.append({'residue_root':root,'derivative_mod5':der,'digits':digits})
    return {'identity':'P(X)=F(1056+1200X,238)/5^40',
            'degree':len(raw)-1,'raw_coefficient_valuations':values,
            'entire_reduction_ascending':residue,'simple_root_controls':controls,
            'uniform_family':'5^-40 F(1056+1200X,238+5^41 V)',
            'positive_V_degree_coefficients_divisible_by_5':True,
            'proof_of_uniformity':'An integer binomial expansion gives factor5^(41*j) in every V^j term; for j>=1 normalization by5^40 leaves a factor5.'}

def affine_residue(coeff,row):
    p,s,t=[row[k] for k in ('p','s','t')]
    values=[evaluate(a,t)%p for a in coeff]
    f=evaluate(values,s)%p
    fs=evaluate([i*values[i] for i in range(1,len(values))],s)%p
    dt=[evaluate([j*a[j] for j in range(1,len(a))],t)%p for a in coeff]
    ft=evaluate(dt,s)%p
    require(f==0,'Affine residue is not a curve point')
    require((fs,ft)==(row['F_s'],row['F_t']),'Affine residue derivative mismatch')
    require(fs or ft,'Singular residue only')
    return {**row,'F':f,'simple_Hensel_coordinate':'s' if fs else 't'}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-receipt',type=Path,
                        default=HERE/'q26_local_solubility_through_29_verification.json')
    args=parser.parse_args()
    data=json.loads((HERE/'q26_local_solubility_through_29_inputs.json').read_text())
    coeff=json.loads((HERE/'q26_high_genus_affine_trace_curve_inputs.json').read_text())['curve']['coefficients_ascending_s']
    p5=direct_p5(coeff,data['p5'])
    table=[affine_residue(coeff,row) for row in data['simple_affine_table']]
    p29=affine_residue(coeff,data['prime29_consistency_point'])
    primes=[p for p in range(2,30) if all(p%d for d in range(2,math.isqrt(p)+1))]
    covered=sorted(set(data['imported_primes'])|{5}|{r['p'] for r in table})
    require(covered==primes==data['all_primes_through_29'],'Incomplete small-prime coverage')
    paths=[Path(__file__).resolve(),HERE/'Q26_LOCAL_SOLUBILITY_THROUGH_PRIME_TWENTY_NINE.md',
           HERE/'q26_local_solubility_through_29_inputs.json',
           HERE/'q26_high_genus_affine_trace_curve_inputs.json',
           HERE/'Q26_HIGH_GENUS_AFFINE_TRACE_CURVE.md',
           HERE/'Q26_SMALL_PRIME_LOCAL_BRANCHES.md',
           HERE/'q26_small_prime_local_branches_verification.json',
           HERE/'ODD_AFFINE_CRITICAL_TRACE_CURVES.md',HERE/'sources/milne-ant300.pdf']
    result={'schema_version':1,'result':'PASS','python':platform.python_version(),
            'python_flint':flint.__version__,'direct_p5_certificate':p5,
            'simple_affine_table':table,'prime29_consistency_point':p29,
            'all_primes_through_29':covered,'real_place_imported':True,
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'scope':'Infinitely many good local points at every prime<=29 and over R. No Q-rational point, all-prime local solubility, Hasse principle, all-critical-real statement or exhaustive local classification.',
            'formal_verification':False,'independent_external_replay':False,
            'external_human_review':False}
    args.write_receipt.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: direct Q5 Hensel family and complete local-prime coverage through29.')

if __name__=='__main__':main()
