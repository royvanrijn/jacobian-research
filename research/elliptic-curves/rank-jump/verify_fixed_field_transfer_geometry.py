#!/usr/bin/env python3
"""Portable polynomial, cubic irreducibility and Hurwitz arithmetic replay."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
import retrospective as r
import fixed_field_transfer_geometry as source
from verify_fresh_symbolic_discriminant import trim,mul,value

OUTPUT=r.OUT/'rank_jump_fixed_field_transfer_geometry_verification_v1.json'


def remainder(a,b):
    a=trim(a);b=trim(b);assert b!=[0]
    while len(a)>=len(b) and a!=[0]:
        k=len(a)-len(b);v=a[-1]/b[-1]
        for j,c in enumerate(b):a[j+k]-=v*c
        a=trim(a)
    return a


def gcd(a,b):
    a,b=trim(a),trim(b)
    while b!=[0]:a,b=b,remainder(a,b)
    return [x/a[-1] for x in a]


def derivative(a):return trim([i*a[i] for i in range(1,len(a))] or [Q(0)])


def compute():
    data=r.read(source.OUTPUT);inputs=r.read(source.INPUT)
    for obj in (data,r.read(source.PROVENANCE)):
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha
    results=[]
    for row in data['rows']:
        assert row['status']=='PASS';inp=next(x for x in inputs['families'] if x['family']==row['family'])
        A=list(map(Q,inp['A']));B=list(map(Q,inp['B']));a3=mul(mul(A,A),A);b2=mul(B,B)
        D=trim([-64*(a3[i] if i<len(a3) else 0)-432*(b2[i] if i<len(b2) else 0) for i in range(max(len(a3),len(b2)))])
        assert D==list(map(Q,row['discriminant_ascending'])) and len(D)==25
        assert len(A)<=9 and len(B)<=13 and gcd(A,D)==[1]
        product=[Q(row['discriminant_unit'])];seen=[];odd=0
        for factor in row['squarefree_factors']:
            q=list(map(Q,factor['coefficients_ascending']));e=factor['multiplicity'];degree=len(q)-1
            assert degree==factor['degree'] and e>0 and gcd(q,derivative(q))==[1]
            assert all(gcd(q,p)==[1] for p in seen);seen.append(q)
            for _ in range(e):product=mul(product,q)
            if e%2:odd+=degree
        assert product==D
        infA=A[8] if len(A)>8 else Q(0);infB=B[12] if len(B)>12 else Q(0)
        infD=-16*(4*infA**3+27*infB**2)
        assert [infA,infB,infD]==list(map(Q,row['infinity_coefficients'])) and infD==D[24]!=0
        witness=row['cubic_irreducibility_witness'];t=Q(inp['irreducibility_witness_parameter'])
        assert t==Q(witness['parameter']) and value(D,t)!=0
        p=witness['prime'];assert p>2 and all(p%d for d in range(2,isqrt(p)+1))
        a,b=value(A,t),value(B,t);ap=a.numerator*pow(a.denominator,-1,p)%p;bp=b.numerator*pow(b.denominator,-1,p)%p
        assert ap==witness['A_mod_p'] and bp==witness['B_mod_p'] and (4*ap**3+27*bp**2)%p
        assert [x for x in range(p) if (x**3+ap*x+bp)%p==0]==witness['roots_mod_p']==[]
        assert odd==row['geometric_transposition_branch_points'] and odd>0 and odd%2==0
        # Inertia transposition: contributions 1,1,3 for degrees 3,2,6.
        for degree,contribution,key in [(3,1,'root_cover_genus'),(2,1,'quadratic_resolvent_genus'),(6,3,'constant_cubic_isomorphism_cover_genus')]:
            assert 2*row[key]-2==-2*degree+odd*contribution
        assert row['constant_field_basechange_genus_lower_bound']==row['constant_cubic_isomorphism_cover_genus']>1
        assert row['constant_field_basechange_degree_lower_bound']==6
        results.append({'family':row['family'],'branch_points':odd,'root_genus':row['root_cover_genus'],
            'resolvent_genus':row['quadratic_resolvent_genus'],'constant_field_carrier_genus':row['constant_cubic_isomorphism_cover_genus'],
            'irreducibility_prime':p})
    files=(Path(__file__),source.INPUT,source.OUTPUT,source.PROVENANCE,Path(__file__).with_name('verify_fresh_symbolic_discriminant.py'))
    return {'schema':'rank-jump.fixed-field-transfer-geometry-verification.v1','status':'PASS','rows':results,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'boundary':'Exact multiplicative support, infinity, specialization irreducibility and Hurwitz arithmetic. Monodromy and constant-field carrier arguments are in the proof note; no rank or solubility assertion is computed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],result['rows'],flush=True)
