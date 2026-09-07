#!/usr/bin/env python3
"""Exact determinant norms and private odd-prime witnesses from the ZIP."""
import argparse
from fractions import Fraction as F
from itertools import combinations,product
from math import gcd,isqrt,lcm,prod
from pathlib import Path
import json
import zipfile
import retrospective as r
import package_three_radical_evidence as archive

OUTPUT=r.OUT/'rank_jump_three_radical_incidence_verification_v1.json'


def normalize(cs):
    d=lcm(*(x.denominator for x in cs));integers=[int(c*d) for c in cs]
    content=gcd(*integers)
    if integers[-1]<0:content=-content
    return [v//content for v in integers]


def norm(A,B,cs):
    columns=[]
    for j in range(3):
        row=[0]*5
        for i,c in enumerate(cs):row[i+j]=c
        for k in (4,3):row[k-2]-=A*row[k];row[k-3]-=B*row[k];row[k]=0
        columns.append(row[:3])
    a,b,c=columns
    return a[0]*(b[1]*c[2]-b[2]*c[1])-b[0]*(a[1]*c[2]-a[2]*c[1])+c[0]*(a[1]*b[2]-a[2]*b[1])


def check_norm(A,B,delta,coefficients,record):
    cs=normalize(coefficients);assert list(map(str,cs))==record['primitive_coefficients']
    value=norm(A,B,cs);assert str(value)==record['norm'] and value!=0
    remaining=abs(value)
    for s in record['discriminant_gcd_steps']:
        g=int(s);assert g==gcd(remaining,abs(delta))>1;remaining//=g
    assert gcd(remaining,delta)==1 and str(remaining)==record['outside_discriminant_norm']
    square=isqrt(remaining)**2==remaining
    assert square==record['outside_norm_is_square']
    if square:assert int(record['outside_norm_square_root'])**2==remaining
    return remaining


def symbolic_norm():
    def mul(a,b):
        result={}
        for x,c in a.items():
            for y,d in b.items():
                e=tuple(i+j for i,j in zip(x,y));result[e]=result.get(e,0)+c*d
        return {e:c for e,c in result.items() if c}
    result={(0,0,0):1}
    for a,b,c in product((-1,1),repeat=3):
        result=mul(result,{(1,0,0):a,(0,1,0):b,(0,0,1):c})
    heron={(4,0,0):1,(0,4,0):1,(0,0,4):1,(2,2,0):-2,(2,0,2):-2,(0,2,2):-2}
    assert result==mul(heron,heron)
    return {'formal_radical_norm_terms':len(result),'identity':'Product over all eight signs of (+/-u +/-v +/-w) equals (u^4+v^4+w^4-2u^2v^2-2u^2w^2-2v^2w^2)^2.'}


def verify():
    manifest=r.read(archive.MANIFEST);assert r.digest(archive.BUNDLE.read_bytes())==manifest['bundle_sha256']
    assert manifest['packager_sha256']==r.digest(Path(archive.__file__).read_bytes())
    with zipfile.ZipFile(archive.BUNDLE) as z:
        raw={entry['path']:z.read(entry['member']) for entry in manifest['members']}
    for entry in manifest['members']:assert len(raw[entry['path']])==entry['bytes'] and r.digest(raw[entry['path']])==entry['sha256']
    objects=[json.loads(raw[str(p.relative_to(r.ROOT))]) for p in archive.PATHS]
    for obj in objects:
        assert obj['status']=='PASS'
        for name,sha in obj['bindings'].items():
            data=raw[name] if name in raw else (r.ROOT/name).read_bytes();assert r.digest(data)==sha,name
    source,private=objects;cases={x['token']:x for x in r.read(archive.first.panel.INPUT)['cases']}
    rows=[];triple_count=control_count=0
    for saved,witness in zip(source['rows'],private['rows']):
        token=saved['token'];assert token==witness['token']
        raw_case=cases[token];a1,a2,a3,a4,a6=map(F,raw_case['model'])
        b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
        c4=b2*b2-24*b4;c6=-b2**3+36*b2*b4-216*b6
        A=-c4/48;B=-c6/864;scale=lcm(A.denominator,B.denominator)
        A=int(A*scale**4);B=int(B*scale**6);delta=-4*A**3-27*B*B
        assert saved['cubic_ascending']==list(map(str,[B,A,0,1])) and str(delta)==saved['discriminant']
        xs=[]
        for px,py in raw_case['generic_sections']:
            x,y=F(px),F(py)
            assert y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
            xs.append((x+b2/12)*scale**2)
        triples=list(combinations(range(len(xs)),3));assert len(triples)==len(saved['rows'])==len(witness['rows'])
        outside=[]
        for triple,record in zip(triples,saved['rows']):
            assert list(triple)==record['generic_indices']
            a,b,c=[xs[i] for i in triple]
            cs=[a*a+b*b+c*c-2*(a*b+a*c+b*c),2*(a+b+c),F(-3)]
            outside.append(check_norm(A,B,delta,cs,record));triple_count+=1
        total=prod(outside);private_bits=[]
        for i,(n,record) in enumerate(zip(outside,witness['rows'])):
            assert record['index']==i and record['generic_indices']==list(triples[i])
            q=int(record['private_norm_factor']);assert q>0 and n%q==0
            assert gcd(q,total//n)==1 and gcd(q,delta)==1
            assert record['private_factor_is_nonsquare']==(isqrt(q)**2!=q)
            # Odd places are essential for the Selmer conclusion. The initial
            # worker strips disc(f), which is odd on case05. Remove two here
            # as well, and certify a nonsquare odd private part on every row.
            while q%2==0:q//=2
            assert isqrt(q)**2!=q
            private_bits.append(q.bit_length())
        assert witness['ramification_rank_lower_bound']==len(triples)
        assert witness['candidate_span_unramified_dimension_upper_bound']==0
        assert witness['new_Selmer_dimension_after_adding_generic_subgroup']==0
        for i,(x,record) in enumerate(zip(xs,saved['generic_controls'])):
            assert record['generic_index']==i
            n=check_norm(A,B,delta,[x,F(-1)],record);assert isqrt(n)**2==n;control_count+=1
        n=check_norm(A,B,delta,[F(A),F(0),F(3)],saved['derivative_false_positive_control'])
        assert isqrt(n)**2==n;control_count+=1
        rows.append({'token':token,'triples':len(triples),'independent_private_odd_ramification_witnesses':len(private_bits),
            'private_odd_factor_bit_range':[min(private_bits),max(private_bits)],
            'candidate_span_new_Selmer_dimension':0})
        print(token,'PASS',len(triples),'private odd witnesses',flush=True)
    assert triple_count==10520 and control_count==285
    return {'schema':'rank-jump.three-radical-incidence-verification.v1','status':'PASS',
        'symbolic_norm':symbolic_norm(),'rows':rows,'verified_triple_norms':triple_count,
        'verified_controls':control_count,'verified_private_odd_ramification_witnesses':triple_count,
        'method':'Independent rational Weierstrass invariants and primitive coefficients, cubic companion-matrix determinants, gcd support checks, and coprimality of every odd private witness with the product of all other norms. No factoring or local number-field computation.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
                    (Path(__file__),archive.BUNDLE,archive.MANIFEST,archive.first.panel.INPUT,Path(archive.__file__),Path(r.__file__))}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=verify()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS all additive norms and complete block ramification ranks')
