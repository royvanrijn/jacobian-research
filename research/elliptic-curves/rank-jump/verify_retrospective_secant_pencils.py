#!/usr/bin/env python3
"""Independent rational polynomial and generic-source checks for fitted secants."""
import argparse
from pathlib import Path
import retrospective as r
import bad_prime_support as bad
import retrospective_secant_pencils as pencils

OUTPUT=r.OUT/'rank_jump_retrospective_secant_verification_v1.json'


def trim(a):
    a=list(a)
    while len(a)>1 and not a[-1]:a.pop()
    return a


def add(*terms):
    out=[r.F(0)]*max(map(len,terms))
    for term in terms:
        for i,a in enumerate(term):out[i]+=a
    return trim(out)


def mul(a,b):
    out=[r.F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)


def scale(a,s):return trim([s*x for x in a])


def compute():
    from sage.all import QQ,CartanMatrix
    data=r.read(pencils.INPUT);assert data['source_sha256']==r.digest(r.INPUT.read_bytes())
    result=r.read(pencils.OUTPUT);assert result['bindings']==pencils.bindings()
    inverse=CartanMatrix(['D',7]).inverse()
    ends=[inverse[i,i] for i in [0,5,6]]
    assert ends==[1,QQ(7)/4,QQ(7)/4]
    rows=[]
    for row in result['rows']:
        assert row['status']=='PASS';i=row['case_index'];old=next(x for x in data['rows'] if x['case_index']==i)
        source=bad.cases()[i];model,points=r.short(source['model'],source['generic_points'][:2])
        assert model==old['short_model'] and points==old['generic_pair']
        profile,_,sigs=r.characterize(source)
        assert sigs[:2]==old['generic_pair_fingerprints'] and r.rank(sigs[:2])==2
        A,B=map(r.F,model[3:]);s=r.F(row['secant_slope']);b=r.F(row['secant_intercept'])
        for x,y in map(lambda p:map(r.F,p),points):assert s*x+b==y
        roots=list(map(r.F,row['residual_roots']));product=[r.F(1)]
        for root in roots:product=mul(product,[-root,r.F(1)])
        assert product==[B-b*b,A-2*s*b,-s*s,r.F(1)]
        assert product==list(map(r.F,row['residual_cubic_coefficients']))
        x0=-b/s;C=x0**3+A*x0+B
        assert str(x0)==row['fixed_section_x'] and str(C)==row['fixed_section_y_squared']
        assert bool(QQ(str(C)).is_square())==row['fixed_section_is_rational']
        a2=[-s*s,s*s];a4=[A-2*s*b,2*s*b];a6=[B-b*b,b*b]
        c4=scale(add(mul(a2,a2),scale(a4,-3)),16)
        delta=scale(add(mul(mul(a2,a2),mul(a4,a4)),scale(mul(mul(a4,a4),a4),-4),
            scale(mul(mul(mul(a2,a2),a2),a6),-4),scale(mul(a6,a6),-27),
            scale(mul(mul(a2,a4),a6),18)),16)
        assert c4==list(map(r.F,row['parent_c4_coefficients']))
        assert delta==list(map(r.F,row['parent_discriminant_coefficients']))
        assert len(delta)==4 and len(c4)==3
        rows.append({'case_index':i,'status':'PASS','generic_pair_source_verified':True,
            'independent_pair_fingerprint_rank':2,'original_quotient_contribution':0,
            'exact_polynomial_invariants_verified':True,'fixed_section_is_rational':row['fixed_section_is_rational'],
            'original_generic_rank':old['original_generic_rank'],
            'original_observed_quotient_rank':old['original_observed_quotient_rank'],
            'apparent_known_jump_over_fitted_base':old['original_known_independent_rank']-row['base_arithmetic_generic_rank']})
    return {'schema':'rank-jump.retrospective-secant-verification.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),pencils.INPUT,pencils.OUTPUT)},
        'D7_nonidentity_simple_component_corrections':list(map(str,ends)),
        'fixed_section_height_lower_bound':'1/4','rows':rows}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();data=compute()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS independent secant-pencil verification')
    else:r.write_new(OUTPUT,data)
