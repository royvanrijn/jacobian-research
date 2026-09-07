#!/usr/bin/env python3
"""Portable rational polynomial and Rabin replay; no Sage dependency."""
import argparse
from fractions import Fraction as Q
from math import gcd,prod,isqrt
from pathlib import Path
import retrospective as r
import fresh_symbolic_discriminant as source
import verify_surface_discriminant_modular as finite

OUTPUT=r.OUT/'rank_jump_fresh_symbolic_discriminant_verification_v1.json'


def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def mul(a,b):
    c=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)


def value(a,t):
    result=Q(0)
    for c in reversed(a):result=result*t+c
    return result


def compute():
    data=r.read(source.OUTPUT);inputs=r.read(source.INPUT)
    for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    rows=[];modular_factor_count=0;degeneration=None
    for row in data['rows']:
        assert row['status']=='PASS';fam=next(x for x in inputs['families'] if x['family']==row['family'])
        A=list(map(Q,fam['A']));B=list(map(Q,fam['B']));a3=mul(mul(A,A),A);b2=mul(B,B)
        D=trim([-64*(a3[i] if i<len(a3) else 0)-432*(b2[i] if i<len(b2) else 0) for i in range(max(len(a3),len(b2)))])
        assert D==list(map(Q,row['discriminant_coefficients_ascending'])) and len(D)-1==row['discriminant_degree']==24
        product=[Q(row['unit'])];factor_degrees=[];irr=[]
        for factor in row['factors']:
            q=list(map(Q,factor['coefficients_ascending']));primitive=list(map(int,factor['primitive_integer_coefficients_ascending']))
            assert all(q[i]*primitive[-1]==q[-1]*primitive[i] for i in range(len(q)))
            assert gcd(*primitive)==1 and len(q)-1==factor['degree']
            for _ in range(factor['exponent']):product=mul(product,q)
            possible=set(range(1,len(q)-1));certificate=[]
            for mod in factor['modular']:
                p=mod['prime'];assert p>1 and all(p%d for d in range(2,isqrt(p)+1))
                pp=[1];factors=mod['factors_ascending'];assert len({tuple(f) for f in factors})==len(factors)
                for f in factors:
                    assert finite.irreducible(f,p);pp=finite.mul(pp,f,p);modular_factor_count+=1
                assert pp==[(c*pow(primitive[-1],-1,p))%p for c in primitive]
                degrees=[len(f)-1 for f in factors];assert degrees==mod['factor_degrees']
                sums={0}
                for d in degrees:sums|={s+d for s in list(sums)}
                possible&=sums-{0,len(q)-1}
                assert sorted(possible)==mod['remaining_possible_proper_factor_degrees']
                certificate.append({'prime':p,'degrees':degrees})
            assert not possible and factor['modular_irreducibility_status']=='PASS'
            factor_degrees.append([factor['degree'],factor['exponent']]);irr.append(certificate)
        assert product==D
        specializations=[]
        for rec in row['specializations']:
            inp=next(x for x in inputs['cases'] if x['token']==rec['token']);t=Q(inp['parameter']);d=t.denominator
            rawA=value(A,t)*d**8;rawB=value(B,t)*d**12;model,_=r.short(inp['model'],[]);model=list(map(Q,model))
            u=Q(rec['scale_to_frozen_short_model']);assert rawA==u**4*model[3] and rawB==u**6*model[4]
            assert value(D,t)*d**24==u**12*(-16*(4*model[3]**3+27*model[4]**2))!=0
            N=int(inp['unresolved_cofactor']);parts=[N]
            for f,claim in zip(row['factors'],rec['factor_value_gcds']):
                q=list(map(Q,f['coefficients_ascending']));v=(value(q,t)*d**f['degree']).numerator
                assert str(v)==claim['homogeneous_numerator'] and str(gcd(N,v))==claim['gcd']
                new=[]
                for n in parts:
                    g=gcd(n,v);new.extend([g,n//g] if 1<g<n else [n])
                parts=new
            assert parts==list(map(int,rec['cofactor_parts'])) and len(parts)-1==rec['proper_split_count']==0
            # Audit constant normalization too, so no omitted content split is hidden.
            content=Q(row['unit']);cg=[gcd(N,content.numerator),gcd(N,content.denominator)]
            assert cg==[1,1]
            specializations.append({'token':rec['token'],'proper_cofactor_splits':0,'unit_gcds':cg})
        if row['family']=='a1-fibration-01':
            assert factor_degrees==[[1,2],[22,1]] and row['factors'][0]['coefficients_ascending']==['2','1']
            t0=Q(-2);a0=value(A,t0);b0=value(B,t0);H=list(map(Q,row['factors'][1]['coefficients_ascending']))
            assert a0!=0 and value(H,t0)!=0 and value(D,t0)==0
            node=-3*b0/(2*a0);assert a0==-3*node**2 and b0==2*node**3
            degeneration={'parameter':'-2','discriminant_order':2,'c4_at_parameter':str(-48*a0),
                'A_at_parameter':str(a0),'B_at_parameter':str(b0),'double_cubic_root':str(node),
                'residual_degree22_factor_value':str(value(H,t0)),
                'linear_homogeneous_values':[{ 'token':x['token'],'n_plus_2d':Q(x['parameter']).numerator+2*Q(x['parameter']).denominator}
                    for x in inputs['cases'] if x['family']==row['family']],
                'interpretation':'Order-two discriminant with c4 a unit: multiplicative I2 fibre in characteristic zero. This singular parameter is absent from both retained smooth fibres.'}
        rows.append({'family':row['family'],'factor_degrees_and_multiplicities':factor_degrees,
            'modular_certificates':irr,'specializations':specializations})
    files=(Path(__file__),source.OUTPUT,source.INPUT,Path(finite.__file__),Path(r.__file__))
    return {'schema':'rank-jump.fresh-symbolic-discriminant-verification.v1','status':'PASS',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'rows':rows,'irreducible_finite_factors_verified':modular_factor_count,'MW16_linear_degeneration':degeneration,
        'boundary':'Exact discriminant divisor and failed cofactor-splitting route. This does not test auxiliary-cover solubility or exclude rank jumps on smooth fibres.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],[(x['family'],x['factor_degrees_and_multiplicities']) for x in result['rows']],flush=True)
