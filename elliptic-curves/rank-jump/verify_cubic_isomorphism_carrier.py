#!/usr/bin/env python3
"""Portable interpolation, norm and discriminant-squareclass verification."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt
import retrospective as r
import cubic_isomorphism_carrier as source
from verify_fresh_symbolic_discriminant import value
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_cubic_isomorphism_carrier_verification_v1.json'


def det(M):
    return M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])


def compute():
    data=r.read(source.OUTPUT);inp=r.read(source.INPUT)
    for obj in (data,inp):
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha
    checks=0
    for row in data['families']:
        family=next(x for x in inp['families'] if x['family']==row['family']);t=Q(row['fixed_field_parameter'])
        a=value(list(map(Q,family['A'])),t);b=value(list(map(Q,family['B'])),t)
        assert list(map(Q,row['fixed_cubic_ascending']))==[b,a,0,1]
        K=Algebra(row['fixed_cubic_ascending']);base=[K.elt([1,0,0]),K.elt([0,1,0]),K.elt([0,0,1])]
        ac=list(map(Q,row['A_equation_v2_vw_w2']));bc=list(map(Q,row['B_equation_v3_v2w_vw2_w3']));ic=list(map(Q,row['index_v3_v2w_vw2_w3']))
        assert ac==[a,3*b,-a*a/3] and bc==[b,-2*a*a/3,-a*b,-2*a**3/27-b*b] and ic==[1,0,a,b]
        # Each asserted identity has separate v,w degree <=6. This 7x7 grid
        # certifies the entire two-variable identity for each fixed field.
        for v in range(7):
            for w in range(7):
                q=K.elt([2*a*w/3,v,w]);M=list(map(list,zip(*[K.mul(q,e) for e in base])))
                assert sum(M[i][i] for i in range(3))==0
                aa=sum(ac[i]*v**(2-i)*w**i for i in range(3));bb=sum(bc[i]*v**(3-i)*w**i for i in range(4))
                assert -sum(M[i][j]*M[j][i] for i in range(3) for j in range(3))/2==aa
                assert -det(M)==-K.norm(q)==bb
                index=sum(ic[i]*v**(3-i)*w**i for i in range(4))
                assert det(list(map(list,zip(base[0],q,K.mul(q,q)))))==index
                assert -4*aa**3-27*bb**2==(-4*a**3-27*b*b)*index**2;checks+=1
        assert row['known_identity_lift']=={'t':str(t),'v':'1','w':'0'}
    models={x['token']:x['model'] for x in inp['models']};discriminants={};witnesses=[]
    for token,model in models.items():
        s,_=r.short(model,[]);a,b=map(Q,s[3:]);D=-4*a**3-27*b*b;assert D
        gal=r.galois(s);assert gal['galois_group']=='S3'
        witnesses.append({'token':token,'irreducibility_prime':gal['irreducibility_prime']});discriminants[token]=D
    assert len(inp['pairs'])==len(data['pairs'])==8
    for pair,row in zip(inp['pairs'],data['pairs']):
        assert all(pair[k]==row[k] for k in ('high','low'))
        ds=[discriminants[pair[k]] for k in ('high','low')];ratio=ds[0]/ds[1]
        assert list(map(Q,row['cubic_discriminants']))==ds and Q(row['discriminant_ratio'])==ratio
        n,d=ratio.numerator,ratio.denominator;sn=isqrt(abs(n));sd=isqrt(d)
        assert str(sn)==row['absolute_numerator_floor_sqrt'] and str(sd)==row['denominator_floor_sqrt']
        assert sn*sn<=abs(n)<(sn+1)**2 and sd*sd<=d<(sd+1)**2
        assert n<0 or sn*sn!=n or sd*sd!=d
        assert row['same_quadratic_resolvent'] is False and row['cubic_fields_isomorphic'] is False
    files=(Path(__file__),source.INPUT,source.OUTPUT,Path(r.__file__),Path(__file__).with_name('verify_unpointed_governing_norm.py'))
    return {'schema':'rank-jump.cubic-isomorphism-carrier-verification.v1','status':'PASS',
        'carrier_polynomial_identity_grid_checks':checks,'numeric_carriers_verified':len(data['families']),
        'nonisomorphic_matched_pairs':8,'individual_cubic_S3_witnesses':witnesses,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'boundary':'Exact identities on all five specified fields and complete eight-pair nonisomorphism test. No rational-point enumeration on the carriers.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['carrier_polynomial_identity_grid_checks'],'identity checks;',result['nonisomorphic_matched_pairs'],'nonisomorphic pairs',flush=True)
