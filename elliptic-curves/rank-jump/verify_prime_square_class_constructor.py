#!/usr/bin/env python3
"""Replay ideal targets and independently recompute every norm over Q."""
import argparse
from fractions import Fraction as Q
from itertools import product
from math import gcd
from pathlib import Path
import retrospective as r
import early_relation_pool as pool
import complete_prime_square_class_constructor as run
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_prime_square_class_constructor_verification_v1.json'


def check_bindings(obj):
    for name,sha in obj.get('bindings',{}).items():
        assert r.digest((r.ROOT/name).read_bytes())==sha,name


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari,matrix,vector
    from sage.env import SAGE_VERSION
    data=r.read(run.INPUT);out=r.read(run.OUTPUT);source=r.read(pool.INPUT)
    for obj in (data,out):check_bindings(obj)
    assert out['status']=='PASS'
    ids=[i for i,c in enumerate(source['columns']) if c['f']==1][:24]
    assert data['prime_ideals']==[{'retained_column':i,**source['columns'][i]} for i in ids]
    directions=[v for v in product(range(-2,3),repeat=3) if gcd(*v)==1 and next(x for x in v if x)>0]
    assert out['directions']==list(map(list,directions)) and len(directions)==49
    A=Algebra(data['cubic_ascending']);basis=list(map(A.elt,out['maximal_order_basis']))
    R=PolynomialRing(QQ,'z');f=R(data['cubic_ascending'])
    nf=pari.nfinit([pari(f),data['S_finite']]);assert str(nf.disc())==data['field_discriminant']
    assert basis==[A.elt([str(pari.lift(b).polcoef(i)) for i in range(3)]) for b in nf.nf_get_zk()]
    B=matrix(QQ,basis).transpose()
    assert B.det()**2*f.discriminant()==QQ(data['field_discriminant'])
    def trace(a):
        M=A.matrix(a);return sum(M[i][i] for i in range(3))
    T=matrix(QQ,[[trace(A.mul(x,y)) for y in basis] for x in basis])
    assert T==matrix(QQ,out['trace_gram']) and T.is_positive_definite()
    assert T.det()==QQ(data['field_discriminant'])
    assert len(out['rows'])==27
    verified=[]
    for i,row in enumerate(out['rows']):
        H=matrix(ZZ,row['ideal_hnf']);U=matrix(ZZ,row['lll_unimodular_matrix'])
        assert abs(U.det())==1 and abs(H.det())==ZZ(row['ideal_norm'])
        if i<24:
            c=data['prime_ideals'][i]
            assert row['id']==f"prime-column-{c['retained_column']}" and row['kind']=='candidate'
            assert row['prime']==c['p'] and row['prime_hnf']==c['hnf']
            P=next(P for P in pari.idealprimedec(nf,c['p']) if str(pari.idealhnf(nf,P))==c['hnf'])
            assert int(P[2])==c['e'] and int(P[3])==1
            I=pari.idealmul(nf,P,P)
            assert row['ideal_norm']==str(c['p']**2)
        else:
            j=i-24;g=data['generic_classes'][j]
            assert row['kind']=='positive_control' and row['generic_index']==j and row['id']==f'generic-{j}'
            beta=A.elt(g['beta_ascending']);assert A.norm(beta)==Q(g['norm'])
            I=pari.idealhnf(nf,pari(R(g['beta_ascending'])))
        assert H==matrix(ZZ,3,3,lambda a,b:ZZ(I[a,b]))
        # Ideal basis inclusion and equal absolute norm suffice to prove (alpha)=I.
        # This replay uses rational multiplication matrices, not PARI element norms.
        reduced=H*U;norms=[];hits=[]
        for v in directions:
            coords=reduced*vector(ZZ,v)
            alpha=tuple(sum(Q(int(coords[j]))*basis[j][k] for j in range(3)) for k in range(3))
            N=A.norm(alpha);assert N.denominator==1;norms.append(str(N))
            if abs(N)!=Q(row['ideal_norm']):continue
            sign=1 if N>0 else -1;alpha=tuple(sign*x for x in alpha)
            assert A.norm(alpha)==Q(row['ideal_norm'])
            hits.append({'direction':list(v),'sign_correction':sign,'alpha_ascending':list(map(str,alpha)),'norm':row['ideal_norm']})
        assert norms==row['tested_norms'] and hits==row['hits']
        verified.append({'id':row['id'],'kind':row['kind'],'norms_verified':len(norms),'generators':len(hits)})
    assert all(x['generators']==0 for x in verified[:24])
    assert all(x['generators']==1 for x in verified[24:])
    return {'schema':'rank-jump.prime-square-class-constructor-verification.v1','status':'PASS',
        'rows':verified,'rational_norm_checks':27*49,'candidate_generators':0,'positive_controls_recovered':3,
        'new_independent_classes_certified':0,'software':{'sage':SAGE_VERSION,'pari':str(pari('version()'))},
        'bindings':run.bindings([Path(__file__),Path(run.__file__),run.INPUT,run.OUTPUT,pool.INPUT,
            Path(__file__).with_name('verify_unpointed_governing_norm.py')]),
        'boundary':'No candidate principal generator found in the frozen box. This is not a nonprincipality certificate. Target ideals use PARI; all 1323 norms and equal-index principal-ideal witnesses use independent rational arithmetic.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args()
    result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['rational_norm_checks'],'norm checks;',result['positive_controls_recovered'],'controls;',result['candidate_generators'],'candidate generators')
