#!/usr/bin/env python3
"""Independent multiplication-matrix replay of the bounded bad-prime test."""
import argparse
from fractions import Fraction as Q
from itertools import product
from math import gcd
from pathlib import Path
import retrospective as r
import bad_prime_principalization as run
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_bad_prime_principalization_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari,matrix,vector
    out=r.read(run.OUTPUT);inp=r.read(run.INPUT)
    for obj in (out,inp):
        for path,sha in obj['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    ds=[v for v in product(range(-2,3),repeat=3) if gcd(*v)==1 and next(x for x in v if x)>0]
    assert len(ds)==49
    results=[]
    for source,data in zip(inp['cases'],out['cases']):
        assert source['token']==data['token'] and data['status']=='PASS'
        A=Algebra(source['cubic']);R=PolynomialRing(QQ,'z');f=R(source['cubic'])
        nf=pari.nfinit([pari(f),[pari(R(b)) for b in source['basis']]])
        assert str(nf.disc())==source['field_discriminant']
        basis=[A.elt(b) for b in data['basis']]
        assert basis==[A.elt([str(pari.lift(b).polcoef(i)) for i in range(3)]) for b in nf.nf_get_zk()]
        B=matrix(QQ,basis).transpose();assert B.det()**2*f.discriminant()==QQ(source['field_discriminant'])
        trace=lambda a:sum(A.matrix(a)[i][i] for i in range(3))
        T=matrix(QQ,[[trace(A.mul(a,b)) for b in basis] for a in basis])
        assert T==matrix(QQ,data['trace_gram'])
        prime_rows=[(p,j,P) for p in source['S_finite'] for j,P in enumerate(pari.idealprimedec(nf,p))]
        candidates=[z for z in data['rows'] if z['kind']=='candidate']
        assert len(candidates)==len(prime_rows)
        for row,(p,j,P) in zip(candidates,prime_rows):
            assert (row['p'],row['prime_index'],row['e'],row['f'])==(p,j,int(P[2]),int(P[3]))
            assert matrix(QQ,row['target_hnf'])==matrix(QQ,pari.idealpow(nf,P,2))
        norms_checked=0;candidate_hits=0;control_hits=0
        for row in data['rows']:
            H=matrix(QQ,row['reduced_hnf']);I=matrix(QQ,row['target_hnf']);U=matrix(ZZ,row['lll_matrix'])
            assert abs(U.det())==1 and abs(H.det())==QQ(row['reduced_norm']) and abs(I.det())==QQ(row['target_norm'])
            multiplier=A.elt(row['multiplier'])
            # Multiplication by the retained element maps the reduced ideal onto the target.
            M=matrix(QQ,A.matrix(multiplier));change=I.inverse()*B.inverse()*M*B*H
            assert all(x.denominator()==1 for x in change.list()) and abs(change.det())==1
            if row['kind']=='control':
                beta=A.elt(source['generic_controls'][row['generic_index']])
                C=I.inverse()*B.inverse()*matrix(QQ,A.matrix(beta))*B
                assert all(x.denominator()==1 for x in C.list()) and abs(C.det())==1
            norms=[];hits=[]
            for v in ds:
                c=B*H*U*vector(ZZ,v);a=tuple(Q(str(x)) for x in c);n=A.norm(a);norms.append(str(n))
                if abs(n)!=Q(row['reduced_norm']):continue
                alpha=A.mul(a,multiplier);N=A.norm(alpha);sign=1 if N>0 else -1;alpha=tuple(sign*x for x in alpha)
                assert A.norm(alpha)==Q(row['target_norm'])
                C=I.inverse()*B.inverse()*matrix(QQ,A.matrix(alpha))*B
                assert all(x.denominator()==1 for x in C.list()) and abs(C.det())==1
                hits.append({'direction':list(v),'sign':sign,'alpha':list(map(str,alpha)),'norm':row['target_norm']})
            assert norms==row['tested_norms'] and hits==row['hits'];norms_checked+=len(norms)
            if row['kind']=='candidate':candidate_hits+=len(hits)
            else:control_hits+=len(hits)
        assert candidate_hits==data['candidate_generators']==0 and control_hits==data['control_generators']==3
        # Equality of retained reductions would be an exact pair principalization witness.
        reductions=[row['reduced_hnf'] for row in candidates]
        collisions=[[i,j] for i in range(len(reductions)) for j in range(i) if reductions[i]==reductions[j]]
        assert not collisions
        results.append({'token':data['token'],'prime_square_targets':len(candidates),'norms_verified':norms_checked,
            'candidate_generators':candidate_hits,'controls_recovered':control_hits,
            'identical_reduced_ideal_pairs':collisions,'principality_of_missed_targets':'UNKNOWN',
            'new_strict_classes_constructed':0})
    return {'schema':'rank-jump.bad-prime-principalization-verification.v1','status':'PASS','cases':results,
        'bindings':run.bind([Path(__file__),run.INPUT,run.OUTPUT,Path(__file__).with_name('verify_unpointed_governing_norm.py')]),
        'boundary':'Certifies the finite generator misses and six exact positive controls, not nonprincipality, class rank or a high/low incidence comparison.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','check']);args=parser.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',sum(c['norms_verified'] for c in result['cases']),'norms; six controls; no candidate generators')
