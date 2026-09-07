#!/usr/bin/env python3
"""Verify a dual incidence obstruction without using the worker's elimination."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import retrospective as r
import relation_root_class as root
import relation_root_affine as run
import complete_relation_root_class as completed
import early_relation_pool as pool
from verify_unpointed_governing_norm import Algebra

OUTPUT=r.OUT/'rank_jump_relation_root_affine_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari,matrix,vector
    d=r.read(pool.INPUT);out=r.read(run.OUTPUT);inp=r.read(root.INPUT);base=r.read(completed.OUTPUT)
    for obj in (out,inp,base):
        for name,sha in obj['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    assert out['status']=='PASS' and out['affine_status']=='INCONSISTENT'
    ids={w['column'] for w in out['dual_witness']};assert len(ids)==len(out['dual_witness'])
    S=set(d['S_finite']);byprime={}
    for i,c in enumerate(d['columns']):byprime.setdefault(c['p'],[]).append(i)
    rational_weights={p:sum(d['columns'][i]['e'] for i in cols if i in ids)%2 for p,cols in byprime.items()}
    # Direct sparse valuation formula; no call to the worker's dictionary(),
    # rank routine, elimination or dual construction.
    pairings=[]
    for row in d['relations']:
        z=0
        for i,v in row['ideal_factorization']:
            c=d['columns'][i]
            if c['p'] in S:continue
            z+=v*(int(i in ids)+c['f']*rational_weights[c['p']])
        pairings.append(z%2)
    assert pairings==out['dual_dictionary_values']==[0]*4134
    A=Algebra(d['cubic_ascending']);B=matrix(QQ,list(map(A.elt,base['maximal_order_basis']))).transpose();Binv=B.inverse()
    R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending']);nf=pari.nfinit([pari(f),d['S_finite']])
    assert B.det()**2*f.discriminant()==QQ(d['field_discriminant'])
    beta=A.elt(base['projection_ascending']);gamma=list(map(A.elt,inp['generic_classes_ascending']))
    elements=[beta]+gamma;coords=[Binv*vector(QQ,g) for g in elements]
    decompositions={};computed=[];generic_pairings=[0]*16;value=0;checks=0
    def mul_ideal(H,J):
        products=[]
        for i in range(3):
            for j in range(3):
                a=A.elt(map(str,B*H.column(i)));b=A.elt(map(str,B*J.column(j)))
                products.append(Binv*vector(QQ,A.mul(a,b)))
        return matrix(ZZ,products).hermite_form(include_zero_rows=False).transpose()
    for witness in out['dual_witness']:
        i=witness['column'];c=d['columns'][i];p=c['p']
        assert {k:witness[k] for k in c}==c and p not in S
        if p not in decompositions:
            assert ZZ(p).is_prime(proof=True) and f.discriminant()%p
            decompositions[p]={str(pari.idealhnf(nf,P)):P for P in pari.idealprimedec(nf,p)}
        P=decompositions[p][c['hnf']];assert int(P[2])==c['e']==1 and int(P[3])==c['f']
        H=matrix(ZZ,pari.idealhnf(nf,P));assert abs(H.det())==p**c['f']
        powers=[matrix(ZZ,3,3,1),H];inverses={0:powers[0],1:H.inverse()}
        def inverse(k):
            while len(powers)<=k:powers.append(mul_ideal(powers[-1],H))
            if k not in inverses:inverses[k]=powers[k].inverse()
            return inverses[k]
        valuations=[witness['projection_valuation']]+witness['generic_valuations']
        for v,x in zip(valuations,coords):
            assert v>=0
            assert all(QQ(y).denominator()%p for y in inverse(v)*x)
            assert any(QQ(y).denominator()%p==0 for y in inverse(v+1)*x)
            checks+=1
        value^=valuations[0]%2
        for j,v in enumerate(valuations[1:]):assert v%2==0;generic_pairings[j]^=v%2
        computed.append({'column':i,'p':p,'residue_degree':c['f'],'projection_valuation':valuations[0]})
    assert value==out['dual_projection_value']==1 and generic_pairings==[0]*16
    return {'schema':'rank-jump.relation-root-affine-verification.v1','status':'PASS',
        'dual_prime_ideal_support':len(ids),'dual_rational_prime_support':len(decompositions),
        'dictionary_pairings_verified':len(pairings),'all_dictionary_pairings_zero':True,
        'generic_pairings':generic_pairings,'root_projection_pairing':value,
        'independent_lattice_valuation_checks':checks,'dual_projection_valuations':computed,
        'conclusion':'The entire affine coset pi(w) times G times the span of all 4134 retained norm projections is disjoint from Selmer. The enlarged linear span adds no strict class.',
        'method':'Direct sparse norm-valuation dual formula; independently multiplied prime ideal lattices and rational p-integrality. No worker elimination or PARI idealval/idealpow is used.',
        'bindings':root.bindings([Path(__file__),run.OUTPUT,pool.INPUT,pool.OUTPUT,root.INPUT,completed.OUTPUT,
            Path(__file__).with_name('verify_unpointed_governing_norm.py')])}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS',result['dictionary_pairings_verified'],'dictionary pairings;',result['independent_lattice_valuation_checks'],'lattice valuation checks; root pairing',result['root_projection_pairing'])
