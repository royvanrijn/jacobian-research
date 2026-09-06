#!/usr/bin/env python3
"""Pure rational group-law and finite-character replay of quartet relations."""
import argparse
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import retrospective as r
from branch_blocks import qrank

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_paired_quartet_relations_inputs_v1.json'
SOURCE=r.OUT/'rank_jump_paired_quartet_relations_v1.json'
OUTPUT=r.OUT/'rank_jump_paired_quartet_relations_verification_v1.json'


def add(A,P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y==-v:return None
    m=(v-y)/(u-x) if x!=u else (3*x*x+A)/(2*y)
    xx=m*m-x-u;return (xx,m*(x-xx)-y)


def mul(A,n,P):
    if n<0:return mul(A,-n,(P[0],-P[1]))
    out=None
    while n:
        if n&1:out=add(A,out,P)
        P=add(A,P,P);n>>=1
    return out


def verify():
    inp=r.read(INPUT);source=r.read(SOURCE)
    for data in (inp,source):
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    rows=[]
    for c,old in zip(inp['cases'],source['rows'],strict=True):
        assert old['execution']['status']=='COMPLETE';result=old['result'];assert result['status']=='PASS'
        assert result['id']==c['id'];model,basis=r.short(c['model'],c['basis'])
        assert all(F(x)==0 for x in model[:3])
        A=F(model[3]);points=[tuple(map(F,P)) for P in basis];lifts=[tuple(map(F,x['point'])) for x in c['lifts']]
        r.short(model,lifts)
        primes=[x['prime'] for x in c['rank_certificate']['signatures']]
        blocks=[(p,r.roots_at(model[3],model[4],p)) for p in primes]
        signatures=[r.point_signature(model,P,blocks) for P in points]
        assert r.rank(signatures)==len(points)==result['basis_rank'] and r.rank(signatures[:17])==17
        coefficients=[list(map(F,row)) for row in result['coordinates_in_witness_basis']]
        for P,coords,relation in zip(lifts,coefficients,result['exact_relations'],strict=True):
            den=relation['lift_multiplier'];ints=relation['basis_coefficients']
            assert [den*x for x in coords]==ints
            value=None
            for n,Q in zip(ints,points,strict=True):value=add(A,value,mul(A,n,Q))
            assert value==mul(A,den,P)
        Q=[v[17:] for v in coefficients];assert qrank(Q)==result['exact_quotient_rank']==3
        for v,g in zip(result['kernel_integer_vectors'],result['kernel_generic_coordinates'],strict=True):
            assert all(sum(v[i]*Q[i][j] for i in range(4))==0 for j in range(len(Q[0])))
            assert list(map(F,g))==[sum(v[i]*coefficients[i][j] for i in range(4)) for j in range(17)]
        assert qrank(result['kernel_integer_vectors'])==1
        pairs=[{'indices':list(ij),'quotient_rank':qrank([Q[i] for i in ij])} for ij in combinations(range(4),2)]
        triples=[{'indices':list(ij),'quotient_rank':qrank([Q[i] for i in ij])} for ij in combinations(range(4),3)]
        rows.append({'id':c['id'],'basis_rank':len(points),'exact_quotient_rank':3,'pair_ranks':pairs,'triple_ranks':triples,
                     'kernel_integer_vectors':result['kernel_integer_vectors'],'kernel_generic_coordinates':result['kernel_generic_coordinates']})
    return {'schema':'rank-jump.paired-quartet-relations-verification.v1','status':'PASS','rows':rows,
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,SOURCE,Path(__file__),HERE/'retrospective.py',HERE/'branch_blocks.py')},
            'boundary':'Exact rational group addition and independent finite Kummer signatures; no numerical-height assumption in this replay.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);mode=p.parse_args().mode
    result=verify()
    if mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print('PASS both exact quartet quotient ranks3; all eight coordinate relations verified')
