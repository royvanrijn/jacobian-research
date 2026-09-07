#!/usr/bin/env python3
"""Exact node factor and residual branch geometry for the MW16 triple block."""
import argparse
from collections import Counter
from itertools import combinations
from pathlib import Path
import retrospective as r
import additive_branch_geometry as first
import complete_additive_branch_geometry as prior

PROTOCOL=Path(__file__).with_name('ADDITIVE_NODE_BRANCH_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_additive_node_branch_v1.json'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,first.INPUT,prior.OUTPUT,Path(prior.__file__),Path(first.__file__),Path(r.__file__))}


def compute():
    from sage.all import QQ,GF,PolynomialRing
    source=next(x for x in r.read(first.INPUT)['families'] if x['family']=='a1-fibration-01')
    old=next(x for x in r.read(prior.OUTPUT)['rows'] if x['family']=='a1-fibration-01')
    R=PolynomialRing(QQ,'t');A=R(source['A']);B=R(source['B']);D=-4*A**3-27*B**2
    L=D.gcd(D.derivative()).monic();assert L.degree()==1
    tau=-L[0];a,b=A(tau),B(tau);assert a!=0
    e=-3*b/(2*a);assert a==-3*e*e and b==2*e**3
    sections=prior.functions(source,R);through=[];values=[]
    for i,(x,y) in enumerate(sections):
        assert x.denominator()(tau) and y.denominator()(tau)
        xx=x.numerator()(tau)/x.denominator()(tau);yy=y.numerator()(tau)/y.denominator()(tau)
        values.append({'index':i,'x':str(xx),'y':str(yy)})
        if xx==e:assert yy==0;through.append(i)
    assert len(through)==10
    expected={tuple(x) for x in combinations(through,3)};records=[];polys=[]
    for row in old['norms']:
        N=R(row['norm_ascending']);q=N;v=0
        while q%L==0:q//=L;v+=1
        assert v==(4 if tuple(row['generic_indices']) in expected else 0)
        branch=first.primitive(q);scalar=q.leading_coefficient()/branch.leading_coefficient()
        assert N==scalar*L**v*branch
        polys.append(branch);records.append({'generic_indices':row['generic_indices'],'norm_degree':row['degree'],
            'node_multiplicity':v,'branch_degree':int(branch.degree()),'branch_ascending':list(map(str,branch.list())),
            'norm_decomposition_scalar':str(scalar),'squarefree_witness_prime':None,'discriminant_coprime_witness_prime':None})
    primes=r.read(PROTOCOL)['limits']['primes'];reductions={}
    for p in primes:
        P=PolynomialRing(GF(p),'t');dp=P(D);vals=[P(q) for q in polys];reductions[p]=vals
        for row,q in zip(records,vals):
            if q.degree()!=row['branch_degree']:continue
            if row['squarefree_witness_prime'] is None and q.gcd(q.derivative())==1:row['squarefree_witness_prime']=p
            if row['discriminant_coprime_witness_prime'] is None and dp.degree()==D.degree() and q.gcd(dp)==1:row['discriminant_coprime_witness_prime']=p
    exceptions=[];unresolved=[];count=0
    for i in range(len(polys)):
        for j in range(i):
            count+=1;w=None
            for p in primes:
                a,b=reductions[p][i],reductions[p][j]
                if a.degree()==records[i]['branch_degree'] and b.degree()==records[j]['branch_degree'] and a.gcd(b)==1:w=p;break
            if w is None:unresolved.append([j,i])
            elif w!=primes[0]:exceptions.append({'pair':[j,i],'prime':w})
    complete=not unresolved and all(x['squarefree_witness_prime'] and x['discriminant_coprime_witness_prime'] for x in records)
    degrees=sorted(x['branch_degree'] for x in records);assert all(d%2==0 for d in degrees)
    genera=[{'classes':k,'genus_range':[1+2**k*(sum(ds)-4)//4 for ds in (degrees[:k],degrees[-k:])]} for k in range(1,5)] if complete else []
    return {'schema':'rank-jump.additive-node-branch.v1','status':'PASS' if complete else 'PARTIAL','bindings':bindings(),
        'family':'a1-fibration-01','node_parameter':str(tau),'node_x':str(e),'node_factor_ascending':list(map(str,L.list())),
        'section_values_at_node':values,'sections_through_node':through,'node_triple_count':len(expected),'node_norm_multiplicity':4,
        'rows':records,'pair_count':count,'default_pair_witness_prime':primes[0],
        'pair_witness_exceptions':exceptions,'unresolved_pairs':unresolved,
        'branch_degree_distribution':{str(k):v for k,v in sorted(Counter(degrees).items())},
        'simultaneous_norm_cover_genus_ranges_first_four':genera,
        'boundary':'A common fourth-power factor at a fixed singular fibre; it creates no shared odd norm branch. Rational points on residual norm covers remain only necessary incidence conditions.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],result['node_parameter'],result['sections_through_node'],result['branch_degree_distribution'],flush=True)
