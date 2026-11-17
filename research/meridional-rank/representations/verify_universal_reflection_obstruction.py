#!/usr/bin/env python3
"""Replay integer Nullstellensatz identities for13a_373 reflection roots.

Standard library only: no Groebner software or symbolic algebra package.
The first-loop reduction and Coxeter transport are written arguments.
Here the normalized trace, reduced root formulas and all8 integer
polynomial certificates are reconstructed and checked exactly.
"""
from copy import deepcopy
from itertools import product
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
N=3


def constant(value):
    return {(0,)*N:value} if value else {}


def variable(i):
    exponent=[0]*N
    exponent[i]=1
    return {tuple(exponent):1}


def add(*polynomials):
    answer={}
    for polynomial in polynomials:
        for exponent,coefficient in polynomial.items():
            answer[exponent]=answer.get(exponent,0)+coefficient
            if not answer[exponent]:
                del answer[exponent]
    return answer


def scale(polynomial,value):
    return {exponent:coefficient*value for exponent,coefficient in polynomial.items()
            if coefficient*value}


def multiply(a,b):
    answer={}
    for ea,ca in a.items():
        for eb,cb in b.items():
            exponent=tuple(x+y for x,y in zip(ea,eb))
            answer[exponent]=answer.get(exponent,0)+ca*cb
            if not answer[exponent]:
                del answer[exponent]
    return answer


def dot(u,v,gram):
    return add(*(multiply(multiply(a,gram[i][j]),b)
                 for i,a in enumerate(u) for j,b in enumerate(v)))


def reflect(u,v,gram):
    factor=dot(u,v,gram)
    return [add(a,scale(multiply(factor,b),-1)) for a,b in zip(v,u)]


def unit(i,size):
    return [constant(int(i==j)) for j in range(size)]


def subtract_vectors(u,v,sign=1):
    return [add(a,scale(b,-sign)) for a,b in zip(u,v)]


def reduced_roots():
    c,b,d=[variable(i) for i in range(3)]
    gram=[[constant(2),c,scale(b,-1)],
          [c,constant(2),d],[scale(b,-1),d,constant(2)]]
    U,V,W=[unit(i,3) for i in range(3)]
    roots={0:U,4:[scale(a,-1) for a in U],6:V,7:V,8:W}
    roots[5]=reflect(W,roots[4],gram)
    roots[9]=reflect(roots[5],W,gram)
    roots[10]=reflect(V,roots[9],gram)
    roots[12]=reflect(roots[10],U,gram)
    roots[11]=reflect(U,roots[10],gram)
    assert all(dot(r,r,gram)==constant(2) for r in roots.values())
    return gram,roots


def equations(gram,roots,signs):
    crossing_triples=[(12,5,6),(11,7,8),(9,11,12)]
    return [f for (o,a,b),sign in zip(crossing_triples,signs)
            for f in subtract_vectors(reflect(roots[o],roots[a],gram),roots[b],sign)]


def read_polynomial(terms):
    answer={}
    for a,b,d,coefficient in terms:
        assert min(a,b,d)>=0 and isinstance(coefficient,int) and coefficient
        exponent=(a,b,d)
        assert exponent not in answer
        answer[exponent]=coefficient
    return answer


def verify_branch(gram,roots,branch):
    fs=equations(gram,roots,branch['signs'])
    assert len(fs)==len(branch['multipliers'])==9
    hs=[read_polynomial(terms) for terms in branch['multipliers']]
    identity=add(*(multiply(h,f) for h,f in zip(hs,fs)))
    assert identity==constant(15), 'integer Nullstellensatz identity fails'
    assert branch['constant']==15
    count=sum(len(h) for h in hs)
    degree=max(sum(e) for h in hs for e in h)
    assert count==branch['monomials']==285
    assert degree==branch['coefficient_degree']==9
    return {'signs':branch['signs'],'constant':15,
            'multiplier_monomials':count,'multiplier_degree':degree,
            'equation_degrees':[max(map(sum,f),default=-1) for f in fs]}


def verify_trace_binding(case):
    assert case['seed_arcs']==[1,2,6,8]
    expected_trace=[[2,6,7],[7,2,3],[3,1,0],[1,3,4],[8,4,5],
                    [5,8,9],[6,9,10],[10,0,12],[4,10,11]]
    assert case['seed_trace']==expected_trace
    components=[{i} for i in range(1,27)]
    for row in case['pd']:
        a=next(c for c in components if row[1] in c)
        b=next(c for c in components if row[3] in c)
        if a is not b:
            components.remove(a)
            components.remove(b)
            components.append(a|b)
    components=sorted(sorted(c) for c in components)
    arc={v:i for i,c in enumerate(components) for v in c}
    crossing={(arc[r[1]],)+tuple(sorted((arc[r[0]],arc[r[2]]))) for r in case['pd']}
    used=set()
    known=set(case['seed_arcs'])
    for o,a,b in expected_trace:
        assert o in known and a in known and b not in known
        assert (o,)+tuple(sorted((a,b))) in crossing
        known.add(b)
        used.add((o,)+tuple(sorted((a,b))))
    assert known==set(range(13))
    left=[(0,1,2),(12,5,6),(11,7,8),(9,11,12)]
    assert crossing-used=={(o,)+tuple(sorted((a,b))) for o,a,b in left}
    return expected_trace


def verify_four_coordinate_reduction(trace,roots3):
    global N
    N=4
    c,b,d,p=[variable(i) for i in range(4)]
    gram=[[constant(2),constant(1),c,p],
          [constant(1),constant(2),{},add(p,b)],
          [c,{},constant(2),d],
          [p,add(p,b),d,constant(2)]]
    roots={a:unit(i,4) for i,a in enumerate([1,2,6,8])}
    for o,a,b in trace:
        roots[b]=reflect(roots[o],roots[a],gram)
    assert all(dot(r,r,gram)==constant(2) for r in roots.values())
    assert reflect(roots[0],roots[1],gram)==roots[2]
    def extend(f):
        return {e+(0,):v for e,v in f.items()}
    for arc,root in roots3.items():
        x,y,z=map(extend,root)
        assert roots[arc]==[x,scale(x,-1),y,z], 'four-to-three root reduction fails'
    N=3


def main():
    inputs=json.loads((HERE/'seven_sheet_cover_inputs.json').read_text())
    case=next(c for c in inputs['cases'] if c['name']=='13a_373')
    trace=verify_trace_binding(case)
    certificates=json.loads((HERE/'universal_reflection_certificates.json').read_text())
    assert certificates['variables']==['c','b','d']
    branches=certificates['branches']
    assert [tuple(b['signs']) for b in branches]==list(product((1,-1),repeat=3))
    gram,roots=reduced_roots()
    results=[verify_branch(gram,roots,b) for b in branches]
    verify_four_coordinate_reduction(trace,roots)
    bad=deepcopy(branches[0])
    bad['multipliers'][0][0][3]+=1
    try:
        verify_branch(gram,roots,bad)
    except AssertionError as error:
        assert 'integer Nullstellensatz identity fails' in str(error)
    else:
        raise AssertionError('corrupt polynomial multiplier accepted')
    files=['verify_universal_reflection_obstruction.py',
           'universal_reflection_certificates.json','seven_sheet_cover_inputs.json',
           'UNIVERSAL_REFLECTION_OBSTRUCTION.tex']
    receipt={'status':'PASS','knot':'13a_373','sign_branches':results,
             'normalized_four_coordinate_trace_replayed':True,
             'all_reduced_root_norms':2,
             'corrupted_multiplier_rejected':True,
             'scope':'Exact integer polynomial identities constant15 for all8 reduced sign branches, independently reconstructed reduced reflection equations and normalized four-coordinate seed trace. The initial nondivisive first-loop reduction, root-sign normalization and transport to all Coxeter targets remain written arguments. Characteristic zero only; no determinant localization or Groebner implementation trusted by replay.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                              for p in files}}
    (HERE/'universal_reflection_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','knot':'13a_373','branches':8,
                      'certificate_constant':15,'multiplier_monomials_per_branch':285,
                      'corrupt_multiplier_rejected':True},indent=2))


if __name__=='__main__':
    main()
