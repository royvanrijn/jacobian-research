#!/usr/bin/env python3
"""Replay exact algebra and diagram transfers for sixteen reflection obstructions.

Standard library only: no Groebner software or symbolic algebra package.
The first-loop reduction and Coxeter transport are written arguments.
Here the normalized trace, reduced root formulas and all8 integer
polynomial certificates are reconstructed and checked exactly. The eight new
A6-diagram transfers also replay their raw local reductions and exact relabelings.
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


def verify_raw_first_loop():
    global N
    N=6
    a01,a02,a03,t,a13,a23=[variable(i) for i in range(6)]
    gram=[[constant(2),a01,a02,a03],
          [a01,constant(2),t,a13],
          [a02,t,constant(2),a23],
          [a03,a13,a23,constant(2)]]
    e0,e1,e2,e3=[unit(i,4) for i in range(4)]
    u7=reflect(e1,e2,gram)
    u3=reflect(u7,e1,gram)
    A=add(constant(1),scale(multiply(t,t),-1))
    assert u7==[{},scale(t,-1),constant(1),{}]
    assert u3==[{},A,t,{}]
    s=dot(u3,e0,gram)
    assert s==add(multiply(A,a01),multiply(t,a02))
    u0=reflect(u3,e0,gram)
    assert u0==subtract_vectors(e0,[multiply(s,a) for a in u3])
    H=add(constant(2),scale(multiply(s,s),-1))
    assert dot(u0,e0,gram)==H
    assert reflect(u0,e0,gram)==[
        add(constant(1),scale(H,-1)),
        multiply(multiply(H,s),A),multiply(multiply(H,s),t),{}]
    assert all(dot(r,r,gram)==constant(2) for r in [u7,u3,u0])
    # Choosing either sign of seed e1 leaves its reflection unchanged.
    # Recomputed trace changes u3's sign but not u7 or u0, so the first
    # closing sign epsilon can be normalized to+1 after proving t=0.
    for epsilon in (1,-1):
        seed1=[scale(a,epsilon) for a in e1]
        new7=reflect(seed1,e2,gram)
        new3=reflect(new7,seed1,gram)
        new0=reflect(new3,e0,gram)
        assert new7==u7
        assert new3==[scale(a,epsilon) for a in u3]
        assert new0==u0
    # The sibling seed changes preserve a four-root basis for every
    # value of the old pair Gram entry, including degenerate pair forms.
    new391=reflect(e3,e2,gram)
    assert new391==[{}, {}, constant(1), scale(a23,-1)]
    newW392=reflect(e2,e3,gram)
    newV392=reflect(newW392,e2,gram)
    assert newW392==[{}, {}, scale(a23,-1), constant(1)]
    assert newV392==[{}, {}, add(constant(1),scale(multiply(a23,a23),-1)),a23]
    assert add(newV392[2],multiply(a23,a23))==constant(1)
    N=3
    return {'variables':['a01','a02','a03','a12','a13','a23'],
            'raw_first_three_steps_verified':True,
            'first_closing_coefficients':['1-H','Hs(1-t^2)','Hst','0'],
            'both_seed_root_sign_normalizations_verified':True}


def verify_short_integer_contradiction():
    """A generated integer identity1 for all4 alpha,beta branches.

    The auxiliary symbols y,h,v,w,p,q are independent. This identity
    avoids every division, including by3 or5, and is separate from the
    optional identities15 in the reduced3-variable ring.
    """
    global N
    N=6
    y,h,v,w,p,q=[variable(i) for i in range(6)]
    one=constant(1)
    def cert_basis(i):
        return [one if i==j else {} for j in range(6)]
    def cert_add(*certificates):
        return [add(*(certificate[i] for certificate in certificates))
                for i in range(6)]
    def cert_times(poly,certificate):
        return [multiply(poly,a) for a in certificate]
    reports=[]
    for alpha,beta in product((1,-1),repeat=2):
        F=[add(h,constant(-1),multiply(y,y)),
           multiply(p,add(y,scale(w,-1))),
           add(one,multiply(p,v)),
           add(scale(multiply(p,h),-1),constant(-beta)),
           add(constant(-1),scale(multiply(q,add(one,scale(multiply(w,y),-1))),-1)),
           add(scale(multiply(multiply(q,w),v),-1),constant(-alpha))]
        def check(certificate,target):
            assert add(*(multiply(a,f) for a,f in zip(certificate,F)))==target
        D=add(w,scale(y,-1))
        Dc=cert_add(cert_times(v,cert_basis(1)),cert_times(D,cert_basis(2)))
        check(Dc,D)
        E=add(v,scale(h,-beta))
        Ec=cert_add(cert_times(scale(h,-beta),cert_basis(2)),
                    cert_times(scale(v,-beta),cert_basis(3)))
        check(Ec,E)
        A=add(constant(-1),scale(multiply(q,h),-1))
        Ac=cert_add(cert_basis(4),cert_times(scale(q,-1),cert_basis(0)),
                    cert_times(scale(multiply(q,y),-1),Dc))
        check(Ac,A)
        B=add(scale(y,beta),constant(-alpha))
        Bc=cert_add(cert_basis(5),cert_times(scale(y,-beta),Ac),
                    cert_times(multiply(q,v),Dc),cert_times(multiply(q,y),Ec))
        check(Bc,B)
        Hc=cert_add(cert_basis(0),cert_times(scale(add(scale(y,beta),constant(alpha)),-1),Bc))
        check(Hc,h)
        C=cert_add(cert_times(constant(-1),Ac),cert_times(scale(q,-1),Hc))
        check(C,one)
        reports.append({'alpha':alpha,'beta':beta,'constant':1,
                        'multiplier_monomials':sum(len(a) for a in C),
                        'multiplier_degree':max(sum(e) for a in C for e in a)})
    N=3
    return reports


def sibling_roots(name,gram):
    U,V,W=[unit(i,3) for i in range(3)]
    roots={0:U,4:[scale(a,-1) for a in U],6:V,7:V}
    if name=='13a_375':
        roots[9]=W
        plan=[(6,9,10),(10,7,8),(4,10,11),(11,4,5),(8,11,12)]
    elif name=='13a_391':
        roots[9]=W
        plan=[(9,6,5),(6,9,10),(4,10,11),(5,11,12),(11,7,8)]
    else:
        assert name=='13a_392'
        roots[12]=W
        plan=[(12,6,5),(5,12,11),(4,11,10),(6,10,9),(10,7,8)]
    for o,a,b in plan:
        roots[b]=reflect(roots[o],roots[a],gram)
    assert all(dot(r,r,gram)==constant(2) for r in roots.values())
    return roots,plan


def verify_siblings(inputs):
    global N
    N=3
    x,y,z=[variable(i) for i in range(3)]
    gram=[[constant(2),x,y],[x,constant(2),z],[y,z,constant(2)]]
    U,V,W=[unit(i,3) for i in range(3)]
    reports=[]
    original={'13a_375':[1,2,6,9],'13a_391':[1,2,5,9],
              '13a_392':[1,2,5,11]}
    for name in original:
        case=next(c for c in inputs['cases'] if c['name']==name)
        assert case['seed_arcs']==original[name]
        roots,plan=sibling_roots(name,gram)
        components=[{i} for i in range(1,27)]
        for row in case['pd']:
            a=next(c for c in components if row[1] in c)
            b=next(c for c in components if row[3] in c)
            if a is not b:
                components.remove(a);components.remove(b);components.append(a|b)
        components=sorted(sorted(c) for c in components)
        arc={v:i for i,c in enumerate(components) for v in c}
        crossings={(arc[r[1]],)+tuple(sorted((arc[r[0]],arc[r[2]]))) for r in case['pd']}
        prefix=[(2,6,7),(7,2,3),(3,1,0),(1,3,4)]
        new_seeds=[1,2,6,12 if name=='13a_392' else 9]
        known=set(new_seeds);used=set()
        for o,a,b in prefix+plan:
            assert o in known and a in known and b not in known
            edge=(o,)+tuple(sorted((a,b)))
            assert edge in crossings
            known.add(b);used.add(edge)
        assert known==set(range(13))
        leftovers=crossings-used
        expected={'13a_375':[(0,1,2),(5,8,9),(12,5,6),(9,12,0)],
                  '13a_391':[(0,1,2),(8,4,5),(12,8,9),(10,12,0)],
                  '13a_392':[(0,1,2),(8,4,5),(11,8,9),(9,12,0)]}[name]
        assert leftovers=={(o,)+tuple(sorted((a,b))) for o,a,b in expected}
        if name in ('13a_375','13a_391'):
            T=roots[10];A=roots[11]
            t=add(y,scale(multiply(x,z),-1))
            assert T==[{},scale(z,-1),constant(1)]
            assert A==[scale(t,-1),scale(z,-1),constant(1)]
            assert roots[12][0]==scale(t,-1)
            if name=='13a_375':
                assert roots[8]==reflect(T,V,gram)
                assert roots[5]==[
                    add(multiply(t,t),constant(-1)),multiply(t,z),scale(t,-1)]
            else:
                B=roots[5]
                r=dot(B,A,gram)
                assert B==[{},constant(1),scale(z,-1)]
                assert roots[12][1]==add(scale(z,-1),scale(r,-1))
                assert roots[8]==reflect(A,V,gram)
                assert reflect(A,roots[8],gram)==V
            # Last unused crossing forces t=-gamma. Substitute that
            # monic linear equation as y=xz-gamma, without dividing.
            for gamma in (1,-1):
                new_y=add(multiply(x,z),constant(-gamma))
                special=[[constant(2),x,new_y],[x,constant(2),z],
                         [new_y,z,constant(2)]]
                rr,_=sibling_roots(name,special)
                if name=='13a_375':
                    assert rr[5]==[scale(a,gamma) for a in rr[10]]
                    assert reflect(rr[5],rr[8],special)==V
                else:
                    target=[scale(a,gamma) for a in reflect(rr[10],U,special)]
                    difference=subtract_vectors(rr[12],target)
                    r=dot(rr[5],rr[11],special)
                    assert difference==[{},scale(r,-1),multiply(r,z)]
        else:
            h=add(constant(1),scale(multiply(z,z),-1))
            t=add(multiply(z,x),multiply(h,y))
            k=dot(V,roots[10],gram)
            assert roots[11]==[{},z,h]
            assert roots[10]==[scale(t,-1),z,h]
            assert roots[9]==[scale(t,-1),add(z,scale(k,-1)),h]
            assert roots[8]==[multiply(k,t),
                              add(constant(1),scale(multiply(k,z),-1)),
                              scale(multiply(k,h),-1)]
        reports.append({'name':name,'original_seed_arcs':original[name],
                        'normalized_seed_arcs':new_seeds,'reduced_trace':plan,
                        'unused_crossings':expected,'root_formula_checks':True})
    # The two seed replacements have determinant1, irrespective of the
    # original pair Gram determinant. Check both identities over Z[z].
    assert add(add(constant(1),scale(multiply(z,z),-1)),multiply(z,z))==constant(1)
    # Division-free auxiliary identity1 for392, using only last/middle
    # crossing coefficients and h=1-z^2.
    N=5
    z,h,t,k,p=[variable(i) for i in range(5)]
    one=constant(1)
    for beta,gamma in product((1,-1),repeat=2):
        F=[add(h,constant(-1),multiply(z,z)),
           add(multiply(p,t),constant(-gamma)),
           multiply(p,add(z,scale(k,-1))),
           add(one,scale(multiply(p,h),-1)),
           multiply(t,add(k,constant(beta)))]
        def cb(i):return [one if i==j else {} for j in range(5)]
        def ca(*cs):return [add(*(c[i] for c in cs)) for i in range(5)]
        def cm(poly,c):return [multiply(poly,a) for a in c]
        def check(c,target):assert add(*(multiply(a,f) for a,f in zip(c,F)))==target
        E=add(z,scale(k,-1));K=add(k,constant(beta))
        Ec=cm(constant(gamma),ca(cm(t,cb(2)),cm(scale(E,-1),cb(1))))
        Kc=cm(constant(gamma),ca(cm(p,cb(4)),cm(scale(K,-1),cb(1))))
        check(Ec,E);check(Kc,K)
        Zc=ca(Ec,Kc);check(Zc,add(z,constant(beta)))
        Hc=ca(cb(0),cm(scale(add(z,constant(-beta)),-1),Zc));check(Hc,h)
        C=ca(cb(3),cm(p,Hc));check(C,one)
    N=3
    return reports


def verify_remaining_alternating(inputs):
    global N
    N=6
    a01,a02,a03,a12,a13,a23=[variable(i) for i in range(6)]
    gram=[[constant(2),a01,a02,a03],
          [a01,constant(2),a12,a13],
          [a02,a12,constant(2),a23],
          [a03,a13,a23,constant(2)]]
    e0,e1,P,Q=[unit(i,4) for i in range(4)]
    u0=reflect(e1,e0,gram)
    u4=reflect(e0,e1,gram)
    u2=reflect(u0,e1,gram)
    assert u0==[constant(1),scale(a01,-1),{},{}]
    assert u4==[scale(a01,-1),constant(1),{},{}]
    assert u2==[a01,add(constant(1),scale(multiply(a01,a01),-1)),{},{}]
    u7=reflect(u2,P,gram)
    assert u7[2]==constant(1) and not u7[3]
    r=dot(u7,e0,gram)
    assert reflect(u7,e0,gram)[2]==scale(r,-1)
    # Raw basis changes, before imposing a01=1 or any orthogonality.
    root12=reflect(Q,u0,gram)
    Vsimple=reflect(root12,P,gram)
    assert Vsimple[2]==constant(1)
    W3266=reflect(P,Q,gram)
    root12_3266=reflect(W3266,u0,gram)
    V3266=reflect(root12_3266,P,gram)
    determinant=add(multiply(V3266[2],W3266[3]),
                    scale(multiply(V3266[3],W3266[2]),-1))
    assert determinant==constant(1)
    V3275=reflect(Q,P,gram)
    assert V3275[2]==constant(1)
    assert all(dot(v,v,gram)==constant(2) for v in [root12,Vsimple,W3266,
                    root12_3266,V3266,V3275])
    # After the alternate first-loop normalization, U,V,W have a free
    #3-variable Gram form despite the unused fourth parameter p.
    N=4
    x,y,z,p=[variable(i) for i in range(4)]
    normalized=[[constant(2),constant(1),{},p],
                [constant(1),constant(2),scale(x,-1),add(p,scale(y,-1))],
                [{},scale(x,-1),constant(2),z],
                [p,add(p,scale(y,-1)),z,constant(2)]]
    e0,e1,V,W=[unit(i,4) for i in range(4)]
    U=subtract_vectors(e0,e1)
    u0=reflect(e1,e0,normalized);u4=reflect(e0,e1,normalized)
    u2=reflect(u0,e1,normalized);u7=reflect(u2,V,normalized)
    assert u0==U and u4==[scale(a,-1) for a in U]
    assert u2==e0 and u7==V
    assert reflect(u7,e0,normalized)==u2
    assert dot(U,V,normalized)==x and dot(U,W,normalized)==y
    assert dot(V,W,normalized)==z
    N=3
    records={
      '13a_3264':{'original':[1,3,5,10],'new':[1,3,6,10],
        'old_third':10,'paired':'13a_373','canonical_third':8,
        'basis_operations':[(10,0,12),(12,5,6)]},
      '13a_3266':{'original':[1,3,5,8],'new':[1,3,6,9],
        'old_third':9,'paired':'13a_375','canonical_third':9,
        'basis_operations':[(5,8,9),(9,0,12),(12,5,6)]},
      '13a_3275':{'original':[1,3,5,9],'new':[1,3,6,9],
        'old_third':9,'paired':'13a_391','canonical_third':9,
        'basis_operations':[(9,5,6)]},
      '13a_3276':{'original':[1,3,5,9],'new':[1,3,6,9],
        'old_third':9,'paired':'13a_392','canonical_third':12,
        'basis_operations':[(9,0,12),(12,5,6)]}}
    plans={
      '13a_373':[(8,4,5),(5,8,9),(6,9,10),(10,0,12),(4,10,11)],
      '13a_375':[(6,9,10),(10,7,8),(4,10,11),(11,4,5),(8,11,12)],
      '13a_391':[(9,6,5),(6,9,10),(4,10,11),(5,11,12),(11,7,8)],
      '13a_392':[(12,6,5),(5,12,11),(4,11,10),(6,10,9),(10,7,8)]}
    by_name={c['name']:c for c in inputs['cases']}
    def crossings(case):
        components=[{i} for i in range(1,27)]
        for row in case['pd']:
            a=next(c for c in components if row[1] in c)
            b=next(c for c in components if row[3] in c)
            if a is not b:
                components.remove(a);components.remove(b);components.append(a|b)
        components=sorted(sorted(c) for c in components)
        arc={v:i for i,c in enumerate(components) for v in c}
        return [(arc[r[1]],arc[r[0]],arc[r[2]]) for r in case['pd']]
    def key(step):
        o,a,b=step
        return (o,)+tuple(sorted((a,b)))
    def fold(a):return 0 if a==4 else 6 if a==7 else a
    def reduced(case):
        return {key(tuple(map(fold,s))) for s in crossings(case)
                if not set(s)&{1,2,3}}
    reports=[]
    local=[(3,1,0),(1,3,4),(0,3,2),(2,6,7),(7,1,2)]
    for name,record in records.items():
        case=by_name[name]
        assert case['seed_arcs']==record['original']
        allcross=set(map(key,crossings(case)))
        assert all(key(step) in allcross for step in local+record['basis_operations'])
        core=reduced(case)
        assert len(core)==8 and core==reduced(by_name[record['paired']])
        known={0,6,record['canonical_third']}
        canonical_plan=plans[record['paired']]
        for step in canonical_plan:
            o,a,b=map(fold,step)
            assert o in known and a in known and b not in known
            assert key((o,a,b)) in core
            known.add(b)
        assert known=={0,5,6,8,9,10,11,12}
        assert record['old_third'] in known
        reports.append({'name':name,'original_seeds':record['original'],
          'changed_seeds':record['new'],'paired_reduced_system':record['paired'],
          'canonical_third_root':record['canonical_third'],
          'unit_determinant_basis_operations':record['basis_operations'],
          'folded_crossings':[list(s) for s in sorted(core)],
          'all_reduced_roots_propagated_from_three':True})
    return reports



def verify_a6_alternating_transfers(old_inputs,new_inputs):
    """Raw local algebra and literal relabeling to the proved reduced systems."""
    global N
    N=6
    a,b,c,d,e,f=[variable(i) for i in range(6)]
    gram=[[constant(2),a,b,d],[a,constant(2),c,e],
          [b,c,constant(2),f],[d,e,f,constant(2)]]
    A,B,C,D=[unit(i,4) for i in range(4)]
    # Both changes of the last two seed roots are triangular with diagonal 1.
    newP=reflect(D,C,gram)
    newQ=reflect(C,D,gram)
    assert newP==[{}, {},constant(1),scale(f,-1)]
    assert newQ==[{}, {},scale(f,-1),constant(1)]
    local_reports=[]
    for kind in ('I','II'):
        u2=reflect(A,B,gram) if kind=='I' else reflect(C,B,gram)
        u3=reflect(C,u2,gram) if kind=='I' else reflect(A,u2,gram)
        t=add(c,scale(multiply(a,b),-1)) if kind=='I' else add(a,scale(multiply(b,c),-1))
        H=dot(u3,A,gram)
        if kind=='I':
            assert u2==[scale(a,-1),constant(1),{},{}]
            assert u3==[scale(a,-1),constant(1),scale(t,-1),{}]
            assert H==add(scale(a,-1),scale(multiply(b,t),-1))
            closing=[add(constant(1),multiply(a,H)),scale(H,-1),multiply(t,H),{}]
        else:
            assert u2==[{},constant(1),scale(c,-1),{}]
            assert u3==[scale(t,-1),constant(1),scale(c,-1),{}]
            assert H==scale(t,-1)
            closing=[add(constant(1),multiply(t,H)),scale(H,-1),multiply(c,H),{}]
        assert reflect(u3,A,gram)==closing
        assert all(dot(u,u,gram)==constant(2) for u in (u2,u3))
        for epsilon in (1,-1):
            F0=closing[0];F1=add(closing[1],constant(-epsilon));F2=closing[2]
            if kind=='I':
                assert add(a,constant(-epsilon))==scale(add(F0,multiply(a,F1)),-epsilon)
                assert t==scale(add(F2,multiply(t,F1)),-epsilon)
            else:
                assert add(t,constant(-epsilon))==scale(add(F0,multiply(t,F1)),-epsilon)
                assert c==scale(add(F2,multiply(c,F1)),-epsilon)
            # Independently substitute the forced Gram restrictions, for both signs.
            aa=constant(epsilon)
            cc=scale(b,epsilon) if kind=='I' else {}
            reduced_gram=[[constant(2),aa,b,d],[aa,constant(2),cc,e],
                          [b,cc,constant(2),f],[d,e,f,constant(2)]]
            v2=reflect(A,B,reduced_gram) if kind=='I' else reflect(C,B,reduced_gram)
            v3=reflect(C,v2,reduced_gram) if kind=='I' else reflect(A,v2,reduced_gram)
            assert reflect(v3,A,reduced_gram)==[scale(x,epsilon) for x in B]
            assert reflect(B,v3,reduced_gram)==[scale(x,-epsilon) for x in A]
            assert reflect(v2,C,reduced_gram)==C
        local_reports.append({'type':kind,'raw_six_variable_coefficients_verified':True,
                              'division_free_branch_reductions':2,
                              'both_normalized_sign_branches_verified':True})
    N=3
    records={
      '13a_304':{'kind':'I','original':[0,1,5,9],'new':[0,1,6,9],
        'C':6,'paired_arc':7,'independent_reduced':[0,6,9],
        'basis_operations':[(9,5,6)],'target':'13a_375',
        'relabel':{5:10,6:11,12:8,11:6,0:9,8:12,9:0,10:5}},
      '13a_362':{'kind':'I','original':[0,1,5,11],'new':[0,1,6,11],
        'C':6,'paired_arc':7,'independent_reduced':[0,6,11],
        'basis_operations':[(11,5,6)],'target':'13a_375',
        'relabel':{0:10,5:11,11:8,10:6,12:9,6:12,8:0,9:5}},
      '13a_623':{'kind':'I','original':[0,1,5,8],'new':[0,1,5,8],
        'C':8,'paired_arc':9,'independent_reduced':[0,5,8],
        'basis_operations':[],'target':'13a_373',
        'relabel':{11:10,10:11,7:0,8:12,0:8,5:6,12:9,6:5}},
      '13a_629':{'kind':'I','original':[0,1,5,7],'new':[0,1,5,8],
        'C':8,'paired_arc':9,'independent_reduced':[0,5,8],
        'basis_operations':[(5,7,8)],'target':'13a_373',
        'relabel':{8:10,7:11,5:0,6:12,11:8,12:6,10:9,0:5}}}
    for name,base in [('13a_3232','13a_304'),('13a_3259','13a_362'),
                      ('13a_3374','13a_623'),('13a_3375','13a_629')]:
        record=deepcopy(records[base]);record['kind']='II';records[name]=record
    def crossings(case):
        components=[{i} for i in range(1,27)]
        for row in case['pd']:
            aa=next(c for c in components if row[1] in c)
            bb=next(c for c in components if row[3] in c)
            if aa is not bb:
                components.remove(aa);components.remove(bb);components.append(aa|bb)
        components=sorted(sorted(c) for c in components)
        arc={v:i for i,c in enumerate(components) for v in c}
        assert len(components)==13
        return [(arc[r[1]],arc[r[0]],arc[r[2]]) for r in case['pd']]
    def key(step):
        o,aa,bb=step
        return (o,)+tuple(sorted((aa,bb)))
    def reduced(case,C,paired):
        def fold(i):return 0 if i==4 else C if i==paired else i
        return {key(tuple(map(fold,step))) for step in crossings(case)
                if not set(step)&{1,2,3}}
    old_by_name={c['name']:c for c in old_inputs['cases']}
    new_by_name={c['name']:c for c in new_inputs['cases']}
    reports=[]
    for name,record in records.items():
        case=new_by_name[name];C=record['C'];paired=record['paired_arc']
        assert case['seed_arcs']==record['original']
        crossing=set(map(key,crossings(case)))
        if record['kind']=='I':local=[(0,1,2),(C,2,3),(3,0,1),(1,3,4),(2,C,paired)]
        else:local=[(C,1,2),(0,2,3),(3,0,1),(1,3,4),(2,C,paired)]
        assert {s for s in crossing if set(s)&{1,2,3}}==set(map(key,local))
        assert all(key(s) in crossing for s in record['basis_operations'])
        # Replay the original four-seed trace independently of the input's claims.
        known=set(record['original'])
        for o,aa,bb in case['seed_trace']:
            assert o in known and aa in known and bb not in known
            assert key((o,aa,bb)) in crossing
            known.add(bb)
        assert known==set(range(13))
        core=reduced(case,C,paired)
        assert len(core)==8
        nodes={i for step in core for i in step}
        relabel=record['relabel']
        assert set(relabel)==nodes and len(set(relabel.values()))==8
        assert set(record['independent_reduced'])<=nodes
        image={key(tuple(relabel[i] for i in step)) for step in core}
        assert image==reduced(old_by_name[record['target']],6,7)
        reports.append({'name':name,'local_type':record['kind'],
            'original_seeds':record['original'],'changed_seeds':record['new'],
            'independent_reduced_roots':record['independent_reduced'],
            'unit_determinant_basis_operations':record['basis_operations'],
            'folded_pair':[C,paired],'paired_reduced_system':record['target'],
            'relabel_new_to_old':relabel,'folded_crossings':[list(s) for s in sorted(core)]})
    return {'raw_local_algebra':local_reports,'transfers':reports}


def singular_form_control():
    # A singular form can still support4 independent reflection directions.
    # This ensures the proof is about direction span, not invertible Gram.
    form=[[2]*4 for _ in range(4)]
    identity=[[int(i==j) for j in range(4)] for i in range(4)]
    def mm(a,b):
        return [[sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    for root in range(4):
        reflection=[row[:] for row in identity]
        reflection[root]=[x-2 for x in reflection[root]]
        assert mm(reflection,reflection)==identity
        displacement=[reflection[i][root]-identity[i][root] for i in range(4)]
        assert displacement==[-2*int(i==root) for i in range(4)]
    assert all(row==form[0] for row in form)
    return {'Gram_rank':1,'independent_root_directions':4,
            'involutions_verified':4,'displacements':'-2e_i'}


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
    raw=verify_raw_first_loop()
    short=verify_short_integer_contradiction()
    singular=singular_form_control()
    siblings=verify_siblings(inputs)
    alternating=verify_remaining_alternating(inputs)
    a6_inputs=json.loads((HERE/'a6_affine_inputs.json').read_text())
    a6_alternating=verify_a6_alternating_transfers(inputs,a6_inputs)
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
           'UNIVERSAL_REFLECTION_OBSTRUCTION.tex','a6_affine_inputs.json']
    receipt={'status':'PASS','knot':'13a_373','sibling_obstructions':siblings,'remaining_alternating_transfers':alternating,'a6_alternating_transfers':a6_alternating,'sign_branches':results,'division_free_short_proof_identities':short,'singular_Gram_positive_control':singular,
             'raw_six_variable_first_loop':raw,'normalized_four_coordinate_trace_replayed':True,
             'all_reduced_root_norms':2,
             'corrupted_multiplier_rejected':True,
             'scope':'Exact integer polynomial identities constant15 for all8 reduced sign branches, independently reconstructed reduced reflection equations and normalized four-coordinate seed trace. The generated auxiliary-symbol integer identities1 verify the division-free short contradiction in every field. The initial first-loop reduction, root-sign normalization and transport to all Coxeter targets remain written arguments; reflection theorem characteristic not2 and Coxeter corollary characteristic zero. Optional reduced identities15 are supporting evidence only, with no role at characteristics3,5. The eight A6-diagram transfers additionally check both raw local reductions, their division-free sign branches, determinant-one seed changes, original traces and exact reduced-system relabelings. No determinant localization or Groebner implementation trusted by replay.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                              for p in files}}
    (HERE/'universal_reflection_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','knot':'13a_373','branches':8,
                      'short_proof_constant':1,'short_proof_sign_branches':4,'optional_certificate_constant':15,
                      'corrupt_multiplier_rejected':True,'sibling_obstructions':len(siblings),'remaining_alternating_transfers':len(alternating),'a6_alternating_transfers':len(a6_alternating['transfers'])},indent=2))


if __name__=='__main__':
    main()
