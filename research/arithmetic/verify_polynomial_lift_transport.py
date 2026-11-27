#!/usr/bin/env python3
"""Exact controls for literal return tuples and coupled polynomial lift transport.

Universal statements are written proofs. Finite controls retain their precise
quotient/inner degrees, full domains, braid words and negative inputs.
"""
import verify_polynomial_passport_components as c
from functools import lru_cache
from collections import deque
from itertools import permutations, combinations, product
from pathlib import Path
from math import comb
import argparse
import hashlib
import platform
import json

if not __debug__:
    raise RuntimeError("lift transport checker requires active assertions")

HERE = Path(__file__).resolve().parent
m=c.m

@lru_cache(None)
def layout(q):
    d=len(q[0]);r=len(q)
    assert all(len(g)==d and sorted(g)==list(range(d)) for g in q)
    assert m.product(q)==m.perm(d,[tuple(range(d))])
    assert sum(c.index(g) for g in q)==d-1
    slots=tuple((i,C) for i,g in enumerate(q) for C in c.cycles(g))
    adjacent=[[] for _ in range(d)]
    for j,(_,C) in enumerate(slots):
        for a in C:adjacent[a].append(j)
    seen={0};done=set();queue=deque([0]);steps=[];breaks=[None]*len(slots)
    while queue:
        a=queue.popleft()
        for j in adjacent[a]:
            if j in done:continue
            done.add(j);i,C=slots[j];prev=a;x=q[i][a]
            while x!=a:
                assert x not in seen
                seen.add(x);queue.append(x);steps.append((i,prev,x));prev=x;x=q[i][x]
            breaks[j]=(i,prev)
    assert len(seen)==d and len(done)==len(slots) and len(steps)==d-1
    by_edge={e:j for j,e in enumerate(breaks)};walk=[];edges=set();a=0
    for _ in range(d):
        for i in range(r-1,-1,-1):
            edge=(i,a);assert edge not in edges;edges.add(edge)
            if edge in by_edge:walk.append(by_edge[edge])
            a=q[i][a]
    assert a==0 and len(edges)==r*d and len(walk)==len(slots)==(r-1)*d+1
    assert len(set(walk))==len(slots)
    return slots,tuple(steps),tuple(breaks),tuple(walk[::-1])

def to_wreath(state,d):
    n=len(state[0]);assert n%d==0;b=n//d;answer=[]
    for g in state:
        assert len(g)==n and sorted(g)==list(range(n))
        q=tuple(g[a]%d for a in range(d));assert sorted(q)==list(range(d));maps=[]
        for a in range(d):
            assert all(g[a+d*t]%d==q[a] for t in range(b))
            A=tuple((g[a+d*t]-q[a])//d for t in range(b))
            assert sorted(A)==list(range(b))
            maps.append(A)
        answer.append((q,tuple(maps)))
    return tuple(answer)

def to_perms(state):
    d=len(state[0][0]);b=len(state[0][1][0]);answer=[]
    for q,maps in state:
        answer.append(tuple(q[x%d]+d*maps[x%d][x//d] for x in range(b*d)))
    return tuple(answer)

def wmul(A,B):
    q,x=A;s,y=B
    return m.compose(q,s),tuple(m.compose(x[s[a]],y[a]) for a in range(len(q)))

def winv(A):
    q,x=A;p=m.inverse(q)
    return p,tuple(m.inverse(x[p[a]]) for a in range(len(q)))

def wproduct(state):
    d=len(state[0][0]);b=len(state[0][1][0]);I=tuple(range(b))
    result=(tuple(range(d)),(I,)*d)
    for A in state:result=wmul(result,A)
    return result

def wbraid(state,move):
    i=abs(move)-1;A,B=state[i:i+2]
    pair=(B,wmul(wmul(winv(B),A),B)) if move>0 else (wmul(wmul(A,B),winv(A)),A)
    return state[:i]+pair+state[i+2:]

def extract(state):
    q=tuple(a[0] for a in state);slots,steps,breaks,order=layout(q)
    d=len(q[0]);b=len(state[0][1][0]);I=tuple(range(b));charts=[None]*d;charts[0]=I
    assert wproduct(state)==(m.perm(d,[tuple(range(d))]),
                            tuple(I if a<d-1 else m.perm(b,[tuple(range(b))])
                                  for a in range(d)))
    for i,a,z in steps:charts[z]=m.compose(state[i][1][a],charts[a])
    inverse=tuple(m.inverse(F) for F in charts)
    returns=[]
    for (i,a) in breaks:
        z=q[i][a]
        returns.append(m.compose(m.compose(inverse[z],state[i][1][a]),charts[a]))
    R=tuple(returns[j] for j in order)
    assert sum(c.index(g) for g in R)==b-1
    assert m.product(R)==m.perm(b,[tuple(range(b))])
    return q,R

def lift(q,R):
    slots,steps,breaks,order=layout(q);d=len(q[0]);b=len(R[0]);I=tuple(range(b))
    assert len(R)==len(slots)
    assert all(len(h)==b and sorted(h)==list(range(b)) for h in R)
    assert sum(c.index(g) for g in R)==b-1
    assert m.product(R)==m.perm(b,[tuple(range(b))])
    maps=[[I]*d for _ in q]
    for j,h in zip(order,R):
        i,a=breaks[j];maps[i][a]=h
    normal=tuple((g,tuple(A)) for g,A in zip(q,maps))
    P,C=wproduct(normal);assert P==m.perm(d,[tuple(range(d))])
    charts=[I]
    for a in range(d-1):charts.append(m.compose(C[a],charts[-1]))
    assert m.compose(C[-1],charts[-1])==m.perm(b,[tuple(range(b))])
    invcharts=[m.inverse(F) for F in charts]
    original=tuple((g,tuple(m.compose(m.compose(invcharts[g[a]],A[a]),charts[a])
                              for a in range(d))) for g,A in normal)
    standard=(P,tuple(I if a<d-1 else m.perm(b,[tuple(range(b))]) for a in range(d)))
    assert wproduct(original)==standard
    assert extract(original)==(q,R)
    return original

def transport(q,R,word):
    state=lift(q,R)
    for move in word:state=wbraid(state,move)
    return extract(state)

def returns_at(q,b,chosen):
    slots,_,_,order=layout(q);I=tuple(range(b))
    return tuple(chosen.get(slots[j],I) for j in order)

def free_seed(outer,b,point=0):
    d=len(outer[0]);q=outer+(tuple(range(d)),)
    R=returns_at(q,b,{(len(outer),(point,)):m.perm(b,[tuple(range(b))])})
    return to_perms(lift(q,R))


def decorated(q,R):
    slots,_,_,order=layout(q);values=[None]*len(slots)
    for j,h in zip(order,R):values[j]=tuple(m.cycle_lengths(h))
    return q,tuple(values)


def canonical_cycle(C):
    at=C.index(min(C));return C[at:]+C[:at]


def decorated_move(data,move):
    q,values=data;i=abs(move)-1;slots=layout(q)[0]
    qp=m.hurwitz(q,i,1 if move>0 else -1);new={}
    for (branch,C),value in zip(slots,values):
        D=C;target=branch
        if move>0:
            if branch==i:
                target=i+1;D=tuple(m.inverse(q[i+1])[x] for x in C)
            elif branch==i+1:target=i
        else:
            if branch==i:target=i+1
            elif branch==i+1:
                target=i;D=tuple(q[i][x] for x in C)
        new[(target,canonical_cycle(D))]=value
    return qp,tuple(new[s] for s in layout(qp)[0])


def decorated_rotation(data,k):
    q,values=data;d=len(q[0]);qr=c.rotated(q,k)
    new={(i,canonical_cycle(tuple((x+k)%d for x in C))):v
         for (i,C),v in zip(layout(q)[0],values)}
    return qr,tuple(new[s] for s in layout(qr)[0])


def all_small_returns(q,b):
    assert b in (2,3)
    M=len(layout(q)[0]);I=tuple(range(b));C=m.perm(b,[tuple(range(b))])
    answer=[]
    for at in range(M):
        R=[I]*M;R[at]=C;answer.append(tuple(R))
    if b==3:
        Ts=[m.perm(b,[e]) for e in combinations(range(b),2)]
        for i,j in combinations(range(M),2):
            for T in Ts:
                R=[I]*M;R[i]=T;R[j]=m.compose(T,C)
                assert c.index(R[j])==1;answer.append(tuple(R))
    return tuple(answer)


def components(vertices,edges):
    remaining=set(vertices);answer=[]
    while remaining:
        seed=min(remaining);seen={seed};queue=deque([seed])
        while queue:
            x=queue.popleft()
            for y in edges[x]:
                assert y in vertices
                if y not in seen:seen.add(y);queue.append(y)
        remaining-=seen;answer.append(seen)
    return answer


def small_block_controls():
    rows=[]
    for outer in [
        (m.perm(2,[(0,1)]),tuple(range(2))),
        (m.perm(3,[(0,1)]),m.perm(3,[(1,2)])),
        (m.perm(4,[(0,1)]),m.perm(4,[(1,2)]),m.perm(4,[(2,3)]))
    ]:
        quotient_orbit,_=c.closure(outer);d=len(outer[0]);r=len(outer)
        moves=tuple(sign*(i+1) for i in range(r-1) for sign in (-1,1))
        for b in (2,3):
            X={(q,R) for q in quotient_orbit for R in all_small_returns(q,b)
               if all(c.index(g) for g in to_perms(lift(q,R)))}
            Y={decorated(q,R) for q,R in X};xe={};ye={};phase_controls=0
            literal_states=set();symmetric_base_control=None
            for q,R in X:
                state=lift(q,R);literal=to_perms(state);c.validate(literal)
                literal_states.add(literal);assert to_wreath(literal,d)==state
                data=decorated(q,R);xe[(q,R)]=[]
                f=3 if b==3 and sum(c.index(h)>0 for h in R)==2 else 1
                sy=sum(decorated_rotation(data,k)==data for k in range(d))
                assert c.stabilizer_size(literal)==b*sy//f
                if d==2 and b==3 and any(m.cycle_lengths(g)==[1]*4+[2] for g in literal):
                    normal=group_from_transposition(literal)
                    assert normal[:2]==(3,2) and c.h.group_closure(literal)==72
                    symmetric_base_control={"normal_block_size":3,"block_count":2,
                                            "exact_group_order":72,"group":"S3 wr C2"}
                for move in moves:
                    new=wbraid(state,move);key=extract(new);assert key in X
                    assert lift(*key)==new
                    assert to_perms(new)==m.hurwitz(literal,abs(move)-1,1 if move>0 else -1)
                    assert decorated(*key)==decorated_move(data,move)
                    xe[(q,R)].append(key)
                if q==outer:
                    twist=tuple(range(1,r))*r*d
                    qp,Rp=transport(q,R,twist)
                    assert qp==q
                    C=m.perm(b,[tuple(range(b))])
                    assert Rp==tuple(m.conjugate(m.inverse(C),h) for h in R)
                    phase_controls+=1
            for data in Y:ye[data]=[decorated_move(data,move) for move in moves]
            xc=components(X,xe);yc=components(Y,ye)
            projected=[frozenset(decorated(q,R) for q,R in comp) for comp in xc]
            assert set(projected)=={frozenset(comp) for comp in yc} and len(xc)==len(yc)
            independent=None
            if r==2:
                # Enumerate every transport of the first original factor, then
                # complete the second from the fixed product. This does not use
                # the return parametrization to generate its candidate domain.
                candidates=set();perms=tuple(permutations(range(b)))
                P=m.perm(b*d,[tuple(range(b*d))])
                for q in quotient_orbit:
                    for maps in product(perms,repeat=d):
                        A=to_perms(((q[0],maps),))[0]
                        B=m.compose(m.inverse(A),P)
                        assert to_wreath((A,B),d)[1][0]==q[1]
                        if c.index(A) and c.index(B) and c.index(A)+c.index(B)==b*d-1:
                            candidates.add((A,B))
                assert candidates==literal_states;independent=len(candidates)
            rows.append({"quotient_degree":d,"block_size":b,"quotient_orbit":len(quotient_orbit),
                         "complete_return_model_states":len(X),"decorated_states":len(Y),
                         "original_components":len(xc),"decorated_components":len(yc),
                         "literal_component_sizes":sorted(map(len,xc)),
                         "root_phase_controls":phase_controls,
                         "independent_wreath_first_factor_completion":independent,
                         "symmetric_base_control":symmetric_base_control})
    return rows


def marked_block(q,R):
    data=decorated(q,R);slots=layout(q)[0]
    marked=[(slot,h) for slot,h in zip(slots,data[1]) if len(h)<sum(h)]
    assert len(marked)==1
    (i,C),_=marked[0];assert len(C)==1 and q[i]==tuple(range(len(q[0])))
    return C[0]


def meridian(r,j):
    gather=list(range(r-1,j+1,-1))
    return tuple(gather+[j+1,j+1]+[-x for x in reversed(gather)])


def group_from_transposition(state):
    n=len(state[0])
    T=next(g for g in state if m.cycle_lengths(g)==[1]*(n-2)+[2])
    edge=next(cy for cy in c.cycles(T,False))
    seen={tuple(sorted(edge))};queue=deque(seen)
    generators=state+tuple(m.inverse(g) for g in state)
    while queue:
        a,b=queue.popleft()
        for g in generators:
            new=tuple(sorted((g[a],g[b])))
            if new not in seen:seen.add(new);queue.append(new)
    cs=c.v.one.v.components(n,list(seen));size=len(cs[0])
    assert all(len(C)==size for C in cs)
    d=n//size
    assert {frozenset(C) for C in cs}=={
        frozenset(range(a,n,d)) for a in range(d)}
    return size,d,len(seen)


def free_cyclic_controls():
    outer3=(m.perm(3,[(0,1)]),m.perm(3,[(1,2)]))
    outer5=(m.perm(5,[(0,1,2)]),m.perm(5,[(2,3,4)]))
    pure5=c.v.one.fixed_product(c.v.one.pure_canonical(2))
    _,pure7=c.v.one.pure_small()
    cases=[((m.perm(2,[(0,1)]),),2),((m.perm(3,[(0,1,2)]),),3),
           (outer3,2),(outer5,2),(outer5,3),(outer5,4),(pure5,2)]
    cases.extend((q,2) for q in pure7);rows=[]
    for outer,b in cases:
        d=len(outer[0]);r=len(outer)+1
        O,_=c.closure(outer);K=c.h.group_closure(outer)
        q0=outer+(tuple(range(d)),);Q,_=c.closure(q0)
        assert len(Q)==r*len(O)
        seed=free_seed(outer,b);orbit,_=c.closure(seed)
        assert len(orbit)==d*r*len(O)
        model=set()
        for q in Q:
            identity,= [i for i,g in enumerate(q) if c.index(g)==0]
            for a in range(d):
                R=returns_at(q,b,{(identity,(a,)):m.perm(b,[tuple(range(b))])})
                model.add((q,R))
        actual=set()
        for state in orbit:
            c.validate(state);q,R=extract(to_wreath(state,d))
            assert (q,R) in model and to_perms(lift(q,R))==state
            assert c.stabilizer_size(state)==b
            actual.add((q,R))
        assert actual==model
        point_generators=[]
        for j in range(r-1):
            word=meridian(r,j);perm=[]
            assert m.braid(q0,word)==q0
            for a in range(d):
                R=returns_at(q0,b,{(r-1,(a,)):m.perm(b,[tuple(range(b))])})
                qp,Rp=transport(q0,R,word);assert qp==q0
                perm.append(marked_block(qp,Rp))
            assert tuple(perm)==m.inverse(outer[j])
            point_generators.append(tuple(perm))
        assert c.h.group_closure(tuple(point_generators))==K
        predicted=b**d*K
        C=m.perm(b,[tuple(range(b))])
        cyclic={m.power(C,k) for k in range(b)}
        assert all(A in cyclic for _,maps in to_wreath(seed,d) for A in maps)
        paths={0:()};queue=deque([0])
        while queue:
            a=queue.popleft()
            for j,g in enumerate(outer):
                z=g[a]
                if z not in paths:paths[z]=(j+1,)+paths[a];queue.append(z)
        assert len(paths)==d
        certificates=[]
        for a,word in sorted(paths.items()):
            witness=word+(r,)+tuple(-x for x in reversed(word))
            value=c.group_word(seed,witness)
            assert value==m.perm(b*d,[tuple(a+d*t for t in range(b))])
            certificates.append({"block":a,"group_word_for_local_cycle":witness})
        exact=c.h.group_closure(seed) if predicted<=70000 else None
        if exact is not None:assert exact==predicted
        normal=None
        if b==2:
            normal=group_from_transposition(seed);assert normal[:2]==(2,d)
        rows.append({"quotient_degree":d,"block_size":b,"quotient_factors":len(outer),
                     "quotient_literal_orbit":len(O),"quotient_group_order":K,
                     "literal_lift_orbit":len(orbit),"absolute_lift_orbit":r*len(O),
                     "cyclic_stabilizer_order":b,"certified_wreath_order":predicted,
                     "exact_generated_group_closure":exact,
                     "group_basis":"cyclic wreath embedding + explicit words for every base generator + exact quotient closure",
                     "cyclic_wreath_embedding_checked":True,
                     "cyclic_base_generator_certificates":certificates,
                     "transposition_normal_closure":normal,
                     "literal_seed":seed,"orbit_sha256":c.digest(orbit),
                     "meridian_point_actions":point_generators})
        print(json.dumps({"free_cyclic_outer_degree":d,"block_size":b,
                          "outer_group":K,"lift_orbit":len(orbit)}),flush=True)
    return rows


def schreier_loop_controls():
    d=5;b=3
    q0=(m.perm(d,[(0,1,2)]),m.perm(d,[(2,3,4)]),tuple(range(d)))
    orbit,parent=c.closure(q0);paths={}
    for q in orbit:
        back=c.word_to_seed(q,parent);paths[q]=tuple(-x for x in reversed(back))
        assert m.braid(q0,paths[q])==q
    all_images=set();generators=[];loops=0
    for q in sorted(orbit):
        for move in (-1,1,-2,2):
            qp=m.braid(q,[move]);word=paths[q]+(move,)+tuple(-x for x in reversed(paths[qp]))
            assert m.braid(q0,word)==q0;loops+=1;image=[]
            for a in range(d):
                R=returns_at(q0,b,{(2,(a,)):m.perm(b,[tuple(range(b))])})
                newq,newR=transport(q0,R,word);assert newq==q0
                image.append(marked_block(newq,newR))
            image=tuple(image)
            if image not in all_images:
                all_images.add(image);generators.append({"braid_loop":word,"point_action":image})
    order=c.h.group_closure(tuple(row["point_action"] for row in generators))
    assert order==60
    return {"quotient_degree":5,"block_size":3,"quotient_orbit":len(orbit),
            "all_directed_schreier_loops":loops,"distinct_loop_actions":generators,
            "complete_stabilizer_image_order":order,"image":"A5 on the five marked blocks"}


def inner_four_obstruction():
    b=4;P=m.perm(b,[tuple(range(b))]);Ts=tuple(m.perm(b,[e]) for e in combinations(range(b),2))
    Tset=set(Ts)
    Rs={(A,B,C) for A in Ts for B in Ts
        if (C:=m.compose(m.inverse(m.compose(A,B)),P)) in Tset}
    assert len(Rs)==16
    inner,_=c.closure(min(Rs));assert inner==Rs
    q=(m.perm(2,[(0,1)]),tuple(range(2)))
    ordered=tuple(sorted(Rs));where={R:i for i,R in enumerate(ordered)};action=[]
    for R in ordered:
        qp,Rp=transport(q,R,[1,1]);assert qp==q
        action.append(where[Rp])
    assert sorted(m.cycle_lengths(tuple(action)))==[8,8]
    assert len({decorated(q,R) for R in Rs})==1
    original={to_perms(lift(q,R)) for R in Rs}
    remaining=set(original);sector=[]
    for _ in range(2):
        seed=min(remaining);orbit,_=c.closure(seed);remaining-=orbit
        assert len(orbit)==16 and c.h.group_closure(seed)==576
        assert c.block_divisors(seed)==[1,2,8] and c.stabilizer_size(seed)==1
        W=to_wreath(seed,2)
        assert all(c.index(g)%2==c.index(part[0])%2 for g,part in zip(seed,W))
        DT=next(g for g in seed if m.cycle_lengths(g)==[1]*4+[2,2])
        distances=tuple(sorted(min((a-b)%8,(b-a)%8) for a,b in c.cycles(DT,False)))
        assert distances in ((2,2),(2,4))
        sector.append({"literal_orbit":len(orbit),"group_order":576,"actual_group":"the parity kernel H",
                       "literal_seed":seed,"DT_cyclic_distances":distances,
                       "orbit_sha256":c.digest(orbit)})
    assert not remaining and len({r["DT_cyclic_distances"] for r in sector})==2
    # Complete the entire passport independently, including its other sector.
    P8=m.perm(8,[tuple(range(8))])
    direct={(A,D) for D in c.h.dt_class(8)
            if m.cycle_lengths(A:=m.compose(P8,D))==[2,2,4]}
    assert len(direct)==c.character_total([[2,2,4],[1]*4+[2,2]])==20
    remaining=set(direct);full=[]
    while remaining:
        seed=min(remaining);orbit,_=c.closure(seed);remaining-=orbit
        order=c.h.group_closure(seed);aut=c.stabilizer_size(seed)
        assert (len(orbit),order,aut) in ((16,576,1),(8,64,2))
        full.append({"literal_orbit":len(orbit),"absolute_orbit_size":len(orbit)*aut//8,
                     "group_order":order,"cyclic_stabilizer_order":aut,"literal_seed":seed,
                     "block_divisors":c.block_divisors(seed),"orbit_sha256":c.digest(orbit)})
    assert sorted(row["literal_orbit"] for row in full)==[8,16,16]
    return {"inner_degree":4,"quotient_degree":2,"all_three_transposition_inner_tuples":16,
            "one_full_inner_braid_component":True,"quotient_loop":[1,1],
            "quotient_loop_action_on_returns":action,"return_action_cycle_lengths":[8,8],
            "same_decoration_sector_components":sector,"complete_ordered_original_passport":20,
            "complete_full_passport_components":full}


def transposition_monodromy_controls():
    rows=[]
    for n,raw,outer in [
        (6,(m.perm(6,[(0,1),(3,4)]),m.perm(6,[(1,2),(4,5)]),m.perm(6,[(2,3)])),
         (m.perm(3,[(0,1)]),m.perm(3,[(1,2)]))),
        (10,(m.perm(10,[(0,1,2),(4,5,6)]),m.perm(10,[(2,3,4),(6,7,8)]),
             m.perm(10,[(8,9)])),
         (m.perm(5,[(0,1,2)]),m.perm(5,[(2,3,4)])))
    ]:
        seed=c.v.one.fixed_product(raw);orbit,_=c.closure(seed)
        a,d,edges=group_from_transposition(seed)
        assert a==n and d==1 and edges==comb(n,2)
        other=free_seed(outer,2);imprimitive,_=c.closure(other)
        assert not orbit&imprimitive
        passport=[m.cycle_lengths(g) for g in seed]
        total=c.character_total(passport)*c.type_placements(passport)
        assert len(orbit)+len(imprimitive)==total
        assert (len(orbit),len(imprimitive),total)==((216,27,243) if n==6 else (1800,75,1875))
        rows.append({"degree":n,"passport":passport,"primitive_literal_orbit":len(orbit),
                     "group":"S"+str(n),"normal_transposition_edges":edges,
                     "imprimitive_literal_orbit":len(imprimitive),"independent_full_total":total,
                     "literal_primitive_seed":seed,"primitive_orbit_sha256":c.digest(orbit)})
    return rows


def expect_rejected(function):
    try:function()
    except AssertionError:return True
    raise AssertionError("invalid control was admitted")


def negative_controls():
    q=(m.perm(2,[(0,1)]),tuple(range(2)));M=len(layout(q)[0])
    I=tuple(range(3));C=m.perm(3,[(0,1,2)]);T=m.perm(3,[(0,1)])
    assert M==3
    expect_rejected(lambda:lift(q,(m.inverse(C),I,I)))
    expect_rejected(lambda:lift(q,(T,T,C)))
    badq=(m.perm(3,[(0,1)]),m.perm(3,[(1,2)]),
          m.perm(3,[(0,1)]),m.perm(3,[(0,1)]))
    expect_rejected(lambda:layout(badq))
    expect_rejected(lambda:to_wreath((m.perm(6,[(0,1,2)]),),3))
    return {"wrong_fixed_return_product_rejected":True,"excess_return_index_rejected":True,
            "nonminimal_quotient_rejected":True,"nonpreserved_block_partition_rejected":True}

def canonical_divisor_controls(obstruction):
    rows=[]
    for n in (6,9,10,12,15):
        state=(m.perm(n,[tuple(range(n))]),)
        divisors=c.block_divisors(state)
        chosen=max(d for d in divisors if 1<d<n)
        for d in divisors:
            if 1<d<n:
                q,R=extract(to_wreath(state,d))
                assert to_perms(lift(q,R))==state
        least=min(k for k in range(2,n+1) if n%k==0)
        assert chosen==n//least
        rows.append({"degree":n,"all_block_divisors":divisors,
                     "canonical_quotient_degree":chosen,"canonical_block_size":least})
    filtered=[]
    for row in obstruction["complete_full_passport_components"]:
        proper=[d for d in row["block_divisors"] if 1<d<8]
        chosen=max(proper)
        assert chosen==(2 if row["group_order"]==576 else 4)
        filtered.append({"literal_orbit":row["literal_orbit"],"group_order":row["group_order"],
                         "canonical_quotient_degree":chosen,"canonical_block_size":8//chosen})
    return {"cyclic_nonunique_decomposition_controls":rows,"degree8_component_filter":filtered}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=HERE/"polynomial_lift_transport_verification.json")
    args=parser.parse_args()
    small=small_block_controls();free=free_cyclic_controls();loops=schreier_loop_controls()
    four=inner_four_obstruction();transpositions=transposition_monodromy_controls()
    negatives=negative_controls();canonical=canonical_divisor_controls(four)
    paths=["POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md","POLYNOMIAL_PASSPORT_COMPONENTS.md",
           "INDEX_TWO_HURWITZ_MIXTURES.md","HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md",
           "verify_polynomial_passport_components.py",
           "verify_index_two_hurwitz_mixtures.py","verify_mixed_hurwitz_isolation.py",
           "verify_terminal_mass_hurwitz_descent.py","verify_terminal_mass_projection_escape.py",
           "verify_paired_hurwitz_leaf_extraction.py","verify_unmarked_polynomial_hurwitz_reduction.py"]
    receipt={"schema_version":1,"kind":"polynomial_lift_return_transport","result":"PASS",
             "python":platform.python_version(),
             "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
             "small_block_return_domains":small,"free_cyclic_lifts":free,
             "schreier_loop_action":loops,"block_size_four_obstruction":four,
             "transposition_monodromy_controls":transpositions,"negative_controls":negatives,
             "canonical_divisor_controls":canonical,
             "scope":"Written all-degree return bijection, coupled transport, small-fibre sufficiency, "
                     "free-cyclic-branch component classification and lifted index-two catalogue. "
                     "Finite replay domains are explicit. Arbitrary primitive and remaining larger-fibre "
                     "component classification remains open; no extrapolation from these controls.",
             "independent_replay":False,"formal_verification":False,"external_review":False}
    args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"result":"PASS","small_return_states":sum(r["complete_return_model_states"] for r in small),
                      "free_cyclic_controls":len(free),"schreier_image_order":loops["complete_stabilizer_image_order"],
                      "block_size_four_return_orbits":[8,8],"receipt":str(args.output)},indent=2))


if __name__=="__main__":
    main()

