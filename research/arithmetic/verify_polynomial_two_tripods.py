#!/usr/bin/env python3
"""Exact controls for the written all-degree two-tripod classification."""
from collections import Counter, deque
from itertools import permutations
from functools import lru_cache
from math import gcd, factorial
from pathlib import Path
import hashlib
import json
import platform
import verify_polynomial_passport_components as c

if not __debug__:
    raise RuntimeError("tripod checker requires active assertions")

HERE=Path(__file__).resolve().parent
m=c.m
WORDS=((2,),(-2,),(1,1),(-1,-1),(1,2,2,-1),(1,-2,-2,-1))
INVERSES=(1,0,3,2,5,4)
CONTROLS=(
 (1,2,3,4,5),(1,2,3,4,7),(1,2,3,4,9),(1,2,3,4,11),
 (1,2,3,4,6),(1,2,3,5,7),(1,2,3,4,8),(1,3,5,7,9),
 (1,2,4,6,8),(2,4,6,8,10),(1,2,4,5,7),(1,2,4,5,9),
 (1,2,3,6,8),(1,2,3,6,7))
REPEATED=((1,1,1,1,2),(1,1,1,1,3),(1,1,1,2,2),
          (1,1,2,2,3),(1,2,2,3,4),(2,2,2,2,2),
          (2,2,4,4,6),(3,3,3,3,3))


@lru_cache(None)
def setup(lengths):
    assert len(lengths)==5 and tuple(sorted(lengths))==lengths
    assert lengths[0]>0
    cycles=[];offset=0
    for a in lengths:
        cycles.append(tuple(range(offset,offset+a)));offset+=a
    return tuple(cycles),m.perm(offset,cycles)


def states(lengths):
    setup(lengths)
    return {(order,delta) for order in permutations(range(5))
            for delta in range(lengths[order[0]])}


def orientation(state):
    order,_=state
    return (-1)**sum(order[i]>order[j] for i in range(5) for j in range(i+1,5))


def decode(lengths,state):
    cycles,g=setup(lengths)
    order,delta=state
    assert tuple(sorted(order))==tuple(range(5))
    center,a,b,d,e=order
    assert 0<=delta<lengths[center]
    u=(cycles[center][0],cycles[a][0],cycles[b][0])
    v=(cycles[center][delta],cycles[d][0],cycles[e][0])
    return g,m.perm(sum(lengths),[u]),m.perm(sum(lengths),[v])


def encode(state):
    g,u,v=state
    cycles=tuple(sorted(c.cycles(g),key=len))
    lengths=tuple(map(len,cycles))
    setup(lengths)
    where={x:i for i,C in enumerate(cycles) for x in C}
    U=[x for x in range(len(g)) if u[x]!=x]
    V=[x for x in range(len(g)) if v[x]!=x]
    iu={where[x] for x in U};iv={where[x] for x in V}
    assert len(U)==len(V)==len(iu)==len(iv)==3 and len(iu&iv)==1
    center=next(iter(iu&iv))
    x=next(x for x in U if where[x]==center)
    y=next(x for x in V if where[x]==center)
    a,b=where[u[x]],where[u[u[x]]]
    d,e=where[v[y]],where[v[v[y]]]
    delta=(cycles[center].index(y)-cycles[center].index(x))%lengths[center]
    return (center,a,b,d,e),delta


def forget_equal_labels(lengths,state):
    order,delta=state
    available={}
    for i,a in enumerate(lengths):available.setdefault(a,[]).append(i)
    canonical=tuple(available[lengths[i]].pop(0) for i in order)
    return canonical,delta


def B(lengths,state,inverse=False):
    (center,a,b,d,e),delta=state
    if delta:
        return (center,d,e,a,b),lengths[center]-delta
    if inverse:
        return (a,d,e,b,center),0
    return (e,center,d,a,b),0


def T(lengths,state,inverse=False):
    (center,a,b,d,e),delta=state
    if inverse:
        if delta>1:
            return (center,a,b,d,e),delta-1
        if delta==0 and lengths[center]>1:
            return (center,a,b,d,e),lengths[center]-1
        return (b,center,a,d,e),0
    if delta:
        return (center,a,b,d,e),(delta+1)%lengths[center]
    return (a,b,center,d,e),1%lengths[a]


def compact_move(lengths,state,index):
    if index<2:return B(lengths,state,index==1)
    if index<4:return T(lengths,state,index==3)
    # Q1 Q2^(+/-2) Q1^-1 = Q2^-1 Q1^(+/-2) Q2.
    return B(lengths,T(lengths,B(lengths,state,True),index==5))


def compact_orbits(lengths,raw=True):
    domain=states(lengths);remaining=set(domain);rows=[]
    edge_tests=0
    while remaining:
        seed=min(remaining);orbit={seed};queue=deque([seed])
        while queue:
            x=queue.popleft()
            literal=decode(lengths,x)
            assert encode(literal)==x
            for index,word in enumerate(WORDS):
                y=compact_move(lengths,x,index)
                assert y in domain and orientation(y)==orientation(x)
                assert compact_move(lengths,y,INVERSES[index])==x
                if raw:
                    assert encode(m.braid(literal,word))==y
                edge_tests+=1
                if y not in orbit:
                    orbit.add(y);queue.append(y)
        assert {orientation(x) for x in orbit}=={orientation(seed)}
        assert len(orbit)==12*sum(lengths)
        remaining-=orbit
        rows.append({"orientation":orientation(seed),"compact_orbit":len(orbit),
                     "coordinate_seed":seed})
    assert len(rows)==2
    return {"lengths":lengths,"complete_compact_domain":len(domain),
            "raw_signed_Schreier_edge_tests":edge_tests if raw else 0,
            "components":rows}


def repeated_controls(lengths):
    M=1
    for multiplicity in Counter(lengths).values():M*=factorial(multiplicity)
    labeled=states(lengths)
    domain={forget_equal_labels(lengths,x) for x in labeled}
    assert len(domain)==24*sum(lengths)//M
    seed=min(domain);seen={seed};queue=deque([seed]);edge_tests=0
    while queue:
        x=queue.popleft();literal=decode(lengths,x)
        assert forget_equal_labels(lengths,encode(literal))==x
        for index,word in enumerate(WORDS):
            y=forget_equal_labels(lengths,compact_move(lengths,x,index))
            assert y in domain
            assert forget_equal_labels(lengths,encode(m.braid(literal,word)))==y
            edge_tests+=1
            if y not in seen:
                seen.add(y);queue.append(y)
    assert seen==domain
    return {"lengths":lengths,"equal_length_label_group_order":M,
            "complete_unlabelled_domain":len(domain),
            "components":1,"raw_signed_Schreier_edge_tests":edge_tests}


def normal_support_certificate(seed):
    N=len(seed[0]);generators=seed+tuple(m.inverse(g) for g in seed)
    small=seed[1:]
    supports={frozenset(x for x in range(N) if g[x]!=x) for g in small}
    assert len(supports)==2 and all(len(s)==3 for s in supports)
    queue=deque(supports)
    while queue:
        support=queue.popleft()
        for g in generators:
            new=frozenset(g[x] for x in support)
            if new not in supports:
                supports.add(new);queue.append(new)
    components=[{x} for x in range(N)]
    for support in supports:
        hits=[S for S in components if S&support]
        components=[S for S in components if not S&support]+[set().union(*hits)]
    D=len(components);n=N//D
    assert all(len(S)==n for S in components)
    expected={frozenset(range(j,N,D)) for j in range(D)}
    assert {frozenset(S) for S in components}==expected
    # Independent order formulas and cycle parity in the written certificate.
    order=(factorial(n)//2)**D*D*(1+(n%2==0))
    return {"normal_three_cycle_supports":len(supports),"normal_support_blocks":D,
            "fibre_degree":n,"actual_group_order_from_written_support_theorem":order}


def first_position(state):
    state=tuple(state)
    position=next(i for i,g in enumerate(state) if len(c.cycles(g))==5)
    while position:
        state=m.braid(state,[position])
        position-=1
    return state


def literal_controls(lengths):
    N=sum(lengths);D=gcd(*lengths)
    coordinate_seed=next(x for x in sorted(states(lengths)) if orientation(x)==1)
    opposite=(coordinate_seed[0][:-2]+coordinate_seed[0][-2:][::-1],0)
    strict=len(set(lengths))==5
    coordinates=(coordinate_seed,opposite) if strict else (coordinate_seed,)
    seeds=[c.v.one.fixed_product(decode(lengths,x)) for x in coordinates]
    passport=[m.cycle_lengths(g) for g in seeds[0]]
    total=c.character_total(passport)*c.type_placements(passport)
    M=1
    for multiplicity in Counter(lengths).values():M*=factorial(multiplicity)
    placements=1 if lengths==(1,1,1,1,3) else 3
    assert total==24*N*N*placements//M
    union=set();rows=[];orientation_tests=0
    expected_blocks=sorted({d for d in range(1,D+1) if D%d==0}|{N})
    for seed in seeds:
        c.validate(seed)
        orbit,_=c.closure(seed)
        assert len(orbit)==total//(2 if strict else 1) and not orbit&union
        union.update(orbit)
        sign=orientation(encode(first_position(seed)))
        if strict:
            for state in orbit:
                assert orientation(encode(first_position(state)))==sign
                orientation_tests+=1
        assert c.block_divisors(seed)==expected_blocks
        assert c.stabilizer_size(seed)==1
        cert=normal_support_certificate(seed)
        assert cert["normal_support_blocks"]==D
        rows.append({"literal_seed":seed,"orientation":sign,"literal_orbit":len(orbit),
                     "absolute_orbit":len(orbit)//N,"deck":1,
                     "all_block_count_divisors":expected_blocks,
                     "normal_support_certificate":cert,"orbit_sha256":c.digest(orbit)})
    assert len(union)==total
    return {"lengths":lengths,"complete_literal_domain":total,
            "all_literal_orientation_tests":orientation_tests,
            "independent_total":"ordinary hook-character formula, all3 type placements",
            "components":rows}


def small_controls():
    def close(gens,N):
        I=tuple(range(N));seen={I};queue=deque([I])
        while queue:
            x=queue.popleft()
            for g in gens:
                y=m.compose(x,g)
                if y not in seen:
                    seen.add(y);queue.append(y)
        return seen
    three=m.perm(5,[(0,1,2)])
    five=(4,0,3,1,2)
    assert len(close((three,five),5))==60
    assert len(close((m.perm(4,[(0,1,2)]),m.perm(4,[(0,1,3)])),4))==12
    assert len(close((three,m.perm(5,[(0,3,4)])),5))==60
    rejected=[]
    for lengths,state in (
        ((1,3,2,4,5),((0,1,2,3,4),0)),
        ((1,2,3,4,5),((0,1,1,3,4),0)),
        ((1,2,3,4,5),((0,1,2,3,4),1))):
        try:decode(lengths,state)
        except AssertionError:rejected.append(True)
        else:raise AssertionError("invalid coordinate accepted")
    bad=m.perm(15,[(0,1,2,3,4),(5,6,7,8),(9,10,11),(12,13)])
    try:encode((bad,m.perm(15,[(0,1,5)]),m.perm(15,[(2,9,14)])))
    except AssertionError:rejected.append(True)
    else:raise AssertionError("non-tree support accepted")
    return {"zero_offset_position_group_order":60,"overlap_orders":[12,60],
            "deliberately_invalid_inputs_rejected":len(rejected)}


def main():
    inputs=[
      "POLYNOMIAL_TWO_TRIPOD_COMPONENTS.md","POLYNOMIAL_PASSPORT_COMPONENTS.md",
      "POLYNOMIAL_INTRINSIC_CLASSES_AND_LIFTING.md",
      "polynomial_lifting_obstructions_verification.json",
      "verify_polynomial_passport_components.py",
      "sources/jones-zvonkin-braid-cacti-2002.pdf"]
    # Record actual helper module paths rather than assuming a historical filename.
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m,c.h,c.v,c.v.one))
    inputs=list(dict.fromkeys(inputs))
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={"schema_version":1,"kind":"five_cycle_two_tripod_components",
             "result":"PASS","python":platform.python_version(),
             "checker_sha256":sha(Path(__file__).resolve()),
             "input_sha256":{p:sha(HERE/p) for p in inputs},
             "small_controls":small_controls(),
             "compact_controls":[compact_orbits(parts) for parts in CONTROLS],
             "repeated_length_controls":[repeated_controls(parts) for parts in REPEATED],
             "literal_controls":[literal_controls(parts) for parts in
                                 ((1,2,3,4,5),(1,2,3,4,6),(2,4,6,8,10))+REPEATED],
             "unrestricted_scope":"E1-E3 written tree, braid-normalization and normal-support proofs; finite domains do not establish the universal statement",
             "lifting_scope":"D7's existing degree15 quadratic values are retained; no uniform lifting-value formula for every strict five-part passport",
             "independent_replay":False,"formal_verification":False,"external_review":False}
    (HERE/"polynomial_two_tripods_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: two-tripod exact controls; "+
          str(sum(row["complete_compact_domain"] for row in receipt["compact_controls"]))+
          " compact states and "+
          str(sum(row["complete_literal_domain"] for row in receipt["literal_controls"]))+
          " literal tuples")


if __name__=="__main__":
    main()
