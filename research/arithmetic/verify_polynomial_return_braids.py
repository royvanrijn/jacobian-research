#!/usr/bin/env python3
"""Free-group return transport and the complete four-sheet gluing criterion.

The universal braid realization uses Artin's theorem as read in the retained
primary source. Exact word tests and complete bounded component controls are
recorded separately from that proof.
"""
from collections import deque
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json
import platform
import argparse
import verify_polynomial_lift_transport as l

if not __debug__:
    raise RuntimeError("return braid checker requires active assertions")

HERE=Path(__file__).resolve().parent
m=l.m
c=l.c


def mul(*words):
    out=[]
    for word in words:
        for x in word:
            assert isinstance(x,int) and x
            if out and out[-1]==-x:
                out.pop()
            else:
                out.append(x)
    return tuple(out)


def inv(word):
    return tuple(-x for x in word[::-1])


def substitute(word,images):
    return mul(*(images[x-1] if x>0 else inv(images[-x-1]) for x in word))


def word_wmul(A,B):
    q,x=A
    s,y=B
    return m.compose(q,s),tuple(mul(x[s[a]],y[a]) for a in range(len(q)))


def word_winv(A):
    q,x=A
    p=m.inverse(q)
    return p,tuple(inv(x[p[a]]) for a in range(len(q)))


def word_wproduct(state):
    d=len(state[0][0])
    value=(tuple(range(d)),((),)*d)
    for A in state:
        value=word_wmul(value,A)
    return value


def word_lift(q,R):
    slots,_,breaks,order=l.layout(q)
    d=len(q[0])
    maps=[[()]*d for _ in q]
    for j,h in zip(order,R):
        i,a=breaks[j]
        maps[i][a]=h
    normal=tuple((g,tuple(A)) for g,A in zip(q,maps))
    P,C=word_wproduct(normal)
    charts=[()]
    for a in range(d-1):
        charts.append(mul(C[a],charts[-1]))
    assert mul(C[-1],charts[-1])==mul(*R)
    W=tuple((g,tuple(mul(inv(charts[g[a]]),A[a],charts[a]) for a in range(d)))
            for g,A in normal)
    assert word_wproduct(W)==(P,tuple(() if a<d-1 else mul(*R) for a in range(d)))
    return W


def word_extract(W):
    q=tuple(x[0] for x in W)
    _,steps,breaks,order=l.layout(q)
    d=len(q[0])
    charts=[None]*d
    charts[0]=()
    for i,a,z in steps:
        charts[z]=mul(W[i][1][a],charts[a])
    returns=[]
    for i,a in breaks:
        returns.append(mul(inv(charts[q[i][a]]),W[i][1][a],charts[a]))
    return q,tuple(returns[j] for j in order)


def word_braid(W,move):
    i=abs(move)-1
    A,B=W[i:i+2]
    if move>0:
        pair=(B,word_wmul(word_wmul(word_winv(B),A),B))
    else:
        pair=(word_wmul(word_wmul(A,B),word_winv(A)),A)
    return W[:i]+pair+W[i+2:]


def peripheral_generator(word):
    # A reduced conjugate of a positive basis generator has a unique
    # one-letter cyclic reduction.
    start=0
    stop=len(word)
    while stop-start>1 and word[start]==-word[stop-1]:
        start+=1
        stop-=1
    assert stop-start==1 and word[start]>0,word
    return word[start]


def formal_transport(q,word):
    M=len(l.layout(q)[0])
    basis=tuple((i,) for i in range(1,M+1))
    W=word_lift(q,basis)
    assert word_extract(W)==(q,basis)
    for move in word:
        W=word_braid(W,move)
    qp,images=word_extract(W)
    assert mul(*images)==tuple(range(1,M+1))
    assert sorted(peripheral_generator(w) for w in images)==list(range(1,M+1))
    assert word_lift(qp,images)==W
    return qp,images


def evaluate(word,values):
    b=len(values[0])
    out=tuple(range(b))
    for x in word:
        value=values[x-1] if x>0 else m.inverse(values[-x-1])
        out=m.compose(out,value)
    return out


def erase(images,active):
    chosen={j:k+1 for k,j in enumerate(active)}
    variables=tuple((chosen[i],) if i in chosen else () for i in range(len(images)))
    values=tuple(substitute(w,variables) for w in images)
    return tuple(w for w in values if w)


def return_rotation(R):
    b=len(R[0])
    return min(c.rotated(R,k) for k in range(b))


def quartic_inner():
    b=4
    P=m.perm(b,[tuple(range(b))])
    Ts=tuple(m.perm(b,[e]) for e in combinations(range(b),2))
    TT=set(Ts)
    triples={(A,B,C) for A in Ts for B in Ts
             if (C:=m.compose(m.inverse(m.compose(A,B)),P)) in TT}
    assert len(triples)==16
    canonical=tuple(sorted({return_rotation(R) for R in triples}))
    assert len(canonical)==4
    assert all(len({c.rotated(R,k) for k in range(4)})==4 for R in triples)
    where={R:i for i,R in enumerate(canonical)}
    images=[]
    for move in (1,2):
        image=tuple(where[return_rotation(m.braid(R,[move]))] for R in canonical)
        assert m.cycle_lengths(image)==[1,3]
        images.append(image)
    assert c.h.group_closure(tuple(images))==12
    return triples,canonical,tuple(images)


def active_positions(data):
    q,types=data
    slots,_,_,order=l.layout(q)
    by_slot=dict(zip(slots,types))
    active=tuple(k for k,j in enumerate(order) if by_slot[slots[j]]!=(1,1,1,1))
    assert len(active)==3
    assert all(by_slot[slots[order[k]]]==(1,1,2) for k in active)
    return active


def fill_triple(q,active,triple):
    R=[tuple(range(4))]*len(l.layout(q)[0])
    for k,T in zip(active,triple):
        R[k]=T
    return tuple(R)


def quartic_edge(data,move,canonical):
    q,_=data
    target=l.decorated_move(data,move)
    active=active_positions(data)
    where={R:i for i,R in enumerate(canonical)}
    image=[]
    qp,words=formal_transport(q,[move])
    assert qp==target[0]
    erasure=erase(words,active)
    assert len(erasure)==3 and mul(*erasure)==(1,2,3)
    assert sorted(peripheral_generator(w) for w in erasure)==[1,2,3]
    for triple in canonical:
        R=fill_triple(q,active,triple)
        newq,newR=l.transport(q,R,[move])
        assert newq==target[0] and l.decorated(newq,newR)==target
        assert tuple(evaluate(w,R) for w in words)==newR
        nonidentity=tuple(T for T in newR if c.index(T))
        assert tuple(evaluate(w,triple) for w in erasure)==nonidentity
        image.append(where[return_rotation(nonidentity)])
    image=tuple(image)
    assert c.index(image)%2==0
    return target,image


def perm_group(generators,n):
    I=tuple(range(n))
    seen={I}
    queue=deque([I])
    while queue:
        x=queue.popleft()
        for g in generators:
            y=m.compose(g,x)
            if y not in seen:
                seen.add(y)
                queue.append(y)
    return seen


def point_orbits(group,n):
    remaining=set(range(n))
    out=[]
    while remaining:
        a=min(remaining)
        orbit={g[a] for g in group}
        remaining-=orbit
        out.append(tuple(sorted(orbit)))
    return tuple(out)


def classify_quartic(data,canonical):
    """Complete four-state criterion above one triple-T decorated component."""
    r=len(data[0])
    moves=tuple(sign*(i+1) for i in range(r-1) for sign in (-1,1))
    paths={data:tuple(range(4))}
    parent={data:None}
    queue=deque([data])
    edge_count=0
    loops=set()
    while queue:
        vertex=queue.popleft()
        for move in moves:
            target,image=quartic_edge(vertex,move,canonical)
            reached=m.compose(image,paths[vertex])
            if target not in paths:
                paths[target]=reached
                parent[target]=(vertex,move)
                queue.append(target)
            loop=m.compose(m.inverse(paths[target]),reached)
            loops.add(loop)
            edge_count+=1
    group=perm_group(loops,4)
    assert len(group) in (1,2,3,4,12)
    orbits=point_orbits(group,4)
    assert sorted(map(len,orbits)) in ([1,1,1,1],[2,2],[1,3],[4])
    return {"decorated_vertices":len(paths),"directed_edges":edge_count,
            "return_label_group_order":len(group),
            "return_label_generators":sorted(loops),
            "label_orbits":orbits,
            "literal_component_sizes":sorted(4*len(O)*len(paths) for O in orbits)},paths


def quotient_seed(d,r):
    q=tuple(m.perm(d,[(a,a+1)]) for a in range(d-1))
    assert len(q)<=r
    return q+(tuple(range(d)),)*(r-len(q))


def word_controls():
    rows=[]
    for d,r in ((2,2),(2,3),(3,2),(3,3),(3,4),(4,3)):
        seed=quotient_seed(d,r)
        orbit,_=c.closure(seed)
        edges=0
        inverses=0
        relations=0
        phase_controls=0
        specializations=0
        for q in sorted(orbit):
            M=len(l.layout(q)[0])
            basis=tuple((i,) for i in range(1,M+1))
            for move in (sign*(i+1) for i in range(r-1) for sign in (-1,1)):
                qp,images=formal_transport(q,[move])
                backq,back=formal_transport(qp,[-move])
                assert backq==q
                assert tuple(substitute(w,images) for w in back)==basis
                assert tuple(substitute(w,back) for w in images)==basis
                edges+=1
                inverses+=2
                for b,nonidentity in (
                    (4,(m.perm(4,[(0,1)]),m.perm(4,[(1,2)]),m.perm(4,[(2,3)]))),
                    (5,(m.perm(5,[(0,1,2)]),m.perm(5,[(2,3)]),m.perm(5,[(3,4)]))),
                    (6,(m.perm(6,[(0,1,2)]),m.perm(6,[(2,3,4)]),m.perm(6,[(4,5)]))),
                ):
                    R=nonidentity+(tuple(range(b)),)*(M-3)
                    q2,R2=l.transport(q,R,[move])
                    assert q2==qp and R2==tuple(evaluate(w,R) for w in images)
                    specializations+=1
            for i in range(1,r-1):
                assert formal_transport(q,[i,i+1,i])==formal_transport(q,[i+1,i,i+1])
                relations+=1
            for i in range(1,r):
                for j in range(i+2,r):
                    assert formal_transport(q,[i,j])==formal_transport(q,[j,i])
                    relations+=1
            twist=tuple(range(1,r))*r*d
            qp,images=formal_transport(q,twist)
            Z=tuple(range(1,M+1))
            assert qp==q and images==tuple(mul(inv(Z),(i,),Z) for i in range(1,M+1))
            phase_controls+=1
        rows.append({"quotient_degree":d,"factor_positions":r,"complete_quotient_orbit":len(orbit),
                     "signed_edges":edges,"two_sided_inverse_tests":inverses,
                     "braid_relation_tests":relations,"full_twist_root_phase_tests":phase_controls,
                     "permutation_specializations_degrees":[4,5,6],
                     "permutation_specialization_tests":specializations})
    q=quotient_seed(2,2)
    concrete=[{"original_word":word,"return_images":formal_transport(q,word)[1]}
              for word in ((1,),(-1,),(1,1))]
    return {"complete_domains":rows,"degree2_literal_substitutions":concrete}


def quartic_partition_controls():
    P=m.perm(4,[tuple(range(4))])
    group=tuple(permutations(range(4)))
    rows=[]
    for index_two,name in (((1,3),"C,T"),((2,2),"D,T")):
        tuples={(A,B) for A in group if tuple(m.cycle_lengths(A))==index_two
                for B in group if m.cycle_lengths(B)==[1,1,2] and m.compose(A,B)==P}
        count=4 if name=="C,T" else 2
        assert len(tuples)==count
        assert len({return_rotation(R) for R in tuples})==1
        rows.append({"ordered_return_pattern":name,"literal_fibre":count,
                     "root_phase_orbits":1,"literal_pairs":sorted(tuples)})
    triples,canonical,inner=quartic_inner()
    A4=perm_group(inner,4)
    subgroups={frozenset([tuple(range(4))])}
    pending=deque(subgroups)
    while pending:
        H=pending.popleft()
        for g in A4-H:
            K=frozenset(perm_group(tuple(H)+(g,),4))
            if K not in subgroups:
                subgroups.add(K)
                pending.append(K)
    assert len(subgroups)==10
    census={}
    for H in subgroups:
        census[len(H)]=census.get(len(H),0)+1
        assert sorted(map(len,point_orbits(H,4))) in ([1,1,1,1],[2,2],[1,3],[4])
    assert census=={1:1,2:3,3:4,4:1,12:1}
    return {"two_return_patterns":rows,"all_triple_T_tuples":sorted(triples),
            "four_rotation_class_representatives":canonical,
            "inner_Q1_Q2_actions":inner,"inner_absolute_image_order":12,
            "complete_A4_subgroup_order_census":census}


def quartic_component_controls(canonical):
    rows=[]
    total_vertices=0
    for d,r in ((2,2),(2,3),(3,2),(3,3),(4,3)):
        q=quotient_seed(d,r)
        M=len(l.layout(q)[0])
        seen=set()
        domains=[]
        for active in combinations(range(M),3):
            R=fill_triple(q,active,canonical[0])
            if not all(c.index(g) for g in l.to_perms(l.lift(q,R))):
                continue
            data=l.decorated(q,R)
            if data in seen:
                continue
            result,paths=classify_quartic(data,canonical)
            seen.update(paths)
            comparisons=[]
            for label_orbit in result["label_orbits"]:
                expected=set()
                for vertex,S in paths.items():
                    for label in label_orbit:
                        for k in range(4):
                            R=fill_triple(vertex[0],active_positions(vertex),
                                          c.rotated(canonical[S[label]],k))
                            literal=l.to_perms(l.lift(vertex[0],R))
                            c.validate(literal)
                            expected.add(literal)
                seed=min(expected)
                actual,_=c.closure(seed)
                assert actual==expected
                aut=c.stabilizer_size(seed)
                assert len(actual)*aut%(4*d)==0
                comparisons.append({"label_orbit":label_orbit,"literal_orbit":len(actual),
                                    "cyclic_stabilizer_order":aut,
                                    "absolute_orbit":len(actual)*aut//(4*d),
                                    "all_block_divisors":c.block_divisors(seed),
                                    "literal_seed":seed,"orbit_sha256":c.digest(actual)})
                total_vertices+=len(actual)
            assert sorted(x["literal_orbit"] for x in comparisons)==result["literal_component_sizes"]
            result.update({"base_active_positions":active,"base_literal_quotient":q,
                           "direct_original_component_controls":comparisons})
            domains.append(result)
        # Completeness at the base quotient: all ways to choose three active
        # return slots were inspected, and every decorated component meets it.
        rows.append({"quotient_degree":d,"factor_positions":r,
                     "all_base_triples_inspected":len(tuple(combinations(range(M),3))),
                     "decorated_components":domains})
        print(json.dumps({"quartic_quotient_degree":d,"factor_positions":r,
                          "decorated_components":len(domains),
                          "label_group_orders":[row["return_label_group_order"] for row in domains]}),
              flush=True)
    assert total_vertices==26576
    degree16=rows[-1]
    chosen=[row for row in degree16["decorated_components"]
            if all(m.cycle_lengths(tuple(g))==[1]*6+[2]*5
                   for g in row["direct_original_component_controls"][0]["literal_seed"])]
    assert len(chosen)==2 and all(row["decorated_vertices"]==64
                                 and row["return_label_group_order"]==3 for row in chosen)
    original=[entry for row in chosen for entry in row["direct_original_component_controls"]]
    assert sorted(row["literal_orbit"] for row in original)==[256,256,768,768]
    assert sorted(row["absolute_orbit"] for row in original)==[16,16,48,48]
    assert all(row["all_block_divisors"]==[1,4,16] and row["cyclic_stabilizer_order"]==1 for row in original)
    assert all(all(m.cycle_lengths(tuple(g))==[1]*6+[2]*5 for g in row["literal_seed"]) for row in original)
    return {"complete_domains":rows,"direct_original_literal_states":total_vertices,
            "degree16_three_2power5_sector":{"decorated_vertices":128,"literal_states":2048,
              "literal_components":[256,256,768,768],"absolute_components":[16,16,48,48],
              "all_block_divisors":[1,4,16]}}


def cyclic_outer_inheritance_controls():
    inputs=json.loads((HERE/"polynomial_passport_components_inputs.json").read_text())
    certificates=inputs["degree11_pair_certificates"]
    assert len(certificates)==6
    lower_labels=[]
    rows=[]
    for item in certificates:
        h=tuple(tuple(g) for g in item["representative"])
        c.validate(h)
        assert m.cycle_lengths(c.group_word(h,item["group_three_cycle_word"]))==[1]*8+[3]
        orbit,_=c.closure(h)
        assert len(orbit)==22
        lower_labels.append(min(orbit))
        for d in (2,3):
            q=(m.perm(d,[tuple(range(d))]),)+(tuple(range(d)),)*len(h)
            R=l.returns_at(q,11,{(i+1,(0,)):g for i,g in enumerate(h)})
            lift=l.to_perms(l.lift(q,R))
            c.validate(lift)
            assert m.cycle_lengths(lift[0])==[d]*11
            for g,old in zip(lift[1:],h):
                assert m.cycle_lengths(g)==[1]*(11*(d-1))+m.cycle_lengths(old)
            q2,R2=l.extract(l.to_wreath(lift,d))
            assert q2==q and tuple(g for g in R2 if c.index(g))==h
            W=l.to_wreath(lift,d)
            assert all(c.index(A)%2==0 for _,maps in W for A in maps)
            # Each old generator is local on root block zero. Conjugation by
            # the outer generator supplies its local conjugates on all blocks.
            local_certificates=[]
            for a in range(d):
                for i,old in enumerate(h,1):
                    value=m.conjugate(m.power(lift[0],a),lift[i])
                    fixed={x for x in range(11*d) if value[x]!=x}
                    assert fixed and {x%d for x in fixed}=={a}
                    local_certificates.append({"block":a,"inner_factor":i,
                                               "group_word":(1,)*a+(i+1,)+(-1,)*a})
            if d%2:
                # The unbounded spin assertion is H5's theorem. This is an
                # exact evaluation of its type formula, not a new Clifford
                # calculation or an independent proof of H5.
                exponent=sum((ell*ell-1)//8 for g in lift for ell in m.cycle_lengths(g))
                exponent+=((11*d)**2-1)//8
                assert exponent%2==0
            rows.append({"outer_degree":d,"inner_degree":11,"literal_seed":lift,
                         "inner_rigid_component_label":min(orbit),
                         "group":"A11 wr C"+str(d),
                         "group_basis":"R5 wreath embedding + H6 A11 certificate + local conjugates",
                         "standard_spin_from_H5":1 if d%2 else None,
                         "local_factor_certificates":local_certificates})
    assert len(set(lower_labels))==6
    return {"six_distinct_inner_components":True,"controls":rows,
            "assurance":"No large lifted group closure was attempted; wreath equality uses R5 and H6."}


def expect_rejected(function):
    try:
        function()
    except AssertionError:
        return True
    raise AssertionError("invalid Artin control was admitted")


def check_artin(images):
    M=len(images)
    assert mul(*images)==tuple(range(1,M+1))
    assert sorted(peripheral_generator(w) for w in images)==list(range(1,M+1))


def negative_controls():
    expect_rejected(lambda:check_artin(((-2,1,2),(2,))))
    expect_rejected(lambda:check_artin(((1,2),(-2,))))
    expect_rejected(lambda:check_artin(((-1,),(2,))))
    return {"peripheral_but_wrong_boundary_rejected":True,
            "boundary_but_nonperipheral_rejected":True,
            "negative_generator_rejected":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=HERE/"polynomial_return_braids_verification.json")
    args=parser.parse_args()
    words=word_controls()
    quartic=quartic_partition_controls()
    _,canonical,_=quartic_inner()
    components=quartic_component_controls(canonical)
    inheritance=cyclic_outer_inheritance_controls()
    negatives=negative_controls()
    paths=["POLYNOMIAL_PASSPORT_QUARTIC_RETURN_CLASSIFICATION.md",
           "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md","POLYNOMIAL_PASSPORT_COMPONENTS.md",
           "verify_polynomial_lift_transport.py","verify_polynomial_passport_components.py",
           "polynomial_passport_components_inputs.json","sources/artin-theory-of-braids-1947.pdf"]
    receipt={"schema_version":1,"kind":"polynomial_return_braids_and_quartic_gluing","result":"PASS",
             "python":platform.python_version(),
             "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
             "free_group_controls":words,"quartic_inner_and_subgroups":quartic,
             "complete_quartic_gluing_controls":components,
             "cyclic_outer_inheritance":inheritance,"negative_controls":negatives,
             "scope":"Uniform R1-R5 written proofs use Artin's actual characterization, the literal "
                     "tree construction and the finite quartic return catalogue. Bounded controls "
                     "do not prove other quotient or larger-fibre classifications. The complete "
                     "arbitrary-passport goal remains open.",
             "independent_replay":False,"formal_verification":False,"external_review":False}
    args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"result":"PASS","direct_quartic_original_states":components["direct_original_literal_states"],
                      "cyclic_outer_inheritance_controls":len(inheritance["controls"]),
                      "receipt":str(args.output)},indent=2))


if __name__=="__main__":
    main()
