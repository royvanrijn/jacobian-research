#!/usr/bin/env python3
"""Intrinsic branch classes and genuine quadratic/cubic lifting obstructions.

Finite permutation-cover and Clifford certificates are separate from the
written all-degree statements. No stable Conway--Parker theorem is applied.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import product
from math import gcd, lcm, factorial
from pathlib import Path
import argparse
import hashlib
import json
import platform
import re
import verify_polynomial_return_braids as r

if not __debug__:
    raise RuntimeError("lifting checker requires active assertions")

HERE=Path(__file__).resolve().parent
c=r.c
m=r.m
l=r.l


def order(g):
    return lcm(*m.cycle_lengths(g))


def parse_mtx(text):
    values=list(map(int,text.split()))
    assert values[:4]==[12,1,45,1]
    perm=tuple(x-1 for x in values[4:])
    assert len(perm)==45 and sorted(perm)==list(range(45))
    return perm


def parse_gap(text):
    assert re.fullmatch(r"\s*a\d+\s*:=\s*(?:\([\d,]+\))+\s*;\s*",text)
    cycles=[tuple(int(x)-1 for x in part.split(","))
            for part in re.findall(r"\(([\d,]+)\)",text)]
    return m.perm(7,cycles)


@lru_cache(None)
def cover_a7():
    A=parse_mtx((HERE/"sources/atlas-3a7-45a-generator-1.mtx").read_text())
    B=parse_mtx((HERE/"sources/atlas-3a7-45a-generator-2.mtx").read_text())
    a=parse_gap((HERE/"sources/atlas-a7-7-generator-1.gap").read_text())
    b=parse_gap((HERE/"sources/atlas-a7-7-generator-2.gap").read_text())
    I,J=tuple(range(45)),tuple(range(7))
    projection={I:J}
    queue=deque([I])
    generators=((A,a),(B,b),(m.inverse(A),m.inverse(a)),(m.inverse(B),m.inverse(b)))
    edge_tests=0
    while queue:
        x=queue.popleft()
        q=projection[x]
        for u,v in generators:
            y=m.compose(x,u)
            image=m.compose(q,v)
            if y in projection:
                assert projection[y]==image
            else:
                projection[y]=image
                queue.append(y)
            edge_tests+=1
    fibres={}
    for x,q in projection.items():
        fibres.setdefault(q,[]).append(x)
    assert len(projection)==7560 and len(fibres)==2520 and edge_tests==30240
    assert all(c.index(q)%2==0 and len(F)==3 for q,F in fibres.items())
    kernel=tuple(sorted(fibres[J]))
    assert len(kernel)==3
    z=next(x for x in kernel if x!=I)
    assert m.power(z,3)==I
    assert all(m.compose(z,g)==m.compose(g,z) for g in (A,B))
    assert order(A)==3 and order(B)==5 and order(m.compose(A,B))==7
    # Orders3,5,7 force trivial abelianization, proving this extension stem.
    indices={I:0,z:1,m.power(z,2):2}
    return fibres,indices,(A,B),edge_tests


def a7_canonical(g,k):
    fibres,_,_,_=cover_a7()
    values=[x for x in fibres[g] if m.power(x,k)==tuple(range(45))]
    assert len(values)==1
    return values[0]


def cubic_a7(state):
    _,indices,_,_=cover_a7()
    out=tuple(range(45))
    for g in state:
        out=m.compose(out,a7_canonical(g,order(g)))
    P=m.product(state)
    out=m.compose(out,m.inverse(a7_canonical(P,order(P))))
    return indices[out]


def clifford_word(roots):
    value={0:1}
    for a,b in roots:
        value=c.multiply_vector(value,a,b)
    return value


def centralizer_reduction_controls():
    x=((0,1),(2,3))
    y=((0,2),(1,3))
    # Four unit roots for the disjoint-transposition commutator; eight for Q8.
    symmetric=clifford_word(((0,1),(2,3),(0,1),(2,3)))
    quaternion=clifford_word(x+y+x[::-1]+y[::-1])
    assert symmetric=={0:-4} and quaternion=={0:-16}
    X,Y=m.perm(7,x),m.perm(7,y)
    assert m.compose(X,Y)==m.compose(Y,X)
    fibres,_,_,_=cover_a7()
    dt=tuple(c.h.dt_class(7))
    checks=0
    for g in dt+(m.perm(7,[tuple(range(7))]),):
        lift=a7_canonical(g,order(g))
        for h,F in fibres.items():
            if m.compose(g,h)==m.compose(h,g):
                assert m.compose(lift,F[0])==m.compose(F[0],lift)
                checks+=1
    return {"Sn_disjoint_transposition_raw_commutator":symmetric,
            "A7_Klein_four_raw_commutator":quaternion,
            "normalization_denominators":[4,16],
            "A7_all_branch_centralizer_tests_in_cubic_cover":checks}


def a7_component_controls():
    P=m.perm(7,[tuple(range(7))])
    dt=tuple(c.h.dt_class(7))
    D=set(dt)
    direct={(A,B,C) for A in dt for B in dt
            if (C:=m.compose(m.inverse(m.compose(A,B)),P)) in D}
    assert len(direct)==c.character_total([[1]*3+[2,2]]*3)==392
    fibres,_,_,_=cover_a7()
    A7=set(fibres)
    assert {m.conjugate(h,dt[0]) for h in A7}==D
    assert len({m.conjugate(h,P) for h in A7})==360
    remaining=set(direct)
    rows=[]
    all_values=Counter()
    edge_tests=0
    while remaining:
        seed=min(remaining)
        orbit,_=c.closure(seed)
        assert orbit<=direct
        remaining-=orbit
        value=cubic_a7(seed)
        for state in orbit:
            assert cubic_a7(state)==value
            for move in (-1,1,-2,2):
                assert cubic_a7(m.braid(state,[move]))==value
                edge_tests+=1
        all_values[value]+=len(orbit)
        group=r.perm_group(seed,7)
        assert (len(group),len(orbit),value) in ((2520,147,1),(2520,147,2),(168,49,0))
        rows.append({"literal_seed":seed,"literal_orbit":len(orbit),
                     "group_order":len(group),"cubic_value":value,
                     "actual_group_sha256":c.digest({tuple([g]) for g in group}),
                     "orbit_sha256":c.digest(orbit)})
    assert dict(all_values)=={0:98,1:147,2:147}
    assert len({row["actual_group_sha256"] for row in rows if row["group_order"]==168})==2
    return {"complete_literal_domain":392,"signed_braid_value_tests":edge_tests,
            "same_A7_DT_class_size":105,"fixed_infinity_A7_class_size":360,
            "cubic_value_counts":dict(all_values),"components":rows}


MATCHINGS=tuple(frozenset((frozenset(p),frozenset(set(range(4))-set(p))))
                for p in ((0,1),(0,2),(0,3)))


def s4_affine(u):
    eta=tuple(MATCHINGS.index(frozenset(frozenset(u[x] for x in e) for e in matching))
              for matching in MATCHINGS)
    t=eta[0]
    eps=(eta[1]-t)%3
    assert eps in (1,2) and eta[2]==(t+2*eps)%3
    return t,eps


def rotate(j,v):
    for _ in range(j%4):
        v=((-v[1])%3,v[0])
    return v


MATRICES=tuple(tuple(rotate(j,(int(k==0),int(k==1))) for k in range(2)) for j in range(4))


def h8_projection(g):
    q,maps=l.to_wreath((g,),2)[0]
    assert c.index(g)%2==c.index(q)%2
    v=[0,0]
    columns=[]
    for a,u in enumerate(maps):
        t,eps=s4_affine(u)
        v[q[a]]=t
        column=[0,0]
        column[q[a]]=eps
        columns.append(tuple(column))
    return tuple(v),MATRICES.index(tuple(columns))


def heisenberg_multiply(X,Y):
    v,j,z=X
    w,k,t=Y
    w=rotate(j,w)
    omega=(v[0]*w[1]-v[1]*w[0])%3
    return tuple((a+b)%3 for a,b in zip(v,w)),(j+k)%4,(z+t+2*omega)%3


def heisenberg_inverse(X):
    v,j,z=X
    return rotate(-j,((-v[0])%3,(-v[1])%3)),(-j)%4,(-z)%3


def heisenberg_power(X,k):
    out=((0,0),0,0)
    for _ in range(k):
        out=heisenberg_multiply(out,X)
    return out


def h8_canonical(g):
    k=order(g)
    assert gcd(k,3)==1
    v,j=h8_projection(g)
    z=heisenberg_power((v,j,0),k)[2]
    lift=(v,j,(-z*pow(k,-1,3))%3)
    assert heisenberg_power(lift,k)==((0,0),0,0)
    return lift


def cubic_h8(state):
    out=((0,0),0,0)
    for g in state:
        out=heisenberg_multiply(out,h8_canonical(g))
    target=h8_canonical(m.product(state))
    assert out[:2]==target[:2]
    return (out[2]-target[2])%3


def h8_controls():
    old=json.loads((HERE/"polynomial_lift_transport_verification.json").read_text())
    inputs=old["block_size_four_obstruction"]["same_decoration_sector_components"]
    states=[tuple(tuple(g) for g in row["literal_seed"]) for row in inputs]
    H=r.perm_group(states[0],8)
    assert len(H)==576 and r.perm_group(states[1],8)==H
    quotient={h8_projection(g) for g in H}
    assert len(quotient)==36
    quotient_counts=Counter(h8_projection(g) for g in H)
    assert set(quotient_counts.values())=={16}
    for g in H:
        for h in states[0]+tuple(m.inverse(h) for h in states[0]):
            q=h8_projection(m.compose(g,h))
            a,b=h8_projection(g),h8_projection(h)
            assert q==heisenberg_multiply(a+(0,),b+(0,))[:2]
    E={q+(z,) for q in quotient for z in range(3)}
    assert len(E)==108
    zero=((0,0),0,0)
    for x in E:
        assert heisenberg_multiply(x,heisenberg_inverse(x))==zero
    cocycle_tests=0
    for a,b,d in product(quotient,repeat=3):
        x,y,z=a+(0,),b+(0,),d+(0,)
        assert heisenberg_multiply(heisenberg_multiply(x,y),z)==heisenberg_multiply(x,heisenberg_multiply(y,z))
        cocycle_tests+=1
    u=m.perm(4,[(0,1,2)])
    W=((tuple(range(2)),(u,tuple(range(4)))),)
    x=l.to_perms(W)[0]
    W=((tuple(range(2)),(tuple(range(4)),u)),)
    y=l.to_perms(W)[0]
    assert x in H and y in H and m.compose(x,y)==m.compose(y,x)
    xq,yq=h8_projection(x)+(0,),h8_projection(y)+(0,)
    comm=heisenberg_multiply(heisenberg_multiply(heisenberg_multiply(xq,yq),heisenberg_inverse(xq)),heisenberg_inverse(yq))
    assert comm[:2]==((0,0),0) and comm[2]
    branch_classes=[]
    rows=[]
    tests=0
    for seed in states:
        orbit,_=c.closure(seed)
        assert len(orbit)==16
        value=cubic_h8(seed)
        for state in orbit:
            assert cubic_h8(state)==value
            for move in (-1,1):
                assert cubic_h8(m.braid(state,[move]))==value
                tests+=1
        A=next(g for g in seed if m.cycle_lengths(g)==[2,2,4])
        D=next(g for g in seed if m.cycle_lengths(g)==[1]*4+[2,2])
        CA={m.conjugate(h,A) for h in H}
        CD={m.conjugate(h,D) for h in H}
        assert len(CA)==72 and len(CD)==36
        branch_classes.append((CA,CD))
        for g in (A,D,m.product(seed)):
            lift=h8_canonical(g)
            for h in H:
                if m.compose(g,h)==m.compose(h,g):
                    hq=h8_projection(h)+(0,)
                    assert heisenberg_multiply(lift,hq)==heisenberg_multiply(hq,lift)
        rows.append({"literal_seed":seed,"cubic_value":value,"literal_orbit":16,
                     "A_H_class_size":72,"DT_H_class_size":36,
                     "orbit_sha256":c.digest(orbit)})
    assert branch_classes[0]==branch_classes[1]
    assert sorted(row["cubic_value"] for row in rows)==[1,2]
    return {"actual_group_order":576,"quotient_order":36,
            "Heisenberg_quotient_cover_order":108,"pullback_cover_order":1728,
            "all_cocycle_associativity_tests":cocycle_tests,
            "commuting_base_three_cycles":[x,y],"stem_commutator_central_coordinate":comm[2],
            "same_intrinsic_branch_classes":True,"signed_braid_value_tests":tests,
            "components":rows}


def intrinsic_class_controls():
    old=json.loads((HERE/"polynomial_involution_components_verification.json").read_text())
    domain=next(row for row in old["exact_component_controls"]["complete_passports"]
                if row["degree"]==10 and row["simple_factors"]==2)
    states=[tuple(tuple(g) for g in row["literal_seed"]) for row in domain["components"]]
    G=r.perm_group(states[0],10)
    assert len(G)==28800 and r.perm_group(states[1],10)==G
    classes=[]
    for state in states:
        beta=next(g for g in state if c.index(g)==2)
        orbit={m.conjugate(h,beta) for h in G}
        classes.append(orbit)
    assert sorted(map(len,classes))==[30,100] and not classes[0]&classes[1]
    return {"passport":"(2^5),(2^2 1^6),T,T","same_literal_group_order":28800,
            "different_beta_G_class_sizes":sorted(map(len,classes)),
            "same_group_is_not_same_Nielsen_data":True}


def cyclic_family_controls():
    _,outer=c.v.one.pure_small()
    rows=[]
    for b in (2,3):
        for q0 in outer:
            seed=l.free_seed(q0,b)
            orbit,_=c.closure(seed)
            value=cubic_a7(q0)
            top_order=c.h.group_closure(q0)
            assert len(orbit)==(4116 if top_order==2520 else 1372)
            for state in orbit:
                q=tuple(tuple(g[a]%7 for a in range(7)) for g in state)
                assert cubic_a7(q)==value
            assert c.stabilizer_size(seed)==b
            class_certificates=[]
            if top_order==2520:
                top=set(cover_a7()[0])
                target_dt=m.perm(7,[(0,1),(2,3)])
                I=tuple(range(b))
                Pc=m.perm(b,[tuple(range(b))])
                powers={m.power(Pc,k) for k in range(b)}
                for factor,g in enumerate(seed):
                    q,maps=l.to_wreath((g,),7)[0]
                    if c.index(q):
                        charts=[None]*7
                        for C in c.cycles(q):
                            charts[C[0]]=I
                            a=C[0]
                            for _ in C:
                                z=q[a]
                                value_chart=m.compose(maps[a],charts[a])
                                if charts[z] is not None:
                                    assert charts[z]==value_chart
                                else:
                                    charts[z]=value_chart
                                a=z
                        assert all(F in powers for F in charts)
                        base=l.to_perms(((tuple(range(7)),tuple(m.inverse(F) for F in charts)),))[0]
                        u=next(u for u in sorted(top) if m.conjugate(u,q)==target_dt)
                        diagonal=l.to_perms(((u,(I,)*7),))[0]
                        conjugator=m.compose(diagonal,base)
                        target=l.to_perms(((target_dt,(I,)*7),))[0]
                    else:
                        a,= [a for a,A in enumerate(maps) if A!=I]
                        u=next(u for u in sorted(top) if u[a]==0)
                        conjugator=l.to_perms(((u,(I,)*7),))[0]
                        target=l.to_perms(((tuple(range(7)),(Pc,)+(I,)*6),))[0]
                    assert m.conjugate(conjugator,g)==target
                    class_certificates.append({"factor":factor,"wreath_conjugator":conjugator,
                                               "common_class_representative":target})
            rows.append({"cyclic_fibre":b,"quotient_group_order":top_order,
                         "literal_orbit":len(orbit),"cubic_value":value,
                         "literal_seed":seed,"orbit_sha256":c.digest(orbit),
                         "intrinsic_class_conjugacy_certificates":class_certificates})
    assert sorted(row["cubic_value"] for row in rows if row["quotient_group_order"]==2520)==[1,1,2,2]
    return {"complete_lift_orbits":rows,"literal_states":sum(row["literal_orbit"] for row in rows),
            "scope":"Finite replay b2,b3; all b>=2 uses D6 and L5-L7."}


REFERENCE15=m.perm(15,[(0,1,2,3,4),(5,6,7,8),(9,10,11),(12,13)])
REFERENCE_CYCLES15=tuple(c.cycles(REFERENCE15))
REFERENCE_ROOTS15=tuple((C[j],C[j+1]) for C in REFERENCE_CYCLES15 for j in range(len(C)-1))


def even_class_roots15(g):
    by_length={len(C):C for C in c.cycles(g)}
    assert sorted(by_length)==[1,2,3,4,5]
    conjugator=[None]*15
    for C in REFERENCE_CYCLES15:
        for x,y in zip(C,by_length[len(C)]):
            conjugator[x]=y
    conjugator=tuple(conjugator)
    if c.index(conjugator)%2:
        conjugator=m.compose(conjugator,m.perm(15,[(12,13)]))
    assert c.index(conjugator)%2==0 and m.conjugate(conjugator,REFERENCE15)==g
    return tuple((conjugator[a],conjugator[b]) for a,b in REFERENCE_ROOTS15)


def quadratic15(state):
    roots=[]
    for g in state:
        if c.index(g)==10:
            roots.extend(even_class_roots15(g))
        else:
            C,=c.cycles(g,False)
            assert len(C)==3
            roots.extend(zip(C[:-1],C[1:]))
    assert len(roots)==14
    value=clifford_word(roots)
    raw,sign=c.raw_cycle_lift(15)
    ratio=1 if value==raw else -1
    assert value=={mask:ratio*coef for mask,coef in raw.items()}
    return ratio*sign


def primitive15_controls():
    raw=(REFERENCE15,m.perm(15,[(4,8,11)]),m.perm(15,[(11,13,14)]))
    first=c.v.one.fixed_product(raw)
    second=c.v.one.fixed_product((raw[0],raw[1],m.inverse(raw[2])))
    passport=[m.cycle_lengths(g) for g in first]
    independent=c.character_total(passport)*c.type_placements(passport)
    assert independent==16200
    all_states=set()
    rows=[]
    spin_edge_tests=0
    for seed in (first,second):
        c.validate(seed)
        orbit,_=c.closure(seed)
        assert len(orbit)==8100 and not orbit&all_states
        all_states.update(orbit)
        assert c.block_divisors(seed)==[1,15]
        assert c.stabilizer_size(seed)==1
        value=quadratic15(seed)
        for move in (-1,1,-2,2):
            assert quadratic15(m.braid(seed,[move]))==value
            spin_edge_tests+=1
        assert quadratic15(c.rotated(seed,1))==value
        rows.append({"literal_seed":seed,"literal_orbit":8100,"absolute_orbit":540,
                     "quadratic_value":value,"monodromy":"A15",
                     "group_certificate":"even generators + complete block-divisor test + a literal3cycle; normal3cycle support proof in D7",
                     "all_block_count_divisors":[1,15],"orbit_sha256":c.digest(orbit)})
    assert len(all_states)==independent and sorted(row["quadratic_value"] for row in rows)==[-1,1]
    # The complete centralizer of the strict partition, including both parities.
    cycles=tuple(C for C in REFERENCE_CYCLES15 if len(C)>1)
    reference=clifford_word(REFERENCE_ROOTS15)
    centralizers=Counter()
    for exponents in product(*(range(len(C)) for C in cycles)):
        g=tuple(range(15))
        for C,k in zip(cycles,exponents):
            g=m.compose(g,m.power(m.perm(15,[C]),k))
        assert m.compose(g,REFERENCE15)==m.compose(REFERENCE15,g)
        sign=(-1)**c.index(g)
        transported=clifford_word(tuple((g[a],g[b]) for a,b in REFERENCE_ROOTS15))
        assert transported=={mask:sign*coef for mask,coef in reference.items()}
        centralizers[sign]+=1
    assert dict(centralizers)=={-1:60,1:60}
    assert len(r.perm_group((m.perm(4,[(0,1,2)]),m.perm(4,[(1,2,3)])),4))==12
    assert len(r.perm_group((m.perm(5,[(0,1,2)]),m.perm(5,[(2,3,4)])),5))==60
    return {"complete_literal_domain":independent,"components":rows,
            "quadratic_signed_braid_tests":spin_edge_tests,
            "quadratic_cycle_rotation_tests":2,
            "strict_partition_centralizer_parities":dict(centralizers),
            "overlapping_three_cycle_base_group_orders":[12,60],
            "spin_scope":"Exact values at the two seeds and their listed edges; constancy on all tuples follows from D1/D7's class-equivariant lift proof."}


def negative_controls():
    source=(HERE/"sources/atlas-3a7-45a-generator-1.mtx").read_text()
    values=source.split()
    values[4]="46"
    try:
        parse_mtx(" ".join(values))
    except AssertionError:
        pass
    else:
        raise AssertionError("corrupted permutation was admitted")
    # A branch image of order3 has three order3 preimages in this cover.
    fibres,_,_,_=cover_a7()
    g=m.perm(7,[(0,1,2)])
    assert all(m.power(x,3)==tuple(range(45)) for x in fibres[g])
    try:
        a7_canonical(g,3)
    except AssertionError:
        pass
    else:
        raise AssertionError("noncoprime lift was treated as unique")
    return {"corrupted_permutation_rejected":True,
            "three_order3_lifts_not_mislabelled_unique":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=HERE/"polynomial_lifting_obstructions_verification.json")
    args=parser.parse_args()
    reduction=centralizer_reduction_controls()
    seven=a7_component_controls()
    eight=h8_controls()
    intrinsic=intrinsic_class_controls()
    family=cyclic_family_controls()
    fifteen=primitive15_controls()
    negative=negative_controls()
    paths=["POLYNOMIAL_INTRINSIC_CLASSES_AND_LIFTING.md",
           "POLYNOMIAL_PASSPORT_COMPONENTS.md","POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
           "POLYNOMIAL_INVOLUTION_COMPONENT_CLASSIFICATION.md",
           "verify_polynomial_passport_components.py","verify_polynomial_lift_transport.py",
           "verify_polynomial_return_braids.py","polynomial_lift_transport_verification.json",
           "polynomial_involution_components_verification.json",
           "sources/atlas-3a7-45a-generator-1.mtx","sources/atlas-3a7-45a-generator-2.mtx",
           "sources/atlas-a7-7-generator-1.gap","sources/atlas-a7-7-generator-2.gap",
           "sources/atlas-a7-group-2026.html",
           "sources/roberts-venkatesh-hurwitz-monodromy-2015.pdf",
           "sources/akita-kawasaki-satoh-schur-covers-2024.pdf"]
    receipt={"schema_version":1,"kind":"intrinsic_classes_and_genuine_polynomial_lifting","result":"PASS",
             "python":platform.python_version(),
             "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
             "cover_group_order":7560,"quotient_order":2520,"projection_signed_edges":30240,
             "centralizer_reduction_controls":reduction,"degree7_cubic_components":seven,
             "degree8_cubic_coupling":eight,"intrinsic_class_refinement":intrinsic,
             "cyclic_wreath_cubic_family":family,"degree15_primitive_quadratic_split":fifteen,
             "negative_controls":negative,
             "scope":"Written D1-D7 establish the exact defined invariants and their stated "
                     "separation scope. Complete finite orbit classifications are only those "
                     "explicitly listed; all-degree cyclic inheritance and Sn reduced-multiplier "
                     "vanishing have separate proofs. No arbitrary-passport classification or "
                     "stable-range theorem is inferred from the finite controls.",
             "independent_replay":False,"formal_verification":False,"external_review":False}
    args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"result":"PASS","cubic_cover_order":7560,"degree7_values":{"0":98,"1":147,"2":147},
                      "degree8_values":[1,2],"primitive_A15_components":[8100,8100],
                      "cyclic_lift_states":family["literal_states"],"receipt":str(args.output)},indent=2))


if __name__=="__main__":
    main()
