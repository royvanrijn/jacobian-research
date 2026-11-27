#!/usr/bin/env python3
"""Exact controls for two bundled involutions and simple polynomial branches.

Component completeness in all degrees imports Bogatyrev--Gendron, published
Theorem 1.2 and Corollary 2. Finite original closures below are independent of
that import and exhaust exact class-character totals in their stated domains.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import factorial, comb
from pathlib import Path
import argparse
import hashlib
import json
import platform
import verify_polynomial_lift_transport as l

if not __debug__:
    raise RuntimeError("involution component checker requires active assertions")

HERE=Path(__file__).resolve().parent
m=l.m
c=l.c


def divisors(n):
    return tuple(d for d in range(1,n+1) if n%d==0)


def mobius(n):
    sign=1
    p=2
    while p*p<=n:
        if n%p==0:
            n//=p
            sign=-sign
            if n%p==0:
                return 0
        p+=1
    return -sign if n>1 else sign


def involution_type(n,a):
    return [1]*(n-2*a)+[2]*a


@lru_cache(None)
def ordered_count(n,a,b,s):
    """B4's factorial formula, including identity heavy entries."""
    if (n<s+1 or min(a,b)<0 or max(a,b)>n//2 or
            a+b+s!=n-1 or s<1):
        return 0
    value=Fraction(n**(s+1)*factorial(n-a-1)*factorial(n-b-1),
                   factorial(a)*factorial(b)*factorial(n-2*a)*factorial(n-2*b))
    assert value.denominator==1
    return int(value)


@lru_cache(None)
def primitive_ordered_count(n,e,f,s):
    assert e>0 and f>0 and e+f==2*s+2 and e%2==f%2==n%2
    value=sum(mobius(k)*k**(s+1)*
              ordered_count(n//k,(n//k-e)//2,(n//k-f)//2,s)
              for k in divisors(n) if k%2 and n//k>=max(e,f))
    assert value>0
    return value


def partitions(n,s):
    return tuple((e,2*s+2-e) for e in range(1,s+2)
                 if e<=2*s+2-e<=n and e%2==n%2)


def imported_component_number(s,n):
    if n<s+1:
        return 0
    t=min(s,n-s-1)
    return t//2+1 if (n+s)%2 else (t+1)//2


def dihedral_order(d):
    return 1 if d==1 else 2 if d==2 else 2*d


def component_descriptors(N,A,B,s):
    assert N>=2 and s>=1 and N//2>=A>=B>=2 and A+B+s==N-1
    placements=(s+2)*(s+1)//(2 if A==B else 1)
    e,f=N-2*A,N-2*B
    result=[]
    for n in divisors(N):
        d=N//n
        if e>0 and d%2 and n>=f:
            degree_partitions=((e,f),)
        elif e==0 and d%2==0:
            degree_partitions=partitions(n,s)
        else:
            degree_partitions=()
        for e0,f0 in degree_partitions:
            primitive=primitive_ordered_count(n,e0,f0,s)
            if d%2:
                fixed_pattern=d**(s+1)*primitive
            else:
                fixed_pattern=(d//2)*d**s*(1 if e0==f0 else 2)*primitive
            literal=placements*fixed_pattern
            aut=2 if n==2 else 1
            assert literal*aut%N==0
            result.append({"inner_degree":n,"outer_degree":d,
                           "primitive_degree_partition":(e0,f0),
                           "primitive_ordered_inner_count":primitive,
                           "fixed_ordered_original_count":fixed_pattern,
                           "literal_component_size":literal,
                           "absolute_component_size":literal*aut//N,
                           "cyclic_stabilizer_order":aut,
                           "group":"S"+str(n)+" wr Delta_"+str(d),
                           "group_order":factorial(n)**d*dihedral_order(d),
                           "all_block_count_divisors":sorted(set(divisors(d)+(N,)))})
    assert sum(row["literal_component_size"] for row in result)==placements*ordered_count(N,A,B,s)
    return result


def color_words(a,b,s,previous=None,prefix=()):
    """All path colorings with A/B matchings; simple-edge tags in path order."""
    if a+b+s==0:
        yield prefix
        return
    for color,count in ((2,s),(0,a),(1,b)):
        if count and (color==2 or color!=previous):
            yield from color_words(a-(color==0),b-(color==1),s-(color==2),
                                   color,prefix+(color,))


def component_seed(row,s):
    n=row["inner_degree"]
    d=row["outer_degree"]
    e,f=row["primitive_degree_partition"]
    a,b=(n-e)//2,(n-f)//2
    I=tuple(range(n))
    if d%2:
        A=tuple((-x)%d for x in range(d))
        B=m.compose(A,m.perm(d,[tuple(range(d))]))
        fixedA=next(x for x in range(d) if A[x]==x)
        fixedB=next(x for x in range(d) if B[x]==x)
        roles={(0,(fixedA,)):0,(1,(fixedB,)):1}
    else:
        A=tuple((1-x)%d for x in range(d))
        B=tuple((-x)%d for x in range(d))
        fixed=sorted(x for x in range(d) if B[x]==x)
        assert len(fixed)==2
        roles={(1,(fixed[0],)):0,(1,(fixed[1],)):1}
    q=(A,B)+(tuple(range(d)),)*s
    roles.update({(i+2,(0,)):i+2 for i in range(s)})
    slots,_,_,order=l.layout(q)
    for coloring in color_words(a,b,s):
        edges=[[],[]]+[[] for _ in range(s)]
        simple=2
        for at,color in enumerate(coloring):
            target=color if color<2 else simple
            edges[target].append((at,at+1))
            simple+=color==2
        raw=[m.perm(n,part) for part in edges]
        R=tuple(raw[roles[slots[j]]] if slots[j] in roles else I for j in order)
        R=c.v.one.fixed_product(R)
        if l.group_from_transposition(R)[:2]!=(n,1):
            continue
        state=l.to_perms(l.lift(q,R))
        c.validate(state)
        assert l.group_from_transposition(state)[:2]==(n,d)
        return state
    raise AssertionError(("no primitive path seed in this finite control",row,s))


def tuple_label(state,s,known_normal_blocks=None):
    N=len(state[0])
    if known_normal_blocks is None:
        n,d,_=l.group_from_transposition(state)
    else:
        n,d=known_normal_blocks
        assert n*d==N
    heavy=[i for i,g in enumerate(state) if c.index(g)>1]
    assert len(heavy)==2
    A,B=sorted((c.index(state[i]) for i in heavy),reverse=True)
    assert A+B+s==N-1
    if d%2:
        e,f=N-2*A,N-2*B
        assert e>0 and f>0
    else:
        q,R=l.extract(l.to_wreath(state,d))
        small=next(i for i in heavy if c.index(state[i])==B)
        slots,_,_,order=l.layout(q)
        by_slot={slots[j]:h for j,h in zip(order,R)}
        values=[by_slot[(small,C)] for C in c.cycles(q[small]) if len(C)==1]
        assert len(values)==2
        assert all(all(ell<=2 for ell in m.cycle_lengths(h)) for h in values)
        e,f=sorted(n-2*c.index(h) for h in values)
        assert e>0 and f>0
    assert (e,f) in partitions(n,s)
    return n,d,e,f


def exact_component_controls():
    cases=((6,2,2,1),(7,3,2,1),(8,3,3,1),(8,4,2,1),
           (9,4,3,1),(10,5,2,2),(10,4,4,1),(12,6,3,2),
           (18,8,8,1),(14,7,4,2))
    results=[]
    total=0
    for N,A,B,s in cases:
        rows=component_descriptors(N,A,B,s)
        expected_passport=[involution_type(N,A),involution_type(N,B)]+[involution_type(N,1)]*s
        independent=c.character_total(expected_passport)*c.type_placements(expected_passport)
        assert independent==sum(row["literal_component_size"] for row in rows)
        found=set()
        controls=[]
        for row in rows:
            seed=component_seed(row,s)
            orbit,_=c.closure(seed)
            assert len(orbit)==row["literal_component_size"] and not found&orbit
            found.update(orbit)
            n,d,edges=l.group_from_transposition(seed)
            assert (n,d)==(row["inner_degree"],row["outer_degree"])
            assert edges==d*comb(n,2)
            q=tuple(part[0] for part in l.to_wreath(seed,d))
            quotient_order=c.h.group_closure(q)
            assert quotient_order==dihedral_order(d)
            expected_label=(n,d)+tuple(row["primitive_degree_partition"])
            assert tuple_label(seed,s)==expected_label
            labels_checked=0
            patterns=Counter()
            for state in orbit:
                pattern=tuple(tuple(m.cycle_lengths(g)) for g in state)
                patterns[pattern]+=1
                assert tuple_label(state,s,(n,d))==expected_label
                labels_checked+=1
            assert all(count==row["fixed_ordered_original_count"] for count in patterns.values())
            assert len(patterns)==c.type_placements(expected_passport)
            assert c.block_divisors(seed)==row["all_block_count_divisors"]
            assert c.stabilizer_size(seed)==row["cyclic_stabilizer_order"]
            group=c.h.group_closure(seed) if row["group_order"]<=70000 else None
            if group is not None:
                assert group==row["group_order"]
            controls.append(dict(row,literal_seed=seed,orbit_sha256=c.digest(orbit),
                                 all_literal_labels_checked=labels_checked,
                                 normal_transposition_support_edges=edges,
                                 exact_quotient_group_order=quotient_order,
                                 exact_generated_group_closure=group,
                                 group_basis="complete normal-transposition graph + exact quotient closure; B2"))
        assert len(found)==independent
        total+=independent
        results.append({"degree":N,"bundled_indices":[A,B],"simple_factors":s,
                        "ordered_character_total":c.character_total(expected_passport),
                        "full_type_placements":c.type_placements(expected_passport),
                        "complete_literal_domain":independent,"components":controls})
        print(json.dumps({"degree":N,"bundled_indices":[A,B],"simple_factors":s,
                          "literal_components":[row["literal_component_size"] for row in rows]}),flush=True)
    return {"complete_passports":results,"complete_literal_states":total}


def factorial_and_partition_controls():
    count=0
    identities=0
    for n in range(2,21):
        for a in range(n//2+1):
            for b in range(n//2+1):
                s=n-1-a-b
                if s<1:
                    continue
                passports=[involution_type(n,a),involution_type(n,b)]+[involution_type(n,1)]*s
                assert ordered_count(n,a,b,s)==c.character_total(passports)
                count+=1
                identities+=a==0 or b==0
    for s in range(1,15):
        for n in range(1,41):
            assert len(partitions(n,s))==imported_component_number(s,n)
    # This larger finite arithmetic-only range checks the proposed divisor
    # accounting, not braid transitivity in these passports.
    sums=0
    for N in range(6,41):
        for A in range(2,N//2+1):
            for B in range(2,A+1):
                s=N-1-A-B
                if s>0:
                    rows=component_descriptors(N,A,B,s)
                    assert len(rows)==(
                        sum(N//n%2==1 and n>=N-2*B for n in divisors(N))
                        if N>2*A else
                        sum(imported_component_number(s,n) for n in divisors(N) if N//n%2==0))
                    sums+=1
    return {"factorial_vs_hook_character_domains":count,"domains_with_identity_heavy_entry":identities,
            "factorial_count_max_degree":20,"primitive_partition_domains":14*40,
            "partition_max_simple_count":14,"partition_max_inner_degree":40,
            "divisor_sum_domains":sums,"divisor_sum_max_degree":40,
            "scope":"Exact arithmetic identities in these finite domains; no additional braid orbit assertions."}


def prime_family_controls():
    results=[]
    for p in (5,7,11,13,17,19):
        rows=component_descriptors(2*p,p,p-3,2)
        assert len(rows)==2 and {tuple(row["primitive_degree_partition"]) for row in rows}=={(1,5),(3,3)}
        counts={tuple(row["primitive_degree_partition"]):row["literal_component_size"] for row in rows}
        assert counts[(1,5)]==p**3*(p*p-9)*(p*p-1)//20
        assert counts[(3,3)]==p**3*(p*p-1)**2//12
        assert all(row["inner_degree"]==p and row["outer_degree"]==2 and
                   row["all_block_count_divisors"]==[1,2,2*p] and
                   row["cyclic_stabilizer_order"]==1 for row in rows)
        results.append({"prime":p,"literal_component_sizes":[counts[(1,5)],counts[(3,3)]],
                        "group_order":2*factorial(p)**2})
    return {"formula_only_prime_controls":results,
            "direct_complete_orbits_for_p5_p7":"Recorded separately among exact_component_controls",
            "scope":"The unbounded prime-family conclusion is B6's proof, not extrapolation from this list."}


def expect_rejected(function):
    try:
        function()
    except AssertionError:
        return True
    raise AssertionError("invalid passport was admitted")


def negative_controls():
    expect_rejected(lambda:component_descriptors(7,3,3,0))
    expect_rejected(lambda:component_descriptors(6,3,2,1))
    expect_rejected(lambda:component_descriptors(8,3,1,3))
    expect_rejected(lambda:primitive_ordered_count(5,0,6,2))
    return {"no_simple_branch_rejected":True,"wrong_index_budget_rejected":True,
            "outside_A_B_at_least_two_scope_rejected":True,
            "nonpositive_primitive_partition_rejected":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=HERE/"polynomial_involution_components_verification.json")
    args=parser.parse_args()
    arithmetic=factorial_and_partition_controls()
    components=exact_component_controls()
    prime=prime_family_controls()
    negatives=negative_controls()
    paths=["POLYNOMIAL_INVOLUTION_COMPONENT_CLASSIFICATION.md",
           "POLYNOMIAL_PASSPORT_COMPONENTS.md","POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
           "verify_polynomial_passport_components.py","verify_polynomial_lift_transport.py",
           "sources/bogatyrev-gendron-pell-abel-2025.pdf"]
    receipt={"schema_version":1,"kind":"two_bundled_involutions_and_simple_branches","result":"PASS",
             "python":platform.python_version(),
             "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
             "arithmetic_controls":arithmetic,"exact_component_controls":components,
             "same_group_prime_family":prime,"negative_controls":negatives,
             "scope":"All-degree component sufficiency is imported explicitly from the published "
                     "Bogatyrev--Gendron Theorem1.2/Corollary2. Monodromy, blocks, literal label "
                     "extraction, factorial/Mobius orbit sizes and prime-family conclusions have "
                     "written derivations B2-B6. Only the listed ten passports are directly "
                     "exhausted by ordinary braid closures and independent hook totals.",
             "independent_replay":False,"formal_verification":False,"external_review":False}
    args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"result":"PASS","complete_passports":10,
                      "complete_literal_states":components["complete_literal_states"],
                      "receipt":str(args.output)},indent=2))


if __name__=="__main__":
    main()
