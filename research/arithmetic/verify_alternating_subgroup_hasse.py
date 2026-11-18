#!/usr/bin/env python3
"""Exact A10 subset/equipartition covers and tame decomposition-group controls.

Run: python3 research/arithmetic/verify_alternating_subgroup_hasse.py
Optional --write-receipt PATH. Uses only the Python standard library.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
from itertools import combinations, permutations
import json
from math import comb
from pathlib import Path
import platform

N=10
IDENTITY=tuple(range(N))
ALL=(1<<N)-1

def require(ok,message):
    if not ok:
        raise AssertionError(message)

def compose(a,b):
    return tuple(a[b[j]] for j in range(N))

def permutation(*cycles):
    out=list(IDENTITY)
    for cycle in cycles:
        for a,b in zip(cycle,cycle[1:]+cycle[:1]):
            out[a]=b
    return tuple(out)

def parity(value):
    return sum(value[i]>value[j] for i in range(N) for j in range(i+1,N))%2

def partitions(n,minimum=1):
    if n==0:
        yield ()
    for first in range(minimum,n+1):
        for rest in partitions(n-first,first):
            yield (first,)+rest

def representative(parts):
    cycles=[];start=0
    for length in parts:
        cycles.append(tuple(range(start,start+length)))
        start+=length
    return permutation(*cycles)

def orbit_masks(generators):
    unseen=set(range(N));out=[]
    while unseen:
        first=min(unseen);seen={first};queue=[first]
        for point in queue:
            for generator in generators:
                image=generator[point]
                if image not in seen:
                    seen.add(image);queue.append(image)
        unseen-=seen
        out.append(sum(1<<point for point in seen))
    return out

def subset_sizes(generators):
    sizes={0}
    for mask in orbit_masks(generators):
        length=mask.bit_count()
        sizes|={length+old for old in tuple(sizes)}
    return sizes

EQUIPARTITIONS=[sum(1<<j for j in subset) for subset in combinations(range(N),5) if 0 in subset]
def image_mask(value,mask):
    return sum(1<<value[j] for j in range(N) if (mask>>j)&1)

def equipartition_fixed(generators):
    for mask in EQUIPARTITIONS:
        if all(image_mask(g,mask) in (mask,ALL^mask) for g in generators):
            return True
    return False

def hits(generators,include_h2=False):
    sizes=subset_sizes(generators)
    result=[]
    for k in ((1,2,3) if include_h2 else (1,3)):
        if k in sizes:
            result.append("H"+str(k))
    if equipartition_fixed(generators):
        result.append("W5")
    return result

def perfect_matchings(points):
    if not points:
        yield IDENTITY
        return
    first=points[0]
    for index in range(1,len(points)):
        second=points[index]
        remaining=points[1:index]+points[index+1:]
        for rest in perfect_matchings(remaining):
            value=list(rest)
            value[first]=second
            value[second]=first
            yield tuple(value)

MATCHINGS=list(perfect_matchings(tuple(range(N))))
def commutes(a,b):
    return all(a[b[j]]==b[a[j]] for j in range(N))

def matching_fixed(generators,matchings=MATCHINGS):
    return any(all(commutes(g,m) for g in generators) for m in matchings)

def matching_cover_hits(generators,matchings=MATCHINGS):
    sizes=subset_sizes(generators)
    out=["H"+str(k) for k in (1,3) if k in sizes]
    if matching_fixed(generators,matchings):
        out.append("W2")
    return out

def generated_group(generators):
    seen={IDENTITY};queue=[IDENTITY]
    for value in queue:
        for generator in generators:
            next_value=compose(value,generator)
            if next_value not in seen:
                seen.add(next_value);queue.append(next_value)
    return seen

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()

    H_coefficients={10:9,8:-864,6:28161,4:-350150,2:1127500}
    s_star=sum(H_coefficients.values())+250000
    derivative_at_1=sum(j*c for j,c in H_coefficients.items())-s_star
    require(s_star==1054656 and derivative_at_1==-38112,
            "exact Hensel basepoint and derivative")
    require(38112==2**5*3*397,
            "Hensel derivative prime valuations")

    require(len(MATCHINGS)==945 and len(set(MATCHINGS))==945,
            "complete perfect-matching action")
    types=[parts for parts in partitions(10) if len(parts)%2==0]
    witnesses={}
    cover_rows=[]
    for parts in types:
        g=representative(parts)
        require(parity(g)==0,"even cycle type")
        require(bool(hits([g])),"three-class element cover")
        require(bool(matching_cover_hits([g])),"matching three-class element cover")
        sums=subset_sizes([g])
        valid=sorted(k for k in range(1,6) if k in sums)
        require(bool(valid),"natural-subset cover")
        cover_rows.append({"cycle_type":list(parts),"natural_subset_sizes":valid,
                           "three_cover_hits":hits([g]),
                           "matching_cover_hits":matching_cover_hits([g])})
    for k in range(1,6):
        parts=tuple(sorted((k,10-k)))
        g=representative(parts)
        valid=[j for j in range(1,6) if j in subset_sizes([g])]
        require(valid==[k],"each subset class forced")
        witnesses[str(k)]=list(parts)

    # Exact finite reduction: C_{S10}(d) is C_{S4}(d) x S6.
    # First factor has 8 elements; impose even total parity.
    d=permutation((0,1),(2,3))
    central4=[]
    for p4 in permutations(range(4)):
        h=tuple(p4)+tuple(range(4,10))
        if compose(h,d)==compose(d,h):
            central4.append(p4)
    require(len(central4)==8,"complete four-letter centralizer")
    centralizer=[]
    for p4 in central4:
        for p6 in permutations(range(4,10)):
            h=tuple(p4)+tuple(p6)
            if parity(h)==0:
                require(compose(h,d)==compose(d,h),"centralizer relation")
                centralizer.append(h)
    require(len(centralizer)==2880 and len(set(centralizer))==2880,
            "complete even centralizer")

    d_matchings=[m for m in MATCHINGS if commutes(d,m)]
    hit_counts=Counter()
    matching_hit_counts=Counter()
    three_failures=0
    for h in centralizer:
        found=hits([d,h],include_h2=True)
        require(bool(found),"four classes cover all DT-cyclic-quotient groups")
        hit_counts[",".join(found)]+=1
        matching_found=matching_cover_hits([d,h],d_matchings)
        require(bool(matching_found),"three matching classes cover every DT local group")
        matching_hit_counts[",".join(matching_found)]+=1
        if not hits([d,h]):
            three_failures+=1
            require(2 in subset_sizes([d,h]),"H2 repairs each missing D")

    bad_h=permutation((0,1),(4,5,6,7,8,9))
    bad_group=generated_group([d,bad_h])
    require(parity(bad_h)==0 and len(bad_group)==12,"C2 x C6 example")
    require(compose(d,bad_h)==compose(bad_h,d),"abelian decomposition group")
    require(sorted(mask.bit_count() for mask in orbit_masks([d,bad_h]))==[2,2,6],
            "counterexample orbit sizes")
    require(not hits([d,bad_h]),"three-cover ramification counterexample")
    require(hits([d,bad_h],True)==["H2"],"new component is essential")
    require(matching_cover_hits([d,bad_h])==["W2"],"matching cover repairs counterexample")

    # Removing each of H1,H3,W5 fails already on a cyclic group.
    essential={
        "H1":permutation(tuple(range(9))),
        "H3":permutation(tuple(range(7)),tuple(range(7,10))),
        "W5":permutation(tuple(range(5)),tuple(range(5,10)))
    }
    for expected,g in essential.items():
        require(hits([g],True)==[expected],"essential element witness")
    require(sum(comb(10,k) for k in range(1,6))==637,"five-subset algebra rank")
    require(10+120+comb(10,5)//2==256,"three-cover rank")
    require(10+45+120+comb(10,5)//2==301,"four-cover rank")
    require(10+120+len(MATCHINGS)==1075,"three matching-cover rank")

    uniform_sanity=[]
    for n in range(8,31,2):
        b=n//2
        element_cases=local_fixed_cycle_cases=0
        for parts in partitions(n):
            if len(parts)%2:
                continue
            sizes={0}
            for length in parts:
                sizes|={old+length for old in tuple(sizes)}
            odd_subset=any(k in sizes for k in range(1,b,2))
            multiplicities=Counter(parts)
            matching=all(multiplicities[length]%2==0 for length in multiplicities if length%2)
            require(odd_subset or matching,"uniform cyclic-cover sanity")
            element_cases+=1
        for parts in partitions(n-4):
            sizes={0}
            for length in parts:
                sizes|={old+length for old in tuple(sizes)}
            odd_subset=any(k in sizes for k in range(1,b,2))
            all_even=all(length%2==0 for length in parts)
            require(odd_subset or all_even,"uniform DT-local-cover sanity")
            local_fixed_cycle_cases+=1
        uniform_sanity.append({"degree":n,"odd_subset_classes":list(range(1,b,2)),
            "cyclic_cycle_type_cases":element_cases,
            "DT_fixed_letter_cycle_type_cases":local_fixed_cycle_cases})

    report={
        "schema_version":1,"kind":"complete_finite_group_certificate",
        "command":"python3 research/arithmetic/verify_alternating_subgroup_hasse.py",
        "python":platform.python_version(),"arithmetic":"integer permutations and finite subsets",
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "hensel_basepoint":{"z":1,"s":s_star,"derivative":derivative_at_1,
            "derivative_factorization":{"2":5,"3":1,"397":1}},
        "even_cycle_types":len(types),"cycle_type_cover":cover_rows,
        "natural_subset_minimum":5,"forced_subset_size_witnesses":witnesses,
        "natural_subset_rank":637,"three_element_cover":["H1","H3","W5"],
        "three_element_cover_rank":256,
        "tame_DT_cover":["H1","H2","H3","W5"],"tame_DT_cover_rank":301,
        "optimal_class_count_DT_cover":["H1","H3","W2"],
        "optimal_class_count_DT_cover_rank":1075,
        "perfect_matching_action_degree":len(MATCHINGS),
        "DT_invariant_perfect_matchings":len(d_matchings),
        "DT_matching_cover_hit_counts":dict(sorted(matching_hit_counts.items())),
        "DT_centralizer_even_elements":len(centralizer),
        "DT_generator_pair_hit_counts":dict(sorted(hit_counts.items())),
        "three_cover_failures_in_DT_centralizer":three_failures,
        "ramified_counterexample":{"inertia":list(d),"frobenius":list(bad_h),
            "group_order":len(bad_group),"orbit_sizes":[2,2,6]},
        "negative_controls":["deleting H1 loses type(9,1)","deleting H3 loses type(7,3)",
            "deleting W5 loses type(5,5)","deleting H2 loses D=C2xC6",
            "H1,H3,W2 covers the ramified counterexample"],
        "scope":"Complete cycle-type and centralizer enumeration for A10. Infinite arithmetic existence additionally uses the written specialization-inertia, local lifting and integral Hilbert argument. Optimality of three subgroup classes for either element or DT-local coverage imports Bubboloni Theorem3.9.",
        "uniform_partition_sanity_degrees8_through30":uniform_sanity,
        "independent_replay":False,"formal_verification":False,"external_review":False,
        "result":"PASS"}
    text=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(text)
    print(text,end="")

if __name__=="__main__":
    main()
