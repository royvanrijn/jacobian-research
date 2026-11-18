#!/usr/bin/env python3
"""Exact matching-walk certificates for all-real paired Laurent covers.

The note proves the uniform classification and realization. This checker
verifies the constructor/spectra in finite controls, exact chamber counts,
and a separate complete DT^4 enumeration in degree eight.
Only Python's standard library is required.
"""
from __future__ import annotations
import argparse
from collections import defaultdict, Counter
from itertools import combinations
import hashlib
import json
from pathlib import Path
import platform

def canonical(edge):
    return tuple(sorted(edge))

def path_edges(n):
    q=n//2-1
    return ([(2*i+1,2*i+2) for i in range(q)],
            [(2*i+2,2*i+3) for i in range(q)])

def compose(p,q):
    return tuple(p[q[i]] for i in range(len(p)))

def inverse(p):
    out=[0]*len(p)
    for i,j in enumerate(p):out[j]=i
    return tuple(out)

def perm(n,edges):
    out=list(range(n))
    for a,b in edges:out[a],out[b]=out[b],out[a]
    return tuple(out)

def reflect_edge(n,edge):
    return canonical(0 if j==0 else n-j for j in edge)

def validate(n,steps,reversible=False):
    assert n>=8 and n%2==0 and len(steps)==n//2
    A,B=path_edges(n)
    current=set(A);states=[set(current)];product=tuple(range(n))
    used=Counter()
    for step in steps:
        assert len(step)==2
        e,f=map(canonical,step)
        assert len(set(e+f))==4
        occupied=set(v for edge in current for v in edge)
        for edge in (e,f):
            assert edge in current or not (set(edge)&occupied)
        current.symmetric_difference_update((e,f))
        assert len(set(v for edge in current for v in edge))==2*len(current)
        used.update((e,f))
        states.append(set(current))
        product=compose(product,perm(n,(e,f)))
    assert current==set(B)
    assert product==compose(perm(n,A),perm(n,B))
    outside={e:c for e,c in used.items() if e not in set(A+B)}
    assert len(outside)==1
    (aux,multiplicity),=outside.items()
    assert 0 in aux and multiplicity==2
    assert all(used[e]==1 for e in A+B)
    seen={0}
    while True:
        enlarged=seen|{v for e in used if set(e)&seen for v in e}
        if enlarged==seen:break
        seen=enlarged
    assert seen==set(range(n))
    if reversible:
        assert aux==(0,n//2)
        assert all({reflect_edge(n,e) for e in states[i]}==states[-1-i]
                   for i in range(len(states)))
    counts=[n-2*len(state) for state in states]
    assert all(r%4==2 for r in counts)
    return states,counts

def symmetric_peak_six(n):
    b=n//2;q=b-1;h=b//2
    A,B=path_edges(n)
    targets=set(range(h-2 if b%2==0 else h-1,q))
    center=(b-1)//2
    prefix=[(A[center],A[q-1])]
    targets.difference_update((center,q-1))
    last=max(targets)
    prefix.append(((0,b),A[last]))
    targets.remove(last)
    for offset,index in enumerate(sorted(targets,reverse=True)):
        prefix.append((A[index],B[q-1-offset]))
    assert len(prefix)==h
    steps=list(prefix)
    if b%2:
        edge=A[h-2]
        steps.append((edge,reflect_edge(n,edge)))
    steps.extend(tuple(reflect_edge(n,e) for e in step)
                 for step in reversed(prefix))
    _,counts=validate(n,steps,True)
    assert max(counts)==6
    return steps

def extend_four(n,steps):
    def shift(edge):return canonical(0 if j==0 else j+2 for j in edge)
    bigger=[((1,2),(n+1,n+2))]
    bigger.extend(tuple(shift(e) for e in step) for step in steps)
    bigger.append(((2,3),(n+2,n+3)))
    _,old=validate(n,steps,True)
    _,new=validate(n+4,bigger,True)
    assert max(new)==max(old)+4
    return bigger

def general_peak(n,peak):
    b=n//2;q=b-1;h=(peak-2)//4
    A,B=path_edges(n)
    steps=[(A[2*i],A[2*i+1])for i in range(h)]
    if 2*h==q:
        steps.extend([((0,1),B[0]),((0,1),B[1])])
        remaining=B[2:]
    else:
        steps.append(((0,1),A[2*h]))
        for i in range(q-2*h-1):
            steps.append((B[i],A[2*h+1+i]))
        first=q-2*h-1
        steps.append(((0,1),B[first]))
        remaining=B[first+1:]
    steps.extend((remaining[2*i],remaining[2*i+1])
                 for i in range(len(remaining)//2))
    _,counts=validate(n,steps)
    assert max(counts)==peak
    return steps

def advance(n,j,states):
    A,B=path_edges(n);q=len(A);aux=(0,j)
    amask=sum(1<<i for i,e in enumerate(A)if j in e)
    bmask=sum(1<<i for i,e in enumerate(B)if j in e)
    deps=[sum(1<<i for i,e in enumerate(A)if set(e)&set(f))for f in B]
    nxt={}
    for (aa,bb,phase),dist in states.items():
        events=[]
        for i,e in enumerate(A):
            if not aa>>i&1:events.append((0,i,e))
        for i,e in enumerate(B):
            if not bb>>i&1 and aa&deps[i]==deps[i] and (not bmask>>i&1 or phase==2):
                events.append((1,i,e))
        if phase==0 and aa&amask==amask:events.append((2,0,aux))
        if phase==1:events.append((3,0,aux))
        for e,f in combinations(events,2):
            if set(e[2])&set(f[2]):continue
            a2,b2,ph=aa,bb,phase
            for typ,i,edge in(e,f):
                if typ==0:a2|=1<<i
                elif typ==1:b2|=1<<i
                elif typ==2:ph=1
                else:ph=2
            out=nxt.setdefault((a2,b2,ph),defaultdict(int))
            real=n-2*(q-a2.bit_count()+b2.bit_count()+(ph==1))
            for peak,count in dist.items():out[max(peak,real)]+=count
    return nxt

def count_all(n):
    b=n//2;q=b-1;full=(1<<q)-1;answer=defaultdict(int)
    for j in range(1,n):
        states={(0,0,0):{2:1}}
        for _ in range(b):states=advance(n,j,states)
        for peak,count in states.get((full,full,2),{}).items():answer[peak]+=count
    return dict(answer)

def count_symmetric(n):
    b=n//2;q=b-1;full=(1<<q)-1
    states={(0,0,0):{2:1}}
    for _ in range(b//2):states=advance(n,b,states)
    A,_=path_edges(n)
    def rev(mask):return int(f"{mask:0{q}b}"[::-1],2)
    answer=defaultdict(int)
    for (aa,bb,phase),dist in states.items():
        if phase!=1:continue
        delta=(full^aa)^rev(bb)
        if b%2==0:
            if delta:continue
        else:
            if delta.bit_count()!=1:continue
            edge=A[delta.bit_length()-1]
            if set(edge)&set(reflect_edge(n,edge)):continue
        for peak,count in dist.items():answer[peak]+=count
    return dict(answer)

def degree_eight_independent():
    n=8;A,B=path_edges(n);ap=perm(n,A);bp=perm(n,B);P=compose(ap,bp)
    C=[perm(n,es) for a,b,c,d in combinations(range(n),4)
       for es in [((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))]]
    assert len(set(C))==210
    pairs=defaultdict(list)
    for a in C:
        for b in C:pairs[compose(a,b)].append((a,b))
    reflections=[];power=tuple(range(n))
    for _ in range(7):
        reflections.append(compose(ap,power));power=compose(power,P)
    total=0;labeled=0;canonical_tuples=[]
    for ab,left in pairs.items():
        for a,b in left:
            for c,d in pairs.get(compose(inverse(ab),P),()):
                total+=1;tup=(a,b,c,d)
                for c0 in reflections:
                    current=c0;good=True
                    for sig in tup:
                        moved=[v for v in range(n)if sig[v]!=v]
                        if any(current[v] not in(v,sig[v])for v in moved):
                            good=False;break
                        current=compose(current,sig)
                        if compose(current,current)!=tuple(range(n)):
                            good=False;break
                    if good and any(sig[0]!=0 for sig in tup):
                        labeled+=1
                        if c0==ap:canonical_tuples.append(tup)
    assert total==93296 and labeled==504 and len(canonical_tuples)==72
    symmetric=0;orders=Counter();r=tuple([0]+list(range(n-1,0,-1)))
    for tup in canonical_tuples:
        seen={tuple(range(n))};queue=list(seen)
        for p in queue:
            for sig in tup:
                candidate=compose(p,sig)
                if candidate not in seen:seen.add(candidate);queue.append(candidate)
        orders[len(seen)]+=1
        cs=[ap]
        for sig in tup:cs.append(compose(cs[-1],sig))
        if all(compose(compose(r,cs[i]),r)==cs[-1-i]for i in range(5)):
            symmetric+=1
    assert orders=={20160:72} and symmetric==4
    return {"raw_DT_quadruples":total,"labeled_all_real_transitive":labeled,
            "absolute_all_real":72,"all_generated_group_orders":dict(orders),
            "reversible":symmetric}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    parser.add_argument("--include-degree26-all",action="store_true",
                        help="Also run the approximately one-minute full degree26 DP")
    args=parser.parse_args()
    controls=[]
    for n in range(8,62,2):
        max_general=n-2 if n%4==0 else n
        max_symmetric=n-2 if n%4==0 else n-4
        for peak in range(6,max_general+1,4):general_peak(n,peak)
        for peak in range(6,max_symmetric+1,4):
            extensions=(peak-6)//4;base=n-4*extensions
            steps=symmetric_peak_six(base)
            for _ in range(extensions):steps=extend_four(base,steps);base+=4
            _,counts=validate(n,steps,True)
            assert max(counts)==peak and set(counts)==set(range(2,peak+1,4))
        controls.append(n)
    all_counts={str(n):count_all(n) for n in (8,10,12)}
    assert all_counts=={"8":{6:72},"10":{10:1008,6:840},"12":{10:55728,6:7848}}
    sym_counts={str(n):count_symmetric(n) for n in (8,10,12,26)}
    assert sym_counts=={
        "8":{6:4},"10":{6:8},"12":{10:60,6:36},
        "26":{22:47628000,18:128943360,14:83956320,10:9108696,6:32296}}
    if args.include_degree26_all:
        all_counts["26"]=count_all(26)
        assert all_counts["26"]=={
            26:36074460723840000,22:700096359090240000,18:1016860269086323200,
            14:159268810077705600,10:1067699019638736,6:9220226200}
    independent=degree_eight_independent()
    bad=symmetric_peak_six(8)
    bad[0]=((0,1),(1,2))
    try:validate(8,bad,True)
    except AssertionError:pass
    else:raise AssertionError("Overlapping-transposition control was accepted")
    report={
        "schema_version":1,"kind":"exact_real_matching_walk_certificate",
        "python":platform.python_version(),
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "constructor_control_degrees":controls,
        "all_real_peak_distributions":all_counts,
        "all_real_totals":{n:sum(d.values())for n,d in all_counts.items()},
        "odd_real_peak_distributions":sym_counts,
        "odd_real_totals":{n:sum(d.values())for n,d in sym_counts.items()},
        "independent_degree8_DT_enumeration":independent,
        "all_real_branch_values_context":"One fixed ordered set of b distinct real values; target fixed, quotient source isomorphism",
        "odd_context":"One fixed negation-symmetric ordered branch set; reversible walks and odd Laurent functions",
        "source_arithmetic":"Exact permutations, bit masks and unbounded integer DP; no numerical cover recovery",
        "limits":"Uniform real/real-algebraic realization and all-degree spectra use the written proofs. No Q coefficients, degree26 rational model, uniform rational parametrization or complex Hurwitz connectivity is certified.",
        "independent_replay":False,"external_review":False,"formal_verification":False,
        "result":"PASS"
    }
    output=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(output)
    print(output,end="")
if __name__=="__main__":main()

