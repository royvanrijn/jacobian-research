#!/usr/bin/env python3
"""Exact controls for the uniform odd-lobe Hurwitz distance obstruction.

The proof is uniform. This checker exhausts every ordinary signed braid
word through the stated finite depths, tracks which original lobes have
received external conjugation, and checks all resulting leaf deletions.
"""
from __future__ import annotations
import argparse,hashlib,json,platform,time
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def cp(p,q):return tuple(p[q[i]]for i in range(len(p)))
def pm(n,cycles):
    p=list(range(n))
    for c in cycles:
        for i,a in enumerate(c):p[a]=c[(i+1)%len(c)]
    return tuple(p)
def inv(p):
    a=[0]*len(p)
    for i,x in enumerate(p):a[x]=i
    return tuple(a)
def conj(p,q):return cp(cp(p,q),inv(p))
def product(t):
    p=tuple(range(len(t[0])))
    for g in t:p=cp(p,g)
    return p
def cycles(p):
    seen=set();out=[]
    for a in range(len(p)):
        if a in seen:continue
        c=[];b=a
        while b not in seen:seen.add(b);c.append(b);b=p[b]
        out.append(tuple(c))
    return out

def bouquet(t,buffer=False):
    n=4*t+1+2*buffer
    out=[];lobes=[]
    for i in range(t):
        a,b,c,d=range(4*i+1,4*i+5)
        out.extend((pm(n,[(0,a),(b,c)]),pm(n,[(a,b),(c,d)])))
        lobes.append(frozenset((0,a,b,c,d)))
    if buffer:out.append(pm(n,[(0,4*t+1),(1,4*t+2)]))
    return tuple(out),lobes

def leaves(t):
    return [a for a in range(len(t[0]))if sum(g[a]!=a for g in t)==1]

def is_wreath_deletion(t,x):
    """Unique pairing of a full 2r-cycle consists of antipodal pairs."""
    n=len(t[0]);r=(n-1)//2
    owner=[i for i,g in enumerate(t)if g[x]!=x]
    assert len(owner)==1
    P=product(t);C=cycles(P)
    assert len(C)==1
    order=[a for a in C[0]if a!=x]
    label={a:i%r for i,a in enumerate(order)}
    at=owner[0]
    lower=list(t);lower[at]=cp(pm(n,[(x,t[at][x])]),t[at])
    assert all(g[x]==x for g in lower)
    for g in lower:
        for i in range(r):
            a,b=order[i],order[i+r]
            if label[g[a]]!=label[g[b]]:return False
    return True

def point_orbit(gens,start):
    orb={start};q=[start]
    for a in q:
        for g in gens:
            b=g[a]
            if b not in orb:orb.add(b);q.append(b)
    return orb

def group_closure(gens):
    out={tuple(range(len(gens[0])))};q=list(out)
    for a in q:
        for g in gens:
            b=cp(a,g)
            if b not in out:out.add(b);q.append(b)
    return out

def unexposed_controls(t,labels,exposed,lobes):
    unex=[i for i in range(len(lobes))if not(exposed>>i&1)]
    for i in unex:
        p,q=t[labels.index(2*i)],t[labels.index(2*i+1)]
        O=lobes[i]
        assert all(p[a]==q[a]==a for a in range(len(p))if a not in O)
        assert point_orbit((p,q),0)==set(O)
        # Two DT involutions whose product has order five generate D5.
        assert sorted(map(len,cycles(cp(p,q))))==[1]*(len(p)-5)+[5]
    return unex

def move(t,labels,exposed,j,sign,oldcolors):
    a,b=t[j:j+2];la,lb=labels[j:j+2]
    target,actor=(la,lb)if sign>0 else(lb,la)
    target_lobe=target//2 if target<oldcolors else-1
    actor_lobe=actor//2 if actor<oldcolors else-1
    if target_lobe>=0 and target_lobe!=actor_lobe:
        exposed|=1<<target_lobe
    pair=(b,conj(b,a))if sign>0 else(conj(a,b),a)
    return t[:j]+pair+t[j+2:],labels[:j]+(lb,la)+labels[j+2:],exposed

def initial_controls(t,buffer):
    tup,lobes=bouquet(t,buffer);n=len(tup[0]);r=len(tup)
    assert n==2*r+1 and len(cycles(product(tup)))==1
    a1,a2=1,5
    assert cp(cp(tup[0],tup[2]),cp(tup[0],tup[2]))==pm(n,[(0,a1,a2)])
    assert all(sorted(map(len,cycles(g)))==[1]*(n-4)+[2,2]for g in tup)
    for i in range(t):
        A,B=tup[2*i:2*i+2]
        assert len(point_orbit((A,B),0))==5
    L=leaves(tup)
    assert all(not is_wreath_deletion(tup,x)for x in L)
    degree=[sum(g[x]!=x for g in tup)for x in range(n)]
    removable=0
    for g in tup:
        edges=[(a,g[a])for a in range(n)if a<g[a]and g[g[a]]==a]
        if all(min(degree[a],degree[b])==1 for a,b in edges):removable+=1
    assert removable==int(buffer)
    return dict(lobes=t,buffer=buffer,degree=n,DT_factors=r,
        leaf_count=len(L),all_initial_leaf_deletions_symmetric=True,
        removable_DT_count=removable,three_cycle_certificate=[0,a1,a2],
        lower_bound_to_any_wreath_leaf=t-1,
        protected_new_leaf_lower_bound=t if buffer else None)

def complete_word_ball(t,buffer,protected=False):
    tup,lobes=bouquet(t,buffer);r=len(tup);P=product(tup)
    depth=t-1 if protected else t-2
    assert not protected or buffer
    fresh=4*t+1
    count=0;unique=set();protected_leaf_cases=0
    def visit(cur,labels,exposed,d):
        nonlocal count,protected_leaf_cases
        count+=1
        assert exposed.bit_count()<=depth-d
        assert product(cur)==P
        unex=unexposed_controls(cur,labels,exposed,lobes)
        L=leaves(cur)
        if protected:
            assert unex and all(fresh not in lobes[i]for i in unex)
            if fresh in L:
                protected_leaf_cases+=1
                assert not is_wreath_deletion(cur,fresh)
        else:
            assert len(unex)>=2 and 0 not in L
            for x in L:
                assert any(x not in lobes[i]for i in unex)
                assert not is_wreath_deletion(cur,x)
        unique.add(cur)
        if d:
            for j in range(r-1):
                for sign in(-1,1):
                    nxt=move(cur,labels,exposed,j,sign,2*t)
                    visit(*nxt,d-1)
    visit(tup,tuple(range(r)),0,depth)
    assert count==sum((2*(r-1))**k for k in range(depth+1))
    return dict(lobes=t,buffer=buffer,degree=len(P),DT_factors=r,
        protected_new_leaf=protected,all_signed_word_prefixes_through_depth=depth,
        word_prefix_count=count,distinct_tuples=len(unique),
        protected_leaf_endpoint_count=protected_leaf_cases if protected else None,
        every_tested_target_deletion_nonwreath=True,
        unexposed_odd_transitive_subgroup_checks=True)

def negative_controls():
    # The r=2 core is D5 rather than a full-alternating tuple.
    t,_=bouquet(1)
    assert len(group_closure(t))==10
    # The retained wreath normal form really is accepted by the detector.
    r=4;n=2*r+1;x=2*r
    canonical=(pm(n,[(x,2),(0,1)]),)+tuple(
        pm(n,[(0,2*j),(1,2*j+1)])for j in range(r-1,0,-1))
    assert x in leaves(canonical)and is_wreath_deletion(canonical,x)
    # Internal Hurwitz moves may change a generator without exposing its lobe.
    t,lobes=bouquet(3)
    u,labels,exposed=move(t,tuple(range(len(t))),0,0,1,6)
    assert exposed==0 and u!=t
    assert len(unexposed_controls(u,labels,exposed,lobes))==3
    return dict(single_lobe_group_order=10,wreath_detector_positive_control=True,
        internal_braid_preserves_unexposed_lobe=True)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write-receipt",type=Path,
        default=HERE/"hurwitz_leaf_distance_verification.json")
    args=ap.parse_args()
    before=time.monotonic()
    initial=[initial_controls(t,b)for t in(2,3,4,5,12,25)for b in(False,True)]
    balls=[complete_word_ball(t,b)for t in(2,3,4,5)for b in(False,True)]
    protected=[complete_word_ball(t,True,True)for t in(2,3,4)]
    controls=negative_controls()
    paths=[HERE/"UNBOUNDED_HURWITZ_LEAF_NORMALIZATION_DISTANCE.md",
        Path(__file__),HERE/"IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md",
        HERE/"FIXED_LEAF_HURWITZ_OBSTRUCTION.md",
        HERE/"sources/jones-primitive-cycle-1209.5169v1.pdf"]
    out=dict(schema_version=1,result="PASS",python=platform.python_version(),
        checker_sha256=sha(Path(__file__)),source_sha256={
            str(p.relative_to(ROOT)):sha(p)for p in paths},
        initial_family_controls=initial,complete_signed_word_balls=balls,
        protected_leaf_word_balls=protected,negative_controls=controls,
        seconds=round(time.monotonic()-before,6),
        uniform_scope="Written proof, not bounded tests: every r>=4 has a full A_(2r+1) tuple with symmetric initial deletions and any Wr-leaf word has length>=floor(r/2)-1.",
        limitations=["No obstruction to O(r) or longer algorithms.",
            "No all-degree Hurwitz orbit distinction.",
            "No assumption that the wreath locus is reachable in unclassified degrees."],
        formal_verification=False,external_review=False,independent_external_replay=False)
    args.write_receipt.write_text(json.dumps(out,indent=2)+"\n")
    print("PASS: uniform lobe-family controls and complete bounded signed-word checks.")
if __name__=="__main__":main()
