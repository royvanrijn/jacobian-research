#!/usr/bin/env python3
"""Exact constructive controls for imprimitive leaf-deletion normal forms."""
import importlib.util
from itertools import combinations, permutations, product
from collections import deque
import argparse, hashlib, heapq, json, platform, random
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("fixed",HERE/"verify_fixed_leaf_hurwitz_obstruction.py")
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
m=f.m;leaf=f.leaf

def full_deletion(t,x):
    at=f.owner(t,x);assert at is not None
    small=list(t);small[at]=m.compose(m.perm(len(t[0]),[(x,t[at][x])]),t[at])
    assert all(g[x]==x for g in small)
    return tuple(small),at

def lower_class(t,x):
    d,at=full_deletion(t,x);e=next(iter(leaf.edges(d[at])))
    orbit={tuple(sorted(e))};todo=list(orbit)
    for a,b in todo:
        for g in d:
            edge=tuple(sorted((g[a],g[b])))
            if edge not in orbit:orbit.add(edge);todo.append(edge)
    n=len(t[0]);adj=[set()for _ in range(n)]
    for a,b in orbit:adj[a].add(b);adj[b].add(a)
    unseen=set(range(n))-{x};blocks=[]
    while unseen:
        todo=[min(unseen)];unseen.remove(todo[0]);C=set(todo)
        for a in todo:
            for b in adj[a]:
                if b not in C:C.add(b);todo.append(b);unseen.remove(b)
        blocks.append(tuple(sorted(C)))
    blocks.sort()
    if len(blocks)==1:return "S",blocks
    assert all(len(C)==2 for C in blocks)
    return "W",blocks

def project(g,blocks):
    labels={a:i for i,C in enumerate(blocks)for a in C};p=[]
    for C in blocks:
        image={labels[g[a]]for a in C};assert len(image)==1
        p.append(next(iter(image)))
    return tuple(p)

def star_tuple(P,c):
    seq=[];a=P[c]
    while a!=c:seq.append(a);a=P[a]
    assert len(seq)==len(P)-1
    return tuple(m.perm(len(P),[(c,a)])for a in seq[::-1])

def delete_fixed(p,l):
    old=[a for a in range(len(p))if a!=l];idx={a:i for i,a in enumerate(old)}
    assert p[l]==l
    return tuple(idx[p[a]]for a in old)

def star_word(t,c):
    P=m.product(t);S=star_tuple(P,c)
    if t==S:return []
    assert len(t)>=2
    deg=[0]*len(P)
    for g in t:
        es=leaf.edges(g);assert len(es)==1
        for a,b in es:deg[a]+=1;deg[b]+=1
    l=next(a for a,d in enumerate(deg)if d==1 and a!=c)
    at=next(i for i,g in enumerate(t)if g[l]!=l)
    target_at=next(i for i,g in enumerate(S)if g[l]!=l)
    bring=list(range(at+1,len(t)));bringS=list(range(target_at+1,len(S)))
    moved=m.braid(t,bring);target=m.braid(S,bringS)
    assert moved[-1]==target[-1]
    old=[a for a in range(len(P))if a!=l];idx={a:i for i,a in enumerate(old)}
    small=tuple(delete_fixed(g,l)for g in moved[:-1])
    targetsmall=tuple(delete_fixed(g,l)for g in target[:-1])
    assert targetsmall==star_tuple(m.product(small),idx[c])
    middle=star_word(small,idx[c])
    word=bring+middle+[-v for v in bringS[::-1]]
    assert m.braid(t,word)==S
    return word

def apply(t,w,x):
    for code in w:
        t=m.hurwitz(t,abs(code)-1,1 if code>0 else-1)
        assert f.owner(t,x)is not None
    return t

def star_upper(t,x,blocks,c):
    ps=tuple(project(g,blocks)for g in t[1:])
    w=star_word(ps,c)
    upper=[(abs(a)+1)*(1 if a>0 else-1)for a in w]
    out=apply(t,upper,x)
    assert tuple(project(g,blocks)for g in out[1:])==star_tuple(m.product(ps),c)
    return out,upper

def canonical(r):
    n=2*r+1;x=2*r
    return (m.perm(n,[(x,2),(0,1)]),)+tuple(m.perm(n,[(0,2*j),(1,2*j+1)])for j in range(r-1,0,-1))

def normal(t,x):
    original=t
    kind,blocks=lower_class(t,x);assert kind=="W"
    r=len(blocks);n=2*r+1
    owner=f.owner(t,x)
    word=[-i for i in range(owner,0,-1)]
    t=apply(t,word,x);assert f.owner(t,x)==0
    d,_=full_deletion(t,x);e=frozenset(next(iter(leaf.edges(d[0]))))
    R=next(i for i,B in enumerate(blocks)if frozenset(B)==e)
    u=t[0][x];U=next(i for i,B in enumerate(blocks)if u in B)
    assert R!=U
    t,w=star_upper(t,x,blocks,R);word+=w
    C=m.product(tuple(project(g,blocks)for g in t[1:]))
    Rp=m.inverse(C)[U]
    push=False
    if Rp!=R:
        at=next(i for i,g in enumerate(t[1:],1)
                if project(g,blocks)==m.perm(r,[(R,Rp)]))
        bring=list(range(at,1,-1))
        t=apply(t,bring,x);word+=bring
        beta=t[1];assert beta[u]==u
        t=apply(t,[1,1],x);word +=[1,1]
        assert f.owner(t,x)==0 and t[0][x]==u
        d,_=full_deletion(t,x)
        assert frozenset(next(iter(leaf.edges(d[0]))))==frozenset(blocks[Rp])
        assert m.product(tuple(project(g,blocks)for g in t[1:]))==C
        push=True
    R=Rp
    t,w=star_upper(t,x,blocks,R);word+=w
    order=[R]
    while len(order)<r:order.append(C[order[-1]])
    assert order[1]==U and C[order[-1]]==R
    orient={R:list(blocks[R])}
    for j in range(r):
        if j==R:continue
        g=next(g for g in t[1:]if project(g,blocks)==m.perm(r,[(R,j)]))
        orient[j]=[g[a]for a in orient[R]]
        assert set(orient[j])==set(blocks[j])
    if orient[U][0]!=u:
        assert orient[U][1]==u
        orient={j:pair[::-1]for j,pair in orient.items()}
    sigma=[-1]*n;sigma[x]=2*r
    for j,k in enumerate(order):
        for bit,a in enumerate(orient[k]):sigma[a]=2*j+bit
    sigma=tuple(sigma);assert sorted(sigma)==list(range(n))
    endpoint=tuple(m.conjugate(sigma,g)for g in t)
    assert endpoint==canonical(r)
    assert m.braid(original,word)==t
    return word,sigma,push


def canonical_certificate(r):
    t=canonical(r);n=2*r+1;x=2*r
    assert m.cycle_lengths(m.product(t))==[n]
    assert lower_class(t,x)[0]=="W"
    A,B1,B2=t[0],t[-1],t[-2]
    C=m.compose(A,B1)
    assert C==m.perm(n,[(0,x,2,1,3)])
    Aprime=m.conjugate(C,A)
    assert Aprime==m.perm(n,[(1,2),(x,3)])
    assert m.power(m.compose(Aprime,B2),2)==m.perm(n,[(1,2,5)])
    degree=leaf.degrees(t)
    removable=[i for i,g in enumerate(t)
        if all(degree[a]==1 or degree[b]==1 for a,b in leaf.edges(g))]
    assert len(removable)==r-2
    smaller=leaf.delete(t,1,[2*r-2,2*r-1])
    assert smaller==canonical(r-1)
    if r==3:
        assert leaf.group_closure(t)==2520
        assert leaf.group_closure(smaller)==10
    return {"r":r,"degree":n,"simultaneous_removable_DTs":r-2,
            "literal_recursive_core":True,"three_cycle_certificate":True}

def complete_degree7():
    alltu=set(f.all_degree7_tuples());unseen=set(alltu)
    counter=tuple(m.perm(7,p)for p in [
        [(0,1),(3,4)],[(1,2),(4,5)],[(2,3),(5,6)]])
    cp=m.product(counter);counter_relabel=[-1]*7;a=0
    for i in range(7):counter_relabel[a]=i;a=cp[a]
    counter_relabel=tuple(counter_relabel)
    normalized_counter=tuple(m.conjugate(counter_relabel,g)for g in counter)
    assert normalized_counter in alltu
    counter_orbit_checked=False
    rows=[];classes={"S":0,"W":0};sigmas={};endpoints={};normal_count=pushes=maxword=0
    while unseen:
        seed=min(unseen);orbit={seed};queue=deque([seed])
        while queue:
            t=queue.popleft()
            for i in range(2):
                for sign in(-1,1):
                    u=m.hurwitz(t,i,sign)
                    assert u in alltu
                    if u not in orbit:orbit.add(u);queue.append(u)
        unseen-=orbit
        group_order=leaf.group_closure(seed)
        counts={"S":0,"W":0};allS=0
        for t in orbit:
            kinds=[]
            for x in range(7):
                if f.owner(t,x)is None:continue
                kind,blocks=lower_class(t,x)
                classes[kind]+=1;counts[kind]+=1;kinds.append(kind)
                if kind=="W":
                    word,sigma,push=normal(t,x)
                    normal_count+=1;pushes+=push;maxword=max(maxword,len(word))
                    endpoint=m.braid(t,word)
                    if x in sigmas:
                        assert sigmas[x]==sigma and endpoints[x]==endpoint
                    else:
                        sigmas[x]=sigma;endpoints[x]=endpoint
            allS+=all(kind=="S"for kind in kinds)
        if normalized_counter in orbit:
            assert len(orbit)==147 and group_order==2520 and counts["W"]==0
            counter_orbit_checked=True
        rows.append({"tuple_count":len(orbit),"upper_group_order":group_order,
                     "symmetric_leaf_choices":counts["S"],
                     "imprimitive_leaf_choices":counts["W"],
                     "tuples_all_leaf_deletions_symmetric":allS})
    rows.sort(key=lambda z:(z["upper_group_order"],z["imprimitive_leaf_choices"]))
    assert len(alltu)==392 and classes=={"S":1008,"W":126}
    assert normal_count==126 and len(sigmas)==7 and counter_orbit_checked
    assert [(z["tuple_count"],z["upper_group_order"],z["imprimitive_leaf_choices"])for z in rows]==[
        (49,168,0),(49,168,0),(147,2520,0),(147,2520,126)]
    counter=tuple(m.perm(7,p)for p in [
        [(0,1),(3,4)],[(1,2),(4,5)],[(2,3),(5,6)]])
    assert leaf.group_closure(counter)==2520
    chosen=[x for x in range(7)if f.owner(counter,x)is not None]
    assert chosen==[0,6] and all(lower_class(counter,x)[0]=="S"for x in chosen)
    assert all(not f.removable_owner(counter,x)for x in chosen)
    return {"complete_fixed_product_domain":392,"orbits":rows,
            "all_leaf_choices_classified":1134,"classification_totals":classes,
            "actual_Wr_inputs_normalized":normal_count,
            "root_adjustments_used":pushes,"maximum_word_length":maxword,
            "fixed_product_and_leaf_same_endpoint_verified":True,
            "explicit_counter_entire_orbit_has_only_symmetric_deletions":True,
            "explicit_all_symmetric_full_A7_counter":counter,
            "counter_full_cycle_normalization":counter_relabel}

def pruefer_tree(seq):
    n=len(seq)+2;deg=[1]*n
    for a in seq:deg[a]+=1
    leaves=[a for a,d in enumerate(deg)if d==1];heapq.heapify(leaves);edges=[]
    for a in seq:
        l=heapq.heappop(leaves);edges.append((l,a));deg[l]-=1;deg[a]-=1
        if deg[a]==1:heapq.heappush(leaves,a)
    edges.append(tuple(leaves))
    return edges

def block_lift(r,edges,twists,root,attachment,bit,ownerpos):
    n=2*r+1;x=2*r
    A=m.perm(n,[(x,2*attachment+bit),(2*root,2*root+1)])
    beta=[]
    for(a,b),twist in zip(edges,twists):
        beta.append(m.perm(n,[(2*a,2*b+twist),(2*a+1,2*b+1-twist)]))
    return tuple(beta[:ownerpos])+(A,)+tuple(beta[ownerpos:])

def all_four_block_controls():
    count=maxword=pushes=0
    for seq in product(range(4),repeat=2):
        tree=pruefer_tree(seq)
        for ordered in permutations(tree):
            for twists in product(range(2),repeat=3):
                for ownerpos in range(4):
                    t=block_lift(4,ordered,twists,0,3,0,ownerpos)
                    assert m.cycle_lengths(m.product(t))==[9]
                    word,sigma,push=normal(t,8)
                    count+=1;maxword=max(maxword,len(word));pushes+=push
    assert count==3072
    return {"all_labelled_block_trees":16,"all_edge_orders_per_tree":6,
            "all_twist_vectors":8,"all_owner_positions":4,
            "root_block":0,"attachment_block":3,
            "actual_inputs_normalized":count,
            "root_adjustments_used":pushes,"maximum_word_length":maxword}

def larger_controls():
    rng=random.Random(21031);report=[]
    for r,samples in((5,12),(6,12),(12,4),(25,2),(50,1)):
        largest=0
        for sample in range(samples):
            seq=[rng.randrange(r)for _ in range(r-2)]
            tree=pruefer_tree(seq);rng.shuffle(tree)
            root,attachment=rng.sample(range(r),2)
            t=block_lift(r,tree,[rng.randrange(2)for _ in tree],
                         root,attachment,rng.randrange(2),rng.randrange(r))
            n=2*r+1
            sigma=list(range(n));rng.shuffle(sigma);sigma=tuple(sigma)
            x=sigma[2*r];t=tuple(m.conjugate(sigma,g)for g in t)
            word,relabel,push=normal(t,x)
            largest=max(largest,len(word))
        report.append({"r":r,"degree":2*r+1,"deterministic_arbitrary_tree_controls":samples,
                       "maximum_word_length":largest})
    return report

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=HERE/"imprimitive_leaf_normal_form_verification.json")
    args=parser.parse_args()
    census=complete_degree7()
    block_controls=all_four_block_controls()
    higher=larger_controls()
    canonical_controls=[canonical_certificate(r)for r in(3,4,5,6,12,25,50)]
    paths=["IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md",
           "FIXED_LEAF_HURWITZ_OBSTRUCTION.md",
           "verify_fixed_leaf_hurwitz_obstruction.py",
           "PAIRED_HURWITZ_LEAF_EXTRACTION.md",
           "verify_paired_hurwitz_leaf_extraction.py",
           "verify_unmarked_polynomial_hurwitz_reduction.py",
           "alternating_hurwitz_source_ledger.json"]
    report={"schema_version":1,"kind":"imprimitive_leaf_deletion_normal_form",
        "result":"PASS","python":platform.python_version(),
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()for p in paths},
        "complete_degree7_dichotomy_and_normal_form":census,
        "ordered_twisted_four_block_inputs":block_controls,
        "higher_degree_controls":higher,"canonical_certificates":canonical_controls,
        "scope":"Written uniform dichotomy and actual paired-block normal form/extraction for every imprimitive chosen-leaf deletion. The symmetric-deletion case and global uniform connectivity remain open.",
        "independent_replay":False,"formal_verification":False,"external_review":False}
    args.output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()

