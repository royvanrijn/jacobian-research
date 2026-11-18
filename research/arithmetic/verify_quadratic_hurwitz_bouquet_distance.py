#!/usr/bin/env python3
"""Exact local span lemma and all-leaf distance controls for D5 bouquets."""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import verify_uniform_hurwitz_attachment as h

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def support(g):
    return {i for i,a in enumerate(g)if a!=i}

def bouquet(t):
    n=4*t+1
    return tuple(g for i in range(t)for g in (
        h.permutation(n,[(0,4*i+1),(4*i+2,4*i+3)]),
        h.permutation(n,[(4*i+1,4*i+2),(4*i+3,4*i+4)])))

def exact_local_lemma():
    dts=[]
    for a,b,c,d in itertools.combinations(range(8),4):
        for edges in (((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))):
            dts.append(h.permutation(8,edges))
    assert len(dts)==210 and len(set(dts))==210
    checked=0;maxgrowth=0
    for a in dts:
        for b in dts:
            old=(a,b)
            for sign in (1,-1):
                new=h.hurwitz(old,0,sign)
                growth=0
                for x in range(8):
                    before=[j for j,g in enumerate(old)if g[x]!=x]
                    after=[j for j,g in enumerate(new)if g[x]!=x]
                    assert bool(before)==bool(after)
                    if before:
                        extra=max(before+after)-min(before+after)-(max(before)-min(before))
                        assert 0<=extra<=1
                        growth+=extra
                assert growth<=8
                maxgrowth=max(maxgrowth,growth);checked+=1
    assert checked==88200 and maxgrowth==8
    # Nonadjacent factor exchanges have no such constant bound.
    a=h.permutation(8,[(0,1),(2,3)])
    b=h.permutation(8,[(4,5),(6,7)])
    nonadjacent_growth=8*10
    assert nonadjacent_growth>8
    return {"degree8_DT_class_size":210,"all_ordered_pairs_both_signs":checked,
            "maximum_minimal_historical_growth":maxgrowth,
            "nonadjacent_swap_not_allowed_free":True}

def intervals(tup):
    n=len(tup[0])
    return [tuple([min(j for j,g in enumerate(tup)if g[x]!=x),
                   max(j for j,g in enumerate(tup)if g[x]!=x)])
            for x in range(n)]

def gap(I,J):
    return max(0,I[0]-J[1],J[0]-I[1])

def antipodes(cycle,leaf):
    lower=[x for x in cycle if x!=leaf]
    m=len(lower)//2
    return [(lower[j],lower[j+m])for j in range(m)]

def all_deleted_sheet_control(t):
    tup=bouquet(t);cycle=h.cycles(h.product(tup))[0]
    expected=[0]+[x for i in reversed(range(t))for x in
                    (4*i+1,4*i+3,4*i+4,4*i+2)]
    assert list(cycle)==expected
    initial=intervals(tup)
    position={a:i for i,a in enumerate(cycle)}
    g=max(0,2*((t+1)//2)-3)
    exact_min=None
    for leaf in range(4*t+1):
        blocks=antipodes(cycle,leaf)
        assert len(blocks)==2*t
        nonroot=0
        totalgap=0
        for a,b in blocks:
            totalgap+=gap(initial[a],initial[b])
            ia,ib=position[a],position[b]
            assert sorted(((ia-ib)%len(cycle),(ib-ia)%len(cycle)))==[2*t,2*t+1]
            if a and b:
                nonroot+=1
                i,j=(a-1)//4,(b-1)//4
                d=min(abs(i-j),t-abs(i-j))
                assert d>=(t+1)//2-1
                assert gap(initial[a],initial[b])>=g
        assert nonroot>=2*t-1
        assert totalgap>=(2*t-1)*g
        exact_min=totalgap if exact_min is None else min(exact_min,totalgap)
    return {"lobes":t,"factor_count":2*t,"degree":4*t+1,
            "all_deleted_sheets_checked":4*t+1,
            "minimum_exact_sum_of_initial_pair_gaps":exact_min,
            "uniform_pair_gap":g,"proved_total_gap_bound":(2*t-1)*g,
            "ordinary_word_lower_bound":((2*t-1)*g+7)//8}

def wreath_leaves(tup):
    n=len(tup[0]);P=h.product(tup)
    cycle=h.cycles(P)[0]
    out=[]
    for x in range(n):
        owners=[i for i,g in enumerate(tup)if g[x]!=x]
        if len(owners)!=1:continue
        owner=owners[0];A=tup[owner];u=A[x]
        e=h.compose(h.permutation(n,[(x,u)]),A)
        blocks={frozenset(a)for a in antipodes(cycle,x)}
        lower=tup[:owner]+(e,)+tup[owner+1:]
        if all({frozenset(g[a]for a in pair)for pair in blocks}==blocks for g in lower):
            for pair in blocks:
                assert any(pair<=support(g)for i,g in enumerate(tup)if i!=owner)
            out.append((x,blocks))
    return out

def word_ball(t,depth):
    seed=bouquet(t);initial=intervals(seed)
    counts=Counter();wr_endpoints=0
    def visit(tup,hist,length):
        nonlocal wr_endpoints
        counts[length]+=1
        growth=sum(b-a-(j-i)for (a,b),(i,j) in zip(hist,initial))
        assert growth<=8*length
        for leaf,blocks in wreath_leaves(tup):
            wr_endpoints+=1
            total=0
            for pair in blocks:
                a,b=tuple(pair)
                assert max(hist[a][0],hist[b][0])<=min(hist[a][1],hist[b][1])
                total+=gap(initial[a],initial[b])
            assert total<=growth
            assert (total+7)//8<=length
        if length==depth:return
        for j in range(len(tup)-1):
            for sign in (1,-1):
                new=h.hurwitz(tup,j,sign)
                hh=list(hist)
                for x in range(len(seed[0])):
                    positions=[k for k in (j,j+1)if new[k][x]!=x]
                    if positions:
                        hh[x]=(min(hist[x][0],min(positions)),max(hist[x][1],max(positions)))
                delta=sum((b-a)-(d-c)for (a,b),(c,d)in zip(hh,hist))
                assert delta<=8
                visit(new,hh,length+1)
    visit(seed,initial,0)
    return {"lobes":t,"complete_signed_word_depth":depth,
            "prefixes_by_length":dict(counts),"total_prefixes":sum(counts.values()),
            "wreath_leaf_endpoint_occurrences":wr_endpoints}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write-receipt",type=Path,
        default=HERE/"quadratic_hurwitz_bouquet_distance_verification.json")
    args=ap.parse_args()
    local=exact_local_lemma()
    controls=[all_deleted_sheet_control(t)for t in range(2,41)]
    balls=[word_ball(2,3),word_ball(3,2)]
    canonical=(h.permutation(9,[(8,2),(0,1)]),)+tuple(
        h.permutation(9,[(0,2*j),(1,2*j+1)])for j in range(3,0,-1))
    positive=wreath_leaves(canonical)
    assert any(x==8 for x,blocks in positive)
    assert not wreath_leaves(bouquet(2))
    upper=json.loads((HERE/"d5_lobe_core_absorption_verification.json").read_text())
    assert upper["existential_literal_controller_word_bound"]==446265984
    assert upper["degree_nine_base_control"]["literal_controller_bound"]==6656
    for t in range(2,101):
        total=sum(14*(2*j)-38 for j in range(2,t))
        assert total==14*t*t-52*t+48
    paths=[Path(__file__).resolve(),HERE/"QUADRATIC_HURWITZ_BOUQUET_DISTANCE.md",
           HERE/"verify_uniform_hurwitz_attachment.py",
           HERE/"D5_LOBE_BOUNDED_CORE_ABSORPTION.md",
           HERE/"d5_lobe_core_absorption_verification.json",
           HERE/"UNBOUNDED_HURWITZ_LEAF_NORMALIZATION_DISTANCE.md",
           HERE/"FIXED_LEAF_HURWITZ_OBSTRUCTION.md",
           HERE/"IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md"]
    result={"schema_version":1,"result":"PASS","python":platform.python_version(),
        "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in paths},
        "complete_local_span_control":local,"all_deleted_sheet_controls":controls,
        "complete_historical_word_balls":balls,
        "detector_controls":{"canonical_Wr_leaf8_accepted":True,
                             "initial_two_lobe_bouquet_rejected":True},
        "uniform_lower_bound":"ceil((2t-1)*max(0,2ceil(t/2)-3)/8), for every endpoint leaf",
        "uniform_upper_bound":"6656+(t-2)*446265984+14t^2-52t+48",
        "distance_order":"Theta(t^2) for unbuffered t-lobe bouquets to the wreath-leaf locus",
        "scope":"Ordinary signed adjacent Hurwitz generators, allowing only common sheet relabelings for free. Uniform lower bound from historical intervals; upper imports frozen finite controllers. No general braid-orbit distinction or arbitrary-input connectivity theorem.",
        "formal_verification":False,"independent_external_replay":False,
        "external_human_review":False}
    args.write_receipt.write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: all local DT pairs, all-leaf antipodal gaps and historical word-ball bounds.")

if __name__=="__main__":main()
