#!/usr/bin/env python3
"""Exact two-sweep/core/endpoint certificate for one-lobe Wr absorption.

The fixed thirteen-sheet controller EXISTS by the retained complete finite
Hurwitz theorem. No recovered literal controller word or its replay is claimed.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import verify_uniform_hurwitz_attachment as h

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CORE_CYCLES=[[(0,1),(2,8)],[(0,6),(10,11)],[(0,4),(11,12)],
             [(0,2),(1,7)],[(0,9),(1,5)],[(9,10),(1,3)]]
ABSOLUTE_CORE_STATES=446265625
CONTROLLER_BOUND=ABSOLUTE_CORE_STATES-1+12*30

def wr(r):
    n=2*r+1
    return (h.permutation(n,[(2*r,2),(0,1)]),)+tuple(
        h.permutation(n,[(0,2*j),(1,2*j+1)])
        for j in range(r-1,0,-1))

def initial(r):
    assert r>=4
    n=2*r+5
    old=tuple(g+tuple(range(2*r+1,n)) for g in wr(r))
    alpha,beta,gamma,delta=range(2*r+1,2*r+5)
    return old+(h.permutation(n,[(0,alpha),(beta,gamma)]),
                h.permutation(n,[(alpha,beta),(gamma,delta)]))

def replay(t,word):
    for code in word:
        t=h.hurwitz(t,abs(code)-1,1 if code>0 else -1)
    return t

def prefix_word(r):
    ans=[]
    for j in range(r,1,-1):ans.extend([-j]*3)
    ans.extend([r+1]*3)
    for j in range(r,2,-1):ans.extend([-j]*3)
    assert len(ans)==6*(r-1)
    # Bring the final three core factors left across the outside factors.
    for idx in range(r-2,2,-1):
        ans.extend(range(idx+1,idx+4))
    assert len(ans)==9*r-18
    return ans

def core_embedding(r):
    ans={i:i for i in (0,1,2,3,5,7)}
    ans.update({4:2*(r-2),6:2*(r-1),8:2*r,
                9:2*r+1,10:2*r+2,11:2*r+3,12:2*r+4})
    assert len(set(ans.values()))==13
    return ans

def embed(g,embedding,n):
    ans=list(range(n))
    for i,j in enumerate(g):ans[embedding[i]]=embedding[j]
    return tuple(ans)

def canonical_target(core):
    C=h.product(core);standard=wr(6);P=h.product(standard)
    q=[None]*13
    a,b=12,9
    for _ in range(13):
        q[a]=b;a=P[a];b=C[b]
    assert a==12 and b==9 and sorted(q)==list(range(13))
    q=tuple(q)
    assert q==(0,4,11,6,12,1,10,3,8,5,2,7,9)
    target=tuple(h.conjugate(q,g) for g in standard)
    assert h.product(target)==C
    assert target[1]==h.permutation(13,[(0,2),(4,7)])
    assert target[0]==h.permutation(13,[(9,11),(0,4)])
    assert q[0]==0
    return target,q

def antipodal_pairs(P,leaf):
    cycle=next(c for c in h.cycles(P) if len(c)>1)
    cycle=[a for a in cycle if a!=leaf]
    assert len(cycle)%2==0
    half=len(cycle)//2
    return {frozenset((cycle[j],cycle[j+half])) for j in range(half)}

def normalize_star(t,leaf):
    """Explicit relabeling only, after the exact final inverse-move slide."""
    n=len(t[0]);r=len(t)
    A=t[0];attachment=A[leaf]
    root=next(e for e in h.edges(A) if leaf not in e)
    zero=t[-1][attachment]
    assert zero in root
    one=next(a for a in root if a!=zero)
    labeling={zero:0,one:1,leaf:2*r}
    for index,g in enumerate(t[1:],1):
        j=r-index
        a,b=g[zero],g[one]
        assert a not in labeling and b not in labeling and a!=b
        labeling[a]=2*j;labeling[b]=2*j+1
    assert sorted(labeling)==list(range(n))
    p=tuple(labeling[i]for i in range(n))
    normalized=tuple(h.conjugate(p,g) for g in t)
    assert normalized==wr(r)
    return p

def case(r,core,target):
    t0=initial(r);n=len(t0[0])
    assert len(t0)==r+2 and len(h.cycles(h.product(t0)))==1
    prefix=prefix_word(r);t=replay(t0,prefix);em=core_embedding(r)
    expected_core=tuple(embed(g,em,n) for g in core)
    expected_outside=tuple(h.permutation(n,[(2,2*(j-2)),(7,2*j+1)])
                           for j in range(r-1,3,-1))
    assert t==expected_core+expected_outside
    assert h.product(t)==h.product(t0)
    assert frozenset((2,7)) in antipodal_pairs(h.product(expected_core),2*r+1)
    # Only endpoints are checked for this imported finite controller.
    # Its literal word existence/bound is proved by the frozen n13 theorem.
    substituted=tuple(embed(g,em,n)for g in target)+expected_outside
    assert h.product(substituted)==h.product(t)
    assert sum(g[2*r+1]!=2*r+1 for g in substituted)==1
    slide=final_word(r)
    assert len(slide)==5*(r-4)
    endpoint=replay(substituted,slide)
    labeling=normalize_star(endpoint,2*r+1)
    assert labeling[0]==0
    assert h.product(endpoint)==h.product(t0)
    assert len(prefix)+len(slide)==14*r-38
    return {"old_factor_count":r,"new_factor_count":r+2,"new_degree":n,
            "outside_DT_count":r-4,
            "prefix_moves":len(prefix),"final_commuting_bring_and_inverse_slide_moves":len(slide),
            "moves_outside_imported_fixed_controller":14*r-38,
            "literal_core_and_outside_formula":True,
            "endpoint_is_canonical_Wr_after_explicit_relabeling":True,
            "marked_common_root0_preserved":True,
            "normalizing_relabeling":list(labeling)}


def final_word(r):
    k=r-4
    ans=[]
    for j in range(k):
        ans.extend(range(6+j,2+j,-1))
    ans.extend(-j for j in range(2,k+2))
    assert len(ans)==5*k
    return ans

def bouquet(t):
    n=4*t+1
    return tuple(g for i in range(t)for g in (
        h.permutation(n,[(0,4*i+1),(4*i+2,4*i+3)]),
        h.permutation(n,[(4*i+1,4*i+2),(4*i+3,4*i+4)])))

def base_target():
    seed=bouquet(2);C=h.product(seed);standard=wr(4);P=h.product(standard)
    leaf=C[0];q=[None]*9
    a,b=8,leaf
    for _ in range(9):
        q[a]=b;a=P[a];b=C[b]
    assert q[0]==0 and sorted(q)==list(range(9))
    target=tuple(h.conjugate(tuple(q),g)for g in standard)
    assert h.product(target)==C
    return target,tuple(q)

def base_control():
    target,q=base_target()
    receipt=json.loads((HERE/"degree_fourteen_hurwitz_verification.json").read_text())
    for n,count in ((9,6561),(13,446265625)):
        row=next(a for a in receipt["cases"]if a["kind"]=="DT" and a["polynomial_degree"]==n)
        assert row["complete"] and row["absolute_orbit"]==row["absolute_total_bound"]==count
        chars=next(a for a in receipt["character_bounds"]if a["polynomial_degree"]==n)
        assert chars["literal_DT_count"]==n*count and chars["absolute_DT_count"]==count
        seed=next(a for a in receipt["seed_controls"]if a["degree"]==n)
        assert seed["DT_seed_full_cycle"] and seed["DT_three_cycle_certificate"]
    return {"degree":9,"absolute_orbit_size":6561,
            "literal_controller_bound":6560+8*12,"target_conjugator":q,
            "target_root0_preserved":True,"literal_base_controller_replayed":False}

def bouquet_control(t,core,target):
    total=bouquet(t);n=len(total[0])
    base,q=base_target()
    total=tuple(g+tuple(range(9,n))for g in base)+total[4:]
    relabel=h.inverse(q)+tuple(range(9,n))
    total=tuple(h.conjugate(relabel,g)for g in total)
    moves=0
    for r in range(4,2*t,2):
        total=replay(total,prefix_word(r));em=core_embedding(r)
        assert total[:6]==tuple(embed(g,em,n)for g in core)
        # Finite-controller endpoint substitution, explicitly not a word replay.
        total=tuple(embed(g,em,n)for g in target)+total[6:]
        total=replay(total,final_word(r))
        active=tuple(g[:2*r+5]for g in total[:r+2])
        relabel=normalize_star(active,2*r+1)+tuple(range(2*r+5,n))
        assert relabel[0]==0
        total=tuple(h.conjugate(relabel,g)for g in total)
        moves+=14*r-38
    assert total==wr(2*t)
    assert moves==14*t*t-52*t+48
    return {"lobes":t,"final_factor_count":2*t,"final_degree":n,
            "moves_outside_finitely_imported_controllers":moves,
            "base_controller_count":1,"thirteen_sheet_controller_count":t-2,
            "marked_root0_preserved_at_every_stage":True,
            "final_tuple_canonical":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write-receipt",type=Path,
        default=HERE/"d5_lobe_core_absorption_verification.json")
    args=ap.parse_args()
    local=(h.permutation(7,[(0,1),(3,4)]),
           h.permutation(7,[(1,2),(5,6)]))
    exchanged=(h.permutation(7,[(0,1),(5,6)]),
               h.permutation(7,[(1,2),(3,4)]))
    assert replay(local,[1]*3)==replay(local,[-1]*3)==exchanged
    core=tuple(h.permutation(13,cs)for cs in CORE_CYCLES)
    assert h.cycles(h.product(core))==[(0,9,11,12,10,8,2,4,6,1,3,5,7)]
    three=h.compose(h.compose(core[0],core[1]),h.compose(core[0],core[1]))
    assert three==h.permutation(13,[(0,1,6)])
    target,q=canonical_target(core)
    pairs=antipodal_pairs(h.product(core),9)
    assert pairs=={frozenset(a)for a in [(0,4),(6,11),(2,7),(1,12),(3,10),(5,8)]}
    assert CONTROLLER_BOUND==446265984
    cases=[case(r,core,target)for r in tuple(range(4,13))+(25,50)]
    # Negative controls: the final target root and the outside anchor differ;
    # no claim that the external factors can simply be left at their old roots.
    assert {frozenset((0,4)),frozenset((2,7)),frozenset((6,11))}.__len__()==3
    assert target[0][9]==11 and target[-1][11]==0
    base=base_control()
    bouquet_controls=[bouquet_control(t,core,target)for t in (2,3,4,5,8,12)]
    paths=[Path(__file__).resolve(),HERE/"D5_LOBE_BOUNDED_CORE_ABSORPTION.md",
           HERE/"verify_uniform_hurwitz_attachment.py",
           HERE/"DEGREE_FOURTEEN_HURWITZ_ORBIT.md",
           HERE/"degree_fourteen_hurwitz_verification.json",
           HERE/"verify_degree_fourteen_hurwitz_orbit.py",
           HERE/"IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md"]
    report={"schema_version":1,"result":"PASS","python":platform.python_version(),
        "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in paths},
        "fixed_core_cycles":CORE_CYCLES,"fixed_core_product_cycle":[0,9,11,12,10,8,2,4,6,1,3,5,7],
        "fixed_core_three_cycle":[0,1,6],"fixed_target_conjugator":q,
        "fixed_target_cycles":[[list(c)for c in h.cycles(g)if len(c)>1]for g in target],
        "fixed_core_antipodal_pairs_after_leaf9":sorted(sorted(a)for a in pairs),
        "imported_absolute_core_orbit_size":ABSOLUTE_CORE_STATES,
        "existential_literal_controller_word_bound":CONTROLLER_BOUND,
        "controller_word_recovered_or_replayed":False,
        "exact_bounded_controls":cases,
        "degree_nine_base_control":base,
        "complete_bouquet_endpoint_controls":bouquet_controls,
        "uniform_bound":"446265984+14r-38 ordinary signed Hurwitz moves, then one explicit sheet relabeling fixing0; canonicalWr_r plus ONE root0 D5lobe, r>=4.",
        "bouquet_uniform_bound":"6656+(t-2)*446265984+14t^2-52t+48 moves, t>=2; common-root0 preserved after every stage.",
        "scope":"Complete root-preserving one-lobe conversion and quadratic total upper bound for the unbuffered bouquet, using frozen finite n9/n13 connectivity. No all-tuple connectivity, arbitrary old-root absorption, linear total bouquet bound or unproved interface substitution.",
        "formal_verification":False,"independent_external_replay":False,
        "external_human_review":False}
    args.write_receipt.write_text(json.dumps(report,indent=2)+"\n")
    print("PASS: exact two sweeps, fixed13core, external anchors and final canonical-star slide.")

if __name__=="__main__":main()
