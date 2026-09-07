#!/usr/bin/env sage -python
"""
Fast final search from the q32-derived A1/MW16 frame to rootless/MW17.

Uses the A1 geometry:
  D.F=2 and rootless requires both old I2/A1 components to have positive
  new-fibre degree, hence the unique dominant label is 1.

For q, the root contribution is 1/2, so only MW vectors on the exact shell

    h(mw) = 2q - 1/2

can contribute.  We therefore skip all other dominant weights/norms.

Search q=2 and q=4 first (smaller shells), then q=6 with an expandable PARI
storage cap.  The existing 10k q6 cache is used as a seen set so larger-cap
runs only test newly exposed candidates.
"""

import argparse
import json
from pathlib import Path
from sage.all import (
    QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix, lcm,
    matrix, pari, vector, xgcd
)

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3/route-scout"
FRAME=LOCAL/"q32-fast-step07-q4-a1-root-adapted-frame.txt"
OLD_CACHE=LOCAL/"q32-fast-step08-q6-mwvectors.json"
OUT=LOCAL/"q32-final-rootless-targeted-shell.json"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--q6-cap",type=int,default=50000)
args=parser.parse_args()

def load_gram(path):
    return matrix(ZZ,[
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])

def bezout_vector_for_pairing(ns,fiber):
    pairings=list(ns*fiber)
    current=ZZ(0)
    result=[ZZ(0)]*ns.nrows()
    for i,value in enumerate(pairings):
        if value==0:
            continue
        divisor,left,right=xgcd(current,ZZ(value))
        result=[left*x for x in result]
        result[i]+=right
        current=divisor
    if abs(current)!=1:
        return None
    if current==-1:
        result=[-x for x in result]
    return vector(ZZ,result)

def child_frame(ns,fiber):
    mate=bezout_vector_for_pairing(ns,fiber)
    if mate is None:
        return None
    sq=ZZ(mate*ns*mate)
    assert sq%2==0
    mate-=(sq//2)*fiber
    kernel=matrix(ZZ,[list(fiber*ns),list(mate*ns)]).right_kernel_matrix()
    child=-(kernel*ns*kernel.transpose())
    assert child.is_positive_definite() and child.det()==948
    return child

def root_rank_count(child):
    qf=pari(child).qfminim(2)
    count=ZZ(qf[0])
    if count==0:
        return 0,0
    return int(matrix(ZZ,qf[2]).rank()),int(count)

if not FRAME.exists():
    raise SystemExit(f"Missing A1/MW16 frame: {FRAME}")

frame=load_gram(FRAME)
assert frame.dimensions()==(17,17)
assert frame.det()==948
cartan=frame[:1,:1]
assert cartan==matrix(ZZ,[[2]])
coupling=frame[:1,1:]
tail=frame[1:,1:]
height=tail-coupling.transpose()*QQ(1)/2*coupling
assert height.is_positive_definite()

U2=matrix(ZZ,((0,1),(1,0)))
ns=block_diagonal_matrix(U2,-frame)

scale=lcm(v.denominator() for v in height.list())
scaled=(scale*height).change_ring(ZZ)

seen_q6=set()
if OLD_CACHE.exists():
    cache=json.loads(OLD_CACHE.read_text())
    for vals in cache.get("vectors",[]):
        vv=tuple(map(int,vals))
        seen_q6.add(vv)
        seen_q6.add(tuple(-x for x in vv))
print(
    "Q32FINAL_SEEN|"
    f"q6_cached_signed={len(seen_q6)}|status=PASS",
    flush=True,
)

def enumerate_mw(q,cap=None):
    target=QQ(2*q)-QQ(1)/2
    bound=ZZ(scale*target)
    result=pari(scaled).qfminim(bound) if cap is None else pari(scaled).qfminim(bound,cap)
    pari_count=int(result[0])
    cols=matrix(ZZ,result[2]).columns()
    candidates=[]
    seen=set()
    for col in cols:
        for sign in (1,-1):
            mw=sign*vector(ZZ,col)
            if mw*height*mw != target:
                continue
            key=tuple(map(int,mw))
            if key in seen:
                continue
            seen.add(key)
            candidates.append(mw)
    return target,pari_count,len(cols),candidates

def search_q(q,cap=None,skip_seen=()):
    target,pari_count,stored,candidates=enumerate_mw(q,cap)
    tested=0
    parity_reject=0
    already=0
    for mw in candidates:
        key=tuple(map(int,mw))
        if key in skip_seen:
            already+=1
            continue

        # A1 dominant label is forced to 1:
        # r*2 + coupling*mw = 1.
        numerator=ZZ(1)-ZZ((coupling*mw)[0])
        if numerator%2:
            parity_reject+=1
            continue
        r=numerator//2
        witness=vector(ZZ,[r]+list(mw))
        assert witness*frame*witness==2*q

        fiber=vector(ZZ,[q//2,2]+list(witness))
        if gcd(tuple(fiber))!=1:
            continue
        child=child_frame(ns,fiber)
        if child is None:
            continue
        tested+=1
        rr,rc=root_rank_count(child)
        if rc!=0:
            continue

        assert rr==0
        print(
            "Q32FINAL_HIT|"
            f"q={q}|mw={','.join(map(str,mw))}|r={r}|"
            f"tested_new_children={tested}|pari_count={pari_count}|stored={stored}|"
            f"fiber={','.join(map(str,fiber))}|"
            "root_data=0,0,1|MW=17|status=PASS_ROOTLESS",
            flush=True,
        )
        return {
            "q":q,
            "mw":list(map(int,mw)),
            "root_coordinate":int(r),
            "witness":list(map(int,witness)),
            "fiber":list(map(int,fiber)),
            "root_data":[0,0,1],
            "mw_rank":17,
            "tested_new_children":tested,
            "pari_count":pari_count,
            "stored_columns":stored,
        }

    complete=(cap is None)
    print(
        "Q32FINAL_SCAN|"
        f"q={q}|target_mw_norm={target}|pari_count={pari_count}|stored={stored}|"
        f"exact_shell_signed={len(candidates)}|already_seen={already}|"
        f"parity_reject={parity_reject}|tested_new_children={tested}|"
        f"complete={int(complete)}|status=NO_ROOTLESS",
        flush=True,
    )
    return None

# Smaller q first: these are exact/full enumerations.
hit=search_q(2)
if hit is None:
    hit=search_q(4)

# q6: reuse the old 10k sample as seen and expand only if needed.
if hit is None:
    hit=search_q(6,args.q6_cap,seen_q6)

status="PASS_Q32_FINAL_ROOTLESS" if hit else "Q32_FINAL_ROOTLESS_PENDING_LARGER_Q6_CAP"
payload={
    "schema":"elkies-k3.h3-q32-final-rootless-targeted-shell.v1",
    "status":status,
    "q6_cap":int(args.q6_cap),
    "forced_A1_label":1,
    "old_q6_seen_signed":len(seen_q6),
    "hit":hit,
    "next":(
        None if hit else
        f"rerun with --q6-cap {2*args.q6_cap}; the old and current samples are "
        "deterministic prefixes, and only newly exposed exact-shell vectors need testing"
    ),
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q32FINAL_RESULT|"
    f"status={status}",
    flush=True,
)
