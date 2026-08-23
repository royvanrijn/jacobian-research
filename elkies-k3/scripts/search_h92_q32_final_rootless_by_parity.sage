#!/usr/bin/env sage -python
import argparse, json
from pathlib import Path
from sage.all import QQ, ZZ, block_diagonal_matrix, lcm, matrix, pari, vector, xgcd

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3/route-scout"
FRAME=LOCAL/"q32-fast-step07-q4-a1-root-adapted-frame.txt"
OUT=LOCAL/"q32-final-rootless-parity-scan.json"

parser=argparse.ArgumentParser()
parser.add_argument("--cap",type=int,default=250000)
parser.add_argument("--pari-gb",type=int,default=4)
args=parser.parse_args()

pari.allocatemem(args.pari_gb*1024**3)

def load_gram(path):
    return matrix(ZZ,[
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])

M=load_gram(FRAME)
assert M.dimensions()==(17,17) and M.det()==948
assert M[:1,:1]==matrix(ZZ,[[2]])

# ------------------------------------------------------------------
# 1. Exact bad parity classes: classes containing a norm-8 vector.
# ------------------------------------------------------------------
q8=pari(M).qfminim(8)
cols=matrix(ZZ,q8[2]).columns()
bad=set()
exact8=0
for z in cols:
    if z*M*z != 8:
        continue
    exact8 += 2
    bad.add(tuple(int(x)%2 for x in z))

print(
    "Q32PARITY_BAD|"
    f"qfminim_count_le8={int(q8[0])}|exact_norm8_signed={exact8}|"
    f"bad_parity_classes={len(bad)}|status=PASS",
    flush=True,
)

# ------------------------------------------------------------------
# 2. q6 MW exact shell.  Test only parity; construct a child once.
# ------------------------------------------------------------------
root=M[:1,:1]
coupling=M[:1,1:]
tail=M[1:,1:]
H=tail-coupling.transpose()*root.inverse()*coupling
scale=lcm(v.denominator() for v in H.list())
B=(scale*H).change_ring(ZZ)
target=QQ(23)/2
assert scale*target in ZZ

res=pari(B).qfminim(ZZ(scale*target),args.cap)
mwcols=matrix(ZZ,res[2]).columns()

def bezout(ns,f):
    cur=ZZ(0); a=[ZZ(0)]*ns.nrows()
    for i,val in enumerate(ns*f):
        if not val: continue
        g,l,r=xgcd(cur,ZZ(val))
        a=[l*x for x in a]
        a[i]+=r
        cur=g
    if abs(cur)!=1: return None
    if cur==-1: a=[-x for x in a]
    return vector(ZZ,a)

def child_root_count(v):
    ns=block_diagonal_matrix(matrix(ZZ,[[0,1],[1,0]]),-M)
    D=vector(ZZ,[3,2]+list(v))
    mate=bezout(ns,D)
    assert mate is not None
    mate-=(mate*ns*mate//2)*D
    ker=matrix(ZZ,[list(D*ns),list(mate*ns)]).right_kernel_matrix()
    child=-(ker*ns*ker.transpose())
    qf=pari(child).qfminim(2)
    return int(qf[0])

exact_shell=0
integral=0
bad_hits=0
good=None
parities_seen=set()

for idx,col in enumerate(mwcols,1):
    for sign in (1,-1):
        mw=sign*vector(ZZ,col)
        if mw*H*mw != target:
            continue
        exact_shell+=1
        num=ZZ(1)-ZZ((coupling*mw)[0])
        if num%2:
            continue
        integral+=1
        r=num//2
        v=vector(ZZ,[r]+list(mw))
        assert v*M*v==12
        parity=tuple(int(x)%2 for x in v)
        parities_seen.add(parity)
        if parity in bad:
            bad_hits+=1
            continue

        roots=child_root_count(v)
        assert roots==0, roots
        good={
            "stored_column_index":idx,
            "mw":list(map(int,mw)),
            "root_coordinate":int(r),
            "witness":list(map(int,v)),
            "parity":list(parity),
            "verified_child_root_count":roots,
        }
        print(
            "Q32PARITY_HIT|"
            f"column={idx}|mw={','.join(map(str,mw))}|r={r}|"
            f"unique_parities_seen={len(parities_seen)}|"
            "child_roots=0|MW=17|status=PASS_ROOTLESS",
            flush=True,
        )
        break
    if good is not None:
        break
    if idx%25000==0:
        print(
            "Q32PARITY_PROGRESS|"
            f"columns={idx}|exact_shell={exact_shell}|"
            f"unique_parities={len(parities_seen)}|bad_hits={bad_hits}|status=RUNNING",
            flush=True,
        )

status="PASS_Q32_ROOTLESS_BY_PARITY" if good else "NO_GOOD_PARITY_IN_CAP"
payload={
    "schema":"elkies-k3.h3-q32-final-rootless-parity.v1",
    "status":status,
    "cap":args.cap,
    "pari_total_count_le_target":int(res[0]),
    "stored_columns":len(mwcols),
    "bad_parity_classes":len(bad),
    "exact_shell_tested":exact_shell,
    "integral_candidates":integral,
    "unique_candidate_parities_seen":len(parities_seen),
    "bad_parity_hits":bad_hits,
    "hit":good,
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q32PARITY_RESULT|"
    f"cap={args.cap}|exact_shell={exact_shell}|unique_parities={len(parities_seen)}|"
    f"status={status}",
    flush=True,
)
