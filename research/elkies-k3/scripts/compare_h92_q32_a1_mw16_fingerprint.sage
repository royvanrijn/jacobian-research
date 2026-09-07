#!/usr/bin/env sage -python
from pathlib import Path
from sage.all import QQ, ZZ, lcm, matrix, pari, vector

ROOT=Path("/Users/royvanrijn/Documents/jacobian-research")
if not (ROOT/"elkies-k3/scripts").is_dir():
    ROOT=Path.home()/"Documents"/"jacobian-research"

NEW=ROOT/"artifacts/local/elkies-k3/route-scout/q32-fast-step07-q4-a1-root-adapted-frame.txt"
HIST=ROOT/"artifacts/generated-results/elkies-k3-h3-mw15-2a1-q4-degree2-first-hit-frames/q4-o0981-r1-n2-d2-4f02793cfc09.txt"

def load_gram(path):
    return matrix(ZZ,[
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])

for p in (NEW,HIST):
    if not p.exists():
        raise SystemExit(f"missing {p}")

N=load_gram(NEW)
H=load_gram(HIST)
assert N.dimensions()==H.dimensions()==(17,17)
assert N[:1,:1]==H[:1,:1]==matrix(ZZ,[[2]])

def mw_height(M):
    root=M[:1,:1]
    c=M[:1,1:]
    t=M[1:,1:]
    return t-c.transpose()*root.inverse()*c

HN=mw_height(N)
HH=mw_height(H)
assert HN.dimensions()==HH.dimensions()==(16,16)

scale=ZZ(1)
for x in list(HN.list())+list(HH.list()):
    scale=lcm(scale,ZZ(QQ(x).denominator()))
IN=(scale*HN).change_ring(ZZ)
IH=(scale*HH).change_ring(ZZ)

def smith_diag(M):
    S,_,_=M.smith_form()
    return tuple(abs(ZZ(S[i,i])) for i in range(S.nrows()))

def theta_counts(M, rational_bounds):
    out={}
    for b in rational_bounds:
        ib=ZZ(scale*QQ(b))
        qf=pari(M).qfminim(ib)
        # PARI count includes both signs.
        out[str(QQ(b))]=int(qf[0])
    return out

bounds=[QQ(2),QQ(3),QQ(4),QQ(5),QQ(6)]
tn=theta_counts(IN,bounds)
th=theta_counts(IH,bounds)

# Historical final q6 witness:
# root coordinate -2, then these 16 MW coordinates.
zh=vector(ZZ,[-5,-1,2,3,-2,1,1,-1,0,2,1,0,-1,1,-1,3])
nh=QQ(zh*HH*zh)
nn=QQ(zh*HN*zh)
assert nh==QQ(23)/2

# Count vectors through the exact final q6 MW shell 23/2. This may be large,
# so use qfminim count only; no child construction.
shell=QQ(23)/2
pari.allocatemem(4*1024**3)
cntN=int(pari(IN).qfminim(ZZ(scale*shell))[0])
cntH=int(pari(IH).qfminim(ZZ(scale*shell))[0])

print(
    "A1MW16_FINGERPRINT|"
    f"frame_equal={int(N==H)}|height_equal={int(HN==HH)}|scale={scale}|"
    f"detN={HN.det()}|detH={HH.det()}|"
    f"smith_equal={int(smith_diag(IN)==smith_diag(IH))}|"
    f"theta_equal={int(tn==th)}|"
    f"q6_shell_count_new={cntN}|q6_shell_count_hist={cntH}|"
    f"hist_final_vector_norm_hist={nh}|hist_final_vector_norm_new={nn}|"
    "status=PASS",
    flush=True,
)
print("A1MW16_SMITH_NEW|"+",".join(map(str,smith_diag(IN))),flush=True)
print("A1MW16_SMITH_HIST|"+",".join(map(str,smith_diag(IH))),flush=True)
print("A1MW16_THETA_NEW|"+json_dump(tn) if False else f"A1MW16_THETA_NEW|{tn}",flush=True)
print(f"A1MW16_THETA_HIST|{th}",flush=True)
