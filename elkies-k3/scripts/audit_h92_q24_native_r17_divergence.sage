#!/usr/bin/env sage -python
"""
Audit where the q24-native D12->A1 suffix diverges from the historical
R17-directed H3 corridor.

For each corresponding stage:
  * reconstruct the native stage and its exact H3->stage NS transport;
  * compare the native positive frame with the historical selected frame by
    exact integral qfisom;
  * pull the pinned R17 fibre into native-stage coordinates and report its
    (a,b), q=ab, and old-fibre degree b.

The first non-isometric corresponding frame is the precise branch divergence.
"""

import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/"elkies-k3"/"scripts"
LOCAL=ROOT/"artifacts/local/elkies-k3/q24-native-suffix"
GEN=ROOT/"artifacts/generated-results"

ENGINE=SCRIPTS/"exact_neighbor_engine.sage"
CLOSE=SCRIPTS/"close_h92_q8_q24_by_q6_translation.sage"
SUFFIX=LOCAL/"q24-native-d12-to-a1.json"
REVERSE=GEN/"elkies-k3-rank17-to-h3-reverse-transport.json"
OUT=LOCAL/"q24-native-divergence-audit.json"

HIST = [
    # stage, source artifact, selected orbit
    ("D12", GEN/"elkies-k3-h3-q6-q8-d13-q24-degree2.json", 85),
    ("A11", GEN/"elkies-k3-h3-d12-o85-q6-degree2.json", 42),
    ("2A5", GEN/"elkies-k3-h3-a11-middle-q8-degree2.json", 922),
    ("3A3", GEN/"elkies-k3-h3-a5a5-c2-q4-degree2.json", 472),
    ("A3+2A2", GEN/"elkies-k3-h3-a3x3-q4-degree2.json", 323),
    ("5A1", GEN/"elkies-k3-h3-mw10-a3a2a2-q4-degree2.json", 207),
    ("4A1", GEN/"elkies-k3-h3-mw12-5a1-q4-degree2-first-hit.json", 52),
    ("3A1", GEN/"elkies-k3-h3-mw13-4a1-q4-degree2-first-hit.json", 114),
    ("2A1", GEN/"elkies-k3-h3-mw14-3a1-q4-degree2-first-hit.json", 498),
    ("A1", GEN/"elkies-k3-h3-mw15-2a1-q4-degree2-first-hit.json", 981),
]

U2=matrix(ZZ,((0,1),(1,0)))

def load_gram(path):
    return matrix(ZZ,[
        [ZZ(x) for x in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])

def ns(frame):
    return block_diagonal_matrix(U2,-frame)

def exact_iso(A,B):
    """Return an integral X with X*A*X^t=B, or None."""
    raw=pari(A).qfisom(pari(B))
    candidates=[]
    if str(raw)!="0":
        S=matrix(ZZ,raw)
        candidates += [S,S.transpose()]
        try:
            I=S.inverse()
            if I.change_ring(ZZ)==I:
                I=I.change_ring(ZZ)
                candidates += [I,I.transpose()]
        except Exception:
            pass
    raw2=pari(B).qfisom(pari(A))
    if str(raw2)!="0":
        S=matrix(ZZ,raw2)
        candidates += [S,S.transpose()]
        try:
            I=S.inverse()
            if I.change_ring(ZZ)==I:
                I=I.change_ring(ZZ)
                candidates += [I,I.transpose()]
        except Exception:
            pass
    seen=set()
    for X in candidates:
        k=tuple(X.list())
        if k in seen: continue
        seen.add(k)
        if X.dimensions()==(17,17) and X*A*X.transpose()==B:
            assert abs(ZZ(X.det()))==1
            return X
    return None

for p in [ENGINE,CLOSE,SUFFIX,REVERSE]+[x[1] for x in HIST]:
    if not p.exists():
        raise SystemExit(f"missing prerequisite: {p}")

# Engine + q24 native D12.
eng={"__name__":"__audit_engine__","__file__":str(ENGINE)}
exec(compile(ENGINE.read_text(),str(ENGINE),"exec"),eng)
splitter=eng["primitive_hyperbolic_split"]
minimize=eng["minimize_child_frame"]

saved=list(sys.argv)
scope={"__name__":"__audit_close__","__file__":str(CLOSE)}
buf=io.StringIO()
try:
    sys.argv=[str(CLOSE)]
    with contextlib.redirect_stdout(buf):
        exec(compile(CLOSE.read_text(),str(CLOSE),"exec"),scope)
finally:
    sys.argv=saved

G_H3=matrix(ZZ,scope["ns"])
D24=vector(ZZ,scope["D24eq"])
sp=splitter(G_H3,D24)
raw=matrix(ZZ,sp["child_frame"])
mi=minimize(raw)
D12=matrix(ZZ,mi["frame"])
T=block_diagonal_matrix(identity_matrix(ZZ,2),matrix(ZZ,mi["basis"]))*matrix(ZZ,sp["transport"])
assert T*G_H3*T.transpose()==ns(D12)

# pinned -> H3 exact basis map
rev=json.loads(REVERSE.read_text())
assert rev["status"]=="PASS_EXACT_PINNED_R17_TO_H3_REVERSE_TRANSPORT"
R=matrix(ZZ,rev["complete_reverse_transport"]["pinned_rank17_to_h3_basis_R"])
assert abs(ZZ(R.det()))==1

suffix=json.loads(SUFFIX.read_text())
assert suffix["status"]=="PASS_Q24_NATIVE_D12_TO_A1"

native_stages=[("D12",D12,T,None)]
current=D12

for item in suffix["steps"]:
    data=json.loads((ROOT/item["search_artifact"]).read_text())
    assert load_gram(ROOT/data["frame"])==current
    rec=next(r for r in data["neighbors"] if int(r["orbit_index"])==int(item["discovered_orbit_index"]))
    adapt=matrix(ZZ,rec["child_root_adapted_basis"])
    neighbor=matrix(ZZ,rec["neighbor_basis"])
    child=matrix(ZZ,rec["child_root_adapted_frame"])
    stepT=block_diagonal_matrix(identity_matrix(ZZ,2),adapt)*neighbor
    assert stepT*ns(current)*stepT.transpose()==ns(child)
    T=stepT*T
    current=child
    native_stages.append((item["target_ade"].replace("A5+A5","2A5").replace("A3+A3+A3","3A3").replace("A2+A2+A3","A3+2A2").replace("A1+A1+A1+A1+A1","5A1").replace("A1+A1+A1+A1","4A1").replace("A1+A1+A1","3A1").replace("A1+A1","2A1"),child,T,item))

assert len(native_stages)==10

records=[]
first_div=None
last_match=None

for idx,((stage,N,Tstage,_),(hname,hpath,horbit)) in enumerate(zip(native_stages,HIST)):
    assert stage==hname,(stage,hname)
    hd=json.loads(hpath.read_text())
    hrec=next(r for r in hd["neighbors"] if int(r["orbit_index"])==horbit)
    H=matrix(ZZ,hrec["child_root_adapted_frame"])
    X=exact_iso(N,H)
    iso=(X is not None)

    if iso:
        last_match=stage
    elif first_div is None:
        first_div=stage

    # pinned basis in this native stage
    C=Tstage*R
    P=C.inverse()
    assert P.change_ring(ZZ)==P
    P=P.change_ring(ZZ)
    D=vector(ZZ,P.row(0))
    Fold=vector(ZZ,[1,0]+[0]*17)
    deg=ZZ(D*ns(N)*Fold)
    if deg<0:
        D=-D
        deg=-deg
    a,b=ZZ(D[0]),ZZ(D[1])
    w=vector(ZZ,D[2:])
    q=a*b
    assert D*ns(N)*D==0
    assert b==deg
    assert w*N*w==2*q

    recout={
        "stage":stage,
        "native_historical_isometric":iso,
        "historical_orbit":horbit,
        "historical_artifact":str(hpath.relative_to(ROOT)),
        "isometry":None if X is None else [[int(v) for v in row] for row in X.rows()],
        "pinned_target":{"a":int(a),"b":int(b),"q":int(q),"degree":int(deg)},
    }
    records.append(recout)

    print(
        "Q24DIV|"
        f"stage={stage}|hist_orbit={horbit}|isometric={int(iso)}|"
        f"target_a={a}|target_b={b}|target_q={q}|degree={deg}|"
        f"status={'MATCH' if iso else 'DIVERGED'}",
        flush=True,
    )

payload={
    "schema":"elkies-k3.h3-q24-native-divergence-audit.v1",
    "status":"PASS_Q24_NATIVE_DIVERGENCE_AUDIT",
    "first_nonisometric_stage":first_div,
    "last_isometric_stage":last_match,
    "stages":records,
    "interpretation":(
        "The first non-isometric stage is the earliest point where the fresh "
        "q24-native first-hit suffix leaves the historical R17-directed positive "
        "frame class. If a previous stage is isometric, its historical next "
        "witness can be transported into native coordinates and checked there "
        "to repair the branch exactly."
    ),
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24DIV_RESULT|"
    f"last_match={last_match}|first_divergence={first_div}|"
    "status=PASS_Q24_NATIVE_DIVERGENCE_AUDIT",
    flush=True,
)
