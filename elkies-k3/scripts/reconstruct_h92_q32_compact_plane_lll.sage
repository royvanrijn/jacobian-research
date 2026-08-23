#!/usr/bin/env sage -python
import json, math
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ, matrix

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

DA,DY,D0=40,38,44
NC=125

records=[]
for spath in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        sig=json.loads(spath.read_text())
        p=ZZ(sig["prime"])
    except Exception:
        continue
    if sig.get("status")!="PASS_Q32_MODP_SIGNATURE":
        continue
    qpath=LOCAL/f"q24-degree46-direct-global-mod-{p}.json"
    if not qpath.exists():
        continue
    try:
        q24=json.loads(qpath.read_text())
    except Exception:
        continue
    if q24.get("status")!="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE":
        continue
    records.append((p,sig,q24))

records.sort(key=lambda t:(0 if t[0]==100003 else 1,int(t[0])))
if len(records)<3:
    raise SystemExit("need >=3 complete modular primes")

def compact_plane(p,sig,q24):
    F=GF(p)
    R=PolynomialRing(F,"U")
    U=R.gen()
    sec=q24["section_mod_p"]
    Z=R([F(v) for v in sec["Z_coefficients_low_to_high"]])
    X=R([F(v) for v in sec["X_coefficients_low_to_high"]])
    Y=R([F(v) for v in sec["Y_coefficients_low_to_high"]])

    P=matrix(F,sig["plane_rref_2x56"])
    rows=[]
    for row in P.rows():
        A=sum(row[i]*U**i for i in range(41))
        B=sum(row[41+i]*U**i for i in range(15))
        C,rem=(-A*X+B*Y).quo_rem(Z**2)
        assert rem==0
        BY=B*Z
        assert A.degree()<=DA and BY.degree()<=DY and C.degree()<=D0
        rows.append(
            [A[i] for i in range(DA+1)]
            +[BY[i] for i in range(DY+1)]
            +[C[i] for i in range(D0+1)]
        )

    M=matrix(F,rows).echelon_form()
    assert M.dimensions()==(2,NC)
    assert tuple(M.pivots())==(0,1)
    return M

planes=[(p,compact_plane(p,sig,q24)) for p,sig,q24 in records]
print(
    "Q32PLANELLL_INPUT|"
    f"primes={','.join(str(p) for p,_ in planes)}|"
    f"count={len(planes)}|coords={NC}|status=PASS",
    flush=True,
)

# Last prime is not used in CRT/LLL at all.
train=planes[:-1]
holdp,H=planes[-1]
mods=[p for p,_ in train]

def crt_scalar(vals,mods):
    x=ZZ(0); M=ZZ(1)
    for rr,p in zip(vals,mods):
        rr=ZZ(rr)%p
        t=((rr-x)%p)*((M%p).inverse_mod(p))%p
        x=(x+M*t)%(M*p)
        M*=p
    if x>M//2:
        x-=M
    return x,M

_,MOD=crt_scalar([0]*len(mods),mods)
print(
    "Q32PLANELLL_MODULUS|"
    f"train={len(train)}|holdout={holdp}|"
    f"bits={MOD.nbits()}|status=PASS",
    flush=True,
)

# CRT-lift the two modular RREF rows.
Rcrt=[[ZZ(0)]*NC for _ in range(2)]
Rcrt[0][0]=ZZ(1)
Rcrt[1][1]=ZZ(1)
for ri in range(2):
    for j in range(2,NC):
        x,M=crt_scalar([int(P[ri,j]) for _,P in train],mods)
        assert M==MOD
        Rcrt[ri][j]=x

# Congruence lattice:
#   row 0 = CRT RREF row 0
#   row 1 = CRT RREF row 1
#   rows j>=2 = M e_j
#
# Every vector in this lattice reduces into the q32 plane modulo every
# training prime. LLL may choose any GL_2-compatible small basis.
B=matrix(ZZ,NC,NC)
for j in range(NC):
    B[0,j]=Rcrt[0][j]
    B[1,j]=Rcrt[1][j]
for j in range(2,NC):
    B[j,j]=MOD

print(
    "Q32PLANELLL_START|"
    f"dimension={NC}|det_bits={123*MOD.nbits()}|status=START",
    flush=True,
)
L=B.LLL(delta=0.99)
print("Q32PLANELLL_REDUCED|dimension=125|status=PASS",flush=True)

def primitive(v):
    vals=[ZZ(x) for x in v]
    g=ZZ(0)
    for x in vals:
        g=ZZ(math.gcd(int(g),abs(int(x))))
    if g>1:
        vals=[x//g for x in vals]
    # deterministic sign
    for x in vals:
        if x:
            if x<0:
                vals=[-a for a in vals]
            break
    return vals

def member_mod(v,p,P):
    F=GF(p)
    vv=[F(ZZ(x)%p) for x in v]
    # RREF pivots 0,1: first two coordinates determine the row combination.
    a,b=vv[0],vv[1]
    return all(
        vv[j] == a*P[0,j]+b*P[1,j]
        for j in range(2,NC)
    )

def height_bits(v):
    return max(abs(ZZ(x)).nbits() for x in v)

tested=[]
full=[]
for bi,row in enumerate(L.rows()):
    raw=[ZZ(x) for x in row]
    v=primitive(raw)

    train_hits=sum(member_mod(v,p,P) for p,P in train)
    hold=member_mod(v,holdp,H)
    h=height_bits(v)
    norm2=sum(x*x for x in v)

    tested.append((hold,train_hits,h,norm2,bi,v))
    if train_hits==len(train) and hold:
        full.append((h,norm2,bi,v))

tested.sort(key=lambda t:(not t[0],-t[1],t[2],t[3]))
for rank,t in enumerate(tested[:12]):
    hold,th,h,n2,bi,v=t
    print(
        "Q32PLANELLL_CANDIDATE|"
        f"rank={rank}|basis_index={bi}|train={th}/{len(train)}|"
        f"heldout={int(hold)}|height_bits={h}|"
        f"v01={v[0]},{v[1]}|"
        f"status={'PASS_ALL_PRIMES' if th==len(train) and hold else 'PARTIAL'}",
        flush=True,
    )

# Extract two QQ-independent candidates that pass the unseen prime.
chosen=[]
for h,n2,bi,v in sorted(full,key=lambda x:(x[0],x[1])):
    if not chosen:
        chosen.append((h,n2,bi,v))
    else:
        M=matrix(ZZ,[chosen[0][3],v])
        if M.rank()==2:
            chosen.append((h,n2,bi,v))
            break

complete=len(chosen)==2

# Strong validation: their reductions must span the full modular q32 plane at
# every available prime, not merely be individual members.
prime_span_ok={}
if complete:
    C=matrix(ZZ,[x[3] for x in chosen])
    assert C.rank()==2
    for p,P in planes:
        F=GF(p)
        Cp=matrix(F,[[F(ZZ(x)%p) for x in row] for row in C.rows()])
        ok=(Cp.rank()==2 and Cp.row_space()==P.row_space())
        prime_span_ok[int(p)]=bool(ok)
    complete=complete and all(prime_span_ok.values())

for idx,item in enumerate(chosen):
    h,n2,bi,v=item
    print(
        "Q32PLANELLL_CHOSEN|"
        f"i={idx}|basis_index={bi}|height_bits={h}|"
        f"v01={v[0]},{v[1]}|status=PASS_ALL_PRIMES",
        flush=True,
    )

if chosen:
    print(
        "Q32PLANELLL_SPAN|"
        f"chosen={len(chosen)}|"
        f"all_prime_spans={sum(prime_span_ok.values())}/{len(prime_span_ok)}|"
        f"status={'PASS' if complete else 'PARTIAL'}",
        flush=True,
    )

payload={
    "schema":"elkies-k3.h3-q32-compact-plane-lll.v1",
    "status":(
        "PASS_Q32_COMPACT_INTRINSIC_INTEGER_PLANE_ALL_PRIMES"
        if complete else
        "PARTIAL_Q32_COMPACT_INTRINSIC_PLANE_LLL"
    ),
    "training_primes":[int(p) for p,_ in train],
    "heldout_prime":int(holdp),
    "modulus_bits":int(MOD.nbits()),
    "coordinate_degrees":[DA,DY,D0],
    "candidate_count_all_primes":len(full),
    "chosen":[
        {
            "basis_index":int(bi),
            "height_bits":int(h),
            "coordinates":[str(x) for x in v],
        }
        for h,n2,bi,v in chosen
    ],
    "prime_span_validation":prime_span_ok,
}
out=LOCAL/"q32-compact-intrinsic-plane-lll.json"
out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
print(
    "Q32PLANELLL_RESULT|"
    f"all_prime_candidates={len(full)}|chosen={len(chosen)}|"
    f"status={payload['status']}",
    flush=True,
)
