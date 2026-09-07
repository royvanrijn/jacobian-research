#!/usr/bin/env sage -python
import json, math
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ, matrix

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

# ---------------------------------------------------------------------------
# Rebuild the compact intrinsic q32 2x125 plane at every available good prime.
# ---------------------------------------------------------------------------
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
    raise SystemExit("need at least 3 complete modular signatures")

DA,DY,D0=40,38,44
NC=125

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
        C,r=(-A*X+B*Y).quo_rem(Z**2)
        assert r==0
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
    "Q32ROWLLL_INPUT|"
    f"primes={','.join(str(p) for p,_ in planes)}|count={len(planes)}|"
    "coords=125|pivots=0,1|status=PASS",
    flush=True,
)

# Keep the final prime completely out of the reconstruction.
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
    # centered representative helps the final lattice row a little
    if x>M//2:
        x-=M
    return x,M

# Same CRT modulus for every coordinate.
_,MOD=crt_scalar([0]*len(mods),mods)
print(
    "Q32ROWLLL_MODULUS|"
    f"train={len(train)}|holdout={holdp}|"
    f"modulus_bits={MOD.nbits()}|status=PASS",
    flush=True,
)

nonpivot=list(range(2,NC))
N=len(nonpivot)
assert N==123

results=[]
for ri in range(2):
    residues=[]
    for j in nonpivot:
        x,M=crt_scalar([int(P[ri,j]) for _,P in train],mods)
        assert M==MOD
        residues.append(x)

    # L = { (a_2,...,a_124,d) : a_j == d*r_j mod MOD }.
    #
    # Rows MOD*e_j plus (r_2,...,r_124,1) generate this lattice.
    dim=N+1
    B=matrix(ZZ,dim,dim)
    for j in range(N):
        B[j,j]=MOD
    for j,r in enumerate(residues):
        B[N,j]=r
    B[N,N]=1

    print(
        "Q32ROWLLL_START|"
        f"row={ri}|dimension={dim}|modulus_bits={MOD.nbits()}|status=START",
        flush=True,
    )
    R=B.LLL(delta=0.99)
    print(
        "Q32ROWLLL_REDUCED|"
        f"row={ri}|dimension={dim}|status=PASS",
        flush=True,
    )

    candidates=[]
    # Test all reduced basis vectors; the desired row should be exceptionally
    # short if its primitive integer coordinates fit beneath the modulus.
    for bi,v in enumerate(R.rows()):
        d=ZZ(v[N])
        if d==0:
            continue
        nums=[ZZ(v[j]) for j in range(N)]
        g=abs(d)
        for a in nums:
            g=math.gcd(int(g),abs(int(a)))
        g=ZZ(g)
        if g>1:
            d//=g
            nums=[a//g for a in nums]
        if d<0:
            d=-d
            nums=[-a for a in nums]

        if d%holdp==0:
            matches=-1
        else:
            dinv=(d%holdp).inverse_mod(holdp)
            matches=sum(
                int((a%holdp)*dinv%holdp)==int(H[ri,j])
                for a,j in zip(nums,nonpivot)
            )

        bits=max(
            [abs(d).nbits()]
            +[abs(a).nbits() for a in nums]
        )
        norm2=sum(a*a for a in nums)+d*d
        candidates.append((matches,bits,norm2,bi,d,nums))

    if not candidates:
        raise RuntimeError(f"LLL row {ri}: no nonzero-denominator candidate")

    # Primary score is held-out prediction. This prime was never used by LLL,
    # so a full 123/123 hit is decisive.
    candidates.sort(key=lambda x:(-x[0],x[1],x[2]))
    best=candidates[0]
    matches,bits,norm2,bi,d,nums=best

    print(
        "Q32ROWLLL_BEST|"
        f"row={ri}|basis_index={bi}|heldout={matches}/{N}|"
        f"height_bits={bits}|den_bits={abs(d).nbits()}|"
        f"status={'PASS_HELDOUT' if matches==N else 'PARTIAL'}",
        flush=True,
    )

    # Also expose shortest vectors, useful if the exact row is not quite in
    # the reduced basis yet.
    bynorm=sorted(candidates,key=lambda x:(x[2],x[1]))[:5]
    for rank,c in enumerate(bynorm):
        mm,bb,nn,ii,dd,aa=c
        print(
            "Q32ROWLLL_SHORT|"
            f"row={ri}|rank={rank}|basis_index={ii}|"
            f"heldout={mm}/{N}|height_bits={bb}|den_bits={abs(dd).nbits()}",
            flush=True,
        )

    results.append({
        "row":ri,
        "heldout_matches":int(matches),
        "height_bits":int(bits),
        "denominator":str(d),
        "numerators":[str(a) for a in nums],
        "pass_heldout":bool(matches==N),
    })

complete=all(r["pass_heldout"] for r in results)

payload={
    "schema":"elkies-k3.h3-q32-compact-row-lll.v1",
    "status":(
        "PASS_Q32_COMPACT_INTRINSIC_ROWS_HELDOUT"
        if complete else
        "PARTIAL_Q32_COMPACT_INTRINSIC_ROWS_LLL"
    ),
    "training_primes":[int(p) for p,_ in train],
    "heldout_prime":int(holdp),
    "modulus_bits":int(MOD.nbits()),
    "coordinate_degrees":[DA,DY,D0],
    "rows":results,
}
if complete:
    # Store the exact rational RREF explicitly.
    qrows=[]
    for ri,res in enumerate(results):
        d=ZZ(res["denominator"])
        nums=[ZZ(x) for x in res["numerators"]]
        row=["0"]*NC
        row[ri]="1"
        for a,j in zip(nums,nonpivot):
            row[j]=str(a/d)
        qrows.append(row)
    payload["qq_plane_rref_2x125"]=qrows

out=LOCAL/"q32-compact-intrinsic-row-lll.json"
out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
print(
    "Q32ROWLLL_RESULT|"
    f"rows_passed={sum(r['pass_heldout'] for r in results)}/2|"
    f"status={payload['status']}",
    flush=True,
)
