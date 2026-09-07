#!/usr/bin/env sage -python
import json
from pathlib import Path
from sage.all import GF, PolynomialRing, QQ, ZZ, Zmod, matrix

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

# Use every currently passing q32 signature for which the corresponding
# horizontal-section artifact exists.
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

if len(records)<3:
    raise SystemExit("need at least 3 complete modular q32 signatures")

# Put 100003 first, otherwise preserve numeric order.
records.sort(key=lambda t:(0 if t[0]==100003 else 1,int(t[0])))

def intrinsic_plane(p,sig,q24):
    F=GF(p)
    R=PolynomialRing(F,"U")
    U=R.gen()

    sec=q24["section_mod_p"]
    Z=R([F(v) for v in sec["Z_coefficients_low_to_high"]])
    X=R([F(v) for v in sec["X_coefficients_low_to_high"]])
    Y=R([F(v) for v in sec["Y_coefficients_low_to_high"]])

    P=matrix(F,sig["plane_rref_2x56"])
    assert P.dimensions()==(2,56)

    triples=[]
    allpolys=[]
    for row in P.rows():
        A=sum(row[i]*U**i for i in range(41))
        B=sum(row[41+i]*U**i for i in range(15))
        Cx=A*Z**2
        Cy=B*Z**3
        C0=-A*X+B*Y
        triples.append((Cx,Cy,C0))
        allpolys += [Cx,Cy,C0]

    # Fixed divisor common to the whole pencil. Normalize monic before divide.
    g=allpolys[0]
    for f in allpolys[1:]:
        g=g.gcd(f)
    if g:
        g=g.monic()
    else:
        g=R.one()

    triples=[tuple(f//g for f in tr) for tr in triples]

    maxdx=max(f.degree() for tr in triples for f in [tr[0]])
    maxdy=max(f.degree() for tr in triples for f in [tr[1]])
    maxd0=max(f.degree() for tr in triples for f in [tr[2]])

    # Use fixed generous slots so the coordinate system is identical at every p.
    DX,DY,D0=88,86,92
    rows=[]
    for Cx,Cy,C0 in triples:
        row=[]
        row += [Cx[i] for i in range(DX+1)]
        row += [Cy[i] for i in range(DY+1)]
        row += [C0[i] for i in range(D0+1)]
        rows.append(row)

    M=matrix(F,rows).echelon_form()
    return {
        "matrix":M,
        "pivots":tuple(M.pivots()),
        "gcd_degree":int(g.degree()),
        "degrees":(int(maxdx),int(maxdy),int(maxd0)),
        "gcd":g,
    }

planes=[]
for p,sig,q24 in records:
    r=intrinsic_plane(p,sig,q24)
    planes.append((p,r))
    print(
        "Q32INTRINSIC_PRIME|"
        f"prime={p}|gcd_degree={r['gcd_degree']}|"
        f"degrees={','.join(map(str,r['degrees']))}|"
        f"pivots={','.join(map(str,r['pivots']))}|status=PASS",
        flush=True,
    )

# Group primes by identical intrinsic pivot/degree/gcd profile; use the largest
# group, because bad-specialization primes may change a gcd degree.
groups={}
for p,r in planes:
    key=(r["gcd_degree"],r["degrees"],r["pivots"])
    groups.setdefault(key,[]).append((p,r))

best=max(groups.values(),key=len)
print(
    "Q32INTRINSIC_GROUP|"
    f"compatible={len(best)}/{len(planes)}|"
    f"primes={','.join(str(p) for p,_ in best)}|"
    f"gcd_degree={best[0][1]['gcd_degree']}|"
    f"pivots={','.join(map(str,best[0][1]['pivots']))}|status=PASS",
    flush=True,
)

if len(best)<3:
    raise SystemExit("intrinsic planes do not yet have a stable modular profile")

def crt_scalar(residues,mods):
    x=ZZ(0); M=ZZ(1)
    for rr,pp in zip(residues,mods):
        rr=ZZ(rr)%pp
        t=((rr-x)%pp)*((M%pp).inverse_mod(pp))%pp
        x=(x+M*t)%(M*pp)
        M*=pp
    return x,M

def rr_scalar(residues,mods):
    x,M=crt_scalar(residues,mods)
    try:
        q=QQ(Zmod(M)(x).rational_reconstruction())
        return q,M
    except Exception:
        return None,M

def red(q,p):
    q=QQ(q); p=ZZ(p)
    d=ZZ(q.denominator())%p
    if not d:
        return None
    return int((ZZ(q.numerator())%p)*d.inverse_mod(p)%p)

# Train on all but last compatible prime, hold the last out.
train=best[:-1]
holdp,holdr=best[-1]
mods=[p for p,_ in train]
mats=[r["matrix"] for _,r in train]
H=holdr["matrix"]
nrows,ncols=H.dimensions()

total=nrows*ncols
recovered=held=0
maxnb=maxdb=0
vals=[]

for i in range(nrows):
    row=[]
    for j in range(ncols):
        residues=[int(M[i,j]) for M in mats]
        q,mod=rr_scalar(residues,mods)
        if q is None:
            row.append(None)
            continue
        recovered+=1
        maxnb=max(maxnb,abs(ZZ(q.numerator())).nbits())
        maxdb=max(maxdb,abs(ZZ(q.denominator())).nbits())
        if red(q,holdp)==int(H[i,j]):
            held+=1
        row.append(q)
    vals.append(row)

print(
    "Q32INTRINSIC_CRT|"
    f"train={len(train)}|train_primes={','.join(str(p) for p,_ in train)}|"
    f"holdout={holdp}|recovered={recovered}/{total}|heldout={held}/{total}|"
    f"max_num_bits={maxnb}|max_den_bits={maxdb}|"
    f"status={'PASS_HELDOUT' if held==total else 'PARTIAL'}",
    flush=True,
)

# More useful diagnostic: how many held-out successes are non-forced entries?
pivot=set(holdr["pivots"])
forced=0
nonforced_held=0
for i in range(nrows):
    for j in range(ncols):
        # In RREF pivot columns are structurally forced.
        if j in pivot:
            forced+=1
            continue
        q=vals[i][j]
        if q is not None and red(q,holdp)==int(H[i,j]):
            nonforced_held+=1

print(
    "Q32INTRINSIC_SIGNAL|"
    f"forced_entries={forced}|nonforced_heldout={nonforced_held}|"
    f"nonforced_total={total-forced}|"
    f"status={'PROMISING' if nonforced_held>0 else 'NO_NONTRIVIAL_SIGNAL'}",
    flush=True,
)

out=LOCAL/"q32-intrinsic-plane-crt-diagnostic.json"
payload={
    "schema":"elkies-k3.h3-q32-intrinsic-plane-crt-diagnostic.v1",
    "compatible_primes":[int(p) for p,_ in best],
    "training_primes":[int(p) for p,_ in train],
    "heldout_prime":int(holdp),
    "profile":{
        "gcd_degree":best[0][1]["gcd_degree"],
        "degrees":list(best[0][1]["degrees"]),
        "pivots":list(best[0][1]["pivots"]),
    },
    "crt":{
        "total":total,
        "recovered":recovered,
        "heldout":held,
        "forced_entries":forced,
        "nonforced_heldout":nonforced_held,
        "max_num_bits":maxnb,
        "max_den_bits":maxdb,
    },
}
out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
