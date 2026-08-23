#!/usr/bin/env sage -python
import json
from pathlib import Path
from sage.all import GF, PolynomialRing, QQ, ZZ, Zmod, matrix

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

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
    raise SystemExit("need >=3 complete primes")

# Tight intrinsic coordinates after dividing the universal Z^2 fixed factor:
#   H = A*x + (B*Z)*y + C,
#   C=(-A*X+B*Y)/Z^2,
# with deg(A)<=40, deg(BZ)<=38, deg(C)<=44.
DA,DY,D0=40,38,44
NC=(DA+1)+(DY+1)+(D0+1)
assert NC==125

def compact_plane(p,sig,q24):
    F=GF(p)
    R=PolynomialRing(F,"U")
    U=R.gen()

    sec=q24["section_mod_p"]
    Z=R([F(v) for v in sec["Z_coefficients_low_to_high"]])
    X=R([F(v) for v in sec["X_coefficients_low_to_high"]])
    Y=R([F(v) for v in sec["Y_coefficients_low_to_high"]])
    assert (Z.degree(),X.degree(),Y.degree())==(24,52,78)

    P=matrix(F,sig["plane_rref_2x56"])
    assert P.dimensions()==(2,56)

    rows=[]
    actual=[]
    for row in P.rows():
        A=sum(row[i]*U**i for i in range(41))
        B=sum(row[41+i]*U**i for i in range(15))
        num=-A*X+B*Y
        z2=Z**2
        q,r=num.quo_rem(z2)
        assert r==0
        C=q
        BY=B*Z
        assert A.degree()<=DA
        assert BY.degree()<=DY
        assert C.degree()<=D0
        actual.append((A,BY,C))
        rows.append(
            [A[i] for i in range(DA+1)]
            +[BY[i] for i in range(DY+1)]
            +[C[i] for i in range(D0+1)]
        )

    M=matrix(F,rows).echelon_form()
    assert M.dimensions()==(2,125)
    return M,tuple(M.pivots()),tuple(
        max(f.degree() for f in slot)
        for slot in zip(*actual)
    )

planes=[]
for p,sig,q24 in records:
    M,piv,deg=compact_plane(p,sig,q24)
    planes.append((p,M,piv,deg))
    print(
        "Q32COMPACT_PRIME|"
        f"prime={p}|degrees={','.join(map(str,deg))}|"
        f"pivots={','.join(map(str,piv))}|status=PASS",
        flush=True,
    )

groups={}
for rec in planes:
    key=(rec[2],rec[3])
    groups.setdefault(key,[]).append(rec)
best=max(groups.values(),key=len)

print(
    "Q32COMPACT_GROUP|"
    f"compatible={len(best)}/{len(planes)}|"
    f"primes={','.join(str(r[0]) for r in best)}|"
    f"degrees={','.join(map(str,best[0][3]))}|"
    f"pivots={','.join(map(str,best[0][2]))}|"
    "coords=125|status=PASS",
    flush=True,
)
if len(best)<3:
    raise SystemExit("no stable compact profile")

def crt_scalar(residues,mods):
    x=ZZ(0); M=ZZ(1)
    for rr,p in zip(residues,mods):
        rr=ZZ(rr)%p
        t=((rr-x)%p)*((M%p).inverse_mod(p))%p
        x=(x+M*t)%(M*p)
        M*=p
    return x,M

def rr_scalar(residues,mods):
    x,M=crt_scalar(residues,mods)
    try:
        return QQ(Zmod(M)(x).rational_reconstruction()),M
    except Exception:
        return None,M

def red(q,p):
    q=QQ(q)
    d=ZZ(q.denominator())%p
    if not d:
        return None
    return int((ZZ(q.numerator())%p)*d.inverse_mod(p)%p)

train=best[:-1]
hold=best[-1]
mods=[r[0] for r in train]
H=hold[1]
nr,nc=H.dimensions()

recovered=held=0
maxnb=maxdb=0
rowvals=[]
for i in range(nr):
    vals=[]
    for j in range(nc):
        q,_=rr_scalar([int(r[1][i,j]) for r in train],mods)
        vals.append(q)
        if q is None:
            continue
        recovered+=1
        maxnb=max(maxnb,abs(ZZ(q.numerator())).nbits())
        maxdb=max(maxdb,abs(ZZ(q.denominator())).nbits())
        if red(q,hold[0])==int(H[i,j]):
            held+=1
    rowvals.append(vals)

piv=set(hold[2])
forced_total=nr*len(piv)
forced_held=0
nonforced_held=0
nonforced_recovered=0
for i in range(nr):
    for j in range(nc):
        q=rowvals[i][j]
        if j in piv:
            if q is not None and red(q,hold[0])==int(H[i,j]):
                forced_held+=1
        else:
            if q is not None:
                nonforced_recovered+=1
                if red(q,hold[0])==int(H[i,j]):
                    nonforced_held+=1

print(
    "Q32COMPACT_CRT|"
    f"train={len(train)}|holdout={hold[0]}|"
    f"recovered={recovered}/{nr*nc}|heldout={held}/{nr*nc}|"
    f"forced_heldout={forced_held}/{forced_total}|"
    f"nonforced_recovered={nonforced_recovered}/{nr*nc-forced_total}|"
    f"nonforced_heldout={nonforced_held}/{nr*nc-forced_total}|"
    f"max_num_bits={maxnb}|max_den_bits={maxdb}|"
    f"status={'PASS_HELDOUT' if held==nr*nc else 'PARTIAL'}",
    flush=True,
)

# Also measure stability across successive training prefixes. A true small
# rational coefficient should reconstruct to the same QQ value before the
# held-out test, unlike a boundary artefact.
stable=0
stable_and_held=0
for i in range(nr):
    for j in range(nc):
        if j in piv:
            continue
        qs=[]
        for k in range(3,len(train)+1):
            q,_=rr_scalar([int(r[1][i,j]) for r in train[:k]],[r[0] for r in train[:k]])
            qs.append(q)
        nonnull=[q for q in qs if q is not None]
        if len(nonnull)>=2 and nonnull[-1]==nonnull[-2]:
            stable+=1
            if red(nonnull[-1],hold[0])==int(H[i,j]):
                stable_and_held+=1

print(
    "Q32COMPACT_STABILITY|"
    f"stable_last_two={stable}|stable_and_heldout={stable_and_held}|"
    f"nonforced_total={nr*nc-forced_total}|"
    f"status={'REAL_SIGNAL' if stable_and_held>0 else 'NO_SIGNAL_YET'}",
    flush=True,
)

out=LOCAL/"q32-compact-intrinsic-crt-diagnostic.json"
out.write_text(json.dumps({
    "schema":"elkies-k3.h3-q32-compact-intrinsic-crt.v1",
    "compatible_primes":[int(r[0]) for r in best],
    "training_primes":[int(r[0]) for r in train],
    "heldout_prime":int(hold[0]),
    "coordinate_degrees":[DA,DY,D0],
    "coordinates_per_row":NC,
    "pivots":list(hold[2]),
    "crt":{
        "recovered":recovered,
        "heldout":held,
        "forced_heldout":forced_held,
        "nonforced_recovered":nonforced_recovered,
        "nonforced_heldout":nonforced_held,
        "stable_last_two":stable,
        "stable_and_heldout":stable_and_held,
        "max_num_bits":maxnb,
        "max_den_bits":maxdb,
    },
},indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
