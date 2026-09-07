#!/usr/bin/env sage -python
import json
from pathlib import Path
from sage.all import PolynomialRing, QQ, ZZ, Zmod

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

# ---------------------------------------------------------------------------
# Collect every valid modular horizontal section, including primes that were
# unsuitable for the later split-spinor q32 divisorial calculation.
# ---------------------------------------------------------------------------
records=[]
for path in sorted(LOCAL.glob("q24-degree46-direct-global-mod-*.json")):
    try:
        d=json.loads(path.read_text())
        p=ZZ(d["prime"])
        sec=d["section_mod_p"]
    except Exception:
        continue
    if d.get("status")!="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE":
        continue
    prof=sec.get("profile",{})
    if (
        prof.get("Z_degree")!=24
        or prof.get("X_degree")!=52
        or prof.get("Y_degree")!=78
    ):
        continue
    records.append((p,d,path))

records.sort(key=lambda t:(0 if t[0]==100003 else 1,int(t[0])))
if len(records)<3:
    raise SystemExit(f"need >=3 modular horizontal sections, found {len(records)}")

print(
    "Q32PCRT_INPUT|"
    f"count={len(records)}|primes={','.join(str(p) for p,_,_ in records)}|"
    "profile=Z24,X52,Y78|status=PASS",
    flush=True,
)

# Last prime is entirely held out.
train=records[:-1]
holdp,hold,_=records[-1]
mods=[p for p,_,_ in train]

def coeffs(rec,key,n):
    vals=[ZZ(v) for v in rec["section_mod_p"][key]]
    vals += [ZZ(0)]*(n-len(vals))
    assert len(vals)==n
    return vals

def crt_scalar(vals,mods):
    x=ZZ(0); M=ZZ(1)
    for rr,p in zip(vals,mods):
        rr=ZZ(rr)%p
        t=((rr-x)%p)*((M%p).inverse_mod(p))%p
        x=(x+M*t)%(M*p)
        M*=p
    return x,M

def rr_scalar(vals,mods):
    x,M=crt_scalar(vals,mods)
    try:
        return QQ(Zmod(M)(x).rational_reconstruction()),M
    except Exception:
        return None,M

def red(q,p):
    q=QQ(q); p=ZZ(p)
    den=ZZ(q.denominator())%p
    if not den:
        return None
    return int((ZZ(q.numerator())%p)*den.inverse_mod(p)%p)

objects=[
    ("Z","Z_coefficients_low_to_high",25),
    ("X","X_coefficients_low_to_high",53),
    ("Y","Y_coefficients_low_to_high",79),
]

results={}
all_complete=True
for name,key,n in objects:
    arrays=[coeffs(d,key,n) for _,d,_ in train]
    hv=coeffs(hold,key,n)

    out=[]
    recovered=held=0
    nonzero_held=0
    nonzero_total=sum(1 for x in hv if x%holdp)
    maxnb=maxdb=0

    for j in range(n):
        q,M=rr_scalar([a[j] for a in arrays],mods)
        out.append(q)
        if q is None:
            continue
        recovered+=1
        maxnb=max(maxnb,abs(ZZ(q.numerator())).nbits())
        maxdb=max(maxdb,abs(ZZ(q.denominator())).nbits())
        ok=(red(q,holdp)==int(hv[j]%holdp))
        if ok:
            held+=1
            if hv[j]%holdp:
                nonzero_held+=1

    complete=(held==n)
    all_complete &= complete
    results[name]=out

    print(
        "Q32PCRT_PROGRESS|"
        f"object={name}|train={len(train)}|holdout={holdp}|"
        f"recovered={recovered}/{n}|heldout={held}/{n}|"
        f"nonzero_heldout={nonzero_held}/{nonzero_total}|"
        f"max_num_bits={maxnb}|max_den_bits={maxdb}|"
        f"status={'PASS_HELDOUT' if complete else 'PARTIAL'}",
        flush=True,
    )

# ---------------------------------------------------------------------------
# If every coefficient reconstructs and predicts the unseen prime, verify the
# section exactly on the characteristic-zero D13 parent.
# ---------------------------------------------------------------------------
identity=False
if all_complete:
    R=PolynomialRing(QQ,"U")
    Z=R(results["Z"])
    X=R(results["X"])
    Y=R(results["Y"])
    assert (Z.degree(),X.degree(),Y.degree())==(24,52,78)
    assert Z.is_monic()

    q8_candidates=[
        LOCAL/"q8-corrected2cover-qq-child.json",
        GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
    ]
    Q8=next((
        p for p in q8_candidates
        if p.exists()
        and json.loads(p.read_text()).get("status")=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"
    ),None)
    if Q8 is None:
        raise SystemExit("missing exact q8/D13 parent artifact")

    q8=json.loads(Q8.read_text())
    A=R([QQ(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]])
    B=R([QQ(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]])

    identity=(Y**2 == X**3 + A*X*Z**4 + B*Z**6)
    print(
        "Q32PCRT_EXACT|"
        f"Z={Z.degree()}|X={X.degree()}|Y={Y.degree()}|"
        f"weierstrass_identity={int(identity)}|"
        f"status={'PASS_EXACT_SECTION' if identity else 'FAIL_IDENTITY'}",
        flush=True,
    )

out=LOCAL/"q32-horizontal-section-crt-qq.json"
payload={
    "schema":"elkies-k3.h3-q32-horizontal-section-crt-qq.v1",
    "status":(
        "PASS_EXACT_Q32_HORIZONTAL_SECTION_FROM_MODULAR_CRT"
        if all_complete and identity else
        "PARTIAL_Q32_HORIZONTAL_SECTION_CRT"
    ),
    "training_primes":[int(p) for p,_,_ in train],
    "heldout_prime":int(holdp),
    "available_primes":[int(p) for p,_,_ in records],
    "profile":{"Z_degree":24,"X_degree":52,"Y_degree":78},
}
if all_complete:
    payload["section"]={
        "Z_coefficients_low_to_high":[str(q) for q in results["Z"]],
        "X_coefficients_low_to_high":[str(q) for q in results["X"]],
        "Y_coefficients_low_to_high":[str(q) for q in results["Y"]],
    }
    payload["exact_weierstrass_identity"]=bool(identity)

out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
print(
    "Q32PCRT_RESULT|"
    f"available={len(records)}|train={len(train)}|holdout={holdp}|"
    f"status={payload['status']}",
    flush=True,
)
