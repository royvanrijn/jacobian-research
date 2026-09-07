#!/usr/bin/env sage -python
"""
Hybrid reconstruction of the q32 horizontal section from:
  - the checkpointed p-adic lift modulo 100003^N, and
  - all independent modular section artifacts.

This reconstructs each coefficient independently.  It therefore avoids the
possibly huge common-denominator vector used by simultaneous LLL.
"""

import json
from pathlib import Path
from sage.all import PolynomialRing, QQ, ZZ, Zmod

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

base=ZZ(100003)
CP=LOCAL/f"q32-horizontal-hensel-checkpoint-p{base}.json"
if not CP.exists():
    raise SystemExit(f"missing checkpoint {CP}")

cp=json.loads(CP.read_text())
assert cp["status"]=="PASS_HENSEL_CHECKPOINT"
digit=int(cp["digit"])
M0=ZZ(cp["modulus"])
assert M0==base**digit

z0=[ZZ(v)%M0 for v in cp["z"]]+[ZZ(1)]
x0=[ZZ(v)%M0 for v in cp["x"]]
y0=[ZZ(v)%M0 for v in cp["y"]]
assert (len(z0),len(x0),len(y0))==(25,53,79)

records=[]
for path in sorted(LOCAL.glob("q24-degree46-direct-global-mod-*.json")):
    try:
        d=json.loads(path.read_text())
        p=ZZ(d["prime"])
    except Exception:
        continue
    if p==base or d.get("status")!="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE":
        continue
    sec=d["section_mod_p"]
    prof=sec.get("profile",{})
    if (prof.get("Z_degree"),prof.get("X_degree"),prof.get("Y_degree"))!=(24,52,78):
        continue
    records.append((p,d,path))

records.sort(key=lambda t:int(t[0]))
if len(records)<2:
    raise SystemExit("need >=2 independent modular primes")

# Final independent prime is not used for CRT.
holdp,hold,_=records[-1]
train=records[:-1]

print(
    "Q32PHYBRID_INPUT|"
    f"hensel_digit={digit}|hensel_bits={M0.nbits()}|"
    f"train_primes={','.join(str(p) for p,_,_ in train)}|"
    f"holdout={holdp}|status=PASS",
    flush=True,
)

def crt_pair(x,M,r,p):
    x=ZZ(x)%M
    r=ZZ(r)%p
    t=((r-x)%p)*((M%p).inverse_mod(p))%p
    return (x+M*t)%(M*p),M*p

def rr(x,M):
    try:
        return QQ(Zmod(M)(ZZ(x)%M).rational_reconstruction())
    except Exception:
        return None

def red(q,p):
    q=QQ(q); p=ZZ(p)
    d=ZZ(q.denominator())%p
    if not d:
        return None
    return int((ZZ(q.numerator())%p)*d.inverse_mod(p)%p)

def arr(d,key,n):
    vals=[ZZ(v) for v in d["section_mod_p"][key]]
    vals += [ZZ(0)]*(n-len(vals))
    assert len(vals)==n
    return vals

objects=[
    ("Z","Z_coefficients_low_to_high",z0,25),
    ("X","X_coefficients_low_to_high",x0,53),
    ("Y","Y_coefficients_low_to_high",y0,79),
]

results={}
complete=True

for name,key,basevals,n in objects:
    train_arrays=[arr(d,key,n) for _,d,_ in train]
    holdvals=arr(hold,key,n)

    qs=[]
    recovered=held=0
    maxnb=maxdb=0
    failed_idx=[]
    held_idx=[]

    for j in range(n):
        x=ZZ(basevals[j])
        M=M0
        for (p,unused_d,unused_path),aa in zip(train,train_arrays):
            x,M=crt_pair(x,M,aa[j],p)

        q=rr(x,M)
        qs.append(q)
        if q is None:
            failed_idx.append(j)
            continue
        recovered+=1
        maxnb=max(maxnb,abs(ZZ(q.numerator())).nbits())
        maxdb=max(maxdb,abs(ZZ(q.denominator())).nbits())
        if red(q,holdp)==int(holdvals[j]%holdp):
            held+=1
            held_idx.append(j)

    obj_complete=(held==n)
    complete &= obj_complete
    results[name]=qs
    print(
        "Q32PHYBRID_PROGRESS|"
        f"object={name}|combined_bits={M.nbits()}|"
        f"recovered={recovered}/{n}|heldout={held}/{n}|"
        f"max_num_bits={maxnb}|max_den_bits={maxdb}|"
        f"failed_indices={','.join(map(str,failed_idx[:16])) or 'none'}|"
        f"status={'PASS_HELDOUT' if obj_complete else 'PARTIAL'}",
        flush=True,
    )

# Exact QQ identity if complete.
identity=False
matches=[]
if complete:
    R=PolynomialRing(QQ,"U")
    Z=R(results["Z"])
    X=R(results["X"])
    Y=R(results["Y"])

    q8_candidates=[
        LOCAL/"q8-corrected2cover-qq-child.json",
        GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
    ]
    Q8=next((
        q for q in q8_candidates
        if q.exists()
        and json.loads(q.read_text()).get("status")=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"
    ),None)
    if Q8 is None:
        raise SystemExit("missing exact D13 parent")
    q8=json.loads(Q8.read_text())
    A=R([QQ(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]])
    B=R([QQ(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]])
    identity=(Y**2==X**3+A*X*Z**4+B*Z**6)

    # Check every unused/independent modular artifact, not only the holdout.
    for p,d,_ in records:
        ok=True
        for poly,key,n in (
            (Z,"Z_coefficients_low_to_high",25),
            (X,"X_coefficients_low_to_high",53),
            (Y,"Y_coefficients_low_to_high",79),
        ):
            target=arr(d,key,n)
            for j,q in enumerate(poly.list()+[QQ(0)]*(n-len(poly.list()))):
                if red(q,p)!=int(target[j]%p):
                    ok=False
                    break
            if not ok:
                break
        if ok:
            matches.append(int(p))

    print(
        "Q32PHYBRID_EXACT|"
        f"identity={int(identity)}|"
        f"independent_matches={','.join(map(str,matches)) or 'none'}|"
        f"status={'PASS_EXACT_AND_INDEPENDENT' if identity and matches else 'PARTIAL'}",
        flush=True,
    )

out=LOCAL/"q32-horizontal-section-hybrid-crt-qq.json"
payload={
    "schema":"elkies-k3.h3-q32-horizontal-section-hybrid-crt.v1",
    "status":(
        "PASS_EXACT_Q32_HORIZONTAL_SECTION_HYBRID"
        if complete and identity and matches else
        "PARTIAL_Q32_HORIZONTAL_SECTION_HYBRID"
    ),
    "base_prime":int(base),
    "hensel_digit":digit,
    "hensel_precision_bits":int(M0.nbits()),
    "training_primes":[int(p) for p,_,_ in train],
    "heldout_prime":int(holdp),
}
if complete:
    payload["section"]={
        "Z_coefficients_low_to_high":[str(q) for q in results["Z"]],
        "X_coefficients_low_to_high":[str(q) for q in results["X"]],
        "Y_coefficients_low_to_high":[str(q) for q in results["Y"]],
    }
    payload["exact_weierstrass_identity"]=bool(identity)
    payload["independent_matches"]=matches

out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
print(
    "Q32PHYBRID_RESULT|"
    f"hensel_digit={digit}|status={payload['status']}",
    flush=True,
)
