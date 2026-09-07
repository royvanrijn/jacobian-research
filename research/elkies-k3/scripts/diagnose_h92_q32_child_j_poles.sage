#!/usr/bin/env sage -python
import json
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

records=[]
for path in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        d=json.loads(path.read_text())
        p=ZZ(d["prime"])
    except Exception:
        continue
    if d.get("status")!="PASS_Q32_MODP_SIGNATURE":
        continue
    records.append((p,d,path))

records.sort(key=lambda t:(0 if t[0]==100003 else 1,int(t[0])))
if not records:
    raise SystemExit("no q32 signatures")

def rat(R,K,rec):
    n=R([K.base_ring()(v) for v in rec["num"]])
    d=R([K.base_ring()(v) for v in rec["den"]])
    return K(n)/K(d)

good=[]
for p,d,path in records:
    F=GF(p)
    R=PolynomialRing(F,"V")
    V=R.gen()
    K=R.fraction_field()

    A=rat(R,K,d["jacobian_A"])
    B=rat(R,K,d["jacobian_B"])
    J=K(ZZ(6912))*A**3/(K(4)*A**3+K(27)*B**2)
    J=K(J)

    den=R(J.denominator())
    num=R(J.numerator())
    fac=den.factor()
    profile=sorted((f.degree(),int(e)) for f,e in fac)

    order8=[f for f,e in fac if int(e)==8 and f.degree()==1]
    simple_linear=[f for f,e in fac if int(e)==1 and f.degree()==1]
    simple_other=[(f.degree(),int(e)) for f,e in fac if not (int(e)==8 and f.degree()==1) and not (int(e)==1 and f.degree()==1)]

    roots8=[]
    roots1=[]
    for f in order8:
        roots8.append(int(-f[0]/f[1]))
    for f in simple_linear:
        roots1.append(int(-f[0]/f[1]))

    ok=(len(order8)==1 and len(simple_linear)==2 and simple_other==[(8,1)])
    print(
        "Q32JPOLES|"
        f"prime={p}|j={num.degree()}/{den.degree()}|"
        f"profile={profile}|I8star_roots={roots8}|"
        f"rational_I1={roots1}|status={'PASS_CANONICAL_TRIPLE' if ok else 'OTHER'}",
        flush=True,
    )
    if ok:
        good.append((p,d,roots8[0],roots1))

print(
    "Q32JPOLES_RESULT|"
    f"canonical={len(good)}/{len(records)}|"
    f"primes={','.join(str(x[0]) for x in good)}|"
    f"status={'PASS_STABLE_CHILD_NORMALIZATION' if len(good)>=3 else 'INSUFFICIENT'}",
    flush=True,
)
