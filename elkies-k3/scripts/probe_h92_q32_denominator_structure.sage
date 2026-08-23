#!/usr/bin/env sage -python
import json
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

pointed=json.loads((LOCAL/"q32-pointed-spinor-weierstrass-anchor.json").read_text())
pby={int(r["prime"]):r for r in pointed["primes"]}

rows=[]
for spath in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        sig=json.loads(spath.read_text())
        p=int(sig["prime"])
    except Exception:
        continue
    if sig.get("status")!="PASS_Q32_MODP_SIGNATURE" or p not in pby:
        continue
    rows.append((p,sig,pby[p]))

for p,sig,pt in rows:
    F=GF(p)
    R=PolynomialRing(F,"V")
    V=R.gen()

    def poly(rec, which):
        return R([F(v) for v in rec[which]])

    A=sig["jacobian_A"]; B=sig["jacobian_B"]; X=pt["marked_section_x"]
    Ad=poly(A,"den"); Bd=poly(B,"den"); Xd=poly(X,"den")
    An=poly(A,"num"); Bn=poly(B,"num"); Xn=poly(X,"num")

    # monic normalization
    Ad/=Ad.leading_coefficient()
    Bd/=Bd.leading_coefficient()
    Xd/=Xd.leading_coefficient()

    gAB=Ad.gcd(Bd).monic()
    gAX=Ad.gcd(Xd).monic()
    gBX=Bd.gcd(Xd).monic()

    # Test for one common polynomial h with denominator exponents 4,6,2:
    # A ~ A0/h^4, B ~ B0/h^6, x ~ X0/h^2.
    def root_if_exact(P,e):
        out=R.one()
        for f,m in P.factor():
            if int(m)%e:
                return None
            out*=f.monic()**(int(m)//e)
        return out.monic()

    hA=root_if_exact(Ad,4)
    hB=root_if_exact(Bd,6)
    hX=root_if_exact(Xd,2)
    common = (hA is not None and hB is not None and hX is not None and hA==hB==hX)

    # Also test weaker proportional power identities directly.
    powAB = (Ad**3 == Bd**2)
    powAX = (Ad == Xd**2)
    powBX = (Bd == Xd**3)

    print(
        "Q32DEN|"
        f"prime={p}|"
        f"degA={Ad.degree()}|degB={Bd.degree()}|degX={Xd.degree()}|"
        f"gAB={gAB.degree()}|gAX={gAX.degree()}|gBX={gBX.degree()}|"
        f"A3eqB2={int(powAB)}|AeqX2={int(powAX)}|BeqX3={int(powBX)}|"
        f"common_h={int(common)}|"
        f"hdeg={(hA.degree() if common else -1)}|"
        "status=PASS",
        flush=True,
    )

    if common:
        h=hA
        A0=An*(h**4)/Ad
        B0=Bn*(h**6)/Bd
        X0=Xn*(h**2)/Xd
        assert A0 in R and B0 in R and X0 in R
        print(
            "Q32POLYMODEL|"
            f"prime={p}|hdeg={h.degree()}|"
            f"A0deg={R(A0).degree()}|B0deg={R(B0).degree()}|X0deg={R(X0).degree()}|"
            f"A0terms={len([c for c in R(A0) if c])}|"
            f"B0terms={len([c for c in R(B0) if c])}|"
            f"X0terms={len([c for c in R(X0) if c])}|"
            "status=PASS_COMMON_DENOMINATOR_MODEL",
            flush=True,
        )

print(f"Q32DEN_RESULT|primes={len(rows)}|status=PASS",flush=True)
