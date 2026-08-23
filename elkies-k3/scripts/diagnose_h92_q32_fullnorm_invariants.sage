#!/usr/bin/env sage -python
"""
Diagnose what remains incoherent after the q32 D12 full normalization.

Rebuild the same normalization:
  * unique I8* -> infinity
  * center finite discriminant roots
  * q=b9/(a6*x4) fixes base scale
  * r=a6^2/x4^3 fixes short-Weierstrass scale

Then report the actual normalized coefficient vectors prime-by-prime and test
three layers independently with a held-out prime:

  X : marked spinor section
  j-data : scale-free c4^3/Delta as a rational function
  Dshape : monic centered degree-10 discriminant

If Dshape/j-data reconstruct but X does not, the curve is coherent and only
the marking is drifting. If Dshape itself does not reconstruct, the q32 base
parameter still has an unremoved projective ambiguity.
"""

import json, math
from pathlib import Path
from sage.all import GF, Integers, PolynomialRing, QQ, ZZ, binomial

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

pointed=json.loads((LOCAL/"q32-pointed-spinor-weierstrass-anchor.json").read_text())
assert pointed["status"]=="PASS_POINTED_Q32_D12_SPINOR_MARKING"
pby={int(r["prime"]):r for r in pointed["primes"]}

records=[]
for spath in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        sig=json.loads(spath.read_text()); p=int(sig["prime"])
    except Exception:
        continue
    if sig.get("status")!="PASS_Q32_MODP_SIGNATURE" or p not in pby:
        continue
    records.append((ZZ(p),sig,pby[p]))
records.sort(key=lambda r:int(r[0]))


def norm(row):
    p,sig,pt=row
    F=GF(p)
    R=PolynomialRing(F,"S")
    S=R.gen()

    def parts(rec):
        n=R([F(v) for v in rec["num"]])
        d=R([F(v) for v in rec["den"]])
        lc=d.leading_coefficient()
        return n/lc,d/lc

    A0,Ad=parts(sig["jacobian_A"])
    B0,Bd=parts(sig["jacobian_B"])
    X0,Xd=parts(pt["marked_section_x"])

    def root(P,e):
        out=R.one()
        for f,m in R(P).factor():
            m=int(m)
            if m%e:return None
            out*=f.monic()**(m//e)
        return out.monic()

    h=root(Ad,4)
    assert h is not None and h==root(Bd,6)==root(Xd,2)
    assert Ad==h**4 and Bd==h**6 and Xd==h**2

    D0=-16*(4*A0**3+27*B0**2)
    stars=[f.monic() for f,e in D0.factor() if int(e)==14 and f.degree()==1]
    assert len(stars)==1
    beta=-stars[0][0]

    def i8(P,w):
        P=R(P); out=R.zero()
        for i,c in enumerate(P.list()):
            for j in range(i+1):
                out += c*F(binomial(i,j))*beta**(i-j)*S**(w-j)
        return R(out)

    A=i8(A0,8); B=i8(B0,12); X=i8(X0,4)
    assert (A.degree(),B.degree(),X.degree())==(6,9,4)

    D=-16*(4*A**3+27*B**2)
    assert D.degree()==10
    mu=-D[9]/(F(10)*D[10])

    def shift(P,z):
        P=R(P); out=R.zero()
        for i,c in enumerate(P.list()):
            for j in range(i+1):
                out += c*F(binomial(i,j))*z**(i-j)*S**j
        return R(out)

    A=shift(A,mu); B=shift(B,mu); X=shift(X,mu)
    D=-16*(4*A**3+27*B**2)
    assert D[9]==0

    a6,b9,x4=A[6],B[9],X[4]
    q=b9/(a6*x4)
    r=a6**2/x4**3
    assert q and r

    def bscale(P,z):
        return R([c*z**i for i,c in enumerate(R(P).list())])

    A=bscale(A,q)/(r**2)
    B=bscale(B,q)/(r**3)
    X=bscale(X,q)/r
    D=-16*(4*A**3+27*B**2)
    assert D.degree()==10 and D[9]==0

    # Monic discriminant shape removes any residual overall equation scalar.
    Dm=D/D[10]

    # c4^3/Delta = (-48 A)^3/Delta, represented canonically by monic denom.
    c4=-48*A
    jnum=c4**3
    jden=D
    g=jnum.gcd(jden)
    jnum//=g; jden//=g
    lc=jden.leading_coefficient()
    jnum/=lc; jden/=lc

    return {
        "X":[int(X[i]) for i in range(5)],
        "A":[int(A[i]) for i in range(7)],
        "B":[int(B[i]) for i in range(10)],
        "Dshape":[int(Dm[i]) for i in range(10)], # leading 1 omitted; d9=0 included
        "jnum":[int(c) for c in jnum.list()],
        "jden":[int(c) for c in jden.list()],
        "lead":[int(A[6]),int(B[9]),int(X[4])],
    }


data=[(p,norm(r)) for p,*_ in records for r in [next(rr for rr in records if rr[0]==p)]]
# Simpler deterministic rebuild without clever comprehension.
data=[]
for row in records:
    data.append((row[0],norm(row)))

for p,d in data:
    print(
        "Q32FULLNORM_VECTOR|"
        f"prime={p}|"
        f"X={','.join(map(str,d['X']))}|"
        f"lead={','.join(map(str,d['lead']))}|"
        f"D0to4={','.join(map(str,d['Dshape'][:5]))}|"
        "status=PASS",
        flush=True,
    )

train=data[:-1]
holdp,hold=data[-1]
primes=[p for p,_ in train]


def crt(vals):
    x=ZZ(0); M=ZZ(1)
    for a,p in zip(vals,primes):
        a=ZZ(a)%p
        t=((a-x)%p)*((M%p).inverse_mod(p))%p
        x=(x+M*t)%(M*p); M*=p
    return x,M


def rr_test(name):
    arrs=[d[name] for _,d in train]
    hv=hold[name]
    if not all(len(a)==len(arrs[0]) for a in arrs) or len(hv)!=len(arrs[0]):
        print(f"Q32FULLNORM_DIAG|object={name}|status=SHAPE_MISMATCH",flush=True)
        return
    got=valid=0
    rows=[]
    for j in range(len(arrs[0])):
        x,M=crt([a[j] for a in arrs])
        try:
            q=QQ(Integers(M)(x).rational_reconstruction())
        except Exception:
            rows.append((j,None,False))
            continue
        got+=1
        den=ZZ(q.denominator())%holdp
        ok=bool(den and
            int((ZZ(q.numerator())%holdp)*den.inverse_mod(holdp)%holdp)==
            int(hv[j]%holdp))
        valid+=int(ok)
        rows.append((j,q,ok))
    print(
        f"Q32FULLNORM_DIAG|object={name}|recovered={got}/{len(arrs[0])}|"
        f"heldout={valid}/{len(arrs[0])}|modulus_bits={M.nbits()}|"
        f"status={'PASS_HELDOUT' if valid==len(arrs[0]) else 'PARTIAL'}",
        flush=True,
    )
    for j,q,ok in rows:
        if ok or (q is not None and abs(ZZ(q.numerator())).nbits()<32 and abs(ZZ(q.denominator())).nbits()<32):
            print(
                f"Q32FULLNORM_COEFF|object={name}|index={j}|value={q}|"
                f"holdout={int(ok)}",
                flush=True,
            )

for name in ("X","Dshape","A","B"):
    rr_test(name)

# j numerator/denominator degree may be stable; test if so.
for key in ("jnum","jden"):
    lengths=sorted(set(len(d[key]) for _,d in data))
    print(f"Q32FULLNORM_SHAPE|object={key}|lengths={lengths}",flush=True)
    if len(lengths)==1:
        rr_test(key)

print(f"Q32FULLNORM_DIAG_RESULT|primes={len(data)}|holdout={holdp}|status=PASS",flush=True)
