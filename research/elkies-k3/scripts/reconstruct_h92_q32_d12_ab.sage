#!/usr/bin/env sage -python
"""
Recover the q32 D12 child from the canonical j-map cube/square roots.

After canonicalizing the child base:
    j = P18/Q10,
    P = c*a6^3,
    P-1728Q = c*b9^2
with monic a6,b9 and the SAME scalar c. Hence
    j = 1728*a^3/(a^3-b^2)
and a canonical model is
    y^2 = x^3 - 3*a*x + 2*b,
    Delta = 1728*(a^3-b^2).

Thus only the monic degree-6 and degree-9 polynomials a,b need reconstruction.
"""

import json
from pathlib import Path
from sage.all import GF, PolynomialRing, QQ, ZZ, Zmod

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

records=[]
for path in sorted(LOCAL.glob("q32-signature-mod-*.json")):
    try:
        d=json.loads(path.read_text())
        p=ZZ(d["prime"])
    except Exception:
        continue
    if d.get("status")=="PASS_Q32_MODP_SIGNATURE":
        records.append((p,d,path))

records.sort(key=lambda t:(0 if t[0]==100003 else 1,int(t[0])))
if len(records)<3:
    raise SystemExit("need >=3 q32 signatures")

def canonical_roots(p,d):
    F=GF(p)
    RV=PolynomialRing(F,"V"); V=RV.gen(); KV=RV.fraction_field()

    def rf(rec):
        return KV(RV([F(v) for v in rec["num"]]))/KV(
            RV([F(v) for v in rec["den"]])
        )

    A=rf(d["jacobian_A"])
    B=rf(d["jacobian_B"])
    J=KV(F(6912))*A**3/(KV(F(4))*A**3+KV(F(27))*B**2)
    N=RV(J.numerator()); D=RV(J.denominator())
    assert N.gcd(D).degree()==0 and N.degree()==D.degree()==18

    e8=[f for f,e in D.factor() if int(e)==8 and f.degree()==1]
    if len(e8)!=1:
        return None
    f=e8[0]
    r=-f[0]/f[1]

    RT=PolynomialRing(F,"T"); T=RT.gen()
    def invpoly(poly):
        return sum(
            F(poly[i])*(F(r)*T+1)**i*T**(18-i)
            for i in range(poly.degree()+1)
        )

    P=invpoly(N); Q=invpoly(D)
    g=P.gcd(Q)
    if g.degree()>0:
        P//=g; Q//=g
    if P.degree()!=18 or Q.degree()!=10:
        return None

    lc=Q.leading_coefficient()
    P/=lc; Q/=lc
    if p in (2,5):
        return None

    # Center q9=0.
    c=-Q[9]/F(10)
    RS=PolynomialRing(F,"S"); S=RS.gen()
    P1=RS(P(S+c)); Q1=RS(Q(S+c))
    assert Q1.is_monic() and Q1[9]==0

    # Scale q8=q7.
    if not Q1[8] or not Q1[7]:
        return None
    scale=Q1[7]/Q1[8]
    RU=PolynomialRing(F,"U"); U=RU.gen()
    P2=RU(P1(scale*U)); Q2=RU(Q1(scale*U))
    lc2=Q2.leading_coefficient()
    P2/=lc2; Q2/=lc2
    assert Q2.is_monic() and Q2[9]==0 and Q2[8]==Q2[7]

    def monic_power_root(poly,e):
        out=RU.one()
        for fac,m in poly.factor():
            if int(m)%e:
                return None
            out*=fac.monic()**(int(m)//e)
        return out.monic()

    a=monic_power_root(P2,3)
    H=P2-F(1728)*Q2
    b=monic_power_root(H,2)
    if a is None or b is None or a.degree()!=6 or b.degree()!=9:
        return None

    cp=P2.leading_coefficient()
    cb=H.leading_coefficient()
    assert cp==cb
    assert P2==cp*a**3
    assert H==cp*b**2

    K=a**3-b**2
    assert K.degree()==10
    assert cp*K==F(1728)*Q2
    assert K[9]==0 and K[8]==K[7]

    return {
        "a":a, "b":b, "K":K,
        "P":P2, "Q":Q2,
        "scalar":cp,
    }

good=[]
for p,d,path in records:
    r=canonical_roots(p,d)
    if r is None:
        print(f"Q32AB_PRIME|prime={p}|status=SKIP",flush=True)
        continue
    print(
        "Q32AB_PRIME|"
        f"prime={p}|adeg={r['a'].degree()}|bdeg={r['b'].degree()}|"
        f"Kdeg={r['K'].degree()}|K9={int(r['K'][9])}|"
        f"K8eqK7={int(r['K'][8]==r['K'][7])}|status=PASS",
        flush=True,
    )
    good.append((p,r))

if len(good)<3:
    raise SystemExit("too few good canonical a,b pairs")

print(
    "Q32AB_GROUP|"
    f"compatible={len(good)}/{len(records)}|"
    f"primes={','.join(str(p) for p,_ in good)}|"
    "profile=a6_monic,b9_monic,K10|status=PASS",
    flush=True,
)

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

train=good[:-1]
holdp,hold=good[-1]
mods=[p for p,_ in train]

results={}
complete=True
for name,n in (("a",7),("b",10)):
    arrays=[[ZZ(r[name][i]) for i in range(n)] for _,r in train]
    hv=[ZZ(hold[name][i]) for i in range(n)]

    vals=[]
    recovered=held=stable=stableheld=0
    maxnb=maxdb=0
    details=[]

    for j in range(n):
        q,M=rr_scalar([x[j] for x in arrays],mods)
        vals.append(q)
        ok=False
        if q is not None:
            recovered+=1
            maxnb=max(maxnb,abs(ZZ(q.numerator())).nbits())
            maxdb=max(maxdb,abs(ZZ(q.denominator())).nbits())
            ok=(red(q,holdp)==int(hv[j]%holdp))
            if ok: held+=1

        qprev=None
        if len(train)>=3:
            qprev,_=rr_scalar([x[j] for x in arrays[:-1]],mods[:-1])
        st=(q is not None and qprev is not None and q==qprev)
        if st:
            stable+=1
            if ok: stableheld+=1

        details.append((j,q,ok,st))

    objok=(held==n)
    complete &= objok
    results[name]=vals
    print(
        "Q32AB_CRT|"
        f"object={name}|train={len(train)}|holdout={holdp}|"
        f"recovered={recovered}/{n}|heldout={held}/{n}|"
        f"stable={stable}|stable_heldout={stableheld}|"
        f"max_num_bits={maxnb}|max_den_bits={maxdb}|"
        f"status={'PASS_HELDOUT' if objok else 'PARTIAL'}",
        flush=True,
    )
    for j,q,ok,st in details:
        if ok or st:
            print(
                "Q32AB_COEFF|"
                f"object={name}|i={j}|value={q}|"
                f"heldout={int(ok)}|stable={int(st)}|status=SIGNAL",
                flush=True,
            )

exact=False
if complete:
    R=PolynomialRing(QQ,"u"); u=R.gen()
    a=R(results["a"]); b=R(results["b"])
    assert a.degree()==6 and b.degree()==9 and a.is_monic() and b.is_monic()
    K=a**3-b**2
    exact=(
        K.degree()==10
        and K[9]==0
        and K[8]==K[7]
    )
    A=-3*a
    B=2*b
    Delta=-16*(4*A**3+27*B**2)
    assert Delta==1728*K

    print(
        "Q32AB_EXACT|"
        f"Kdeg={K.degree()}|K9={K[9]}|K8eqK7={int(K[8]==K[7])}|"
        f"Adeg={A.degree()}|Bdeg={B.degree()}|Deltadeg={Delta.degree()}|"
        f"status={'PASS_D12_NORMAL_FORM' if exact else 'FAIL'}",
        flush=True,
    )

out=LOCAL/"q32-d12-ab-crt.json"
payload={
    "schema":"elkies-k3.h3-q32-d12-ab-crt.v1",
    "status":(
        "PASS_EXACT_Q32_D12_AB_NORMAL_FORM"
        if complete and exact else
        "PARTIAL_Q32_D12_AB_CRT"
    ),
    "primes":[int(p) for p,_ in good],
    "training_primes":[int(p) for p,_ in train],
    "heldout_prime":int(holdp),
    "normal_form":"y^2=x^3-3*a(u)*x+2*b(u)",
}
if complete:
    payload["a_coefficients_low_to_high"]=[str(q) for q in results["a"]]
    payload["b_coefficients_low_to_high"]=[str(q) for q in results["b"]]

out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
print(
    "Q32AB_RESULT|"
    f"good={len(good)}|train={len(train)}|holdout={holdp}|"
    f"status={payload['status']}",
    flush=True,
)
