#!/usr/bin/env sage -python
"""
Canonicalize the modular q32/D12 j-map using only the unique I8* pole and the
degree-10 simple-pole divisor.

For each prime:
  1. find the unique linear denominator factor of multiplicity 8;
  2. send it to infinity via t = 1/(V-r);
  3. obtain j=P18/Q10 and make Q10 monic;
  4. translate t=s+c so coeff(s^9,Q)=0;
  5. scale s=a*u with a=q7/q8, so the monic denominator satisfies q8=q7.

Steps 4-5 canonically remove the full affine freedom remaining after the
I8* pole is at infinity, using rational operations only (no chosen I1 roots).

Then CRT/rational-reconstruct this canonical j-map across primes and hold one
prime out entirely.
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
    if d.get("status")!="PASS_Q32_MODP_SIGNATURE":
        continue
    records.append((p,d,path))

records.sort(key=lambda t:(0 if t[0]==100003 else 1,int(t[0])))
if len(records)<3:
    raise SystemExit("need >=3 q32 signatures")

def canon_j(p,d):
    F=GF(p)
    RV=PolynomialRing(F,"V")
    V=RV.gen()
    KV=RV.fraction_field()

    def rf(rec):
        n=RV([F(v) for v in rec["num"]])
        den=RV([F(v) for v in rec["den"]])
        return KV(n)/KV(den)

    A=rf(d["jacobian_A"])
    B=rf(d["jacobian_B"])
    J=KV(F(6912))*A**3/(KV(F(4))*A**3+KV(F(27))*B**2)
    N=RV(J.numerator())
    D=RV(J.denominator())
    assert N.gcd(D).degree()==0
    assert N.degree()==D.degree()==18

    fac=D.factor()
    e8=[f for f,e in fac if int(e)==8 and f.degree()==1]
    if len(e8)!=1:
        return None, f"I8star_count_{len(e8)}"
    f=e8[0]
    r=-f[0]/f[1]

    RT=PolynomialRing(F,"T")
    T=RT.gen()

    # T^18 * poly(r+1/T)
    def invert_poly(poly):
        return sum(
            F(poly[i])*(F(r)*T+1)**i*T**(18-i)
            for i in range(poly.degree()+1)
        )

    P=invert_poly(N)
    Q=invert_poly(D)
    g=P.gcd(Q)
    if g.degree()>0:
        P//=g
        Q//=g

    if P.degree()!=18 or Q.degree()!=10:
        return None, f"inverted_degrees_{P.degree()}_{Q.degree()}"

    # Normalize denominator monic.
    lc=Q.leading_coefficient()
    P/=lc
    Q/=lc
    assert Q.is_monic()

    # Kill t^9 coefficient: t = s+c, c=-q9/10.
    if p in (2,5):
        return None,"bad_characteristic_for_centering"
    c=-Q[9]/F(10)

    RS=PolynomialRing(F,"S")
    S=RS.gen()
    P1=RS(P(S+c))
    Q1=RS(Q(S+c))
    assert Q1.degree()==10 and Q1.is_monic() and Q1[9]==0

    q8=Q1[8]
    q7=Q1[7]
    if not q8 or not q7:
        return None,"canonical_scale_coefficient_zero"

    # s=a*u, a=q7/q8 => after monic renormalization q8'=q7'.
    a=q7/q8

    RU=PolynomialRing(F,"U")
    U=RU.gen()
    P2=RU(P1(a*U))
    Q2=RU(Q1(a*U))
    lc2=Q2.leading_coefficient()
    P2/=lc2
    Q2/=lc2

    assert Q2.degree()==10 and Q2.is_monic()
    assert Q2[9]==0
    assert Q2[8]==Q2[7]

    # Intrinsic Kodaira signatures visible in j:
    # numerator should be constant * cube (j=0 ramification multiples of 3),
    # P-1728Q constant * square (j=1728 ramification multiples of 2).
    cube_ok=all(int(e)%3==0 for unused_f,e in P2.factor())
    H=P2-F(1728)*Q2
    square_ok=all(int(e)%2==0 for unused_f,e in H.factor())

    return {
        "P":P2,
        "Q":Q2,
        "r":int(r),
        "center":int(c),
        "scale":int(a),
        "cube_ok":cube_ok,
        "square_ok":square_ok,
    },None

good=[]
for p,d,path in records:
    res,reason=canon_j(p,d)
    if res is None:
        print(
            f"Q32JCANON_PRIME|prime={p}|reason={reason}|status=SKIP",
            flush=True,
        )
        continue
    print(
        "Q32JCANON_PRIME|"
        f"prime={p}|Pdeg={res['P'].degree()}|Qdeg={res['Q'].degree()}|"
        f"q9={int(res['Q'][9])}|"
        f"q8eq7={int(res['Q'][8]==res['Q'][7])}|"
        f"cube={int(res['cube_ok'])}|square={int(res['square_ok'])}|"
        "status=PASS",
        flush=True,
    )
    good.append((p,res))

if len(good)<3:
    raise SystemExit("too few canonical j maps")

print(
    "Q32JCANON_GROUP|"
    f"compatible={len(good)}/{len(records)}|"
    f"primes={','.join(str(p) for p,_ in good)}|"
    "profile=P18/Q10,q9=0,q8=q7|status=PASS",
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

# Last good prime held out.
train=good[:-1]
holdp,hold=good[-1]
mods=[p for p,_ in train]

results={}
complete=True
for name,n in (("P",19),("Q",11)):
    arrays=[]
    for unused_p,r in train:
        poly=r[name]
        arrays.append([ZZ(poly[i]) for i in range(n)])
    hv=[ZZ(hold[name][i]) for i in range(n)]

    qs=[]
    recovered=held=0
    maxnb=maxdb=0
    stable=0
    stableheld=0

    for j in range(n):
        q,M=rr_scalar([a[j] for a in arrays],mods)
        qs.append(q)
        if q is not None:
            recovered+=1
            maxnb=max(maxnb,abs(ZZ(q.numerator())).nbits())
            maxdb=max(maxdb,abs(ZZ(q.denominator())).nbits())
            if red(q,holdp)==int(hv[j]%holdp):
                held+=1

        # last-two-prefix stability
        if len(train)>=4:
            q1,_=rr_scalar([a[j] for a in arrays[:-1]],mods[:-1])
            q2=q
            if q1 is not None and q2 is not None and q1==q2:
                stable+=1
                if red(q2,holdp)==int(hv[j]%holdp):
                    stableheld+=1

    objok=(held==n)
    complete &= objok
    results[name]=qs
    print(
        "Q32JCANON_CRT|"
        f"object={name}|train={len(train)}|holdout={holdp}|"
        f"recovered={recovered}/{n}|heldout={held}/{n}|"
        f"stable={stable}|stable_heldout={stableheld}|"
        f"max_num_bits={maxnb}|max_den_bits={maxdb}|"
        f"status={'PASS_HELDOUT' if objok else 'PARTIAL'}",
        flush=True,
    )

# Exact characteristic-zero sanity checks if fully reconstructed.
exact_cube=False
exact_square=False
if complete:
    RQ=PolynomialRing(QQ,"u")
    P=RQ(results["P"])
    Q=RQ(results["Q"])
    assert P.degree()==18 and Q.degree()==10
    assert Q.is_monic() and Q[9]==0 and Q[8]==Q[7]

    exact_cube=all(int(e)%3==0 for unused_f,e in P.factor())
    H=P-QQ(1728)*Q
    exact_square=all(int(e)%2==0 for unused_f,e in H.factor())

    print(
        "Q32JCANON_EXACT|"
        f"cube={int(exact_cube)}|square={int(exact_square)}|"
        f"status={'PASS_J_BELYI_SHAPE' if exact_cube and exact_square else 'FAIL_SHAPE'}",
        flush=True,
    )

out=LOCAL/"q32-d12-canonical-j-crt.json"
payload={
    "schema":"elkies-k3.h3-q32-d12-canonical-j-crt.v1",
    "status":(
        "PASS_EXACT_Q32_D12_CANONICAL_J"
        if complete and exact_cube and exact_square else
        "PARTIAL_Q32_D12_CANONICAL_J_CRT"
    ),
    "primes":[int(p) for p,_ in good],
    "training_primes":[int(p) for p,_ in train],
    "heldout_prime":int(holdp),
    "normalization":{
        "I8star":"infinity",
        "denominator":"monic_degree_10",
        "q9":"0",
        "scale_condition":"q8=q7",
    },
}
if complete:
    payload["j_numerator_coefficients_low_to_high"]=[str(q) for q in results["P"]]
    payload["j_denominator_coefficients_low_to_high"]=[str(q) for q in results["Q"]]
    payload["cube_shape"]=bool(exact_cube)
    payload["square_shape"]=bool(exact_square)

out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
print(
    "Q32JCANON_RESULT|"
    f"good={len(good)}|train={len(train)}|holdout={holdp}|"
    f"status={payload['status']}",
    flush=True,
)
