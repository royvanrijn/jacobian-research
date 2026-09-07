#!/usr/bin/env sage -python
"""
Normalize the modular q32 D12 polynomial model geometrically by moving its
unique I8* fibre to infinity, then reconstruct that normalized model over QQ.

For each good prime:
  A=A0/h^4, B=B0/h^6, xP=X0/h^2.
Let beta be the root of the unique multiplicity-14 discriminant factor.
Set T=1/(V-beta) and transform line-bundle sections by weights 8,12,4:
  Ainf(T)=T^8 A0(beta+1/T)
  Binf(T)=T^12 B0(beta+1/T)
  Xinf(T)=T^4 X0(beta+1/T).

For I8* at infinity, expected degrees are Ainf<=6, Binf<=9.
This bypasses reconstruction of h and of the original V-coordinate.
"""

import json, math
from itertools import combinations
from pathlib import Path
from sage.all import GF, PolynomialRing, QQ, ZZ, binomial, matrix

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
if len(records)<5:
    raise SystemExit("need at least five good primes")


def modular_model(row):
    p,sig,pt=row
    F=GF(p)
    R=PolynomialRing(F,"V")
    V=R.gen()

    def parts(rec):
        n=R([F(v) for v in rec["num"]])
        d=R([F(v) for v in rec["den"]])
        lc=d.leading_coefficient()
        return n/lc,d/lc

    A0,Ad=parts(sig["jacobian_A"])
    B0,Bd=parts(sig["jacobian_B"])
    X0,Xd=parts(pt["marked_section_x"])

    def exact_root(P,e):
        out=R.one()
        for f,m in R(P).factor():
            m=int(m)
            if m%e:
                return None
            out*=f.monic()**(m//e)
        return out.monic()

    h=exact_root(Ad,4)
    assert h is not None and h==exact_root(Bd,6)==exact_root(Xd,2)
    assert Ad==h**4 and Bd==h**6 and Xd==h**2

    Delta=-16*(4*A0**3+27*B0**2)
    fac=Delta.factor()
    stars=[(f,e) for f,e in fac if int(e)==14 and f.degree()==1]
    if len(stars)!=1:
        raise ArithmeticError(f"prime {p}: expected unique linear I8* factor^14, got {fac}")
    fstar=stars[0][0].monic()
    beta=-fstar[0]

    T=R.gen()  # just reuse polynomial symbol algebraically

    def infinity_transform(P,w):
        P=R(P)
        out=R.zero()
        for i,c in enumerate(P.list()):
            if not c:
                continue
            for j in range(i+1):
                exponent=w-j
                assert exponent>=0
                out += c*F(binomial(i,j))*beta**(i-j)*T**exponent
        return R(out)

    Ai=infinity_transform(A0,8)
    Bi=infinity_transform(B0,12)
    Xi=infinity_transform(X0,4)

    # Cancellation at beta encodes the I8* orders.
    if Ai.degree()>6 or Bi.degree()>9:
        raise ArithmeticError(
            f"prime {p}: I8* infinity degree failure A={Ai.degree()} B={Bi.degree()}"
        )

    Di=-16*(4*Ai**3+27*Bi**2)
    # Finite discriminant should consist of ten I1 fibres generically.
    sqfree=(Di.gcd(Di.derivative()).degree()==0)

    print(
        "Q32I8INF_MOD|"
        f"prime={p}|beta={int(beta)}|"
        f"Adeg={Ai.degree()}|Bdeg={Bi.degree()}|Xdeg={Xi.degree()}|"
        f"Ddeg={Di.degree()}|Dsqfree={int(sqfree)}|"
        f"status={'PASS' if Ai.degree()==6 and Bi.degree()==9 else 'PASS_BOUND'}",
        flush=True,
    )

    return {
        "beta":ZZ(int(beta)),
        "A":[ZZ(int(Ai[i])) for i in range(7)],
        "B":[ZZ(int(Bi[i])) for i in range(10)],
        "X":[ZZ(int(Xi[i])) for i in range(5)],
    }


mods_data=[(r[0],modular_model(r)) for r in records]
train=mods_data[:-1]
holdp,hold=mods_data[-1]
primes=[p for p,_ in train]

print(
    "Q32I8INF_INPUT|"
    f"train={','.join(map(str,primes))}|holdout={holdp}|count={len(records)}|status=PASS",
    flush=True,
)


def crt_scalar(vals,mods):
    x=ZZ(0); M=ZZ(1)
    for a,p in zip(vals,mods):
        a=ZZ(a)%p
        t=((a-x)%p)*((M%p).inverse_mod(p))%p
        x=(x+M*t)%(M*p)
        M*=p
    if x>M//2:
        x-=M
    return x,M

_,MOD=crt_scalar([0]*len(primes),primes)
print(f"Q32I8INF_MODULUS|bits={MOD.nbits()}|status=PASS",flush=True)


def primitive(v):
    v=[ZZ(x) for x in v]
    g=ZZ(0)
    for x in v:
        g=ZZ(math.gcd(int(g),abs(int(x))))
    if g>1:
        v=[x//g for x in v]
    if v[-1]<0:
        v=[-x for x in v]
    return v


def candidates(Rred):
    rows=[list(v) for v in Rred.rows()]
    seen=set()
    def emit(v):
        v=primitive(v); k=tuple(v)
        if v[-1] and k not in seen:
            seen.add(k); return v
        return None
    for v in rows:
        q=emit(v)
        if q is not None: yield q
    k=min(12,len(rows))
    for i,j in combinations(range(k),2):
        for s in (1,-1):
            q=emit([rows[i][c]+s*rows[j][c] for c in range(len(rows[i]))])
            if q is not None: yield q
    if len(rows)<=12:
        k=min(7,len(rows))
        for i,j,l in combinations(range(k),3):
            for s1 in (1,-1):
                for s2 in (1,-1):
                    q=emit([
                        rows[i][c]+s1*rows[j][c]+s2*rows[l][c]
                        for c in range(len(rows[i]))
                    ])
                    if q is not None: yield q


def reconstruct(name):
    arrays=[d[name] for _,d in train]
    hv=hold[name]
    N=len(arrays[0])
    residues=[]
    for j in range(N):
        x,M=crt_scalar([a[j] for a in arrays],primes)
        assert M==MOD
        residues.append(x)

    B=matrix(ZZ,N+1,N+1)
    for j in range(N):
        B[j,j]=MOD
    for j,r in enumerate(residues):
        B[N,j]=r
    B[N,N]=1

    print(
        f"Q32I8INF_LLL_START|object={name}|dimension={N+1}|"
        f"modulus_bits={MOD.nbits()}|status=START",
        flush=True,
    )
    RR=B.LLL(delta=0.99)

    scored=[]
    for v in candidates(RR):
        d=v[-1]
        if d%holdp==0:
            m=-1
        else:
            inv=(d%holdp).inverse_mod(holdp)
            m=sum(
                int((v[j]%holdp)*inv%holdp)==int(hv[j]%holdp)
                for j in range(N)
            )
        bits=max(abs(x).nbits() for x in v)
        scored.append((m,bits,sum(x*x for x in v),v))
    scored.sort(key=lambda z:(-z[0],z[1],z[2]))
    m,bits,_,v=scored[0]

    print(
        f"Q32I8INF_LLL_BEST|object={name}|heldout={m}/{N}|"
        f"height_bits={bits}|common_den_bits={abs(v[-1]).nbits()}|"
        f"status={'PASS_HELDOUT' if m==N else 'PARTIAL'}",
        flush=True,
    )
    for rank,(mm,bb,nn,vv) in enumerate(scored[:5]):
        print(
            f"Q32I8INF_LLL_SHORT|object={name}|rank={rank}|heldout={mm}/{N}|"
            f"height_bits={bb}|den_bits={abs(vv[-1]).nbits()}",
            flush=True,
        )

    if m!=N:
        return None
    d=QQ(v[-1])
    return [QQ(x)/d for x in v[:-1]],v


exact={}; vectors={}
# X first: cheapest diagnostic.
for name in ("X","A","B"):
    ans=reconstruct(name)
    if ans is None:
        print(
            f"Q32I8INF_RESULT|failed={name}|modulus_bits={MOD.nbits()}|"
            "status=NEED_MORE_PRIMES_OR_DIFFERENT_NORMALIZATION",
            flush=True,
        )
        raise SystemExit(0)
    exact[name],vectors[name]=ans

RQ=PolynomialRing(QQ,"T")
A=RQ(exact["A"])
B=RQ(exact["B"])
X=RQ(exact["X"])
assert A.degree()<=6 and B.degree()<=9 and X.degree()<=4

rhs=X**3+A*X+B

def qqsqrt(q):
    q=QQ(q)
    if q<0:return None
    n=ZZ(q.numerator());d=ZZ(q.denominator())
    if not n.is_square() or not d.is_square():return None
    return QQ(n.sqrt())/QQ(d.sqrt())

def psqrt(P):
    P=RQ(P)
    if not P:return RQ.zero()
    fac=P.factor()
    u=qqsqrt(fac.unit())
    if u is None:return None
    out=RQ(u)
    for f,e in fac:
        if int(e)%2:return None
        out*=f**(int(e)//2)
    assert out**2==P
    return out

Y=psqrt(rhs)
if Y is None:
    print("Q32I8INF_SECTION|square=0|status=MODEL_RECOVERED_SECTION_FAILED",flush=True)
    raise SystemExit(0)

Delta=-16*(4*A**3+27*B**2)
print(
    "Q32I8INF_EXACT|"
    f"Adeg={A.degree()}|Bdeg={B.degree()}|Xdeg={X.degree()}|Ydeg={Y.degree()}|"
    f"Ddeg={Delta.degree()}|identity=1|status=PASS",
    flush=True,
)

def serial(P):
    return {
        "degree":int(P.degree()),
        "coefficients_low_to_high":[str(c) for c in P.list()],
    }

payload={
    "schema":"elkies-k3.h3-q32-d12-i8star-infinity-qq.v1",
    "status":"PASS_EXACT_Q32_D12_I8STAR_INFINITY_QQ_HELDOUT",
    "training_primes":[int(p) for p in primes],
    "heldout_prime":int(holdp),
    "crt_modulus_bits":int(MOD.nbits()),
    "base_normalization":"T=1/(V-beta), beta=unique I8* root",
    "A":serial(A),"B":serial(B),
    "spinor_X":serial(X),"spinor_Y":serial(Y),
    "discriminant":serial(Delta),
    "exact_weierstrass_identity":True,
}
out=LOCAL/"q32-d12-i8star-infinity-qq.json"
out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
print(
    "Q32I8INF_RESULT|"
    f"modulus_bits={MOD.nbits()}|holdout={holdp}|"
    f"status={payload['status']}",
    flush=True,
)
