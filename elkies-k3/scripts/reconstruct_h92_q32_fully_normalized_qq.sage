#!/usr/bin/env sage -python
"""
Fully normalize the q32 D12 model across characteristics.

Stage 1: move unique I8* fibre to infinity, giving degrees A=6,B=9,X=4 and
degree-10 squarefree finite discriminant.

Residual equivalences are
    T -> a*T + b
    x -> u^2*x, y -> u^3*y.

Fix translation by centering the ten finite discriminant roots:
    mu = -d9/(10*d10),  z=T-mu.

For leading coefficients a6,b9,x4, the residual scalings are
    a6 -> u^4 a^-6 a6
    b9 -> u^6 a^-9 b9
    x4 -> u^2 a^-4 x4.

Hence
    q = b9/(a6*x4) -> a*q
    r = a6^2/x4^3 -> u^2*r.

Set canonical base S=z/q, i.e. z=q*S, then divide
    X by r, A by r^2, B by r^3.

The result is invariant under the residual affine base and Weierstrass scaling,
so coefficients from all good primes should be reductions of one QQ model.
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


def normalize(row):
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

    def exact_root(P,e):
        out=R.one()
        for f,m in R(P).factor():
            m=int(m)
            if m%e:return None
            out*=f.monic()**(m//e)
        return out.monic()

    h=exact_root(Ad,4)
    assert h is not None and h==exact_root(Bd,6)==exact_root(Xd,2)
    assert Ad==h**4 and Bd==h**6 and Xd==h**2

    Delta=-16*(4*A0**3+27*B0**2)
    stars=[(f,e) for f,e in Delta.factor() if int(e)==14 and f.degree()==1]
    if len(stars)!=1:
        raise ArithmeticError(f"prime {p}: no unique I8*")
    beta=-stars[0][0].monic()[0]

    def i8_transform(P,w):
        P=R(P)
        out=R.zero()
        for i,c in enumerate(P.list()):
            for j in range(i+1):
                out += c*F(binomial(i,j))*beta**(i-j)*S**(w-j)
        return R(out)

    A=i8_transform(A0,8)
    B=i8_transform(B0,12)
    X=i8_transform(X0,4)
    assert A.degree()==6 and B.degree()==9 and X.degree()==4

    D=-16*(4*A**3+27*B**2)
    assert D.degree()==10 and D.gcd(D.derivative()).degree()==0

    # Center finite discriminant roots.
    d10=D[10]
    d9=D[9]
    mu=-d9/(F(10)*d10)

    def translate(P,shift):
        P=R(P)
        out=R.zero()
        for i,c in enumerate(P.list()):
            for j in range(i+1):
                out += c*F(binomial(i,j))*shift**(i-j)*S**j
        return R(out)

    Ac=translate(A,mu)
    Bc=translate(B,mu)
    Xc=translate(X,mu)
    Dc=-16*(4*Ac**3+27*Bc**2)
    assert Dc.degree()==10 and Dc[9]==0

    a6=Ac[6]
    b9=Bc[9]
    x4=Xc[4]
    if not a6 or not b9 or not x4:
        raise ArithmeticError(f"prime {p}: leading coefficient degeneration")

    q=b9/(a6*x4)
    r=a6**2/x4**3
    assert q and r

    # Canonical base z=q*S.
    def scale_base(P,q):
        P=R(P)
        return R([c*q**i for i,c in enumerate(P.list())])

    An=scale_base(Ac,q)/(r**2)
    Bn=scale_base(Bc,q)/(r**3)
    Xn=scale_base(Xc,q)/r

    Dn=-16*(4*An**3+27*Bn**2)
    assert Dn.degree()==10 and Dn[9]==0
    assert Dn.gcd(Dn.derivative()).degree()==0

    # These become characteristic-independent normalizations.
    # q,r were chosen so B9/(A6*X4)=1 and A6^2/X4^3=1.
    aa=An[6]; bb=Bn[9]; xx=Xn[4]
    assert bb==aa*xx
    assert aa**2==xx**3

    print(
        "Q32FULLNORM_MOD|"
        f"prime={p}|beta={int(beta)}|mu={int(mu)}|"
        f"Adeg={An.degree()}|Bdeg={Bn.degree()}|Xdeg={Xn.degree()}|Ddeg={Dn.degree()}|"
        f"D9zero={int(Dn[9]==0)}|lead_rel=1|status=PASS",
        flush=True,
    )

    return {
        "A":[ZZ(int(An[i])) for i in range(7)],
        "B":[ZZ(int(Bn[i])) for i in range(10)],
        "X":[ZZ(int(Xn[i])) for i in range(5)],
    }


data=[(r[0],normalize(r)) for r in records]
train=data[:-1]
holdp,hold=data[-1]
primes=[p for p,_ in train]

print(
    "Q32FULLNORM_INPUT|"
    f"train={','.join(map(str,primes))}|holdout={holdp}|count={len(data)}|status=PASS",
    flush=True,
)


def crt(vals,mods):
    x=ZZ(0);M=ZZ(1)
    for a,p in zip(vals,mods):
        a=ZZ(a)%p
        t=((a-x)%p)*((M%p).inverse_mod(p))%p
        x=(x+M*t)%(M*p);M*=p
    if x>M//2:x-=M
    return x,M

_,MOD=crt([0]*len(primes),primes)
print(f"Q32FULLNORM_MODULUS|bits={MOD.nbits()}|status=PASS",flush=True)


def primitive(v):
    v=[ZZ(x) for x in v]
    g=ZZ(0)
    for x in v:g=ZZ(math.gcd(int(g),abs(int(x))))
    if g>1:v=[x//g for x in v]
    if v[-1]<0:v=[-x for x in v]
    return v


def candidate_vectors(RR):
    rows=[list(v) for v in RR.rows()]
    seen=set()
    def emit(v):
        v=primitive(v);k=tuple(v)
        if v[-1] and k not in seen:
            seen.add(k);return v
        return None
    for v in rows:
        q=emit(v)
        if q is not None:yield q
    k=min(12,len(rows))
    for i,j in combinations(range(k),2):
        for s in (1,-1):
            q=emit([rows[i][c]+s*rows[j][c] for c in range(len(rows[i]))])
            if q is not None:yield q
    if len(rows)<=12:
        k=min(7,len(rows))
        for i,j,l in combinations(range(k),3):
            for s1 in (1,-1):
                for s2 in (1,-1):
                    q=emit([
                        rows[i][c]+s1*rows[j][c]+s2*rows[l][c]
                        for c in range(len(rows[i]))
                    ])
                    if q is not None:yield q


def scalar_rr_object(name):
    arrs=[d[name] for _,d in train]
    hv=hold[name]
    vals=[]
    ok=0
    maxnb=maxdb=0
    for j in range(len(arrs[0])):
        x,M=crt([a[j] for a in arrs],primes)
        try:
            q=QQ(Integers(M)(x).rational_reconstruction())
        except Exception:
            vals.append(None);continue
        vals.append(q)
        maxnb=max(maxnb,abs(ZZ(q.numerator())).nbits())
        maxdb=max(maxdb,abs(ZZ(q.denominator())).nbits())
        den=ZZ(q.denominator())%holdp
        if den and int((ZZ(q.numerator())%holdp)*den.inverse_mod(holdp)%holdp)==int(hv[j]%holdp):
            ok+=1
    return vals,ok,maxnb,maxdb


# Need Integers only here; import lazily.
from sage.all import Integers

def reconstruct(name):
    arrs=[d[name] for _,d in train]
    hv=hold[name]
    N=len(arrs[0])

    # First try ordinary RR now that normalization should make coefficients small.
    vals,ok,nb,db=scalar_rr_object(name)
    print(
        f"Q32FULLNORM_RR|object={name}|heldout={ok}/{N}|"
        f"max_num_bits={nb}|max_den_bits={db}|"
        f"status={'PASS_HELDOUT' if ok==N else 'PARTIAL'}",
        flush=True,
    )
    if ok==N and all(v is not None for v in vals):
        return vals,None

    residues=[]
    for j in range(N):
        x,M=crt([a[j] for a in arrs],primes)
        assert M==MOD
        residues.append(x)

    Bmat=matrix(ZZ,N+1,N+1)
    for j in range(N):Bmat[j,j]=MOD
    for j,z in enumerate(residues):Bmat[N,j]=z
    Bmat[N,N]=1

    RR=Bmat.LLL(delta=0.99)
    scored=[]
    for v in candidate_vectors(RR):
        d=v[-1]
        if d%holdp==0:m=-1
        else:
            inv=(d%holdp).inverse_mod(holdp)
            m=sum(int((v[j]%holdp)*inv%holdp)==int(hv[j]%holdp) for j in range(N))
        bits=max(abs(x).nbits() for x in v)
        scored.append((m,bits,sum(x*x for x in v),v))
    scored.sort(key=lambda z:(-z[0],z[1],z[2]))
    m,bits,_,v=scored[0]
    print(
        f"Q32FULLNORM_LLL|object={name}|heldout={m}/{N}|"
        f"height_bits={bits}|den_bits={abs(v[-1]).nbits()}|"
        f"status={'PASS_HELDOUT' if m==N else 'PARTIAL'}",
        flush=True,
    )
    if m!=N:return None,None
    d=QQ(v[-1])
    return [QQ(x)/d for x in v[:-1]],v


exact={};vecs={}
for name in ("X","A","B"):
    vals,v=reconstruct(name)
    if vals is None:
        print(
            f"Q32FULLNORM_RESULT|failed={name}|modulus_bits={MOD.nbits()}|"
            "status=NEED_MORE_PRIMES",
            flush=True,
        )
        raise SystemExit(0)
    exact[name]=vals
    if v is not None:vecs[name]=[str(x) for x in v]

RQ=PolynomialRing(QQ,"S")
X=RQ(exact["X"]);A=RQ(exact["A"]);B=RQ(exact["B"])
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
    fac=P.factor();u=qqsqrt(fac.unit())
    if u is None:return None
    out=RQ(u)
    for f,e in fac:
        if int(e)%2:return None
        out*=f**(int(e)//2)
    assert out**2==P
    return out

Y=psqrt(rhs)
if Y is None:
    print("Q32FULLNORM_SECTION|square=0|status=MODEL_PASS_SECTION_FAIL",flush=True)
    raise SystemExit(0)

D=-16*(4*A**3+27*B**2)
assert D.degree()==10 and D[9]==0
print(
    "Q32FULLNORM_EXACT|"
    f"Adeg={A.degree()}|Bdeg={B.degree()}|Xdeg={X.degree()}|Ydeg={Y.degree()}|"
    f"Ddeg={D.degree()}|D9zero=1|identity=1|status=PASS",
    flush=True,
)

def serial(P):
    return {"degree":int(P.degree()),"coefficients_low_to_high":[str(c) for c in P.list()]}

payload={
    "schema":"elkies-k3.h3-q32-d12-fully-normalized-qq.v1",
    "status":"PASS_EXACT_Q32_D12_FULLY_NORMALIZED_QQ_HELDOUT",
    "training_primes":[int(p) for p in primes],
    "heldout_prime":int(holdp),
    "crt_modulus_bits":int(MOD.nbits()),
    "normalization":{
        "I8star":"infinity",
        "finite_discriminant_centroid":"0",
        "base_scale":"q=b9/(a6*x4)",
        "weierstrass_x_scale":"r=a6^2/x4^3",
    },
    "A":serial(A),"B":serial(B),"spinor_X":serial(X),"spinor_Y":serial(Y),
    "discriminant":serial(D),"exact_weierstrass_identity":True,
    "lll_vectors":vecs,
}
out=LOCAL/"q32-d12-fully-normalized-qq.json"
out.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
print(
    "Q32FULLNORM_RESULT|"
    f"modulus_bits={MOD.nbits()}|holdout={holdp}|status={payload['status']}",
    flush=True,
)
