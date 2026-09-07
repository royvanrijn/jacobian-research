#!/usr/bin/env sage -python
"""
Enumerate polynomial sections on the canonical q32/D12 child over GF(p).

Target model:
    E: y^2 = x^3 - 3*a6(u)*x + 2*b9(u)

Zero-pole sections satisfy deg x<=4, deg y<=6.  The D12/MW5 lattice predicts
exactly 26 nonzero such sections spanning rank 4.

We solve the two natural degree profiles:
  A) deg x=3.  Cancellation of u^9 forces leading x coefficient c to satisfy
         c^3 - 3c + 2 = (c-1)^2(c+2)=0,
     so c in {1,-2}; then deg y<=4.
  B) deg x=4.  Write lc(x)=s^2, lc(y)=s^3 and recursively solve the top
     coefficients of y from the square-root-at-infinity expansion.

The remaining low coefficients form small zero-dimensional polynomial ideals.
"""

import argparse, json
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ, EllipticCurve, vector

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

parser=argparse.ArgumentParser()
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()
p=ZZ(args.prime)

SIG=LOCAL/f"q32-signature-mod-{p}.json"
if not SIG.exists():
    raise SystemExit(f"missing {SIG}")
d=json.loads(SIG.read_text())
assert d["status"]=="PASS_Q32_MODP_SIGNATURE"

F=GF(p)

# -------------------------------------------------------------------------
# Canonical child a6,b9, exactly as in the passing canonical-j/a,b probes.
# -------------------------------------------------------------------------
RV=PolynomialRing(F,"V"); V=RV.gen(); KV=RV.fraction_field()

def rf(rec):
    return KV(RV([F(v) for v in rec["num"]]))/KV(
        RV([F(v) for v in rec["den"]])
    )

A0=rf(d["jacobian_A"])
B0=rf(d["jacobian_B"])
J=KV(F(6912))*A0**3/(KV(F(4))*A0**3+KV(F(27))*B0**2)
N=RV(J.numerator()); D=RV(J.denominator())
assert N.gcd(D).degree()==0 and N.degree()==D.degree()==18

e8=[f for f,e in D.factor() if int(e)==8 and f.degree()==1]
assert len(e8)==1
f=e8[0]; r=-f[0]/f[1]

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
assert P.degree()==18 and Q.degree()==10
lc=Q.leading_coefficient(); P/=lc; Q/=lc
c=-Q[9]/F(10)

RS=PolynomialRing(F,"S"); S0=RS.gen()
P1=RS(P(S0+c)); Q1=RS(Q(S0+c))
assert Q1.is_monic() and Q1[9]==0
scale=Q1[7]/Q1[8]
RU=PolynomialRing(F,"u"); u=RU.gen()
P2=RU(P1(scale*u)); Q2=RU(Q1(scale*u))
lc2=Q2.leading_coefficient(); P2/=lc2; Q2/=lc2
assert Q2.is_monic() and Q2[9]==0 and Q2[8]==Q2[7]

def monic_power_root(poly,e):
    out=RU.one()
    for fac,m in poly.factor():
        assert int(m)%e==0
        out*=fac.monic()**(int(m)//e)
    return out.monic()

a=monic_power_root(P2,3)
H2=P2-F(1728)*Q2
b=monic_power_root(H2,2)
assert a.degree()==6 and b.degree()==9
assert (a**3-b**2).degree()==10

print(
    "Q32POLYSEC_MODEL|"
    f"prime={p}|adeg={a.degree()}|bdeg={b.degree()}|"
    f"Kdeg={(a**3-b**2).degree()}|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# Recursive square-root solver.
# -------------------------------------------------------------------------
def solve_branch_dx3(clead):
    names=("ell","x2","x1","x0","inv")
    S=PolynomialRing(F,names=names,order="degrevlex")
    ell,x2,x1,x0,inv=S.gens()
    K=S.fraction_field()
    U=PolynomialRing(K,"z"); z=U.gen()
    aa=U([K(v) for v in a.list()])
    bb=U([K(v) for v in b.list()])
    x=K(F(clead))*z**3+K(x2)*z**2+K(x1)*z+K(x0)
    R=x**3-K(3)*aa*x+K(2)*bb
    assert R.degree()<=8

    ys={4:K(ell)}
    equations=[S((K(ell)**2-K(R[8])).numerator())]
    for k in range(7,3,-1):
        j=k-4
        known=K(0)
        for i in range(5):
            q=k-i
            if q<0 or q>4 or i==4 or q==4:
                continue
            if i in ys and q in ys:
                known+=ys[i]*ys[q]
        # known loop above excludes terms involving y4; equivalently:
        known=sum(
            ys[i]*ys[k-i]
            for i in ys
            if (k-i) in ys and i!=4 and (k-i)!=4
        )
        yj=(K(R[k])-known)/(K(2)*ys[4])
        ys[j]=yj

    y=sum(ys[i]*z**i for i in range(5))
    residual=y**2-R
    for k in range(4):
        equations.append(S(K(residual[k]).numerator()))
    equations.append(inv*ell-1)

    I=S.ideal(equations)
    print(
        "Q32POLYSEC_BRANCH|"
        f"profile=dx3|lead={int(F(clead))}|vars={S.ngens()}|"
        f"equations={len(equations)}|status=GROEBNER_START",
        flush=True,
    )
    sols=I.variety()
    print(
        "Q32POLYSEC_BRANCH|"
        f"profile=dx3|lead={int(F(clead))}|solutions={len(sols)}|"
        "status=GROEBNER_PASS",
        flush=True,
    )

    out=[]
    for sol in sols:
        vals={g:F(sol[g]) for g in S.gens()}
        if not vals[ell]:
            continue
        xx=RU([
            vals[x0],vals[x1],vals[x2],F(clead)
        ])
        # Re-evaluate recursive y coefficients numerically from RHS.
        RR=xx**3-F(3)*a*xx+F(2)*b
        yy=[F(0)]*5
        yy[4]=vals[ell]
        for k in range(7,3,-1):
            j=k-4
            known=sum(
                yy[i]*yy[k-i]
                for i in range(5)
                if 0<=k-i<5 and i!=4 and (k-i)!=4
            )
            yy[j]=(RR[k]-known)/(F(2)*yy[4])
        yy=RU(yy)
        assert yy**2==RR
        out.append((xx,yy,"dx3"))
    return out

def solve_branch_dx4():
    names=("s","x3","x2","x1","x0","inv")
    S=PolynomialRing(F,names=names,order="degrevlex")
    ss,x3,x2,x1,x0,inv=S.gens()
    K=S.fraction_field()
    U=PolynomialRing(K,"z"); z=U.gen()
    aa=U([K(v) for v in a.list()])
    bb=U([K(v) for v in b.list()])
    x=K(ss)**2*z**4+K(x3)*z**3+K(x2)*z**2+K(x1)*z+K(x0)
    R=x**3-K(3)*aa*x+K(2)*bb
    assert R.degree()<=12

    ys={6:K(ss)**3}
    for k in range(11,5,-1):
        j=k-6
        known=sum(
            ys[i]*ys[k-i]
            for i in ys
            if (k-i) in ys and i!=6 and (k-i)!=6
        )
        ys[j]=(K(R[k])-known)/(K(2)*ys[6])

    y=sum(ys[i]*z**i for i in range(7))
    residual=y**2-R
    equations=[S(K(residual[k]).numerator()) for k in range(6)]
    equations.append(inv*ss-1)
    I=S.ideal(equations)

    print(
        "Q32POLYSEC_BRANCH|profile=dx4|vars=6|equations=7|"
        "status=GROEBNER_START",
        flush=True,
    )
    sols=I.variety()
    print(
        f"Q32POLYSEC_BRANCH|profile=dx4|solutions={len(sols)}|"
        "status=GROEBNER_PASS",
        flush=True,
    )

    out=[]
    for sol in sols:
        vals={g:F(sol[g]) for g in S.gens()}
        if not vals[ss]:
            continue
        xx=RU([
            vals[x0],vals[x1],vals[x2],vals[x3],vals[ss]**2
        ])
        RR=xx**3-F(3)*a*xx+F(2)*b
        yy=[F(0)]*7
        yy[6]=vals[ss]**3
        for k in range(11,5,-1):
            j=k-6
            known=sum(
                yy[i]*yy[k-i]
                for i in range(7)
                if 0<=k-i<7 and i!=6 and (k-i)!=6
            )
            yy[j]=(RR[k]-known)/(F(2)*yy[6])
        yy=RU(yy)
        assert yy**2==RR
        out.append((xx,yy,"dx4"))
    return out

sections=[]
sections += solve_branch_dx3(F(1))
sections += solve_branch_dx3(F(-2))
sections += solve_branch_dx4()

# Deduplicate exact polynomial points.
unique={}
for x,y,profile in sections:
    key=(
        tuple(int(x[i]) for i in range(5)),
        tuple(int(y[i]) for i in range(7)),
    )
    unique[key]=(x,y,profile)
sections=list(unique.values())

print(
    "Q32POLYSEC_COUNT|"
    f"prime={p}|nonzero={len(sections)}|expected=26|"
    f"profiles="+",".join(
        f"{name}:{sum(1 for _,_,q in sections if q==name)}"
        for name in ("dx3","dx4")
    )+
    f"|status={'PASS_26' if len(sections)==26 else 'PARTIAL'}",
    flush=True,
)

# -------------------------------------------------------------------------
# Identify the unique height-1 pair intrinsically: P and -P are the only
# zero-pole sections whose doubles are again among the 26 zero-pole sections.
# -------------------------------------------------------------------------
Kfun=RU.fraction_field()
E=EllipticCurve(Kfun,[0,0,0,Kfun(-3*a),Kfun(2*b)])

def coeff_key_from_point(P):
    if P.is_zero():
        return None
    xx,yy=P.xy()
    if RU(xx.denominator()).degree()>0 or RU(yy.denominator()).degree()>0:
        return None
    xp=RU(xx.numerator())/F(xx.denominator())
    yp=RU(yy.numerator())/F(yy.denominator())
    if xp.degree()>4 or yp.degree()>6:
        return None
    return (
        tuple(int(xp[i]) for i in range(5)),
        tuple(int(yp[i]) for i in range(7)),
    )

keys=set(unique)
double_closed=[]
for x,y,profile in sections:
    Pnt=E(Kfun(x),Kfun(y))
    k2=coeff_key_from_point(2*Pnt)
    if k2 in keys:
        double_closed.append((x,y,profile,k2))

print(
    "Q32POLYSEC_DOUBLE|"
    f"count={len(double_closed)}|"
    f"status={'PASS_HEIGHT1_PAIR' if len(double_closed)==2 else 'PARTIAL'}",
    flush=True,
)

for i,(x,y,profile,k2) in enumerate(double_closed):
    print(
        "Q32POLYSEC_HEIGHT1|"
        f"index={i}|profile={profile}|"
        f"xdeg={x.degree()}|ydeg={y.degree()}|"
        f"x={','.join(str(int(x[j])) for j in range(5))}|"
        f"y={','.join(str(int(y[j])) for j in range(7))}|"
        "status=PASS",
        flush=True,
    )

# -------------------------------------------------------------------------
# Point-adapted normalization from the unique +/- height-1 pair.
# Both signs produce the same normalized model.
# -------------------------------------------------------------------------
adapted=None
if len(double_closed)==2:
    xP,yP,_profile,_=double_closed[0]

    # Leading degree may be dx3 or dx4.  First use the section's own leading
    # coefficients to fix Weierstrass scaling whenever their weights allow it.
    dx=int(xP.degree()); dy=int(yP.degree())
    lx=xP[dx]; ly=yP[dy]

    # For the expected spinor height-1 profile dx=3,dy=4 the local infinity
    # normalization is already rigid in the canonical short model; use xP
    # itself to fix the affine base coordinate.
    W=PolynomialRing(F,"w"); w=W.gen()

    # Kill xP's next-to-leading coefficient by an affine translation.
    # For monic leading coefficient c, x=c*u^d+q*u^(d-1)+..., use
    # u=v-q/(d*c).
    shift=-lx**(-1)*xP[dx-1]/F(dx)
    V2=PolynomialRing(F,"v"); v=V2.gen()
    xp1=V2(xP(v+shift))
    yp1=V2(yP(v+shift))
    a1=V2(a(v+shift)); b1=V2(b(v+shift))
    assert xp1[dx-1]==0

    # Rational scale from the next two nonzero coefficients of xP:
    # under v=lambda*w and Weierstrass line-bundle scaling, their ratio fixes
    # lambda without extracting roots.  Choose lambda so normalized
    # coefficients of w^(d-2), w^(d-3) are equal when both are nonzero.
    q2=xp1[dx-2] if dx>=2 else F(0)
    q3=xp1[dx-3] if dx>=3 else F(0)
    if q2 and q3:
        lam=q3/q2
    else:
        lam=F(1)

    # Substitute base.  Do not yet apply an additional x/y line-bundle
    # scaling; recording raw substituted coefficients is enough to test
    # cross-prime arithmetic simplification in the next probe.
    xp2=W(xp1(lam*w))
    yp2=W(yp1(lam*w))
    aa2=W(a1(lam*w))
    bb2=W(b1(lam*w))

    adapted={
        "profile":[dx,dy],
        "base_shift":int(shift),
        "base_scale":int(lam),
        "xP":[int(xp2[i]) for i in range(dx+1)],
        "yP":[int(yp2[i]) for i in range(dy+1)],
        "a":[int(aa2[i]) for i in range(aa2.degree()+1)],
        "b":[int(bb2[i]) for i in range(bb2.degree()+1)],
    }
    print(
        "Q32POLYSEC_ADAPTED|"
        f"profile={dx},{dy}|shift={int(shift)}|scale={int(lam)}|"
        f"xP_next={int(xp2[dx-1])}|"
        f"xP_eq={int(xp2[dx-2]==xp2[dx-3]) if dx>=3 and q2 and q3 else -1}|"
        "status=PASS",
        flush=True,
    )

OUT=args.output.resolve() if args.output else LOCAL/f"q32-d12-polynomial-sections-mod-{p}.json"
OUT.write_text(json.dumps({
    "schema":"elkies-k3.h3-q32-d12-polynomial-sections-modp.v1",
    "status":"PASS_Q32_D12_POLYNOMIAL_SECTIONS_MODP" if len(sections)==26 else "PARTIAL_Q32_D12_POLYNOMIAL_SECTIONS_MODP",
    "prime":int(p),
    "count":len(sections),
    "canonical_model":{
        "a":[int(a[i]) for i in range(7)],
        "b":[int(b[i]) for i in range(10)],
    },
    "sections":[
        {
            "profile":profile,
            "x":[int(x[i]) for i in range(5)],
            "y":[int(y[i]) for i in range(7)],
        }
        for x,y,profile in sections
    ],
    "height1_candidates":[
        {
            "profile":profile,
            "x":[int(x[i]) for i in range(5)],
            "y":[int(y[i]) for i in range(7)],
        }
        for x,y,profile,_ in double_closed
    ],
    "point_adapted":adapted,
},indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q32POLYSEC_RESULT|"
    f"prime={p}|sections={len(sections)}|height1={len(double_closed)}|"
    f"status={'PASS_26_AND_HEIGHT1' if len(sections)==26 and len(double_closed)==2 else 'PARTIAL'}",
    flush=True,
)
