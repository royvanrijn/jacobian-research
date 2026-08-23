#!/usr/bin/env sage -python
"""
Recover the *actual-twist* low-degree polynomial sections of the q32/D12 child.

Important correction to the previous probe:
the canonical j-map model y^2=x^3-3*a*x+2*b is defined only up to quadratic
twist.  The q32 child itself determines the twist.  We recover it directly
from the original modular Jacobian A(V),B(V).

After the same canonical base transformation used for j, write
    A_actual / (-3a) = w^2,
    B_actual / ( 2b) = w^3.
Then w = D*g(u)^2.  Removing the rational square g gives the actual minimal
constant twist
    y^2 = x^3 - 3 D^2 a x + 2 D^3 b.

For this model the nontrivial zero-pole D12 component classes must occur in
the deg(x)=3 branch.  The lattice predicts exactly 8 such sections:
  2 of height 1 (spinor correction 3),
  6 of height 3 (vector correction 1).
The unique +/- height-1 pair is characterized by its double also being a
zero-pole polynomial section, so no MW-coordinate identification is needed.
"""

import argparse, json
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ, EllipticCurve

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
RV=PolynomialRing(F,"V"); V=RV.gen(); KV=RV.fraction_field()

def rf(rec):
    return KV(RV([F(v) for v in rec["num"]]))/KV(
        RV([F(v) for v in rec["den"]])
    )

Aorig=rf(d["jacobian_A"])
Borig=rf(d["jacobian_B"])

# -------------------------------------------------------------------------
# 1. Canonical base coordinate from j.
# -------------------------------------------------------------------------
J=KV(F(6912))*Aorig**3/(KV(F(4))*Aorig**3+KV(F(27))*Borig**2)
N=RV(J.numerator()); Den=RV(J.denominator())
assert N.gcd(Den).degree()==0 and N.degree()==Den.degree()==18

e8=[f for f,e in Den.factor() if int(e)==8 and f.degree()==1]
assert len(e8)==1
f=e8[0]
r=-f[0]/f[1]

RT=PolynomialRing(F,"T"); T=RT.gen()
def invpoly(poly):
    return sum(
        F(poly[i])*(F(r)*T+1)**i*T**(18-i)
        for i in range(poly.degree()+1)
    )

P=invpoly(N)
Q=invpoly(Den)
g=P.gcd(Q)
if g.degree()>0:
    P//=g
    Q//=g
assert P.degree()==18 and Q.degree()==10

lc=Q.leading_coefficient()
P/=lc
Q/=lc
assert Q.is_monic()

center=-Q[9]/F(10)
RS=PolynomialRing(F,"S"); S=RS.gen()
P1=RS(P(S+center))
Q1=RS(Q(S+center))
assert Q1.is_monic() and Q1[9]==0

assert Q1[8] and Q1[7]
base_scale=Q1[7]/Q1[8]

RU=PolynomialRing(F,"u"); u=RU.gen(); KU=RU.fraction_field()
P2=RU(P1(base_scale*u))
Q2=RU(Q1(base_scale*u))
lc2=Q2.leading_coefficient()
P2/=lc2
Q2/=lc2
assert Q2.degree()==10 and Q2.is_monic()
assert Q2[9]==0 and Q2[8]==Q2[7]

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
assert a.is_monic() and b.is_monic()
assert (a**3-b**2).degree()==10

print(
    "Q32TWIST_BASE|"
    f"prime={p}|I8star={int(r)}|center={int(center)}|"
    f"scale={int(base_scale)}|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 2. Recover the quadratic twist from the *actual* q32 Jacobian.
# -------------------------------------------------------------------------
vmap=KU(F(r)) + KU(1)/(KU(F(base_scale))*KU(u)+KU(F(center)))

def eval_rational_at(value,arg):
    num=RV(value.numerator())
    den=RV(value.denominator())
    return KU(num(arg))/KU(den(arg))

Aeval=eval_rational_at(Aorig,vmap)
Beval=eval_rational_at(Borig,vmap)

Acan=KU(-F(3)*a)
Bcan=KU(F(2)*b)
cA=Aeval/Acan
cB=Beval/Bcan
w=cB/cA

assert cA==w**2
assert cB==w**3

wn=RU(w.numerator())
wd=RU(w.denominator())
nf=wn.factor()
df=wd.factor()

square_part=KU.one()
for fac,e in nf:
    assert int(e)%2==0, ("nonconstant twist numerator",fac,e)
    square_part*=KU(fac.monic())**(int(e)//2)
for fac,e in df:
    assert int(e)%2==0, ("nonconstant twist denominator",fac,e)
    square_part/=KU(fac.monic())**(int(e)//2)

Dfun=w/(square_part**2)
Dnum=RU(Dfun.numerator())
Dden=RU(Dfun.denominator())
assert Dnum.degree()<=0 and Dden.degree()<=0
twist=F(Dnum[0])/F(Dden[0])
assert twist

At=-F(3)*twist**2*a
Bt= F(2)*twist**3*b

# Verify actual curve differs from the twist model only by the rational square
# scaling square_part.
assert Aeval == KU(square_part**4)*KU(At)
assert Beval == KU(square_part**6)*KU(Bt)

print(
    "Q32TWIST_RECOVER|"
    f"prime={p}|twist={int(twist)}|"
    f"square_num_deg={RU(square_part.numerator()).degree()}|"
    f"square_den_deg={RU(square_part.denominator()).degree()}|"
    f"twist_square={int(twist.is_square())}|status=PASS_ACTUAL_TWIST",
    flush=True,
)

# -------------------------------------------------------------------------
# 3. Solve the deg(x)=3 polynomial-section branches on the actual twist.
# -------------------------------------------------------------------------
def solve_dx3(clead):
    names=("ell","x2","x1","x0","inv")
    Sring=PolynomialRing(F,names=names,order="degrevlex")
    ell,x2,x1,x0,inv=Sring.gens()
    K=Sring.fraction_field()
    U=PolynomialRing(K,"z"); z=U.gen()

    AA=U([K(v) for v in At.list()])
    BB=U([K(v) for v in Bt.list()])
    x=K(F(clead))*z**3+K(x2)*z**2+K(x1)*z+K(x0)
    R=x**3+AA*x+BB
    assert R[9]==0
    assert R.degree()<=8

    ys={4:K(ell)}
    equations=[Sring((K(ell)**2-K(R[8])).numerator())]

    # Match coefficients 7,6,5,4 recursively.
    for k in range(7,3,-1):
        j=k-4
        known=sum(
            ys[i]*ys[k-i]
            for i in ys
            if (k-i) in ys and i!=4 and (k-i)!=4
        )
        ys[j]=(K(R[k])-known)/(K(2)*ys[4])

    y=sum(ys[i]*z**i for i in range(5))
    residual=y**2-R
    equations += [
        Sring(K(residual[k]).numerator())
        for k in range(4)
    ]
    equations.append(inv*ell-1)

    I=Sring.ideal(equations)
    print(
        "Q32TWIST_DX3|"
        f"lead={int(F(clead))}|vars=5|eqs={len(equations)}|"
        "status=GROEBNER_START",
        flush=True,
    )
    sols=I.variety()
    print(
        "Q32TWIST_DX3|"
        f"lead={int(F(clead))}|solutions={len(sols)}|"
        "status=GROEBNER_PASS",
        flush=True,
    )

    out=[]
    for sol in sols:
        vals={g:F(sol[g]) for g in Sring.gens()}
        if not vals[ell]:
            continue
        xx=RU([vals[x0],vals[x1],vals[x2],F(clead)])
        RR=xx**3+At*xx+Bt

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
        out.append((xx,yy))
    return out

# Leading cancellation:
# c^3 - 3 D^2 c + 2 D^3 = (c-D)^2(c+2D).
sections=[]
sections += solve_dx3(twist)
sections += solve_dx3(-F(2)*twist)

unique={}
for x,y in sections:
    key=(
        tuple(int(x[i]) for i in range(4)),
        tuple(int(y[i]) for i in range(5)),
    )
    unique[key]=(x,y)
sections=list(unique.values())

print(
    "Q32TWIST_DX3_COUNT|"
    f"prime={p}|count={len(sections)}|expected=8|"
    f"status={'PASS_EXPECTED_8' if len(sections)==8 else 'PARTIAL'}",
    flush=True,
)

# -------------------------------------------------------------------------
# 4. The unique +/- height-1 pair is the pair whose double remains integral
#    with deg x<=4, deg y<=6.
# -------------------------------------------------------------------------
E=EllipticCurve(KU,[0,0,0,KU(At),KU(Bt)])

def polynomial_zero_pole(Pnt):
    if Pnt.is_zero():
        return None
    xx,yy=Pnt.xy()
    xd=RU(xx.denominator())
    yd=RU(yy.denominator())
    if xd.degree()>0 or yd.degree()>0:
        return None
    xp=RU(xx.numerator())/F(xd[0])
    yp=RU(yy.numerator())/F(yd[0])
    if xp.degree()>4 or yp.degree()>6:
        return None
    return xp,yp

height1=[]
for x,y in sections:
    Pnt=E(KU(x),KU(y))
    doubled=polynomial_zero_pole(2*Pnt)
    if doubled is not None:
        height1.append((x,y,doubled))

print(
    "Q32TWIST_HEIGHT1|"
    f"candidates={len(height1)}|expected=2|"
    f"status={'PASS_UNIQUE_PAIR' if len(height1)==2 else 'PARTIAL'}",
    flush=True,
)

for i,(x,y,(x2,y2)) in enumerate(height1):
    print(
        "Q32TWIST_HEIGHT1_POINT|"
        f"index={i}|x={','.join(str(int(x[j])) for j in range(4))}|"
        f"y={','.join(str(int(y[j])) for j in range(5))}|"
        f"double_xdeg={x2.degree()}|double_ydeg={y2.degree()}|status=PASS",
        flush=True,
    )

# -------------------------------------------------------------------------
# 5. Sign-independent section-adapted base normalization using x(P).
# -------------------------------------------------------------------------
adapted=None
if len(height1)==2:
    xP=height1[0][0]
    assert xP==height1[1][0]
    assert xP.degree()==3

    # u = v + shift kills the quadratic coefficient of xP.
    shift=-xP[2]/(F(3)*xP[3])

    VR=PolynomialRing(F,"v"); v=VR.gen()
    xp1=VR(xP(v+shift))
    A1=VR(At(v+shift))
    B1=VR(Bt(v+shift))
    assert xp1[2]==0

    # v = lambda*w, choose coeff(w)=constant when possible.
    assert xp1[1] and xp1[0]
    lam=xp1[0]/xp1[1]

    WR=PolynomialRing(F,"w"); ww=WR.gen()
    xp=WR(xp1(lam*ww))
    AA=WR(A1(lam*ww))
    BB=WR(B1(lam*ww))
    assert xp[2]==0 and xp[1]==xp[0]

    # Translate x by the distinguished x-coordinate. This is sign-independent:
    # y^2 = X^3 + c2 X^2 + c4 X + c6.
    c2=F(3)*xp
    c4=F(3)*xp**2+AA
    c6=xp**3+AA*xp+BB
    assert c6.degree()<=8

    adapted={
        "base_shift":int(shift),
        "base_scale":int(lam),
        "x_height1":[int(xp[i]) for i in range(4)],
        "c2":[int(c2[i]) for i in range(c2.degree()+1)],
        "c4":[int(c4[i]) for i in range(c4.degree()+1)],
        "c6":[int(c6[i]) for i in range(c6.degree()+1)],
    }

    print(
        "Q32TWIST_ADAPTED|"
        f"xPdeg={xp.degree()}|xP2={int(xp[2])}|"
        f"xP1eq0={int(xp[1]==xp[0])}|"
        f"c2deg={c2.degree()}|c4deg={c4.degree()}|c6deg={c6.degree()}|"
        "status=PASS_SECTION_ADAPTED",
        flush=True,
    )

OUT=args.output.resolve() if args.output else LOCAL/f"q32-d12-actual-twist-dx3-mod-{p}.json"
OUT.write_text(json.dumps({
    "schema":"elkies-k3.h3-q32-d12-actual-twist-dx3-modp.v1",
    "status":(
        "PASS_Q32_D12_ACTUAL_TWIST_DX3_HEIGHT1"
        if len(sections)==8 and len(height1)==2 and adapted is not None
        else "PARTIAL_Q32_D12_ACTUAL_TWIST_DX3"
    ),
    "prime":int(p),
    "twist":int(twist),
    "canonical_base":{
        "I8star_root":int(r),
        "center":int(center),
        "scale":int(base_scale),
    },
    "twisted_model":{
        "a4":[int(At[i]) for i in range(7)],
        "a6":[int(Bt[i]) for i in range(10)],
    },
    "dx3_sections":[
        {
            "x":[int(x[i]) for i in range(4)],
            "y":[int(y[i]) for i in range(5)],
        }
        for x,y in sections
    ],
    "height1":[
        {
            "x":[int(x[i]) for i in range(4)],
            "y":[int(y[i]) for i in range(5)],
        }
        for x,y,_ in height1
    ],
    "section_adapted":adapted,
},indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q32TWIST_RESULT|"
    f"prime={p}|dx3={len(sections)}|height1={len(height1)}|"
    f"status={'PASS_ACTUAL_TWIST_SECTION_ANCHOR' if len(sections)==8 and len(height1)==2 and adapted is not None else 'PARTIAL'}",
    flush=True,
)
