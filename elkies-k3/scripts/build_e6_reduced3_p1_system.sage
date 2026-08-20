from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="E6 P1 reduced3: eliminate b2,b3,b4,b5,b6,b7 analytically.")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--export",default=None)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("choose good prime")

names=["a1","a2","a3","a4",
       "lam","mu","s0","s1","sl","sm",
       "x10","x11","x12",
       "y10","y11","y12","y13","y14"]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
RF=FractionField(R)
Rt=PolynomialRing(RF,"t")
t=Rt.gen()

a1,a2,a3,a4=[RF(d[f"a{i}"]) for i in range(1,5)]
lam,mu,s0,s1,sl,sm=[RF(d[n]) for n in ["lam","mu","s0","s1","sl","sm"]]

# I4 at t=0.
a0=-3*s0**2
b0=2*s0**3
b1=-s0*a1
b2=a1**2/(12*s0)-s0*a2
b3=(a1**3 + 36*a1*a2*s0**2 - 216*a3*s0**4)/(216*s0**3)

# A(1)=-3*s1^2.
a5=-3*s1**2-(a0+a1+a2+a3+a4)

# Temporary variables B4,B5,B6 until all I4@1 conditions are solved.
S=PolynomialRing(RF,["B4","B5","B6"])
B4,B5,B6=S.gens()
St=PolynomialRing(S,"t")
tt=St.gen()

def lift(x):
    return S(RF(x))

AA=sum(lift(c)*tt**i for i,c in enumerate([a0,a1,a2,a3,a4,a5]))

# B(1)=2*s1^3 solves B7.
B7=(2*lift(s1)**3
    -(lift(b0)+lift(b1)+lift(b2)+lift(b3)+B4+B5+B6)
    -1)

BB=(lift(b0)+lift(b1)*tt+lift(b2)*tt**2+lift(b3)*tt**3
    +B4*tt**4+B5*tt**5+B6*tt**6+B7*tt**7+tt**8)

# B'(1)+s1*A'(1)=0 solves B6 linearly in B4,B5.
Efirst=S(BB.derivative(tt)(1)+lift(s1)*AA.derivative(tt)(1))
c6=S(Efirst.derivative(B6))
if c6==0:
    raise RuntimeError("expected nonzero B6 coefficient")

B6_expr=S(-Efirst.subs({B6:0})/c6)
B7_expr=S(B7.subs({B6:B6_expr}))

# Substitute B6 and enforce Delta''(1)=Delta'''(1)=0.
BB2=(lift(b0)+lift(b1)*tt+lift(b2)*tt**2+lift(b3)*tt**3
     +B4*tt**4+B5*tt**5+B6_expr*tt**6+B7_expr*tt**7+tt**8)

DeltaS=-16*(4*AA**3+27*BB2**2)
D1S=DeltaS.derivative(tt)
D2S=D1S.derivative(tt)
D3S=D2S.derivative(tt)
E1=S(D2S(1))
E2=S(D3S(1))

c11=S(E1.derivative(B4)); c12=S(E1.derivative(B5))
c21=S(E2.derivative(B4)); c22=S(E2.derivative(B5))

for nm,c in [("c11",c11),("c12",c12),("c21",c21),("c22",c22)]:
    if c.degree(B4)>0 or c.degree(B5)>0:
        raise RuntimeError(f"{nm} unexpectedly depends on B4/B5")

c11=RF(c11); c12=RF(c12); c21=RF(c21); c22=RF(c22)
k1=RF(E1.subs({B4:0,B5:0}))
k2=RF(E2.subs({B4:0,B5:0}))

detlin=c11*c22-c12*c21
if detlin==0:
    raise RuntimeError("I4@1 B4/B5 system singular identically")

b4=RF((-k1*c22+c12*k2)/detlin)
b5=RF((-c11*k2+c21*k1)/detlin)
b6=RF(B6_expr.subs({B4:S(b4),B5:S(b5)}))
b7=RF(B7_expr.subs({B4:S(b4),B5:S(b5)}))

A=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
B=(b0+b1*t+b2*t**2+b3*t**3+b4*t**4+b5*t**5+b6*t**6+b7*t**7+t**8)
Ap=A.derivative(t)
Bp=B.derivative(t)

# P1 reductions.
x10,x11,x12=[RF(d[n]) for n in ["x10","x11","x12"]]
y10,y11,y12,y13,y14=[RF(d[n]) for n in ["y10","y11","y12","y13","y14"]]
x13=s1-x10-x11-x12
y15=-(y10+y11+y12+y13+y14)
X=x10+x11*t+x12*t**2+x13*t**3
Y=y10+y11*t+y12*t**2+y13*t**3+y14*t**4+y15*t**5

eqs=[]
tags=[]

def add(tag,e):
    e=RF(e)
    if e!=0:
        num=R(e.numerator())
        if num!=0:
            tags.append(tag)
            eqs.append(num)

for pref,a,s in [("I2_lam",lam,sl),("I2_mu",mu,sm)]:
    add(pref+"_A",A(a)+3*s**2)
    add(pref+"_B",B(a)-2*s**3)
    add(pref+"_d1",Bp(a)+s*Ap(a))

S1=Y**2-X**3-A*X-B
for k in range(S1.degree()+1):
    add(f"P1_{k}",S1[k])

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"E6RED3|field={field}|vars={R.ngens()}|eqs={len(eqs)}|naive_dim={R.ngens()-len(eqs)}",flush=True)
print("E6RED3|eliminated=a0,b0,b1,b2,b3,a5,b4,b5,b6,b7,x13,x14,y15,y16|count=14",flush=True)
print(f"E6RED3|b4_num_terms={len(R(b4.numerator()).monomials())}|b4_den_terms={len(R(b4.denominator()).monomials())}",flush=True)
print(f"E6RED3|b5_num_terms={len(R(b5.numerator()).monomials())}|b5_den_terms={len(R(b5.denominator()).monomials())}",flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6RED3_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

print("E6RED3|expected_geometry_dim=3; remaining equations have dependencies",flush=True)
print("E6RED3|saturate=all substitution denominators and s0*s1*sl*sm*lam*(lam-1)*mu*(mu-1)*(lam-mu)",flush=True)

if args.export:
    if args.p==0:
        raise SystemExit("--export requires --p")
    out=Path(args.export)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"E6RED3|export={out}",flush=True)
