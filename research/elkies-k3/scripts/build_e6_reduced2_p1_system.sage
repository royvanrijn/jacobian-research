from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Further reduced E6 P1 system using closed-form I4 conditions at t=0.")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--export",default=None)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("choose good prime")

names=["a1","a2","a3","a4",
       "b4","b5",
       "lam","mu","s0","s1","sl","sm",
       "x10","x11","x12",
       "y10","y11","y12","y13","y14"]
R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t")
t=Rt.gen()

a1,a2,a3,a4=[d[f"a{i}"] for i in range(1,5)]
b4,b5=d["b4"],d["b5"]
lam,mu,s0,s1,sl,sm=[d[n] for n in ["lam","mu","s0","s1","sl","sm"]]

a0=-3*s0**2
b0=2*s0**3
b1=-s0*a1

b2=a1**2/(12*s0)-s0*a2
b3=(a1**3+36*a1*a2*s0**2-216*a3*s0**4)/(216*s0**3)

a5=-3*s1**2-(a0+a1+a2+a3+a4)

RF=FractionField(R)
S=PolynomialRing(RF,["B6"])
B6=S.gen()
St=PolynomialRing(S,"t")
tt=St.gen()
lift=lambda x:S(RF(x))

AA=sum(lift(c)*tt**i for i,c in enumerate([a0,a1,a2,a3,a4,a5]))
b7=2*lift(s1)**3-(lift(b0)+lift(b1)+lift(b2)+lift(b3)+lift(b4)+lift(b5)+B6)-1
BB=(lift(b0)+lift(b1)*tt+lift(b2)*tt**2+lift(b3)*tt**3+
    lift(b4)*tt**4+lift(b5)*tt**5+B6*tt**6+b7*tt**7+tt**8)

E=S(BB.derivative(tt)(1)+lift(s1)*AA.derivative(tt)(1))
coef=E.monomial_coefficient(B6)
if coef==0:
    raise RuntimeError("b6 coefficient vanished")
b6=RF(-E.subs({B6:0})/coef)
b7=RF(b7.subs({B6:S(b6)}))

A=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
B=b0+b1*t+b2*t**2+b3*t**3+b4*t**4+b5*t**5+b6*t**6+b7*t**7+t**8

Ap=A.derivative(t)
Bp=B.derivative(t)
Delta=-16*(4*A**3+27*B**2)
D1=Delta.derivative(t)
D2=D1.derivative(t)
D3=D2.derivative(t)

x10,x11,x12=[d[n] for n in ["x10","x11","x12"]]
y10,y11,y12,y13,y14=[d[n] for n in ["y10","y11","y12","y13","y14"]]

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

add("I4_1_D2",D2(1))
add("I4_1_D3",D3(1))

for pref,a,s in [("I2_lam",lam,sl),("I2_mu",mu,sm)]:
    add(pref+"_A",A(a)+3*s**2)
    add(pref+"_B",B(a)-2*s**3)
    add(pref+"_d1",Bp(a)+s*Ap(a))

S1=Y**2-X**3-A*X-B
for k in range(S1.degree()+1):
    add(f"P1_{k}",S1[k])

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"E6RED2|field={field}|vars={R.ngens()}|eqs={len(eqs)}|naive_dim={R.ngens()-len(eqs)}",flush=True)
print("E6RED2|eliminated_total=12|new=b2,b3",flush=True)
print(f"E6RED2|b2={b2}",flush=True)
print(f"E6RED2|b3={b3}",flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6RED2_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

for tag,e in zip(tags[:2],eqs[:2]):
    for var in (b4,b5):
        print(f"E6RED2_DIAG|eq={tag}|var={var}|degree={e.degree(var)}",flush=True)

print("E6RED2|expected_geometry_dim=3; equation count has algebraic dependencies",flush=True)
print("E6RED2|saturate=s0*s1*sl*sm*lam*(lam-1)*mu*(mu-1)*(lam-mu)",flush=True)

if args.export:
    if args.p==0:
        raise SystemExit("--export requires --p")
    out=Path(args.export)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            num=e.numerator()
            h.write(str(num).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"E6RED2|export={out}",flush=True)
