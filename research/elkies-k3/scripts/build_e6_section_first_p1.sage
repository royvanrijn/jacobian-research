from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="E6 P1 section-first reduction; avoid b4/b5 expression swell.")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--sign",type=int,choices=[-1,1],default=1)
ap.add_argument("--export",default=None)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("choose good prime")

eps=K(args.sign)

# RED2 fiber variables, but reduced section variables.
names=["a1","a2","a3","a4","b4","b5",
       "lam","mu","s0","s1","sl","sm",
       "x0","x1","y0","y1","y2"]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
RF=FractionField(R)
Rt=PolynomialRing(RF,"t"); t=Rt.gen()

a1,a2,a3,a4,b4,b5=[RF(d[n]) for n in ["a1","a2","a3","a4","b4","b5"]]
lam,mu,s0,s1,sl,sm=[RF(d[n]) for n in ["lam","mu","s0","s1","sl","sm"]]
x0,x1,y0,y1,y2=[RF(d[n]) for n in ["x0","x1","y0","y1","y2"]]

# I4@0 substitutions.
a0=-3*s0**2
b0= 2*s0**3
b1=-s0*a1
b2=a1**2/(12*s0)-s0*a2
b3=(a1**3+36*a1*a2*s0**2-216*a3*s0**4)/(216*s0**3)

# I4@1 A-condition.
a5=-3*s1**2-(a0+a1+a2+a3+a4)

# Keep b4,b5 free; solve b6,b7 from B(1), first derivative at 1.
S=PolynomialRing(RF,["B6"])
B6=S.gen()
St=PolynomialRing(S,"t"); tt=St.gen()
lift=lambda z:S(RF(z))
AA=sum(lift(c)*tt**i for i,c in enumerate([a0,a1,a2,a3,a4,a5]))
B7=2*lift(s1)**3-(lift(b0)+lift(b1)+lift(b2)+lift(b3)+lift(b4)+lift(b5)+B6)-1
BB=(lift(b0)+lift(b1)*tt+lift(b2)*tt**2+lift(b3)*tt**3+
    lift(b4)*tt**4+lift(b5)*tt**5+B6*tt**6+B7*tt**7+tt**8)
E=S(BB.derivative(tt)(1)+lift(s1)*AA.derivative(tt)(1))
c6=E.derivative(B6)
if c6==0: raise RuntimeError("B6 coefficient vanished")
b6=RF((-E.subs({B6:0}))/c6)
b7=RF(B7.subs({B6:S(b6)}))

A=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
B=b0+b1*t+b2*t**2+b3*t**3+b4*t**4+b5*t**5+b6*t**6+b7*t**7+t**8
Ap=A.derivative(t); Bp=B.derivative(t)
Delta=-16*(4*A**3+27*B**2)
D1=Delta.derivative(t); D2=D1.derivative(t); D3=D2.derivative(t)

# Top section equations:
# P1_10 => y5=0
# P1_9  => x3=0
# P1_8  => y4^2=1; choose eps=+/-1.
# X(1)=s1 and Y(1)=0 then give x2,y3.
x2=s1-x0-x1
y4=RF(eps)
y3=-y4-y0-y1-y2

X=x0+x1*t+x2*t**2
Y=y0+y1*t+y2*t**2+y3*t**3+y4*t**4

eqs=[]; tags=[]
def add(tag,e):
    e=RF(e)
    if e!=0:
        num=R(e.numerator())
        if num!=0:
            tags.append(tag); eqs.append(num)

# Remaining exact I4@1 conditions.
add("I4_1_D2",D2(1))
add("I4_1_D3",D3(1))

# I2 fibers.
for pref,a,s in [("I2_lam",lam,sl),("I2_mu",mu,sm)]:
    add(pref+"_A",A(a)+3*s**2)
    add(pref+"_B",B(a)-2*s**3)
    add(pref+"_d1",Bp(a)+s*Ap(a))

# Remaining P1 coefficients now only degree <=8.
S1=Y**2-X**3-A*X-B
for k in range(S1.degree()+1):
    add(f"P1_{k}",S1[k])

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"E6SF|field={field}|sign={args.sign}|vars={R.ngens()}|eqs={len(eqs)}|naive_dim={R.ngens()-len(eqs)}",flush=True)
print("E6SF|section_reduction=P1_10=>y5=0;P1_9=>x3=0;P1_8=>y4=sign",flush=True)
print(f"E6SF|X={X}",flush=True)
print(f"E6SF|Y={Y}",flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6SF_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

# Diagnose triangular variables in high remaining P1 coefficients.
for tag,e in zip(tags,eqs):
    if tag in ("P1_7","P1_6","P1_5"):
        print("E6SF_DIAG|eq="+tag+"|"+
              "|".join(f"{v}:deg{e.degree(d[v])}" for v in ["y2","b4","b5","x1","y1"]),
              flush=True)

print("E6SF|expected_geometry_dim=3",flush=True)
print("E6SF|saturate=s0*s1*sl*sm*lam*(lam-1)*mu*(mu-1)*(lam-mu)",flush=True)

if args.export:
    if args.p==0: raise SystemExit("--export requires --p")
    out=Path(args.export); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"E6SF|export={out}",flush=True)
