from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Triangularly reduced E6 MW3 P1 construction.")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--export",default=None)
ap.add_argument("--show",action="store_true")
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("choose good prime !=2,3,79")

# Keep only genuinely free variables after obvious local/Tate and section substitutions.
names=["a1","a2","a3","a4",
       "b2","b3","b4","b5",
       "lam","mu","s0","s1","sl","sm",
       "x10","x11","x12",
       "y10","y11","y12","y13","y14"]
R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
Rt=PolynomialRing(R,"t"); t=Rt.gen()

a1,a2,a3,a4=[d[f"a{i}"] for i in range(1,5)]
b2,b3,b4,b5=[d[f"b{i}"] for i in range(2,6)]
lam,mu,s0,s1,sl,sm=[d[n] for n in ["lam","mu","s0","s1","sl","sm"]]

# Immediate eliminations from I4@0:
a0=-3*s0**2
b0= 2*s0**3
b1=-s0*a1

# I4@1 A-condition solves a5.
a5=-3*s1**2-(a0+a1+a2+a3+a4)

# Temporarily introduce b6 and derive b7 from B(1)=2s1^3.
S=PolynomialRing(R,["B6"])
B6=S.gen()
St=PolynomialRing(S,"t"); tt=St.gen()
lift=lambda x:S(x)

AA=lift(a0)+lift(a1)*tt+lift(a2)*tt**2+lift(a3)*tt**3+lift(a4)*tt**4+lift(a5)*tt**5
b7=2*lift(s1)**3-(lift(b0)+lift(b1)+lift(b2)+lift(b3)+lift(b4)+lift(b5)+B6)-1
BB=lift(b0)+lift(b1)*tt+lift(b2)*tt**2+lift(b3)*tt**3+lift(b4)*tt**4+lift(b5)*tt**5+B6*tt**6+b7*tt**7+tt**8

# I4@1 first derivative condition solves b6 linearly.
E=S(BB.derivative(tt)(1)+lift(s1)*AA.derivative(tt)(1))
coef=E.monomial_coefficient(B6)
const=E.subs({B6:0})
if coef==0:
    raise RuntimeError("expected linear b6 coefficient")
b6_R=R(-const/coef)
b7_R=R(b7.subs({B6:S(b6_R)}))

# Final A,B in reduced ring.
A=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
B=b0+b1*t+b2*t**2+b3*t**3+b4*t**4+b5*t**5+b6_R*t**6+b7_R*t**7+t**8
Ap=A.derivative(t); Bp=B.derivative(t)
Delta=-16*(4*A**3+27*B**2)
D1=Delta.derivative(t); D2=D1.derivative(t); D3=D2.derivative(t)

# P1: x4=0,y6=0 and t=1 singular-hit conditions solve x3,y5.
x10,x11,x12=[d[n] for n in ["x10","x11","x12"]]
y10,y11,y12,y13,y14=[d[n] for n in ["y10","y11","y12","y13","y14"]]
x13=s1-x10-x11-x12
y15=-(y10+y11+y12+y13+y14)

X=x10+x11*t+x12*t**2+x13*t**3
Y=y10+y11*t+y12*t**2+y13*t**3+y14*t**4+y15*t**5

eqs=[]; tags=[]
def add(tag,e):
    e=R(e)
    if e!=0:
        tags.append(tag); eqs.append(e)

# Remaining I4 exactness at 0 and 1.
# A,B,first derivative equations are already built into substitutions.
add("I4_0_D2",D2(0))
add("I4_0_D3",D3(0))
add("I4_1_D2",D2(1))
add("I4_1_D3",D3(1))

# I2 local equations at lambda and mu.
for pref,a,s in [("I2_lam",lam,sl),("I2_mu",mu,sm)]:
    add(pref+"_A",A(a)+3*s**2)
    add(pref+"_B",B(a)-2*s**3)
    add(pref+"_d1",Bp(a)+s*Ap(a))

# P1 section identity.
S1=Y**2-X**3-A*X-B
for k in range(S1.degree()+1):
    add(f"P1_{k}",S1[k])

# P1 IV* leading constraints are already x4=0,y6=0.
# P1 I4@1 singular hit already X(1)=s1,Y(1)=0.
# Its finite identity at I4@0 and I2 fibers is left unconstrained, as desired.

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"E6RED|stage=p1|field={field}|vars={R.ngens()}|eqs={len(eqs)}|naive_dim={R.ngens()-len(eqs)}",flush=True)
print(f"E6RED|eliminated=a0,b0,b1,a5,b6,b7,x13,x14,y15,y16|count=10",flush=True)
print(f"E6RED|A={A}",flush=True)
print(f"E6RED|B={B}",flush=True)
print(f"E6RED|X={X}",flush=True)
print(f"E6RED|Y={Y}",flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6RED_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)
    if args.show:
        print(f"E6RED_FORMULA|i={i}|eq={e}",flush=True)

print("E6RED|expected_geometry_dim=3 before deeper component labels",flush=True)
print("E6RED|saturate=lam*(lam-1)*mu*(mu-1)*(lam-mu)*s0*s1*sl*sm and exact fiber residuals",flush=True)

if args.export:
    if args.p==0:
        raise SystemExit("--export requires --p")
    out=Path(args.export); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"E6RED|export={out}",flush=True)
