from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="E6 P1 triangular reduction through P1_7 and P1_6.")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--export",default=None)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("choose good prime")

# Start from the section-first variables.
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

a0=-3*s0**2
b0=2*s0**3
b1=-s0*a1
b2=a1**2/(12*s0)-s0*a2
b3=(a1**3+36*a1*a2*s0**2-216*a3*s0**4)/(216*s0**3)
a5=-3*s1**2-(a0+a1+a2+a3+a4)

# Solve b6,b7 from B(1), B'(1)+s1*A'(1)=0.
S=PolynomialRing(RF,["B6"]); B6=S.gen()
St=PolynomialRing(S,"t"); tt=St.gen()
lift=lambda z:S(RF(z))

AA=sum(lift(c)*tt**i for i,c in enumerate([a0,a1,a2,a3,a4,a5]))
B7=2*lift(s1)**3-(lift(b0)+lift(b1)+lift(b2)+lift(b3)+lift(b4)+lift(b5)+B6)-1
BB=(lift(b0)+lift(b1)*tt+lift(b2)*tt**2+lift(b3)*tt**3+
    lift(b4)*tt**4+lift(b5)*tt**5+B6*tt**6+B7*tt**7+tt**8)
E=S(BB.derivative(tt)(1)+lift(s1)*AA.derivative(tt)(1))
c6=E.derivative(B6)
b6=RF((-E.subs({B6:0}))/c6)
b7=RF(B7.subs({B6:S(b6)}))

A=a0+a1*t+a2*t**2+a3*t**3+a4*t**4+a5*t**5
B=b0+b1*t+b2*t**2+b3*t**3+b4*t**4+b5*t**5+b6*t**6+b7*t**7+t**8
Ap=A.derivative(t); Bp=B.derivative(t)

x2=s1-x0-x1
y3=-1-y0-y1-y2
X=x0+x1*t+x2*t**2
Y=y0+y1*t+y2*t**2+y3*t**3+t**4

S1=Y**2-X**3-A*X-B

# First triangular solve: P1_7 for b4.
E7=RF(S1[7])
c7=E7.derivative(b4)
if c7==0 or c7.derivative(b4)!=0:
    raise RuntimeError("P1_7 not affine-linear in b4")
b4_expr=RF(-E7.subs({b4:0})/c7)

def subpoly(poly,subs):
    return sum(RF(poly[k].subs(subs))*t**k for k in range(poly.degree()+1))

subs1={b4:b4_expr}
A1=subpoly(A,subs1); B1=subpoly(B,subs1)
X1=subpoly(X,subs1); Y1=subpoly(Y,subs1)
S1a=Y1**2-X1**3-A1*X1-B1

# Second triangular solve: P1_6 for b5.
E6=RF(S1a[6])
c6p=E6.derivative(b5)
if c6p==0 or c6p.derivative(b5)!=0:
    raise RuntimeError("P1_6 not affine-linear in b5 after b4 elimination")
b5_expr=RF(-E6.subs({b5:0})/c6p)

subs2={b5:b5_expr}
A2=subpoly(A1,subs2); B2=subpoly(B1,subs2)
X2=subpoly(X1,subs2); Y2=subpoly(Y1,subs2)
Sfinal=Y2**2-X2**3-A2*X2-B2

Ap=A2.derivative(t); Bp=B2.derivative(t)
Delta=-16*(4*A2**3+27*B2**2)
D1=Delta.derivative(t); D2=D1.derivative(t); D3=D2.derivative(t)

# New ring without b4,b5.
new_names=[n for n in names if n not in ("b4","b5")]
RR=PolynomialRing(K,new_names,order="degrevlex")
rrd=RR.gens_dict()

def to_RR(fr):
    fr=RF(fr)
    mp={d[n]:rrd[n] for n in new_names}
    num=RR(fr.numerator().subs(mp))
    den=RR(fr.denominator().subs(mp))
    return num,den

eqs=[]; tags=[]
def add(tag,e):
    num,den=to_RR(e)
    if num!=0:
        tags.append(tag); eqs.append(num)

# Fiber conditions.
add("I4_1_D2",D2(1))
add("I4_1_D3",D3(1))
for pref,a,s in [("I2_lam",lam,sl),("I2_mu",mu,sm)]:
    aa=RF(a); ss=RF(s)
    add(pref+"_A",A2(aa)+3*ss**2)
    add(pref+"_B",B2(aa)-2*ss**3)
    add(pref+"_d1",Bp(aa)+ss*Ap(aa))

# P1_7,P1_6 vanish by construction. Keep 0..5.
for k in range(6):
    add(f"P1_{k}",Sfinal[k])

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"E6TRI2|field={field}|vars={RR.ngens()}|eqs={len(eqs)}|naive_dim={RR.ngens()-len(eqs)}",flush=True)
print(f"E6TRI2|b4_num_terms={len(R(b4_expr.numerator()).monomials())}|b4_den_terms={len(R(b4_expr.denominator()).monomials())}",flush=True)
print(f"E6TRI2|b5_num_terms={len(R(b5_expr.numerator()).monomials())}|b5_den_terms={len(R(b5_expr.denominator()).monomials())}",flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6TRI2_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

# Find next affine-linear variable in P1_5, then P1_4.
for target in ("P1_5","P1_4"):
    if target not in tags: continue
    e=eqs[tags.index(target)]
    vals=[]
    for nm in ["y2","x1","y1","a4","a3","a2","x0","y0"]:
        if nm not in rrd: continue
        v=rrd[nm]
        if e.degree(v)==1:
            coeff=e.derivative(v)
            vals.append((nm,len(coeff.monomials()),len(e.monomials())))
    print(f"E6TRI2_NEXT|eq={target}|linear={vals}",flush=True)

print("E6TRI2|expected_geometry_dim=3",flush=True)

if args.export:
    if args.p==0: raise SystemExit("--export requires --p")
    out=Path(args.export); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(new_names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"E6TRI2|export={out}",flush=True)
