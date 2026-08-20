from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="E6 section-first P1 with automatic triangular elimination from P1_7.")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--export",default=None)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("choose good prime")

# Start from same variables as section-first, then eliminate one selected variable.
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

# Solve b6,b7 from B(1), first derivative, retaining b4,b5.
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
Delta=-16*(4*A**3+27*B**2)
D1=Delta.derivative(t); D2=D1.derivative(t); D3=D2.derivative(t)

x2=s1-x0-x1
y3=-1-y0-y1-y2
X=x0+x1*t+x2*t**2
Y=y0+y1*t+y2*t**2+y3*t**3+t**4

S1=Y**2-X**3-A*X-B
P17=RF(S1[7])

# Choose the cleanest affine-linear variable among candidates.
candidates=["y2","b4","b5","x1","y1"]
choices=[]
for nm in candidates:
    v=RF(d[nm])
    # derivative constant wrt v iff affine-linear
    coeff=P17.derivative(v)
    if coeff==0:
        continue
    if coeff.derivative(v)!=0:
        continue
    const=P17.subs({v:0})
    score=(len(R(coeff.numerator()).monomials())+len(R(coeff.denominator()).monomials()),
           len(R(const.numerator()).monomials()),
           nm)
    choices.append((score,nm,v,coeff,const))

if not choices:
    raise RuntimeError("no affine-linear elimination variable found in P1_7")

choices.sort(key=lambda z:z[0])
score,elim_name,elim_var,coeff,const=choices[0]
elim_expr=RF(-const/coeff)

print(f"E6TRI|elim_from=P1_7|variable={elim_name}|score={score}",flush=True)
print(f"E6TRI|elim_expr_num_terms={len(R(elim_expr.numerator()).monomials())}|den_terms={len(R(elim_expr.denominator()).monomials())}",flush=True)

# Substitute eliminated variable coefficient-wise.  A,B,X,Y live in RF[t],
# while elim_var lives in RF; direct Polynomial.subs() rejects that parent mix.
subs={elim_var:elim_expr}

def sub_coeffs(poly):
    return sum(
        RF(poly[k].subs(subs))*t**k
        for k in range(poly.degree()+1)
    )

A=sub_coeffs(A)
B=sub_coeffs(B)
X=sub_coeffs(X)
Y=sub_coeffs(Y)

Ap=A.derivative(t)
Bp=B.derivative(t)
Delta=-16*(4*A**3+27*B**2)
D1=Delta.derivative(t)
D2=D1.derivative(t)
D3=D2.derivative(t)
S1=Y**2-X**3-A*X-B

# New polynomial ring without eliminated variable.
new_names=[n for n in names if n!=elim_name]
RR=PolynomialRing(K,new_names,order="degrevlex")
rrd=RR.gens_dict()

# Coercion map by name.
def to_RR(fr):
    fr=RF(fr)
    num=fr.numerator()
    den=fr.denominator()
    # eliminated variable is absent after substitution
    mp={d[n]:rrd[n] for n in new_names}
    num2=RR(num.subs(mp))
    den2=RR(den.subs(mp))
    return num2,den2

eqs=[]; tags=[]
def add(tag,e):
    num,den=to_RR(e)
    if num!=0:
        tags.append(tag); eqs.append(num)

# remaining fiber equations
add("I4_1_D2",D2(1))
add("I4_1_D3",D3(1))
for pref,a,s in [("I2_lam",RF(lam).subs(subs),RF(sl).subs(subs)),
                 ("I2_mu",RF(mu).subs(subs),RF(sm).subs(subs))]:
    add(pref+"_A",A(a)+3*s**2)
    add(pref+"_B",B(a)-2*s**3)
    add(pref+"_d1",Bp(a)+s*Ap(a))

# P1_7 is identically zero after substitution; keep lower coefficients.
for k in range(0,7):
    add(f"P1_{k}",S1[k])

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"E6TRI|field={field}|vars={RR.ngens()}|eqs={len(eqs)}|naive_dim={RR.ngens()-len(eqs)}",flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6TRI_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

# Diagnose next affine-linear possibilities in P1_6 and P1_5.
for target in ("P1_6","P1_5"):
    if target not in tags:
        continue
    e=eqs[tags.index(target)]
    vals=[]
    for nm in ["y2","b4","b5","x1","y1","a4","a3"]:
        if nm not in rrd: continue
        v=rrd[nm]
        deg=e.degree(v)
        if deg==1:
            coeff=e.derivative(v)
            vals.append((nm,len(coeff.monomials()),e.degree()))
    print(f"E6TRI_NEXT|eq={target}|linear={vals}",flush=True)

print("E6TRI|expected_geometry_dim=3",flush=True)

if args.export:
    if args.p==0: raise SystemExit("--export requires --p")
    out=Path(args.export); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(new_names)+"\n"+str(args.p)+"\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(eqs) else "\n")
    print(f"E6TRI|export={out}",flush=True)
