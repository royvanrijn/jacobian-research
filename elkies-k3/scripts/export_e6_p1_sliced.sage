from sage.all import *
from pathlib import Path
import argparse, random

ap=argparse.ArgumentParser(description="Frozen E6 P1 system after four safe triangular eliminations.")
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--slices",type=int,default=3)
ap.add_argument("--out",required=True)
ap.add_argument("--show",action="store_true")
args=ap.parse_args()

if not is_prime(args.p) or args.p in (2,3,79):
    raise SystemExit("choose good prime != 2,3,79")

K=GF(args.p)
random.seed(int(args.seed))
set_random_seed(args.seed)

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

S=PolynomialRing(RF,["B6"]); B6=S.gen()
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
x2=s1-x0-x1
y3=-1-y0-y1-y2
X=x0+x1*t+x2*t**2
Y=y0+y1*t+y2*t**2+y3*t**3+t**4
S1=Y**2-X**3-A*X-B

def subpoly(poly,subs):
    return sum(RF(poly[k].subs(subs))*t**k for k in range(poly.degree()+1))

E7=RF(S1[7]); c7=E7.derivative(b4)
if c7==0 or c7.derivative(b4)!=0: raise RuntimeError("P1_7 not affine-linear in b4")
b4_expr=RF(-E7.subs({b4:0})/c7)
A1=subpoly(A,{b4:b4_expr}); B1=subpoly(B,{b4:b4_expr})
X1=subpoly(X,{b4:b4_expr}); Y1=subpoly(Y,{b4:b4_expr})
S1a=Y1**2-X1**3-A1*X1-B1

E6=RF(S1a[6]); c6p=E6.derivative(b5)
if c6p==0 or c6p.derivative(b5)!=0: raise RuntimeError("P1_6 not affine-linear in b5")
b5_expr=RF(-E6.subs({b5:0})/c6p)
A2=subpoly(A1,{b5:b5_expr}); B2=subpoly(B1,{b5:b5_expr})
X2=subpoly(X1,{b5:b5_expr}); Y2=subpoly(Y1,{b5:b5_expr})
S2=Y2**2-X2**3-A2*X2-B2

E5=RF(S2[5]); c5=E5.derivative(y2)
if c5==0 or c5.derivative(y2)!=0: raise RuntimeError("P1_5 not affine-linear in y2")
y2_expr=RF(-E5.subs({y2:0})/c5)
A3=subpoly(A2,{y2:y2_expr}); B3=subpoly(B2,{y2:y2_expr})
X3=subpoly(X2,{y2:y2_expr}); Y3=subpoly(Y2,{y2:y2_expr})
S3=Y3**2-X3**3-A3*X3-B3

E4=RF(S3[4]); c4=E4.derivative(a3)
if c4==0 or c4.derivative(a3)!=0: raise RuntimeError("P1_4 not affine-linear in a3")
a3_expr=RF(-E4.subs({a3:0})/c4)
A4=subpoly(A3,{a3:a3_expr}); B4f=subpoly(B3,{a3:a3_expr})
X4=subpoly(X3,{a3:a3_expr}); Y4=subpoly(Y3,{a3:a3_expr})
S4=Y4**2-X4**3-A4*X4-B4f

Ap=A4.derivative(t); Bp=B4f.derivative(t)
Delta=-16*(4*A4**3+27*B4f**2)
D1=Delta.derivative(t); D2=D1.derivative(t); D3=D2.derivative(t)

remaining=[n for n in names if n not in ("b4","b5","y2","a3")]
RR=PolynomialRing(K,remaining,order="degrevlex"); rrd=RR.gens_dict()
all_subs={b4:b4_expr,b5:b5_expr,y2:y2_expr,a3:a3_expr}

def settle(fr):
    fr=RF(fr)
    for _ in range(8):
        old=fr; fr=RF(fr.subs(all_subs))
        if fr==old: break
    return fr

def convert(fr):
    fr=settle(fr); num=R(fr.numerator()); out=RR(0)
    for mon,coef in num.dict().items():
        term=RR(coef)
        for i,exp in enumerate(mon):
            if exp:
                nm=names[i]
                if nm not in rrd: raise RuntimeError(f"eliminated variable survived: {nm}")
                term*=rrd[nm]**exp
        out+=term
    return out

eqs=[]; tags=[]
def add(tag,e):
    poly=convert(e)
    if poly!=0: tags.append(tag); eqs.append(poly)

add("I4_1_D2",D2(1)); add("I4_1_D3",D3(1))
for pref,a,s in [("I2_lam",settle(lam),settle(sl)),("I2_mu",settle(mu),settle(sm))]:
    add(pref+"_A",A4(a)+3*s**2)
    add(pref+"_B",B4f(a)-2*s**3)
    add(pref+"_d1",Bp(a)+s*Ap(a))
for k in range(4): add(f"P1_{k}",S4[k])
base_eq_count=len(eqs)

for j in range(args.slices):
    form=RR(K.random_element())
    for v in RR.gens():
        c=K.random_element()
        if c: form += c*v
    eqs.append(form); tags.append(f"SLICE_{j+1}")

print(f"E6SLICE|stage=build|p={args.p}|seed={args.seed}|vars={RR.ngens()}|base_eqs={base_eq_count}|slices={args.slices}|eqs={len(eqs)}",flush=True)
print("E6SLICE|history=P1_7->b4,P1_6->b5,P1_5->y2,P1_4->a3",flush=True)
for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6SLICE_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

meta=Path(args.out).with_suffix(".meta.txt"); meta.parent.mkdir(parents=True,exist_ok=True)
meta.write_text("remaining="+repr(remaining)+"\n"+"b4="+str(b4_expr)+"\n"+"b5="+str(b5_expr)+"\n"+"y2="+str(y2_expr)+"\n"+"a3="+str(a3_expr)+"\n")
out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write(",".join(remaining)+"\n"+str(args.p)+"\n")
    for i,e in enumerate(eqs):
        h.write(str(e).replace("**","^")); h.write(",\n" if i+1<len(eqs) else "\n")
print(f"E6SLICE|stage=export|out={out}|meta={meta}",flush=True)
