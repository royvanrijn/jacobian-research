from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Iterative triangular eliminator for the E6 P1 construction.")
ap.add_argument("--p",type=int,default=0)
ap.add_argument("--max-steps",type=int,default=8)
ap.add_argument("--max-expr-terms",type=int,default=250)
ap.add_argument("--export",default=None)
args=ap.parse_args()

K=QQ if args.p==0 else GF(args.p)
if args.p and (not is_prime(args.p) or args.p in (2,3,79)):
    raise SystemExit("choose good prime")

# Base section-first system variables.
names=["a1","a2","a3","a4","b4","b5",
       "lam","mu","s0","s1","sl","sm",
       "x0","x1","y0","y1","y2"]

R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()
RF=FractionField(R)
Rt=PolynomialRing(RF,"t")
t=Rt.gen()

a1,a2,a3,a4,b4,b5=[RF(d[n]) for n in ["a1","a2","a3","a4","b4","b5"]]
lam,mu,s0,s1,sl,sm=[RF(d[n]) for n in ["lam","mu","s0","s1","sl","sm"]]
x0,x1,y0,y1,y2=[RF(d[n]) for n in ["x0","x1","y0","y1","y2"]]

# Fiber reductions at t=0 and t=1.
a0=-3*s0**2
b0=2*s0**3
b1=-s0*a1
b2=a1**2/(12*s0)-s0*a2
b3=(a1**3+36*a1*a2*s0**2-216*a3*s0**4)/(216*s0**3)
a5=-3*s1**2-(a0+a1+a2+a3+a4)

S=PolynomialRing(RF,["B6"])
B6=S.gen()
St=PolynomialRing(S,"t")
tt=St.gen()
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

x2=s1-x0-x1
y3=-1-y0-y1-y2
X=x0+x1*t+x2*t**2
Y=y0+y1*t+y2*t**2+y3*t**3+t**4

# Current substitution dictionary on base coefficient field.
subs={}

def sub_coeff(fr):
    fr=RF(fr)
    # Apply substitutions repeatedly to allow chained expressions.
    for _ in range(4):
        old=fr
        if subs:
            fr=RF(fr.subs(subs))
        if fr==old:
            break
    return fr

def subpoly(poly):
    return sum(sub_coeff(poly[k])*t**k for k in range(poly.degree()+1))

def current_objects():
    AA=subpoly(A)
    BB=subpoly(B)
    XX=subpoly(X)
    YY=subpoly(Y)
    Ap=AA.derivative(t)
    Bp=BB.derivative(t)
    Delta=-16*(4*AA**3+27*BB**2)
    D1=Delta.derivative(t); D2=D1.derivative(t); D3=D2.derivative(t)
    Sec=YY**2-XX**3-AA*XX-BB
    return AA,BB,XX,YY,Ap,Bp,D2,D3,Sec

def term_count(fr):
    fr=RF(fr)
    return len(R(fr.numerator()).monomials()) + len(R(fr.denominator()).monomials())

def live_names():
    eliminated={str(v) for v in subs.keys()}
    return [n for n in names if n not in eliminated]

def pick_elimination(Sec):
    # Search from highest remaining P1 coefficient downward.
    maxk=Sec.degree()
    candidates=[]
    priority=["b4","b5","y2","a3","a2","a4","x1","y1","a1","x0","y0"]
    rank={n:i for i,n in enumerate(priority)}
    for k in range(maxk,-1,-1):
        e=sub_coeff(Sec[k])
        if e==0:
            continue
        for nm in live_names():
            if nm not in d:
                continue
            v=RF(d[nm])
            coeff=e.derivative(v)
            if coeff==0:
                continue
            if coeff.derivative(v)!=0:
                continue
            const=e.subs({v:0})
            expr=RF(-const/coeff)
            score=(
                k*-1,  # prefer highest coefficient
                term_count(coeff),
                term_count(expr),
                rank.get(nm,99)
            )
            candidates.append((score,k,nm,v,expr,coeff,e))
        # If highest nonzero coefficient has any affine-linear candidate, use only it.
        if candidates:
            topk=max(c[1] for c in candidates)
            candidates=[c for c in candidates if c[1]==topk]
            break
    if not candidates:
        return None
    candidates.sort(key=lambda z:z[0])
    return candidates[0]

print(f"E6LOOP|stage=start|vars={len(names)}|max_steps={args.max_steps}|max_expr_terms={args.max_expr_terms}",flush=True)

history=[]
for step in range(1,args.max_steps+1):
    AA,BB,XX,YY,Ap,Bp,D2,D3,Sec=current_objects()
    choice=pick_elimination(Sec)
    if choice is None:
        print(f"E6LOOP|stage=stop|reason=no_affine_linear_section_variable|step={step}",flush=True)
        break
    score,k,nm,v,expr,coeff,e=choice
    terms=term_count(expr)
    print(f"E6LOOP|step={step}|eq=P1_{k}|variable={nm}|expr_terms={terms}|coeff_terms={term_count(coeff)}",flush=True)
    if terms>args.max_expr_terms:
        print(f"E6LOOP|stage=stop|reason=expression_threshold|variable={nm}|terms={terms}",flush=True)
        break
    subs[v]=expr
    history.append((k,nm,expr))
else:
    print(f"E6LOOP|stage=stop|reason=max_steps",flush=True)

# Build final equations after all triangular section eliminations.
AA,BB,XX,YY,Ap,Bp,D2,D3,Sec=current_objects()

remaining=live_names()
RR=PolynomialRing(K,remaining,order="degrevlex")
rrd=RR.gens_dict()

def convert(fr):
    fr=sub_coeff(fr)
    num=fr.numerator()
    den=fr.denominator()
    result_num=RR(0)
    result_den=RR(0)
    def cv(poly):
        out=RR(0)
        for mon,coef in R(poly).dict().items():
            term=RR(coef)
            for i,exp in enumerate(mon):
                if exp==0:
                    continue
                nm=names[i]
                if nm not in rrd:
                    raise RuntimeError(f"eliminated variable {nm} survived conversion")
                term*=rrd[nm]**exp
            out+=term
        return out
    return cv(num),cv(den)

eqs=[]; tags=[]
def add(tag,e):
    num,den=convert(e)
    if num!=0:
        tags.append(tag); eqs.append(num)

# Fiber equations.
add("I4_1_D2",D2(1))
add("I4_1_D3",D3(1))
for pref,a,s in [("I2_lam",sub_coeff(lam),sub_coeff(sl)),
                 ("I2_mu",sub_coeff(mu),sub_coeff(sm))]:
    add(pref+"_A",AA(a)+3*s**2)
    add(pref+"_B",BB(a)-2*s**3)
    add(pref+"_d1",Bp(a)+s*Ap(a))

# Remaining section coefficient equations.
for k in range(Sec.degree()+1):
    e=sub_coeff(Sec[k])
    if e!=0:
        add(f"P1_{k}",e)

field="QQ" if args.p==0 else f"GF({args.p})"
print(f"E6LOOP|stage=final|field={field}|vars={RR.ngens()}|eqs={len(eqs)}|naive_dim={RR.ngens()-len(eqs)}",flush=True)
print("E6LOOP|history="+",".join(f"P1_{k}->{nm}" for k,nm,_ in history),flush=True)

for i,(tag,e) in enumerate(zip(tags,eqs)):
    print(f"E6LOOP_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

if args.export:
    if args.p==0:
        raise SystemExit("--export requires --p")
    out=Path(args.export)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w") as h:
        h.write(",".join(remaining)+"\\n"+str(args.p)+"\\n")
        for i,e in enumerate(eqs):
            h.write(str(e).replace("**","^"))
            h.write(",\\n" if i+1<len(eqs) else "\\n")
    print(f"E6LOOP|export={out}",flush=True)
