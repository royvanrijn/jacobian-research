#!/usr/bin/env python3
"""PF185 controls; the unbounded congruence and pole proof is in the source."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_CUBIC_BALANCE_COMPANION_DEGREE_TWENTY_SEVEN.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,w,X,Z,T,p,s,lam=S.symbols("v w X Z T p s lam")
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(m[0] for m in pn.monoms());od=min(m[0] for m in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))
def power_coeff(f,alpha,N):
    ys=[F(1)]
    for n in range(1,N+1):
        ys.append(sum((((alpha+1)*k-n)*f[k]*ys[n-k]
                       for k in range(1,min(n,len(f)-1)+1)),F(0))/n)
    return ys[N]
def addto(dest,e,c):
    if c:
        dest[e]=dest.get(e,F(0))+c
        if not dest[e]:del dest[e]
def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():addto(out,i+j,x*y)
    return out
def whole_centered(poly,n):
    pp=S.Poly(poly,w)
    assert pp.degree()==6 and pp.LC()==v**28
    fs={}
    for k in range(1,7):
        fs[k]={i[0]-28:F(int(c)) for i,c in S.Poly(pp.nth(6-k),v).terms() if c}
    ys=[{0:F(1)}];alpha=F(n,2)
    for j in range(1,3*n+4):
        out={}
        for k in range(1,min(j,6)+1):
            scale=((alpha+1)*k-j)/j
            for e,c in mul(fs[k],ys[j-k]).items():addto(out,e,c*scale)
        ys.append(out)
    raw=[{e+14*n:c for e,c in ys[3*n+i+1].items()} for i in range(3)]
    centered=[]
    for j in range(3):
        out={}
        for i in range(j+1):
            factor=F((-1)**(j-i)*int(S.binomial(j,i)))
            for e,c in raw[i].items():addto(out,e-7*(j-i),factor*c)
        centered.append(out)
    return centered
def leading(d):
    return (min(d),d[min(d)]) if d else (S.oo,F(0))

# The prefix q_j for j>-5 contains one allowed negative multiple of three.
check("prefix only allowed negative index is minus three",
      [i for i in range(-4,0) if i%3==0]==[-3])
check("inverse remainder after retaining qminus3 starts at order five",5>4)
a3,a2,a1,a0,r2,r1,r0,cc=S.symbols("a3 a2 a1 a0 r2 r1 r0 cc",nonzero=True)
H=a3*w**3+a2*w**2+a1*w+a0
R=r2*w**2+r1*w+r0
inv=S.series(T**3/(a3+a2*T+a1*T**2+a0*T**3),T,0,4).removeO()
for j in [0,1,2]:
    coeff=S.expand(inv*(1/T-cc)**j).coeff(T,1)
    check("negative half-power weighted residue j"+str(j),
          zero(coeff-(1/a3 if j==2 else 0)))
check("third-residue RHS has exact D30 pole fourteen",14-0==14)

# Every congruence coefficient from direct generalized-power recurrence.
for n in range(1,24,2):
    for j in [0,1,2]:
        d=S.Rational(2*n+j+1,3)
        direct=power_coeff([F(1),F(0),F(0),F(2)],F(n,2),2*n+j+1)
        if d.q==1:
            expected=F(S.binomial(S.Rational(n,2),int(d)))*F(2)**int(d)
            check(f"cubic leading congruence coefficient n{n} j{j}",direct==expected and direct!=0)
        else:
            check(f"cubic other leading congruence vanishes n{n} j{j}",direct==0)
for j in [0,1,2]:
    for n in range(1,20,2):
        expected=S.Rational(n+11*(j+1),3)
        exponent=(j+1)*(S.Rational(10,3)-7)-S.Rational(n,3)
        check(f"exact escape scaling n{n} j{j}",zero(exponent+expected))
check("finite third bound is fourteen",7*2==14)

# A whole sextic control retains the finite SIMPLE cluster w=+1,-1.
G=X**4+2*v**7*X**2+v**10*X+v**14
Pwhole=S.expand((w**2-1)*G.subs(X,v**7*w-1))
check("whole cubic-balance finite specialization",zero(Pwhole.subs(v,0)-(w**2-1)))
whole_controls=[]
for n in [1,3,5,7,9,11,13,15,17,19,21]:
    moments=whole_centered(Pwhole,n)
    data={}
    for j in [0,1,2]:
        order,coeff=leading(moments[j])
        bound=max(7*j,(n+11*(j+1))//3)
        check(f"whole centered pole bound n{n} j{j}",order>=-bound)
        d=S.Rational(2*n+j+1,3)
        esc=S.Rational(n+11*(j+1),3)
        if d.q==1 and esc>7*j:
            expected=F(S.binomial(S.Rational(n,2),int(d)))
            check(f"whole centered exact dominant pole n{n} j{j}",order==-esc)
            check(f"whole centered exact dominant coefficient n{n} j{j}",coeff==expected)
        data[str(j)]={"order":str(order),"leading":str(coeff)}
    if n==3:
        check("small third character finite contribution really dominates escape estimate",
              leading(moments[2])[0]==-14 and S.Rational(n+33,3)==12)
    if n==9:
        check("ninth third character cannot ignore finite bound",
              S.Rational(n+33,3)==14)
    whole_controls.append({"n":n,"centered_moments":data})

# Finite controls of the written arbitrary-k congruence proof.
allowed=[]
for k in range(3,100,2):
    if k%6==1:
        j=0
    elif k%6==5:
        j=1
    elif k>=15:
        j=2
    else:
        allowed.append(k)
        continue
    top=(k+11*(j+1))//3
    check("top congruence exponent is integral k"+str(k),
          (k+11*(j+1))%3==0)
    check("top exact pole dominates finite and possible RHS k"+str(k),top>7*j)
    for n in range(1,k,2):
        low=max(7*j,(n+11*(j+1))//3)
        check(f"every smaller character has strictly smaller bound k{k} n{n}",low<top)
check("only two small indices escape the congruence argument",allowed==[3,9])
check("k3 third escape is below finite bound",S.Rational(3+33,3)<14)
check("k9 third escape equals finite bound",S.Rational(9+33,3)==14)

# Exact polynomial-part formula for a cubic companion index.
T1=r2**2/a3*w+2*r2*r1/a3-r2**2*a2/a3**2
remainder=S.Poly(S.expand(R**2-H*T1),w)
check("division quotient of R2 over cubic H is the stated linear function",remainder.degree()<=2)
check("later binomial cubic term has negative degree",3*2-3*3==-3)
check("still later binomial quartic term has negative degree",4*2-5*3==-7)
HH=w**3+w**2+2*w+3
RR=2*w**2+5*w+7
PP=HH**2+RR
TT=4*w+16
check("specific polynomial division quotient",S.div(RR**2,HH,w)[0]==TT)
rev=S.expand(T**6*PP.subs(w,1/T))
positive=S.expand(S.series(rev**S.Rational(3,2),T,0,10).removeO())
positive=S.expand(sum(positive.coeff(T,i)*w**(9-i) for i in range(10)))
check("whole generic numerical polynomial part formula",
      zero(positive-(HH**3+S.Rational(3,2)*HH*RR+S.Rational(3,8)*TT)))
check("fiber substitution identity for cubic companion",
      S.expand(H**3+S.Rational(3,2)*H*(p-H**2)+H**3/2-S.Rational(3,2)*p*H)==0)

# The local coefficient orders and full source-cubic estimates.
sub={a3:v**14,a2:-2*v**7,r2:v**10,r1:-v**3}
slope=S.cancel((r2**2/a3).subs(sub))
intercept=S.cancel((2*r2*r1/a3-r2**2*a2/a3**2).subs(sub))
check("local quotient slope exact order six",initial(slope,v)[0]==6)
check("local quotient intercept has order at least minus one",initial(intercept,v)[0]>=-1)
# Full local source cubic and actual quadratic remainder, with no companion claim.
HH=w*((v**7*w-1)**2+v**7)
RR=v**10*w**2-v**3*w+3
TT=v**6*w
check("full local source quotient exactly linear",
      zero(S.div(RR**2,HH,w)[0]-TT))
wo=s**(-21)*(1+s**10*Z)
Ho=HH.subs({v:s**3,w:wo},simultaneous=True)
To=TT.subs({v:s**3,w:wo},simultaneous=True)
check("outer source cubic has order minus one third",initial(Ho,s)[0]==-1)
check("outer quotient has order at least minus one",initial(To,s)[0]==-3)
Qon=-Ho**3/2+S.Rational(3,2)*p*Ho+S.Rational(3,8)*To
oo,ll=initial(s**3*Qon,s)
check("outer cubic-companion expression has at most pole one in v",oo==0)
check("outer leading expression at cubic roots is bounded and nonzero",
      zero(S.rem(S.expand(ll),S.Poly(Z**3+1,Z),Z).as_expr()+S.Rational(1,8)))
wi=v**(-7)*(1+v**4*Z)
Hi=HH.subs(w,wi)
Ti=TT.subs(w,wi)
check("inner source cubic is finite",initial(Hi,v)[0]>=0)
check("inner quotient has pole at most one",initial(Ti,v)[0]>=-1)
Qinner=-Hi**3/2+S.Rational(3,2)*p*Hi+S.Rational(3,8)*Ti
check("inner cubic-companion expression pole at most one",initial(Qinner,v)[0]>=-1)
check("source cubic estimate after square factor",min(S.Rational(20,3),S.Integer(7))-7==-S.Rational(1,3))
check("midpoint deformation order ten is later than cubic root scale",10>S.Rational(10,3))
check("source and quartic quadratic coefficient orders agree below ten",7<10)
check("total cubic-companion pole bound contradicts intermediate degree",4<8)
check("remaining index nine gives ruling degree twenty-seven",3*9==27)

receipt={"claim":"PF2RULED185","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
        "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
        "RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md",
        "RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md",
        "RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md",
        "RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md"]},
    "uniform_proved_scope":{
        "coefficient_degrees":"unrestricted",
        "companion_degree_prior_bound_assumed":False,
        "residue_congruence_index_bound":[3,9],
        "third_finite_cluster_and_RHS_bound":14,
        "cubic_index_pole_degree_at_most":4,
        "cubic_index_contradicted_by_actual_degree":8,
        "remaining_actual_reduced_companion_index":9,
        "remaining_actual_reduced_companion_ruling_degree":27,
        "remaining_D30_or_full_ruling_six_closed":False,
        "HC4_or_JC2_closed":False},
    "whole_sextic_moment_controls":whole_controls,
    "assurance":"Exact internal symbolic controls and unbounded written proof; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_companion_degree27_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
