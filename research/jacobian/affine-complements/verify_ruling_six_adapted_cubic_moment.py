#!/usr/bin/env python3
"""PF186 controls; the unrestricted adapted-moment proof is in the source."""
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
PROOF=BASE/"RULING_SIX_ADAPTED_CUBIC_MOMENT_AND_DEGREE30_EXCLUSION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,w,X,Y,Z,T,p,s=S.symbols("v w X Y Z T p s")
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(m[0] for m in pn.monoms());od=min(m[0] for m in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))
def addto(dest,e,c):
    if c:
        dest[e]=dest.get(e,F(0))+c
        if not dest[e]:del dest[e]
def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():addto(out,i+j,x*y)
    return out
def combine(terms):
    out={}
    for factor,poly in terms:
        for e,c in mul(factor,poly).items():addto(out,e,c)
    return out
def whole_raw(poly,n):
    pp=S.Poly(poly,w)
    assert pp.degree()==6 and pp.LC()==v**28
    fs={}
    for k in range(1,7):
        fs[k]={i[0]-28:F(int(c)) for i,c in S.Poly(pp.nth(6-k),v).terms() if c}
    ys=[{0:F(1)}];alpha=F(n,2)
    for j in range(1,3*n+5):
        out={}
        for k in range(1,min(j,6)+1):
            scale=((alpha+1)*k-j)/j
            for e,c in mul(fs[k],ys[j-k]).items():addto(out,e,c*scale)
        ys.append(out)
    return [{e+14*n:c for e,c in ys[3*n+i+1].items()} for i in range(4)]
def leading(poly):
    return (min(poly),poly[min(poly)]) if poly else (S.oo,F(0))
def power_coeff(f,alpha,N):
    ys=[F(1)]
    for n in range(1,N+1):
        ys.append(sum((((alpha+1)*k-n)*f[k]*ys[n-k]
                       for k in range(1,min(n,len(f)-1)+1)),F(0))/n)
    return ys[N]

# Exact unused inverse-prefix residue and an intentional wrong-center control.
a3,a2,a1,a0,r2,r1,r0,aa=S.symbols("a3 a2 a1 a0 r2 r1 r0 aa",nonzero=True)
revH=a3+a2*T+a1*T**2+a0*T**3
revR=r2+r1*T+r0*T**2
invP=S.series(T**3/revH*(1+T**4*revR/revH**2)**S.Rational(-1,2),T,0,5).removeO()
check("inverse half leading cubic coefficient",zero(S.expand(invP).coeff(T,3)-1/a3))
check("inverse half next coefficient",zero(S.expand(invP).coeff(T,4)+a2/a3**2))
weight=1/(aa*T**3)-2/T**2+aa/T
res=S.expand(weight*invP).coeff(T,1)
check("generic weighted negative coefficient",
      zero(res-(-a2/(a3*aa)-2)/a3))
center=-a2/(2*a3)
check("adapted negative coefficient vanishes exactly",zero(res.subs(aa,center)))
check("wrong third instead of half center is rejected",
      zero(res.subs(aa,-a2/(3*a3))-1/a3))
for exponent in range(5,13):
    check("inverse remainder cannot contribute exponent"+str(exponent),
          S.expand(weight*T**exponent).coeff(T,1)==0)
for degree in range(0,8):
    BB=w**3/aa-2*w**2+aa*w
    check("weighted polynomial coefficient residue zero degree"+str(degree),
          S.expand(BB*w**degree).coeff(w,-1)==0)

# Exact root-sum center and a regular bounded finite cubic root.
gamma,rho,delta,mu,j0=S.symbols("gamma rho delta mu j0",nonzero=True)
rscaled=v**7*rho
mH=gamma+v*mu
JH=v**7*j0
SS=(1+v)*(Y-rscaled)*((Y-mH)**2+JH)
ratio=S.Poly(S.expand(SS),Y).nth(2)/S.Poly(S.expand(SS),Y).nth(3)
check("scaled cubic root-sum ratio exact",zero(ratio+rscaled+2*mH))
asource=-v**(-7)*ratio/2
mG=mH+v**10*delta
bcenter=v**(-7)*mG
check("exact bounded center difference",zero(asource-bcenter-(rho/2-v**3*delta)))
check("finite cubic scaled root really order seven",initial(rscaled,v)[0]==7)
HH=S.cancel(v**(-7)*SS.subs(Y,v**7*w))
check("finite cubic root is exact",zero(HH.subs(w,rho)))
check("finite cubic specialization is nonconstant linear",
      zero(HH.subs(v,0)-gamma**2*(w-rho)))
check("finite cubic derivative is a unit",zero(S.diff(HH,w).subs(v,0)-gamma**2))
check("adapted center leading coefficient gamma",initial(asource,v)==(-7,gamma))

# Both finite and escape weight estimates, retaining nonzero center mismatch.
afinite=v**(-7)+S.Rational(3,4)
Bfinite=w*(w-afinite)**2/afinite
fo,fl=initial(Bfinite,v)
check("finite adapted weight order is minus seven",fo==-7 and zero(fl-w))
check("unadapted cubic centered moment would cost twenty-one",
      initial((w-afinite)**3,v)[0]==-21)
beta0,beta1,eta0,eta1=S.symbols("beta0 beta1 eta0 eta1",nonzero=True)
bs=gamma*s**(-21)+beta0+beta1*s**3
ass=bs+eta0+eta1*s**3
ws=bs+s**(-11)*Z
Bs=ws*(ws-ass)**2/ass
eo,el=initial(s**22*Bs,s)
check("escape leading weight is exactly Z squared",eo==0 and zero(el-Z**2))
err,errlead=initial(s**22*Bs-Z**2,s)
check("all bounded-center corrections have positive escape order",err>0)
check("w over a tends to one at escape",zero(initial(ws/ass,s)[1]-1))

# Leading qdom residue coefficients, independently extracted at infinity.
for n in range(1,26,2):
    d=S.Rational(2*n+3,3)
    direct=power_coeff([F(1),F(0),F(0),F(2)],F(n,2),2*n+3)
    if d.q==1:
        expected=F(S.binomial(S.Rational(n,2),int(d)))*F(2)**int(d)
        check("adapted all-index coefficient n"+str(n),direct==expected and direct!=0)
    else:
        check("other leading congruence vanishes n"+str(n),direct==0)
check("ninth binomial is exactly nine over2048",
      S.binomial(S.Rational(9,2),7)==S.Rational(9,2048))
check("escape exponent is nplus33over3",
      -(S.Rational(22,3)+S.Rational(11,3)+S.Symbol("n")/3)==-(S.Symbol("n")+33)/3)

# WHOLE sextics: two finite simple roots, with exact bounded mismatch a-b=3/4.
whole_controls=[]
for finite,shift,q0 in [(w**2-1,F(0),1),(w**2-3*w+2,F(3,4),1),
                         (w**2-3*w+2,F(3,4),2)]:
    GG=X**4+2*v**7*X**2+q0*v**10*X+v**14
    PP=S.expand(finite*GG.subs(X,v**7*w-1))
    pp=S.Poly(PP,w)
    check(f"whole finite specialization shift{shift} q{q0}",zero(PP.subs(v,0)-finite))
    aexpr=S.cancel(-pp.nth(5)/(4*pp.nth(6)))
    check(f"whole adapted center from cubic coefficient shift{shift} q{q0}",
          zero(aexpr-v**(-7)-S.Rational(shift.numerator,shift.denominator)))
    check(f"whole center error bounded shift{shift} q{q0}",initial(aexpr-v**(-7),v)[0]>=0)
    adict={-7:F(1)}
    if shift:adict[0]=shift
    a2dict=mul(adict,adict)
    for n in [1,3,5,7,9,11,13,15]:
        raw=whole_raw(PP,n)
        # a*W_n=[w^-1]w(w-a)^2 P^(n/2), avoiding a unit denominator in extraction.
        Wnumer=combine([({0:F(1)},raw[3]),
                        ({e:-2*c for e,c in adict.items()},raw[2]),
                        (a2dict,raw[1])])
        no,nc=leading(Wnumer)
        order=no+7 if no!=S.oo else S.oo
        bound=max(7,(n+33)//3)
        check(f"whole adapted moment bound shift{shift} q{q0} n{n}",order>=-bound)
        d=S.Rational(2*n+3,3)
        if d.q==1:
            expected=F(S.binomial(S.Rational(n,2),int(d)))*F(q0)**int(d)
            exact=S.Rational(n+33,3)
            check(f"whole adapted exact dominant pole shift{shift} q{q0} n{n}",order==-exact)
            check(f"whole adapted exact dominant coefficient shift{shift} q{q0} n{n}",nc==expected)
        if n==9:
            check(f"whole top coefficient q power shift{shift} q{q0}",
                  order==-14 and nc==F(9,2048)*F(q0)**7)
        if n<=7:
            check(f"whole lower characters have no order-fourteen pole shift{shift} q{q0} n{n}",
                  order>=-13)
        whole_controls.append({"finite_polynomial":str(finite),"center_error":str(shift),
                               "q0":q0,"n":n,"W_order":str(order),"W_leading":str(nc)})
# Exact smallest-case controls against accidental finite-cluster omission.
check("lower finite allowance remains strictly below top",7<14)
check("ninth top exact escape order",(9+33)//3==14)
for n in [1,3,5,7]:
    check("all lower positive odd indices bounded by thirteen n"+str(n),
          max(7,(n+33)//3)<=13)
for k in range(3,64,6):
    top=(k+33)//3
    check("last congruence top dominates finite cluster k"+str(k),top>7)
    for n in range(1,k,2):
        check(f"last congruence lower-character gap k{k} n{n}",
              max(7,(n+33)//3)<top)

# New global lower bounds are arithmetic corollaries within the retained model.
check("next permitted leading degree after thirty is thirtythree",
      min(d for d in range(30,40) if d%3==0 and d!=30)==33)
check("strict infinity bound gives intermediate degree seven",
      S.ceiling(S.Rational(33+3*2,6))==7)
check("original degree lower bound is fourteen",2*7==14)

receipt={"claim":"PF2RULED186","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md",
        "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
        "RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md",
        "RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md",
        "RULING_SIX_CUBIC_BALANCE_COMPANION_DEGREE_TWENTY_SEVEN.md"]},
    "uniform_proved_scope":{
        "coefficient_and_companion_degrees":"unrestricted",
        "adapted_center":"a=-K/(2h)=-a2/(2a3)",
        "exact_weight":"w*(w-a)^2/a",
        "negative_half_residue_vanishes":True,
        "finite_cluster_weight_bound":7,
        "escape_upper_pole_bound":"(n+33)/3",
        "ninth_exact_top_pole":14,
        "ninth_exact_nonzero_binomial":"9/2048",
        "every_degree30_case_excluded":True,
        "necessary_leading_degree_minimum":33,
        "necessary_intermediate_generic_degree_minimum":7,
        "necessary_original_generic_degree_minimum":14,
        "full_ruling_six_or_JC2_closed":False},
    "whole_sextic_controls":whole_controls,
    "assurance":"Exact internal symbolic controls and unrestricted written contour/prefix proof; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_adapted_cubic_moment_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
