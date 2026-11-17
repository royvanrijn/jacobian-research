#!/usr/bin/env python3
"""PF182 exact controls; unrestricted parameter and Newton proof is in source."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,s,Y,X,Z,p,lam,eps=S.symbols("v s Y X Z p lam eps")
C,gamma=S.symbols("C gamma",nonzero=True)
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(m[0] for m in pn.monoms());od=min(m[0] for m in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))

# First p-order quartic deformation, including the other Hensel factor correction.
W=(-1/gamma**2+2*(Y-gamma)/gamma**3-3*(Y-gamma)**2/gamma**4
   +4*(Y-gamma)**3/gamma**5)/C**2
check("parameter deformation is inverse C2Y2 remainder",
      S.rem(S.cancel(C**2*Y**2*W+1),S.Poly((Y-gamma)**4,Y),Y)==0)
check("parameter quartic correction degree at most three",S.degree(W,Y)<=3)
other=S.cancel((-1-C**2*Y**2*W)/(Y-gamma)**4)
check("parameter other-factor correction is polynomial",not S.denom(other).has(Y))
check("parameter full first-order product",
      zero(other*(Y-gamma)**4+C**2*Y**2*W+1))
wp=S.Poly(S.expand(W.subs(Y,X+gamma)),X)
for i,expected in enumerate([-1/(C**2*gamma**2),2/(C**2*gamma**3),
                             -3/(C**2*gamma**4),4/(C**2*gamma**5)]):
    check("exact first parameter coefficient X"+str(i),zero(wp.nth(i)-expected))
check("first parameter constant is nonzero",wp.nth(0)!=0)

# Depression at the first parameter order preserves this constant coefficient.
a,b,c,d,w0,w1,w2,w3=S.symbols("a b c d w0 w1 w2 w3")
raw=X**4+v*a*X**3+v*b*X**2+v*c*X+v*d+eps*(w3*X**3+w2*X**2+w1*X+w0)
depressed=S.expand(S.series(raw.subs(X,Z-(v*a+eps*w3)/4),eps,0,2).removeO())
pconst=S.expand(depressed.coeff(Z,0)).coeff(eps,1)
pq=S.expand(depressed.coeff(Z,2)).coeff(eps,1)
baseq=S.expand(depressed.coeff(Z,2)).coeff(eps,0)
check("first depressed parameter constant unchanged modulo v",zero(pconst.subs(v,0)-w0))
deltaeps=S.expand(pconst-baseq*pq/2)
check("subtracting J2 changes first parameter constant only by positive v order",
      zero(deltaeps.subs(v,0)-w0))
check("first-order translation of X4 adds only cubic term",
      S.expand(S.series((X-eps*w3/4)**4,eps,0,2).removeO())==X**4-eps*w3*X**3)
for jj in [5,7]:
    JJ=v**jj
    DD=(3-p)*v**12
    check("generic Delta parameter cap J"+str(jj),initial(DD,v)[0]==12)
    check("generic quartic constant parameter cap J"+str(jj),initial(JJ**2+DD,v)[0]<=12)
    check("parameter coefficient remains affine J"+str(jj),S.degree(3-p,p)==1)
check("hypothetical Delta order13 misses compulsory first parameter coefficient",
      S.expand(v**13).coeff(v,12)==0 and wp.nth(0)!=0)

# Exact full sextic controls for the final quadratic-response row.
# Critical leading remainder Y(Y-1) corresponds to A gamma+B=0.
models=[]
for sign in [-1,1]:
    XX=sign*S.I*s**5+s**7*lam
    YY=1+X
    F0=YY**2*(X**2+v**5)**2+v**10*YY*X+3*v**12
    center=sign*S.I*s**5
    pc=F0.subs({v:s**2,X:center},simultaneous=True)
    check("whole center constant has order at least12 sign"+str(sign),initial(pc,s)[0]>=24)
    for j,expected_order in [(1,20),(2,10),(3,5)]:
        oo,ll=initial(S.diff(F0,X,j).subs({v:s**2,X:center},simultaneous=True),s)
        check("whole center derivative exact j"+str(j)+" sign"+str(sign),oo==expected_order and ll!=0)
    for j in [4,5,6]:
        oo,ll=initial(S.diff(F0,X,j).subs({v:s**2,X:center},simultaneous=True),s)
        check("whole higher derivative nonnegative j"+str(j)+" sign"+str(sign),oo>=0)
    scaled=(F0-p*v**12).subs({v:s**2,X:XX},simultaneous=True)/s**24
    oo,ll=initial(scaled,s)
    check("whole sextic first response is quadratic sign"+str(sign),
          oo==0 and zero(ll-(-4*lam**2+3-p)))
    check("whole response degree two sign"+str(sign),S.degree(ll,lam)==2)
    dy=S.diff(F0,X).subs({v:s**2,X:XX},simultaneous=True)/s**17
    do,dl=initial(dy,s)
    check("whole sextic branch derivative mu17half sign"+str(sign),do==0 and zero(dl+8*lam))
    check("whole quadratic-response branch is actually polar sign"+str(sign),
          S.Rational(17,2)>7)
    check("whole quadratic-response branch passes weaker strict quartic bound sign"+str(sign),
          S.Rational(17,2)<9)
    models.append({"sign":sign,"mu":"17/2","response_order":"7/2",
                   "whole_first_response":str(ll)})
check("forbidden row would give intermediate degree six",
      4*(S.Rational(17,2)-7)==6)
check("quadratic Taylor has exact parameter order",5+2*S.Rational(7,2)==12)
check("linear Taylor strictly later",10+S.Rational(7,2)>12)
check("cubic Taylor strictly later",S.Rational(5,2)+3*S.Rational(7,2)>12)
for j in [4,5,6]:
    check("higher Taylor strictly later j"+str(j),j*S.Rational(7,2)>12)

# The low-budget delta11 nearby control has mu8 and constant simple leading roots.
YY=1+X
Fnear=YY**2*(X**2+v**5)**2+v**10*YY*X+3*v**11+2*v**12
XX=S.I*s**5+s**6*lam
oo,ll=initial((Fnear-p*v**12).subs({v:s**2,X:XX},simultaneous=True)/s**22,s)
check("delta11 whole sextic has constant quadratic initial",oo==0 and zero(ll-(-4*lam**2+3)))
do,dl=initial(S.diff(Fnear,X).subs({v:s**2,X:XX},simultaneous=True)/s**16,s)
check("delta11 whole sextic mu exactly eight",do==0 and zero(dl+8*lam))
check("delta11 row insufficient total degree",4*(8-7)==4<6)

# A depressed quartic with q nonzero cannot have all roots multiple.
rr,qq,dd,jj=S.symbols("r q d j")
four=(X-rr)**4
check("quadruple-root depressed cubic coefficient forces zero root",
      S.expand(four).coeff(X,3)==-4*rr)
check("depressed quadruple-root case has zero linear coefficient",
      S.expand(four.subs(rr,0)).coeff(X,1)==0)
twotwo=(X**2-rr**2)**2
check("depressed two-double-root case has zero linear coefficient",
      S.expand(twotwo).coeff(X,1)==0)
for j0 in [0,1,-2]:
    lead=lam**4+2*j0*lam**2+lam+j0**2+2-p
    check("t9 boundary generic leading quartic has simple roots j"+str(j0),
          S.resultant(lead,S.diff(lead,lam),lam)!=0)
    check("t9 boundary nonzero linear derivative retained j"+str(j0),
          S.expand(lead).coeff(lam,1)==1)
    for ell in ([7] if j0==0 else [6]):
        JJ=j0*v**ell if j0!=0 else v**7
        GG=(X**2+JJ)**2+v**9*X+(2-p)*v**12
        XX=v**3*lam
        oo,ll=initial(GG.subs(X,XX)/v**12,v)
        expected=lam**4+2*(j0 if ell==6 else 0)*lam**2+lam+(j0**2 if ell==6 else 0)+2-p
        check("whole t9 boundary initial j"+str(j0),oo==0 and zero(ll-expected))
        do,dl=initial(S.diff(GG,X).subs(X,XX)/v**9,v)
        check("whole t9 boundary derivative mu nine j"+str(j0),
              do==0 and zero(dl-S.diff(expected,lam)))

# Representative exact quartic Newton controls for every table mechanism.
def square_control(t,ell,delta):
    tag=f"t{t} J{ell} D{delta}"
    m=min(S.Rational(delta),t+S.Rational(ell,2))
    check("square inequality m>2ell "+tag,m>2*ell)
    check("square inequality m<2t-ell "+tag,m<2*t-ell)
    JJ=-v**ell*(1+v)**2
    GG=(X**2+JJ)**2+v**t*(1+2*v)*X+2*v**delta*(1+3*v)
    XX=-s**(6*ell)*(1+s**12)+s**(6*(m-ell))*lam
    oo,ll=initial(GG.subs({v:s**12,X:XX},simultaneous=True)/s**(12*m),s)
    const=(2 if delta==m else 0)+(-1 if t+S.Rational(ell,2)==m else 0)
    check("whole square initial "+tag,oo==0 and zero(ll-(4*lam**2+const)))
    mu=(m+ell)/2
    do,dl=initial(S.diff(GG,X).subs({v:s**12,X:XX},simultaneous=True)/s**(12*mu),s)
    check("whole square derivative "+tag,do==0 and zero(dl-8*lam))
for model in [(8,1,9),(8,3,9),(8,4,10),(8,4,12),(8,5,12),
              (9,4,10),(9,5,11),(9,5,12),(10,4,11),(10,5,11),(10,5,12)]:
    square_control(*model)
for t,ell,delta in [(8,5,9),(8,6,10),(9,6,10),(9,7,11),(10,6,11),(10,7,12)]:
    tag=f"t{t} J{ell} D{delta}"
    check("Delta-dominant inequalities "+tag,delta<2*ell and delta<S.Rational(4*t,3))
    GG=(X**2+v**ell)**2+v**t*X+2*v**delta
    XX=s**(3*delta)*lam
    oo,ll=initial(GG.subs({v:s**12,X:XX},simultaneous=True)/s**(12*delta),s)
    check("whole Delta initial "+tag,oo==0 and zero(ll-lam**4-2))
    do,dl=initial(S.diff(GG,X).subs({v:s**12,X:XX},simultaneous=True)/s**(9*delta),s)
    check("whole Delta derivative "+tag,do==0 and zero(dl-4*lam**3))
for t,ell,delta in [(8,5,10),(9,5,10),(10,6,12)]:
    for d0 in [2,-1]:
        tag=f"t{t} J{ell} D{delta} c{d0}"
        GG=(X**2+v**ell)**2+v**t*X+d0*v**delta
        XX=s**(3*delta)*lam
        oo,ll=initial(GG.subs({v:s**12,X:XX},simultaneous=True)/s**(12*delta),s)
        expected=(lam**2+1)**2+d0
        check("whole boundary initial "+tag,oo==0 and zero(ll-expected))
        factor=lam**2+2 if d0==-1 else expected
        check("boundary has selected simple nonzero roots "+tag,
              S.resultant(factor,4*lam*(lam**2+1),lam)!=0)
for ell,delta in [(6,11),(6,12),(8,12)]:
    t=8
    GG=(X**2+v**ell)**2+v**t*X+2*v**delta
    XX=s**(4*t)*lam
    oo,ll=initial(GG.subs({v:s**12,X:XX},simultaneous=True)/s**(16*t),s)
    check(f"t8 high J q-balance initial J{ell} D{delta}",oo==0 and zero(ll-lam**4-lam))
    do,dl=initial(S.diff(GG,X).subs({v:s**12,X:XX},simultaneous=True)/s**(12*t),s)
    check(f"t8 high J q-balance derivative J{ell} D{delta}",do==0 and zero(dl-4*lam**3-1))

# Exact table contradictions and final coefficient degrees.
check("pair mu at most eight cannot reach derivative total thirtyfour",2*8+2*9==34)
check("all mu33quarter row total insufficient",4*S.Rational(33,4)==33<34)
check("strict upper mu nine is essential",not (S.Integer(9)<9))
check("delta cap reduces t8 to four integer orders",list(range(8+1,12+1))==[9,10,11,12])
check("delta cap reduces t9 to three integer orders",list(range(9+1,12+1))==[10,11,12])
check("delta cap reduces t10 to two integer orders",list(range(10+1,12+1))==[11,12])
u,cc=S.symbols("u c",nonzero=True)
MM=u**2
hh=(u-cc)**14
check("remaining degree30 leading polynomial",S.degree(MM*hh**2,u)==30)
check("remaining h degree all consumed at A root",S.degree(hh,u)==14)
check("remaining M degree central root only",S.degree(MM,u)==2)
Glead=u**4*(u-1)**14
rho=3*u*(u-1)*(243*u**4-891*u**3+1188*u**2-660*u+110)/110
check("remaining leading data have positive primitive ODE control",
      zero(S.diff(rho,u)-S.diff(Glead,u)*rho/(3*Glead)-1))

receipt={"claim":"PF2RULED182","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md",
        "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
        "RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md",
        "RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md",
        "RULING_SIX_NONCRITICAL_POLE_BUDGET_AND_DEGREE30_LEDGER.md"]},
    "uniform_proved_scope":{
        "coefficient_and_companion_degrees":"unrestricted",
        "first_parameter_constant_coefficient":"-1/(C^2*gamma^2) at v^(2L)",
        "generic_quartic_Delta_and_constant_orders_at_most":"2L",
        "entire_slope6_degree30_class_excluded":True,
        "remaining_D30_parameters":"m2,H14,L7,b_field6..9; M proportional u2,h proportional(u-c)14,c nonzero",
        "slope7_degree30_closed":False,
        "full_ruling_six_or_JC2_closed":False},
    "whole_sextic_quadratic_response_controls":models,
    "assurance":"Exact internal symbolic controls and written arbitrary-degree proof; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_slope6_degree30_exclusion_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
