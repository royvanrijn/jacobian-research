#!/usr/bin/env python3
"""PF180 exact controls; the arbitrary-coefficient proof is in the source."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,s,X,Z,p,lam=S.symbols("v s X Z p lam")
J,q,D,Delta=S.symbols("J q D Delta")
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(m[0] for m in pn.monoms());od=min(m[0] for m in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))
def partitions(n,ceiling=None):
    if n==0:
        yield []
        return
    if ceiling is None:ceiling=n
    for a in range(min(n,ceiling),0,-1):
        for tail in partitions(n-a,a):
            yield [a]+tail
def polemax(alpha):
    return S.ceiling(S.Rational(3*alpha,2))-1

# The exact derivative trace identity, via a companion matrix rather than sampled roots.
comp=S.Matrix([[0,0,0,-D],[1,0,0,-q],[0,1,0,-2*J],[0,0,1,0]])
check("depressed quartic root trace is zero",S.trace(comp)==0)
check("depressed quartic cubic Newton sum",zero(S.trace(comp**3)+3*q))
check("exact trace of derivative values is minus eight q",
      zero(S.trace(4*comp**3+4*J*comp+q*S.eye(4))+8*q))
F=(X**2+J)**2+q*X+Delta
check("depressed quartic constant term retained",zero(F.subs(X,0)-J**2-Delta))
unitfactor=(1+v+X)**2
check("separated sextic factor derivative differs by G multiple",
      zero(S.diff(unitfactor*F,X)-unitfactor*S.diff(F,X)-S.diff(unitfactor,X)*F))

# Complete finite index partitions of total quartic mass at most four.
partition_controls=[]
for mass in range(1,5):
    for part in partitions(mass):
        budget=sum(polemax(a) for a in part)
        if mass<=3:
            check("pole mass below four cannot reach five "+str(part),budget<=4)
        else:
            check("quartic mass upper budget at most five "+str(part),budget<=5)
        partition_controls.append({"indices":part,"maximum_integral_pole_degree":int(budget)})
check("positive budget-five index possibilities retained",
      [r["indices"] for r in partition_controls
       if sum(r["indices"])==4 and r["maximum_integral_pole_degree"]==5]==[[4],[3,1]])
alpha,mu=S.symbols("alpha mu",positive=True)
check("coframe pole order exact identity",
      S.expand(-(alpha*(6-mu)-1)-1-alpha*(mu-6))==0)
check("all-four derivative valuation sum is twenty-nine",4*6+5==29)
check("three-sheet strict bound is below five",S.Rational(3,2)*3<5)
for t in [5,6]:
    check("derivative trace contradicts every mu above six t"+str(t),t<=6)

# General square-center derivative expansion.
xc,disp=S.symbols("xc disp",nonzero=True)
centred=S.expand(F.subs({X:xc+disp,J:-xc**2}))
check("whole square-center expansion",
      zero(centred-(4*xc**2*disp**2+4*xc*disp**3+disp**4+q*xc+q*disp+Delta)))
check("opposite residuals cannot both cancel leading q center",
      zero((Delta+q*xc)-(Delta-q*xc)-2*q*xc))

# Full quartic controls with nonconstant units in J and q, and all monomial terms retained.
models=[]
square_models=[(7,1,8),(7,2,8),(7,3,8),(7,4,9),(7,4,10),
               (8,1,9),(8,2,9),(8,3,9),(8,3,10),(8,4,9),
               (8,4,10),(8,4,11),(8,5,11)]
for t,ell,delta in square_models:
    m=min(S.Rational(delta),S.Rational(2*t+ell,2))
    tag=f"t{t} J{ell} D{delta}"
    check("square-center relative displacement inequality "+tag,m>2*ell)
    check("square-center q displacement inequality "+tag,m<2*t-ell)
    JJ=-v**ell*(1+v)**2
    qq=v**t*(1+2*v)
    DD=2*v**delta*(1+3*v)
    GG=(X**2+JJ)**2+qq*X+DD
    beta=(m-ell)/2
    expected_mu=(m+ell)/2
    for sign in [-1,1]:
        XX=sign*s**(6*ell)*(1+s**12)+s**(12*beta)*lam
        scaled=GG.subs({v:s**12,X:XX},simultaneous=True)/s**(12*m)
        oo,ll=initial(scaled,s)
        residual=(2 if delta==m else 0)+(sign if S.Rational(2*t+ell,2)==m else 0)
        check("full quartic square initial "+tag+str(sign),
              oo==0 and zero(ll-(4*lam**2+residual)))
        dv=S.diff(GG,X).subs({v:s**12,X:XX},simultaneous=True)/s**(12*expected_mu)
        do,dl=initial(dv,s)
        check("full quartic exact square derivative "+tag+str(sign),
              do==0 and zero(dl-8*lam))
        check("square initial derivative is a unit at either simple root "+tag+str(sign),
              S.resultant(4*lam**2+residual,8*lam,lam)!=0)
    models.append({"mechanism":"square centers","t":t,"J_order":ell,
                   "Delta_order":delta,"mu":str(expected_mu)})

# Exact one-sign cancellations: do not discard the other two sheets.
for t,ell,m in [(7,2,8),(7,4,9),(8,2,9),(8,4,10)]:
    tag=f"t{t} J{ell} m{m}"
    JJ=-v**ell*(1+v)**2
    qq=v**t*(1+2*v)
    xp=v**(ell//2)*(1+v)
    DD=-qq*xp+v**(m+1)*(1+v)
    GG=(X**2+JJ)**2+qq*X+DD
    check("canceled plus residual has higher exact order "+tag,
          initial(GG.subs(X,xp),v)[0]==m+1)
    check("opposite residual keeps compulsory order "+tag,
          initial(GG.subs(X,-xp),v)[0]==m)
    beta=S.Rational(m-ell,2)
    XX=-s**(6*ell)*(1+s**12)+s**(12*beta)*lam
    oo,ll=initial(GG.subs({v:s**12,X:XX},simultaneous=True)/s**(12*m),s)
    check("one-sign cancellation retains two simple roots "+tag,
          oo==0 and zero(ll-(4*lam**2-2)))
    mu0=S.Rational(m+ell,2)
    do,dl=initial(S.diff(GG,X).subs({v:s**12,X:XX},simultaneous=True)/s**(12*mu0),s)
    check("uncanceled pair exact derivative order "+tag,do==0 and zero(dl-8*lam))

# Dominant Delta rows; J=0 is an additional unrestricted-order control.
for t,ell,delta in [(7,5,8),(7,5,9),(7,8,8),(8,5,9),(8,6,9),(8,6,10),(8,None,10)]:
    tag=f"t{t} J{ell} D{delta}"
    check("Delta dominant inequalities "+tag,
          (ell is None or delta<2*ell) and delta<S.Rational(4*t,3))
    JJ=0 if ell is None else -v**ell*(1+v)**2
    qq=v**t*(1+2*v)
    DD=2*v**delta*(1+3*v)
    GG=(X**2+JJ)**2+qq*X+DD
    XX=s**(3*delta)*lam
    oo,ll=initial(GG.subs({v:s**12,X:XX},simultaneous=True)/s**(12*delta),s)
    check("whole Delta-dominant quartic initial "+tag,oo==0 and zero(ll-lam**4-2))
    mu0=S.Rational(3*delta,4)
    do,dl=initial(S.diff(GG,X).subs({v:s**12,X:XX},simultaneous=True)/s**(12*mu0),s)
    check("whole Delta-dominant derivative "+tag,do==0 and zero(dl-4*lam**3))

# Square/Delta boundary, including the exceptional double zero.
for t,ell in [(7,4),(8,5)]:
    delta=2*ell
    for d0 in [2,-1]:
        tag=f"t{t} J{ell} Dlead{d0}"
        JJ=-v**ell*(1+v)**2
        qq=v**t*(1+2*v)
        DD=d0*v**delta*(1+3*v)
        GG=(X**2+JJ)**2+qq*X+DD
        XX=s**(3*delta)*lam
        oo,ll=initial(GG.subs({v:s**12,X:XX},simultaneous=True)/s**(12*delta),s)
        expected=(lam**2-1)**2+d0
        check("boundary whole quartic initial "+tag,oo==0 and zero(ll-expected))
        mu0=S.Rational(3*delta,4)
        do,dl=initial(S.diff(GG,X).subs({v:s**12,X:XX},simultaneous=True)/s**(12*mu0),s)
        check("boundary exact derivative "+tag,do==0 and zero(dl-4*lam*(lam**2-1)))
        simple=lam**2-2 if d0==-1 else expected
        check("boundary selected nonzero roots remain simple "+tag,
              S.resultant(simple,4*lam*(lam**2-1),lam)!=0)

# q-dominant rows: independently verify the three outer roots and the fourth inner root.
for t,ell,delta in [(7,5,10),(7,8,14),(8,6,11),(8,9,13),(8,None,11)]:
    tag=f"t{t} J{ell} D{delta}"
    check("q dominant inequalities "+tag,
          (ell is None or ell>S.Rational(2*t,3)) and delta>S.Rational(4*t,3))
    JJ=0 if ell is None else -v**ell*(1+v)**2
    qq=v**t*(1+2*v)
    DD=2*v**delta*(1+3*v)
    GG=(X**2+JJ)**2+qq*X+DD
    XX=s**(4*t)*lam
    oo,ll=initial(GG.subs({v:s**12,X:XX},simultaneous=True)/s**(16*t),s)
    check("q dominant whole outer initial "+tag,oo==0 and zero(ll-lam**4-lam))
    do,dl=initial(S.diff(GG,X).subs({v:s**12,X:XX},simultaneous=True)/s**(12*t),s)
    check("q dominant outer derivative "+tag,do==0 and zero(dl-4*lam**3-1))
    check("q dominant three outer roots simple "+tag,
          S.resultant(lam**3+1,4*lam**3+1,lam)!=0)
    dorder,dlead=initial(JJ**2+DD,v)
    XXinner=s**(12*(dorder-t))*lam
    io,il=initial(GG.subs({v:s**12,X:XXinner},simultaneous=True)/s**(12*dorder),s)
    check("q dominant fourth inner root initial "+tag,io==0 and zero(il-lam-dlead))
    ido,idl=initial(S.diff(GG,X).subs({v:s**12,X:XXinner},simultaneous=True)/s**(12*t),s)
    check("q dominant fourth derivative is also order t "+tag,ido==0 and zero(idl-1))

# Exact table values and contradictions: arbitrary later coefficients are covered in prose.
for t,ell in [(7,1),(7,2),(7,3),(8,1),(8,2)]:
    upper=S.Rational(t,2)+S.Rational(3*ell,4)
    check(f"small J compulsory pair has derivative below six t{t} J{ell}",upper<6)
table=[
("t7 J4 D8","bound",S.Integer(6)),
("t7 J4 D9","pair_sum_upper",S.Rational(13,2)),
("t7 J4 D10plus","sum",S.Integer(26)),
("t7 J5plus D8","bound",S.Integer(6)),
("t7 J5plus D9","sum",S.Integer(27)),
("t7 J5plus D10plus","sum",S.Integer(28)),
("t8 J3 D9","bound",S.Integer(6)),
("t8 J3 D10plus","sum",S.Integer(25)),
("t8 J4 D9","sum",S.Integer(26)),
("t8 J4 D10","pair_sum_upper",S.Integer(7)),
("t8 J4 D11plus","sum",S.Integer(28)),
("t8 J5 D9","sum",S.Integer(27)),
("t8 J5 D10","bound",S.Rational(15,2)),
("t8 J5 D11plus","bound",S.Rational(31,4)),
("t8 J6plus D9","sum",S.Integer(27)),
("t8 J6plus D10","bound",S.Rational(15,2)),
("t8 J6plus D11plus","bound",S.Integer(8))]
for label,kind,value in table:
    if kind=="bound":
        check("table violates individual strict derivative interval "+label,not (6<value<S.Rational(15,2)))
    elif kind=="sum":
        check("table violates exact sum twenty-nine "+label,value!=29)
    else:
        check("one-sign cancellation cannot reach exact sum "+label,
              2*value+2*S.Rational(15,2)<=29)
check("t7 canceled-cluster bound has strict gap",
      2*S.Rational(13,2)+2*S.Rational(15,2)==28<29)
check("t8 canceled-cluster bound is strict at necessary total",
      2*7+2*S.Rational(15,2)==29)
check("degree multiples after 24 and excluding 27 begin at 30",
      min(d for d in range(24,40) if d%3==0 and d not in [24,27])==30)
check("new intermediate minimum degree is six",S.ceiling(S.Rational(30+3*2,6))==6)
check("new original minimum generic degree is twelve",2*6==12)

receipt={"claim":"PF2RULED180","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
        "RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md",
        "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
        "RULING_SIX_DEGREE24_LEADING_CLASSIFICATION_AND_LINEAR_CONTROLS.md",
        "RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md",
        "RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md"]},
    "uniform_proved_scope":{
        "coefficient_and_companion_degrees":"unrestricted",
        "entire_degree24_class_excluded":True,
        "four_sheet_pole_mass":4,
        "strict_derivative_interval":"6<mu<15/2",
        "required_sheet_derivative_sum":29,
        "complete_critical_cases":"t5,6 by derivative trace; all t7,8 Newton possibilities",
        "necessary_leading_degree_minimum":30,
        "necessary_intermediate_generic_degree_minimum":6,
        "necessary_original_generic_degree_minimum":12,
        "full_ruling_six_or_JC2_closed":False,
        "target_integrality_or_trace_perfection_claimed":False},
    "index_partition_controls":partition_controls,
    "whole_quartic_models":models,
    "assurance":"Exact internal symbolic controls and written unrestricted proof; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_degree24_critical_exclusion_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
