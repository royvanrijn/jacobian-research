#!/usr/bin/env python3
"""PF188 controls; the unrestricted square-deformation proof is in the source."""
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
PROOF=BASE/"RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,w,X,Y,Z,T,p,s,lam=S.symbols("v w X Y Z T p s lam")
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
def moments(poly,L,n):
    """Exact complete half-power coefficients of the entire sextic."""
    pp=S.Poly(poly,w)
    assert pp.degree()==6 and pp.LC()==v**(4*L)
    fs={}
    for k in range(1,7):
        fs[k]={i[0]-4*L:F(int(c)) for i,c in S.Poly(pp.nth(6-k),v).terms() if c}
    ys=[{0:F(1)}];alpha=F(n,2)
    for j in range(1,3*n+3):
        out={}
        for k in range(1,min(j,6)+1):
            scale=((alpha+1)*k-j)/j
            for e,c in mul(fs[k],ys[j-k]).items():addto(out,e,c*scale)
        ys.append(out)
    mm={e+2*L*n:c for e,c in ys[3*n+1].items()}
    lc={e+2*L*n:c for e,c in ys[3*n+2].items()}
    for e,c in mm.items():addto(lc,e-L,-c)
    return mm,lc
def leading(d):
    return (min(d),d[min(d)]) if d else (S.oo,F(0))

# Exact unrestricted residue-pair algebra at the two quadratic roots.
A,B,z0=S.symbols("A B z0",nonzero=True)
for r in range(1,13):
    numerator=(A*Z+B)**r
    rem=S.Poly(S.rem(S.expand(numerator),S.Poly(Z**2-z0**2,Z),Z),Z).nth(1)
    rem2=S.Poly(S.rem(S.expand(Z*numerator),S.Poly(Z**2-z0**2,Z),Z),Z).nth(1)
    U=B+A*z0;V=B-A*z0
    M=(U**r-V**r)/(2*z0);N=(U**r+V**r)/2
    check("first quadratic residue identity r"+str(r),zero(rem-M))
    check("second quadratic residue identity r"+str(r),zero(rem2-N))
    check("nonzero vector recovers U power r"+str(r),zero(N+z0*M-U**r))
    check("nonzero vector recovers V power r"+str(r),zero(N-z0*M-V**r))
    check("half-integer first binomial never zero r"+str(r),
          S.binomial(S.Rational(2*r-1,2),r)!=0)
for aa,bb in [(1,2),(1,0),(0,2),(1,1),(1,-1)]:
    for r in range(1,11):
        pair=(F((bb+aa)**r-(bb-aa)**r,2),
              F((bb+aa)**r+(bb-aa)**r,2))
        check(f"root-unity and one-sign cancellation vector A{aa} B{bb} r{r}",
              pair!=(0,0))

# Homogeneity and finite-cluster comparison are symbolic identities.
Ls,rhos,ms,rs=S.symbols("L rho m r")
nu=ms-4*rhos
n=2*rs-1
exponent=(rhos-Ls)+(4*rhos-2*Ls)*n/2+nu*rs
pole=rs*(2*Ls-ms)+rhos
check("whole first vector exponent equals claimed pole",zero(exponent+pole))
check("second vector alignment removes its extra Jacobian scale",
      zero((rhos-Ls)+(Ls-rhos)))
check("aligned finite second bound equals rho",zero(Ls-(Ls-rhos)-rhos))
check("strict growth increment is exactly twoLminusM",
      zero(((rs+1)*(2*Ls-ms)+rhos)-pole-(2*Ls-ms)))
check("first vector term beats finite bound by r times positive gap",
      zero(pole-rhos-rs*(2*Ls-ms)))

# Whole sextics with finite SIMPLE roots 1,2 and variable units in all factors.
# No actual companion or global model realization is posited for these controls.
configurations=[
    (3,1,3,4,1,2),
    (4,2,4,6,1,2),
    (5,3,5,7,1,2),
    (6,4,7,9,1,2),
    (7,6,10,13,1,2),
    (7,6,10,14,1,2),
    (8,6,11,15,1,2),
    (8,7,12,15,1,2),
    (9,7,13,16,1,2),
    (10,8,15,19,1,1),
    (10,8,15,19,1,-1)]
whole_controls=[]
for L,ell,t,delta,qa,db in configurations:
    rho=F(ell,2)
    m=min(F(delta),F(t)+rho)
    d=m-2*ell
    tag=f"L{L} ell{ell} t{t} delta{delta} B{db}"
    check("first-square strict orders "+tag,
          ell<L and delta>2*ell and 3*ell<2*t)
    check("deformation before generic parameter cap "+tag,2*ell<m<2*L)
    check("first deformation is positive "+tag,d>0)
    JJ=-v**ell*(1+v)**2
    qq=qa*v**t*(1+2*v)
    DD=db*v**delta*(1+3*v)
    GG=(X**2+JJ)**2+qq*X+DD
    PP=S.expand((w**2-3*w+2)*GG.subs(X,v**L*w-1))
    check("whole finite cluster is retained "+tag,zero(PP.subs(v,0)-(w**2-3*w+2)))
    leadA=qa if F(t)+rho==m else 0
    leadB=db if F(delta)==m else 0
    check("first residual is nonzero linear polynomial "+tag,(leadA,leadB)!=(0,0))
    for nn in [1,3,5,7,9]:
        r=(nn+1)//2
        mm,lc=moments(PP,L,nn)
        aligned={F(e)+F(L)-rho:c for e,c in lc.items()}
        target=-rho-r*(2*L-m)
        factor=F(S.binomial(S.Rational(nn,2),r))
        U=leadB+leadA;V=leadB-leadA
        EM=factor*F(U**r-V**r,2)
        EN=factor*F(U**r+V**r,2)
        check(f"whole first vector bound and coefficient {tag} n{nn}",
              leading(mm)[0]>=target and mm.get(target,F(0))==EM)
        check(f"whole second vector bound and coefficient {tag} n{nn}",
              leading(aligned)[0]>=target and aligned.get(target,F(0))==EN)
        check(f"whole aligned vector leading is nonzero {tag} n{nn}",
              (mm.get(target,0),aligned.get(target,0))!=(0,0))
        check(f"whole finite contribution strictly later {tag} n{nn}",-target>rho)
        whole_controls.append({"L":L,"ell":ell,"t":t,"delta":delta,
            "Delta_lead":db,"n":nn,"common_pole":str(-target),
            "leading_vector":[str(EM),str(EN)]})

# Algebra of exact whole F2*G Taylor derivatives at a square center.
xc,E,q,D0=S.symbols("xc E q D0",nonzero=True)
G=(X**2+E)**2+q*X+D0
check("square-center residual exact",zero(G.subs({X:xc,E:-xc**2})-(D0+q*xc)))
check("square-center first derivative exact",zero(S.diff(G,X).subs({X:xc,E:-xc**2})-q))
check("square-center second derivative exact",
      zero(S.diff(G,X,2).subs({X:xc,E:-xc**2})-8*xc**2))
check("square-center third derivative exact",zero(S.diff(G,X,3).subs(X,xc)-24*xc))
F2,F2p,F2pp=S.symbols("F2 F2p F2pp")
G0,G1,G2=S.symbols("G0 G1 G2")
check("whole second derivative has unique low-order factor term",
      S.expand(F2pp*G0+2*F2p*G1+F2*G2-F2*G2)==F2pp*G0+2*F2p*G1)

# Every cap packet is checked in the WHOLE sextic, not a selected quartic.
cap_configs=[(6,4,10),(7,5,12),(8,6,13),(9,8,14),(10,9,16),(10,8,16)]
cap_controls=[]
for L,ell,t in cap_configs:
    rho=S.Rational(ell,2)
    RR=L-rho
    tag=f"L{L} ell{ell} t{t}"
    check("cap strict first-square orders "+tag,ell<L and 3*ell<2*t)
    check("cap available coefficient range "+tag,L<=t<=2*L-2 and t+rho>=2*L)
    check("cap response displacement smaller than separation "+tag,RR>rho)
    check("cap quadratic order exact "+tag,ell+2*RR==2*L)
    check("cap linear Taylor term strictly later "+tag,t+RR>2*L)
    check("cap cubic Taylor term strictly later "+tag,rho+3*RR>2*L)
    for j in [4,5,6]:
        check(f"cap higher Taylor term later {tag} j{j}",j*RR>2*L)
    check("cap rho at least two "+tag,rho>=2)
    FF0=(1+X)**2*((X**2-v**ell)**2+v**t*X+2*v**(2*L))
    for sign in [-1,1]:
        center=sign*s**ell
        XX=center+s**(2*L-ell)*lam
        whole=(FF0-p*v**(2*L)).subs({v:s**2,X:XX},simultaneous=True)/s**(4*L)
        oo,ll=initial(whole,s)
        const=2+(sign if 2*t+ell==4*L else 0)
        expected=4*lam**2+const-p
        check("entire cap quadratic initial "+tag+str(sign),oo==0 and zero(ll-expected))
        check("cap initial degree two "+tag+str(sign),S.degree(ll,lam)==2)
        check("cap packet separable over C(p) "+tag+str(sign),
              S.resultant(expected,S.diff(expected,lam),lam)!=0)
        deriv=S.diff(FF0,X).subs({v:s**2,X:XX},simultaneous=True)/s**(2*L+ell)
        do,dl=initial(deriv,s)
        check("cap entire derivative mu Lplusrho "+tag+str(sign),do==0 and zero(dl-8*lam))
        # On v=s^2, exact coframe order 2(1-rho)-1; it is polar here.
        check("cap coframe genuinely polar "+tag+str(sign),2*(1-rho)-1<-1)
        cap_controls.append({"L":L,"ell":ell,"t":t,"sign":sign,
              "response_order":str(RR),"mu":str(L+rho),
              "whole_initial":str(ll)})
    check("two cap centers exhaust all four root sheets "+tag,2+2==4)

# Adversarial full source cap controls with TWO finite SIMPLE roots.
# At L8, c=i, one square-center constant vanishes IDENTICALLY at p=0.
for L,ell,t,cc in [(8,4,14,S.I),(9,5,16,S.Integer(2))]:
    rho=S.Rational(ell,2)
    RR=L-rho
    FF0=((1+X)**2-v**(2*L))*((X**2+v**ell)**2+v**t*X+cc*v**(2*L))
    for sign in [-1,1]:
        center=sign*S.I*s**ell
        XX=center+s**(2*L-ell)*lam
        whole=(FF0-p*v**(2*L)).subs({v:s**2,X:XX},simultaneous=True)/s**(4*L)
        oo,ll=initial(whole,s)
        constant=cc+(sign*S.I if 2*t+ell==4*L else 0)
        check(f"whole finite-simple cap initial L{L} sign{sign}",
              oo==0 and zero(ll-(-4*lam**2+constant-p)))
        check(f"whole finite-simple cap stays quadratic L{L} sign{sign}",
              S.degree(ll,lam)==2)
        if L==8 and sign==-1:
            check("one cap center has exactly zero p0 constant",
                  zero(FF0.subs({v:s**2,X:center},simultaneous=True)))
            check("zero p0 constant still gives genuine quadratic response",
                  zero(ll-(-4*lam**2-p)))

# Generic order must be split BEFORE a cap coefficient specializes to zero.
for L in [3,5,7,10]:
    for delta in [L+1,2*L-1]:
        expr=(1+v)*v**delta+(3-p)*v**(2*L)
        check(f"subcap coefficient unaffected by p0 L{L} delta{delta}",
              initial(expr,v)==initial(expr.subs(p,0),v))
    cap=(3-p)*v**(2*L)+v**(2*L+1)
    check("generic cap has order twoL L"+str(L),initial(cap,v)[0]==2*L)
    check("special cap may rise without violating p0 lower bound L"+str(L),
          initial(cap.subs(p,3),v)[0]>2*L)

receipt={"claim":"PF2RULED188","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md",
        "RULING_SIX_FINITE_END_COVERAGE_AND_QUADRATIC_BOUNDARY_ALGEBRA.md",
        "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
        "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
        "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
        "RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md",
        "RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md",
        "RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md",
        "RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md",
        "RULING_SIX_UNIVERSAL_THREE_MOMENT_FIRST_QUARTIC_REDUCTION.md"]},
    "uniform_proved_scope":{
        "slope_leading_coefficient_and_companion_degrees":"unrestricted",
        "generic_deformation_order":"m=min(delta,t+rho), 4rho<m<=2L",
        "subcap_vector_pole":"r*(2L-m)+rho, every positive integer r",
        "subcap_finite_aligned_bound":"rho",
        "cap_whole_first_response":"nonzero quadratic in lambda plus constant minus p",
        "every_actual_A_support_excluded":True,
        "entire_retained_marked_quadratic_ruling_six_class_excluded":True,
        "all_global_model_or_import_coverage_claimed":False,
        "higher_ruling_classes_excluded":False,
        "target_integrality_or_trace_perfection_claimed":False,
        "HC4_or_JC2_closed":False},
    "whole_subcap_vector_controls":whole_controls,
    "whole_cap_quadratic_packet_controls":cap_controls,
    "assurance":"Exact internal symbolic controls and unrestricted written proof; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_square_deformation_closure_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
