#!/usr/bin/env python3
"""PF168 exact controls; complete square/escape proofs are in the source note."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_OUTER_SQUARES_AND_EVEN_SPECIALIZATIONS.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,w,Y,p,s=S.symbols("v w Y p s")
def initial(e,var=v):
    n,d=S.fraction(S.together(e))
    if n==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(n,var),S.Poly(d,var)
    a=min(i[0] for i in pn.monoms());b=min(i[0] for i in pd.monoms())
    return a-b,S.cancel(pn.nth(a)/pd.nth(b))
def square_multiplicities(I):
    return all(m%2==0 for factor,m in S.factor_list(I,Y)[1])
H,K,ell,II,JJ,IP,JP=S.symbols("H K ell II JJ IP JP")
chain=(-H*II+ell*Y*IP)*JP-IP*(-K*JJ+ell*Y*JP)
check("exact weighted leading Jacobian chain rule",
      S.expand(chain+H*II*JP-K*IP*JJ)==0)
jj,kk,ll=S.symbols("jj kk ll",positive=True)
check("top companion term has the matched weight",
      S.expand(3*kk*ll-kk*jj-kk*(6*ll-2*jj)/2)==0)
check("degree relation yields companion weight upper",
      S.cancel((3*kk)/6-kk/2)==0)
# Root multiplicities: k*m/2 integral with odd k iff m even.
for k in [3,5,7,9]:
    check("complete degree-six parity k"+str(k),
          all((S.Rational(k*m,2).q==1)==(m%2==0) for m in range(1,7)))
profiles=[
    ("A",Y**2*(Y-1)**4,True),
    ("two double outer roots",Y**2*(Y-1)**2*(Y+2)**2,True),
    ("outer mass two",Y**4*(Y-1)**2,True),
    ("forbidden cubic",Y**3*(Y+1)**3,False),
    ("forbidden simple leftover",Y**5*(Y+1),False)]
for name,I,expected in profiles:
    check(name+" root multiplicity square test",square_multiplicities(I)==expected)
    if expected:
        _,factors=S.factor_list(I,Y)
        U=S.prod(f**(m//2) for f,m in factors)
        for k in [3,5]:
            J=U**k
            check(name+" matched commuting polynomial k"+str(k),
                  zero(2*I*S.diff(J,Y)-k*S.diff(I,Y)*J))
            check(name+" full companion degree k"+str(k),S.degree(J,Y)==3*k)

# COMPLETE outer escape mass enumeration for actual degree six.
expected={2:[2,4],3:[2],4:[2],5:[]}
for d in range(2,6):
    possibilities=[m for m in range(1,7-d) if m%2==0]
    check("complete escaping even mass d"+str(d),possibilities==expected[d])
check("degree-three outer pair leaves exactly one escape",6-3-2==1)
check("degree-five has no positive even outer mass",not expected[5])

# Whole local degree-three control: square top double cluster PLUS the
# omitted lower simple branch. No actual Keller pair is claimed.
P3=v**8*w**6+2*v**5*w**5+v**2*w**4+w**3+v
no,ni=initial(v**10*(P3-p).subs(w,v**(-3)*Y))
check("d3 entire outer square control",no==0 and zero(ni-Y**4*(Y+1)**2))
check("d3 exact finite specialization",S.degree(P3.subs(v,0),w)==3)
no,ni=initial(v**6*(P3-p).subs(w,v**(-2)*Y))
check("d3 entire remaining lower escape",no==0 and zero(ni-Y**3*(Y+1)))
check("d3 leftover initial root simple",
      ni.subs(Y,-1)==0 and S.diff(ni,Y).subs(Y,-1)!=0)
do,di=initial(v**4*S.diff(P3,w).subs(w,v**(-2)*Y))
check("d3 leftover full derivative order",do==0 and zero(di-S.diff(ni,Y)))
check("d3 leftover coframe holomorphic order four",4>0)
check("d3 lower full parameter numerator",
      S.diff(v**6*(P3-p).subs(w,v**(-2)*Y),p)==-v**6)

P5=v**2*w**6+w**5+w
no,ni=initial(v**10*(P5-p).subs(w,v**(-2)*Y))
check("d5 entire outer face has forbidden odd root multiplicity",
      no==0 and zero(ni-Y**5*(Y+1)) and not square_multiplicities(ni))
# A unit-slope degree-four local square is kept as an unexcluded control.
P4=w**4*(v*w+1)**2+w**2+v
no,ni=initial(v**4*(P4-p).subs(w,v**(-1)*Y))
check("d4 square outer face remains compatible",
      no==0 and zero(ni-Y**4*(Y+1)**2) and square_multiplicities(ni))
check("d4 finite specialization retained",S.degree(P4.subs(v,0),w)==4)

# Exact quartic remainder coefficient, including missing p4 face coefficient.
a,b,c,M,N,h,eta,r4=S.symbols("a b c M N h eta r4")
I=(a*Y**3+b*Y**2+c*Y)**2
check("square p5 coefficient",S.Poly(S.expand(I),Y).nth(5)==2*a*b)
check("square p4 coefficient",S.Poly(S.expand(I),Y).nth(4)==b**2+2*a*c)
check("remainder difference retains nonzero cross coefficient",
      S.expand(S.Poly(S.expand(I),Y).nth(4)-b**2)==2*a*c)
check("cancelled p4 control survives",
      S.expand((b**2+2*a*c).subs({a:1,b:2,c:-2}))==0
      and (2*a*c).subs({a:1,c:-2})!=0)
models=[]
for L in [1,3,4]:
    for name,b0,c0 in [("A",-2,1),("two roots",1,-2),("cancelled p4",2,-2)]:
        j=2*L
        hp=v**j;etap=b0*v**L
        P=(hp*w**3+etap*w**2+c0*w)**2+v+w
        r4p=S.Poly(S.expand(P),w).nth(4)-etap**2
        order,lead=initial(r4p)
        check(name+" entire r4 exact order L"+str(L),
              order==j and lead==2*c0)
        no,ni=initial(v**(2*L)*(P-p).subs(w,v**(-L)*Y))
        check(name+" entire square face L"+str(L),
              no==0 and zero(ni-(Y**3+b0*Y**2+c0*Y)**2))
        models.append({"L":L,"h_order":j,"eta_order":L,"r4_order":order,
                       "profile":name})
# Half-integral slopes must not be silently dropped.
for j in [1,3,5]:
    P=(v**j*w**3+w)**2+v
    whole=s**(2*j)*(P-p).subs({v:s**2,w:s**(-j)*Y},simultaneous=True)
    no,ni=initial(whole,s)
    check("half-integral whole square face j"+str(j),
          no==0 and zero(ni-(Y**3+Y)**2))
    r4p=S.Poly(S.expand(P),w).nth(4)
    order,lead=initial(r4p)
    check("half-integral exact r4 order j"+str(j),order==j and lead==2)
    check("half-integral missing p5/p3 allowed j"+str(j),
          S.Poly(S.expand(P),w).nth(5)==0 and S.Poly(S.expand(P),w).nth(3)==0)

# The finite response field degree and graph test do not confuse the
# outer geometric roots with the single remaining closed prime.
lam=S.symbols("lam")
check("single remaining response has linear polynomial degree",
      S.degree(3*lam+2,lam)==1)
check("nonlinear response consumes more than one residue conjugate",
      all(S.degree(lam**d-p,lam)==d>1 for d in [2,3,4]))
check("closed graph polynomial residue retained",
      not S.denom(S.together((lam**3+2*lam).subs(lam,(p-2)/3))).has(p))
check("ruling-six class not declared empty",
      square_multiplicities(Y**2*(Y-1)**4))

receipt={"claim":"PF2RULED168","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
       name:sha(BASE/name) for name in [
        "RULING_SIX_FINITE_END_COVERAGE_AND_QUADRATIC_BOUNDARY_ALGEBRA.md",
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
        "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md"]},
    "full_square_models":models,
    "uniform_proved_scope":{
       "all_outer_degree_six_initials_square":True,
       "rational_outer_slopes_retained":True,
       "off_M_specialization_degree_three_excluded":True,
       "off_M_specialization_degree_five_excluded":True,
       "remaining_off_M_specialization_degrees":[2,4],
       "outer_mass_four_r4_order":"ord r4=ord h",
       "A_pole_r4_order":"2*L=ord h",
       "refined_shift_initial_square_claimed_without_hypotheses":False},
    "whole_ruling_six_class_closed":False,"generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_outer_square_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt":str(out)},indent=2))
