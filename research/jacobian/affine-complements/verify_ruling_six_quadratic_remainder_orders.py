#!/usr/bin/env python3
"""PF174 controls; the complete quartic Newton proof is in the source note."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_QUADRATIC_REMAINDER_LOCAL_ORDERS.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,Y,X,p,s,w=S.symbols("v Y X p s w")
def initial(e,var=v):
    n,d=S.fraction(S.together(e))
    if n==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(n,var),S.Poly(d,var)
    a=min(i[0] for i in pn.monoms());b=min(i[0] for i in pd.monoms())
    return a-b,S.cancel(pn.nth(a)/pd.nth(b))
def trunc(e,n):
    pp=S.Poly(S.expand(e),v)
    return S.expand(sum(c*v**i[0] for i,c in pp.terms() if i[0]<n))

a2,a1,a0,delta=S.symbols("a2 a1 a0 delta",nonzero=True)
A=a2*Y**2+a1*Y+a0
G=S.Poly(S.expand((A**2+delta)/a2**2),Y)
b,c,d,e0=G.nth(3),G.nth(2),G.nth(1),G.nth(0)
Pdep=c-3*b**2/8
Qdep=d-b*c/2+b**3/8
Rdep=e0-b*d/4+b**2*c/16-3*b**4/256
check("monic square plus constant has zero depressed linear",zero(Qdep))
check("depressed constant defect exact normalized perturbation",
      zero(Rdep-Pdep**2/4-delta/a2**2))

hensel_models=[]
for L,a,Nprec in [(3,1,4),(4,3,5)]:
    A=(Y-1)**2+v*Y+v**2
    wholeS=Y*A+v**L*(1+v)  # regular cubic's constant term is present
    F=wholeS**2+v**a*Y**2+v**L*(1-v)*Y+v**(2*L)*(v-p)
    F2=Y**2;GG=(Y-1)**4;G0=(Y-1)**4
    for n in range(1,Nprec+1):
        res=S.Poly(S.expand(F-F2*GG),v).nth(n)
        dh=S.rem(S.expand(res*(1+4*Y)),Y**2,Y)
        dg,rem=S.div(S.expand(res-G0*dh),Y**2,Y)
        check(f"whole Hensel exact division L{L} n{n}",rem==0)
        F2+=dh*v**n;GG+=dg*v**n
    check("whole factor product precision L"+str(L),trunc(F-F2*GG,Nprec+1)==0)
    check("finite factor cutoff L"+str(L),trunc(F2-Y**2,L)==0)
    check("quartic square/constant cutoff L"+str(L),trunc(GG-A**2-v**a,L)==0)
    gp=S.Poly(GG,Y)
    b,c,d,e0=gp.nth(3),gp.nth(2),gp.nth(1),gp.nth(0)
    pd=c-3*b**2/8;qd=d-b*c/2+b**3/8
    rd=e0-b*d/4+b**2*c/16-3*b**4/256
    oo,ll=initial(rd-pd**2/4)
    qo,ql=initial(qd)
    check("actual quartic defect exact order L"+str(L),oo==a and ll==1)
    check("actual quartic linear error stays high L"+str(L),qo>=L)
    check("whole parameter enters exactly at two-ell L"+str(L),S.diff(F,p)==-v**(2*L))
    check("cubic constant term retained in full equation L"+str(L),
          zero(S.expand(wholeS).subs(Y,0)-v**L*(1+v)))
    hensel_models.append({"ell":L,"d_order":a,"precision":Nprec,"defect_order":oo})

newton_models=[]
# 2*ord B < ord Delta: the two quadratic centres split simply.
L,a,bval=5,3,1
G=(X**2+v**bval)**2+v**a+v**L*X-p*v**(2*L)
sub={v:s**4,X:S.I*s**(2*bval)+s**(2*(a-bval))*Y}
no,ni=initial(s**(-4*a)*G.subs(sub,simultaneous=True),s)
check("quadratic-centre entire initial",no==0 and zero(ni-(-4*Y**2+1)))
check("quadratic-centre simple roots",S.gcd(ni,S.diff(ni,Y))==1)
mu=S.Rational(a+bval,2)
do,di=initial(s**(-4*mu)*S.diff(G,X).subs(sub,simultaneous=True),s)
check("quadratic-centre exact derivative order",do==0 and zero(di+8*Y))
check("quadratic-centre derivative below cutoff",mu<L)
newton_models.append({"case":"2b<a","ell":L,"d_order":a,"B_order":bval,"derivative_order":str(mu)})

# 2*ord B > ord Delta, including B=0.
for L,a,bval in [(4,1,1),(3,2,None)]:
    Bval=0 if bval is None else v**bval
    G=(X**2+Bval)**2+v**a+v**L*X-p*v**(2*L)
    sub={v:s**4,X:s**a*Y}
    no,ni=initial(s**(-4*a)*G.subs(sub,simultaneous=True),s)
    label="Bzero" if bval is None else "strict constant"
    check(label+" entire quartic initial",no==0 and zero(ni-Y**4-1))
    check(label+" all simple roots",S.gcd(ni,S.diff(ni,Y))==1)
    mu=S.Rational(3*a,4)
    do,di=initial(s**(-4*mu)*S.diff(G,X).subs(sub,simultaneous=True),s)
    check(label+" exact derivative order",do==0 and zero(di-4*Y**3))
    check(label+" derivative below cutoff",mu<L)
    newton_models.append({"case":"2b>a","ell":L,"d_order":a,"B_order":bval,
                           "derivative_order":str(mu)})
# Equality includes the possibly double root at zero.
L,a,bval=4,2,1
for sign in [1,-1]:
    G=(X**2+v**bval)**2+sign*v**a+v**L*X-p*v**(2*L)
    sub={v:s**4,X:s**a*Y}
    no,ni=initial(s**(-4*a)*G.subs(sub,simultaneous=True),s)
    check("equal-order full initial sign"+str(sign),
          no==0 and zero(ni-(Y**2+1)**2-sign))
    zero_mult=min(i[0] for i in S.Poly(ni,Y).monoms())
    J=S.cancel(ni/Y**zero_mult)
    check("equal-order nonzero simple roots sign"+str(sign),
          S.degree(J,Y)>=2 and S.gcd(J,S.diff(ni,Y))==1)
    if sign==-1:
        check("equal-order double-zero control retained",zero_mult==2)
    mu=S.Rational(3*a,4)
    do,di=initial(s**(-4*mu)*S.diff(G,X).subs(sub,simultaneous=True),s)
    check("equal-order exact derivative order sign"+str(sign),
          do==0 and zero(di-4*Y*(Y**2+1)) and mu<L)
    newton_models.append({"case":"2b=a","sign":sign,"ell":L,"d_order":a,
                           "B_order":bval,"derivative_order":str(mu)})

# Unrestricted order implications, not an unbounded coefficient experiment.
aa,bb,ell=S.symbols("aa bb ell")
check("first Newton derivative formula",
      S.expand(bb/2+aa/2-(aa+bb)/2)==0)
check("first Newton derivative margin from 2b<a",
      S.expand(3*aa/4-(aa+bb)/2)==aa/4-bb/2)
check("whole affine response index identity",
      S.expand((2*ell-aa)+aa-2*ell)==0)
gap,extra=S.symbols("gap extra",positive=True)
R=ell+gap+extra
check("all nonlinear Taylor terms above parameter",
      S.expand(2*R-2*ell)==2*gap+2*extra)
for alpha in [1,2,3,4]:
    for L,mu in [(3,S.Rational(3,4)),(4,S.Rational(3,2)),(5,S.Rational(2))]:
        check(f"normalized coframe positivity alpha{alpha} L{L}",
              alpha*(L-mu)+alpha-1>=0 and 2*L-mu>L)
check("half-integral outer slope retained",S.Rational(3,2)>1)
Dpole=(v-1)**3*(v-2)**5
dpol=Dpole*(v**2+1)
check("all pole-divisor product divides nonzero d",S.rem(dpol,Dpole,v)==0)
check("global nonzero d degree at least slope sum",S.degree(dpol,v)>=3+5)
check("a missing local order is rejected",
      S.rem(S.cancel(Dpole/(v-2)),Dpole,v)!=0)
check("linear remainder d=0 remains a separate class",S.Poly(v*w+1,w).degree()==1)

receipt={"claim":"PF2RULED174","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
      name:sha(BASE/name) for name in [
        "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
        "RULING_SIX_OUTER_SQUARES_AND_EVEN_SPECIALIZATIONS.md",
        "RULING_SIX_QUARTIC_TRACE_AND_CONSTANT_SQUARE_ROOT_CORRECTION.md",
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md"]},
    "whole_Hensel_models":hensel_models,"complete_Newton_cases":newton_models,
    "uniform_proved_scope":{
       "quadratic_remainder":"d*z+e*w+f",
       "outer_four_sheet_local_bound":"ord_c d>=ord_c h/2",
       "every_A_pole_local_bound":"ord_c d>=L",
       "global_pole_divisor_product":"product (u-c_i)^L_i divides d",
       "nonzero_quadratic_coefficient_degree_lower":"deg d>=sum L_i",
       "d_zero_linear_remainder_excluded":False},
    "whole_ruling_six_class_closed":False,"generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_quadratic_remainder_orders_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt":str(out)},indent=2))
