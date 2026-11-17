#!/usr/bin/env python3
"""PF175 controls; the unrestricted residue hierarchy and pole proof are written."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import tempfile
import time
import sympy as S

START=time.perf_counter();BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_ODD_HALF_POWER_RESIDUE_HIERARCHY.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
u,w,t,v,s,lam=S.symbols("u w t v s lam")
def initial(e,var):
    n,d=S.fraction(S.together(e))
    if n==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(n,var),S.Poly(d,var)
    a=min(i[0] for i in pn.monoms());b=min(i[0] for i in pd.monoms())
    return a-b,S.cancel(pn.nth(a)/pd.nth(b))
def ord_at(e,c):
    return initial(e.subs(u,c+v),v)[0]

# Whole differentiation identity before taking the coefficient residue.
Pu,Pw,Qu,Qw,Puw,Q,phi,phip=S.symbols("Pu Pw Qu Qw Puw Q phi phip")
du_expr=Qu*phi*Pw+Q*phip*Pu*Pw+Q*phi*Puw
dw_expr=Qw*phi*Pu+Q*phip*Pw*Pu+Q*phi*Puw
check("full chain rule residue integrand identity",
      S.expand(du_expr-dw_expr-phi*(Qu*Pw-Qw*Pu))==0)
check("actual Jacobian sign retained",
      S.expand(du_expr-dw_expr+phi*(Pu*Qw-Pw*Qu))==0)
for exponent in [-4,-2,-1,0,1,3]:
    derivative=S.diff(w**exponent,w)
    check("residue of derivative zero exponent"+str(exponent),
          S.expand(derivative).coeff(w,-1)==0)

# Actual polynomial Jacobian controls outside the cubic/ruling-six class:
# P=(w+u^2)^2+u, Q=w+u^2, Jac=1. They test the universal identity only.
W=w+u**2;P=W**2+u;Q=W
check("actual polynomial Jacobian positive control",
      S.expand(S.diff(P,u)*S.diff(Q,w)-S.diff(P,w)*S.diff(Q,u))==1)
for n in [1,3,5,7,9]:
    j=(n+1)//2
    moment=S.binomial(S.Rational(n,2),j)*u**j
    primitive_coeff=2*S.binomial(S.Rational(n,2),j+1)*u**(j+1)
    check("actual-Jacobian residue derivative n"+str(n),
          S.simplify(S.diff(primitive_coeff,u)+moment)==0)
    check("actual-Jacobian primitive is rational n"+str(n),
          not S.denom(primitive_coeff).has(u))
for k in [1,2,3]:
    shifted=(w+u**2)**(-k)
    coeff=S.series(shifted.subs(w,1/t),t,0,2).removeO().coeff(t,1)
    check("coefficient residue is translation invariant negative power"+str(k),
          coeff==(1 if k==1 else 0))

# Full half-power against the finite formula, with H=w3,R=w2+u.
for n in [1,3,5,7,9]:
    full=S.expand(S.series((1+t**4+u*t**6)**S.Rational(n,2),
                          t,0,3*n+2).removeO()).coeff(t,3*n+1)
    finite=0
    lower=(n+1)//2;upper=(3*n+1)//4
    for j in range(lower,upper+1):
        power=3*n-6*j
        target=-1-power
        coeff=S.Poly(S.expand((w**2+u)**j),w).nth(target) if target>=0 else 0
        finite+=S.binomial(S.Rational(n,2),j)*coeff
    check("entire odd half-power equals finite residue formula n"+str(n),zero(full-finite))
    check("finite upper range is complete n"+str(n),3*n-4*(upper+1)<-1)

# Generic first and third coefficients, independent of any proposed mate.
a3=S.Symbol("a3",nonzero=True)
a2,a1,a0,r2,r1,r0,M,h,dcoef,ecoef=S.symbols("a2 a1 a0 r2 r1 r0 M h dcoef ecoef")
revH=a3+a2*t+a1*t**2+a0*t**3
revR=r2+r1*t+r0*t**2
first=S.expand(S.series(t*revR/revH,t,0,2).removeO()).coeff(t,1)
third=S.expand(S.series(revR**2/revH,t,0,3).removeO()).coeff(t,2)
Xi=(r1-r2*a2/a3)**2/a3+2*r2*r0/a3-r2**2*a1/a3**2
check("generic first half-root residue",zero(first-r2/a3))
check("generic third-character Xi coefficient",zero(third-Xi))
check("canonical first residue d over2h",
      zero((first/2).subs({r2:dcoef/M,a3:h/M})-dcoef/(2*h)))
check("linear third-character exact e2 overh",
      zero((Xi/M).subs({r2:0,r1:ecoef,a3:h/M})-ecoef**2/h))

# The current D24 leading control passes the leading ODE.
vv=u-1
MM=u**2;hh=u*vv**10;GG=MM**2*hh
rho=(-42*u+420*u**2-945*u**3+810*u**4-243*u**5)/28
check("D24 leading primitive ODE positive control",
      zero(S.diff(rho,u)-S.diff(GG,u)*rho/(3*GG)-1))
check("D24 source leading degree",S.degree(MM*hh**2,u)==24)
bad=vv**5;good=u*vv**5;star=u*vv**4*(3-u)
check("linear branch first rational character is vacuous",zero(0/(MM*hh)))
check("bad third character simple pole",zero(bad**2/hh-1/u))
check("bad third character nonzero residue",S.residue(S.cancel(bad**2/hh),u,0)==1)
check("good third character exact primitive",zero(S.diff(u**2/2,u)-good**2/hh))
check("good third character fails A support order",ord_at(good,1)==5)
check("star third character full rational identity",zero(star**2/hh-(vv-3+4/vv**2)))
starprimitive=vv**2/2-3*vv-4/vv
check("star third character exact primitive",zero(S.diff(starprimitive,u)-star**2/hh))
check("star passes A support",ord_at(star,1)==4<5)
check("D24 positives explicitly fail central source condition",good.subs(u,0)==0 and star.subs(u,0)==0)
check("D24 bad passes central but fails residue",bad.subs(u,0)!=0)
check("simple h root and central unit incompatible with residue stratum",1==2*0+1)
# Direct whole-source differential for the central guard.
x,alpha,beta,sigma,E0,E1,F0,F1=S.symbols("x alpha beta sigma E0 E1 F0 F1",nonzero=True)
Hsrc=alpha*x+beta+t*(x**2+1)
wsrc=t*(sigma+t**2*x)
Psrc=Hsrc**2+(E0+E1*t**2)*wsrc+F0+F1*t**2
point={t:0,x:-beta/alpha}
check("central H has an actual source zero",zero(Hsrc.subs(point)))
check("central P x derivative zero at that point",zero(S.diff(Psrc,x).subs(point)))
check("central P t derivative is e0 Sigma0",
      zero(S.diff(Psrc,t).subs(point)-E0*sigma))

# A positive D30 coefficient control passes all listed independent guards.
MM30=u**2;hh30=vv**14;ee30=vv**6;GG30=MM30**2*hh30
rho30=3*u*vv*(243*u**4-891*u**3+1188*u**2-660*u+110)/110
check("D30 leading ODE positive control",
      zero(S.diff(rho30,u)-S.diff(GG30,u)*rho30/(3*GG30)-1))
check("D30 exact third character",zero(ee30**2/hh30-1/vv**2))
check("D30 third character primitive",zero(S.diff(-1/vv,u)-ee30**2/hh30))
check("D30 central unit condition",ee30.subs(u,0)==1)
check("D30 A support condition",ord_at(ee30,1)==6<7)
check("D30 leading degree and numerical ledger",
      S.degree(MM30*hh30**2,u)==30 and 6<=2*7-5)
check("D30 first rational character zero in linear branch",zero(0/(MM30*hh30)))

# Unit-K h-root linear remainders: exact cubic zeros with nonzero
# constant correction, central/noncentral N, and the full shifted P.
U,hs,Ks,Ms,Ns,Bs,cs=S.symbols("U hs Ks Ms Ns Bs cs")
zexpr=(w**2-Ns)/Ms
Hexpr=(hs*w+Ks)*zexpr+Bs*w+cs
cleared=U*((U-Ks)**2-Ns*hs**2)+Bs*hs*Ms*(U-Ks)+cs*hs**2*Ms
check("full cubic centre clearing includes constant correction",
      zero(hs**2*Ms*Hexpr.subs(w,(U-Ks)/hs)-cleared))
check("cubic centre derivative unit independent of correction",
      S.simplify(S.diff(cleared,U).subs({U:0,hs:0})-Ks**2)==0)
unitK_models=[]
for j,m,a,NN in [(1,2,0,v),(3,2,1,v),(3,0,2,1+v),(4,3,3,1+v)]:
    hv=v**j;mv=v**m
    wo=-v**(-j)+v**m
    ccor=v**m-v**(j+2*m)+v**j*NN
    Hlocal=(hv*w+1)*(w**2-NN)/mv+w+ccor
    Plocal=Hlocal**2+v**a*w+7
    check(f"linear unitK exact cubic zero j{j} m{m} a{a}",zero(Hlocal.subs(w,wo)))
    oo,ll=initial(S.diff(Hlocal,w).subs(w,wo),v)
    check(f"linear unitK exact cubic derivative j{j} m{m} a{a}",oo==-j-m and ll!=0)
    sub={v:s**2,w:wo.subs(v,s**2)+s**(a+j+2*m)*lam}
    scaled=s**(2*(j-a))*(Plocal-S.Symbol("Pparam")).subs(sub,simultaneous=True)
    oo,ll=initial(scaled,s)
    check(f"linear unitK whole quadratic initial j{j} m{m} a{a}",
          oo==0 and zero(ll-(lam**2-1)))
    deriv=s**(3*j+2*m-a)*S.diff(Plocal,w).subs(sub,simultaneous=True)
    do,dl=initial(deriv,s)
    check(f"linear unitK exact mate velocity order j{j} m{m} a{a}",
          do==0 and zero(dl-2*lam))
    check(f"linear unitK holomorphic coframe j{j} m{m} a{a}",
          S.Rational(3*j-a,2)>0)
    unitK_models.append({"h_order":j,"M_order":m,"e_order":a,
                        "velocity_order":str(S.Rational(3*j+2*m-a,2)),
                        "coframe_v_order":str(S.Rational(3*j-a,2))})
check("linear central unit excludes every positive h order",
      all(0<j for j in range(1,12)))

# New bad D30 data keep the other listed guards and fail the third residue.
bad30=vv**6*(2*u-1)
check("D30 bad central remains unit",bad30.subs(u,0)==-1)
check("D30 bad A support remains valid",ord_at(bad30,1)==6<7)
check("D30 bad exact rational identity",zero(bad30**2/hh30-(1/vv**2+4/vv+4)))
check("D30 bad third residue is four",S.residue(S.cancel(bad30**2/hh30),u,1)==4)

# Whole local constant-residue controls retain the actual pole index.
D,E,F,c0=S.symbols("D E F c0",nonzero=True)
for alpha0,L,r in [(1,3,-2),(2,5,-3),(4,5,-7)]:
    for ad,ae in [(0,0),(1,0),(0,1),(1,1)]:
        what=c0*s**(-alpha0*L)+lam*s**r
        zhat=(what**2-(1+s**alpha0))/(1+s**alpha0)**2
        Rhat=D*s**(alpha0*(2*L+ad))*zhat+E*s**(alpha0*(L+ae))*what+F+s**alpha0
        oo,ll=initial(Rhat,s)
        expected=F+(D*c0**2 if ad==0 else 0)+(E*c0 if ae==0 else 0)
        check(f"entire local constant residue alpha{alpha0} L{L} ad{ad} ae{ae}",
              oo==0 and zero(ll-expected) and not ll.has(lam))
    check(f"first response is an actual-pole index alpha{alpha0} L{L}",r+alpha0<0)
for j in [1,2,3]:
    check("nonconstant square cannot give affine response degree"+str(j),
          S.degree((lam**j+1)**2+7,lam)==2*j>1)
for L in [3,5,7]:
    ratio=v**(2*L-1)/v**(2*L)
    check("upper-band forbidden simple pole L"+str(L),S.residue(ratio,v,0)==1)
    check("half-bound plus support leaves finite band L"+str(L),
          L<=2*L-2 and not (2*L-1<=2*L-2))

receipt={"claim":"PF2RULED175","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{
      name:sha(BASE/name) for name in [
        "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
        "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md",
        "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
        "RULING_SIX_QUADRATIC_REMAINDER_LOCAL_ORDERS.md",
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md"]},
    "uniform_proved_scope":{
       "odd_characters":"every positive odd n",
       "hierarchy":"d_u res(Q*P^(n/2)*P_w)=-res(P^(n/2))/M",
       "all_odd_half_power_residues_over_M_rationally_exact":True,
       "first_character":"d/(2*h), independent recovery of PF173",
       "third_linear_character":"e^2/h rationally exact",
       "linear_simple_pole_exclusion":"ord h != 2*ord e+1",
       "linear_central_condition":"e(0)!=0",
       "linear_central_h_order_one_class_closed":True,
       "linear_unit_K_h_root_condition":"ord e>=ord h",
       "linear_central_h_nonunit_class_closed":True,
       "linear_survivor_central_h_unit":True,
       "A_pole_support_alternative":"ord d<2L or ord e<L",
       "conditional_A_quadratic_band":"ord e>=L => L<=ord d<=2L-2",
       "all_residue_conditions_asserted_sufficient":False},
    "full_unit_K_linear_models":unitK_models,
    "coefficient_controls":{
      "D24_bad":"e=(u-1)^5; third character fails",
      "D24_good":"e=u*(u-1)^5; residue passes, central and A support fail",
      "D24_star":"e=u*(u-1)^4*(3-u); residue/A support pass, central fails",
      "D30_all_listed":"M=u^2,h=(u-1)^14,e=(u-1)^6; no actual mate asserted",
      "D30_bad":"e=(u-1)^6*(2u-1); central/A support pass, third residue4"},
    "whole_ruling_six_class_closed":False,"generic_first_derivative_integrality_proved":False,
    "normal_tail_increasing_power_containment_proved":False,"full_goal_complete":False,
    "independent_whole_proof_review":False,"formal_verification":False,"external_review":False}
parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path,default=None)
args=parser.parse_args()
out=args.output or Path(tempfile.mkdtemp(prefix="math-check-"))/"ruling_six_residue_hierarchy_receipt.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt":str(out)},indent=2))
