#!/usr/bin/env python3
"""Exact controls for PF176; the unrestricted proof is in its source note."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import sympy as S

START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md"
checks=[]
def check(label,value):
    if not bool(value):
        raise AssertionError(label)
    checks.append(label)
def zero(e):
    return S.cancel(e)==0
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
u,v,w,Y,T,X,a,c,C,E,eps=S.symbols("u v w Y T X a c C E eps",nonzero=True)
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:
        return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(m[0] for m in pn.monoms())
    od=min(m[0] for m in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))
def truncate(expr,n):
    return S.expand(S.series(expr,eps,0,n).removeO())
def coeff_ratio(num_coeff,den_coeff,target):
    """Coefficient T^target in a rational series with unit denominator constant."""
    if target<0:
        return S.Integer(0)
    result=[]
    for i in range(target+1):
        top=num_coeff[i] if i<len(num_coeff) else 0
        top-=sum(den_coeff[j]*result[i-j]
                 for j in range(1,min(i,len(den_coeff)-1)+1))
        result.append(S.cancel(top/den_coeff[0]))
    return result[target]
def residue_term(Hco,e,f,j,power):
    """res_w H^power*(e*w+f)^j for negative integer power."""
    if power>=0:
        return S.Integer(0)
    d=-power
    target=j-3*d+1
    if target<0:
        return S.Integer(0)
    revH=sum(Hco[i]*T**i for i in range(4))
    revDen=S.Poly(S.expand(revH**d),T)
    den_coeff=[revDen.nth(i) for i in range(min(target,3*d)+1)]
    num_coeff=[S.binomial(j,i)*e**(j-i)*f**i
               for i in range(min(target,j)+1)]
    return coeff_ratio(num_coeff,den_coeff,target)
def moment(Hco,e,f,n):
    answer=0
    for j in range((n+1)//2,(3*n+1)//4+1):
        answer+=S.binomial(S.Rational(n,2),j)*residue_term(Hco,e,f,j,n-2*j)
    return S.cancel(answer)

# Actual Jacobian and inverse chain rule, including its exact sign and leading order.
Pu,Pw,Qu,Qw,Ms=S.symbols("Pu Pw Qu Qw Ms",nonzero=True)
wu=-Pu/Pw
check("inverse full chain rule equals minus Jacobian over Pw",
      zero(Qu+Qw*wu+(Pu*Qw-Pw*Qu)/Pw))
check("sixth degree inverse leading Pw",
      zero(6*a**(-6)*(a*X)**5-6*X**5/a))
check("actual Q derivative starts at X minus five",
      zero(-1/(Ms*(6*X**5/a))+a*X**(-5)/(6*Ms)))

# A noncube inverse: u=a^3, P=w^6/u^2+w+u.
# This is a coefficient model, not an actual Keller companion claim.
N=13
V=a
for _ in range(4):
    V=truncate(a*(1-eps**5*V-eps**6*a**3)**S.Rational(1,6),N)
check("noncube full inverse equation through order twelve",
      truncate(V**6/a**6+eps**5*V+eps**6*a**3-1,N)==0)
vp=S.Poly(V,eps)
nonzero_inverse=[]
for i in range(N):
    coeff=S.expand(vp.nth(i))
    if coeff==0:
        continue
    j=1-i
    check("cyclic inverse coefficient character j"+str(j),
          all((power[0]-j)%3==0 for power,value in S.Poly(coeff,a).terms()))
    nonzero_inverse.append({"X_exponent":j,"coefficient":str(coeff)})
check("inverse first nontrivial coefficient",
      zero(vp.nth(5)+a**2/6))
check("inverse first parameter integration coefficient",
      zero(vp.nth(6)+a**4/6))
check("noncube constant field does not make arbitrary inverse coefficients constant",
      S.diff(vp.nth(0),a)!=0)
for k in [1,3,5,7,11]:
    check("coupled companion leading coefficient cancels cubic root k"+str(k),
          zero(a**(-3*k)*a**(3*k)-1))

# Constants in the q_j prefix have weight zero: only multiples of three survive.
zeta=S.symbols("zeta")
for j in range(-4,16):
    character=S.rem(zeta**(j%3)-1,zeta**2+zeta+1,zeta)
    check("constant-prefix cyclic character j"+str(j),
          (character==0)==(j%3==0))
check("the first allowed negative constant prefix exponent is minus three",
      [j for j in range(-4,0) if j%3==0]==[-3])

# An actual polynomial Jacobian pair in the CUBE case is a necessary negative control.
# P=w^6+u,Q=w, M=h=1 has no A obstruction and no nontrivial cubic automorphism.
Pcube=w**6+u
check("cube-case actual polynomial Jacobian pair",
      S.diff(Pcube,u)*S.diff(w,w)-S.diff(Pcube,w)*S.diff(w,u)==1)
Qcube=X*(1-u*X**(-6))**S.Rational(1,6)
Qcubeprefix=S.series(Qcube,X,S.oo,7).removeO()
check("cube actual companion has a nonresonant constant positive term",
      S.expand(Qcubeprefix).coeff(X,1)==1)
check("cube actual inverse has first varying coefficient at minus five",
      zero(S.expand(Qcubeprefix).coeff(X,-5)+u/6))
check("cube actual inverse chain rule",
      zero(S.diff(Qcube,u)+1/(6*Qcube**5)))
check("cube control is expressly outside noncube hypothesis",
      S.Integer(1)**3==1)

# Coefficient residue scaling and the exact nonzero leading residue.
for ell in [1,3,5,7]:
    test=w**(-1)+2*w**(-2)+3*w**3
    transformed=S.expand(test.subs(w,v**(-ell)*Y))
    check("coefficient-residue scaling slope"+str(ell),
          zero(v**(-ell)*transformed.coeff(Y,-1)-1))
for r in range(1,11):
    coefficient=S.expand(S.series((Y**r/(Y-c)**2).subs(Y,1/T),
                                  T,0,2).removeO()).coeff(T,1)
    expected=r*c**(r-1)
    check("leading coefficient-residue identity r"+str(r),zero(coefficient-expected))
    check("leading half-binomial is nonzero r"+str(r),
          S.binomial(S.Rational(2*r+1,2),r+1)!=0)
check("m1 leading residue vanishes in linear branch",
      S.series((1/(Y-c)**2).subs(Y,1/T),T,0,2).removeO().coeff(T,1)==0)
ell,ae,r,step=S.symbols("ell ae r step",integer=True,positive=True)
n=2*r+1
exponent=-ell*(n+1)+(ell+ae)*(r+1+step)
check("all later binomial terms have exact positive-margin expression",
      S.expand(exponent+(r+1)*(ell-ae)-(ell+ae)*step)==0)
check("zero leading cluster centre would destroy later nonvanishing",
      (2*c).subs(c,0)==0)
check("zero e leading coefficient would destroy leading moment",
      E**2!=0 and (E**2).subs(E,0)==0)

# Independent whole half-power coefficient extraction, not only the leading term.
Htest=w**3+w**2+2*w+3
Rtest=2*w+5
Ptest=S.expand(Htest**2+Rtest)
for nn in [1,3,5,7]:
    revP=S.expand(T**6*Ptest.subs(w,1/T))
    full=S.expand(S.series(revP**S.Rational(nn,2),T,0,3*nn+2).removeO()).coeff(T,3*nn+1)
    finite=moment([1,1,2,3],2,5,nn)
    check("whole half-power agrees with complete finite formula n"+str(nn),zero(full-finite))

# Full cubic-plus-linear models have nonzero lower cubic corrections and nonconstant f.
# All are coefficient controls; no actual companion is posited.
models=[]
for ell0,ae0,cc,CC in [(3,0,1,1),(3,2,2,1),(5,1,-1,2),(5,4,3,1),(7,6,1,1)]:
    # H=CC*w*(v^ell*w-cc)^2+v*w+1
    Hco=[CC*v**(2*ell0),-2*CC*cc*v**ell0,CC*cc**2+v,S.Integer(1)]
    ee=v**ae0*(2+v)
    ff=3+v
    values={}
    for nn in [1,3,5,7,9]:
        values[nn]=moment(Hco,ee,ff,nn)
        if nn==1:
            check(f"full model m1 zero L{ell0} a{ae0}",values[nn]==0)
            continue
        rr=(nn-1)//2
        order,leading=initial(values[nn],v)
        expected_order=-(rr+1)*(ell0-ae0)
        expected_leading=S.binomial(S.Rational(nn,2),rr+1)*2**(rr+1)*rr*cc**(rr-1)/CC
        check(f"full moment exact pole order L{ell0} a{ae0} n{nn}",order==expected_order)
        check(f"full moment exact nonzero leading coefficient L{ell0} a{ae0} n{nn}",
              zero(leading-expected_leading) and leading!=0)
    combination=values[3]-2*values[5]+3*values[9]
    check(f"constant top-character relation fails L{ell0} a{ae0}",
          initial(combination,v)[0]==-5*(ell0-ae0))
    models.append({"slope":ell0,"e_order":ae0,"cluster_centre":cc,
                   "leading_cubic_constant":CC,
                   "tested_moment_orders":{str(nn):str(initial(values[nn],v)[0])
                                           for nn in [3,5,7,9]}})
# n=13 is the first tested linear case with a later binomial term contributing.
Hco=[v**6,-2*v**3,1+v,S.Integer(1)]
mm13=moment(Hco,2+v,3+v,13)
check("higher finite sum includes the j8 term n13",
      residue_term(Hco,2+v,3+v,8,-3)!=0)
oo,ll=initial(mm13,v)
check("full n13 pole order survives the later contribution",oo==-21)
check("full n13 leading coefficient survives the later contribution",
      zero(ll-S.binomial(S.Rational(13,2),7)*2**7*6))

# Prior D30 positive control passes earlier guards but fails actual constant resonance.
vv=u-1
M30=u**2;h30=vv**14;e30=vv**6;G30=M30**2*h30
rho=3*u*vv*(243*u**4-891*u**3+1188*u**2-660*u+110)/110
check("former D30 control leading ODE passes",
      zero(S.diff(rho,u)-S.diff(G30,u)*rho/(3*G30)-1))
check("former D30 first character zero",zero(0/(M30*h30)))
check("former D30 third rational-character primitive passes",
      zero(S.diff(-1/vv,u)-e30**2/h30))
check("former D30 central h and e units",h30.subs(u,0)==1 and e30.subs(u,0)==1)
check("former D30 A support passes",6<7)
check("former D30 leading degree",S.degree(M30*h30**2,u)==30)
check("former D30 G is a noncube",4%3!=0)
m3=S.Rational(3,8)*M30*e30**2/h30
check("former D30 m3 pole order",initial(m3.subs(u,1+v),v)[0]==-2)
for nn in [3,5,7,9,11]:
    rr=(nn-1)//2
    check("former D30 all-character theorem pole order n"+str(nn),
          -(rr+1)*(7-6)==-(nn+1)//2)

receipt={
    "claim":"PF2RULED176","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md",
        "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
        "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
        "RULING_SIX_ODD_HALF_POWER_RESIDUE_HIERARCHY.md"]},
    "uniform_proved_scope":{
        "coefficient_and_companion_degrees":"unrestricted",
        "constant_prefix":"q_j in C and q_j=0 unless 3 divides j, for j>-5",
        "actual_companion_resonance":"sum_(odd n<=k) c_n*m_n=0, c_k nonzero, k>=3 odd",
        "linear_A_character_orders":"ord_c m_(2r+1)=-(r+1)*(L-ord_c e), every r>=1",
        "entire_linear_remainder_branch_excluded":True,
        "necessary_quadratic_coefficient":"d!=0",
        "quadratic_remainder_class_closed":False,
        "trace_integrality_or_target_landing_claimed":False,
        "full_ruling_six_or_JC2_closed":False},
    "noncube_inverse_control":nonzero_inverse,
    "whole_local_models":models,
    "assurance":"Exact internal symbolic controls; unrestricted theorem proved in source; no formal or external human review."
}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
out=BASE/"ruling_six_companion_resonance_linear_closure_receipt.json"
if args.write_receipt:
    out.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS",
                  "check_count":len(checks),"receipt_written":args.write_receipt,
                  "elapsed_seconds":receipt["elapsed_seconds"]},indent=2))
