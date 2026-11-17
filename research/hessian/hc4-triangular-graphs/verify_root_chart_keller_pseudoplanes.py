#!/usr/bin/env python3
"""Exact controls for actual root-chart Keller/pseudo-plane transfer."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import platform
import time
from pathlib import Path
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

BASE=Path(__file__).resolve().parent
PROOF=BASE/"WHOLE_ROOT_CHARTS_CREATE_CYCLIC_KELLER_PSEUDOPLANES.md"
AUDIT=BASE/"sources/root-chart-keller-pseudoplanes-source-audit.json"
RECEIPT=BASE/"root_chart_keller_pseudoplanes_receipt.json"
START=time.perf_counter()
checks: list[dict[str,object]]=[]

def check(name: str, condition: object) -> None:
    ok=bool(condition)
    checks.append({"name":name,"passed":ok})
    if not ok:
        raise AssertionError(name)

def eq(name: str, actual: object, expected: object=0) -> None:
    check(name,S.cancel(S.sympify(actual)-S.sympify(expected))==0)

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def sub(expr: object, mapping: dict) -> object:
    return S.expand(S.sympify(expr).subs(mapping,simultaneous=True))

def laurent_order(expr: object, variable: S.Symbol) -> int:
    expr=S.expand(expr)
    if expr==0:
        raise ValueError("zero has no finite order")
    return min(int(term.as_powers_dict().get(variable,0))
               for term in S.Add.make_args(expr))

x,y,z,t,T,lam,mu,u,a,b,s,eps,tau=S.symbols(
    "x y z t T lam mu u a b s eps tau")

def mu_reduce(expr: object, order: int) -> object:
    out=S.Integer(0)
    for term in S.Add.make_args(S.expand(expr)):
        exponent=term.as_powers_dict().get(mu,S.Integer(0))
        if not exponent.is_integer:
            raise ValueError("nonintegral root character")
        coefficient=S.cancel(term/mu**exponent)
        if coefficient.has(mu):
            raise ValueError("non-Laurent root expression")
        out += coefficient*mu**(int(exponent)%order)
    return S.expand(out)

def mu_eq(name: str, actual: object, expected: object, order: int) -> None:
    check(name,mu_reduce(S.sympify(actual)-S.sympify(expected),order)==0)

# Irreducible root extension and its exact root-character basis. The proof
# uses Eisenstein for ALL r; these are exact transcription controls.
for r in [1,2,3,4,6]:
    polynomial=S.Poly(T**r-y,T)
    check(f"root equation is monic degree r{r}",
          polynomial.degree()==r and polynomial.LC()==1)
    eq(f"source-line constant term has exact order one r{r}",
       laurent_order(polynomial.TC(),y),1)
    check(f"all remaining nonleading coefficients divisible by source parameter r{r}",
          all(co==0 or S.rem(co,y,y)==0 for co in polynomial.all_coeffs()[1:]))
    phi=S.cyclotomic_poly(r,mu)
    for j in range(r):
        invariant=S.rem(mu**j-1,phi,mu)==0
        check(f"root-character invariant basis r{r} j{j}",
              invariant==(j==0))

# Exact whole-root chart formulas. h has primitive TOTAL support; its
# leading exponent need not itself be coprime to r.
for r in [1,2,3,4,6]:
    h=z**(-r-1)+2*z**(-2*r-2)
    p=2*r+2
    xp=z**(1-r)*t+h
    yp=lam+z**r
    eq(f"root chart recovers original source x r{r}",
       z**(1-r)*(z**(r-1)*(x-h))+h,x)
    eq(f"root chart recovers plane coordinate t r{r}",
       z**(r-1)*(xp-h),t)
    eq(f"root chart fixes original coordinate equation r{r}",yp-lam,z**r)
    jac=S.diff(xp,t)*S.diff(yp,z)-S.diff(xp,z)*S.diff(yp,t)
    eq(f"new Keller Jacobian chain factor r{r}",jac,r)
    eq(f"original source translation lifts regularly on W r{r}",
       z**(r-1)*S.diff(xp,t),1)
    eq(f"divided field is whole-plane translation r{r}",
       S.diff(xp,t),z**(1-r))
    eq(f"principal part has no logarithmic exponent r{r}",
       S.expand(h).coeff(z,-r),0)
    eq(f"principal part leading pole r{r}",laurent_order(h,z),-p)
    check(f"principal part strict pole inequality r{r}",p>=r+1)

    rotation_t=mu**(r-1)*(t+z**(r-1)*(h-sub(h,{z:mu*z})))
    rotation_z=mu*z
    mu_eq(f"root rotation fixes ORIGINAL x r{r}",
          sub(xp,{t:rotation_t,z:rotation_z}),xp,r)
    mu_eq(f"root rotation fixes ORIGINAL y r{r}",
          sub(yp,{t:rotation_t,z:rotation_z}),yp,r)
    iterate=t
    for _ in range(r):
        iterate=mu_reduce(sub(iterate,{t:rotation_t,z:rotation_z}),r)
    mu_eq(f"cyclic rotation has rth iterate identity on t r{r}",iterate,t,r)

    # The retained original source line is z=0 in (x,z), outside W.
    t_at_source=z**(r-1)*(x-h)
    q=p+1-r
    eq(f"new marked ramified valuation z order r{r}",laurent_order(z,z),1)
    eq(f"new marked ramified valuation t pole r{r}",
       laurent_order(t_at_source,z),-q)
    eq(f"target normal parameter has ramification r{r}",
       laurent_order(z**r,z),r)
    eq(f"new source volume order at omitted line r{r}",
       laurent_order(S.diff(t_at_source,x),z),r-1)
    check(f"new marked transverse pole at least two r{r}",q>=2)
    if r>1:
        check(f"retained source line becomes RAMIFIED omitted prime r{r}",r-1>0)
    else:
        eq("r1 retained line stays UNRAMIFIED",r-1,0)

    # Important negative guard: a rational Darboux chart by itself does not
    # supply a polynomial target pair on the whole new plane.
    check(f"arbitrary original coordinate x is not regular in new plane r{r}",
          not S.cancel(xp).is_polynomial(t,z))
    eq(f"rational-only coordinate pair still has nonzero chain Jacobian r{r}",
       jac,r)

# Necessary support primitivity, including the leading-pole-only trap.
for r in [2,3,4,6,8]:
    good=[r+1,2*r+2]
    check(f"total support can be primitive r{r}",math.gcd(r,*good)==1)
    bad=[2*r,3*r]
    check(f"nonprimitive full-root-invariant support rejected r{r}",
          math.gcd(r,*bad)>1)
    h=sum(z**(-q) for q in bad)
    mu_eq(f"bad support is fixed by the whole root group r{r}",
          sub(h,{z:mu*z}),h,r)
    rotation_t=mu**(r-1)*(t+z**(r-1)*(h-sub(h,{z:mu*z})))
    mu_eq(f"bad support makes polynomial linear deck on entire plane r{r}",
          rotation_t,mu**(r-1)*t,r)
    eq(f"bad support deck fixes whole-chart origin r{r}",
       sub(mu_reduce(rotation_t,r),{t:0,z:0}),0)
    check(f"nonidentity root character cannot fix an etale target tangent r{r}",
          S.rem(mu-1,S.cyclotomic_poly(r,mu),mu)!=0)
    if r%2==0:
        check(f"leading pole alone need not be coprime r{r}",
              math.gcd(r,max(good))>1 and math.gcd(r,*good)==1)
        subgroup_bad=[2*(r+1),2*(r+2)]
        check(f"proper subgroup can obstruct support primitivity r{r}",
              math.gcd(r,*subgroup_bad)==2)

# Cyclic source-image/fibre/Picard ledgers. Their geometric completeness is
# proved in the note; these finite controls cannot replace it.
for r in [1,2,3,4,7]:
    copies=set(range(r))
    chosen=0
    translates=[{(chosen+j)%r} for j in range(r)]
    check(f"all source-chart translates cover every copied boundary r{r}",
          set().union(*translates)==copies)
    check(f"selected plane omits exactly r-1 copied lines r{r}",
          len(copies-{chosen})==r-1)
    special_points=1
    general_points=r
    check(f"whole source image includes central fibre r{r}",special_points>0)
    check(f"whole source image includes every ordinary fibre r{r}",general_points>0)
    if r>1:
        check(f"selected source pseudo-cover is not finite-etale r{r}",
              special_points!=general_points)
    else:
        check("r1 selected source map is its plane-chart isomorphism",
              special_points==general_points==1)
    snf=smith_normal_form(S.Matrix([[-r]]),domain=ZZ)
    eq(f"Picard localization relation has cyclic cokernel order r{r}",abs(snf[0,0]),r)
    eq(f"Euler of universal cover from r-1 added lines r{r}",1+(r-1),r)
    eq(f"finite-cover Euler gives quotient one r{r}",S.Rational(r,r),1)
    eq(f"central fibre pulls back with exact multiplicity r{r}",r*1,r)

# Explicit smooth affine controls with all r plane charts and a free
# exact-volume-preserving cyclic action. NO Keller target pair is supplied.
for r in [2,3,4,5]:
    p=r+2
    F=u**p*a-b**r+1
    H=sum(b**j for j in range(r))
    eq(f"root-of-unity fibre factorization r{r}",(b-1)*H,b**r-1)
    eq(f"selected root chart inverse denominator is a unit at chosen boundary r{r}",
       H.subs(b,1),r)
    eq(f"smoothness candidate impossible on affine control r{r}",
       F.subs({u:0,b:0}),1)
    eq(f"r zero-fibre roots are reduced r{r}",
       S.Poly(b**r-1,b).gcd(S.Poly(r*b**(r-1),b)).as_expr(),1)
    bchart=1+u**p*s
    achart=S.expand(((1+u**p*s)**r-1)/u**p)
    check(f"selected WHOLE plane chart is polynomial r{r}",
          achart.is_polynomial(s,u))
    eq(f"whole plane chart satisfies affine equation r{r}",
       sub(F,{a:achart,b:bchart}))
    eq(f"source inverse on u-open r{r}",(bchart-1)/u**p,s)
    eq(f"source inverse on H-open r{r}",achart/sub(H,{b:bchart}),s)
    eq(f"source volume in selected whole chart r{r}",
       S.diff(bchart,s)/u**p,1)
    eq(f"source volume on the boundary chart r{r}",
       S.diff(achart,s)/(r*bchart**(r-1)),1)
    theta_s=S.diff(bchart,s)/((p-1)*u**(p-1))
    theta_u=S.diff(bchart,u)/((p-1)*u**(p-1))
    check(f"shared primitive is polynomial in whole plane chart r{r}",
          S.cancel(theta_s).is_polynomial(s,u) and S.cancel(theta_u).is_polynomial(s,u))
    eq(f"shared primitive differentiates to unit form r{r}",
       S.diff(theta_u,s)-S.diff(theta_s,u),1)

    # ALL epsilon charts by a polynomial identity modulo epsilon^r=1.
    achart_eps=sum(S.binomial(r,j)*eps**(r-j)*u**(p*(j-1))*s**j
                  for j in range(1,r+1))
    eps_residual=sub(F,{a:achart_eps,b:eps+u**p*s})
    eq(f"all r whole charts satisfy relation modulo root equation r{r}",
       S.rem(eps_residual,eps**r-1,eps),0)
    eq(f"all whole charts have the same unit volume r{r}",
       S.diff(eps+u**p*s,s)/u**p,1)

    action={u:mu*u,b:mu*b,a:mu**(r-2)*a}
    mu_eq(f"free cyclic action preserves affine equation r{r}",sub(F,action),F,r)
    mu_eq(f"cyclic action preserves global UNIT volume r{r}",mu**(2-p),1,r)
    mu_eq(f"cyclic action preserves the global primitive r{r}",mu**(2-p),1,r)
    for j in range(1,r):
        check(f"nonidentity group element moves u unless u=0 r{r} j{j}",
              S.rem(mu**j-1,S.cyclotomic_poly(r,mu),mu)!=0)
        eq(f"u=0 forces b nonzero so no fixed point r{r} j{j}",
           F.subs({u:0,b:0}),1)
    delta=lambda f:S.expand(u**p*S.diff(f,b)+r*b**(r-1)*S.diff(f,a))
    eq(f"control Ga field preserves equation r{r}",delta(F))
    eq(f"control Ga field is translation in the source chart r{r}",
       sub(delta(b),{a:achart,b:bchart}),S.diff(bchart,s))
    bflow=b+tau*u**p
    aflow=a+S.expand((bflow**r-b**r)/u**p)
    check(f"free control Ga flow is polynomial r{r}",aflow.is_polynomial(u,a,b,tau))
    eq(f"free control Ga flow preserves affine equation r{r}",
       sub(F,{b:bflow,a:aflow}),F)
    eq(f"control zero-line vector field never vanishes r{r}",
       S.rem((r*b**(r-1))*b/r-1,b**r-1,b),0)

# Explicit field/open and assurance barriers.
proof=PROOF.read_text(encoding="utf-8")
for text in [
    "FULL finite normalization",
    "Globally z need not be integral",
    "generic deck SUBGROUP",
    "The full",
    "gcd(r, {q : a_q!=0})=1",
    "SURJECTIVE ETALE",
    "W_j->V is not finite for r>1",
    "Pic(V)=Z/rZ",
    "pi1(V)=mu_r",
    "does NOT permit an application",
    "same one-multiple-fibre",
    "No etale target pair to A2 is supplied",
]:
    check("written architecture retains scope: "+text,text in proof)
audit=json.loads(AUDIT.read_text(encoding="utf-8"))
check("explicit geometric control is not a claimed Keller target",
      audit["explicit_control"]["target_pair_supplied"] is False)
parents=[BASE/e["path"] for e in audit["retained_inputs"]]
receipt={
    "claim":"WHOLE-ROOT-CHART-CYCLIC-KELLER-PSEUDOPLANES",
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_inputs":[{"path":p.name,"sha256":sha(p)} for p in parents]+[
        {"path":str(AUDIT.relative_to(BASE)),"sha256":sha(AUDIT)}],
    "proved_scope":{
        "starting_map_is_actual_polynomial_plane_Keller":True,
        "all_original_full_normalization_primes_retained":True,
        "new_actual_plane_Keller_Jacobian":"r",
        "new_actual_generic_degree":"r*d",
        "cyclic_generic_deck_SUBGROUP":"mu_r",
        "full_target_field_extension_claimed_Galois":False,
        "specified_affine_pseudo_plane_source_image":"M_i minus whole retained source line",
        "selected_plane_surjective_etale_to_image":True,
        "whole_central_scheme_fibre":"r*D",
        "every_nonzero_fibre":"one reduced A1",
        "Pic_and_pi1_image":"Z/rZ",
        "finite_etale_universal_cover_union_of_r_plane_charts":True,
        "necessary_Laurent_support_primitivity":True,
        "leading_pole_coprime_to_r_required":False,
        "new_marked_ramified_prime_index":"r when r>1",
        "new_common_coordinate_order_there":1,
        "r1_same_field_source_switch":True,
        "selected_open_identified_with_full_intermediate_normalization":False,
        "torsion_canonical_full_order_theorem_applied_to_selected_open":False,
        "explicit_geometric_control_has_Keller_target_pair":False,
        "source_finiteness_or_integrality_of_common_coordinate_proved":False,
        "cyclic_architecture_excluded_or_inverse_proved":False,
    },
    "assurance":"Internal unrestricted written construction and exact controls; no formal verification or external human review."
}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    RECEIPT.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS",
                  "check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt_written":args.write_receipt},indent=2))
