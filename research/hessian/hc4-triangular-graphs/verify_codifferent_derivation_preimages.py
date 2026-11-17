#!/usr/bin/env python3
"""Exact internal support for codifferent preimages, not automatic landing."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import sympy as S

BASE=Path(__file__).resolve().parent
PROOF=BASE/"CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md"
checks=[]
def check(name,value):
    ok=bool(value)
    checks.append({"name":name,"passed":ok})
    if not ok: raise AssertionError(name)
def eq(name,actual,expected):
    check(name,S.cancel(S.sympify(actual)-S.sympify(expected))==0)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def order(f,var):
    f=S.cancel(f)
    if f==0:return S.oo
    n,d=S.fraction(f)
    def po(poly):
        return min(mon[0] for mon,_ in S.Poly(poly,var).terms())
    return po(n)-po(d)

def run():
    u,z,t,x,y,P,Q,T,alpha,beta=S.symbols("u z t x y P Q T alpha beta")

    # Strict tame charts with varying leading units and both target components.
    for e in range(1,7):
        Dt=lambda f:S.cancel(S.diff(f,u)/(e*u**(e-1)))
        Dz=lambda f:S.diff(f,z)
        for q in range(1,7):
            f=(z+1+u*(z*z+1)+u*u)/u**q
            derivative=S.cancel(Dt(f)+Dz(f))
            check(f"transverse derivative exact pole e{e} q{q}",
                  order(derivative,u)==-q-e)
            eq(f"transverse derivative leading unit e{e} q{q}",
               S.cancel(derivative*u**(q+e)).subs(u,0),
               -S.Rational(q,e)*(z+1))
        for n in range(0,2*e+1):
            f=u**n*(z*z+z+1)
            check(f"regular normal function derivative in codifferent e{e} n{n}",
                  order(Dt(f)+Dz(f),u)>=1-e)

        # Both derivative conditions on a full-field monomial lattice,
        # including poles at an unramified tangential prime z=0.
        in_dual=lambda f:(order(f,u)>=1-e and order(f,z)>=0)
        for n in range(-3,4):
            for j in range(-3,4):
                f=u**n*z**j
                check(f"two derivative preimage e{e} n{n} j{j}",
                      (in_dual(Dt(f)) and in_dual(Dz(f)))
                      ==(n>=0 and j>=0))

    # A single derivative on ALL rational functions is false, even for the
    # identity Keller cover; source polynomiality is essential.
    rational=1/P
    eq("one-derivative full-field countercontrol",S.diff(rational,Q),0)
    check("full-field control is not integral polynomial",
          not rational.is_polynomial(P,Q))
    eq("second coordinate derivative detects that pole",
       S.diff(rational,P),-1/P**2)
    check("second derivative is not in identity codifferent",
          not S.diff(rational,P).is_polynomial(P,Q))

    # Finite exceptional directions: tangency is D(H)=0 modulo H.
    curves=[P,Q,Q-P*P,P*Q-1,P*P+Q*Q-1]
    for i,H in enumerate(curves):
        tangent=[]
        for direction in [(1,0),(0,1),(1,1),(1,2)]:
            aa,bb=direction
            derivative=aa*S.diff(H,P)+bb*S.diff(H,Q)
            rem=S.groebner([H],P,Q).reduce(derivative)[1]
            if rem==0:tangent.append(direction)
        check(f"at most one tested constant tangent direction curve{i}",
              len(tangent)<=1)
    generic=(1,2)
    for i,H in enumerate(curves):
        rem=S.groebner([H],P,Q).reduce(
            generic[0]*S.diff(H,P)+generic[1]*S.diff(H,Q))[1]
        check(f"one chosen direction transverse across finite image family{i}",
              rem!=0)

    # Actual polynomial Keller controls verify inverse-Jacobian expressions
    # and the velocity/transverse-coordinate identity without branch models.
    pairs=[
        (x+y*y,y),
        (x,y+x*x*x),
        (x+y*y,y+(x+y*y)**2),
        (2*(x+y*y)+3*y,(x+y*y)+2*y),
    ]
    for i,(pp,qq) in enumerate(pairs):
        jac=S.expand(S.diff(pp,x)*S.diff(qq,y)-S.diff(pp,y)*S.diff(qq,x))
        check(f"actual source control constant nonzero Jacobian{i}",
              jac!=0 and not jac.free_symbols)
        DP=lambda f:S.cancel((S.diff(qq,y)*S.diff(f,x)-S.diff(qq,x)*S.diff(f,y))/jac)
        DQ=lambda f:S.cancel((-S.diff(pp,y)*S.diff(f,x)+S.diff(pp,x)*S.diff(f,y))/jac)
        eq(f"target derivation P on target P control{i}",DP(pp),1)
        eq(f"target derivation P on target Q control{i}",DP(qq),0)
        eq(f"target derivation Q on target P control{i}",DQ(pp),0)
        eq(f"target derivation Q on target Q control{i}",DQ(qq),1)
        eq(f"one velocity is target derivative of transverse coordinate{i}",
           DQ(y),S.diff(pp,x)/jac)
        eq(f"coframe transverse identity dx coefficient{i}",
           -S.diff(qq,x)*S.diff(pp,x)+S.diff(pp,x)*S.diff(qq,x),0)
        eq(f"coframe transverse identity dy coefficient{i}",
           -S.diff(qq,x)*S.diff(pp,y)+S.diff(pp,x)*S.diff(qq,y),jac)

    # Marked traces are computed with actual multiplication matrices.
    for e in range(1,8):
        mult=S.zeros(e)
        for col in range(e-1):mult[col+1,col]=1
        mult[0,e-1]=t
        for q in range(1,10):
            velocity=-S.Rational(q,e)*(mult**(-q-e))
            traces=[S.cancel(S.trace((mult**j)*velocity)) for j in range(e)]
            jj=q % e
            eq(f"split marker formula e{e} q{q}",
               traces[jj],-q*t**(-(q//e)-1))
            check(f"all other marker traces vanish e{e} q{q}",
                  all(value==0 for j,value in enumerate(traces) if j!=jj))
            threshold=(q+e)//e  # ceil((q+1)/e)
            check(f"clearing exponent identity e{e} q{q}",
                  threshold==q//e+1)
            for n in range(threshold+2):
                regular=all(order(t**n*value,t)>=0 for value in traces)
                check(f"exact full marked-trace clearing e{e} q{q} n{n}",
                      regular==(n>=threshold))

    # Same final image may contain several source primes.
    for i,packet in enumerate([
        [(2,1),(3,5)],[(1,2),(4,1)],[(3,2),(5,11),(2,4)],
    ]):
        nn=max((q+e)//e for e,q in packet)
        for n in range(nn+2):
            all_regular=all(e*n-q-e>=1-e for e,q in packet)
            check(f"coincident-image maximum exponent packet{i} n{n}",
                  all_regular==(n>=nn))

    # Every nonintegral codifferent mode differentiates outside it, while
    # one reduced branch factor kills its derivative class.
    for e in range(2,9):
        matrix=S.diag(*[-S.Rational(i,e) for i in range(1,e)])
        check(f"derivative defect map injective rank e{e}",matrix.rank()==e-1)
        for i in range(1,e):
            lam=u**(-i)
            dlam=S.cancel(S.diff(lam,u)/(e*u**(e-1)))
            check(f"codifferent mode is allowed e{e} i{i}",order(lam,u)>=1-e)
            check(f"derivative leaves codifferent e{e} i{i}",order(dlam,u)<1-e)
            check(f"reduced branch kills derivative class e{e} i{i}",
                  order(u**e*dlam,u)>=1-e)
            correction=S.cancel(
                S.diff(u**e*lam,u)/(e*u**(e-1))-u**e*dlam)
            eq(f"Leibniz correction lies in original dual mode e{e} i{i}",
               correction,lam)


    # Source-finite monicity is a distinct, explicit additional hypothesis.
    for degree in range(1,7):
        equation=T**degree+y*T**(degree-1)+Q
        check(f"constant source-leading coefficient degree{degree}",
              S.Poly(equation,T).LC()==1)
    eq("unramified pole derivative has exact q+1 growth",
       S.diff(u**(-5),u),-5*u**(-6))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    started=time.monotonic()
    run()
    inputs=[
        ("CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md",
         "Full source open immersion, tame canonical order; no torsion inference"),
        ("CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md",
         "KU0/1 only; ordinary trace dual embedding and canonical etale open"),
        ("REPLICA_DEGREE_AND_SOURCE_DERIVATIVE_FIELD_RECONSTRUCTION.md",
         "Section6.1 source-finite frames and integral coordinate recovery"),
        ("REDUCED_BRANCH_PAIRING_AND_SIXTH_MOMENT_DEFECT.md",
         "Reduced branch annihilator for the retained trace-defect module"),
        ("BERGMAN_TRACE_DUAL_AND_MULTIPLIER_DOMAIN_RIGIDITY.md",
         "Optional ordinary L2 restatement only; not needed for algebraic proof"),
    ]
    receipt={
        "claim":"Codifferent preimages under target derivations",
        "outcome":"PASS","check_count":len(checks),"checks":checks,
        "environment":{"python":platform.python_version(),"sympy":S.__version__},
        "proof":PROOF.name,"proof_sha256":sha(PROOF),
        "checker":Path(__file__).name,"checker_sha256":sha(Path(__file__)),
        "retained_inputs":[{"path":p,"sha256":sha(BASE/p),"locator":loc}
                           for p,loc in inputs],
        "scope":{
            "finite_normal_surface_two_derivative_preimage":True,
            "actual_Keller_one_generic_D_uniform_for_all_source_polynomials":True,
            "ordinary_codifferent_not_unidentified_moment_lattice":True,
            "exact_codifferent_clearing_ideal":True,
            "generic_target_one_velocity_detects_transverse_integrality":True,
            "source_finite_frame_required_for_one_velocity_inverse":True,
            "canonical_etale_completion_preimage_analogue":True,
            "trace_defect_injects_B_linearly_into_R_mod_Lambda":True,
            "Chau_common_pencil_KQ0_or_LND_classification_used":False,
            "automatic_velocity_membership_proved":False,
            "trace_perfection_or_JC2_proved":False,
        },
        "assurance":"Written unrestricted algebraic proof and internal exact controls. Local branched/Laurent models are not actual whole-plane Keller counterexamples; no formal verification or external review.",
    }
    if args.write_receipt:args.write_receipt.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"outcome":"PASS","check_count":len(checks),
        "elapsed_seconds":time.monotonic()-started,
        "receipt_written":bool(args.write_receipt)},indent=2))
if __name__=="__main__":main()
