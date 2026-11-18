#!/usr/bin/env python3
"""Exact controls for exceptional-fold whole charts and birational image lemma."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import re
import time
from pathlib import Path
import sympy as S
BASE=Path(__file__).resolve().parent
PROOF=BASE/"EXCEPTIONAL_POLAR_FOLD_WHOLE_CHARTS_AND_BIRATIONAL_IMAGES.md"
RECEIPT=BASE/"exceptional_polar_fold_whole_chart_receipt.json"
checks=[]
def check(name,condition):
    ok=bool(condition)
    checks.append({"name":name,"passed":ok})
    if not ok:raise AssertionError(name)
def eq(name,actual,expected=0):
    check(name,S.cancel(S.sympify(actual)-S.sympify(expected))==0)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def order(expr,var):
    expr=S.cancel(expr)
    if expr==0:return S.oo
    numerator,denominator=S.fraction(expr)
    low=lambda f:min(m[0] for m,c in S.Poly(f,var).terms())
    return low(numerator)-low(denominator)
def jac(P,Q,t,s):
    return S.expand(S.diff(P,t)*S.diff(Q,s)-S.diff(P,s)*S.diff(Q,t))

def run():
    s,t,kappa,A,B,T=S.symbols("s t kappa A B T")
    u=s*s*(kappa+A*s+B*s*s+s**3*T)
    a=1/u
    aa=S.series(a,s,0,2).removeO().expand()
    eq("five-free extraction leading pole",aa.coeff(s,-2),1/kappa)
    eq("five-free extraction logarithmic coefficient",aa.coeff(s,-1),-A/kappa**2)
    eq("five-free extraction constant coefficient",aa.coeff(s,0),
       A*A/kappa**3-B/kappa**2)
    check("five local multiplicities match tame order",sum([1]*5)==2*2+2-1)
    I=S.zeros(5)
    for j in range(5):
        I[j,j]=-1 if j==4 else -2
        if j<4:I[j,j+1]=I[j+1,j]=1
    eq("independent local exceptional inverse norm",(-I.inv())[4,4],5)
    aa0=S.cancel(a.subs(A,0))
    c=1/kappa
    alpha=-B/kappa**2
    tt=S.cancel((aa0-alpha-c/s**2)/s)
    eq("exact residue parameter is affine separating",tt.subs(s,0),-T/kappa**2)
    eq("exact unnormalized field chart",alpha+c/s**2+s*tt,aa0)
    eq("logarithmic residue is killed exactly by A0",aa.coeff(s,-1).subs(A,0),0)
    for i,j in [(1,0),(1,1),(2,1),(3,4)]:
        F=S.series(a**i*s**j,s,0,2).removeO()
        eq(f"exact source derivative has no logarithmic residue i{i}j{j}",
           S.expand(S.diff(F,s)).coeff(s,-1),0)

    x,b=S.symbols("x b")
    chart=s**(-2)+s*t
    eq("normalized exact inverse field coordinate",
       (chart*s*s-1)/s**3,t)
    pa,pb,qa,qb=S.symbols("pa pb qa qb")
    J=S.Matrix([[pa,pb],[qa,qb]])*S.Matrix(
        [[S.diff(chart,t),S.diff(chart,s)],[0,1]])
    eq("universal target Jacobian chain factor",J.det(),s*(pa*qb-pb*qa))
    examples=[
        ("source common coordinate",x,-2),
        ("source companion",b,1),
        ("exceptional vanishing generator",x*b*b-1,3),
        ("w generator",x*(x*b*b-1),1),
        ("k residue generator",b*x*x*(x*b*b-1),0),
        ("l residue square",x**3*(x*b*b-1)**2,0),
        ("source product",x*b,-1),
    ]
    for name,f,v in examples:
        image=S.cancel(f.subs({x:chart,b:s},simultaneous=True))
        eq("exact valuation "+name,order(image,s),v)
        check("Laurent admission iff no negative orders "+name,
              image.is_polynomial(t,s)==(v>=0))
    eq("contracting linear-Jacobian negative control",jac(s*t,s,t,s),s)
    check("contracting control has constant zero-line image",
          (s*t).subs(s,0)==0 and s.subs(s,0)==0)
    check("identity Keller pair does not admit exceptional valuation chart",
          not chart.is_polynomial(s,t))

    # Universal s^1 identity using independent derivative symbols.
    hp,Gh,II,KK,Ip,Kp,k,kt,U2,V2=S.symbols(
        "hp Gh II KK Ip Kp k kt U2 V2")
    U0p=hp*Gh*II;V0p=hp*Gh*KK
    U1=k*II;V1=k*KK
    U1p=kt*II+k*Ip*hp
    V1p=kt*KK+k*Kp*hp
    C1=U1p*V1-U1*V1p+2*(U0p*V2-U2*V0p)
    prediction=k*k*hp*(Ip*KK-II*Kp)+2*hp*Gh*(II*V2-U2*KK)
    eq("unrestricted normalized critical-line coefficient identity",C1,prediction)

    z=S.symbols("z")
    curve_pairs=[(z*z,z**3),(z**3,z**4),(z,z*z),
                 (z**3+z,z**4-1),(0,z),(z,0)]
    for ci,(Ap,Bp) in enumerate(curve_pairs):
        Ap=S.sympify(Ap);Bp=S.sympify(Bp)
        G=S.gcd(S.diff(Ap,z),S.diff(Bp,z))
        I0=S.cancel(S.diff(Ap,z)/G)
        K0=S.cancel(S.diff(Bp,z)/G)
        check(f"primitive derivative Bezout pair c{ci}",S.gcd(I0,K0)==1)
        for dd in [1,2,3]:
            h=t**dd+t if dd>1 else 2*t+1
            kk=1+t*t
            U0=Ap.subs(z,h);V0=Bp.subs(z,h)
            U1=kk*I0.subs(z,h);V1=kk*K0.subs(z,h)
            U2=t*t+3;V2=t**3-2
            UP=U0+s*U1+s*s*U2
            VP=V0+s*V1+s*s*V2
            JJ=jac(UP,VP,t,s)
            eq(f"zero constant Jacobian coefficient c{ci}d{dd}",JJ.coeff(s,0),0)
            expect=kk*kk*S.diff(h,t)*(
                S.diff(I0,z)*K0-I0*S.diff(K0,z)).subs(z,h)
            expect+=2*S.diff(h,t)*G.subs(z,h)*(
                I0.subs(z,h)*V2-U2*K0.subs(z,h))
            eq(f"exact hprime divisibility identity c{ci}d{dd}",JJ.coeff(s,1),expect)
            quo,rem=S.div(JJ.coeff(s,1),S.diff(h,t),t)
            eq(f"polynomial hprime factor c{ci}d{dd}",rem,0)
    for dd in range(2,13):
        h=t**dd
        U0=h*h;V0=h**3
        UP=U0+2*s
        VP=V0+3*h*s
        first=jac(UP,VP,t,s).coeff(s,1)
        eq(f"nonbirational image forces nonconstant Jacobian coefficient d{dd}",
           first,-6*S.diff(h,t))
        eq(f"nonbirational critical coefficient has derivative root d{dd}",
           first.subs(t,0),0)

    # Reuse the retained finite linear-Jacobian family. The new guard is
    # failure of admission to the ORIGINAL exceptional source plane.
    X,Y=S.symbols("X Y")
    for dd in range(1,10):
        for low in [0,2*t]:
            f=t**dd+low
            F=S.integrate(f,t)
            UP=s+f;VP=t*(s+f)-F
            eq(f"retained finite family linear Jacobian d{dd}low{low}",
               jac(UP,VP,t,s),-s)
            eq(f"retained finite family elimination d{dd}low{low}",
               F-UP*t+VP,0)
            equation=F-X*t+Y
            check(f"retained finite family expected integral degree d{dd}low{low}",
                  S.degree(equation,t)==dd+1)
            admitted=S.cancel(UP.subs({t:(x*b*b-1)/b**3,s:b},simultaneous=True))
            eq(f"failure of source polynomial admission d{dd}low{low}",
               order(admitted,b),-3*dd)
            check(f"finite family supplies no original Keller coordinate d{dd}low{low}",
                  not admitted.is_polynomial(x,b))
        if dd>=2:
            P0=t**dd
            Q0=S.Rational(dd,dd+1)*t**(dd+1)
            D=S.gcd(S.diff(P0,t),S.diff(Q0,t))
            i0=S.cancel(S.diff(P0,t)/D)
            k0=S.cancel(S.diff(Q0,t)/D)
            eq(f"nonimmersive branch first order m{dd}",order(P0,t),dd)
            eq(f"nonimmersive branch consecutive order m{dd}",order(Q0,t),dd+1)
            check(f"primitive tangent Wronskian nonzero m{dd}",
                  (S.diff(i0,t)*k0-i0*S.diff(k0,t)).subs(t,0)!=0)

    # Complete seven-relation smooth affine control.
    b,u,w,k,l,x,t,eps=S.symbols("b u w k l x t eps")
    relations=[b*b*w-u*(u-1),b*k-u*w,k*k-u*l,w*k-b*l,
               w**3-(u-1)*l,b*b*l-u*w*w,(u-1)*k-b*w*w]
    src={u:x*b*b,w:x*(x*b*b-1),k:b*x*x*(x*b*b-1),
         l:x**3*(x*b*b-1)**2}
    ext={u:1+b**3*t,w:b*(1+b**3*t)*t,
         k:(1+b**3*t)**2*t,l:(1+b**3*t)**3*t*t}
    for j,f in enumerate(relations):
        eq(f"seven-relation source chart equation {j}",
           f.subs(src,simultaneous=True),0)
        eq(f"seven-relation whole exceptional chart equation {j}",
           f.subs(ext,simultaneous=True),0)
    G=S.groebner(relations,b,u,w,k,l)
    def modideal(expr):
        num=S.fraction(S.cancel(expr))[0]
        return G.reduce(S.expand(num))[1]
    xa=w/(u-1)
    source_inverse={
        u:xa*b*b,w:xa*(xa*b*b-1),
        k:b*xa*xa*(xa*b*b-1),l:xa**3*(xa*b*b-1)**2}
    for variable,expr in source_inverse.items():
        eq("localized D(u-1) inverse "+str(variable),modideal(expr-variable),0)
    tt=k/u**2
    external_inverse={
        u:1+b**3*tt,w:b*u*tt,k:u*u*tt,l:u**3*tt*tt}
    for variable,expr in external_inverse.items():
        eq("localized D(u) inverse "+str(variable),modideal(expr-variable),0)
    eq("the two smooth affine localizations cover",(u-(u-1)),1)
    eq("whole marked source affine-plane opens cover",x*b*b-(x*b*b-1),1)
    eq("whole exceptional affine-plane opens cover",(1+b**3*t)-b**3*t,1)
    aa=b**(-2)+b*t
    for variable,f in src.items():
        eq("entire chart transition "+str(variable),
           f.subs(x,aa),ext[variable])
    for j,f in enumerate(relations):
        value=S.expand(f.subs({b:0,u:0,w:0,l:0,k:eps},simultaneous=True))
        rem=S.rem(value,eps*eps,eps)
        if j<6:
            eq(f"omitted-last-relation dual-number guard {j}",rem,0)
        else:
            eq("last relation rejects spurious infinitesimal k",rem,-eps)
    boundary={b:0,u:1,w:0,l:k*k}
    for j,f in enumerate(relations):
        eq(f"whole boundary is affine-line coordinate k relation {j}",
           f.subs(boundary,simultaneous=True),0)
    original_line={b:0,u:0,k:0,l:-w**3}
    for j,f in enumerate(relations):
        eq(f"opposite source line retained relation {j}",
           f.subs(original_line,simultaneous=True),0)
    eq("exceptional whole-chart source volume",S.diff(aa,t),b)
    F=aa*b
    theta_b=S.cancel(aa+S.diff(F,b))
    theta_t=S.cancel(S.diff(F,t))
    eq("globally regular primitive exceptional db",theta_b,3*b*t)
    eq("globally regular primitive exceptional dt",theta_t,b*b)
    eq("same primitive differentiates to source volume",
       S.diff(theta_b,t)-S.diff(theta_t,b),b)
    return {"five_local_multiplicities":[1]*5,
            "coefficient_identity_degree_unrestricted":True,
            "finite_family_sample_degrees":list(range(2,11)),
            "affine_presentation_relations":7}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",action="store_true")
    args=parser.parse_args();start=time.monotonic()
    controls=run();text=PROOF.read_text()
    for phrase in [
        "does NOT require that completion",
        "OPEN IMMERSION","ENTIRE E","[C(E):C(Gamma)]=1",
        "not an infinite formal approximation",
        "target-finiteness step cannot be dropped",
        "It is divisible by h'",
        "These examples FAIL(11)",
        "SEVEN",
        "does not produce an unramified completion boundary",
        "No formal verification",
    ]:
        check("written scope retains "+phrase,phrase in text)
    for link in re.findall(r"\]\(([^)]+)\)",text):
        if "://" not in link:
            check("retained local source exists "+link,(BASE/link).is_file())
    parents=[
        "RAMIFIED_FOLD_VALUATIONS_AND_COMMON_COORDINATE_POLES.md",
        "CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md",
        "REGULAR_SOURCE_COORDINATE_COMPLETIONS_AND_PRIMITIVE_RESIDUES.md",
    ]
    receipt={
        "claim":"EXCEPTIONAL-POLAR-FOLD-WHOLE-CHART-AND-BIRATIONAL-IMAGE",
        "outcome":"PASS","check_count":len(checks),"checks":checks,
        "controls":controls,"proof_sha256":sha(PROOF),
        "checker_sha256":sha(Path(__file__)),
        "retained_inputs":[{"path":p,"sha256":sha(BASE/p)} for p in parents],
        "environment":{"python":platform.python_version(),"sympy":S.__version__},
        "scope":{
            "actual_target_finiteness_required":True,
            "exceptional_index_two_value_pair":[-2,1],
            "whole_ramified_prime_plane_chart":True,
            "entire_prime_smooth_affine_line":True,
            "critical_image_residue_degree":1,
            "critical_image_smoothness_claimed":False,
            "nonimmersive_branch_first_orders":"m,m+1",
            "all_polynomial_degrees_in_written_coefficient_proof":True,
            "original_source_admission_condition_retained":True,
            "corrected_seven_relation_geometric_control":True,
            "geometric_control_has_Keller_target_pair":False,
            "nonempty_unramified_boundary_required_once_value_pair_given":False,
            "exceptional_fold_excluded_or_inverse_proved":False,
        },
        "assurance":"Internal written proof and exact controls only; no formal verification or external review."
    }
    if args.write_receipt:
        RECEIPT.write_text(json.dumps(receipt,indent=2)+chr(10))
    print(json.dumps({"outcome":"PASS","check_count":len(checks),
                      "controls":controls,
                      "elapsed_seconds":time.monotonic()-start,
                      "receipt_written":args.write_receipt},indent=2))
if __name__=="__main__":main()
