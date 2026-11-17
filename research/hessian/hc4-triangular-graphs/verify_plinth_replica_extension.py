#!/usr/bin/env python3
"""Exact controls for plinth ideals and full-normalization replica extension."""
from __future__ import annotations
import argparse, hashlib, json, platform, time
from pathlib import Path
import sympy as S

BASE=Path(__file__).resolve().parent
PROOF=BASE/"PLINTH_IDEALS_AND_INVARIANT_REPLICA_EXTENSION.md"
checks=[]
def check(name,value):
    ok=bool(value);checks.append({"name":name,"passed":ok})
    if not ok:raise AssertionError(name)
def eq(name,a,b=0):check(name,S.cancel(S.sympify(a)-S.sympify(b))==0)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def order(f,z):
    f=S.cancel(f)
    if f==0:return S.oo
    n,d=S.fraction(f)
    return (min(m[0] for m,c in S.Poly(n,z).terms())
            -min(m[0] for m,c in S.Poly(d,z).terms()))

def run():
    x,y,z,t,u,T,lam=S.symbols("x y z t u T lam")
    # All these are local whole-chart controls, not finite plane targets.
    for r in range(1,6):
        for p in range(r+1,r+6):
            threshold=(p+r-1)//r
            xx=z**(1-r)*t+z**(-p)
            yy=lam+z**r
            eq(f"unit root-chart volume r{r}p{p}",
               S.diff(xx,t)*S.diff(yy,z),r)
            check(f"strict primitive forces at least double plinth r{r}p{p}",
                  threshold>=2)
            for n in range(0,threshold+2):
                q=(y-lam)**n
                image=S.cancel(q.subs(y,yy)*xx)
                check(f"exact plinth threshold r{r}p{p}n{n}",
                      (order(image,z)>=0)==(n>=threshold))
                if n<threshold:
                    eq(f"regular invariant cannot cancel a pole r{r}p{p}n{n}",
                       order(image+yy**3+1,z),r*n-p)
    # Group all boundary components with the same source fibre.
    packets=[
        [(0,2,1),(0,5,2),(3,4,2)],
        [(1,3,2),(-2,7,3)],
        [(0,6,1),(0,8,3),(2,5,4)],
    ]
    for k,packet in enumerate(packets):
        exponents={}
        for ll,p,r in packet:
            exponents[ll]=max(exponents.get(ll,0),(p+r-1)//r)
        g=S.prod((y-ll)**nn for ll,nn in exponents.items())
        for ll,p,r in packet:
            value=order(g.subs(y,ll+z**r),z)
            check(f"same-fibre max clears every component packet{k}p{p}r{r}",
                  value>=p)
        for ll,nn in exponents.items():
            smaller=S.cancel(g/(y-ll))
            check(f"same-fibre exponent is minimal packet{k}lambda{ll}",
                  any(order(smaller.subs(y,aa+z**rr),z)<pp
                      for aa,pp,rr in packet if aa==ll))

    # Genuine smooth affine plane-completion controls M_p, with no actual
    # finite etale plane target supplied. Their translation is complete.
    for p in range(2,7):
        b=x*y**p
        a=x*(b-1)
        functions=[b,a,b*b-y**p*a,a*a+b, a*b+b**3]
        for j,F in enumerate(functions):
            deg=int(S.degree(F,x))
            lc=S.Poly(F,x).LC()
            eq(f"leading coefficient is last Taylor derivative p{p}j{j}",
               S.diff(F,x,deg),S.factorial(deg)*lc)
            check(f"leading coefficient divisible by exact plinth p{p}j{j}",
                  S.rem(lc,y**p,y)==0)
        eq(f"plinth generator has degree-one antiderivative p{p}",
           S.diff(b,x),y**p)
        eq(f"nontrivial leading cancellation remains in plinth p{p}",
           b*b-y**p*a,b)
        aa,bb,s=S.symbols("aa bb s")
        rel=y**p*aa-bb*(bb-1)
        flowed=rel.subs({bb:bb+s*y**p,
                        aa:aa+s*(2*bb-1)+s*s*y**p},simultaneous=True)
        eq(f"complete source translation preserves affine relation p{p}",flowed,rel)

    # An invariant source multiplier cannot clear a polar coordinate.
    for e in range(2,7):
        for pole in range(1,6):
            yy=u**(-pole)*(1+u*z)
            velocity=u**(-pole-e)*(z+1)
            for degree in range(0,5):
                q=(S.Integer(1) if degree==0
                   else T+1 if degree==1 else T**degree+T+1)
                multiplied=S.cancel(q.subs(T,yy)*velocity)
                eq(f"invariant multiplier worsens exact pole e{e}q{pole}d{degree}",
                   order(multiplied,u),-(degree+1)*pole-e)
                check(f"no invariant multiplier gives ordinary dual landing e{e}q{pole}d{degree}",
                      order(multiplied,u)<1-e)

    # Denominator clearing when the invariant coordinate IS regular:
    # the retained regular-coordinate completion n=2,q=2.
    zz=x*y*y
    a1=x*(zz-1)
    a2=x*x*(zz-1)
    eq("regular-coordinate clearing on z",y*y*S.diff(zz,x),y**4)
    eq("regular-coordinate clearing on a1",y*y*S.diff(a1,x),y*y*(2*zz-1))
    eq("regular-coordinate clearing on a2",y*y*S.diff(a2,x),3*y*y*a1+zz)
    D=lambda f:S.expand(y*y*S.diff(f,x))
    for generator in [y,zz,a1,a2]:
        ff=generator
        for _ in range(4):ff=D(ff)
        eq("cleared invariant replica is nilpotent on generator "+str(generator),ff,0)

    # The Lie/commutator identities used for the ordinary canonical module.
    q=y*y+y+1
    ss=x*x+y
    ll=x**3+x*y+y*y
    eq("multiplier-module commutator",
       q*S.diff(ss*ll,x)-ss*q*S.diff(ll,x),q*S.diff(ss,x)*ll)
    eq("invariant replica preserves source volume divergence",S.diff(q,x),0)
    for n in range(1,5):
        base=x**4+x*y+y*y
        plain=S.diff(base,x,n)
        ff=base
        for _ in range(n):ff=q*S.diff(ff,x)
        eq(f"invariant replica finite Taylor identity n{n}",ff,q**n*plain)

    # Derived Taylor obstruction from the linked noncontracting affine N.
    # M below is its canonical coefficient module, NOT claimed trace-dual
    # for any finite Keller plane target.
    xx=u*u+u**3*t
    yy=1/u
    gens=[x,x*y,x*y*y-1,x*y**3-y]
    for j,f in enumerate(gens):
        check(f"linked N generator regular at entire generic E j{j}",
              order(f.subs({x:xx,y:yy}),u)>=0)
    delta_b=S.diff(x*y,x)
    delta_c=S.diff(x*y*y-1,x)
    eq("Taylor invariant from regular b is y",delta_b,y)
    eq("Taylor derivative from regular c is y squared",delta_c,y*y)
    eq("b Taylor chain terminates",S.diff(x*y,x,2),0)
    eq("c Taylor chain terminates",S.diff(x*y*y-1,x,2),0)
    check("invariant y is in canonical M but not regular N",
          order(delta_b.subs({x:xx,y:yy}),u)==-1)
    check("derivative y squared leaves canonical M",
          order(delta_c.subs({x:xx,y:yy}),u)==-2)
    eq("canonical local generator y has exactly simple pole",u*yy,1)
    jac=S.diff(xx,u)*S.diff(yy,t)-S.diff(xx,t)*S.diff(yy,u)
    eq("linked control volume has one boundary zero",jac,u)
    theta_u=S.cancel(xx*S.diff(yy,u))
    eq("linked global primitive is regular on other chart",theta_u,-(1+u*t))
    eq("linked primitive differentiates to volume",-S.diff(theta_u,t),u)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    started=time.monotonic();run()
    inputs=[
        ("EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md",
         "Common actual source frame and strict p_i>r_i boundary inequality"),
        ("COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md",
         "Whole root charts and same-volume primitive"),
        ("HAMILTONIAN_COMMON_COORDINATE_AND_ADDITIVE_COMPLETION.md",
         "Actual complete source translation and kernel"),
        ("CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md",
         "Generic-target first-velocity pole formula; ordinary trace dual only"),
        ("CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md",
         "Full finite source open and canonical coefficient identity"),
        ("RAMIFIED_FOLD_VALUATIONS_AND_COMMON_COORDINATE_POLES.md",
         "Linked affine noncontracting survivor control only"),
        ("REGULAR_SOURCE_COORDINATE_COMPLETIONS_AND_PRIMITIVE_RESIDUES.md",
         "Regular-coordinate denominator-clearing control"),
    ]
    receipt={
        "claim":"Plinth ideals and invariant replica extension",
        "outcome":"PASS","check_count":len(checks),"checks":checks,
        "environment":{"python":platform.python_version(),"sympy":S.__version__},
        "proof":PROOF.name,"proof_sha256":sha(PROOF),
        "checker":Path(__file__).name,"checker_sha256":sha(Path(__file__)),
        "retained_inputs":[{"path":p,"sha256":sha(BASE/p),"locator":loc}
                           for p,loc in inputs],
        "scope":{
            "ordinary_full_finite_S_and_canonical_etale_A_distinct":True,
            "exact_A_plinth_and_leading_coefficient_ideal":True,
            "proper_common_frame_never_source_finite":True,
            "invariant_replica_preserves_S_or_Lambda_iff_y_in_S":True,
            "given_source_coordinate_only_not_arbitrary_LND_classification":True,
            "linked_control_canonical_module_not_claimed_trace_dual":True,
            "automatic_velocity_or_coordinate_landing_proved":False,
            "HC4_JC2_inverse_proved":False,
        },
        "assurance":"Internal exact controls with written unrestricted proof; geometric controls have no supplied finite Keller plane target. No formal verification or external review.",
    }
    if args.write_receipt:args.write_receipt.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"outcome":"PASS","check_count":len(checks),
      "elapsed_seconds":time.monotonic()-started,
      "receipt_written":bool(args.write_receipt)},indent=2))
if __name__=="__main__":main()
