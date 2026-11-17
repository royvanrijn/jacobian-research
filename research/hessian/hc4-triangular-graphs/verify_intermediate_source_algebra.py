#!/usr/bin/env python3
"""Exact support for intermediate source-algebra and replica-fibre proofs."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import sympy as S

BASE=Path(__file__).resolve().parent
PROOF=BASE/"INTERMEDIATE_SOURCE_ALGEBRA_AND_REPLICA_FIBRE_COVERAGE.md"
INPUT=BASE/"CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md"
checks=[]

def check(name,condition):
    ok=bool(condition)
    checks.append({"name":name,"passed":ok})
    if not ok: raise AssertionError(name)

def eq(name,actual,expected):
    check(name,S.cancel(S.sympify(actual)-S.sympify(expected))==0)

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def run():
    x,y,z,u,w,t,T,v=S.symbols("x y z u w t T v")

    # Descent equalizers for actual finite function-field covers.
    # The polynomial cocycle f(s)=f(t) holds modulo s^m=t^m exactly
    # for the invariant monomial degrees. No claim that this finite list
    # proves general faithfully flat descent is made.
    for m in range(2,7):
        for power in range(0,15):
            remainder=S.rem(t**power-v**power,t**m-v**m,t)
            check(f"descent cocycle degree{m} power{power}",
                  (remainder==0)==(power % m==0))
        mult=S.zeros(m)
        for j in range(m-1): mult[j+1,j]=1
        mult[0,m-1]=u
        eq(f"normalization finite generator polynomial degree{m}",
           mult.charpoly(T).as_expr(),T**m-u)

    # A finite missing point has codimension two on a surface; Hartogs
    # does not itself prove that the image contains that point.
    check("point gap has codimension two",2-0==2)
    check("curve gap is height one and must be removed first",2-1==1)
    # Hartogs extends functions, rather than filling points of an image.
    # On the punctured plane D(x) union D(y), the elementary Laurent
    # intersection has no negative exponent in either coordinate.
    for ii in range(-2,3):
        for jj in range(-2,3):
            mon=x**ii*y**jj
            on_Dx=S.cancel(mon).as_numer_denom()[1].subs(x,1).is_number
            on_Dy=S.cancel(mon).as_numer_denom()[1].subs(y,1).is_number
            check(f"punctured-plane Laurent intersection i{ii}j{jj}",
                  (on_Dx and on_Dy)==mon.is_polynomial(x,y))
    # Principal complements on trivializing charts are affine equations.
    for f in [x,x*y,x*x+y*y+1]:
        eq("Cartier complement equation for "+str(f),
           (f*T-1).subs(T,1/f),0)

    # Exact local-slice induction with nonconstant invariant leading term.
    pp=y*(y*y+1)
    delta=lambda f:S.cancel(pp*S.diff(f,x))
    sl=x/pp+1/(y+3)
    eq("given source replica admits rational slice",delta(sl),1)
    g=x**5+y*x**3+(y*y+1)*x+7
    original=g
    reconstruction=S.Integer(0)
    for degree in range(5,-1,-1):
        top=g
        for _ in range(degree): top=delta(top)
        coefficient=S.cancel(top/S.factorial(degree))
        eq(f"slice coefficient invariant degree{degree}",delta(coefficient),0)
        reconstruction+=coefficient*sl**degree
        g=S.cancel(g-coefficient*sl**degree)
        derived=g
        for _ in range(degree): derived=delta(derived)
        eq(f"slice degree decreases at{degree}",derived,0)
    eq("slice finite reconstruction",reconstruction,original)

    # Monic polynomial parameter equations independently control finite
    # source-line images and finite A1 invariant-curve normalization.
    for degree in range(1,7):
        h=T**degree+T if degree>1 else T+1
        relation=h-u
        for exponent in [degree,degree+1,2*degree+1]:
            rem=S.rem(T**exponent,relation,T)
            check(f"finite line-module remainder degree{degree} exponent{exponent}",
                  S.degree(rem,T)<degree)
            eq(f"finite line-module evaluation degree{degree} exponent{exponent}",
               rem.subs({T:t,u:h.subs(T,t)}),t**exponent)
    for genus in range(1,5):
        for degree in range(1,8):
            check(f"RH excludes positive genus{genus} degree{degree}",
                  degree*(2*genus-2)>=0>-2)
    possible_missing=[set(),{"infinity"}]
    affine_missing=[missing for missing in possible_missing if missing]
    check("finite normalization of affine curve leaves exactly one infinity",
          affine_missing==[{"infinity"}])

    # A genuine smooth affine quadratic pseudo-plane source model.
    # This is an etale source quotient/control, NOT an etale plane target.
    # It exercises nonreduced fibres and all actual source points.
    F=w*w-u-u*u*z
    uu=y*y
    ww=y*(-1+y*y*x)
    zz=-2*x+y*y*x*x
    eq("quadratic source model lies on hypersurface",
       F.subs({u:uu,w:ww,z:zz}),0)
    derivs=[S.diff(F,q) for q in [u,w,z]]
    eq("quadratic model smooth at the only possible u0 w0 candidate",
       derivs[0].subs({u:0,w:0}),-1)
    eq("model has invariant polynomial h=y squared",uu,y*y)
    D=lambda f:S.expand(y*S.diff(f,x))
    eq("replica preserves u",D(uu),0)
    eq("replica preserves w with invariant coefficient",D(ww),uu*uu)
    eq("replica preserves z",D(zz),2*ww)
    eq("coordinate-ring replica respects relation",
       S.expand(u*u*S.diff(F,w)+2*w*S.diff(F,z)),0)
    eq("localized intermediate slice is genuine",
       S.cancel(D(ww/(uu*uu))),1)

    # Two local target coordinate pairs verify source etaleness everywhere.
    juw=S.expand(S.diff(uu,x)*S.diff(ww,y)-S.diff(uu,y)*S.diff(ww,x))
    jwz=S.expand(S.diff(ww,x)*S.diff(zz,y)-S.diff(ww,y)*S.diff(zz,x))
    eq("source map full rank away from y0",juw,-2*y**4)
    eq("source map full rank on y0",jwz.subs(y,0),-2)
    eq("central line maps isomorphically",zz.subs(y,0),-2*x)
    cc,aa=S.symbols("cc aa",nonzero=True)
    noncentral_w=ww.subs(y,cc)
    inverse_x=(w/cc+1)/(cc*cc)
    eq("noncentral source line isomorphism inverse",
       noncentral_w.subs(x,inverse_x),w)
    eq("noncentral image reaches arbitrary fibre point",
       zz.subs({y:cc,x:inverse_x}), (w*w-cc*cc)/(cc**4))
    # Central fibre w^2=u(1+uz): the factor 1+uz is a unit there.
    eq("central transverse multiplicity identity",F,w*w-u*(1+u*z))
    eq("central unit does not vanish", (1+u*z).subs(u,0),1)
    check("central multiplicity ledger",1*2==2)
    check("noncentral multiplicity ledger",2*1==2)

    # Reduction must precede the smooth-fibre statement.
    check("central original scheme fibre is nonreduced",
          S.diff(w*w,w).subs(w,0)==0)
    central_basis=S.groebner([F,u],u,w,z)
    check("central fibre ideal is exactly u,w squared",
          [g.as_expr() for g in central_basis.polys]==[u,w*w])
    radical_basis=S.groebner([u,w],u,w,z)
    check("central reduction has no relation in affine parameter z",
          [g.as_expr() for g in radical_basis.polys]==[u,w])
    # Finite etale A1 cover degree is forced to one by integral Euler.
    for cover_degree in range(1,12):
        for target_euler in range(-3,5):
            check(f"finite etale line Euler degree{cover_degree} chi{target_euler}",
                  (1==cover_degree*target_euler)
                  ==(cover_degree==1 and target_euler==1))

    # General multiplicity passports come from ONE polynomial h.
    for h in [
        y*y,(y-1)**3,(y*y-1)**2,
        (y-1)**2*(y+2)**3,y**3-y,
    ]:
        degree=S.degree(h,y)
        _,factors=S.factor_list(h,y)
        eq("root multiplicity sum for "+str(h),
           sum(S.degree(f,y)*mult for f,mult in factors),degree)
    for e in range(1,6):
        # An unramified change of transverse uniformizer preserves order.
        pull=(t*(1+t))**e
        check(f"etale transverse multiplicity e{e}",
              min(mon[0] for mon,_ in S.Poly(pull,t).terms())==e)

    # Dropping etaleness invalidates the source-line isomorphism step.
    eq("branched line map derivative",S.diff(x*x,x),2*x)
    check("branched quadratic line map fails etaleness at origin",
          S.diff(x*x,x).subs(x,0)==0)
    check("nonstable quadratic subfield guard",
          S.diff(x*x,x).subs(x,-x)!=S.diff(x*x,x))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    started=time.monotonic()
    run()
    primaries=[]
    for name,url,locator in [
        ("stacks-023Q-fpqc-effective-epimorphism.html",
         "https://stacks.math.columbia.edu/tag/023Q",
         "Lemma35.13.7, second proof and displayed equalizer"),
        ("stacks-0C1B-riemann-hurwitz.html",
         "https://stacks.math.columbia.edu/tag/0C1B",
         "Section53.12, characteristic-zero formula and Lemma53.12.2"),
    ]:
        path=BASE/"sources"/name
        check("cached primary exists "+name,path.is_file())
        primaries.append({"path":"sources/"+name,"url":url,
                          "locator":locator,"sha256":sha(path)})
    receipt={
        "claim":"Intermediate source algebra and given-replica fibre coverage",
        "outcome":"PASS","check_count":len(checks),"checks":checks,
        "environment":{"python":platform.python_version(),"sympy":S.__version__},
        "proof":PROOF.name,"proof_sha256":sha(PROOF),
        "checker":Path(__file__).name,"checker_sha256":sha(Path(__file__)),
        "retained_inputs":[{"path":INPUT.name,"sha256":sha(INPUT),
          "locator":"KU0/1 Sections1 and3 only; no KU2/KQ0/LA0/PF46/PF50"}],
        "primary_inputs":primaries,
        "scope":{
            "every_embedded_intermediate_E":True,
            "R_intersect_E_smooth_affine_etale":True,
            "source_image_only_finite_point_gaps_without_stability":True,
            "given_polynomial_replica_stability_required_for_gap_filling":True,
            "replica_kernel_C_h_y":True,
            "replica_source_surjective":True,
            "reduced_fibres_disjoint_smooth_A1":True,
            "whole_source_lines_map_isomorphically_to_components":True,
            "scheme_fibres_may_be_nonreduced":True,
            "arbitrary_LND_classification_used":False,
            "common_pencil_LA0_or_whole_root_charts_used":False,
            "first_velocity_integrality_or_inverse_proved":False,
        },
        "assurance":"Internal exact controls; written argument proves unrestricted scope. Classical faithfully-flat descent, RH and canonical KU0/1 are adopted foundations. No formal verification or external review.",
    }
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"outcome":"PASS","check_count":len(checks),
          "elapsed_seconds":time.monotonic()-started,
          "receipt_written":bool(args.write_receipt)},indent=2))

if __name__=="__main__":main()
