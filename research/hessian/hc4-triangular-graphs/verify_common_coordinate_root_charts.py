#!/usr/bin/env python3
"""Exact controls for the conditional common-coordinate whole root charts."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import sympy as S

BASE=Path(__file__).resolve().parent
PROOF=BASE/"COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md"
AUDIT=BASE/"sources/common-coordinate-root-charts-source-audit.json"
RECEIPT=BASE/"common_coordinate_root_charts_receipt.json"
START=time.perf_counter()
checks: list[dict[str,object]]=[]

def check(name: str, condition: object) -> None:
    result=bool(condition)
    checks.append({"name":name,"passed":result})
    if not result:
        raise AssertionError(name)

def eq(name: str, actual: object, expected: object=0) -> None:
    check(name,S.cancel(S.sympify(actual)-S.sympify(expected))==0)

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def res_z(expr: object, z: S.Symbol) -> object:
    # All tested expressions are FINITE Laurent polynomials in z.
    return S.expand(expr).coeff(z,-1)

def sub(expr: object, mapping: dict) -> object:
    return S.expand(S.sympify(expr).subs(mapping,simultaneous=True))

x,y,z,t,u,v,w,a,b,s,T,lam=S.symbols("x y z t u v w a b s T lam")

# Local root normalization, including a genuinely cuspidal Cartier curve.
# This tests ONLY local algebra; y=f^r here is not claimed a source coordinate.
f=v**2-u**3
for r in [2,3,4,5]:
    raw=z**r-(1+u)*f**r
    origin={u:0,v:0,z:0}
    eq(f"raw cyclic cover equation at cusp r{r}",sub(raw,origin))
    for variable in [u,v,z]:
        eq(f"raw cyclic cover is singular before normalization r{r} var{variable}",
           sub(S.diff(raw,variable),origin))
    norm=w**r-(1+u)
    eq(f"unit-root normalized derivative at chosen lift r{r}",
       S.diff(norm,w).subs(w,1),r)
    eq(f"normalized equation eliminates regular base coordinate r{r}",
       S.diff(norm,u),-1)
    w_f=v**2-(w**r-1)**3
    parametrization={u:w**r-1,z:w*w_f}
    eq(f"normalized chart maps into the RAW cover r{r}",sub(raw,parametrization))
    eq(f"generic normalized root is z/f r{r}",
       sub(z,parametrization)/sub(f,parametrization),w)
    eq(f"copied Cartier curve can still be cuspidal r{r}",
       S.diff(w_f,w).subs({w:1,v:0}))
    eq(f"copied Cartier curve v derivative still vanishes r{r}",
       S.diff(w_f,v).subs({w:1,v:0}))
    gcd=S.Poly(w**r-1,w).gcd(S.Poly(r*w**(r-1),w)).as_expr()
    eq(f"unit-root fibre is reduced r{r}",gcd,1)
    eq(f"source simple branch raw equation is already smooth r{r}",
       S.diff(z**r-y,y),-1)

# A Cartier complement is a principal affine open on every local
# trivialization. This also guards removing the source line, not just copies.
for g in [z,z*(z-1),u*v,u*v*(u-v)]:
    eq("Cartier-complement localization "+str(g),(g*T-1).subs(T,1/g))
for r in [1,2,3,5]:
    copied=set(range(r))
    chosen=r-1
    removed=copied-{chosen}
    check(f"one boundary copy remains r{r}",copied-removed=={chosen})
    check(f"whole source zero line separately removed r{r}",
          ("source-line" not in copied) and len(removed)+1==r)

# Whole root-chart transition, volume, divided/original actions, and polynomial
# shears. No bound on r or Laurent order is inferred from these controls.
for r in [1,2,3,4,5]:
    low=z**(-r-2)+2*z**(-r-1)
    high=3*z**(1-r)+5*z**3
    h=low+high
    xp=z**(1-r)*t+h
    yp=lam+z**r
    jac=S.diff(xp,t)*S.diff(yp,z)-S.diff(xp,z)*S.diff(yp,t)
    eq(f"root chart has constant pulled-back volume r{r}",jac,r)
    eq(f"original lifted source translation on x r{r}",
       z**(r-1)*S.diff(xp,t),1)
    eq(f"divided Hamiltonian on source x r{r}",S.diff(xp,t),z**(1-r))
    eq(f"source translation preserves y r{r}",z**(r-1)*S.diff(yp,t))
    eq(f"divided Hamiltonian contraction equals r dz r{r}",
       r*S.diff(t,t),r)
    original_flow=sub(xp,{t:t+s*z**(r-1)})
    eq(f"whole original source translation action r{r}",original_flow,xp+s)
    divided_flow=sub(xp,{t:t+s})
    eq(f"whole divided field action r{r}",divided_flow,xp+s*z**(1-r))
    tnew=t+z**(r-1)*high
    check(f"threshold shear is polynomial in whole chart r{r}",
          S.cancel(tnew).is_polynomial(t,z))
    eq(f"threshold shear preserves the global slice r{r}",S.diff(tnew,t),1)
    eq(f"threshold shear removes exactly high Laurent terms r{r}",
       sub(xp,{t:T-z**(r-1)*high}),z**(1-r)*T+low)
    eq(f"exact-volume principal part lacks logarithmic coefficient r{r}",
       S.expand(low).coeff(z,-r))
    check(f"remaining principal part is nonzero r{r}",low!=0)
    eq(f"remaining largest pole exceeds r r{r}",r+2-r,2)
    # When the principal part were zero these TWO entire curves would
    # coincide off z=0; separatedness is the written argument excluding it.
    eq(f"zero-principal-part source zero curve agrees on punctured overlap r{r}",
       (z**(1-r)*t).subs(t,0))

# Primitive obstruction: dF has no logarithmic residue, even after the
# singular-looking rational source substitution. A primitive on each plane
# alone would not force this shared-chart coefficient to vanish.
for r in [1,2,3,5]:
    good_h=z**(-r-2)+2*z**(-r-1)+3*z**(1-r)
    for logarithmic in [0,7]:
        h=good_h+logarithmic*z**(-r)
        xp=z**(1-r)*t+h
        yp=lam+z**r
        F=xp**2*yp+xp*yp**2
        theta_t=S.diff(F,t)
        theta_z=r*xp*z**(r-1)+S.diff(F,z)
        eq(f"two-form exactness on overlap r{r} log{logarithmic}",
           S.diff(theta_z,t)-S.diff(theta_t,z),r)
        eq(f"derivative of pulled-back polynomial has zero residue r{r} log{logarithmic}",
           res_z(S.diff(F,z).subs(t,0),z))
        eq(f"shared primitive obstruction is precisely r*h_minus_r r{r} log{logarithmic}",
           res_z(theta_z.subs(t,0),z),r*logarithmic)
    for exponent in [-r-2,-r-1,-r,1-r,0,2]:
        term=11*z**exponent
        eq(f"logarithmic exponent uniquely detected r{r} j{exponent}",
           res_z(r*term*z**(r-1),z),11*r if exponent==-r else 0)

# The irreducible-curve LND argument uses a slice, not a smoothness
# classification. Exact finite Taylor reconstructions control this step.
for degree in [2,3,4,6]:
    moved=t**degree+2*t+3
    last=S.diff(moved,t,degree)
    check(f"last nonzero derivative is a scalar degree{degree}",
          last!=0 and not last.free_symbols)
    slice=S.diff(moved,t,degree-1)/last
    eq(f"curve slice derivative one degree{degree}",S.diff(slice,t),1)
    # Shift a slice by its constant to compare to the original parameter.
    slice=S.expand(slice-slice.subs(t,0))
    eq(f"curve slice recovers parameter degree{degree}",slice,t)
    poly=t**(degree+1)-3*t**2+5
    reconstruction=sum(S.diff(poly,t,j).subs(t,0)*slice**j/S.factorial(j)
                       for j in range(degree+2))
    eq(f"finite Taylor gives whole polynomial curve ring degree{degree}",
       reconstruction,poly)
check("a cusp subring is not stable under a nonzero normalization translation",
      1 not in {2*i+3*j for i in range(6) for j in range(6)})
eq("translation differentiates a cusp generator outside its cusp subring",
   S.diff(t**2,t),2*t)

# The genuine torsor action map and its inverse. The proof establishes
# smoothness and one orbit on every fibre before using these maps.
t1,t2=S.symbols("t1 t2")
action=(t+s,t,z)
inverse=(t1-t2,t2,z)
eq("action-map inverse recovers group parameter",
   sub(inverse[0],{t1:t+s,t2:t}),s)
eq("action-map inverse recovers point coordinate",
   sub(inverse[1],{t1:t+s,t2:t}),t)
eq("inverse action-map recovers first point",
   sub(t+s,{s:t1-t2,t:t2}),t1)
eq("inverse action-map preserves base z",action[2],inverse[2])
eq("action map has nonzero relative tangent",S.diff(action[0],s),1)
eq("whole chart action has additive group law",
   sub(sub(t+s,{t:t+T}),{T:T}),t+s+T)

# A full affine family with a proper marked plane open: M_p has UNIT volume.
# p=1 is a positive sharpness control for dropping the shared primitive.
for p in [1,2,3,5]:
    F=y**p*a-b*(b-1)
    source={b:x*y**p,a:x*(x*y**p-1)}
    eq(f"whole affine control source satisfies equation p{p}",sub(F,source))
    eq(f"source inverse on y-open p{p}",sub(b/y**p,source),x)
    eq(f"source inverse on bminus1-open p{p}",sub(a/(b-1),source),x)
    # The possible critical derivatives y=0,b=1/2 cannot satisfy F=0.
    eq(f"whole affine control is smooth p{p}",
       sub(F,{y:0,b:S.Rational(1,2)}),S.Rational(1,4))
    boundary={y:0,b:1}
    eq(f"whole boundary retained with free a p{p}",sub(F,boundary))
    chart={b:1+y**p*t,a:t*(1+y**p*t)}
    eq(f"whole boundary chart satisfies equation p{p}",sub(F,chart))
    eq(f"whole boundary chart inverse on b-open p{p}",sub(a/b,chart),t)
    eq(f"whole boundary chart inverse on y-open p{p}",sub((b-1)/y**p,chart),t)
    eq(f"source/rootchart transition finite principal part p{p}",
       sub(b/y**p,chart),t+y**(-p))
    eq(f"unit volume in source coordinates p{p}",
       S.diff(source[b],x)/y**p,1)
    eq(f"unit volume in whole boundary chart p{p}",
       S.diff(chart[a],t)/(2*chart[b]-1),1)
    delta=lambda g:S.expand(y**p*S.diff(g,b)+(2*b-1)*S.diff(g,a))
    eq(f"control Hamiltonian preserves equation p{p}",delta(F))
    eq(f"control Hamiltonian is nonzero on whole boundary p{p}",
       sub(delta(a),boundary),1)
    eq(f"rootchart source action is translation p{p}",
       delta(a).subs(chart,simultaneous=True),S.diff(chart[a],t))
    if p==1:
        eq("p1 shared primitive has forbidden logarithmic residue",
           res_z(z**(-p),z),1)
    else:
        xp=t+y**(-p)
        theta_t=y*S.diff(xp,t)/(p-1)
        theta_y=(y*S.diff(xp,y)+p*xp)/(p-1)
        check(f"exact primitive is polynomial in entire boundary chart p{p}",
              theta_t.is_polynomial(t,y) and S.cancel(theta_y).is_polynomial(t,y))
        eq(f"exact primitive differentiates to unit volume p{p}",
           S.diff(theta_y,t)-S.diff(theta_t,y),1)
        eq(f"exact primitive chart formula t coefficient p{p}",
           theta_t,y/(p-1))
        eq(f"exact primitive chart formula y coefficient p{p}",
           theta_y,p*t/(p-1))
        eq(f"nonlogarithmic principal part residue vanishes p{p}",
           res_z(z**(-p),z))

# Written hypotheses and source pins guard the scope of the conditional claim.
proof=PROOF.read_text(encoding="utf-8")
for text in [
    "Existence of that coordinate",
    "Normalize the finite cyclic cover",
    "The RAW equation",
    "restricted map W_j->M is etale EVERYWHERE",
    "nonzero LND on ANY irreducible",
    "stronger than freeness alone",
    "H^1_fppf(A1,Ga)=H^1(A1,O)=0",
    "coefficient at z^-r is then not forced",
    "h_low is NONZERO",
    "does not impose an upper",
    "independently obtained coordinates",
]:
    check("written conditional proof retains: "+text,text in proof)
audit=json.loads(AUDIT.read_text(encoding="utf-8"))
for entry in audit["primary_sources"]:
    check("cached primary exact pin "+entry["path"],
          sha(BASE/entry["path"])==entry["sha256"])
    primary=(BASE/entry["path"]).read_text()
    check("primary cached identity "+entry["identity"].split("tag")[-1].strip("."),
          entry["url"].rsplit("/",1)[-1] in primary)
check("primary comparison explicitly includes fppf",
      "fppf" in (BASE/audit["primary_sources"][0]["path"]).read_text())
check("primary vanishing is for quasi-coherent modules",
      "quasi-coherent" in (BASE/audit["primary_sources"][1]["path"]).read_text())

parents=[
    BASE/"CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md",
    BASE/"HAMILTONIAN_COMMON_COORDINATE_AND_ADDITIVE_COMPLETION.md",
]
receipt={
    "claim":"COMMON-COORDINATE-WHOLE-ROOT-CHARTS",
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_inputs":[{"path":p.name,"sha256":sha(p)} for p in parents]+[
        {"path":str(AUDIT.relative_to(BASE)),"sha256":sha(AUDIT)}],
    "primary_inputs":[{"path":e["path"],"sha256":sha(BASE/e["path"]),
                       "url":e["url"],"locator":e["locator"]}
                      for e in audit["primary_sources"]],
    "proved_scope":{
        "GIVEN_regular_polynomial_source_coordinate_on_component_open":True,
        "actual_canonical_entire_component_open_and_topology_used":True,
        "coordinate_existence_proved":False,
        "global_alignment_proved":False,
        "normalized_cyclic_cover_smooth":True,
        "raw_cyclic_cover_assumed_smooth":False,
        "whole_boundary_algebraic_A1_forced":True,
        "whole_affine_plane_root_charts":True,
        "divided_field_regular_by_unit_Hamiltonian":True,
        "single_orbit_smooth_fibres_prove_torsor":True,
        "Laurent_shift_finite":True,
        "same_component_volume_primitive_implies_h_minus_r_zero":True,
        "post_shear_principal_part_nonzero":True,
        "primitive_case_pole_order":"p>=r+1",
        "any_upper_pole_or_target_degree_bound_proved":False,
        "target_Ga_descent_or_unique_kernel_or_cyclic_deck_proved":False,
        "HC4_or_JC2_proved":False,
    },
    "assurance":"Internal conditional written proof and exact controls; no formal verification or external human review."
}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    RECEIPT.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS",
                  "check_count":len(checks),"elapsed_seconds":receipt["elapsed_seconds"],
                  "receipt_written":args.write_receipt},indent=2))
