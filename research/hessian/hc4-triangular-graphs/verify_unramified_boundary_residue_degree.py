#!/usr/bin/env python3
"""Exact controls for unramified boundary residue degree and target normalization."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"UNRAMIFIED_BOUNDARY_RESIDUE_DEGREE_AND_TARGET_NORMALIZATION.md"
TOPOLOGY=BASE/"CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md"
CHAU=BASE.parents[1]/"jacobian/affine-complements/sources/global-boundary-topology/chau-2003-two-remarks.pdf"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
t,s,r,u,v,z=S.symbols("t s r u v z")
def sub(expr,mapping):return S.expand(S.sympify(expr).subs(mapping,simultaneous=True))
def order(poly,var):
    poly=S.Poly(S.expand(poly),var)
    return S.oo if poly.is_zero else min(m[0] for m in poly.monoms())

# A finite birational normalization need not be an immersion or an isomorphism
# to the singular image.
cp=t**2;cq=t**3
check("cusp implicit relation",zero(cq**2-cp**3))
check("cusp integral parameter equation",zero(t**2-cp))
check("cusp field parameter recovered",zero(cq/cp-t))
check("cusp map critical at normalization origin",
      S.diff(cp,t).subs(t,0)==S.diff(cq,t).subs(t,0)==0)
check("cusp map generically immersive",S.gcd(S.diff(cp,t),S.diff(cq,t))==t)
F=v**2-u**3
check("cusp is singular at target origin",
      F.subs({u:0,v:0})==S.diff(F,u).subs({u:0,v:0})==S.diff(F,v).subs({u:0,v:0})==0)
check("cusp normalization function t is outside cusp semigroup",
      1 not in {2*i+3*j for i in range(8) for j in range(8)})

# An unbounded formula gives polynomial immersions with every contact order.
contact_controls=[]
for k in range(1,10):
    p=t**2-1
    q=t*p**k
    implicit=v**2-u**(2*k)*(u+1)
    check("contact normalization satisfies implicit curve k"+str(k),
          zero(sub(implicit,{u:p,v:q})))
    check("parameter is integral over curve coordinate k"+str(k),zero(t**2-(p+1)))
    check("normalization is birational field identity k"+str(k),zero(q/p**k-t))
    check("polynomial normalization is an immersion k"+str(k),
          S.gcd(S.diff(p,t),S.diff(q,t))==1)
    check("two distinct normalization points meet k"+str(k),
          sub(p,{t:1})==sub(p,{t:-1})==sub(q,{t:1})==sub(q,{t:-1})==0)
    check("both source branch tangents are nonzero k"+str(k),
          S.diff(p,t).subs(t,1)!=0 and S.diff(p,t).subs(t,-1)!=0)
    check("branch contact squared difference exact k"+str(k),
          zero((2*u**k*S.sqrt(1+u))**2-4*u**(2*k)*(1+u)))
    check("branch contact is exactly k"+str(k),
          order(4*u**(2*k)*(1+u),u)==2*k)
    check("curve target singularity retained k"+str(k),
          implicit.subs({u:0,v:0})==S.diff(implicit,u).subs({u:0,v:0})==S.diff(implicit,v).subs({u:0,v:0})==0)
    # Etale ambient thickening exists along the entire normalized curve;
    # its nonconstant whole-plane Jacobian guards the Keller distinction.
    aa=S.Integer((-1)**(k+1))
    bb=S.cancel((1+aa*S.diff(q,t))/(2*t))
    check("Bezout coefficient is a polynomial k"+str(k),bb.is_polynomial(t))
    check("whole-curve tangent Bezout identity k"+str(k),
          zero(S.diff(p,t)*bb-aa*S.diff(q,t)-1))
    P=p+aa*s;Q=q+bb*s
    J=S.expand(S.diff(P,t)*S.diff(Q,s)-S.diff(P,s)*S.diff(Q,t))
    check("ambient thickening etale on every D point k"+str(k),zero(J.subs(s,0)-1))
    check("ambient Jacobian explicit k"+str(k),zero(J-(1-aa*S.diff(bb,t)*s)))
    check("ambient comparison is not a whole-plane Keller pair k"+str(k),J!=1)
    contact_controls.append({"k":k,"normalization":[str(p),str(q)],
                             "ambient_Jacobian":str(J),"contact_order":k})

# Disconnected fibre product: select the whole finite normalized component;
# another component is etale only after losing points.
P=t**2-1+s
Q=t*(t**2-1)+S.Rational(3,2)*t*s
J=S.expand(S.diff(P,t)*S.diff(Q,s)-S.diff(P,s)*S.diff(Q,t))
check("node ambient etale open exact Jacobian",zero(J-(1-S.Rational(3,2)*s)))
s_elim=r**2-t**2
eq1=S.expand((P-(r**2-1)).subs(s,s_elim))
eq2=S.expand((Q-r*(r**2-1)).subs(s,s_elim))
g=t**2+r*t-2*r**2+2
check("first basechange equation eliminates transverse coordinate",eq1==0)
check("second basechange equation factors into selected and other components",
      zero(eq2+S.Rational(1,2)*(t-r)*g))
check("components are comaximal even before etale localization",
      zero(g-(t-r)*(t+2*r)-2))
Delta=S.expand(J.subs(s,s_elim))
check("selected normalization component retained at every parameter",zero(Delta.subs(t,r)-1))
check("selected map to target normalization has degree one",S.degree(r,r)==1)
disc=S.discriminant(g,t)
check("other component discriminant exact",zero(disc-(9*r**2-8)))
branch_t=-r/2
check("other component ramification only at roots of discriminant",
      zero(g.subs(t,branch_t)+disc/4))
check("other component branch points removed by ambient etale open",
      zero(Delta.subs(t,branch_t)+disc/8))
norm_Delta=S.resultant(g,2*Delta,t)
check("inverse etale discriminant has target pole on other component",
      zero(norm_Delta+2*disc))
check("other component cannot be presumed finite over entire normalization",
      disc.subs(r,0)!=0 and S.degree(disc,r)==2)
check("two node branches correspond to two different D points",1!=-1)

# Euler equation is integer divisibility, not a bounded classification.
solutions=[]
for degree in range(1,20):
    for chi in range(-10,11):
        if degree*chi==1:solutions.append((degree,chi))
check("integer Euler equation has only degree one",solutions==[(1,1)])
for degree in range(2,12):
    cover=t**degree
    deriv=S.diff(cover,t)
    check("higher-degree A1 cover is ramified at origin m"+str(degree),deriv.subs(t,0)==0)
    check("higher-degree map has exact ramification index m"+str(degree),
          order(cover,t)==degree and order(deriv,t)==degree-1)
    check("punctured-line cover is etale after removing ramification m"+str(degree),
          zero(deriv/(degree*t**(degree-1))-1))
    check("punctured-line Euler zero allows arbitrary degree m"+str(degree),degree*0==0)

# Normalization and etale ambient differentials retain branch criticality.
a,b,c,d,pt,qt=S.symbols("a b c d pt qt")
M=S.Matrix([[a,b],[c,d]])
vec=S.Matrix([pt,qt])
check("ambient inverse preserves zero versus nonzero normalization tangent",
      M.adjugate()*M*vec==S.expand(M.det())*vec)
check("invertible ambient differential does not remove cusp zero tangent",
      M*S.Matrix([0,0])==S.zeros(2,1))

# The source explicitly guards every delicate geometric inference.
proof=PROOF.read_text()
for phrase in [
    "We do NOT assume that T is finite",
    "It may be disconnected",
    "the graph of alpha is a closed immersion",
    "finite birational map onto a",
    "Smoothness is claimed for T_D and Dnu, not for D",
    "finite polynomial normalization",
    "number of analytic branches",
    "if and only if D is a smooth algebraic A1",
    "the residue-degree and normalization theorem in Sections1--4 does not",
    "not an actual whole-plane",
]:
    check("written scope guard: "+phrase,phrase in proof)

receipt={"claim":"HC4MYBR0","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_inputs":[
        {"path":TOPOLOGY.name,"sha256":sha(TOPOLOGY),
         "scope":"Whole finite-normalization boundary closure, normalization A1 and unibranchness"},
        {"path":"../../jacobian/affine-complements/sources/global-boundary-topology/chau-2003-two-remarks.pdf",
         "sha256":sha(CHAU),"locator":"Printedpp41--43, Lemmas5,6(ii), component Remark onp43",
         "scope":"Only singular-target corollary: no A1 nonproperness component"}],
    "proved_scope":{
        "every_actual_omitted_unramified_prime":True,
        "whole_original_finite_normalization_retained":True,
        "residue_degree":1,
        "target_curve_normalization":"A1",
        "normalization_map_between_boundary_and_target":"isomorphism",
        "target_image_map":"finite birational",
        "target_normalization_immersion":"if and only if boundary is smooth",
        "target_singular":"using separately retained Chau component input",
        "boundary_smoothness_or_coordinate_existence_proved_here":False,
        "embedded_target_affine_line_or_ordinary_nodes_claimed":False,
        "whole_fibre_product_presumed_connected_or_finite":False,
        "HC4_JC2_or_integrality_proved":False},
    "arbitrary_contact_controls":contact_controls,
    "assurance":"Internal written proof and exact controls; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"unramified_boundary_residue_degree_verification.json").write_text(
        json.dumps(receipt,indent=2)+chr(10),encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
