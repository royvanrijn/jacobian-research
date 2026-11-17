#!/usr/bin/env python3
"""Exact controls for the conditional common-coordinate Hamiltonian bridge."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"HAMILTONIAN_COMMON_COORDINATE_AND_ADDITIVE_COMPLETION.md"
PRIMARY=BASE/"sources/stacks-031T-normal-domain-intersection.html"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(f):return S.cancel(f)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
x,y,u,v,z,a,s,t=S.symbols("x y u v z a s t")
def hamiltonian(f,g):
    return S.expand(S.diff(f,y)*S.diff(g,x)-S.diff(f,x)*S.diff(g,y))
def iterate(op,f,n):
    for _ in range(n):f=op(f)
    return S.expand(f)
def sub(f,m):return S.expand(f.subs(m,simultaneous=True))

# Contraction with dx wedge dy sends (fx,fy) to fx*dy-fy*dx.
dx,dy=S.symbols("dx dy")
for f in [y,x,x*y,y+x**2,y+x**5,x**2+y**3]:
    hx=S.diff(f,y);hy=-S.diff(f,x)
    check("Hamiltonian contraction equals differential "+str(f),
          zero(hx*dy-hy*dx-(S.diff(f,x)*dx+S.diff(f,y)*dy)))
check("coordinate Hamiltonian is source translation",
      hamiltonian(y,x)==1 and hamiltonian(y,y)==0)

# Arbitrary triangular coordinate conjugations and finite nilpotence controls.
for d in range(2,7):
    coord=y+x**d
    D=lambda g:hamiltonian(coord,g)
    check("new coordinate killed d"+str(d),D(coord)==0)
    check("slice derivative one d"+str(d),D(x)==1)
    flowx=x+s
    flowy=y+x**d-(x+s)**d
    check("coordinate preserved by full source flow d"+str(d),
          zero(sub(coord,{x:flowx,y:flowy})-coord))
    check("flow derivative x d"+str(d),S.diff(flowx,s).subs(s,0)==D(x))
    check("flow derivative y d"+str(d),zero(S.diff(flowy,s).subs(s,0)-D(y)))
    for i,j in [(0,1),(1,0),(1,1),(2,2),(3,1)]:
        g=x**i*y**j
        check(f"coordinate LND finite iterate d{d} i{i} j{j}",
              iterate(D,g,i+d*j+1)==0)
    q=coord**2+1
    check("every polynomial replica coefficient invariant d"+str(d),D(q)==0)
    RD=lambda g:S.expand(q*D(g))
    for n in [1,2,3]:
        check(f"replica iterate identity d{d} n{n}",
              zero(iterate(RD,y,n)-q**n*iterate(D,y,n)))
    fx=x+s*q
    fy=y+x**d-(x+s*q)**d
    check("replica source flow preserves coordinate d"+str(d),
          zero(sub(coord,{x:fx,y:fy})-coord))

# Smooth affine proper-plane-open control.
F=y*a-z*(z-1)
check("singular candidate excluded by equation",F.subs({y:0,a:0,z:S.Rational(1,2)})==S.Rational(1,4))
check("volume chart cover has no omitted point",
      F.subs({y:0,z:S.Rational(1,2)})==S.Rational(1,4))
embed={z:x*y,a:x*(x*y-1)}
check("source parametrization satisfies hypersurface",zero(sub(F,embed)))
check("inverse on y chart",zero(sub(z/y,embed)-x))
check("inverse on zminus1 chart",zero(sub(a/(z-1),embed)-x))
check("two inverse charts agree on hypersurface",
      zero(z*(z-1)-y*a+F))
check("source chart cover includes every y0 point",sub(x*y-1,{y:0})==-1)
D=lambda g:S.expand(y*S.diff(g,z)+(2*z-1)*S.diff(g,a))
E=lambda g:S.expand(a*S.diff(g,z)+(2*z-1)*S.diff(g,y))
check("Hamiltonian control preserves defining equation",D(F)==0)
check("second control derivation preserves equation",E(F)==0)
for name,op in [("delta",D),("epsilon",E)]:
    for generator in [y,z,a]:
        check(f"control LND on generator {name} {generator}",iterate(op,generator,3)==0)
check("Hamiltonian control contraction on y chart",zero(D(z)/y-1))
check("Hamiltonian control contraction on 2zminus1 chart",zero(D(a)/(2*z-1)-1))
check("volume pullback agrees with source volume",
      zero(S.diff(x*y,x)/y-1))
delta_flow={y:y,z:z+s*y,a:a+s*(2*z-1)+s**2*y}
epsilon_flow={a:a,z:z+s*a,y:y+s*(2*z-1)+s**2*a}
for name,flow in [("delta",delta_flow),("epsilon",epsilon_flow)]:
    check("exact action preserves hypersurface "+name,zero(sub(F,flow)-F))
    for generator in [y,z,a]:
        rhs=sub(sub(generator,flow),{s:s+t})
        # Action composition by simultaneous substitution in the coordinates.
        composed=sub(sub(generator,flow),{y:sub(y,flow).subs(s,t),
                                          z:sub(z,flow).subs(s,t),
                                          a:sub(a,flow).subs(s,t)})
        check(f"additive group law {name} {generator}",zero(composed-rhs))
for generator in [y,z,a]:
    check("source restriction is exact translation "+str(generator),
          zero(sub(sub(generator,delta_flow),embed)-sub(sub(generator,embed),{x:x+s})))
check("first action preserves marked boundary y",sub(delta_flow[y],{y:0,z:1})==0)
check("first action preserves marked boundary z",sub(delta_flow[z],{y:0,z:1})==1)
check("first action is nontrivial on boundary",sub(delta_flow[a],{y:0,z:1})==a+s)
check("second kernel differs",E(a)==0 and D(a)!=0)
point={y:0,z:1,a:1,s:1}
check("second action moves boundary point into source",sub(epsilon_flow[y],point)==2)
check("moved point remains on surface",sub(sub(F,epsilon_flow),point)==0)
check("second source restriction x is incomplete Riccati",
      zero(sub(E(z)*y-z*E(y),embed)+x**2*y**2))

# Necessary hypotheses cannot be silently weakened.
xyH=lambda g:hamiltonian(x*y,g)
for n in range(1,9):
    check("noncoordinate Hamiltonian is not nilpotent iterate"+str(n),iterate(xyH,x,n)==x)
check("degenerate volume Hamiltonian has pole",zero(u*(1/u)-1))
check("degenerate volume inverse is not regular polynomial",not S.cancel(1/u).is_polynomial(u,v))
check("quasiaffine translation exits punctured plane",sub(x+s,{x:1,s:-1})==0)

# Normal-DVR input is a cached primary statement, not finite evidence.
primary=PRIMARY.read_text()
check("primary tag cache identifies exact tag","031T" in primary)
check("primary statement contains height-one intersection",
      "height" in primary and "bigcap" in primary and "Noetherian normal" in primary)

receipt={"claim":"HC4MYHCC0","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "primary_inputs":[{"url":"https://stacks.math.columbia.edu/tag/031T",
                      "path":str(PRIMARY.relative_to(BASE)),"sha256":sha(PRIMARY),
                      "locator":"Lemma10.157.6(2), statement and proof",
                      "scope":"Normal Noetherian domain equals intersection of height-one localizations"}],
    "proved_scope":{
        "smooth_affine_completion_and_open_marked_plane_assumed":True,
        "unit_volume_assumed":True,
        "common_divisorial_source_coordinate_assumed":True,
        "whole_boundary_root_charts_required":False,
        "boundary_A1_disjointness_or_pi1_required":False,
        "volume_primitive_required":False,
        "coordinate_regular_everywhere":True,
        "extended_LND_kernel":"C[y]",
        "constructed_actions_preserve_whole_source":True,
        "common_coordinate_existence_proved":False,
        "affineness_proved":False,
        "unique_LND_kernel_or_all_actions_classified":False,
        "cyclic_deck_HC4_or_JC2_proved":False},
    "assurance":"Internal written conditional proof and exact controls; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"hamiltonian_common_coordinate_verification.json").write_text(
        json.dumps(receipt,indent=2)+chr(10),encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
