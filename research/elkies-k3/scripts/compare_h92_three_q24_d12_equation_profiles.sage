#!/usr/bin/env sage -python
"""
Compare the THREE direct H3 D13 q24 -> D12/MW5 orbits in the ACTUAL current
canonical D13 equation frame.

This is the equation-complexity tie-breaker missing from the original lattice
search.  All three candidates are transported by the exact same isometries:

  pinned dominant D13
    -> inverse nef-to-dominant q6-simple reflections
    -> source H3
    -> 102 actual E6/E8 component reflections
    -> raw component-nef q8 presentation
    -> 22 q6 Weyl reflections
    -> q6 Eichler translation tau
    -> current equation D13.

Then, using exactly the same equation zero / deterministic D13 adaptation as
close_h92_q8_q24_by_q6_translation.sage, each divisor is decomposed as

    D = O8 + P + V + k F8.

Report:
  * equation D13 MW coordinates of P;
  * height, D13 correction and P.O (collision degree);
  * degree of P on the preceding q6 fibration;
  * vertical fibre coefficient, root L1 and root support;
  * child root data (must remain D12);
  * simple cost ordering.

Orbit 85 is a mandatory regression:
    MW=(-2,1,-1,1), height=52, corr=0, P.O=24,
    vertical_F=-7, root_L1=69, root_support=13.
"""

import argparse
import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, gcd, matrix, vector


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    h=Path.home()
    candidates += [
        h/"Documents"/"jacobian-research",
        h/"jacobian-research",
        h/"src"/"jacobian-research",
        h/"git"/"jacobian-research",
        h/"projects"/"jacobian-research",
    ]
    seen=set()
    for c in candidates:
        try:
            c=c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate jacobian-research")


def run_scope(path):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    capture=io.StringIO()
    try:
        sys.argv=[str(path)]
        with contextlib.redirect_stdout(capture):
            exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope, capture.getvalue()


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
GEN=ROOT/"artifacts/generated-results"
LOCAL=ROOT/"artifacts/local/elkies-k3"

Q24SEARCH=GEN/"elkies-k3-h3-q6-q8-d13-q24-degree2.json"
PULLBACK=ROOT/"elkies-k3/scripts/audit_h3_d13_q24_section_pullback_v4.sage"
CLOSE=ROOT/"elkies-k3/scripts/close_h92_q8_q24_by_q6_translation.sage"

for path in (Q24SEARCH,PULLBACK,CLOSE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

search=json.loads(Q24SEARCH.read_text())
assert search["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
records=[
    r for r in search["neighbors"]
    if tuple(r["child_root_data"])==(12,264,4)
    and int(r["child_mw_rank"])==5
]
assert sorted(int(r["orbit_index"]) for r in records)==[73,85,184]

print("Q24ALL3EQ|stage=load_transport_scopes|status=START",flush=True)
pb,pb_log=run_scope(PULLBACK)
cl,cl_log=run_scope(CLOSE)
print("Q24ALL3EQ|stage=load_transport_scopes|status=PASS",flush=True)

# Required exact transport infrastructure from the passing scripts.
pb_need=(
    "pinned_d13","Bdom_simple","q6_simple_ns","simple_ns_roots",
    "nef_to_dominant","Bsimple","source_ns","actual_roots",
    "physical_reflections","reflect_row","F_component",
)
cl_need=(
    "ns","F6","F8eq","O8","Badapt","adapted","coords",
    "minimal_section_for_mw","vertical_basis","weyl_transport","tau",
    "isotropic_mate","roots_and_data",
)
missing_pb=[k for k in pb_need if k not in pb]
missing_cl=[k for k in cl_need if k not in cl]
if missing_pb or missing_cl:
    raise SystemExit(
        "embedded scripts missing variables: pullback={} close={}".format(
            missing_pb,missing_cl
        )
    )

pinned_d13=pb["pinned_d13"]
Gpin=matrix(ZZ,[[0,1]+[0]*17,[1,0]+[0]*17]+[
    [0,0]+[-pinned_d13[i,j] for j in range(17)] for i in range(17)
])
# Simpler exact pinned NS from the pullback scope when available.
if "G_d13" in pb:
    Gpin=pb["G_d13"]

Bdom_simple=pb["Bdom_simple"]
q6_simple_ns=pb["q6_simple_ns"]
simple_ns_roots=pb["simple_ns_roots"]
nef_to_dominant=pb["nef_to_dominant"]
Bsimple=pb["Bsimple"]
source_ns=pb["source_ns"]
actual_roots=pb["actual_roots"]
physical_reflections=pb["physical_reflections"]
reflect_row=pb["reflect_row"]
Fraw=vector(ZZ,pb["F_component"])

ns=cl["ns"]
assert ns==source_ns
F6=vector(ZZ,cl["F6"])
F8eq=vector(ZZ,cl["F8eq"])
O8=vector(ZZ,cl["O8"])
Badapt=cl["Badapt"]
adapted=cl["adapted"]
coords=cl["coords"]
minimal_section_for_mw=cl["minimal_section_for_mw"]
vertical_basis=cl["vertical_basis"]
weyl_transport=cl["weyl_transport"]
tau=cl["tau"]
isotropic_mate=cl["isotropic_mate"]
roots_and_data=cl["roots_and_data"]

def pinned_to_raw_component(Dpin):
    D=vector(ZZ,Dpin)*Bdom_simple

    # Dominant -> nef: inverse of stored nef->dominant reflections.
    for index1,unused_label in reversed(nef_to_dominant):
        root=simple_ns_roots[index1-1]
        D=reflect_row(D,q6_simple_ns,root)

    # q6 simple -> source H3.
    D=vector(ZZ,D*Bsimple)

    # Same 102 actual physical component reflections used on the q8 fibre.
    for i,unused_pairing in physical_reflections:
        D=reflect_row(D,source_ns,actual_roots[i])

    assert D*source_ns*D==0
    assert D*source_ns*Fraw==2
    return D

def raw_to_equation(Draw):
    Dphysical=vector(ZZ,weyl_transport(vector(ZZ,Draw)))
    Deq=vector(ZZ,tau(Dphysical))
    assert Deq*ns*Deq==0
    assert Deq*ns*F8eq==2
    assert gcd(tuple(Deq))==1
    return Deq

profiles=[]
for rec in sorted(records,key=lambda r:int(r["orbit_index"])):
    orbit=int(rec["orbit_index"])
    witness=vector(ZZ,rec["witness"])
    Dpin=vector(ZZ,rec["fiber"])

    # Original pinned diagnostics.
    Ppin=vector(ZZ,[23,1]+list(witness))
    assert Dpin==vector(ZZ,[-1,1]+[0]*17)+Ppin-10*vector(ZZ,[1,0]+[0]*17)
    assert Ppin*Gpin*Ppin==-2
    assert Ppin*Gpin*vector(ZZ,[1,0]+[0]*17)==1
    assert Ppin*Gpin*vector(ZZ,[-1,1]+[0]*17)==22

    Draw=pinned_to_raw_component(Dpin)
    Deq=raw_to_equation(Draw)

    # Child regression.
    mate=isotropic_mate(ns,Deq)
    orth=matrix(ZZ,[list(Deq*ns),list(mate*ns)]).right_kernel_matrix()
    child=-(orth*ns*orth.transpose())
    root_data=tuple(map(int,roots_and_data(child)[2]))
    assert root_data==(12,264,4)

    cD=coords(Deq)
    assert cD[1]==2
    z=vector(ZZ,cD[-4:])
    prof=minimal_section_for_mw(adapted,z)
    P=vector(ZZ,prof["section"]*Badapt)

    assert P*ns*P==-2
    assert P*ns*F8eq==1

    vertical=Deq-O8-P
    assert vertical*ns*F8eq==0
    coeff=vertical_basis.solve_left(vector(QQ,vertical))
    assert all(v in ZZ for v in coeff)
    coeff=vector(ZZ,coeff)
    vf=ZZ(coeff[0])
    vr=vector(ZZ,coeff[1:])

    q6_degree=ZZ(P*ns*F6)
    q6_zero=ZZ(P*ns*cl["Ostd"]) if "Ostd" in cl else None

    item={
        "orbit":orbit,
        "original_mw_projection":list(map(int,rec.get("mw_projection",[]))),
        "original_witness":list(map(int,witness)),
        "original_witness_l1":sum(abs(int(v)) for v in witness),
        "original_witness_linf":max(abs(int(v)) for v in witness),
        "equation_mw":list(map(int,z)),
        "height":str(prof["height"]),
        "D13_correction":str(prof["correction"]),
        "P_dot_O":int(prof["pole"]),
        "P_q6_degree":int(q6_degree),
        "P_q6_zero_intersection":None if q6_zero is None else int(q6_zero),
        "vertical_fibre_coefficient":int(vf),
        "vertical_root_coefficients":list(map(int,vr)),
        "vertical_root_L1":sum(abs(int(v)) for v in vr),
        "vertical_root_support":sum(bool(v) for v in vr),
        "equation_divisor":list(map(int,Deq)),
        "equation_horizontal_section":list(map(int,P)),
        "child_root_data":list(root_data),
    }
    profiles.append(item)

    print(
        "Q24ALL3EQ_PROFILE|"
        f"orbit={orbit}|"
        f"eq_mw={','.join(map(str,z))}|"
        f"height={prof['height']}|corr={prof['correction']}|PdotO={prof['pole']}|"
        f"P_q6_degree={q6_degree}|"
        f"vertical_F={vf}|root_L1={item['vertical_root_L1']}|"
        f"root_support={item['vertical_root_support']}|"
        f"wL1={item['original_witness_l1']}|wLinf={item['original_witness_linf']}|"
        "child=D12|status=PASS",
        flush=True,
    )

# Exact orbit-85 regression against the already passing transport.
p85=next(x for x in profiles if x["orbit"]==85)
assert p85["equation_mw"]==[-2,1,-1,1],p85
assert p85["height"]=="52"
assert p85["D13_correction"]=="0"
assert p85["P_dot_O"]==24
assert p85["vertical_fibre_coefficient"]==-7
assert p85["vertical_root_L1"]==69
assert p85["vertical_root_support"]==13
print(
    "Q24ALL3EQ_REGRESSION|orbit=85|mw=-2,1,-1,1|height=52|PdotO=24|"
    "vertical_F=-7|root_L1=69|support=13|status=PASS",
    flush=True,
)

# Rank by the quantities that actually inflate equation-level work:
# collision degree first, then preceding-fibration degree, then vertical support/L1.
ranking=sorted(
    profiles,
    key=lambda x:(
        x["P_dot_O"],
        x["P_q6_degree"],
        x["vertical_root_support"],
        x["vertical_root_L1"],
        abs(x["vertical_fibre_coefficient"]),
        x["original_witness_l1"],
    )
)
print(
    "Q24ALL3EQ_RANKING|"
    + "|".join(
        "orbit{}:PO={},q6deg={},VF={},VL1={},Vsupp={}".format(
            x["orbit"],x["P_dot_O"],x["P_q6_degree"],
            x["vertical_fibre_coefficient"],x["vertical_root_L1"],
            x["vertical_root_support"]
        )
        for x in ranking
    )
    + "|status=PASS",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-three-d12-equation-profiles.v1",
    "status":"PASS_EXACT_THREE_Q24_D12_EQUATION_PROFILES",
    "profiles":profiles,
    "ranking":[x["orbit"] for x in ranking],
    "ranking_rule":[
        "P_dot_O / smooth-collision degree",
        "horizontal P degree on preceding q6 fibration",
        "vertical root support",
        "vertical root L1",
        "absolute vertical fibre coefficient",
        "original pinned witness L1",
    ],
    "boundary":(
        "This is an exact NS/equation-frame complexity comparison. It does not "
        "construct rational coordinates for the orbit-73 or orbit-184 horizontal "
        "sections or their D12 pencils."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/"q24-three-d12-equation-profiles.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24ALL3EQ_RESULT|best_orbit={}|ranking={}|status={}".format(
        ranking[0]["orbit"],
        ",".join(str(x["orbit"]) for x in ranking),
        payload["status"],
    ),
    flush=True,
)
