#!/usr/bin/env sage -python
"""
Compare all three q24 -> D12/MW5 candidates in the CURRENT canonical D13
equation frame, using only the current passing transport layers.

Authoritative transport:
  pinned dominant D13
    -> current audit_h92_q8_q24_effective_zero_choices.sage
       (dominant -> nef -> source H3 -> 102 physical component reflections)
    -> current q6 Weyl transport
    -> current q6 Eichler translation
    -> canonical equation D13.

Orbit 85 must reproduce the already-certified equation profile:
  MW=(-2,1,-1,1), height=52, corr=0, P.O=24,
  vertical_F=-7, root_L1=69, support=13.
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


def run_scope(path, argv=()):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    capture=io.StringIO()
    try:
        sys.argv=[str(path)]+list(argv)
        with contextlib.redirect_stdout(capture):
            exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope,capture.getvalue()


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
GEN=ROOT/"artifacts/generated-results"
LOCAL=ROOT/"artifacts/local/elkies-k3"

Q24SEARCH=GEN/"elkies-k3-h3-q6-q8-d13-q24-degree2.json"
CURRENT_AUDIT=ROOT/"elkies-k3/scripts/audit_h92_q8_q24_effective_zero_choices.sage"
CURRENT_CLOSE=ROOT/"elkies-k3/scripts/close_h92_q8_q24_by_q6_translation.sage"
TMP=LOCAL/"q24-all3-effective-zero-temp.json"

for path in (Q24SEARCH,CURRENT_AUDIT,CURRENT_CLOSE):
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

print("Q24ALL3CURRENT|stage=current_effective_zero_audit|status=START",flush=True)
qa,qa_log=run_scope(CURRENT_AUDIT,("--output",str(TMP)))
print("Q24ALL3CURRENT|stage=current_effective_zero_audit|status=PASS",flush=True)

print("Q24ALL3CURRENT|stage=current_equation_close|status=START",flush=True)
cl,cl_log=run_scope(CURRENT_CLOSE)
print("Q24ALL3CURRENT|stage=current_equation_close|status=PASS",flush=True)

need_qa=(
    "Gpinned","Bpinned_to_simple","q6_simple_ns","simple_root_classes",
    "bridge","Bsimple","source_ns","actual_roots","physical_reflections",
    "reflect","F_actual",
)
need_cl=(
    "ns","F6","Ostd","F8eq","O8","Badapt","adapted","coords",
    "minimal_section_for_mw","vertical_basis","weyl_transport","tau",
    "isotropic_mate","roots_and_data",
)
mqa=[x for x in need_qa if x not in qa]
mcl=[x for x in need_cl if x not in cl]
if mqa or mcl:
    raise SystemExit(f"current scopes missing variables: audit={mqa} close={mcl}")

Gpinned=qa["Gpinned"]
Bpinned_to_simple=qa["Bpinned_to_simple"]
q6_simple_ns=qa["q6_simple_ns"]
simple_root_classes=qa["simple_root_classes"]
bridge=qa["bridge"]
Bsimple=qa["Bsimple"]
source_ns=qa["source_ns"]
actual_roots=qa["actual_roots"]
physical_reflections=qa["physical_reflections"]
reflect=qa["reflect"]
Fraw=vector(ZZ,qa["F_actual"])

ns=cl["ns"]
assert ns==source_ns
F6=vector(ZZ,cl["F6"])
Ostd=vector(ZZ,cl["Ostd"])
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

Fpin=vector(ZZ,[1,0]+[0]*17)
Opin=vector(ZZ,[-1,1]+[0]*17)


def pinned_to_raw(Dpin):
    # pinned dominant D13 -> q6 simple ambient
    D=vector(ZZ,Dpin)*Bpinned_to_simple

    # Invert current stored nef -> dominant reflection bridge.
    for i,unused_pairing in reversed(bridge):
        D=reflect(D,q6_simple_ns,simple_root_classes[i])

    # q6 simple -> source H3
    D=vector(ZZ,D*Bsimple)

    # Replay the exact 102 component reflections used by the current audit.
    for i,unused_pairing in physical_reflections:
        D=reflect(D,source_ns,actual_roots[i])

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

    assert Dpin==vector(ZZ,[12,2]+list(witness))
    assert Dpin*Gpinned*Dpin==0
    assert Dpin*Gpinned*Fpin==2
    assert Dpin*Gpinned*Opin==10

    # Pinned horizontal section invariants.
    Ppin=Opin+Dpin+10*Fpin-Opin  # simplify next line explicitly
    Ppin=vector(ZZ,[23,1]+list(witness))
    assert Dpin==Opin+Ppin-10*Fpin
    assert Ppin*Gpinned*Ppin==-2
    assert Ppin*Gpinned*Opin==22

    Draw=pinned_to_raw(Dpin)
    Deq=raw_to_equation(Draw)

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
    assert P*ns*P==-2 and P*ns*F8eq==1

    vertical=Deq-O8-P
    assert vertical*ns*F8eq==0
    coeff=vertical_basis.solve_left(vector(QQ,vertical))
    assert all(v in ZZ for v in coeff)
    coeff=vector(ZZ,coeff)
    vf=ZZ(coeff[0])
    vr=vector(ZZ,coeff[1:])

    item={
        "orbit":orbit,
        "pinned_mw_projection":list(map(int,rec.get("mw_projection",[]))),
        "pinned_witness":list(map(int,witness)),
        "pinned_witness_l1":sum(abs(int(v)) for v in witness),
        "pinned_witness_linf":max(abs(int(v)) for v in witness),
        "equation_mw":list(map(int,z)),
        "height":str(prof["height"]),
        "D13_correction":str(prof["correction"]),
        "P_dot_O":int(prof["pole"]),
        "P_q6_degree":int(P*ns*F6),
        "P_q6_standard_zero_intersection":int(P*ns*Ostd),
        "vertical_fibre_coefficient":int(vf),
        "vertical_root_coefficients":list(map(int,vr)),
        "vertical_root_L1":sum(abs(int(v)) for v in vr),
        "vertical_root_support":sum(bool(v) for v in vr),
        "child_root_data":list(root_data),
        "equation_divisor":list(map(int,Deq)),
        "equation_horizontal_section":list(map(int,P)),
    }
    profiles.append(item)

    print(
        "Q24ALL3CURRENT_PROFILE|"
        f"orbit={orbit}|eq_mw={','.join(map(str,z))}|"
        f"height={prof['height']}|corr={prof['correction']}|PdotO={prof['pole']}|"
        f"P_q6_degree={item['P_q6_degree']}|P_q6_O={item['P_q6_standard_zero_intersection']}|"
        f"vertical_F={vf}|root_L1={item['vertical_root_L1']}|"
        f"root_support={item['vertical_root_support']}|"
        f"wL1={item['pinned_witness_l1']}|wLinf={item['pinned_witness_linf']}|"
        "child=D12|status=PASS",
        flush=True,
    )

# Hard current orbit-85 regression.
p85=next(x for x in profiles if x["orbit"]==85)
assert p85["equation_mw"]==[-2,1,-1,1],p85
assert p85["height"]=="52"
assert p85["D13_correction"]=="0"
assert p85["P_dot_O"]==24
assert p85["vertical_fibre_coefficient"]==-7
assert p85["vertical_root_L1"]==69
assert p85["vertical_root_support"]==13

print(
    "Q24ALL3CURRENT_REGRESSION|orbit=85|mw=-2,1,-1,1|height=52|corr=0|"
    "PdotO=24|vertical_F=-7|root_L1=69|support=13|status=PASS",
    flush=True,
)

ranking=sorted(
    profiles,
    key=lambda x:(
        x["P_dot_O"],
        x["P_q6_degree"],
        x["vertical_root_support"],
        x["vertical_root_L1"],
        abs(x["vertical_fibre_coefficient"]),
        x["pinned_witness_l1"],
    )
)

print(
    "Q24ALL3CURRENT_RANKING|"
    +"|".join(
        "orbit{}:PO={},q6deg={},q6O={},VF={},VL1={},Vsupp={}".format(
            x["orbit"],x["P_dot_O"],x["P_q6_degree"],
            x["P_q6_standard_zero_intersection"],
            x["vertical_fibre_coefficient"],x["vertical_root_L1"],
            x["vertical_root_support"],
        )
        for x in ranking
    )
    +"|status=PASS",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-three-d12-current-equation-profiles.v1",
    "status":"PASS_EXACT_THREE_Q24_D12_CURRENT_EQUATION_PROFILES",
    "profiles":profiles,
    "ranking":[x["orbit"] for x in ranking],
    "ranking_rule":[
        "equation P.O / collision degree",
        "horizontal section degree on preceding q6 fibration",
        "vertical root support",
        "vertical root L1",
        "absolute vertical fibre coefficient",
        "pinned witness L1",
    ],
}
OUT=args.output.resolve() if args.output else LOCAL/"q24-three-d12-current-equation-profiles.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24ALL3CURRENT_RESULT|best_orbit={}|ranking={}|status={}".format(
        ranking[0]["orbit"],
        ",".join(str(x["orbit"]) for x in ranking),
        payload["status"],
    ),
    flush=True,
)
