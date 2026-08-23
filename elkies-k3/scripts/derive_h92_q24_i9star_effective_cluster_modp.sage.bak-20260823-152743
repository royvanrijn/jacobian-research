#!/usr/bin/env sage -python
"""
Derive the effective infinitely-near-point cluster for the H3 q24 divisor on
the exact D13 parent.

Status: ACTIVE_PROOF once all assertions pass.

The component-graph certificate matches the deterministic D13 root *graph* to
the geometric I9* resolution, but graph data alone does not fix the global
sign of the simple-root classes.  The independent D13 chamber certificate
fixes that sign: in the pinned root-adapted frame the effective simple roots
are the NEGATIVES of the first thirteen frame basis vectors.

Therefore, if

    D24 = O + P - 7 F + sum_i vr[i] * r_i

is the equation-frame decomposition produced by
close_h92_q8_q24_by_q6_translation.sage, then the coefficient of the effective
component E_i is -vr[i].  When the global -7F twist is represented at infinity
(as in probe_h92_q24_d12_degree2_rr_modp.sage), the finite I9* local divisor is

    - sum_i vr[i] * E_i.

A section represented in the singular affine chart must consequently have
resolved divisorial valuation at least vr[i] on E_i.

For a sequence of point blow-ups these valuation lower bounds translate into
multiplicity conditions at infinitely-near centres.  If a new exceptional E
is created at a centre through old components E_j, then the already-forced
old valuations contribute sum_j t(E_j), so the new centre needs only

    max(0, t(E) - sum_j t(E_j))

additional order.  Rank-two tangent cones are treated branchwise; if the two
branches would require different positive residual orders this script refuses
to collapse them into one point condition.

For the selected q24 orbit the derived nonzero cluster is expected to be the
long-chain centres C01,C02,C04,C06 with extra orders 2,2,2,3.  This statement
is derived and checked, not inserted as the definition.

Output:
  artifacts/local/elkies-k3/q24-i9star-effective-cluster-mod-<p>.json
"""

import argparse
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    home=Path.home()
    candidates += [
        home/"Documents"/"jacobian-research",
        home/"jacobian-research",
        home/"src"/"jacobian-research",
        home/"git"/"jacobian-research",
        home/"projects"/"jacobian-research",
    ]
    seen=set()
    for candidate in candidates:
        try:
            candidate=candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if ((candidate/"elkies-k3/scripts").is_dir()
                and (candidate/"artifacts/generated-results").is_dir()):
            return candidate
    raise SystemExit("Could not locate jacobian-research")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
PINNED=ROOT/"elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
GRAPH=LOCAL/f"q24-i9star-component-graph-mod-{args.prime}.json"
RESOLUTION=LOCAL/f"q24-i9star-resolution-mod-{args.prime}.json"
TRANSLATION=LOCAL/"q8-q24-physical-to-equation-translation.json"

for path in (PINNED,GRAPH,RESOLUTION,TRANSLATION):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

frame=load_gram(PINNED)
graph=json.loads(GRAPH.read_text())
resolution=json.loads(RESOLUTION.read_text())
translation=json.loads(TRANSLATION.read_text())

assert graph["status"]=="PASS_H3_Q24_AFFINE_D13_COMPONENT_GRAPH"
assert resolution["status"]=="PASS_EXPLICIT_MODP_I9STAR_D13_COMPONENT_RESOLUTION"
assert translation["status"]=="PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION"
assert frame.nrows()==17 and frame.det()==948
assert matrix(ZZ,graph["lattice_graph"]["root_gram"])==frame[:13,:13]

# ---------------------------------------------------------------------------
# Re-certify the effective sign in the pinned D13 chamber.
# ---------------------------------------------------------------------------
Q24_WITNESS=vector(ZZ,(
    0,5,0,1,2,1,2,2,2,2,4,8,2,0,-1,1,1,
))
U2=matrix(ZZ,((0,1),(1,0)))
ns=block_diagonal_matrix(U2,-frame)
q24=vector(ZZ,[12,2]+list(Q24_WITNESS))
effective_simple=tuple(
    vector(ZZ,[0,0]+[-ZZ(index==node) for index in range(17)])
    for node in range(13)
)
pairings=tuple(q24*ns*curve for curve in effective_simple)
assert pairings==(0,0,1,0,0,0,0,0,0,0,0,0,0)

vr=list(map(int,translation["q24_equation"]["vertical_root_coefficients"]))
assert vr==list(map(int,graph["lattice_graph"]["vertical_root_coefficients"]))
assert all(value>0 for value in vr)

print(
    "Q24D13EFFECTIVE_SIGN|effective_simple=-frame_basis|"
    f"q24_pairings={','.join(map(str,pairings))}|"
    "status=PASS_EFFECTIVE_SIGN_MINUS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Push the REQUIRED vanishing orders to geometric components.
# ---------------------------------------------------------------------------
profiles=[]
for raw_profile in graph["spinor_orientation"]["candidate_profiles"]:
    orientation=int(raw_profile["orientation"])
    mapping={str(k):str(v) for k,v in raw_profile["lattice_to_geometric"].items()}
    required={"F0":0}
    effective_divisor_coefficients={"F0":0}
    for i,value in enumerate(vr):
        component=mapping[str(i)]
        # r_i = -E_i, hence vr_i*r_i = -vr_i*E_i.
        effective_divisor_coefficients[component]=-int(value)
        required[component]=int(value)
    assert set(required)==set(graph["geometric_graph"]["affine_vertices"])
    profiles.append({
        "orientation":orientation,
        "lattice_to_geometric":mapping,
        "effective_divisor_coefficients_at_I9star":effective_divisor_coefficients,
        "required_divisorial_vanishing":required,
    })

# ---------------------------------------------------------------------------
# Convert divisorial thresholds into additional centre multiplicities.
# ---------------------------------------------------------------------------
chronology=graph["geometric_graph"]["chronology"]
for profile in profiles:
    required=profile["required_divisorial_vanishing"]
    centre_plan=[]
    unresolved=[]
    for step in chronology:
        center=str(step["center"])
        active=[str(v) for v in step["active_geometric_components"]]
        new=[str(v) for v in step["new_geometric_components"]]
        baseline=sum(int(required[v]) for v in active)
        targets=[int(required[v]) for v in new]
        if len(new)==1:
            extra=max(0,targets[0]-baseline)
            branch_extra={new[0]:extra}
        elif len(new)==2:
            raw_extra=[max(0,target-baseline) for target in targets]
            # A common point-order condition applies equally to both tangent
            # branches.  If positive residual requirements differ, actual
            # branch quotient data are required and we refuse to over-impose.
            if raw_extra[0]!=raw_extra[1] and max(raw_extra)>0:
                unresolved.append({
                    "center":center,
                    "active":active,
                    "new":new,
                    "baseline":baseline,
                    "targets":targets,
                    "branch_extra":dict(zip(new,raw_extra)),
                })
                extra=None
                branch_extra=dict(zip(new,raw_extra))
            else:
                extra=raw_extra[0]
                branch_extra=dict(zip(new,raw_extra))
        else:
            raise ArithmeticError(f"unexpected number of new components at {center}: {new}")
        centre_plan.append({
            "center":center,
            "active":active,
            "new":new,
            "baseline_required_from_active":int(baseline),
            "new_required":targets,
            "additional_point_order":None if extra is None else int(extra),
            "branch_additional_orders":branch_extra,
        })
    if unresolved:
        raise ArithmeticError(
            "effective q24 cluster has unresolved asymmetric positive split conditions: "
            + repr(unresolved)
        )
    profile["centre_plan"]=centre_plan

# The two spinor orientations differ only on C10a/C10b, and their thresholds
# are already dominated by C08.  Hence the derived infinitely-near point
# conditions must coincide.
def nonzero_plan(profile):
    return [
        (row["center"],int(row["additional_point_order"]))
        for row in profile["centre_plan"]
        if row["additional_point_order"] not in (None,0)
    ]

plans=[nonzero_plan(profile) for profile in profiles]
assert plans[0]==plans[1]
assert plans[0]==[("C01",2),("C02",2),("C04",2),("C06",3)]

print(
    "Q24D13EFFECTIVE_CLUSTER|"
    + "|".join(f"{name}={order}" for name,order in plans[0])
    + "|orientations_agree=1|status=PASS_CLUSTER_2_2_2_3",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-i9star-effective-cluster-modp.v1",
    "status":"PASS_H3_Q24_EFFECTIVE_I9STAR_CLUSTER",
    "prime":int(args.prime),
    "effective_root_sign":-1,
    "q24_effective_simple_pairings":list(map(int,pairings)),
    "global_fibre_twist_location":"infinity",
    "local_I9star_fibre_twist":0,
    "raw_vertical_root_coefficients":vr,
    "profiles":profiles,
    "common_nonzero_centre_plan":[
        {"center":name,"additional_point_order":order}
        for name,order in plans[0]
    ],
    "boundary":(
        "This certifies the effective sign and converts the q24 vertical D13 "
        "cycle into the infinitely-near point multiplicity cluster.  It does "
        "not yet calculate the linear condition matrix on the ten-dimensional "
        "post-collision global space."
    ),
    "next":(
        "Pull the ten post-collision numerator functions successively through "
        "the C01->C02->C04->C06 blow-up charts, impose additional orders "
        "2,2,2,3 in the local surface rings, and test whether the resulting "
        "kernel has dimension two and compiles to a degree-four D12 Jacobian."
    ),
}
OUT=(args.output.resolve() if args.output else
     LOCAL/f"q24-i9star-effective-cluster-mod-{args.prime}.json")
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24D13EFFECTIVE_RESULT|sign=-1|cluster=C01:2,C02:2,C04:2,C06:3|"
    "status=PASS_H3_Q24_EFFECTIVE_I9STAR_CLUSTER",
    flush=True,
)
