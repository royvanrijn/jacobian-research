#!/usr/bin/env sage -python
"""
Build an exact backward lift manifest for the preferred H3 corridor

  D13 --q24/85--> D12 --q6/42--> A11 --...--> rootless/MW17 -> pinned R17.

The purpose is equation lifting, not neighbor search.

For every selected neighbor arrow, recover in the PARENT NS marking:

  * the exact child fibre F_child = row_0(step transition);
  * the exact child zero O_child = -row_0 + row_1;
  * the exact MW projection of the horizontal parent section;
  * its Shioda height, exact local root correction, and P.O;
  * a closest root-lattice representative of that section class;
  * the exact abstract decomposition

        F_child = O_parent + P + V + k F_parent,

    where V is root-lattice vertical and k = q/2 - P.O.

This is the structural data used by the original resolved-RR compilers.
It deliberately does not claim that V is already expressed in the effective
geometric fibre-component basis; that is the next equation-specific conversion.

The output is printed in REVERSE order (R17/rootless backwards to D13), while
the JSON also retains forward order.
"""

import argparse, json
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, identity_matrix,
    matrix, vector
)

ROOT=Path(__file__).resolve().parents[2]
GEN=ROOT/"artifacts/generated-results"
LOCAL=ROOT/"artifacts/local/elkies-k3"
R17CERT=LOCAL/"q24-equation-d13-to-pinned-r17.json"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

U2=matrix(ZZ,((0,1),(1,0)))

STEPS=[
    # parent, child, artifact, orbit, q
    ("D13/MW4","D12/MW5",
     "elkies-k3-h3-q6-q8-d13-q24-degree2.json",85,24),
    ("D12/MW5","A11/MW6",
     "elkies-k3-h3-d12-o85-q6-degree2.json",42,6),
    ("A11/MW6","2A5/MW7",
     "elkies-k3-h3-a11-middle-q8-degree2.json",922,8),
    ("2A5/MW7","3A3/MW8",
     "elkies-k3-h3-a5a5-c2-q4-degree2.json",472,4),
    ("3A3/MW8","A3+2A2/MW10",
     "elkies-k3-h3-a3x3-q4-degree2.json",323,4),
    ("A3+2A2/MW10","5A1/MW12",
     "elkies-k3-h3-mw10-a3a2a2-q4-degree2.json",207,4),
    ("5A1/MW12","4A1/MW13",
     "elkies-k3-h3-mw12-5a1-q4-degree2-first-hit.json",52,4),
    ("4A1/MW13","3A1/MW14",
     "elkies-k3-h3-mw13-4a1-q4-degree2-first-hit.json",114,4),
    ("3A1/MW14","2A1/MW15",
     "elkies-k3-h3-mw14-3a1-q4-degree2-first-hit.json",498,4),
    ("2A1/MW15","A1/MW16",
     "elkies-k3-h3-mw15-2a1-q4-degree2-first-hit.json",981,4),
    ("A1/MW16","rootless/MW17",
     "elkies-k3-h3-mw16-a1-q6-degree2-cap10000-stream-chunk001.json",2247,6),
]

def load_gram(path):
    return matrix(ZZ,[
        [ZZ(x) for x in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])

def ns(frame):
    return block_diagonal_matrix(U2,-frame)

def root_rank_from_record(rec):
    # Parent frame is root-adapted: derive its actual root rank intrinsically
    # from the artifact frame using the fact that the leading block is ADE.
    # We instead use witness/mw_projection lengths when available; this is
    # exact and avoids re-running root enumeration.
    if "mw_projection" in rec:
        return 17-len(rec["mw_projection"])
    raise KeyError("neighbor record has no mw_projection")

def step_transition(rec):
    child=matrix(ZZ,rec["child_root_adapted_frame"])
    A=matrix(ZZ,rec["child_root_adapted_basis"])
    N=matrix(ZZ,rec["neighbor_basis"])
    T=block_diagonal_matrix(identity_matrix(ZZ,2),A)*N
    return child,T

def exact_section_profile(frame,rec,q):
    r=root_rank_from_record(rec)
    assert 0 < r <= 17
    G=matrix(ZZ,frame)
    R=G[:r,:r]
    C=G[:r,r:]
    Tail=G[r:,r:]
    H=Tail-C.transpose()*R.inverse()*C

    z=vector(ZZ,rec["mw_projection"])
    witness=vector(ZZ,rec["witness"])
    assert len(z)==17-r
    assert vector(ZZ,witness[r:])==z

    h=QQ(z*H*z)

    base=vector(ZZ,[0]*r+list(z))
    pair=vector(QQ,base*G[:,:r])
    dual=pair*R.inverse()

    # Exact closest-vector enumeration in the root lattice.  The first shell
    # gives the discriminant-class correction.  Keep all vectors in that shell
    # and choose the representative making the vertical remainder simplest.
    L=IntegralLattice(R)
    it=L.enumerate_close_vectors(-dual)

    shell=[]
    min_norm=None
    for _ in range(200000):
        shift=vector(ZZ,next(it))
        pf=base+vector(ZZ,list(shift)+[0]*(17-r))
        total=QQ(pf*G*pf)
        if min_norm is None:
            min_norm=total
        if total>min_norm:
            break
        if total<min_norm:
            min_norm=total
            shell=[]
        if total==min_norm:
            vr=vector(ZZ,witness[:r])-vector(ZZ,pf[:r])
            shell.append((sum(abs(int(x)) for x in vr),
                          sum(bool(x) for x in vr),
                          tuple(map(int,shift)),pf,vr))

    if not shell:
        raise ArithmeticError("empty closest-vector shell")

    shell.sort(key=lambda row:(row[0],row[1],row[2]))
    l1,support,unused_shift,pframe,vroot=shell[0]

    corr=QQ(min_norm-h)
    po=QQ((h+corr-4)/2)
    if po not in ZZ or po<0:
        raise ArithmeticError(
            f"invalid section profile h={h}, corr={corr}, P.O={po}"
        )
    po=ZZ(po)

    # In NS=U-G, a section with P.O=po has U coordinates [po+1,1].
    P=vector(ZZ,[po+1,1]+list(pframe))
    F=vector(ZZ,[1,0]+[0]*17)
    O=vector(ZZ,[-1,1]+[0]*17)
    D=vector(ZZ,rec["fiber"])

    assert D==vector(ZZ,[q//2,2]+list(witness))
    Gns=ns(G)
    assert P*Gns*P==-2
    assert P*Gns*F==1
    assert P*Gns*O==po
    assert D*Gns*D==0
    assert D*Gns*F==2

    k=ZZ(q//2-po)
    V=vector(ZZ,[0,0]+list(vroot)+[0]*(17-r))
    assert O+P+V+k*F==D

    return {
        "root_rank":r,
        "mw_projection":list(map(int,z)),
        "height":str(h),
        "local_correction":str(corr),
        "P_dot_O":int(po),
        "fibre_twist":int(k),
        "section_frame":list(map(int,pframe)),
        "section_class":list(map(int,P)),
        "vertical_root_coefficients":list(map(int,vroot)),
        "vertical_root_L1":int(l1),
        "vertical_root_support":int(support),
        "closest_shell_size":len(shell),
    }

if not R17CERT.exists():
    raise SystemExit(f"missing R17 certificate: {R17CERT}")
cert=json.loads(R17CERT.read_text())
assert cert["status"]=="PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"

forward=[]

for parent_name,child_name,artifact_name,orbit,q in STEPS:
    path=GEN/artifact_name
    if not path.exists():
        raise SystemExit(f"missing {path}")
    data=json.loads(path.read_text())
    assert data["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"

    frame=load_gram(ROOT/data["frame"])
    hits=[r for r in data["neighbors"] if int(r["orbit_index"])==orbit]
    assert len(hits)==1
    rec=hits[0]
    assert int(rec["q"])==q
    assert int(rec["old_fiber_degree"])==2

    child,T=step_transition(rec)
    assert T*ns(frame)*T.transpose()==ns(child)
    assert abs(ZZ(T.det()))==1

    Fchild=vector(ZZ,T.row(0))
    Schild=vector(ZZ,T.row(1))
    Ochild=Schild-Fchild
    Fold=vector(ZZ,[1,0]+[0]*17)
    Oold=vector(ZZ,[-1,1]+[0]*17)
    Gns=ns(frame)

    assert Fchild==vector(ZZ,rec["fiber"])
    assert Fchild*Gns*Fchild==0
    assert Fchild*Gns*Fold==2
    assert Ochild*Gns*Ochild==-2
    assert Ochild*Gns*Fchild==1

    prof=exact_section_profile(frame,rec,q)

    item={
        "parent":parent_name,
        "child":child_name,
        "q":q,
        "orbit":orbit,
        "artifact":str(path.relative_to(ROOT)),
        "child_root_data":list(map(int,rec["child_root_data"])),
        "child_mw_rank":int(rec["child_mw_rank"]),
        "new_fibre_in_parent":list(map(int,Fchild)),
        "new_zero_in_parent":list(map(int,Ochild)),
        "new_zero_old_fibre_degree":int(Ochild*Gns*Fold),
        "new_zero_old_O_intersection":int(Ochild*Gns*Oold),
        "horizontal":prof,
        "transition":[list(map(int,row)) for row in T.rows()],
    }
    forward.append(item)

# Hard regression for the corrected historical orbit42 profile.
o42=next(x for x in forward if x["orbit"]==42)
assert o42["horizontal"]["height"]=="7"
assert o42["horizontal"]["local_correction"]=="3"
assert o42["horizontal"]["P_dot_O"]==3
assert o42["horizontal"]["fibre_twist"]==0

print(
    "H3BACKWARD_REGRESSION|"
    "orbit42=1|height=7|corr=3|PdotO=3|fiber_twist=0|status=PASS",
    flush=True,
)

# Print from the pinned endpoint backwards, because this is the intended
# equation-lift planning direction.
for rev_index,item in enumerate(reversed(forward),1):
    h=item["horizontal"]
    print(
        "H3BACKWARD_STEP|"
        f"reverse={rev_index}|parent={item['parent']}|child={item['child']}|"
        f"q={item['q']}|orbit={item['orbit']}|"
        f"mw={','.join(map(str,h['mw_projection']))}|"
        f"height={h['height']}|corr={h['local_correction']}|"
        f"PdotO={h['P_dot_O']}|fiber_twist={h['fibre_twist']}|"
        f"vertical_L1={h['vertical_root_L1']}|"
        f"vertical_support={h['vertical_root_support']}|"
        f"child_zero_old_degree={item['new_zero_old_fibre_degree']}|"
        "status=PASS_EXACT_PARENT_STRUCTURE",
        flush=True,
    )

# A compact forward lift classification.
for item in forward:
    h=item["horizontal"]
    # This is not a theorem about runtime; it is a structural heuristic.
    if h["P_dot_O"]==0:
        route="NODE_INTERPOLATION_OR_POLYNOMIAL_SECTION"
    elif h["P_dot_O"]<=1:
        route="LOW_DENOMINATOR_SECTION"
    elif h["P_dot_O"]<=3:
        route="RESOLVED_RR_OR_SMALL_DENOMINATOR"
    else:
        route="RESOLVED_RR_PREFERRED"
    item["equation_lift_hint"]=route

payload={
    "schema":"elkies-k3.h3-r17-backward-exact-lift-manifest.v1",
    "status":"PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST",
    "source_certificate":str(R17CERT.relative_to(ROOT)),
    "forward_steps":forward,
    "reverse_order":[
        {
            "parent":x["parent"],
            "child":x["child"],
            "q":x["q"],
            "orbit":x["orbit"],
        }
        for x in reversed(forward)
    ],
    "proof_boundary":(
        "Exact integral NS/lattice structure for equation lifting.  The child "
        "fibre, child zero, horizontal MW class, local correction, P.O, and "
        "abstract root-lattice vertical decomposition are exact.  Effective "
        "geometric fibre-component realization and rational-function RR "
        "construction remain equation-specific steps."
    ),
}

OUT=args.output.resolve() if args.output else LOCAL/"h3-r17-backward-exact-lift-manifest.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "H3BACKWARD_RESULT|steps={}|status={}".format(
        len(forward),payload["status"]
    ),
    flush=True,
)
