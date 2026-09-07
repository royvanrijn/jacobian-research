#!/usr/bin/env sage -python
"""
Identify the 13 resolved geometric I9* components with the current equation
D13 root basis, and derive the exact local q24 divisor/order targets.

Inputs:
  * passing explicit I9* resolution mod p=100003 with 12 centers / 13
    irreducible exceptional curves;
  * current passing q24 physical->equation translation.

Outputs:
  * exact D13 resolution graph;
  * graph isomorphism(s) equation-root <-> geometric exceptional component;
  * q24-deleted leaf (must be geometric E12);
  * explicit D12 affine bisection Theta and its geometric component hit
    (expected E03);
  * the q24 divisor representative coefficients on the actual old fibre;
  * generic L(D24) minimum orders and the stronger special-Theta orders on
    every geometric component.

Only the unavoidable D13 spinor-leaf swap remains ambiguous.
"""

import argparse
import contextlib
import io
import json
import sys
from collections import deque
from pathlib import Path

from sage.all import (
    PolynomialRing, QQ, ZZ, matrix, pari, sage_eval, vector
)


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
    buf=io.StringIO()
    try:
        sys.argv=[str(path)]
        with contextlib.redirect_stdout(buf):
            exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope


def add_edge(graph,a,b):
    if a==b:
        return
    graph.setdefault(a,set()).add(b)
    graph.setdefault(b,set()).add(a)


def remove_edge(graph,a,b):
    graph.setdefault(a,set()).discard(b)
    graph.setdefault(b,set()).discard(a)


def path_between(graph,start,end):
    todo=deque([(start,[start])])
    seen={start}
    while todo:
        node,path=todo.popleft()
        if node==end:
            return path
        for nxt in graph[node]:
            if nxt in seen:
                continue
            seen.add(nxt)
            todo.append((nxt,path+[nxt]))
    raise ArithmeticError("graph path missing")


def d13_shape(graph):
    branch=[v for v in graph if len(graph[v])==3]
    assert len(branch)==1,branch
    branch=branch[0]
    leaves=[v for v in graph if len(graph[v])==1]
    assert len(leaves)==3,leaves
    spin=[v for v in leaves if branch in graph[v]]
    assert len(spin)==2,spin
    far=next(v for v in leaves if v not in spin)
    chain=path_between(graph,far,branch)
    assert len(chain)==11,chain
    return far,branch,tuple(chain),tuple(sorted(spin))


def highest_root(cartan):
    qf=pari(cartan).qfminim(2)
    half=matrix(ZZ,qf[2]).transpose().rows()
    roots=[vector(ZZ,r) for r in half]+[-vector(ZZ,r) for r in half]
    positive=[r for r in roots if all(v>=0 for v in r)]
    assert positive
    return max(positive,key=lambda r:sum(r))


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
RES=LOCAL/f"q24-i9star-resolution-mod-{args.prime}.json"
CLOSE=ROOT/"elkies-k3/scripts/close_h92_q8_q24_by_q6_translation.sage"

for path in (RES,CLOSE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

resolution=json.loads(RES.read_text())
assert resolution["status"]=="PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS"
assert resolution["actual_blowup_centers"]==12
assert resolution["actual_exceptional_irreducible_components"]==13

cl=run_scope(CLOSE)
need=("ns","F8eq","O8","D24eq","root_source","vertical_coeffs")
missing=[x for x in need if x not in cl]
if missing:
    raise SystemExit("close script missing: "+",".join(missing))

ns=cl["ns"]
F8=vector(ZZ,cl["F8eq"])
O8=vector(ZZ,cl["O8"])
D=vector(ZZ,cl["D24eq"])
roots=[vector(ZZ,r) for r in cl["root_source"]]
vertical_coeffs=vector(ZZ,cl["vertical_coeffs"])
assert len(roots)==13 and len(vertical_coeffs)==14

# ---------------------------------------------------------------------------
# 1. Reconstruct final geometric exceptional graph from blow-up incidence.
# ---------------------------------------------------------------------------

centers=sorted(resolution["centers"],key=lambda r:int(r["label"][1:]))
split_parts={}
split_records={}

graph={"F0":set()}

# Polynomial ring only for evaluating homogeneous tangent factors at stored
# projective child directions.
Fp = None
from sage.all import GF
Fp=GF(ZZ(args.prime))
S=PolynomialRing(Fp,names=("u","x","y"))
u,x,y=S.gens()

def active_nodes(label):
    if label not in split_parts:
        return [label]
    rec=split_records[label]
    # In the present resolution a split center has one remaining singular
    # child. Evaluate each tangent line at that projective direction; all
    # branches through that direction are active at the child.
    children=rec["children"]
    assert len(children)==1
    direction=tuple(Fp(v) for v in children[0]["direction"])
    out=[]
    for node,factor_text in split_parts[label]:
        factor=S(sage_eval(factor_text,locals={"u":u,"x":x,"y":y}))
        if factor(*direction)==0:
            out.append(node)
    assert out
    return out

for rec in centers:
    label=rec["label"]
    parts=int(rec["new_exceptional_irreducible_components"])
    factors=rec["tangent_factors"]

    if parts==1:
        new_nodes=[label]
        graph.setdefault(label,set())
    else:
        assert parts==len(factors)==2
        new_nodes=[f"{label}a",f"{label}b"]
        for node in new_nodes:
            graph.setdefault(node,set())
        split_parts[label]=[
            (node,factors[i]["factor"]) for i,node in enumerate(new_nodes)
        ]
        split_records[label]=rec
        # The two projective tangent lines meet before the next blow-up.
        add_edge(graph,new_nodes[0],new_nodes[1])

    active=[]
    for old in rec["active_components"]:
        if old=="F0":
            active.append("F0")
        elif old in split_parts:
            active.extend(active_nodes(old))
        else:
            active.append(old)

    # Blowing up a common point separates all old branches at that point.
    for i in range(len(active)):
        for j in range(i+1,len(active)):
            remove_edge(graph,active[i],active[j])

    # A single exceptional component meets every active branch.  At a split
    # tangent cone, both exceptional lines pass through the active tangent
    # direction in this I9* tree; the following singular blow-up separates
    # the triple point when necessary.
    for n in new_nodes:
        for a in active:
            add_edge(graph,n,a)

# Remove F0 to get the D13 nonidentity root graph.
exceptional={k:set(v) for k,v in graph.items() if k!="F0"}
for k in exceptional:
    exceptional[k].discard("F0")

assert len(exceptional)==13
assert sum(len(v) for v in exceptional.values())//2==12
geo_far,geo_branch,geo_chain,geo_spin=d13_shape(exceptional)

print(
    "Q24I9GRAPH|"
    f"far={geo_far}|branch={geo_branch}|"
    f"chain={'-'.join(geo_chain)}|spinors={','.join(geo_spin)}|"
    "vertices=13|edges=12|status=PASS_D13_GRAPH",
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Current equation D13 root graph and graph isomorphism.
# ---------------------------------------------------------------------------

Gns=matrix(ZZ,[[int(a*ns*b) for b in roots] for a in roots])
assert all(Gns[i,i]==-2 for i in range(13))

eqgraph={f"R{i+1}":set() for i in range(13)}
for i in range(13):
    for j in range(i+1,13):
        assert Gns[i,j] in (0,1)
        if Gns[i,j]==1:
            add_edge(eqgraph,f"R{i+1}",f"R{j+1}")

eq_far,eq_branch,eq_chain,eq_spin=d13_shape(eqgraph)

# Effective components are C_i=-root_i in the positive-frame convention.
C=[-r for r in roots]
Dpair=[int(D*ns*c) for c in C]
assert sorted(Dpair)==[0]*12+[1]
deleted_index=Dpair.index(1)
deleted=f"R{deleted_index+1}"
assert deleted==eq_far,(deleted,eq_far,Dpair)

# Chain identification is forced. Only the two spinor leaves can swap.
base_map={eq_chain[i]:geo_chain[i] for i in range(11)}
maps=[]
for swap in (False,True):
    m=dict(base_map)
    gs=(geo_spin[1],geo_spin[0]) if swap else geo_spin
    m[eq_spin[0]]=gs[0]
    m[eq_spin[1]]=gs[1]
    assert len(m)==13
    maps.append(m)

assert maps[0][deleted]==geo_far

print(
    "Q24I9MAP|"
    f"deleted_root={deleted}|geometric_deleted={maps[0][deleted]}|"
    f"eq_branch={eq_branch}|geo_branch={geo_branch}|"
    f"spinor_ambiguity=2|status=PASS_GRAPH_ISOMETRIES",
    flush=True,
)

for mi,m in enumerate(maps,1):
    print(
        "Q24I9MAP_DETAIL|"
        f"map={mi}|"
        +";".join(f"{r}->{m[r]}" for r in sorted(m,key=lambda s:int(s[1:])))
        +"|status=PASS",
        flush=True,
    )

# ---------------------------------------------------------------------------
# 3. Explicit old affine component and special D12 fibre member.
# ---------------------------------------------------------------------------

Cartan=-Gns
high13=highest_root(Cartan)
affine=F8+sum((high13[i]*roots[i] for i in range(13)),vector(ZZ,[0]*len(F8)))
assert affine*ns*affine==-2
assert affine*ns*F8==0
assert affine*ns*O8==1
assert D*ns*affine==1

keep=[i for i in range(13) if i!=deleted_index]
C12=Cartan.matrix_from_rows_and_columns(keep,keep)
qf12=pari(C12).qfminim(2)
assert (C12.rank(),ZZ(qf12[0]),abs(ZZ(C12.det())))==(12,264,4)
high12local=highest_root(C12)
m12=[ZZ(0)]*13
for i,v in zip(keep,high12local):
    m12[i]=ZZ(v)

# D = Theta + sum m_i*C_i, C_i=-root_i.
Theta=D+sum((m12[i]*roots[i] for i in range(13)),vector(ZZ,[0]*len(D)))
assert Theta*ns*Theta==-2
assert Theta*ns*F8==2
assert Theta*ns*D==0
assert Theta*ns*affine==0

theta_hits=[
    (i,int(Theta*ns*C[i]))
    for i in range(13) if Theta*ns*C[i]
]
assert len(theta_hits)==1 and theta_hits[0][1]==1
theta_hit_index=theta_hits[0][0]
theta_hit=f"R{theta_hit_index+1}"

print(
    "Q24D12GEOM_THETA|"
    f"deleted={deleted}->{maps[0][deleted]}|"
    f"theta_hit={theta_hit}->{maps[0][theta_hit]}|"
    f"square={Theta*ns*Theta}|old_degree={Theta*ns*F8}|"
    f"new_degree={Theta*ns*D}|old_affine={Theta*ns*affine}|"
    f"D12_highest={','.join(map(str,m12))}|"
    "status=PASS_EXPLICIT_D12_FIBRE",
    flush=True,
)

# ---------------------------------------------------------------------------
# 4. Exact local line-bundle orders on actual resolved components.
#
# Current representative:
#     Drep = O + P + vf*F + sum vr_i * root_i
#
# with effective C_i=-root_i and
#     F = C_aff + sum high13_i*C_i.
#
# Therefore divisor coefficients along old-fibre components are
#     d_aff = vf
#     d_i   = vf*high13_i - vr_i.
#
# Any f in L(Drep) needs ord_C(f) >= -d_C.
# For the special D12 fibre generator f_Theta,
#     div(f_Theta)+Drep = Theta + sum m12_i*C_i,
# so its exact generic component orders are
#     m12_i-d_i,  and -d_aff on C_aff.
# ---------------------------------------------------------------------------

vf=ZZ(vertical_coeffs[0])
vr=vector(ZZ,vertical_coeffs[1:])
assert vf==-7
assert len(vr)==13

d_aff=vf
d_root=[vf*high13[i]-vr[i] for i in range(13)]
generic_min_aff=-d_aff
generic_min_root=[-v for v in d_root]
theta_order_aff=-d_aff
theta_order_root=[m12[i]-d_root[i] for i in range(13)]

assert all(v>=0 for v in generic_min_root)
assert generic_min_aff>=0
assert all(theta_order_root[i]>=generic_min_root[i] for i in range(13))

print(
    "Q24I9LOCAL_AFFINE|"
    f"component=F0|Drep_coeff={d_aff}|generic_min_order={generic_min_aff}|"
    f"theta_order={theta_order_aff}|status=PASS",
    flush=True,
)

for mi,mapping in enumerate(maps,1):
    print(f"Q24I9LOCAL_MAP|map={mi}|status=BEGIN",flush=True)
    rows=[]
    for i in range(13):
        eq=f"R{i+1}"
        geo=mapping[eq]
        rows.append({
            "equation_root":eq,
            "geometric_component":geo,
            "D_pairing":Dpair[i],
            "D13_highest_multiplicity":int(high13[i]),
            "D12_special_multiplicity":int(m12[i]),
            "Drep_component_coefficient":int(d_root[i]),
            "generic_min_order":int(generic_min_root[i]),
            "theta_exact_order":int(theta_order_root[i]),
        })
        print(
            "Q24I9LOCAL|"
            f"map={mi}|eq={eq}|geo={geo}|Dpair={Dpair[i]}|"
            f"mD13={high13[i]}|mD12={m12[i]}|"
            f"Drep={d_root[i]}|minord={generic_min_root[i]}|"
            f"thetaord={theta_order_root[i]}|status=PASS",
            flush=True,
        )
    print(f"Q24I9LOCAL_MAP|map={mi}|status=END",flush=True)

# JSON output.
map_payload=[]
for mapping in maps:
    map_payload.append({
        "equation_to_geometric":mapping,
        "local_orders":[
            {
                "equation_root":f"R{i+1}",
                "geometric_component":mapping[f"R{i+1}"],
                "D_pairing":Dpair[i],
                "D13_highest_multiplicity":int(high13[i]),
                "D12_special_multiplicity":int(m12[i]),
                "Drep_component_coefficient":int(d_root[i]),
                "generic_min_order":int(generic_min_root[i]),
                "theta_exact_order":int(theta_order_root[i]),
            }
            for i in range(13)
        ],
    })

payload={
    "schema":"elkies-k3.h3-q24-i9star-geometric-component-map.v1",
    "status":"PASS_EXACT_Q24_I9STAR_GEOMETRIC_COMPONENT_MAP",
    "prime":int(args.prime),
    "geometric_graph":{
        k:sorted(v) for k,v in sorted(exceptional.items())
    },
    "geometric_D13":{
        "far_leaf":geo_far,
        "branch":geo_branch,
        "long_chain":list(geo_chain),
        "spinor_leaves":list(geo_spin),
    },
    "equation_D13":{
        "far_leaf":eq_far,
        "branch":eq_branch,
        "long_chain":list(eq_chain),
        "spinor_leaves":list(eq_spin),
        "q24_pairings_with_effective_components":Dpair,
        "q24_deleted_root":deleted,
    },
    "maps":map_payload,
    "old_affine_component_class":list(map(int,affine)),
    "theta":{
        "class":list(map(int,Theta)),
        "square":-2,
        "old_fibre_degree":2,
        "new_fibre_degree":0,
        "old_affine_intersection":0,
        "equation_root_hit":theta_hit,
        "geometric_hit_under_map1":maps[0][theta_hit],
        "D12_highest_multiplicities":list(map(int,m12)),
    },
    "line_bundle_orders":{
        "vertical_fibre_coefficient":int(vf),
        "vertical_root_coefficients":list(map(int,vr)),
        "affine":{
            "Drep_component_coefficient":int(d_aff),
            "generic_min_order":int(generic_min_aff),
            "theta_exact_order":int(theta_order_aff),
        },
    },
    "next":(
        "Pull the 10-dimensional post-collision q24 basis through the stored "
        "resolution charts. For each of the two spinor assignments, impose "
        "the generic_min_order divisibility along the actual exceptional "
        "components. The correct resolved module should have rank 8 and leave "
        "kernel dimension 2. Independently impose theta_exact_order to isolate "
        "the one-dimensional special D12 fibre generator."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q24-i9star-geometric-map-mod-{args.prime}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24I9GEOM_RESULT|"
    f"deleted={deleted}->{maps[0][deleted]}|"
    f"theta_hit={theta_hit}->{maps[0][theta_hit]}|"
    f"spinor_maps=2|status={payload['status']}",
    flush=True,
)
