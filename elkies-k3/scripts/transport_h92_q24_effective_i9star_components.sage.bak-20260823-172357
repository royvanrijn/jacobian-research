#!/usr/bin/env sage -python
"""
Transport the ACTUAL EFFECTIVE pinned D13 fibre components to the current
canonical q8 equation, then identify them with the explicit resolved I9*
components.

This avoids the invalid assumption that the deterministic D13 root basis is
itself the effective fibre-component basis.

Pipeline:
  pinned orbit-85 D13 effective components
    -> current effective-zero audit transport
    -> raw physical q8
    -> current q6 Weyl
    -> current Eichler translation tau
    -> canonical equation D13.

Then:
  * identify which transported component is the current Weierstrass identity;
  * map the full affine-D13 graph to the resolved chart graph (two spinor
    choices only);
  * certify q24 pairings still equal R3=1, affine=1, others=0;
  * transport the exact D12 special fibre and Theta;
  * solve the CURRENT q24 divisor representative on the actual effective
    components;
  * emit generic L(D24) minimum orders and exact Theta-member orders on every
    resolved component.
"""

import argparse
import contextlib
import io
import json
import sys
from collections import deque
from pathlib import Path

from sage.all import (
    GF, PolynomialRing, QQ, ZZ, matrix, pari, sage_eval, vector
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


def run_scope(path, argv=()):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    buf=io.StringIO()
    try:
        sys.argv=[str(path)]+list(argv)
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
    raise ArithmeticError("no graph path")


def finite_d13_shape(full_graph, identity):
    graph={
        k:set(v)-{identity}
        for k,v in full_graph.items()
        if k!=identity
    }
    assert len(graph)==13
    assert sum(len(v) for v in graph.values())//2==12
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
    return graph,far,branch,tuple(chain),tuple(sorted(spin))


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
AUDIT=ROOT/"elkies-k3/scripts/audit_h92_q8_q24_effective_zero_choices.sage"
CLOSE=ROOT/"elkies-k3/scripts/close_h92_q8_q24_by_q6_translation.sage"
RES=LOCAL/f"q24-i9star-resolution-mod-{args.prime}.json"
TMP=LOCAL/"q24-effective-components-temp.json"

for p in (AUDIT,CLOSE,RES):
    if not p.exists():
        raise SystemExit(f"Missing prerequisite: {p}")

print("Q24EFFCOMP|stage=load_current_transport|status=START",flush=True)
qa=run_scope(AUDIT,("--output",str(TMP)))
cl=run_scope(CLOSE)
print("Q24EFFCOMP|stage=load_current_transport|status=PASS",flush=True)

resolution=json.loads(RES.read_text())
assert resolution["status"]=="PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS"

need_qa=(
    "Gpinned","Bpinned_to_simple","q6_simple_ns","simple_root_classes",
    "bridge","Bsimple","source_ns","actual_roots","physical_reflections",
    "reflect","F_actual",
)
need_cl=(
    "ns","F8eq","O8","D24eq","P24","weyl_transport","tau",
)
mqa=[x for x in need_qa if x not in qa]
mcl=[x for x in need_cl if x not in cl]
if mqa or mcl:
    raise SystemExit(f"missing transport variables audit={mqa} close={mcl}")

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
F8=vector(ZZ,cl["F8eq"])
O8=vector(ZZ,cl["O8"])
D=vector(ZZ,cl["D24eq"])
P24=vector(ZZ,cl["P24"])
weyl_transport=cl["weyl_transport"]
tau=cl["tau"]

Fpin=vector(ZZ,[1,0]+[0]*17)
Opin=vector(ZZ,[-1,1]+[0]*17)

# Pinned D13 highest root and effective components.
Rpin=Gpinned[2:,2:]
# Gpinned is U + (-positive frame), so the positive D13 Cartan is pinned[:13].
CartanPin=-Rpin[:13,:13]
high13=highest_root(CartanPin)
assert tuple(high13)==(2,2,1,1,2,1,2,2,2,2,2,2,2)

pinned={}
for i in range(13):
    pinned[f"R{i+1}"]=vector(
        ZZ,[0,0]+[-ZZ(i==j) for j in range(17)]
    )
pinned["A0"]=Fpin+vector(ZZ,[0,0]+list(high13)+[0]*4)

for name,c in pinned.items():
    assert c*Gpinned*c==-2
    assert c*Gpinned*Fpin==0


def pinned_to_raw(Cpin):
    C=vector(ZZ,Cpin)*Bpinned_to_simple
    for i,unused_pairing in reversed(bridge):
        C=reflect(C,q6_simple_ns,simple_root_classes[i])
    C=vector(ZZ,C*Bsimple)
    for i,unused_pairing in physical_reflections:
        C=reflect(C,source_ns,actual_roots[i])
    assert C*source_ns*C==-2
    assert C*source_ns*Fraw==0
    return C


def raw_to_equation(Craw):
    C=vector(ZZ,weyl_transport(vector(ZZ,Craw)))
    C=vector(ZZ,tau(C))
    assert C*ns*C==-2
    assert C*ns*F8==0
    return C


effective={
    name:raw_to_equation(pinned_to_raw(c))
    for name,c in pinned.items()
}

# Verify full fibre relation survives transport.
fibre_sum=effective["A0"]
for i in range(13):
    fibre_sum += high13[i]*effective[f"R{i+1}"]
assert fibre_sum==F8

# Transported pinned q24 incidence pattern must survive exactly.
Dpairs={name:int(D*ns*c) for name,c in effective.items()}
assert [Dpairs[f"R{i+1}"] for i in range(13)]==[0,0,1,0,0,0,0,0,0,0,0,0,0]
assert Dpairs["A0"]==1

# Current equation zero selects the identity component among these 14.
Opairs={name:int(O8*ns*c) for name,c in effective.items()}
identity_hits=[name for name,v in Opairs.items() if v==1]
assert len(identity_hits)==1,(identity_hits,Opairs)
identity=identity_hits[0]
assert all(v in (0,1) for v in Opairs.values())

print(
    "Q24EFFCOMP_IDENTITY|"
    f"current_zero_meets={identity}|"
    f"pinned_affine_maps_to_identity={int(identity=='A0')}|"
    f"q24_degree1={','.join(name for name,v in Dpairs.items() if v==1)}|"
    "status=PASS_EFFECTIVE_COMPONENT_TRANSPORT",
    flush=True,
)

# Full transported affine-D13 graph.
tgraph={name:set() for name in effective}
names=list(effective)
for i,a in enumerate(names):
    for b in names[i+1:]:
        v=int(effective[a]*ns*effective[b])
        assert v in (0,1),(a,b,v)
        if v==1:
            add_edge(tgraph,a,b)
assert sum(len(v) for v in tgraph.values())//2==13

tfinite,tfar,tbranch,tchain,tspin=finite_d13_shape(tgraph,identity)

print(
    "Q24EFFCOMP_GRAPH|"
    f"identity={identity}|far={tfar}|branch={tbranch}|"
    f"chain={'-'.join(tchain)}|spinors={','.join(tspin)}|"
    "vertices=14|edges=13|status=PASS_AFFINE_D13",
    flush=True,
)

# ---------------------------------------------------------------------------
# Reconstruct explicit resolved geometric graph, including F0.
# ---------------------------------------------------------------------------

Fp=GF(ZZ(args.prime))
S=PolynomialRing(Fp,names=("u","x","y"))
u,x,y=S.gens()

ggraph={"F0":set()}
split_parts={}
split_records={}

def active_split_nodes(label):
    rec=split_records[label]
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

for rec in sorted(resolution["centers"],key=lambda r:int(r["label"][1:])):
    label=rec["label"]
    parts=int(rec["new_exceptional_irreducible_components"])
    factors=rec["tangent_factors"]
    if parts==1:
        new_nodes=[label]
        ggraph.setdefault(label,set())
    else:
        assert parts==2
        new_nodes=[label+"a",label+"b"]
        for n in new_nodes:
            ggraph.setdefault(n,set())
        split_parts[label]=[
            (new_nodes[i],factors[i]["factor"]) for i in range(2)
        ]
        split_records[label]=rec
        add_edge(ggraph,new_nodes[0],new_nodes[1])

    active=[]
    for old in rec["active_components"]:
        if old=="F0":
            active.append("F0")
        elif old in split_parts:
            active.extend(active_split_nodes(old))
        else:
            active.append(old)

    for i in range(len(active)):
        for j in range(i+1,len(active)):
            remove_edge(ggraph,active[i],active[j])
    for n in new_nodes:
        for a in active:
            add_edge(ggraph,n,a)

assert len(ggraph)==14
assert sum(len(v) for v in ggraph.values())//2==13
gfinite,gfar,gbranch,gchain,gspin=finite_d13_shape(ggraph,"F0")
assert gfar=="E12" and gbranch=="E11"

# Current identity -> F0 fixes the long chain; only spinor swap remains.
maps=[]
for swap in (False,True):
    m={identity:"F0"}
    for a,b in zip(tchain,gchain):
        m[a]=b
    gs=(gspin[1],gspin[0]) if swap else gspin
    m[tspin[0]]=gs[0]
    m[tspin[1]]=gs[1]
    assert len(m)==14
    # Verify every graph edge/nonedge.
    for a in tgraph:
        for b in tgraph:
            assert ((b in tgraph[a]) == (m[b] in ggraph[m[a]]))
    maps.append(m)

for mi,m in enumerate(maps,1):
    print(
        "Q24EFFCOMP_MAP|"
        f"map={mi}|"
        +";".join(f"{name}->{m[name]}" for name in ["A0"]+[f"R{i}" for i in range(1,14)])
        +"|status=PASS",
        flush=True,
    )

# ---------------------------------------------------------------------------
# Exact D12 special fibre and Theta.
# ---------------------------------------------------------------------------

keep=[i for i in range(13) if i!=2]  # pinned R3 is deleted.
C12=CartanPin.matrix_from_rows_and_columns(keep,keep)
assert (C12.rank(),ZZ(pari(C12).qfminim(2)[0]),abs(ZZ(C12.det())))==(12,264,4)
high12local=highest_root(C12)
m12=[ZZ(0)]*13
for i,v in zip(keep,high12local):
    m12[i]=ZZ(v)

Theta=D
for i in range(13):
    Theta -= m12[i]*effective[f"R{i+1}"]

assert Theta*ns*Theta==-2
assert Theta*ns*F8==2
assert Theta*ns*D==0

theta_pairs={name:int(Theta*ns*c) for name,c in effective.items()}
theta_hits=[name for name,v in theta_pairs.items() if v]
assert theta_hits==["R12"],theta_pairs
assert theta_pairs["R12"]==1

print(
    "Q24EFFCOMP_THETA|"
    f"deleted=R3->{maps[0]['R3']}|"
    f"theta_hit=R12->{maps[0]['R12']}|"
    f"pinned_affine=A0->{maps[0]['A0']}|"
    f"current_identity={identity}->F0|"
    f"square={Theta*ns*Theta}|old_degree={Theta*ns*F8}|new_degree={Theta*ns*D}|"
    "status=PASS_EXACT_D12_SPECIAL_FIBRE",
    flush=True,
)

# ---------------------------------------------------------------------------
# Current q24 divisor representative on ACTUAL effective components.
# D = O8 + P24 + V.
# Solve V = k F8 + sum(nonidentity r_C C), using the current identity.
# ---------------------------------------------------------------------------

V=D-O8-P24
assert V*ns*F8==0

nonidentity=[name for name in effective if name!=identity]
vbasis=matrix(QQ,[list(F8)]+[list(effective[name]) for name in nonidentity])
coeff=vbasis.solve_left(vector(QQ,V))
assert all(v in ZZ for v in coeff),coeff
coeff=vector(ZZ,coeff)
assert vector(QQ,coeff)*vbasis==vector(QQ,V)
k=ZZ(coeff[0])
r={name:ZZ(0) for name in effective}
for name,v in zip(nonidentity,coeff[1:]):
    r[name]=ZZ(v)

multiplicity={"A0":ZZ(1)}
for i in range(13):
    multiplicity[f"R{i+1}"]=ZZ(high13[i])

Drep_coeff={
    name:k*multiplicity[name]+r[name]
    for name in effective
}

# Special effective D12 fibre multiplicities on old fibre components.
special_mult={"A0":ZZ(0)}
for i in range(13):
    special_mult[f"R{i+1}"]=ZZ(m12[i])

generic_min={name:-Drep_coeff[name] for name in effective}
theta_order={
    name:special_mult[name]-Drep_coeff[name]
    for name in effective
}

print(
    "Q24EFFCOMP_VERTICAL|"
    f"F={k}|"
    f"nonidentity={','.join(name+':'+str(r[name]) for name in nonidentity)}|"
    f"support={sum(bool(r[name]) for name in nonidentity)}|"
    f"L1={sum(abs(int(r[name])) for name in nonidentity)}|status=PASS",
    flush=True,
)

for mi,m in enumerate(maps,1):
    for name in ["A0"]+[f"R{i}" for i in range(1,14)]:
        print(
            "Q24EFFCOMP_LOCAL|"
            f"map={mi}|pinned={name}|geo={m[name]}|"
            f"identity={int(name==identity)}|"
            f"fiber_mult={multiplicity[name]}|Dpair={Dpairs[name]}|"
            f"Drep={Drep_coeff[name]}|minord={generic_min[name]}|"
            f"special_mult={special_mult[name]}|thetaord={theta_order[name]}|"
            "status=PASS",
            flush=True,
        )

payload={
    "schema":"elkies-k3.h3-q24-effective-i9star-component-transport.v1",
    "status":"PASS_EXACT_Q24_EFFECTIVE_I9STAR_COMPONENT_TRANSPORT",
    "prime":int(args.prime),
    "current_identity_pinned_label":identity,
    "q24_pairings":Dpairs,
    "current_zero_pairings":Opairs,
    "transported_component_classes":{
        name:list(map(int,c)) for name,c in effective.items()
    },
    "transported_graph":{
        name:sorted(v) for name,v in tgraph.items()
    },
    "geometric_graph":{
        name:sorted(v) for name,v in ggraph.items()
    },
    "maps":maps,
    "D12_special_fibre":{
        "deleted_pinned_component":"R3",
        "D12_highest_multiplicities":{
            f"R{i+1}":int(m12[i]) for i in range(13)
        },
        "theta_class":list(map(int,Theta)),
        "theta_hit_pinned":"R12",
        "theta_hit_geometric_map1":maps[0]["R12"],
    },
    "current_divisor_representative":{
        "vertical_fibre_coefficient":int(k),
        "vertical_nonidentity_coefficients":{
            name:int(r[name]) for name in nonidentity
        },
        "component_coefficients":{
            name:int(Drep_coeff[name]) for name in effective
        },
        "generic_min_orders":{
            name:int(generic_min[name]) for name in effective
        },
        "theta_exact_orders":{
            name:int(theta_order[name]) for name in effective
        },
    },
    "next":(
        "Use the actual geometric map(s) and these component orders to pull the "
        "10-dimensional post-collision q24 space through the explicit blow-up "
        "charts. Generic min orders should cut it to H0 dimension 2; theta "
        "exact orders should isolate the special D12 fibre member."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q24-effective-i9star-components-mod-{args.prime}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24EFFCOMP_RESULT|"
    f"identity={identity}->F0|R3->{maps[0]['R3']}|R12->{maps[0]['R12']}|"
    f"A0->{maps[0]['A0']}|spinor_maps=2|status={payload['status']}",
    flush=True,
)
