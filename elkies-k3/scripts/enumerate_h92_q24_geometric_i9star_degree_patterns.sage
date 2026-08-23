#!/usr/bin/env sage -python
"""
Recover the q24 degree pattern on the ACTUAL resolved I9* fibre without
transporting effective components through Weyl reflections.

Facts used:
  * current canonical D13 old fibre is I9*;
  * explicit resolution gives its affine 14-component graph;
  * O8 and P24 both meet the identity component F0
      (P24 has D13 correction 0);
  * O8.P24 = 24;
  * D24.F8 = 2 and D24 is the D12/MW5 neighbour;
  * current exact decomposition has invariant fibre twist k = -7.

For every nonnegative component-degree pattern d_C = D24.C satisfying
    sum multiplicity(C) * d_C = 2,
enumerate those for which the zero-degree old components form D12.

For each candidate solve uniquely
    D24 = O + P24 + A + k F8,
    A = sum_{nonidentity C} a_C C,
directly in the affine-I9* intersection lattice.

The correct geometric pattern must have k=-7.  This identifies the actual
components hit by q24 and gives the divisor coefficients / resolved local
orders needed for the RR compiler, with no NS chamber guess.
"""

import argparse
import json
from collections import deque
from itertools import combinations_with_replacement
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, pari, sage_eval, vector


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


def add_edge(g,a,b):
    if a==b:
        return
    g.setdefault(a,set()).add(b)
    g.setdefault(b,set()).add(a)


def remove_edge(g,a,b):
    g.setdefault(a,set()).discard(b)
    g.setdefault(b,set()).discard(a)


def root_data(C):
    if C.nrows()==0:
        return (0,0,1)
    if C.rank()!=C.nrows():
        return (int(C.rank()),-1,0)
    qf=pari(C).qfminim(2)
    count=ZZ(qf[0])
    half=[vector(ZZ,c) for c in matrix(ZZ,qf[2]).columns()]
    roots=half+[-r for r in half]
    if not roots:
        return (0,0,1)
    rb=matrix(ZZ,[list(r) for r in roots]).row_module().basis_matrix()
    rg=rb*C*rb.transpose()
    return (int(rb.rank()),int(count),abs(int(rg.det())))


def highest_root(C):
    qf=pari(C).qfminim(2)
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
PROFILES=LOCAL/"q24-three-d12-current-equation-profiles.json"

for p in (RES,PROFILES):
    if not p.exists():
        raise SystemExit(f"Missing prerequisite: {p}")

res=json.loads(RES.read_text())
profiles=json.loads(PROFILES.read_text())
assert res["status"]=="PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS"
assert profiles["status"]=="PASS_EXACT_THREE_Q24_D12_CURRENT_EQUATION_PROFILES"
p85=next(r for r in profiles["profiles"] if int(r["orbit"])==85)
assert p85["height"]=="52"
assert p85["D13_correction"]=="0"
assert int(p85["P_dot_O"])==24
assert int(p85["vertical_fibre_coefficient"])==-7

# ---------------------------------------------------------------------------
# 1. Rebuild the resolved affine I9* graph.
# ---------------------------------------------------------------------------

Fp=GF(ZZ(args.prime))
S=PolynomialRing(Fp,names=("u","x","y"))
u,x,y=S.gens()

graph={"F0":set()}
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

for rec in sorted(res["centers"],key=lambda r:int(r["label"][1:])):
    label=rec["label"]
    parts=int(rec["new_exceptional_irreducible_components"])
    factors=rec["tangent_factors"]

    if parts==1:
        new_nodes=[label]
        graph.setdefault(label,set())
    else:
        assert parts==2 and len(factors)==2
        new_nodes=[label+"a",label+"b"]
        for node in new_nodes:
            graph.setdefault(node,set())
        split_parts[label]=[
            (new_nodes[i],factors[i]["factor"]) for i in range(2)
        ]
        split_records[label]=rec
        add_edge(graph,new_nodes[0],new_nodes[1])

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
            remove_edge(graph,active[i],active[j])
    for node in new_nodes:
        for old in active:
            add_edge(graph,node,old)

vertices=["F0"]+sorted(
    [v for v in graph if v!="F0"],
    key=lambda z:(int(z[1:3]),z)
)
assert len(vertices)==14
assert sum(len(v) for v in graph.values())//2==13

M=matrix(ZZ,14,14)
for i in range(14):
    M[i,i]=-2
for i,a in enumerate(vertices):
    for b in graph[a]:
        j=vertices.index(b)
        M[i,j]=1
assert M==M.transpose()
assert M.rank()==13

# Fibre multiplicities = primitive positive kernel of affine intersection matrix.
ker=M.right_kernel_matrix()
assert ker.nrows()==1
mult=vector(ZZ,ker.row(0))
if mult[0]<0:
    mult=-mult
g=ZZ(0)
for v in mult:
    from sage.all import gcd
    g=gcd(g,abs(ZZ(v)))
mult=vector(ZZ,[v//g for v in mult])
assert all(v>0 for v in mult)
assert mult[0]==1,(vertices,mult)
assert M*mult==0

print(
    "Q24GEODEG_GRAPH|"
    f"vertices={','.join(vertices)}|"
    f"multiplicities={','.join(map(str,mult))}|"
    f"mult1={','.join(vertices[i] for i,v in enumerate(mult) if v==1)}|"
    "status=PASS_AFFINE_I9STAR",
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Enumerate all nef old-component degree patterns of total fibre degree 2.
# ---------------------------------------------------------------------------

patterns=[]

# Enumerate allocations using recursion; degree entries are tiny.
def rec_alloc(i,remaining,current):
    if i==len(vertices):
        if remaining==0:
            patterns.append(tuple(current))
        return
    m=int(mult[i])
    for d in range(remaining//m+1):
        current.append(d)
        rec_alloc(i+1,remaining-m*d,current)
        current.pop()

rec_alloc(0,2,[])
assert patterns

nonid=list(range(1,14))
Mnon=M.matrix_from_rows_and_columns(nonid,nonid)
assert Mnon.det()!=0

B2=ZZ(44)  # (O+P)^2 = -2-2+2*(O.P), O.P=24.

records=[]
for dtuple in patterns:
    d=vector(ZZ,dtuple)
    zero=[i for i,v in enumerate(d) if v==0]
    Cz=(-M).matrix_from_rows_and_columns(zero,zero)
    rd=root_data(Cz)
    if rd!=(12,264,4):
        continue

    # A=sum a_i C_i. O+P meets only identity F0, twice.
    rhs_non=vector(QQ,[d[i] for i in nonid])
    a=vector(QQ,Mnon.solve_right(rhs_non))
    if not all(v in ZZ for v in a):
        continue
    a=vector(ZZ,a)

    # Affine-row consistency includes (O+P).F0 = 2.
    affine_contrib=sum(a[j-1]*M[0,j] for j in nonid)
    if affine_contrib != d[0]-2:
        continue

    A2=ZZ(a*Mnon*a)
    kval=QQ(-(B2+A2))/4
    if kval not in ZZ:
        continue
    kval=ZZ(kval)

    coeff=[]
    for i in range(14):
        ai=ZZ(0) if i==0 else a[i-1]
        coeff.append(kval*mult[i]+ai)

    positive=[vertices[i] for i,v in enumerate(d) if v]
    rec={
        "degrees":list(map(int,d)),
        "positive_components":positive,
        "zero_components":[vertices[i] for i in zero],
        "root_data":list(rd),
        "A_coefficients_nonidentity":list(map(int,a)),
        "A_square":int(A2),
        "fibre_twist_k":int(kval),
        "Drep_component_coefficients":list(map(int,coeff)),
        "generic_min_orders":[int(-v) for v in coeff],
    }
    records.append(rec)

records.sort(key=lambda r:(
    abs(r["fibre_twist_k"]+7),
    r["positive_components"],
))

print(
    "Q24GEODEG_ENUM|"
    f"all_degree_patterns={len(patterns)}|D12_integral_candidates={len(records)}|"
    "status=PASS",
    flush=True,
)

for idx,r in enumerate(records,1):
    print(
        "Q24GEODEG_CAND|"
        f"i={idx}|positive={','.join(r['positive_components'])}|"
        f"k={r['fibre_twist_k']}|A2={r['A_square']}|"
        f"A={','.join(map(str,r['A_coefficients_nonidentity']))}|"
        f"Drep={','.join(map(str,r['Drep_component_coefficients']))}|"
        "status=CANDIDATE",
        flush=True,
    )

selected=[r for r in records if r["fibre_twist_k"]==-7]
assert selected, "no geometric D12 degree pattern has the certified fibre twist k=-7"

print(
    "Q24GEODEG_KMINUS7|"
    f"count={len(selected)}|"
    f"patterns={';'.join(','.join(r['positive_components']) for r in selected)}|"
    "status=PASS_CERTIFIED_TWIST_FILTER",
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. For each surviving pattern, compute special D12-fibre multiplicities
#    and exact local orders of its affine bisection member.
# ---------------------------------------------------------------------------

for si,r in enumerate(selected,1):
    d=vector(ZZ,r["degrees"])
    zero=[i for i,v in enumerate(d) if v==0]
    C12=(-M).matrix_from_rows_and_columns(zero,zero)
    assert root_data(C12)==(12,264,4)
    high=highest_root(C12)

    special=[ZZ(0)]*14
    for i,v in zip(zero,high):
        special[i]=ZZ(v)

    coeff=vector(ZZ,r["Drep_component_coefficients"])
    theta_orders=[special[i]-coeff[i] for i in range(14)]

    r["D12_highest_multiplicities_full"]=list(map(int,special))
    r["theta_exact_orders"]=list(map(int,theta_orders))

    print(
        "Q24GEODEG_LOCALSET|"
        f"candidate={si}|positive={','.join(r['positive_components'])}|"
        f"zeroD12={','.join(r['zero_components'])}|status=BEGIN",
        flush=True,
    )
    for i,name in enumerate(vertices):
        print(
            "Q24GEODEG_LOCAL|"
            f"candidate={si}|component={name}|fiber_mult={mult[i]}|"
            f"Ddegree={d[i]}|Drep={coeff[i]}|"
            f"minord={-coeff[i]}|special_mult={special[i]}|"
            f"thetaord={theta_orders[i]}|status=PASS",
            flush=True,
        )
    print(
        "Q24GEODEG_LOCALSET|"
        f"candidate={si}|status=END",
        flush=True,
    )

payload={
    "schema":"elkies-k3.h3-q24-geometric-i9star-degree-pattern.v1",
    "status":"PASS_Q24_GEOMETRIC_I9STAR_DEGREE_PATTERN_ENUMERATION",
    "prime":int(args.prime),
    "vertices":vertices,
    "graph":{k:sorted(v) for k,v in graph.items()},
    "fibre_multiplicities":list(map(int,mult)),
    "O_and_P_meet":"F0",
    "O_dot_P":24,
    "certified_fibre_twist":-7,
    "all_D12_integral_candidates":records,
    "selected_k_minus_7":selected,
    "next":(
        "For each k=-7 survivor, pull the 10-dimensional post-collision basis "
        "through the explicit resolution and impose generic_min_orders. The "
        "correct geometric degree pattern must give resolved codimension 8 / "
        "kernel dimension 2. theta_exact_orders independently isolate the "
        "special D12 fibre generator."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q24-geometric-i9star-degree-pattern-mod-{args.prime}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24GEODEG_RESULT|"
    f"kminus7_candidates={len(selected)}|"
    f"patterns={';'.join(','.join(r['positive_components']) for r in selected)}|"
    f"status={payload['status']}",
    flush=True,
)
