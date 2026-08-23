#!/usr/bin/env sage -python
"""
Scout symbolically cheaper exits from the H3 D13/MW4 frame.

A. Compare all three direct q24 -> D12/MW5 neighbors.
B. Construct the certified lateral q4 neighbor
       D13/MW4 --q4--> A12+A1/MW4
   and export its root-adapted frame for a small-q growth search.
"""

import json
from pathlib import Path

from sage.all import (
    ZZ, QQ, block_diagonal_matrix, matrix, pari, vector
)

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
SCOUT = LOCAL / "route-scout"
SCOUT.mkdir(parents=True, exist_ok=True)

FRAME = ROOT / "elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
Q24_ART = GEN / "elkies-k3-h3-q6-q8-d13-q24-degree2.json"
ENGINE = ROOT / "elkies-k3/scripts/exact_neighbor_engine.sage"

exec(compile(ENGINE.read_text(), str(ENGINE), "exec"))

U2 = matrix(ZZ, ((0,1),(1,0)))
Q4_WITNESS = vector(ZZ, (
    3,-1,1,4,8,5,7,6,5,4,-1,-3,3,-1,1,-1,0,
))


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def write_gram(path, gram, header):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {line}" for line in header]
    lines += [" ".join(map(str,row)) for row in matrix(ZZ,gram).rows()]
    path.write_text("\n".join(lines)+"\n")


def connected_components(cartan):
    unseen=set(range(cartan.nrows()))
    out=[]
    while unseen:
        first=min(unseen)
        unseen.remove(first)
        todo=[first]
        comp=[]
        while todo:
            i=todo.pop()
            comp.append(i)
            adj=[j for j in list(unseen) if cartan[i,j] != 0]
            for j in adj:
                unseen.remove(j)
                todo.append(j)
        out.append(tuple(sorted(comp)))
    return tuple(out)


def component_name(cartan, comp):
    block=cartan.matrix_from_rows_and_columns(comp,comp)
    r=block.nrows()
    det=abs(ZZ(block.det()))
    count=ZZ(pari(block).qfminim(2)[0])
    if det==r+1 and count==r*(r+1):
        return f"A{r}"
    if r>=4 and det==4 and count==2*r*(r-1):
        return f"D{r}"
    exc={(6,3,72):"E6",(7,2,126):"E7",(8,1,240):"E8"}
    return exc.get((r,det,count),f"R{r}d{det}n{count}")


def ade_name(cartan):
    if cartan.nrows()==0:
        return "rootless"
    return "+".join(component_name(cartan,c) for c in connected_components(cartan))


def highest_root(cartan):
    qf=pari(cartan).qfminim(2)
    half=matrix(ZZ,qf[2]).transpose().rows()
    roots=[vector(ZZ,r) for r in half]+[-vector(ZZ,r) for r in half]
    positive=[r for r in roots if all(v>=0 for v in r)]
    assert positive
    return max(positive,key=lambda r:sum(r))


parent=load_gram(FRAME)
assert roots_and_data(parent)[2] == (13,312,4)
ns=block_diagonal_matrix(U2,-parent)
Fold=vector(ZZ,[1,0]+[0]*17)
Oold=vector(ZZ,[-1,1]+[0]*17)

root=parent[:13,:13]
assert ade_name(root)=="D13"

simple_components=tuple(
    vector(ZZ,[0,0]+[-ZZ(i==j) for j in range(17)])
    for i in range(13)
)
hr13=highest_root(root)
affine=Fold+vector(ZZ,[0,0]+list(hr13)+[0]*4)
assert affine*ns*affine==-2

if not Q24_ART.exists():
    raise SystemExit(f"Missing q24 search artifact: {Q24_ART}")

q24=json.loads(Q24_ART.read_text())
assert q24["status"]=="PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"

d12_records=[
    r for r in q24["neighbors"]
    if tuple(r["child_root_data"])==(12,264,4)
    and int(r["child_mw_rank"])==5
]
assert len(d12_records)==3, len(d12_records)

direct=[]
for rec in sorted(d12_records,key=lambda r:int(r["orbit_index"])):
    D=vector(ZZ,rec["fiber"])
    assert D*ns*D==0 and D*ns*Fold==2
    simple_pair=[int(D*ns*c) for c in simple_components]
    affine_pair=int(D*ns*affine)
    zero_idx=[i for i,v in enumerate(simple_pair) if v==0]
    inherited="rootless"
    if zero_idx:
        inherited=ade_name(root.matrix_from_rows_and_columns(zero_idx,zero_idx))

    H=matrix(QQ,rec["child_mw_height"])
    diag=[QQ(H[i,i]) for i in range(H.nrows())]
    shortest_diag=min(diag)

    dominant_labels=rec.get("dominant_labels", [])
    mw=rec.get("mw_projection", [])
    witness=rec["witness"]

    item={
        "orbit":int(rec["orbit_index"]),
        "mw_projection":list(map(int,mw)),
        "dominant_labels":list(map(int,dominant_labels)),
        "witness":list(map(int,witness)),
        "old_zero_intersection":int(D*ns*Oold),
        "simple_pairings":simple_pair,
        "affine_pairing":affine_pair,
        "inherited_ade":inherited,
        "degree_one_simple":[i+1 for i,v in enumerate(simple_pair) if v==1],
        "degree_one_affine":bool(affine_pair==1),
        "child_mw_height":[[str(v) for v in row] for row in H.rows()],
        "child_mw_basis_diagonal":[str(v) for v in diag],
        "shortest_basis_diagonal":str(shortest_diag),
        "witness_l1":int(sum(abs(ZZ(v)) for v in witness)),
        "witness_linf":int(max(abs(ZZ(v)) for v in witness)),
    }
    direct.append(item)

    print(
        "H3EXIT_Q24|"
        f"orbit={item['orbit']}|mw={','.join(map(str,item['mw_projection']))}|"
        f"labels={','.join(map(str,item['dominant_labels']))}|"
        f"O={item['old_zero_intersection']}|"
        f"simple={','.join(map(str,simple_pair))}|affine={affine_pair}|"
        f"inherited={inherited}|"
        f"degree1_simple={item['degree_one_simple']}|"
        f"shortest_MWdiag={shortest_diag}|"
        f"wL1={item['witness_l1']}|wLinf={item['witness_linf']}|"
        "status=PASS_D12_CANDIDATE",
        flush=True,
    )

# Lateral q4.
D4=vector(ZZ,[2,2]+list(Q4_WITNESS))
assert D4*ns*D4==0
assert D4*ns*Fold==2
assert D4*ns*Oold==0

q4_pair=[int(D4*ns*c) for c in simple_components]
assert q4_pair == [0,0,0,0,0,1,0,0,0,0,0,0,0]
assert int(D4*ns*affine)==1

split=primitive_hyperbolic_split(ns,D4)
raw_child=split["child_frame"]
minimized=minimize_child_frame(raw_child)
lat_frame=matrix(ZZ,minimized["frame"])
lat_data=tuple(map(int,minimized["root_data"]))
assert lat_data==(13,158,26),lat_data
lat_ade=ade_name(lat_frame[:13,:13])
assert set(lat_ade.split("+"))=={"A12","A1"}

lat_path=SCOUT/"d13-q4-lateral-a12a1-mw4-frame.txt"
write_gram(
    lat_path,
    lat_frame,
    [
        "H3 route scout lateral q4 child",
        "source D13/MW4",
        "q=4 factor order (2,2), old-fibre degree 2",
        "expected ADE A12+A1, MW4",
    ],
)

lat_H=minimized["mw_height"]
print(
    "H3EXIT_LATERAL_Q4|"
    "q=4|O=0|simple=0,0,0,0,0,1,0,0,0,0,0,0,0|affine=1|"
    f"child={lat_ade}|root_data={','.join(map(str,lat_data))}|MW=4|"
    f"height_diag={','.join(str(lat_H[i,i]) for i in range(lat_H.nrows()))}|"
    f"frame={lat_path}|status=PASS_CHEAP_LATERAL_EXIT",
    flush=True,
)

manifest={
    "schema":"elkies-k3.h3-d13-symbolic-exit-scout.v1",
    "status":"PASS_H3_D13_SYMBOLIC_EXIT_SCOUT",
    "direct_q24_d12":direct,
    "lateral_q4":{
        "q":4,
        "old_zero_intersection":0,
        "simple_pairings":q4_pair,
        "affine_pairing":1,
        "child_ade":lat_ade,
        "child_root_data":list(lat_data),
        "child_mw_rank":4,
        "child_mw_height":[[str(v) for v in row] for row in lat_H.rows()],
        "frame":str(lat_path.relative_to(ROOT)),
    },
}
manifest_path=SCOUT/"d13-symbolic-exit-scout.json"
manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{manifest_path}",flush=True)
print(
    "H3EXIT_SCOUT_RESULT|direct_D12=3|lateral_q4=PASS|"
    "next=SEARCH_LATERAL_CHILD_Q4_Q6_Q8_FOR_MW_GROWTH|status=PASS",
    flush=True,
)
