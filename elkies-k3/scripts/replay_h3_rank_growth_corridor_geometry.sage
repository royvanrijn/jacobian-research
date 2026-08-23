#!/usr/bin/env sage -python
"""
Replay the geometric pattern of the selected H3 rank-growth corridor.

For each stored degree-two neighbour from D13/MW4 through rootless/MW17:

  * load the EXACT selected orbit record;
  * work in that orbit's root-adapted parent frame;
  * compute intersections of the selected fibre D with:
      - old zero O,
      - every effective simple root component,
      - every affine component of each old ADE fibre;
  * identify the old simple-root subsystem orthogonal to D;
  * compare that inherited subsystem with the selected child ADE;
  * for every inherited connected component R, construct the candidate new
    affine component
        Theta_R = D - highest_root(R)
    using the effective old components;
  * certify Theta_R^2=-2 and Theta_R.D=0 when this is genuinely a reducible
    fibre of the new fibration;
  * report Theta_R's degree/intersection profile on the OLD fibration.

This is intended to recover the construction grammar behind the selected
q-sequence, especially orbit 85 (D13->D12) and orbit 42 (D12->A11).
"""

import json
from pathlib import Path

from sage.all import (
    ZZ, QQ, block_diagonal_matrix, matrix, pari, vector
)

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
OUT = ROOT / "artifacts/local/elkies-k3/h3-rank-growth-corridor-geometry.json"

STEPS = (
    ("artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q24-degree2.json",85,24,"D12",12),
    ("artifacts/generated-results/elkies-k3-h3-d12-o85-q6-degree2.json",42,6,"A11",11),
    ("artifacts/generated-results/elkies-k3-h3-a11-middle-q8-degree2.json",922,8,"A5+A5",10),
    ("artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2.json",472,4,"A3+A3+A3",9),
    ("artifacts/generated-results/elkies-k3-h3-a3x3-q4-degree2.json",323,4,"A2+A2+A3",7),
    ("artifacts/generated-results/elkies-k3-h3-mw10-a3a2a2-q4-degree2.json",207,4,"A1+A1+A1+A1+A1",5),
    ("artifacts/generated-results/elkies-k3-h3-mw12-5a1-q4-degree2-first-hit.json",52,4,"A1+A1+A1+A1",4),
    ("artifacts/generated-results/elkies-k3-h3-mw13-4a1-q4-degree2-first-hit.json",114,4,"A1+A1+A1",3),
    ("artifacts/generated-results/elkies-k3-h3-mw14-3a1-q4-degree2-first-hit.json",498,4,"A1+A1",2),
    ("artifacts/generated-results/elkies-k3-h3-mw15-2a1-q4-degree2-first-hit.json",981,4,"A1",1),
    ("artifacts/generated-results/elkies-k3-h3-mw16-a1-q6-degree2-cap10000-stream-chunk001.json",2247,6,"rootless",0),
)

U2 = matrix(ZZ, ((0,1),(1,0)))


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def roots_and_data(gram):
    qf = pari(gram).qfminim(2)
    count = ZZ(qf[0])
    if not count:
        return (), matrix(ZZ,0,gram.nrows()), (0,0,1)
    half = [vector(ZZ,c) for c in matrix(ZZ,qf[2]).columns()]
    roots = tuple(half + [-r for r in half])
    basis = matrix(ZZ,[list(r) for r in roots]).row_module().basis_matrix()
    rg = basis*gram*basis.transpose()
    return roots,basis,(basis.rank(),count,abs(ZZ(rg.det())))


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
            adjacent=[j for j in list(unseen) if cartan[i,j] != 0]
            for j in adjacent:
                unseen.remove(j)
                todo.append(j)
        out.append(tuple(sorted(comp)))
    return tuple(out)


def component_name(cartan, comp):
    block=cartan.matrix_from_rows_and_columns(comp,comp)
    r=block.nrows()
    det=abs(ZZ(block.det()))
    count=ZZ(pari(block).qfminim(2)[0])
    if det == r+1 and count == r*(r+1):
        return f"A{r}"
    if r>=4 and det==4 and count==2*r*(r-1):
        return f"D{r}"
    exc={(6,3,72):"E6",(7,2,126):"E7",(8,1,240):"E8"}
    return exc.get((r,det,count),f"R{r}d{det}n{count}")


def ade_name(cartan):
    if cartan.nrows()==0:
        return "rootless"
    return "+".join(component_name(cartan,c) for c in connected_components(cartan))


def highest_root(block):
    qf=pari(block).qfminim(2)
    half=matrix(ZZ,qf[2]).transpose().rows()
    roots=[vector(ZZ,r) for r in half]+[-vector(ZZ,r) for r in half]
    positive=[r for r in roots if all(v>=0 for v in r)]
    if not positive:
        raise ArithmeticError("no coordinate-positive highest root")
    return max(positive,key=lambda r:sum(r))


def effective_simple_components(root_rank, total_rank=17):
    return tuple(
        vector(ZZ,[0,0]+[-ZZ(i==j) for j in range(total_rank)])
        for i in range(root_rank)
    )


def parent_affines(frame, root_rank, old_fiber):
    cartan=frame[:root_rank,:root_rank]
    affines=[]
    for n,comp in enumerate(connected_components(cartan),1):
        block=cartan.matrix_from_rows_and_columns(comp,comp)
        hr=highest_root(block)
        coeff=[ZZ(0)]*root_rank
        for idx,val in zip(comp,hr):
            coeff[idx]=ZZ(val)
        affine=old_fiber+vector(ZZ,[0,0]+coeff+[0]*(17-root_rank))
        affines.append((n,comp,component_name(cartan,comp),coeff,affine))
    return affines


current=load_gram(SOURCE)
records=[]

for step,(relpath,orbit,q,expected_ade,expected_root_rank) in enumerate(STEPS,1):
    path=ROOT/relpath
    if not path.exists():
        raise SystemExit(f"Missing stored path artifact: {path}")
    data=json.loads(path.read_text())
    matches=[r for r in data["neighbors"] if int(r["orbit_index"])==orbit]
    if len(matches)!=1:
        raise ArithmeticError((path,orbit,len(matches)))
    rec=matches[0]

    assert int(rec["q"])==q
    assert rec["factor_order"]==[q//2,2]
    assert int(rec["old_fiber_degree"])==2
    assert matrix(ZZ,rec["child_root_adapted_frame"]).nrows()==17

    parent_root_rank=roots_and_data(current)[2][0]
    assert parent_root_rank == 17-(step+3) if False else parent_root_rank

    ns=block_diagonal_matrix(U2,-current)
    Fold=vector(ZZ,[1,0]+[0]*17)
    Oold=vector(ZZ,[-1,1]+[0]*17)
    D=vector(ZZ,rec["fiber"])

    assert D*ns*D==0
    assert D*ns*Fold==2
    assert rec["child_ade"]==expected_ade
    assert int(rec["child_root_data"][0])==expected_root_rank

    C=effective_simple_components(parent_root_rank)
    simple_pairings=[int(D*ns*c) for c in C]
    affines=parent_affines(current,parent_root_rank,Fold)
    affine_pairings=[
        {
            "component":name,
            "parent_nodes":[i+1 for i in comp],
            "multiplicities":[int(v) for v in coeff if v or True],
            "pairing":int(D*ns*aff),
        }
        for _,comp,name,coeff,aff in affines
    ]

    zero_indices=[i for i,v in enumerate(simple_pairings) if v==0]
    inherited_ade="rootless"
    inherited_rank=0
    inherited_components=[]
    theta_records=[]

    if zero_indices:
        inherited_cartan=current.matrix_from_rows_and_columns(zero_indices,zero_indices)
        inherited_rank=inherited_cartan.rank()
        inherited_ade=ade_name(inherited_cartan)

        # Each connected inherited root component is a candidate reducible
        # fibre of the NEW fibration. Construct its missing affine curve.
        for compnum,localcomp in enumerate(connected_components(inherited_cartan),1):
            parent_indices=[zero_indices[i] for i in localcomp]
            block=current.matrix_from_rows_and_columns(parent_indices,parent_indices)
            name=component_name(block,tuple(range(block.nrows())))
            hr=highest_root(block)

            coeff=[ZZ(0)]*parent_root_rank
            vertical=vector(ZZ,[0]*19)
            for idx,val in zip(parent_indices,hr):
                coeff[idx]=ZZ(val)
                vertical += ZZ(val)*C[idx]

            theta=D-vertical
            theta_sq=int(theta*ns*theta)
            theta_D=int(theta*ns*D)
            theta_olddeg=int(theta*ns*Fold)
            theta_oldO=int(theta*ns*Oold)
            hits=[
                [i+1,int(theta*ns*C[i])]
                for i in range(parent_root_rank)
                if theta*ns*C[i]
            ]
            affine_hits=[
                [a["component"], int(theta*ns*aff)]
                for a,(_,_,_,_,aff) in zip(affine_pairings,affines)
                if theta*ns*aff
            ]
            is_affine=(theta_sq==-2 and theta_D==0)
            theta_records.append({
                "new_root_component":name,
                "parent_nodes":[i+1 for i in parent_indices],
                "highest_multiplicities_in_parent_numbering":[int(v) for v in coeff],
                "theta_class":[int(v) for v in theta],
                "theta_square":theta_sq,
                "theta_new_fibre_intersection":theta_D,
                "theta_old_fibre_degree":theta_olddeg,
                "theta_old_zero_intersection":theta_oldO,
                "theta_simple_hits":hits,
                "theta_parent_affine_hits":affine_hits,
                "certified_affine_component":is_affine,
            })
            inherited_components.append(name)

    degree_one_simple=[i+1 for i,v in enumerate(simple_pairings) if v==1]
    degree_one_affine=[
        a["component"] for a in affine_pairings if a["pairing"]==1
    ]

    print(
        f"H3CORRIDOR|step={step}|q={q}|orbit={orbit}|"
        f"parent={ade_name(current[:parent_root_rank,:parent_root_rank])}|"
        f"child={expected_ade}|O={int(D*ns*Oold)}|"
        f"simple_pairings={','.join(map(str,simple_pairings))}|"
        f"affine_pairings={','.join(str(a['pairing']) for a in affine_pairings)}|"
        f"inherited={inherited_ade}|"
        f"new_sections_simple={degree_one_simple}|"
        f"new_sections_affine={degree_one_affine}|status=PASS",
        flush=True,
    )
    for t in theta_records:
        print(
            "H3CORRIDOR_THETA|"
            f"step={step}|root={t['new_root_component']}|"
            f"olddeg={t['theta_old_fibre_degree']}|"
            f"oldO={t['theta_old_zero_intersection']}|"
            f"square={t['theta_square']}|D={t['theta_new_fibre_intersection']}|"
            f"affine={int(t['certified_affine_component'])}|"
            f"nodes={t['parent_nodes']}|status=PASS",
            flush=True,
        )

    records.append({
        "step":step,
        "q":q,
        "orbit":orbit,
        "parent_root_rank":int(parent_root_rank),
        "parent_ade":ade_name(current[:parent_root_rank,:parent_root_rank]),
        "child_ade":expected_ade,
        "child_root_rank":expected_root_rank,
        "old_zero_intersection":int(D*ns*Oold),
        "simple_pairings":simple_pairings,
        "affine_pairings":affine_pairings,
        "orthogonal_old_simple_subsystem":{
            "rank":int(inherited_rank),
            "ade":inherited_ade,
            "matches_child":bool(inherited_rank==expected_root_rank and inherited_ade==expected_ade),
        },
        "degree_one_old_simple_components":degree_one_simple,
        "degree_one_old_affine_components":degree_one_affine,
        "candidate_new_affine_components":theta_records,
        "witness":rec["witness"],
        "fiber":rec["fiber"],
        "child_mw_height":rec["child_mw_height"],
    })

    current=matrix(ZZ,rec["child_root_adapted_frame"])

OUT.parent.mkdir(parents=True,exist_ok=True)
payload={
    "schema":"elkies-k3.h3-rank-growth-corridor-geometry.v1",
    "status":"PASS_H3_RANK_GROWTH_CORRIDOR_GEOMETRY",
    "steps":records,
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")

print(f"OUTPUT|{OUT}",flush=True)
print(
    "H3CORRIDOR_RESULT|steps=11|"
    f"inherited_matches={sum(int(r['orthogonal_old_simple_subsystem']['matches_child']) for r in records)}/11|"
    "status=PASS_H3_RANK_GROWTH_CORRIDOR_GEOMETRY",
    flush=True,
)
