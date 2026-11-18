#!/usr/bin/env sage -python
"""GAP/PrimGrp classification-backed minimum-degree certificate for A10 local covers.

Run: sage -python research/arithmetic/verify_alternating_hasse_minimal_degree.py
The class-list completeness is imported from GAP's maximal-subgroup machinery
and PrimGrp classification; exact coset actions and weighted covers are replayed.
"""
from pathlib import Path
import argparse,hashlib,json,platform,sys
from itertools import combinations
from sage.all import libgap
from sage.version import version as sage_version

sys.path.insert(0,str(Path(__file__).resolve().parent))
from verify_alternating_subgroup_hasse import permutation

def require(ok,message):
    if not ok:
        raise AssertionError(message)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    G=libgap.AlternatingGroup(10)
    require(int(G.Size())==1814400,"A10 order")
    maxima=list(libgap.MaximalSubgroupClassReps(G))
    require(len(maxima)==7,"seven classified maximal subgroup classes")
    labels={10:"H1",45:"H2",120:"H3",210:"H4",126:"W5",945:"W2",2520:"M10"}
    witnesses={
      "cycle_9_1":[permutation(tuple(range(9)))],
      "cycle_7_3":[permutation(tuple(range(7)),tuple(range(7,10)))],
      "cycle_5_5":[permutation(tuple(range(5)),tuple(range(5,10)))],
      "DT_C2xC6":[permutation((0,1),(2,3)),permutation((0,1),(4,5,6,7,8,9))]
    }
    rows=[]
    for H in maxima:
        degree=int(G.Index(H))
        action=libgap.ActionHomomorphism(G,libgap.RightCosets(G,H),libgap.OnRight)
        hits={}
        for label,gens in witnesses.items():
            images=[]
            for g in gens:
                pg=libgap.PermList([x+1 for x in g])
                image=libgap.Image(action,pg)
                images.append([int(x)-1 for x in libgap.ListPerm(image,degree)])
            fixed=[j for j in range(degree) if all(g[j]==j for g in images)]
            hits[label]=len(fixed)
        rows.append(dict(label=labels[degree],index=degree,order=int(H.Size()),
            generators=[[int(x) for x in libgap.ListPerm(g,10)]
                        for g in libgap.GeneratorsOfGroup(H)],
            fixed_coset_counts=hits))
    rows.sort(key=lambda r:r["index"])
    require([r["index"] for r in rows]==[10,45,120,126,210,945,2520],
            "complete maximal-index list")
    all_witnesses=set(witnesses)
    candidates=[]
    for count in range(1,len(rows)+1):
        for inds in combinations(range(len(rows)),count):
            covered=set().union(*(set(k for k,v in rows[j]["fixed_coset_counts"].items() if v)
                                  for j in inds))
            if covered==all_witnesses:
                candidates.append(dict(rank=sum(rows[j]["index"] for j in inds),
                    classes=[rows[j]["label"] for j in inds],count=count))
    minimum=min(candidates,key=lambda x:x["rank"])
    triples=[x for x in candidates if x["count"]==3]
    best_triple=min(triples,key=lambda x:x["rank"])
    require(minimum==dict(rank=301,classes=["H1","H2","H3","W5"],count=4),
            "global rank lower bound")
    require(best_triple==dict(rank=1075,classes=["H1","H3","W2"],count=3),
            "three-class rank lower bound")
    require(min(x["count"] for x in candidates)==3,"no two-class cover of witnesses")
    primitive=[]
    for j in range(1,int(libgap.NrPrimitiveGroups(10))+1):
        P=libgap.PrimitiveGroup(10,j)
        primitive.append(dict(id=j,order=int(P.Size()),contained_in_A10=bool(G.IsSubgroup(P))))
    require(list(libgap.eval("ANMAXPRIMS[10]"))==[libgap(6)],
            "primitive-maximal table entry")
    source=Path(str(libgap.FilenameFunc(libgap.MaximalSubgroupsSymmAlt)))
    report=dict(schema_version=1,result="PASS",kind="classification_backed_finite_certificate",
        command="sage -python research/arithmetic/verify_alternating_hasse_minimal_degree.py",
        python=platform.python_version(),sage=str(sage_version),gap=str(libgap.eval("GAPInfo.Version")),
        primgrp=str(libgap.eval('PackageInfo("primgrp")[1].Version')),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        maximal_subgroup_method_source=dict(filename=source.name,
            sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            global_function="MaximalSubgroupsSymmAlt",start_line=int(libgap.StartlineFunc(libgap.MaximalSubgroupsSymmAlt)),
            primitive_maximal_ids_degree10=[6]),
        classified_maximal_subgroups=rows,primitive_groups_degree10=primitive,
        four_witness_minimum=minimum,three_class_minimum=best_triple,
        witness_scope="All cyclic subgroups and DT-inertia cyclic-quotient subgroups must in particular cover these four witnesses. The separate elementary proof supplies the matching upper bounds for the entire domain.",
        classification_assurance="Completeness imports GAP MaximalSubgroupClassReps and the classified PrimGrp library. This is an exact computational result using those inputs, not a classification-free proof, external replay, or formal verification.",
        independent_replay=False,formal_verification=False,external_review=False)
    output=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(output)
    print(output,end="")

if __name__=="__main__":
    main()
