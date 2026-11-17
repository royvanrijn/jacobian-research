#!/usr/bin/env python3
"""Finite Hall-capacity catalogue, independent of knot targets."""

from collections import Counter
import hashlib
import json
from pathlib import Path

from joint_a5 import group_model,generated,all_subgroups,tuple_orbits,canonical,certify_distinct_same_source

HERE=Path(__file__).resolve().parent


def main():
    elements,table,inverses,orders,automorphisms,core=group_model()
    subgroups,generators=all_subgroups(table)
    assert len(subgroups)==59
    records=[]
    for subgroup in subgroups:
        commutators={table[table[table[a][b]][inverses[a]]][inverses[b]] for a in subgroup for b in subgroup}
        derived=generated(commutators,table)
        normal=all(frozenset(table[table[g][h]][inverses[g]] for h in subgroup)==subgroup for g in range(60))
        if 1<len(subgroup)<60:
            assert len(derived)<len(subgroup)
        records.append({"elements":sorted(subgroup),"generators":generators[subgroup],
                        "derived_order":len(derived),"normal":normal})
    assert [len(r["elements"]) for r in records if r["normal"]]==[1,60]
    a4=[h for h in subgroups if len(h)==12]
    assert len(a4)==5
    action=[]
    for g in range(60):
        action.append(tuple(a4.index(frozenset(table[table[g][h]][inverses[g]] for h in subgroup)) for subgroup in a4))
    assert len(set(action))==60
    # The faithful action on the five order12 subgroups bounds Aut(A5)
    # by5!=120; the120 explicit S5 conjugations attain that bound.
    eulerian={}
    for d in [1,2,3,4]:
        counts={}
        for h in subgroups:
            counts[h]=len(h)**d-sum(counts[k] for k in counts if k<h)
        eulerian[d]=counts[subgroups[-1]]
    pairs=tuple_orbits(2,table,orders,automorphisms)
    triples=tuple_orbits(3,table,orders,automorphisms)
    assert len(pairs)==19 and len(triples)==1668
    assert eulerian[2]==2280==120*len(pairs)
    assert eulerian[3]==120*len(triples)
    assert len({canonical(p,automorphisms) for p in pairs})==19
    assert len({canonical(p,automorphisms) for p in triples})==1668
    controls={}
    controls["free_group_F2_kernel_capacity"]=len(pairs)
    controls["free_group_F3_kernel_capacity"]=len(triples)
    controls["two_trefoil_sum_A5_kernels"]=sum(all(orders[x] in [1,3] for x in p) for p in pairs)
    controls["three_trefoil_sum_A5_kernels"]=sum(all(orders[x] in [1,3] for x in p) for p in triples)
    assert controls["two_trefoil_sum_A5_kernels"]==1
    assert controls["three_trefoil_sum_A5_kernels"]==49
    twenty=[{"source_group":"F3-control","generator_images":list(p)+( [0] )} for p in pairs]
    twenty.append({"source_group":"F3-control","generator_images":[0,*pairs[0]]})
    assert certify_distinct_same_source(twenty,"F3-control",table,automorphisms)==20
    one=triples[0]
    duplicates=[{"source_group":"same-H","generator_images":list(one)},
                {"source_group":"same-H","generator_images":[automorphisms[1][x] for x in one]}]
    try:
        certify_distinct_same_source(duplicates,"same-H",table,automorphisms)
    except AssertionError as error:
        assert "same kernel" in str(error)
    else:
        raise AssertionError("duplicate kernels accepted")
    mixed=[{"source_group":"knot-A","generator_images":list(triples[0])},
           {"source_group":"knot-B","generator_images":list(triples[1])}]
    try:
        certify_distinct_same_source(mixed,"knot-A",table,automorphisms)
    except AssertionError as error:
        assert "different source" in str(error)
    else:
        raise AssertionError("cross-knot pooling accepted")
    output={"status":"PASS","group":"A5 even permutations of five letters","elements":elements,
            "subgroups":records,"faithful_action_on_five_order12_subgroups":action,
            "automorphism_order":120,"generating_tuple_counts":eulerian,
            "generating_pair_orbits":pairs,"generating_triple_orbits":triples,
            "kernel_capacities":{d:eulerian[d]//120 for d in eulerian},
            "controls":controls,"twenty_distinct_free_group_maps":twenty,
            "duplicate_kernel_rejected":True,"cross_source_pooling_rejected":True,
            "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ["joint_a5.py","build_joint_a5_catalogue.py"]},
            "scope":"Classical Hall capacity; no knot maps combined across sources. All subgroup/automorphism facts are finite exact checks, not a formal proof of the general kernel-count theorem."}
    (HERE/"joint_a5_catalogue.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"status":"PASS","subgroups":len(subgroups),"Aut_order":120,
                      "generating_tuple_counts":eulerian,"kernel_capacities":output["kernel_capacities"],"controls":controls},indent=2))


if __name__=="__main__":
    main()
