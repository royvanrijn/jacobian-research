#!/usr/bin/env python3
"""Independent PSL2/PGL2(F5) replay of per-knot quotient-kernel counts.

Tuple orbits are formed by exhaustive orbit marking, not the producer's
earliest-entry/centralizer normalization. Diagram propagation is by
crossing closure, not the producer's fixed trace. PD parsing is shared.
"""

from collections import Counter
import hashlib
from itertools import permutations,product
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def matrix_model():
    def multiply(a,b):
        return ((a[0]*b[0]+a[1]*b[2])%5,(a[0]*b[1]+a[1]*b[3])%5,
                (a[2]*b[0]+a[3]*b[2])%5,(a[2]*b[1]+a[3]*b[3])%5)
    def determinant(a):
        return (a[0]*a[3]-a[1]*a[2])%5
    def inverse(a):
        scale=pow(determinant(a),-1,5)
        return tuple(scale*x%5 for x in (a[3],-a[1],-a[2],a[0]))
    def special(a):
        return min(a,tuple(-x%5 for x in a))
    def projective(a):
        scale=pow(next(x for x in a if x),-1,5)
        return tuple(scale*x%5 for x in a)
    identity=special((1,0,0,1))
    elements=sorted({special(a) for a in product(range(5),repeat=4) if determinant(a)==1})
    elements.remove(identity)
    elements.insert(0,identity)
    index={a:i for i,a in enumerate(elements)}
    table=[[index[special(multiply(a,b))] for b in elements] for a in elements]
    inverses=[index[special(inverse(a))] for a in elements]
    larger=sorted({projective(a) for a in product(range(5),repeat=4) if determinant(a)})
    automorphisms=[tuple(index[special(multiply(multiply(a,b),inverse(a)))] for b in elements) for a in larger]
    assert len(elements)==60 and len(larger)==len(set(automorphisms))==120
    core=[[table[table[b][inverses[a]]][b] for a in range(60)] for b in range(60)]
    orders=[]
    for i in range(60):
        current,n=i,1
        while current:
            current=table[current][i]
            n+=1
        orders.append(n)
    assert Counter(orders)=={1:1,2:15,3:20,5:24}
    return table,orders,automorphisms,core


def subgroup_elements(values,table):
    elements,queue={0},[0]
    while queue:
        a=queue.pop()
        for b in values:
            c=table[a][b]
            if c not in elements:
                elements.add(c)
                queue.append(c)
    return frozenset(elements)


def subgroup_size(values,table):
    return len(subgroup_elements(values,table))


def marked_orbits(arity,table,automorphisms):
    visited=set()
    representatives=[]
    generating=0
    for values in product(range(60),repeat=arity):
        if values in visited:
            continue
        orbit={tuple(a[i] for i in values) for a in automorphisms}
        assert not (orbit&visited)
        visited.update(orbit)
        if subgroup_size(values,table)==60:
            assert len(orbit)==120
            representatives.append(min(orbit))
            generating+=len(orbit)
    assert len(visited)==60**arity
    return sorted(representatives),generating


def crossing_closure(crossings,seeds,values,core):
    labels={s:a for s,a in zip(seeds,(0,)+tuple(values))}
    changed=True
    while changed:
        changed=False
        for over,a,b in crossings:
            if over not in labels:
                continue
            if a not in labels and b in labels:
                a,b=b,a
            if a not in labels:
                continue
            value=core[labels[over]][labels[a]]
            if b in labels:
                if labels[b]!=value:
                    return False
            else:
                labels[b]=value
                changed=True
    assert len(labels)==len(crossings)
    return True


def permutation_controls(catalogue):
    elements=[tuple(p) for p in catalogue["elements"]]
    assert elements==[p for p in permutations(range(5)) if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0]
    index={p:i for i,p in enumerate(elements)}
    def multiply(a,b):
        return tuple(a[b[i]] for i in range(5))
    def inverse(a):
        return tuple(a.index(i) for i in range(5))
    automorphisms=[]
    for a in permutations(range(5)):
        automorphisms.append(tuple(index[multiply(multiply(a,b),inverse(a))] for b in elements))
    table=[[index[multiply(a,b)] for b in elements] for a in elements]
    subgroups=[frozenset(r["elements"]) for r in catalogue["subgroups"]]
    assert len(set(subgroups))==59 and frozenset(range(60)) in subgroups and frozenset([0]) in subgroups
    for h,record in zip(subgroups,catalogue["subgroups"]):
        assert all(table[a][b] in h for a in h for b in h)
        assert set(record["generators"])<=h
        assert subgroup_elements(record["generators"],table)==h
        for value in range(60):
            # Any subgroup chain can be built by adjoining one element.
            # Closure within the supplied full lattice certifies completeness.
            if value in h:
                continue
            generated=set(h)
            generated.add(value)
            changed=True
            while changed:
                before=len(generated)
                generated.update(table[a][b] for a in list(generated) for b in [value,*record["generators"]])
                changed=len(generated)>before
            assert frozenset(generated) in subgroups
    order12=[h for h in subgroups if len(h)==12]
    assert len(order12)==5
    action=[]
    for g in range(60):
        p=elements[g]
        pinv=inverse(p)
        action.append(tuple(order12.index(frozenset(index[multiply(multiply(p,elements[h]),pinv)] for h in subgroup)) for subgroup in order12))
    assert len(set(action))==60
    assert [list(a) for a in action]==catalogue["faithful_action_on_five_order12_subgroups"]
    twenty=catalogue["twenty_distinct_free_group_maps"]
    assert len(twenty)==20 and {m["source_group"] for m in twenty}=={"F3-control"}
    canonical=lambda values:min(tuple(a[i] for i in values) for a in automorphisms)
    assert len({canonical(m["generator_images"]) for m in twenty})==20
    assert all(subgroup_size(m["generator_images"],table)==60 for m in twenty)
    first=twenty[0]["generator_images"]
    assert len({canonical([a[i] for i in first]) for a in automorphisms})==1
    return elements,table,automorphisms,canonical


def additional_controls(catalogue,table,orders,core,triples,permutation_table,canonical):
    orbifolds=json.loads((HERE/"joint_a5_orbifold_results.json").read_text())
    plan=json.loads((HERE/"joint_a5_orbifold_plan.json").read_text())
    pretzel=json.loads((HERE/"joint_a5_pretzel_control.json").read_text())
    geometry=json.loads((HERE/"joint_a5_geometry_control.json").read_text())
    for receipt in [orbifolds,plan,pretzel,geometry]:
        for path,digest in receipt["input_sha256"].items():
            assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    inverses=[next(j for j in range(60) if table[i][j]==0) for i in range(60)]
    permutation_inverses=[next(j for j in range(60) if permutation_table[i][j]==0) for i in range(60)]
    permutation_orders=[]
    for i in range(60):
        current,n=i,1
        while current:
            current=permutation_table[current][i]
            n+=1
        permutation_orders.append(n)
    assert orbifolds["status"]=="COMPLETE" and len(orbifolds["results"])==len(plan["patterns"])==40
    summaries=[]
    for pattern,row in zip(plan["patterns"],orbifolds["results"]):
        assert pattern==row["cone_orders"]
        a,b,c,d=pattern
        count=sum(a%orders[x]==b%orders[y]==c%orders[z]==0
                  and d%orders[inverses[table[table[x][y]][z]]]==0 for x,y,z in triples)
        keys=set()
        for mapping in row["representatives"]:
            assert mapping["source_group"]==row["source_group"]
            values=mapping["generator_images"]
            assert values==catalogue["generating_triple_orbits"][mapping["orbit_index"]]
            assert subgroup_size(values,permutation_table)==60
            x,y,z=values
            fourth=permutation_inverses[permutation_table[permutation_table[x][y]][z]]
            assert mapping["fourth_image"]==fourth
            assert all(exponent%permutation_orders[value]==0 for exponent,value in zip(pattern,[x,y,z,fourth]))
            key=canonical(values)
            assert key not in keys
            keys.add(key)
        assert count==len(keys)==row["A5_quotient_kernels"]
        assert row["rank_at_least3"]==(count>=20)
        summaries.append({"cone_orders":pattern,"A5_quotient_kernels":count})
    assert sum(r["A5_quotient_kernels"]>=20 for r in summaries)==39
    assert [r for r in summaries if r["cone_orders"]==[2,3,3,5]][0]["A5_quotient_kernels"]==23
    from verify_diagram_audit import map_witness,plane_map,walk
    pd=pretzel["pd"]
    plane_map(pd)
    assert len(pd)==13 and len(walk(pd))==26
    witness=map_witness(geometry["pd"],pd)
    raw_witness=map_witness(pretzel["raw_tangle_pd"],pd)
    crossings=[c[:3] for c in oriented_diagram(pd)]
    seeds=pretzel["seeds"]
    assert len(set(seeds))==4 and pretzel["normalized_arc"]==seeds[0]
    seen=set(seeds)
    for over,a,b in pretzel["seed_trace"]:
        assert over in seen and a in seen and b not in seen
        assert (over,a,b) in crossings or (over,b,a) in crossings
        seen.add(b)
    assert len(seen)==13
    count=sum(crossing_closure(crossings,seeds,values,core) for values in triples)
    assert count==pretzel["A5_quotient_kernels"]==23
    keys=set()
    for mapping in pretzel["epimorphism_representatives"]:
        assert mapping["source_group"]==pretzel["source_group"]
        values=mapping["generator_images"]
        assert values==catalogue["generating_triple_orbits"][mapping["orbit_index"]]
        assert subgroup_size(values,permutation_table)==60
        key=canonical(values)
        assert key not in keys
        keys.add(key)
        labels=mapping["arc_images"]
        assert labels[seeds[0]]==0 and [labels[s] for s in seeds[1:]]==values
        for over,a,b in crossings:
            assert permutation_table[permutation_table[labels[over]][permutation_inverses[labels[a]]]][labels[over]]==labels[b]
    assert len(keys)==23
    raw_count=sum(crossing_closure(crossings,seeds,values,core) for values in product(range(60),repeat=3))
    image_orders=Counter()
    allowed=[[i for i in range(60) if exponent%orders[i]==0] for exponent in [2,3,3,5]]
    for x,y,z in product(*allowed[:3]):
        fourth=inverses[table[table[x][y]][z]]
        if fourth in allowed[3]:
            image_orders[subgroup_size((x,y,z),table)]+=1
    assert raw_count==sum(image_orders.values())==pretzel["Hom_to_A5"]==2901
    assert image_orders=={1:1,3:20,12:120,60:2760}
    assert {int(k):v for k,v in pretzel["Hom_image_orders"].items()}==image_orders
    assert pretzel["raw_Hom_meridional_lower_bound"]==3
    assert pretzel["joint_kernel_meridional_lower_bound"]==pretzel["bridge_upper_bound"]==4
    return {"orbifold_patterns":summaries,"orbifold_patterns_reaching_twenty":39,
            "known_pretzel_control":{"knot":"P(2,3,3,5)","A5_quotient_kernels":23,"Hom_to_A5":2901,
                                    "Hom_image_orders":dict(image_orders),"raw_lower":3,"joint_lower":4,"bridge_upper":4,
                                    "spherogram_to_paired_endpoint_diagram":witness,"raw_to_oriented_diagram":raw_witness}}


def main():
    report=json.loads((HERE/"joint_a5_results.json").read_text())
    panel=json.loads((HERE/"joint_a5_panel.json").read_text())
    catalogue=json.loads((HERE/"joint_a5_catalogue.json").read_text())
    for receipt in [report,panel,catalogue]:
        for path,digest in receipt["input_sha256"].items():
            assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    elements,permutation_table,permutation_aut,canonical=permutation_controls(catalogue)
    table,orders,automorphisms,core=matrix_model()
    pairs,pair_count=marked_orbits(2,table,automorphisms)
    triples,triple_count=marked_orbits(3,table,automorphisms)
    assert pair_count==2280 and len(pairs)==19
    assert triple_count==200160 and len(triples)==1668
    assert len(catalogue["generating_pair_orbits"])==19 and len(catalogue["generating_triple_orbits"])==1668
    assert catalogue["generating_tuple_counts"]["2"]==pair_count
    assert catalogue["generating_tuple_counts"]["3"]==triple_count
    assert sum(all(orders[x] in [1,3] for x in p) for p in pairs)==1
    assert sum(all(orders[x] in [1,3] for x in p) for p in triples)==49
    assert report["status"]=="COMPLETE" and len(panel["cases"])==len(report["results"])==293
    results=[]
    for case,result in zip(panel["cases"],report["results"]):
        assert case["name"]==result["name"] and result["status"]=="COMPLETE"
        crossings=[c[:3] for c in oriented_diagram(case["pd"])]
        seeds=result["seeds"]
        seen=set(seeds)
        assert len(seen)==4 and result["normalized_arc"]==seeds[0]
        for over,a,b in result["seed_trace"]:
            assert over in seen and a in seen and b not in seen
            assert (over,a,b) in crossings or (over,b,a) in crossings
            seen.add(b)
        assert len(seen)==13
        source_hash=hashlib.sha256(json.dumps({"name":case["name"],"pd":case["pd"],"normalized_arc":seeds[0]},sort_keys=True).encode()).hexdigest()
        source_group=case["name"]+":double-cover:"+source_hash
        assert source_group==result["source_group"]
        keys=set()
        for mapping in result["epimorphism_representatives"]:
            assert mapping["source_group"]==source_group
            values=mapping["generator_images"]
            assert values==catalogue["generating_triple_orbits"][mapping["orbit_index"]]
            assert subgroup_size(values,permutation_table)==60
            key=canonical(values)
            assert key not in keys
            keys.add(key)
            labels=mapping["arc_images"]
            assert labels[seeds[0]]==0 and [labels[s] for s in seeds[1:]]==values
            for over,a,b in crossings:
                p,q,r=elements[labels[over]],elements[labels[a]],elements[labels[b]]
                qinv=tuple(q.index(i) for i in range(5))
                assert tuple(p[qinv[p[i]]] for i in range(5))==r
        count=sum(crossing_closure(crossings,seeds,values,core) for values in triples)
        assert count==len(keys)==result["A5_quotient_kernels"]
        assert result["orbit_trials"]==1668 and result["sharp_four"]==(count>=20)
        results.append({"name":case["name"],"source_group":source_group,"A5_quotient_kernels":count,
                        "sharp_four":count>=20,"four_seed_upper_verified":True})
    histogram=dict(sorted(Counter(r["A5_quotient_kernels"] for r in results).items()))
    assert histogram=={0:151,1:118,2:22,3:1,5:1}
    assert max(r["A5_quotient_kernels"] for r in results)==5
    assert not any(r["sharp_four"] for r in results)
    assert [r["A5_quotient_kernels"] for r in results[:16]]==[1,2,1,2,1,2,1,2,0,1,0,1,0,1,0,1]
    controls=additional_controls(catalogue,table,orders,core,triples,permutation_table,canonical)
    inputs=["verify_joint_a5.py","joint_a5.py","build_joint_a5_catalogue.py","joint_a5_catalogue.json",
            "freeze_joint_a5_panel.py","joint_a5_panel.json","search_joint_a5_kernels.py","joint_a5_results.json","parabolic_lifts.py",
            "audit_joint_a5_orbifolds.py","joint_a5_orbifold_plan.json","joint_a5_orbifold_results.json",
            "build_joint_a5_pretzel_control.py","joint_a5_pretzel_control.json",
            "build_joint_a5_geometry_control.py","joint_a5_geometry_control.json","verify_diagram_audit.py"]
    receipt={"status":"PASS","cases":293,"generating_pairs":2280,"automorphisms":120,
             "pair_kernel_capacity":19,"triple_kernel_capacity":1668,
             "same_knot_kernel_histogram":histogram,"results":results,
             "first16_kernel_counts":[r["A5_quotient_kernels"] for r in results[:16]],
             "three_trefoil_sum_control_kernels":49,"two_trefoil_sum_control_kernels":1,
             "twenty_free_group_maps_checked":True,"automorphic_duplicates_count_as_one":True,
             **controls,
             "scope":"Independent PSL2/PGL2(F5) tuple-orbit catalogue and crossing-closure enumeration reproduce every per-knot count. Permutation relators, onto images, Aut inequivalence and four-seed traces checked directly. PD parser shared. No maps pooled across knots; no new sharp target case.",
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/"joint_a5_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","cases":293,"generating_pairs":pair_count,"Aut":120,
                      "same_knot_kernel_histogram":histogram,"three_trefoil_sum_control_kernels":49},indent=2))


if __name__=="__main__":
    main()
