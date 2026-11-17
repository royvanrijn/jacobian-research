#!/usr/bin/env python3
"""Check index-five cover homology by evaluated Fox prefixes and Smith replay.

The prefix-derivative matrix does not import the producer's edge-walk
builder. The PD parser and the existing elementary Smith replay are shared.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram
from verify_integer_lifts import nullspace, replay

HERE = Path(__file__).resolve().parent


def multiply(a,b):
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a):
    return tuple(a.index(i) for i in range(len(a)))


def fox_cover_matrix(images, words, tree):
    r, s = len(images[0]), len(images)
    identity = tuple(range(r))
    # Verify the supplied deleted edges form a spanning tree.
    components = [{i} for i in range(r)]
    assert len(tree) == len(set(tree)) == r-1
    for edge in tree:
        generator, start = divmod(edge,r)
        end = inverse(images[generator])[start]
        a = next(c for c in components if start in c)
        b = next(c for c in components if end in c)
        assert a is not b
        components = [c for c in components if c is not a and c is not b]+[a|b]
    assert components == [set(range(r))]
    retained = [i for i in range(r*s) if i not in tree]
    rows = []
    for word in words:
        prefix, derivatives = identity, []
        for letter in word:
            generator = abs(letter)-1
            if letter < 0:
                prefix = multiply(prefix,inverse(images[generator]))
            derivatives.append((generator,1 if letter>0 else -1,inverse(prefix)))
            if letter > 0:
                prefix = multiply(prefix,images[generator])
        assert prefix == identity, "not a relator for the stated representation"
        block = [[0]*(r*s) for _ in range(r)]
        for generator,sign,coefficient in derivatives:
            for sheet in range(r):
                block[sheet][r*generator+coefficient[sheet]] += sign
        for row in block:
            boundary = [0]*r
            for edge,value in enumerate(row):
                generator,start = divmod(edge,r)
                boundary[start] -= value
                boundary[inverse(images[generator])[start]] += value
            assert not any(boundary)
        rows.extend([row[i] for i in retained] for row in block)
    return rows, len(retained)


def generated_group(images):
    group,queue = {tuple(range(5))},[tuple(range(5))]
    while queue:
        element = queue.pop()
        for generator in images:
            product = multiply(element,generator)
            if product not in group:
                group.add(product)
                queue.append(product)
    return group


def controls():
    from audit_core_subgroup_homology import covering_matrix
    from integer_lifts import smith_certificate
    images = [(1,2,0,3,4),(0,1,3,4,2)]
    assert len(generated_group(images)) == 60
    checked = []
    for label,words,expected_free,expected_torsion in [
        ("F2 index5 Schreier control",[],6,[]),
        ("C3*C3 index5 Kurosh control",[[1,1,1],[2,2,2]],0,[3]*4),
    ]:
        matrix,tree = covering_matrix(images,words)
        independent,columns = fox_cover_matrix(images,words,tree)
        assert matrix == independent and columns == 6
        smith = smith_certificate(matrix,columns)
        replay(independent,smith)
        assert columns-len(smith["diagonal"]) == expected_free
        assert [d for d in smith["diagonal"] if d>1] == expected_torsion
        checked.append(label)
    try:
        fox_cover_matrix(images,[[1]],tree)
    except AssertionError as error:
        assert "not a relator" in str(error)
    else:
        raise AssertionError("invalid relator accepted")
    return checked


def main():
    report = json.loads((HERE/"core_subgroup_homology_results.json").read_text())
    for path,digest in report["input_sha256"].items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest() == digest
    panel = json.loads((HERE/"parabolic_panel.json").read_text())
    prior = json.loads((HERE/"nonabelian_core_results.json").read_text())
    # Independent literal enumeration of A5, without the producer's core table.
    from itertools import permutations
    a5 = [p for p in permutations(range(5)) if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0]
    results,skipped = [],[]
    assert len(report["results"]) == len(prior["results"]) == len(panel["cases"]) == 16
    for case,old,result in zip(panel["cases"],prior["results"],report["results"]):
        assert case["name"] == old["name"] == result["name"]
        if old["first_nonconstant_colouring"] is None:
            assert old["normalized_colourings"] == 1 and result["status"] == "NO_NONCONSTANT_A5_MAP"
            skipped.append(case["name"])
            continue
        images = [a5[i] for i in old["first_nonconstant_colouring"]["arc_labels"]]
        assert [list(p) for p in images] == result["generator_images"]
        assert len(generated_group(images)) == 60
        crossings = oriented_diagram(case["pd"])
        for over,incoming,outgoing,_ in crossings:
            assert images[outgoing] == multiply(multiply(images[over],inverse(images[incoming])),images[over])
        root = old["seeds"][0]
        assert images[root] == tuple(range(5))
        words = [[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in crossings]+[[root+1]]
        assert words == result["presentation_words"]
        matrix,columns = fox_cover_matrix(images,words,result["tree_edges"])
        diagonal = replay(matrix,result["smith"])
        assert columns == result["columns"] == 61 and len(matrix) == result["rows"] == 70
        free,torsion = columns-len(diagonal),[d for d in diagonal if d>1]
        assert free == result["H1_free_rank"] == 0
        assert torsion == result["H1_torsion_factors"] and len(torsion) == 1
        assert result["maximum_H1_dimension_all_fields"] == 1
        assert not result["sharp_four_via_schreier"]
        for prime in [2,3,5,7,11]:
            assert len(nullspace(matrix,columns,prime)) == free+sum(d%prime==0 for d in torsion)
        results.append({"name":case["name"],"A5_image_order":60,"degree":5,
                        "H1_free_rank":free,"H1_torsion_factors":torsion,
                        "maximum_H1_dimension_all_fields":1,"sharp_four":False})
    assert len(results) == 12 and len(skipped) == 4
    inputs = ["verify_core_subgroup_homology.py","audit_core_subgroup_homology.py",
              "core_subgroup_homology_plan.json","core_subgroup_homology_results.json",
              "parabolic_panel.json","nonabelian_core_results.json","parabolic_lifts.py",
              "verify_integer_lifts.py","integer_lifts.py"]
    receipt = {"status":"PASS","covers":len(results),"results":results,"skipped_trivial":skipped,
               "controls":controls(),"invalid_relator_rejected":True,
               "scope":"Independent permutation relations, Fox-prefix chain matrix and elementary Smith replay on exactly the retained first A5 witnesses. PD parser and Smith producer/replay shared; no exclusion of other witnesses or covers. No new sharp case.",
               "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/"core_subgroup_homology_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","covers":len(results),"maximum_H1_dimension_all_fields":1,
                      "torsion_factor_counts":dict(Counter(r["H1_torsion_factors"][0] for r in results)),
                      "new_sharp_cases":0},indent=2))


if __name__ == "__main__":
    main()
