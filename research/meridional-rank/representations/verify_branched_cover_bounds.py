#!/usr/bin/env python3
"""Compare full cover-edge homology with root lifts and known controls."""

from collections import Counter
from itertools import combinations, product
import hashlib
import json
from math import factorial
from pathlib import Path

from branched_cover_homology import cover_homology, betti_in_prime
from integer_lifts import relator_matrix, smith_certificate
from parabolic_lifts import oriented_diagram
from verify_integer_lifts import replay

HERE = Path(__file__).resolve().parent


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def all_assignments(degree):
    permutations = []
    for a, b in combinations(range(degree), 2):
        p = list(range(degree))
        p[a], p[b] = b, a
        permutations.append(tuple(p))
    identity = tuple(range(degree))
    assignments = []
    for tail in product(range(len(permutations)), repeat=3):
        assignment = (0,) + tail
        gens = [permutations[i] for i in assignment]
        group, todo = {identity}, [identity]
        while todo:
            x = todo.pop()
            for y in gens:
                value = compose(x, y)
                if value not in group:
                    group.add(value)
                    todo.append(value)
        if len(group) == factorial(degree):
            assignments.append(assignment)
    conjugations = [[permutations.index(compose(compose(a,b),a)) for b in permutations] for a in permutations]
    return assignments, conjugations


def lift_rank_check(data, degree, transpositions, words):
    cover = cover_homology(degree, transpositions, words)
    replay(cover["presentation"], cover["smith"])
    matrix = relator_matrix(degree, transpositions, words)
    scalar = smith_certificate(matrix, len(transpositions))
    replay(matrix, scalar)
    torsion = [v for v in scalar["diagonal"] if v > 1]
    assert cover["H1_torsion_factors"] == torsion
    assert scalar["characteristic_zero_dimension"] == degree - 1 + cover["H1_free_rank"]
    for p in [2,3,5,7,11]:
        dimension = len(transpositions) - sum(d % p != 0 for d in scalar["diagonal"])
        assert dimension == degree - 1 + betti_in_prime(cover,p)
    return {"name": data, **cover, "scalar_diagonal": scalar["diagonal"]}


def main():
    controls = []
    for t in range(1,5):
        # Standard amalgamated presentations of connected sums of t trefoils.
        words = [[1,j+2,1,-(j+2),-1,-(j+2)] for j in range(t)]
        two = lift_rank_check(f"sum_{t}_trefoils_double",2,[(0,1)]*(t+1),words)
        assert two["H1_free_rank"] == 0 and two["H1_torsion_factors"] == [3]*t
        three = lift_rank_check(f"sum_{t}_trefoils_simple_three",3,[(0,1)]+[(0,2)]*t,words)
        assert three["H1_free_rank"] == t-1 and three["H1_torsion_factors"] == []
        controls.extend([two,three])
    sharp = json.loads((HERE/"sharp_five_certificates.json").read_text())["certificates"]
    from search_lift_panel import diagram
    case = sharp[0]
    _, cr = diagram(case["gauss_code"])
    words = [[o+1,a+1,o+1,b+1] for o,a,b in cr]
    four = lift_rank_check(case["name"],4,case["arc_transpositions"],words)
    assert four["H1_free_rank"] == 1 and four["H1_torsion_factors"] == [5]
    controls.append(four)

    panel = json.loads((HERE/"parabolic_panel.json").read_text())
    report = json.loads((HERE/"all_symmetric_degrees_results.json").read_text())
    for path,digest in report["input_sha256"].items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    tables = {r:all_assignments(r) for r in [2,3,4,5]}
    results = []
    assert len(panel["cases"]) == len(report["results"]) == 16
    for case, result in zip(panel["cases"],report["results"]):
        assert case["name"] == result["name"]
        cr = [c[:3] for c in oriented_diagram(case["pd"])]
        seen = set(result["seeds"])
        assert len(seen) == 4
        for o,old,new in result["seed_trace"]:
            assert o in seen and old in seen and new not in seen
            assert (o,old,new) in cr or (o,new,old) in cr
            seen.add(new)
        assert len(seen) == len(cr)
        degrees = []
        for degree_result in result["degrees"]:
            r = degree_result["degree"]
            assignments, conjugations = tables[r]
            accepted = []
            for assignment in assignments:
                labels = [None]*len(cr)
                for seed,label in zip(result["seeds"],assignment):
                    labels[seed] = label
                for o,old,new in result["seed_trace"]:
                    labels[new] = conjugations[labels[o]][labels[old]]
                if all(labels[b] == conjugations[labels[o]][labels[a]] for o,a,b in cr):
                    accepted.append((assignment,labels))
            assert [(tuple(c["assignment"]),c["arc_labels"]) for c in degree_result["certificates"]] == accepted
            assert len(accepted) == (1 if r==2 else 0)
            if r==2:
                words = [[o+1,a+1,o+1,b+1] for o,a,b in cr]
                data = lift_rank_check(case["name"],2,[(0,1)]*len(cr),words)
                assert data["H1_free_rank"] == 0 and len(data["H1_torsion_factors"]) == 1
                assert data["H1_torsion_factors"][0] % 3 != 0
                cert = degree_result["certificates"][0]["smith"]
                rows = []
                for o,a,b in cr:
                    row = [0]*len(cr)
                    row[o] += 2
                    row[a] -= 1
                    row[b] -= 1
                    if any(row):
                        rows.extend([row,[-v for v in row]])
                replay(rows,cert)
                assert cert["maximum_dimension_all_characteristics"] == 2
                degrees.append({"degree":r,"H1_free_rank":0,"H1_torsion_factors":data["H1_torsion_factors"],"maximum_lift_dimension":2})
            else:
                degrees.append({"degree":r,"surjective_bases":0})
        assert result["maximum_scalar_lift_bound_all_degrees_fields"] == 2
        results.append({"name":case["name"],"degrees":degrees,
                        "all_dimensions_binary_affine_criterion":"EXCLUDED by bridge<=4, k+1 bound, and no S3 transposition base"})
    inputs = ["verify_branched_cover_bounds.py","branched_cover_homology.py","all_symmetric_degrees_results.json",
              "all_symmetric_degrees_plan.json","parabolic_panel.json","sharp_five_certificates.json",
              "integer_lifts.py","verify_integer_lifts.py","parabolic_lifts.py"]
    output = {"status":"PASS", "all_degree_field_panel":results,"homology_controls":controls,
              "uniform_claim_assurance":"Written proofs, not inferred from finite checks; no formal or external review or priority claim",
              "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/"branched_cover_verification.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"status":"PASS","homology_controls":len(controls),"panel_cases":len(results),
                      "all_degree_field_maximum_bound":2,"new_sharp_knot_cases":0},indent=2))


if __name__ == "__main__":
    main()
