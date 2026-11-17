#!/usr/bin/env python3
"""Fixed index-five cover homology audit from the retained A5 witnesses.

Freeze first. Each existing first nonconstant witness is used once;
no witness, group, field or knot replacement follows a negative result.
"""

import hashlib
import json
from pathlib import Path
import time

from integer_lifts import smith_certificate
from parabolic_lifts import oriented_diagram
from search_nonabelian_core import a5_tables

HERE = Path(__file__).resolve().parent


def presentation(case, result, group):
    labels = result["first_nonconstant_colouring"]["arc_labels"]
    root = result["seeds"][0]
    assert labels[root] == 0
    images = [group[i] for i in labels]
    words = [[o+1, -(a+1), o+1, -(b+1)]
             for o, a, b, _ in oriented_diagram(case["pd"])]
    words.append([root+1])
    return images, words


def covering_matrix(images, words):
    r, s = len(images[0]), len(images)
    inverses = [tuple(p.index(i) for i in range(r)) for p in images]
    # Positive edges realize the RIGHT monodromy action j.g=phi(g)^(-1)j.
    parent = list(range(r))

    def root(v):
        while v != parent[v]:
            v = parent[v]
        return v

    tree = []
    for i, p in enumerate(inverses):
        for v, w in enumerate(p):
            if root(v) != root(w):
                tree.append(i*r+v)
                parent[root(v)] = root(w)
    assert len(tree) == r-1
    rows = []
    for word in words:
        for start in range(r):
            v = start
            row = [0]*(r*s)
            for letter in word:
                i = abs(letter)-1
                if letter > 0:
                    row[r*i+v] += 1
                    v = inverses[i][v]
                else:
                    v = images[i][v]
                    row[r*i+v] -= 1
            assert v == start
            rows.append(row)
    retained = [i for i in range(r*s) if i not in tree]
    matrix = [[row[i] for i in retained] for row in rows]
    return matrix, tree


def main():
    panel = json.loads((HERE/"parabolic_panel.json").read_text())
    prior = json.loads((HERE/"nonabelian_core_results.json").read_text())
    plan = {
        "targets": [c["name"] for c in panel["cases"]],
        "source_panel_sha256": hashlib.sha256((HERE/"parabolic_panel.json").read_bytes()).hexdigest(),
        "witness_source_sha256": hashlib.sha256((HERE/"nonabelian_core_results.json").read_bytes()).hexdigest(),
        "selection": "Exactly the first nonconstant A5 witness already retained per target; skip targets whose complete Hom count is1. No alternate witnesses or target replacements.",
        "cover_degree": 5, "maximum_covers": 12,
        "maximum_matrix_rows": 70, "maximum_matrix_columns": 61,
        "characteristics": "All, by integer Smith factors; no prime sweep",
        "sharp_four_homology_threshold": 7,
        "motivation": "For H=pi1(double branched cover) and J of index5, rank J<=5(rank H-1)+1. Thus dim H1(J;F)>=7 forces rank H>=3 and mr(K)>=4, matching the fixed four-seed upper certificate. The A5 maps supply exact index5 subgroups even when raw Hom counts are too small.",
        "claim_boundary": "Only the retained first A5 witnesses are tested. No-hit does not exclude other witnesses, covers, finite quotients, or nonabelian rank arguments.",
        "no_adaptive_replacements": True,
    }
    plan_path = HERE/"core_subgroup_homology_plan.json"
    if not plan_path.exists():
        plan_path.write_text(json.dumps(plan, indent=2)+"\n")
        print("FROZEN core_subgroup_homology_plan.json; rerun to execute")
        return
    assert json.loads(plan_path.read_text()) == plan
    inputs = ["audit_core_subgroup_homology.py", "core_subgroup_homology_plan.json",
              "parabolic_panel.json", "nonabelian_core_results.json", "search_nonabelian_core.py",
              "integer_lifts.py", "parabolic_lifts.py"]
    hashes = {p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    out = HERE/"core_subgroup_homology_results.json"
    if out.exists():
        report = json.loads(out.read_text())
        assert report["input_sha256"] == hashes
    else:
        report = {"input_sha256": hashes, "results": [], "claim_boundary": plan["claim_boundary"]}
    group, _ = a5_tables()
    done = {r["name"] for r in report["results"]}
    for case, result in zip(panel["cases"], prior["results"]):
        assert case["name"] == result["name"] and result["status"] == "COMPLETE"
        if case["name"] in done:
            continue
        if result["first_nonconstant_colouring"] is None:
            assert result["normalized_colourings"] == 1
            data = {"name": case["name"], "status": "NO_NONCONSTANT_A5_MAP"}
        else:
            started = time.monotonic()
            images, words = presentation(case, result, group)
            matrix, tree = covering_matrix(images, words)
            columns = len(matrix[0])
            assert len(matrix) <= 70 and columns <= 61
            smith = smith_certificate(matrix, columns)
            free = columns-len(smith["diagonal"])
            torsion = [d for d in smith["diagonal"] if d > 1]
            maximum = free+len(torsion)
            data = {"name": case["name"], "status": "COMPLETE", "degree":5,
                    "generator_images": images, "presentation_words":words,
                    "tree_edges":tree, "rows":len(matrix), "columns":columns,
                    "smith":smith, "H1_free_rank":free, "H1_torsion_factors":torsion,
                    "maximum_H1_dimension_all_fields":maximum,
                    "sharp_four_via_schreier":maximum>=7,
                    "seconds":round(time.monotonic()-started,6)}
        report["results"].append(data)
        out.write_text(json.dumps(report, indent=2)+"\n")
        print(json.dumps({k:data[k] for k in ["name", "status", "maximum_H1_dimension_all_fields", "H1_torsion_factors", "sharp_four_via_schreier"] if k in data}), flush=True)


if __name__ == "__main__":
    main()
