#!/usr/bin/env python3
"""Freeze all293 retained uncertainty names for same-knot A5 kernel counts."""

import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    names_path="sources/current13/data/still_open.json"
    pd_path="sources/knotinfo/2026.11.2_selected.json"
    names=json.loads((HERE/names_path).read_text())
    source=json.loads((HERE/pd_path).read_text())
    catalogue=json.loads((HERE/"joint_a5_catalogue.json").read_text())
    assert len(names)==len(set(names))==293
    assert catalogue["status"]=="PASS" and len(catalogue["generating_triple_orbits"])==1668
    rows={r["name"]:r for r in source["records"]}
    cases=[]
    for name in names:
        row=rows[name]
        assert row["bridge_index"]=="[3,4]" and len(row["pd_notation"])==13
        cases.append({"name":name,"pd":row["pd_notation"],"reported_bridge_interval":[3,4]})
    plan={"status":"FROZEN_BEFORE_KERNEL_TEST","selection":"All293 names in the pinned still_open list, unchanged order; no selected successes or replacements",
          "method":"Normalize one core seed to identity. Test all1668 Aut(A5) classes of generating triples on the other three meridional seeds. Each accepted class is exactly one quotient kernel of this one normalized double-cover group.",
          "threshold":{"A5_generating_pair_Aut_orbits":19,"distinct_kernels_for_rank3":20,"resulting_meridional_lower_bound":4},
          "motivation":"The user's same-knot joint-kernel criterion strengthens the earlier raw Hom count when multiple independent A5 quotients coexist. Full Aut orbit reduction makes an exact293-case audit bounded and cheap. Quotients of different knots are never combined.",
          "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in [names_path,pd_path,"joint_a5_catalogue.json","freeze_joint_a5_panel.py"]},
          "limits":{"cases":293,"seed_count":4,"seed_trials_per_case":715,"epimorphism_orbits_per_case":1668,
                    "maximum_orbit_trials":293*1668,"seconds_per_case":10,"seconds_per_run":300,"checkpoint_every_cases":1,"no_adaptive_replacements":True},
          "acceptance":"A positive requires20 onto maps of one fixed core presentation, pairwise inequivalent under all120 automorphisms, and a four-seed upper trace. Independent group and relation replay plus named-diagram and prior-family checks remain required.",
          "claim_boundary":"A count below20 only fails this particular A5 rank3 certificate. It does not prove rank2, bridge3, universal finite-quotient detection or full MRC. The source uncertainty list is not a novelty proof.",
          "cases":cases}
    path=HERE/"joint_a5_panel.json"
    if path.exists():
        assert json.loads(path.read_text())==plan
    else:
        path.write_text(json.dumps(plan,indent=2)+"\n")
    print(json.dumps({"status":"FROZEN","cases":len(cases),"orbit_trials":plan["limits"]["maximum_orbit_trials"],
                      "name_prefixes":{p:sum(c["name"].startswith(p) for c in cases) for p in ["13a","13n"]}},indent=2))


if __name__=="__main__":
    main()
