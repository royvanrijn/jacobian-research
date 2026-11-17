#!/usr/bin/env python3
"""Freeze a crossing-budget Montesinos identity audit on the unchanged16.

Candidate construction uses source determinants, never named-diagram
matches. A miss is not an exclusion of Montesinos membership.
"""

from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import prod
from pathlib import Path

HERE=Path(__file__).resolve().parent


def continued_fractions(cost):
    yield [cost]
    for first in range(1,cost-1):
        for rest in continued_fractions(cost-first):
            yield [first]+rest


def slope(cf):
    value=Fraction(cf[-1])
    for a in reversed(cf[:-1]):
        value=a+1/value
    return 1/value


def cost_tuples(total,count):
    if count==1:
        if total>=2:
            yield [total]
        return
    for first in range(2,total-2*(count-1)+1):
        for tail in cost_tuples(total-first,count-1):
            yield [first]+tail


def build_candidates(determinants):
    fractions={}
    for cost in range(2,10):
        terms=[]
        for cf in continued_fractions(cost):
            value=slope(cf)
            assert 0<value<1 and cf[-1]>=2 and sum(cf)==cost
            terms.append({"slope":[value.numerator,value.denominator],"continued_fraction":[0]+cf,"cost":cost})
        assert len(terms)==len({tuple(t["slope"]) for t in terms})==2**(cost-2)
        fractions[cost]=sorted(terms,key=lambda t:t["slope"])
    candidates=[]
    totals=Counter()
    for branches in [3,4]:
        for integer in range(14-2*branches):
            for costs in cost_tuples(13-integer,branches):
                for terms in product(*(fractions[c] for c in costs)):
                    totals[branches]+=1
                    values=[Fraction(*t["slope"]) for t in terms]
                    determinant=prod(x.denominator for x in values)*(integer+sum(values,Fraction(0)))
                    assert determinant.denominator==1
                    if int(determinant) in determinants:
                        candidates.append({"branches":branches,"integer_twists":integer,
                                           "fractions":[t["slope"] for t in terms],
                                           "continued_fractions":[t["continued_fraction"] for t in terms],
                                           "tangle_costs":costs,"determinant_screen":int(determinant)})
    assert totals=={3:7423,4:2561}
    return candidates,dict(totals)


def main():
    source="parabolic_panel.json"
    determinants_source="branched_cover_verification.json"
    panel=json.loads((HERE/source).read_text())
    lower=json.loads((HERE/determinants_source).read_text())
    determinants={c["name"]:c["degrees"][0]["H1_torsion_factors"][0] for c in lower["all_degree_field_panel"]}
    assert set(determinants)=={c["name"] for c in panel["cases"]}
    candidates,totals=build_candidates(set(determinants.values()))
    plan={"status":"FROZEN_BEFORE_IDENTITY_TEST","targets":[c["name"] for c in panel["cases"]],
          "source_panel":source,"source_determinants":determinants_source,
          "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in [source,determinants_source,"freeze_current13_montesinos.py"]},
          "target_determinants":determinants,"crossing_budget":13,"branches":[3,4],
          "candidate_rule":"Ordered positive proper rational tangles with unique positive continued fractions ending>=2, every tangle cost>=2, nonnegative integer twist e, e+sum(costs)=13. Keep tuples whose product(denominators)*(e+sum(slopes)) equals a target determinant. Test all individual planar half turns. No flype or Reidemeister search.",
          "pre_screen_tuple_counts":totals,"candidate_count":len(candidates),"candidates":candidates,
          "diagram_count":sum(2**c["branches"] for c in candidates),
          "limits":{"maximum_pre_screen_tuples":9984,"maximum_diagrams":100360,"maximum_run_seconds":300,"checkpoint_every_candidates":25,"no_adaptive_replacements":True},
          "motivation":"The same16 current targets retain lower3/upper4 bounds after all-field symmetric and fixed H3 lift failures. Exact three- or four-tangle Montesinos identification would give matching mr=bridge=3 or4 via the classical Boileau-Zieschang theorem. The crossing-budget and determinant screen give a finite geometric test before changing target knots or coefficient fields.",
          "acceptance":"Exact full crossing-preserving planar-map identity up to global ambient reflection; retain PD, fraction data and witness. Positive matches need separate paired-endpoint construction and seed/theorem-scope verification.",
          "claim_boundary":"A negative excludes only these marked candidate diagrams. It does not exclude Montesinos membership after flypes, isotopies, mixed signs, nonminimal diagrams or other presentations. An identification uses a known MRC theorem, so no new family theorem or novelty claim follows."}
    path=HERE/"current13_montesinos_plan.json"
    if path.exists():
        assert json.loads(path.read_text())==plan
    else:
        path.write_text(json.dumps(plan,indent=2)+"\n")
    print(json.dumps({"status":"FROZEN","targets":len(plan["targets"]),"pre_screen_tuples":totals,
                      "selected_tuples":len(candidates),"diagrams":plan["diagram_count"],
                      "selected_by_branches":dict(Counter(c["branches"] for c in candidates))},indent=2))


if __name__=="__main__":
    main()
