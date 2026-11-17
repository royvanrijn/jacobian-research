#!/usr/bin/env python3
"""Check retained open-list membership and the original Dutra gap dependency.

This is a scoped source/dependency audit, not a replay of either full paper.
The four-letter common word is reduced directly in the stated free group.
"""

import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def reduce_word(word):
    stack=[]
    for letter in word:
        assert letter in [-2,-1,1,2]
        if stack and stack[-1]==-letter:
            stack.pop()
        else:
            stack.append(letter)
    return stack


def main():
    observations=json.loads((HERE/"source_status_observations.json").read_text())
    blob=(HERE/observations["open_list"]["retained_bytes"]).read_bytes()
    git_sha=hashlib.sha1(b"blob "+str(len(blob)).encode()+b"\0"+blob).hexdigest()
    assert git_sha==observations["open_list"]["git_blob_sha"]
    open_names=json.loads(blob)
    assert len(open_names)==len(set(open_names))==293
    panel=json.loads((HERE/"parabolic_panel.json").read_text())
    names=[c["name"] for c in panel["cases"]]
    certificates=json.loads((HERE/"sources/current13/data/certificates.json").read_text())
    closed={c["name"] for c in certificates}
    assert len(closed)==243 and not (closed&set(open_names))
    assert len(names)==16 and set(names)<=set(open_names)
    # Dutra p7: pi1(E0)=F(x1,x2), with the two displayed embeddings.
    alpha=[[2,-1,-2],[1]]
    omega=[[2,-1,-2,1,-2],[2]]
    alpha_product=reduce_word(alpha[0]+alpha[1])
    omega_product=reduce_word(omega[0]+omega[1])
    assert alpha_product==omega_product==[2,-1,-2,1]
    assert reduce_word(alpha_product)==alpha_product and len(alpha_product)==4
    registry=json.loads((ROOT/"MATH_STATUS.json").read_text())
    dependencies=[]
    for entry in registry["entries"]:
        for dependency in entry.get("dependencies",[]):
            if any(token in str(dependency).lower() for token in ["2212.13081","dut22","dutra2022","dutra 2022"]):
                dependencies.append({"id":entry["id"],"state":entry["state"],"dependency":dependency})
    own_dependencies=[r for r in dependencies if r["id"].startswith("MRC-")]
    assert not own_dependencies
    sources=json.loads((HERE/"sources.json").read_text())
    dutra=next(s for s in sources["sources"] if s["id"]=="DUT22")
    assert "gap" in dutra["status"].lower()
    inputs=["audit_whitehead_dependencies.py","source_status_observations.json","parabolic_panel.json",
            "sources/current13/data/still_open.json","sources/current13/data/certificates.json",
            "sources/dutra-whitehead-v1.pdf","sources/whitehead-v1.html","sources.json"]
    result={"status":"PASS scoped source and dependency audit","current_open_list_blob":git_sha,
            "all16_targets_still_in_open_list":True,"open_names":293,"separate_claimed_closures":243,
            "Dutra_original_v1_gap":{"free_group_rank":2,"alpha_generator_words":alpha,"omega_generator_words":omega,
                "common_nontrivial_reduced_word":alpha_product,"Lemma3_1_counterexample_conjugator":[],
                "source_pages_visually_checked":[7,12],"dependency":"Lemma3.8(2) explicitly invokes Lemma3.1 on p12",
                "conclusion":"Lemma3.1 as stated is false; this does not disprove the main MRC conclusion"},
            "direct_DUT22_dependencies":dependencies,"own_proved_entries_using_DUT22":own_dependencies,
            "retained_family_disposition":"DUT22 remains a preprint with a checked proof gap, not an established family. The2017 graph-manifold/braid-pattern theorem is a different earlier source and is retained. No source assumption silently upgraded.",
            "GGLW_scope_read":["Introduction footnote and main statements","Section2.3 self-dual doubling and nonvanishing conditions","Section2.4 actual peripheral correspondence","Definition4.1, selected valuation steps of Theorem4.2 and Corollary4.5","Earlier complete Lemma6.2 proof retained"],
            "GGLW_remaining":"Full44-page microlocal, character-curve and two-bridge existence proof not independently audited. No certificate for one of the16 targets is inferred.",
            "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/"whitehead_dependency_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"all16_still_open":True,"original_Dutra_common_word":alpha_product,
                      "own_DUT22_dependencies":len(own_dependencies)},indent=2))


if __name__=="__main__":
    main()
