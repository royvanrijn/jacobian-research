#!/usr/bin/env python3
"""Check the bounded generator probe's inputs and certificate mechanism.

No failed heuristic run is promoted to a rank or non-primitivity theorem.
"""

import hashlib
import json
from pathlib import Path
import time

from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def inverse(word):
    return [-i for i in word[::-1]]


def reduce_word(word):
    result=list(word)
    changed=True
    while changed:
        changed=False
        for i in range(len(result)-1):
            if result[i]==-result[i+1]:
                del result[i:i+2]
                changed=True
                break
    return result


def cyclic(word):
    result=reduce_word(word)
    while len(result)>1 and result[0]==-result[-1]:
        result=result[1:-1]
    return result


def substitute(word,images):
    expanded=[]
    for letter in word:
        expanded.extend(images[letter-1] if letter>0 else inverse(images[-letter-1]))
    return reduce_word(expanded)


def main():
    report=json.loads((HERE/"core_generator_results.json").read_text())
    plan=json.loads((HERE/"core_generator_plan.json").read_text())
    for receipt in [report,plan]:
        for path,digest in receipt["input_sha256"].items():
            assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    panel=json.loads((HERE/"joint_a5_panel.json").read_text())["cases"][:16]
    sources=json.loads((HERE/"joint_a5_results.json").read_text())["results"][:16]
    assert report["status"]=="COMPLETE" and len(report["results"])==16
    autos=report["automorphisms"]
    assert len(autos)==96
    for auto in autos:
        for i in range(3):
            assert substitute(auto["images"][i],auto["inverse_images"])==[i+1]
            assert substitute(auto["inverse_images"][i],auto["images"])==[i+1]
    for case,source,result in zip(panel,sources,report["results"]):
        assert case["name"]==source["name"]==result["name"]
        assert source["seeds"]==result["seeds"] and source["seed_trace"]==result["seed_trace"]
        words={source["seeds"][0]:[]}
        words.update({g:[i+1] for i,g in enumerate(source["seeds"][1:])})
        for over,a,b in source["seed_trace"]:
            words[b]=reduce_word(words[over]+inverse(words[a])+words[over])
        assert [words[i] for i in range(13)]==result["arc_words"]
        relators=[reduce_word(words[o]+inverse(words[a])+words[o]+inverse(words[b])) for o,a,b,_ in oriented_diagram(case["pd"])]
        assert relators==result["three_generator_relators"]
        assert result["status"]=="BOUNDED_MISS" and result["candidates"]==444 and result["certificate"] is None
    from audit_core_generators import primitive_reduction
    positive=[]
    for word in [[1,2,1,3],[1,2,3,2,1]]:
        certificate,_=primitive_reduction(word,autos,time.monotonic()+2)
        assert certificate is not None
        value=cyclic(word)
        for index in certificate["automorphism_indices"]:
            value=cyclic(substitute(value,autos[index]["images"]))
        assert len(value)==1 and value==certificate["final_relator"]
        positive.append({"word":word,"certificate":certificate})
    for word in [[1,1,1],[2,2,2],[3,3,3],[1,2,-1,-2],[1,1,1,2,2,2]]:
        certificate,_=primitive_reduction(word,autos,time.monotonic()+2)
        assert certificate is None
    inputs=list(report["input_sha256"])+["core_generator_results.json","core_generator_plan.json","verify_core_generator_probe.py",
                                            "audit_core_generators.py","parabolic_lifts.py"]
    receipt={"status":"PASS scoped exploratory-input and control check","cases":16,"bounded_candidates_per_case":444,
             "positive_controls":positive,"negative_controls":5,"actual_two_generator_certificates":0,
             "rank_conclusions":[],"scope":"Independent free-word substitution and inverse-automorphism checks. The failed candidate search is not independently rerun or promoted to an impossibility theorem. The group ranks remain unknown.",
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in dict.fromkeys(inputs)}}
    (HERE/"core_generator_probe_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":receipt["status"],"cases":16,"two_generator_certificates":0,"rank_conclusions":[]},indent=2))


if __name__=="__main__":
    main()
