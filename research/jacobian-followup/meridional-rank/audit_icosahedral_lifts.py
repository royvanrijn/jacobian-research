#!/usr/bin/env python3
"""Freeze and execute all-field H3 root lifts of existing first A5 maps."""

import hashlib
import json
from pathlib import Path

from integer_lifts import smith_certificate
from parabolic_lifts import oriented_diagram
from quadratic_integer_lifts import h3_model,root_matrix,integral_expansion,factor_exact,residue_fields,kernel

HERE=Path(__file__).resolve().parent


def main():
    panel=json.loads((HERE/"parabolic_panel.json").read_text())
    audit=json.loads((HERE/"core_involution_audit.json").read_text())
    plan={"targets":[c["name"] for c in panel["cases"]],
          "source_panel_sha256":hashlib.sha256((HERE/"parabolic_panel.json").read_bytes()).hexdigest(),
          "involution_audit_sha256":hashlib.sha256((HERE/"core_involution_audit.json").read_bytes()).hexdigest(),
          "selection":"Exactly the existing first A5 witnesses and their unique deck implementers, no new maps or target replacement",
          "coefficient_ring":"Z[tau], tau^2=tau+1", "base":"H3 reflection matrices with Cartan edges5,3",
          "maximum_bases":12,"maximum_integer_matrix_rows":78,"maximum_integer_matrix_columns":26,
          "characteristics":"Every residue characteristic: expanded integer Smith factors isolate all rank-drop primes away from discriminant5; include5 separately and test every embedding",
          "factor_trial_cap":1000000,"residue_root_search_cap":1000000,
          "motivation":"All twelve retained first A5 core maps admit a deck involution implemented by a double transposition, giving rank3 H3 reflection bases. Restricted affine-root nullity4 would prove mr=bridge=4. This changes the obstructed symmetric base and uses an exact rank-drop criterion before coefficient tests.",
          "claim_boundary":"Only these fixed H3 bases are covered. A positive needs explicit independent field-kernel, meridian and upper-bound replay plus prior-family audit; a negative does not exclude other H3 colourings or other representations.",
          "no_adaptive_replacements":True}
    plan_path=HERE/"icosahedral_lift_plan.json"
    if not plan_path.exists():
        plan_path.write_text(json.dumps(plan,indent=2)+"\n")
        print("FROZEN icosahedral_lift_plan.json; rerun to execute")
        return
    assert json.loads(plan_path.read_text())==plan
    inputs=["audit_icosahedral_lifts.py","quadratic_integer_lifts.py","icosahedral_lift_plan.json",
            "core_involution_audit.json","parabolic_panel.json","integer_lifts.py","parabolic_lifts.py"]
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    out=HERE/"icosahedral_lift_results.json"
    model,root_vectors,group_certificate=h3_model()
    if out.exists():
        report=json.loads(out.read_text())
        assert report["input_sha256"]==hashes
    else:
        report={"input_sha256":hashes,"H3_group":group_certificate,"results":[],"claim_boundary":plan["claim_boundary"]}
    done={r["name"] for r in report["results"]}
    for case,data in zip(panel["cases"],audit["results"]):
        assert case["name"]==data["name"]
        if case["name"] in done:
            continue
        if not data.get("implementers"):
            result={"name":case["name"],"status":"NO_RETAINED_H3_BASE"}
        else:
            assert len(data["implementers"])==1
            images=data["implementers"][0]["meridian_A5_components"]
            keys=[(tuple(p),1) for p in images]
            matrices,roots=[model[k] for k in keys],[root_vectors[k] for k in keys]
            words=[[o+1,a+1,o+1,b+1] for o,a,b,_ in oriented_diagram(case["pd"])]
            scalar=root_matrix(matrices,roots,words)
            expanded=integral_expansion(scalar)
            n=len(case["pd"])
            assert len(expanded)<=78 and 2*n<=26
            smith=smith_certificate(expanded,2*n)
            assert len(smith["diagonal"])%2==0
            generic=n-len(smith["diagonal"])//2
            candidates=sorted(set([5]+factor_exact(smith["diagonal"][-1],plan["factor_trial_cap"])))
            if generic>=4:
                candidates=sorted(set(candidates+[2]))
            fields=[]
            for p in candidates:
                for field in residue_fields(p,plan["residue_root_search_cap"]):
                    basis=kernel(scalar,n,field)
                    fields.append({"prime":p,"tau_root":field.root,"field_degree":2 if field.root is None else 1,
                                   "dimension":len(basis),"basis":basis})
            maximum=max([generic]+[f["dimension"] for f in fields])
            result={"name":case["name"],"status":"COMPLETE","meridian_matrices":matrices,"root_vectors":roots,
                    "root_matrix":scalar,"expanded_smith":smith,"generic_dimension":generic,
                    "candidate_primes":candidates,"residue_fields":fields,"maximum_dimension_all_fields":maximum,
                    "sharp_four":maximum>=4}
        report["results"].append(result)
        out.write_text(json.dumps(report,indent=2)+"\n")
        print(json.dumps({k:result[k] for k in ["name","status","generic_dimension","candidate_primes","maximum_dimension_all_fields","sharp_four"] if k in result}),flush=True)


if __name__=="__main__":
    main()
