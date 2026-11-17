"""Scoped authority and evidence checks for the36-term integral twelve-variable map."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=HERE/"integration12.json"
    if not out.exists():out.write_text('{"status":"pending"}\n')
    known={e["id"]:e for e in json.loads((ROOT/"MATH_STATUS.json").read_text())["entries"]}
    entry=known["KCUB3"]
    assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
    assert entry["state"]=="proved" and entry["programme_status"]=="archived"
    assert entry["dependencies"]==["BCRDB1"] and known["BCRDB1"]["state"]==known["BCR2"]["state"]=="proved"
    assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
    assert all((ROOT/p).is_file() for p in [entry["canonical_source"],entry["checker"],*entry["software_lock"]])
    assert entry["independent_replay"] and not entry["formal_verification"] and not entry["external_review"]
    assert not entry["closes_problems"] and not entry["narrows_problems"]
    for key in ["KCUB1","KCUB2"]:
        assert known[key]["state"]=="proved" and known[key]["replaced_by"]==["KCUB3"]
        assert "BCRDB1" in known[key]["scope"] and "dominat" in known[key]["scope"]
    v=json.loads((HERE/"integer-validation12.json").read_text())
    assert v["status"]=="pass" and (v["dimension"],v["degree"],v["nonlinear_total"])==(12,3,36)
    assert v["baseline_nonlinear_total"]==40
    assert v["source_sha256"]==sha(HERE/"verify_integer12.py")
    assert v["arithmetic_sha256"]==sha(HERE/"verify.py")
    assert v["rational_checker_sha256"]==sha(HERE/"verify_sparse12.py")
    assert v["primary_source_sha256"]==sha(HERE/"construct_integer12.py")
    assert v["witness_sha256"]==sha(HERE/"cubic12-integer-witness.json")
    assert v["rational_validation_sha256"]==sha(HERE/"sparse-validation12.json")
    assert all(v[k] for k in ["integral_coefficients","integral_collision_points","determinant_one",
        "identity_derivative","full_polynomial_determinant","all_variables_symbolic","source_and_target_conjugation_replayed"])
    q=json.loads((HERE/"sparse-validation12.json").read_text())
    assert q["source_sha256"]==sha(HERE/"verify_sparse12.py")
    assert q["primary_source_sha256"]==sha(HERE/"construct_sparse12.py")
    assert q["witness_sha256"]==sha(HERE/"cubic12-sparse-witness.json")
    base=json.loads((HERE/"de-bondt12-witness.json").read_text())
    assert q["base_sha256"]==sha(HERE/"de-bondt12-witness.json")
    assert q["base_source_sha256"]==base["source_sha256"]==sha(ROOT/base["source"])
    descent=json.loads((HERE/"linear-descent12.json").read_text())
    assert descent["source_sha256"]==sha(HERE/"optimize12.py")
    assert descent["input_sha256"]==sha(HERE/"de-bondt12-witness.json")
    assert descent["initial_support"]==40 and descent["final_support"]==36
    for filename,key in [("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")]:
        for item in json.loads((ROOT/filename).read_text())[key]:
            for identifier in ["KCUB1","KCUB2","KCUB3","BCR2","BCRDB1","BCRDBGRAPH1"]:
                if identifier in item.get("reviewed_claims",{}):
                    assert item["reviewed_claims"][identifier]==claim_fingerprint(known[identifier])
    for name in ["README.md","PROOF12.md","BASELINE_CORRECTION.md","PROOF13.md","PROOF.md","EXPLORATION.md"]:
        for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
            if not target.startswith(("http:","https:","#")):assert (HERE/target.split("#")[0]).is_file(),(name,target)
    seen=set()
    def visit(identifier,path):
        assert identifier not in path
        if identifier in seen:return
        for dep in known[identifier]["dependencies"]:
            if dep in known:visit(dep,path|{identifier})
        seen.add(identifier)
    visit("KCUB3",set())
    assert known["OP-GLDC"]["state"]=="open"
    interval=json.loads((HERE.parent/"collision-frame-search/validation.json").read_text())["normalized_support_interval"];assert interval==[10,48]
    receipt=dict(status="pass",claim="KCUB3",dimension=12,degree=3,baseline_nonlinear_occurrences=40,nonlinear_occurrences=36,
        integral_coefficients=True,integral_collision_points=True,independent_full_polynomial_verification=True,
        new_dimension_record=False,dimension_optimality=False,support_optimality=False,
        known_twelve_variable_baseline_restored=True,dominated_experiments=["KCUB1","KCUB2"],
        normalized_three_variable_support_interval=interval,umbrella_status="open",
        schema=True,source_hashes=True,links=True,reviewed_fingerprints=True,dependency_graph_acyclic=True,
        validation_sha256=sha(HERE/"integer-validation12.json"),
        navigation_errors=["render_status.py: unrelated HC4MYGS1 invalid kind",
           "research.py render: METHOD-EC-CENTRE-POLICY missing source artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json"])
    out.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
if __name__=="__main__":main()

