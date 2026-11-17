"""Scoped checks for KCUB4's nine arbitrary-degree target obstructions."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=HERE/"graph36-integration.json"
    if not out.exists():out.write_text('{"status":"pending"}\n')
    known={e["id"]:e for e in json.loads((ROOT/"MATH_STATUS.json").read_text())["entries"]};entry=known["KCUB4"]
    assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
    assert entry["state"]=="proved" and entry["programme_status"]=="archived" and entry["dependencies"]==["KCUB3"]
    assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
    assert all((ROOT/p).is_file() for p in [entry["canonical_source"],entry["checker"],*entry["software_lock"]])
    assert entry["independent_replay"] and not entry["formal_verification"] and not entry["external_review"]
    assert not entry["closes_problems"] and not entry["narrows_problems"]
    v=json.loads((HERE/"graph36-validation.json").read_text())
    assert v["status"]=="pass" and v["source_sha256"]==sha(HERE/"verify_graph36.py")
    assert v["arithmetic_sha256"]==sha(HERE/"verify.py") and v["input_sha256"]==sha(HERE/"cubic12-integer-witness.json")
    assert all(v[k] for k in ["all_nine_arbitrary_degree_certificates","normalized_raw_source_jets_checked",
        "dual_annihilation_and_nonzero_defect","commutative_target_monomial_coverage_checked",
        "zero_branch_pruning_coverage_checked","higher_degree_terms_excluded_by_positive_order",
        "known_BCR2_square_completion_control"])
    assert [c["pivot"] for c in v["arbitrary_target_degree_cases"]]==list(range(4,13))
    assert [c["maximum_jet"] for c in v["arbitrary_target_degree_cases"]]==[6,7,4,4,6,7,6,7,5]
    for j in range(4,13):
        assert v["case_hashes"][str(j)]["raw_sha256"]==sha(HERE/"graph36"/("pivot-%02d.json"%j))
        assert v["case_hashes"][str(j)]["arbitrary_degree_sha256"]==sha(HERE/"adaptive36-wide"/("pivot-%02d.json"%j))
    for filename,key in [("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")]:
        for item in json.loads((ROOT/filename).read_text())[key]:
            if "KCUB4" in item.get("reviewed_claims",{}):assert item["reviewed_claims"]["KCUB4"]==claim_fingerprint(entry)
    for name in ["GRAPH_DELETION_OBSTRUCTION.md","README.md"]:
        for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
            if not target.startswith(("http:","https:","#")):assert (HERE/target.split("#")[0]).is_file(),(name,target)
    assert known["OP-GLDC"]["state"]=="open"
    receipt=dict(status="pass",claim="KCUB4",fixed_graphs=9,arbitrary_target_polynomial_degree=True,
       one_raw_output_stage_only=True,general_dimension_exclusion=False,umbrella_status="open",
       schema=True,source_hashes=True,links=True,reviewed_fingerprints=True,
       validation_sha256=sha(HERE/"graph36-validation.json"))
    out.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
if __name__=="__main__":main()

