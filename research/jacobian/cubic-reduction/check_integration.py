"""Scoped metadata, source-hash and navigation audit for KCUB1."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=HERE/"integration.json"
    if not out.exists():out.write_text('{"status":"pending"}\n')
    ledger=json.loads((ROOT/"MATH_STATUS.json").read_text());known={e["id"]:e for e in ledger["entries"]}
    entry=known["KCUB1"]
    assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
    assert entry["state"]=="proved" and entry["programme_status"]=="archived"
    assert entry["checker"]=="jacobian/cubic-reduction/verify.py"
    assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
    assert entry["dependencies"]==["KFRAME1"]
    assert not entry["closes_problems"] and not entry["narrows_problems"]
    assert not entry["formal_verification"] and not entry["external_review"] and entry["independent_replay"]
    assert all((ROOT/p).is_file() for p in [entry["canonical_source"],entry["checker"],*entry["software_lock"]])
    validation=json.loads((HERE/"validation.json").read_text())
    assert validation["status"]=="pass" and validation["dimension"]==17 and validation["degree"]==3
    assert validation["source_sha256"]==sha(HERE/"verify.py")
    assert validation["primary_source_sha256"]==sha(HERE/"construct17.py")
    assert validation["witness_sha256"]==sha(HERE/"cubic17-witness.json")
    assert validation["nonlinear_total"]==58 and validation["collision_points"]==3
    assert all(validation[k] for k in ["full_polynomial_schur_determinant","determinant_one","identity_derivative",
        "retained_rows_independent_of_removed_directions","removed_rows_are_identity_plus_retained_polynomials","quadratic_midpoint_identity"])
    for filename,key in (("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")):
        for item in json.loads((ROOT/filename).read_text())[key]:
            if "KCUB1" in item.get("reviewed_claims",{}):
                assert item["reviewed_claims"]["KCUB1"]==claim_fingerprint(entry)
    for name in ["README.md","PROOF.md","EXPLORATION.md"]:
        for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
            if not target.startswith(("http:","https:","#")):assert (HERE/target.split("#")[0]).is_file(),(name,target)
    for name,source in [("single-slices.json","slice_external.py"),("quadratic-slice-repair.json","repair_slices.py"),("search-checkpoint.json","search.py")]:
        data=json.loads((HERE/name).read_text());assert data["source_sha256"]==sha(HERE/source)
    seen=set()
    def visit(identifier,path):
        assert identifier not in path
        if identifier in seen:return
        for dep in known[identifier]["dependencies"]:
            if dep in known:visit(dep,path|{identifier})
        seen.add(identifier)
    visit("KCUB1",set())
    assert known["OP-GLDC"]["state"]=="open"
    old=json.loads((HERE.parent/"collision-frame-search/validation.json").read_text())
    assert old["normalized_support_interval"]==[10,48]
    receipt=dict(status="pass",claim="KCUB1",dimension=17,degree=3,nonlinear_occurrences=58,
        ordinary_degree_collision_minimum=3,dimension_optimality=False,support_optimality=False,
        normalized_three_variable_support_interval=[10,48],umbrella_status="open",
        schema=True,source_hashes=True,local_links=True,claim_fingerprints=True,dependency_graph_acyclic=True,
        validation_sha256=sha(HERE/"validation.json"),
        navigation_regeneration="Unavailable: pre-existing missing KDSQ6 consumer and METHOD-EC-CENTRE-POLICY source; no bypass.")
    out.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
if __name__=="__main__":main()

