"""Scoped ledger, source, receipt and link checks for GLDC4/KFRAME1."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ledger=json.loads((ROOT/"MATH_STATUS.json").read_text())
    known={e["id"]:e for e in ledger["entries"]}
    assert known["OP-GLDC"]["state"]=="open"
    for identifier in ("GLDC4","KFRAME1"):
        entry=known[identifier]
        assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
        assert entry["state"]=="proved" and entry["programme_status"]=="archived"
        assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
        assert all((ROOT/p).is_file() for p in
                   [entry["checker"],entry["canonical_source"],*entry["software_lock"]])
        assert all(d in known or d.startswith("external: ") for d in entry["dependencies"])
        assert not entry["closes_problems"]
    assert known["GLDC4"]["narrows_problems"]==["OP-GLDC"]
    assert known["KFRAME1"]["narrows_problems"]==[]
    visited=set()
    def visit(identifier,active):
        assert identifier not in active
        if identifier in visited:return
        for dep in known[identifier]["dependencies"]:
            if dep in known:visit(dep,active|{identifier})
        visited.add(identifier)
    visit("GLDC4",set())
    for filename,key in (("knowledge/lessons.json","lessons"),
                         ("knowledge/work_items.json","actions")):
        for row in json.loads((ROOT/filename).read_text())[key]:
            for identifier in ("GLDC3","GLDC4","KFRAME1","OP-GLDC"):
                if identifier in row.get("reviewed_claims",{}):
                    assert row["reviewed_claims"][identifier]==claim_fingerprint(known[identifier])
    for filename in ("README.md","PROOF.md","NORMALIZED_WITNESS.md"):
        for target in re.findall(r"\]\(([^)]+)\)",(HERE/filename).read_text()):
            if not target.startswith(("http:","https:","#")):
                assert (HERE/target.split("#")[0]).is_file(),(filename,target)
    validation=json.loads((HERE/"validation.json").read_text())
    assert validation["status"]=="pass" and validation["controls_rerun"]
    assert validation["source_sha256"]==sha(HERE/"verify.py")
    assert validation["manifest_sha256"]==sha(HERE/"results/manifest.json")
    assert validation["coverage_sha256"]==sha(HERE/"coverage-validation.json")
    assert validation["controls_sha256"]==sha(HERE/"controls-validation.json")
    assert validation["retained_candidates"]==0 and validation["excluded_nine_term_ideals"]==0
    assert validation["normalized_witness_independent"]["nonlinear_counts"]==[5,18,26]
    assert validation["frame_certificate_independent"]["independent_coefficient_inputs"]
    receipt=dict(
        status="pass",claims=["GLDC4","KFRAME1"],umbrella_status="open",
        normalized_support_interval=[10,49],degree_bound=7,arity=3,
        production_visits=6040767,independent_visits=7228226,
        retained_candidates=0,coefficient_ideals_for_nine_term_proof=0,
        independent_rational_witness_and_family=True,
        scoped_schema=True,scoped_links=True,reviewed_fingerprints=True,dependency_graph_acyclic=True,
        validation_sha256=sha(HERE/"validation.json"),
        navigation_regeneration="Blocked by pre-existing missing artifacts; no renderer bypass.",
        navigation_errors=[
            "render_status.py: KDSQ6 missing consumer artifacts/generated-results/README.md",
            "research.py render: METHOD-EC-CENTRE-POLICY missing source artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json"])
    (HERE/"integration.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

