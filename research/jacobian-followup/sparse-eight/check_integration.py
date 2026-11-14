"""Scoped integrity and navigation checks for GLDC3 and KPAIR1."""
import hashlib
import json
from pathlib import Path
import re
import sys
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ledger=json.loads((ROOT/"MATH_STATUS.json").read_text())
    known={e["id"]:e for e in ledger["entries"]}
    assert known["OP-GLDC"]["state"]=="open"
    for identifier in ("GLDC3","KPAIR1"):
        entry=known[identifier]
        assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
        assert entry["state"]=="proved" and entry["programme_status"]=="archived"
        assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
        assert all((ROOT/p).is_file() for p in
                   [entry["checker"],entry["canonical_source"],*entry["software_lock"]])
        assert all(dep in known or dep.startswith("external: ") for dep in entry["dependencies"])
        assert not entry["closes_problems"] and entry["narrows_problems"]==["OP-GLDC"]
        assert identifier in known["OP-GLDC"]["dependencies"]
    for filename,key in (("knowledge/lessons.json","lessons"),
                         ("knowledge/work_items.json","actions")):
        for row in json.loads((ROOT/filename).read_text())[key]:
            for identifier in ("GLDC2","GLDC3","KPAIR1","OP-GLDC"):
                if identifier in row.get("reviewed_claims",{}):
                    assert row["reviewed_claims"][identifier]==claim_fingerprint(known[identifier])
    for filename in ("README.md","PROOF.md","COORDINATE_GATES.md"):
        for target in re.findall(r"\]\(([^)]+)\)",(HERE/filename).read_text()):
            if not target.startswith(("https:","http:","#")):
                assert (HERE/target.split("#")[0]).is_file(),(filename,target)
    validation=json.loads((HERE/"validation.json").read_text())
    assert validation["status"]=="pass" and validation["algebra_rerun"]
    assert validation["source_sha256"]==sha(HERE/"verify.py")
    assert validation["manifest_sha256"]==sha(HERE/"results/manifest.json")
    assert validation["coverage_sha256"]==sha(HERE/"coverage-validation.json")
    assert validation["independent_input_compiler_sha256"]==sha(HERE.parent/"sparse-seven/verify.py")
    assert validation["structural_sha256"]==sha(HERE/"structural-validation.json")
    gates=json.loads((HERE/"coordinate-gates-validation.json").read_text())
    assert gates["status"]=="pass"
    assert gates["source_sha256"]==sha(HERE/"verify_coordinate_gates.py")
    assert gates["scalar_pairing_exclusion_max_three_variable_degree"]==10
    assert gates["four_variable_degrees"]==[4,6,7,1]
    receipt=dict(
        status="pass",claims=["GLDC3","KPAIR1"],umbrella_status="open",
        lower_bound=9,degree_bound=7,arity=3,
        eight_term_orbits=17482,exact_rational_backends=2,
        scoped_schema=True,scoped_links=True,reviewed_fingerprints=True,
        validation_sha256=sha(HERE/"validation.json"),
        coordinate_gates_validation_sha256=sha(HERE/"coordinate-gates-validation.json"),
        navigation_regeneration="Blocked by existing unrelated missing artifacts; renderer failures were not bypassed.",
        navigation_errors=[
            "render_status.py: KDSQ6 missing consumer artifacts/generated-results/README.md",
            "research.py render: METHOD-EC-CENTRE-POLICY missing source artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json"],
    )
    (HERE/"integration.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":
    main()

