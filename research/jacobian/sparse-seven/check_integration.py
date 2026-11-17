"""Scoped metadata and evidence integrity checks for the new result."""
import hashlib
import json
from pathlib import Path
import re
import sys
ROOT=Path("/home/royvanrijn/src/jacobian-research/research")
HERE=ROOT/"jacobian/sparse-seven"
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
ledger=json.loads((ROOT/"MATH_STATUS.json").read_text())
by_id={e["id"]:e for e in ledger["entries"]}
claim=by_id["GLDC2"]
assert claim["state"]=="proved" and by_id["OP-GLDC"]["state"]=="open"
assert claim["artifact_hash"]=="sha256:"+hashlib.sha256((ROOT/claim["checker"]).read_bytes()).hexdigest()
assert all((ROOT/p).is_file() for p in [claim["canonical_source"],claim["checker"],*claim["software_lock"]])
assert claim["narrows_problems"]==["OP-GLDC"] and not claim["closes_problems"]
assert "GLDC2" in by_id["OP-GLDC"]["dependencies"]
for filename,key in (("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")):
    data=json.loads((ROOT/filename).read_text())
    for row in data[key]:
        for identifier in ("GLDC2","OP-GLDC"):
            if identifier in row.get("reviewed_claims",{}):
                assert row["reviewed_claims"][identifier]==claim_fingerprint(by_id[identifier]),(filename,row["id"])
for filename in ("README.md","PROOF.md"):
    for target in re.findall(r"\]\(([^)]+)\)",(HERE/filename).read_text()):
        if not target.startswith(("http:","https:","#")):
            assert (HERE/target.split("#")[0]).exists(),(filename,target)
manifest=json.loads((HERE/"results/manifest.json").read_text())
validation=json.loads((HERE/"validation.json").read_text())
assert validation["status"]=="pass" and validation["algebra_rerun"]
assert validation["manifest_sha256"]==hashlib.sha256((HERE/"results/manifest.json").read_bytes()).hexdigest()
assert validation["coverage_sha256"]==hashlib.sha256((HERE/"coverage-validation.json").read_bytes()).hexdigest()
assert validation["source_sha256"]==hashlib.sha256((HERE/"verify.py").read_bytes()).hexdigest()
supports=json.loads((HERE/"results/supports.json").read_text())["supports"]
assert all((any(e[1]==e[2]==0 for e in rows[1])) !=
           (any(e[1]==e[2]==0 for e in rows[2])) for rows in supports)
receipt=dict(
    status="pass",claim="GLDC2",umbrella_status="open",
    coefficient_replay="2326 rational unit ideals; both backends rerun",
    scoped_metadata=True,scoped_links=True,reviewed_fingerprints=True,
    validation_sha256=hashlib.sha256((HERE/"validation.json").read_bytes()).hexdigest(),
    navigation_regeneration="blocked by existing missing artifacts; not bypassed",
    navigation_errors=[
        "render_status.py: KDSQ6 missing consumer artifacts/generated-results/README.md",
        "research.py render: METHOD-EC-CENTRE-POLICY missing source artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json"],
)
(HERE/"integration.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
print(json.dumps(receipt))

