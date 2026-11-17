#!/usr/bin/env python3
"""Scoped integration checks for KCUB6 and its bounded discovery experiments."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
known={e["id"]:e for e in json.loads((ROOT/"MATH_STATUS.json").read_text())["entries"]}
e=known["KCUB6"]
assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
assert e["state"]=="proved" and e["dependencies"]==["BCRDB1","KCUB5"]
assert e["independent_replay"] and not e["formal_verification"] and not e["external_review"]
assert not e["closes_problems"] and not e["narrows_problems"] and e["programme_status"]=="archived"
assert e["artifact_hash"]=="sha256:"+sha(ROOT/e["checker"])
for p in [e["canonical_source"],e["checker"],*e["software_lock"]]:assert (ROOT/p).is_file(),p
v=json.loads((HERE/"fixed-tensor-bound-validation.json").read_text());p=json.loads((HERE/"fixed-tensor-bound-primary.json").read_text())
assert v["status"]=="pass" and v["rank"]==p["rank"]==34
assert v["source_sha256"]==sha(HERE/"verify_fixed_tensor_bound.py")
assert v["primary_diagnostic_sha256"]==sha(HERE/"fixed-tensor-bound-primary.json")
assert p["source_sha256"]==sha(HERE/"compute_fixed_tensor_bound.py") and v["tensor_sha256"]==p["tensor_sha256"]
for folder,script,size,steps in [("power-basis-from89","power_basis_from89.py",178,16980),("power-basis-weighted","power_basis_weighted.py",858,1966)]:
 r=json.loads((HERE/folder/"checkpoint.json").read_text())
 assert r["status"]=="bounded-no-improvement" and r["best"]==89 and r["rank"]==121
 assert r["candidate_forms"]==size and r["steps"]==steps and r["source_sha256"]==sha(HERE/script)
 assert r["seed_candidate_sha256"]==sha(HERE/"power-basis-anneal/candidate.json")
 for name,h in r["source_inputs"].items():assert h==sha(HERE.parent/"de-bondt-2610.24889"/name)
for filename,key in [("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")]:
 for item in json.loads((ROOT/filename).read_text())[key]:
  if "KCUB6" in item.get("reviewed_claims",{}):assert item["reviewed_claims"]["KCUB6"]==claim_fingerprint(e)
for name in ["FIXED_TENSOR_BOUND.md","README.md"]:
 for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
  if not target.startswith(("http:","https:","#")):assert (HERE/target.split("#")[0]).is_file(),target
assert known["OP-GLDC"]["state"]=="open"
receipt=dict(status="pass",claim="KCUB6",fixed_tensor_lower_bound=34,fixed_tensor_upper_bound=89,
 global_dimension_lower_bound=False,both_searches_terminal=True,source_hashes=True,schema=True,reviewed_fingerprints=True,links=True,
 validation_sha256=sha(HERE/"fixed-tensor-bound-validation.json"))
(HERE/"fixed-tensor-bound-integration.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
