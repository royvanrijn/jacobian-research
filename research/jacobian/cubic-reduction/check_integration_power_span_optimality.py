#!/usr/bin/env python3
"""Scoped authority, classification and finite two-exchange checks for KCUB8."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 known={e["id"]:e for e in load(ROOT/"MATH_STATUS.json")["entries"]};entry=known["KCUB8"]
 assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
 assert entry["state"]=="proved" and entry["dependencies"]==["KCUB5"] and entry["programme_status"]=="archived"
 assert entry["independent_replay"] and not entry["formal_verification"] and not entry["external_review"]
 assert not entry["closes_problems"] and not entry["narrows_problems"]
 assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
 for p in [entry["canonical_source"],entry["checker"],*entry["software_lock"]]:assert (ROOT/p).is_file(),p
 out=HERE/"power-span-split";v=load(out/"validation.json")
 assert v["status"]=="pass" and v["source_sha256"]==sha(HERE/"verify_power_span_optimality.py")
 assert v["arithmetic_sha256"]==sha(HERE/"verify.py") and v["helper_sha256"]==sha(HERE/"verify_power_single_compress.py")
 assert v["input_sha256"]==sha(HERE/"power-basis-anneal/candidate.json")
 assert v["primary_checkpoint_sha256"]==sha(out/"checkpoint.json") and v["classification_sha256"]==sha(out/"classification.json")
 assert (v["tree_nodes"],v["factor_identities"],v["isolated_projective_points"],v["projective_lines"])==(408,236,65,8)
 assert (v["cube_span_dimension"],v["line_cube_span_dimension"],v["annihilator_constraint_rank"],v["special_target_rank"],v["zero_locus_special_span_bound"])==(89,24,87,7,5)
 assert all(v[k] for k in ["all_factor_trees_checked","all_cube_supports_covered","all_line_cubic_spans_checked","annihilator_restrictions_checked"])
 for case in v["cases"]:assert case["tree_sha256"]==sha(out/("case-%02d.json"%case["case"]))
 q=load(HERE/"power-two-exchange/validation.json")
 assert q["status"]=="pass" and q["source_sha256"]==sha(HERE/"verify_power_two_exchange.py")
 assert q["fixed_span_validation_sha256"]==sha(out/"validation.json")
 assert q["primary_sha256"]==sha(HERE/"power-two-exchange/checkpoint.json") and q["pool_sha256"]==sha(HERE/"power-two-exchange/pool.json")
 assert (q["candidate_forms"],q["cube_pool_rank"],q["inside_saved_span"],q["nonzero_residual_directions"],q["pairs_tested"],q["max_canceled"])==(930,154,151,445,1774,2)
 primary=load(HERE/"power-two-exchange/checkpoint.json");assert primary["source_sha256"]==sha(HERE/"power_two_exchange.py")
 for name,digest in primary["source_inputs"].items():assert digest==sha(HERE.parent/"de-bondt-2610.24889"/name)
 for filename,key in [("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")]:
  for item in load(ROOT/filename)[key]:
   if "KCUB8" in item.get("reviewed_claims",{}):assert item["reviewed_claims"]["KCUB8"]==claim_fingerprint(entry)
 for name in ["SPAN_OPTIMALITY.md","README.md","SINGLE_CUBE_OBSTRUCTION.md"]:
  for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
   if not target.startswith(("http:","https:","#")):assert (HERE/target.split("#")[0]).is_file(),(name,target)
 assert known["OP-GLDC"]["state"]=="open"
 receipt=dict(status="pass",claim="KCUB8",fixed_scalar_span_minimum=89,global_tensor_rank_minimum=False,global_Keller_dimension_minimum=False,
  classification=True,annihilator_premises=True,finite_two_exchange=True,schema=True,source_hashes=True,reviewed_fingerprints=True,links=True,
  validation_sha256=sha(out/"validation.json"),pair_validation_sha256=sha(HERE/"power-two-exchange/validation.json"),umbrella_status="open")
 (out/"integration.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
if __name__=="__main__":main()
