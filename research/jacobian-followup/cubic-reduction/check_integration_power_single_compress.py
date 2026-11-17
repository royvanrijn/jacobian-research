#!/usr/bin/env python3
"""Scoped authority, certificate and discovery-boundary checks for KCUB7."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 known={e["id"]:e for e in json.loads((ROOT/"MATH_STATUS.json").read_text())["entries"]};e=known["KCUB7"]
 assert REQUIRED_FIELDS<=set(e)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
 assert e["state"]=="proved" and e["dependencies"]==["KCUB5"] and e["programme_status"]=="archived"
 assert e["independent_replay"] and not e["formal_verification"] and not e["external_review"]
 assert not e["closes_problems"] and not e["narrows_problems"]
 assert e["artifact_hash"]=="sha256:"+sha(ROOT/e["checker"])
 for p in [e["canonical_source"],e["checker"],*e["software_lock"]]:assert (ROOT/p).is_file(),p
 out=HERE/"power-single-compress";v=json.loads((out/"validation.json").read_text())
 assert v["status"]=="pass" and v["source_sha256"]==sha(HERE/"verify_power_single_compress.py")
 assert v["arithmetic_sha256"]==sha(HERE/"verify.py") and v["input_sha256"]==sha(HERE/"power-basis-anneal/candidate.json")
 assert v["certification_sha256"]==sha(out/"certification.json") and v["primary_checkpoint_sha256"]==sha(out/"checkpoint.json")
 assert v["construction_validation_sha256"]==sha(HERE/"power89/validation.json")
 assert (v["unit_identities"],v["coordinate_support_cases"],v["proportional_coefficient_pairs"],v["covered_nonempty_cliques"])==(86,15,134,225)
 assert all(v[k] for k in ["all_membership_constraints_reconstructed","all_coordinate_maps_verified_mod_membership","all_coordinate_supports_covered","all_cancellation_pairs_covered","positive_compression_control"])
 for case in v["cases"]:assert case["certificate_sha256"]==sha(out/("cert-%02d.json"%case["case"]))
 line=json.loads((HERE/"power-line-compress/validation.json").read_text());assert line["status"]=="pass" and line["pairs"]==3916
 assert line["source_sha256"]==sha(HERE/"verify_power_line_compress.py")
 pencils=json.loads((HERE/"power-basis-pencils/checkpoint.json").read_text())
 assert pencils["status"]=="bounded-no-improvement" and (pencils["best"],pencils["candidate_forms"],pencils["rank"],pencils["steps"])==(89,256,154,10554)
 assert pencils["source_sha256"]==sha(HERE/"power_basis_pencils.py") and pencils["pool_sha256"]==sha(HERE/"power-basis-pencils/pool.json")
 beam=json.loads((HERE/"linear-beam12/checkpoint.json").read_text());assert beam["status"]=="bounded-no-improvement" and beam["best"]==36 and beam["tested"]==69747
 assert beam["source_sha256"]==sha(HERE/"linear_beam12.py") and beam["input_sha256"]==sha(HERE/"cubic12-integer-witness.json")
 translations=json.loads((HERE/"translation-lines12/checkpoint.json").read_text());assert translations["status"]=="completed-lines-no-improvement" and translations["best"]==36 and len(translations["lines"])==12
 assert translations["source_sha256"]==sha(HERE/"translation_lines12.py")
 for filename,key in [("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")]:
  for item in json.loads((ROOT/filename).read_text())[key]:
   if "KCUB7" in item.get("reviewed_claims",{}):assert item["reviewed_claims"]["KCUB7"]==claim_fingerprint(e)
 for name in ["SINGLE_CUBE_OBSTRUCTION.md","README.md"]:
  for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
   if not target.startswith(("http:","https:","#")):assert (HERE/target.split("#")[0]).is_file(),(name,target)
 assert known["OP-GLDC"]["state"]=="open" and known["KCUB5"]["state"]=="proved"
 receipt=dict(status="pass",claim="KCUB7",unit_identities=86,all_single_new_forms_excluded=True,
   global_dimension_optimality=False,umbrella_status="open",schema=True,source_hashes=True,reviewed_fingerprints=True,links=True,
   validation_sha256=sha(out/"validation.json"))
 (out/"integration.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
if __name__=="__main__":main()
