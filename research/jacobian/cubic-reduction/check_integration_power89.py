"""Scoped authority, identity and artifact checks for the89-variable cubic-linear witness."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=HERE/"power89/integration.json"
    if not out.exists():out.write_text('{"status":"pending"}\n')
    known={e["id"]:e for e in json.loads((ROOT/"MATH_STATUS.json").read_text())["entries"]};entry=known["KCUB5"]
    assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
    assert entry["state"]=="proved" and entry["programme_status"]=="archived" and entry["dependencies"]==["BCRDB1"]
    assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
    assert all((ROOT/p).is_file() for p in [entry["canonical_source"],entry["checker"],*entry["software_lock"]])
    assert entry["independent_replay"] and not entry["formal_verification"] and not entry["external_review"]
    assert not entry["closes_problems"] and not entry["narrows_problems"]
    v=json.loads((HERE/"power89/validation.json").read_text())
    assert v["status"]=="pass" and v["dimension"]==89 and v["coefficient_ring"]=="Z" and v["degree"]==3
    assert v["source_sha256"]==sha(HERE/"verify_power89.py") and v["witness_sha256"]==sha(HERE/"power89/witness.json")
    assert v["primary_source_sha256"]==sha(HERE/"build_power89.py")
    assert v["arithmetic_sha256"]==sha(HERE/"verify.py") and v["base_determinant_checker_sha256"]==sha(HERE/"verify_sparse12.py")
    assert v["candidate_sha256"]==sha(HERE/"power-basis-anneal/candidate.json")
    assert all(v[k] for k in ["all18_cube_identities","all18_projection_identities","integer_matrix_factorization",
        "right_inverse_verified","base12_full_polynomial_determinant","determinant_via_companion_and_sylvester",
        "source_and_target_shear_retained","full_rational_collision_verified","distinct_points_verified","identity_linear_part"])
    assert not v["expanded89_variable_determinant"]
    run=json.loads((HERE/"power-basis-anneal/checkpoint.json").read_text())
    assert run["status"]=="candidate" and run["best"]==89 and run["source_sha256"]==sha(HERE/"power_basis_anneal.py")
    for name,digest in run["source_inputs"].items():assert digest==sha(HERE.parent/"de-bondt-2610.24889"/name)
    timeout=json.loads((HERE/"power-search/checkpoint.json").read_text())
    assert timeout["status"]=="unknown" and timeout["reason"]=="timeout" and timeout["source_sha256"]==sha(HERE/"power_search.py")
    plateau=json.loads((HERE/"power-basis/checkpoint.json").read_text())
    assert plateau["status"]=="bounded-no-improvement" and plateau["best"]==90
    for filename,key in [("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")]:
        for item in json.loads((ROOT/filename).read_text())[key]:
            if "KCUB5" in item.get("reviewed_claims",{}):assert item["reviewed_claims"]["KCUB5"]==claim_fingerprint(entry)
    for name in ["POWER89.md","README.md"]:
        for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
            if not target.startswith(("http:","https:","#")):assert (HERE/target.split("#")[0]).is_file(),(name,target)
    visited=set()
    def visit(key,path):
        assert key not in path
        if key in visited:return
        for dep in known[key]["dependencies"]:
            if dep in known:visit(dep,path|{key})
        visited.add(key)
    visit("KCUB5",set())
    assert known["OP-GLDC"]["state"]=="open"
    receipt=dict(status="pass",claim="KCUB5",dimension=89,previous_audited_dimension=90,form="Y+(AY)^{circ3}",
        coefficient_ring="Z",explicit_rational_collision=True,determinant_via_exact_factorization=True,
        dimension_optimality=False,umbrella_status="open",schema=True,source_hashes=True,links=True,
        reviewed_fingerprints=True,dependency_graph_acyclic=True,
        validation_sha256=sha(HERE/"power89/validation.json"))
    out.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
if __name__=="__main__":main()

