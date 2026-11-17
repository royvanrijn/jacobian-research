"""Scoped metadata, source and certificate checks for the thirteen-variable cubic."""
from pathlib import Path
import hashlib,itertools,json,re,sys
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
from research import claim_fingerprint
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=HERE/"integration13.json"
    if not out.exists():out.write_text('{"status":"pending"}\n')
    known={e["id"]:e for e in json.loads((ROOT/"MATH_STATUS.json").read_text())["entries"]}
    entry=known["KCUB2"]
    assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
    assert entry["state"]=="proved" and entry["programme_status"]=="archived"
    assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
    assert entry["dependencies"]==["KFRAME1","KCUB1"]
    assert not entry["closes_problems"] and not entry["narrows_problems"] and not entry["formal_verification"] and not entry["external_review"]
    assert entry["independent_replay"]
    assert all((ROOT/p).is_file() for p in [entry["canonical_source"],entry["checker"],*entry["software_lock"]])
    d=json.loads((HERE/"sparse-validation13.json").read_text())
    assert d["status"]=="pass" and d["dimension"]==13 and d["degree"]==3 and d["nonlinear_total"]==42
    assert d["source_sha256"]==sha(HERE/"verify_sparse13.py")
    assert d["arithmetic_sha256"]==sha(HERE/"verify.py") and d["replay_sha256"]==sha(HERE/"verify13.py")
    assert d["primary_source_sha256"]==sha(HERE/"construct_sparse13.py")
    assert d["witness_sha256"]==sha(HERE/"cubic13-sparse-witness.json")
    assert d["base_sha256"]==sha(HERE/"cubic13-witness.json") and d["discovery_sha256"]==sha(HERE/"reuse-witness.json")
    assert all(d[k] for k in ["determinant_one","identity_derivative","all_variables_symbolic",
        "full_polynomial_schur_determinant","base_determinant_verified","entire_transformations_replayed"])
    assert d["collision_points"]==3 and d["elementary_base_moves"]=={"carrier":10,"product":13,"linear":4}
    for filename,key in [("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")]:
        for item in json.loads((ROOT/filename).read_text())[key]:
            if "KCUB2" in item.get("reviewed_claims",{}):assert item["reviewed_claims"]["KCUB2"]==claim_fingerprint(entry)
    for name in ["README.md","PROOF13.md","PROOF.md","EXPLORATION.md"]:
        for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
            if not target.startswith(("http:","https:","#")):assert (HERE/target.split("#")[0]).is_file(),(name,target)
    manifest=json.loads((HERE/"orders-manifest.json").read_text())
    assert manifest["status"]=="complete" and len(manifest["cases"])==18
    assert manifest["source_sha256"]==sha(HERE/"search_reuse_orders.py")
    assert {row["name"] for row in manifest["cases"]}=={p+"-"+o for p,o in itertools.product(["012","021","102","120","201","210"],["grlex","grevlex","aux-first"])}
    for row in manifest["cases"]:
        f=HERE/row["checkpoint"];q=json.loads(f.read_text())
        assert row["sha256"]==sha(f) and q["status"]=="bounded-search-exhausted"
        assert q["source_sha256"]==manifest["source_sha256"]
        assert q["parameters"]["beam"]==800 and q["parameters"]["max_dimension"]==12
    for name,source in [("quadratic-slice-repair13.json","repair_cubic13.py"),("linear-descent13.json","optimize13.py")]:
        q=json.loads((HERE/name).read_text());assert q["source_sha256"]==sha(HERE/source)
    seen=set()
    def visit(identifier,path):
        assert identifier not in path
        if identifier in seen:return
        for dep in known[identifier]["dependencies"]:
            if dep in known:visit(dep,path|{identifier})
        seen.add(identifier)
    visit("KCUB2",set())
    assert known["OP-GLDC"]["state"]=="open"
    interval=json.loads((HERE.parent/"collision-frame-search/validation.json").read_text())["normalized_support_interval"]
    assert interval==[10,48]
    receipt=dict(status="pass",claim="KCUB2",dimension=13,degree=3,nonlinear_occurrences=42,
       dimension_optimality=False,support_optimality=False,umbrella_status="open",
       normalized_three_variable_support_interval=interval,
       schema=True,source_hashes=True,local_links=True,reviewed_fingerprints=True,dependency_graph_acyclic=True,
       validation_sha256=sha(HERE/"sparse-validation13.json"),bounded_order_runs=18,
       navigation_regeneration="Global regeneration fails at unrelated metadata and missing artifacts; no bypass.",
       navigation_errors=["render_status.py: HC4MYGS1 invalid kind", "research.py render: METHOD-EC-CENTRE-POLICY missing source artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json"])
    out.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps(receipt))
if __name__=="__main__":main()

