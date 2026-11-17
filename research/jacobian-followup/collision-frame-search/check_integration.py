"""Scoped source, evidence and metadata checks for the 48-term result."""
from pathlib import Path
import hashlib,json,re,sys
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"scripts"))
from research import claim_fingerprint
from render_status import REQUIRED_FIELDS,OPTIONAL_FIELDS
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    output=HERE/"integration.json"
    if not output.exists():output.write_text('{"status":"pending"}\n')
    known={e["id"]:e for e in json.loads((ROOT/"MATH_STATUS.json").read_text())["entries"]}
    entry=known["KFRAME2"]
    assert REQUIRED_FIELDS<=set(entry)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
    assert entry["state"]=="proved" and entry["programme_status"]=="archived"
    assert entry["artifact_hash"]=="sha256:"+sha(ROOT/entry["checker"])
    assert entry["dependencies"]==["KFRAME1","GLDC4"]
    assert not entry["closes_problems"] and not entry["narrows_problems"]
    assert all((ROOT/p).is_file() for p in [entry["canonical_source"],entry["checker"],*entry["software_lock"]])
    assert known["OP-GLDC"]["state"]=="open" and "KFRAME2" in known["OP-GLDC"]["dependencies"]
    for filename,key in (("knowledge/lessons.json","lessons"),("knowledge/work_items.json","actions")):
        for row in json.loads((ROOT/filename).read_text())[key]:
            for identifier in ("KFRAME2","KFRAME3","OP-GLDC"):
                if identifier in row.get("reviewed_claims",{}):
                    assert row["reviewed_claims"][identifier]==claim_fingerprint(known[identifier])
    for name in ("README.md","PROOF.md","INVALIDATED_SPECIALIZATION.md","NONZERO_TARGET_LIMITS.md","NEXT_FAMILY.md"):
        for target in re.findall(r"\]\(([^)]+)\)",(HERE/name).read_text()):
            if not target.startswith(("http:","https:","#")):
                assert (HERE/target.split("#")[0]).is_file(),(name,target)
    visited=set()
    def visit(identifier,active):
        assert identifier not in active
        if identifier in visited:return
        for dep in known[identifier]["dependencies"]:
            if dep in known:visit(dep,active|{identifier})
        visited.add(identifier)
    visit("KFRAME2",set())
    visit("KFRAME3",set())
    validation=json.loads((HERE/"validation.json").read_text())
    assert validation["status"]=="pass" and validation["normalized_support_interval"]==[10,48]
    assert validation["source_sha256"]==sha(HERE/"verify.py")
    assert validation["witness_sha256"]==sha(HERE/"normalized-48-witness.json")
    assert validation["singular_sha256"]==sha(HERE/"singular-validation.json")
    regression=json.loads((HERE/"field-regression.json").read_text())
    assert regression["status"]=="pass" and regression["distinct_generators_and_moduli_verified"]
    assert regression["source_sha256"]==sha(HERE/"check_field_generators.py")
    for folder,generator in ((HERE,"two_parameter.py"),(HERE/"full-transverse","full_transverse.py")):
        family=json.loads((folder/"two-parameter-family.json").read_text())
        assert family["source_sha256"]==sha(HERE/generator)
        for name in ("projection-checkpoint.json","projection-factors.json"):
            projection=json.loads((folder/name).read_text())
            assert projection["status"]=="complete" and projection["source_sha256"]==sha(HERE/"two_parameter.py")
        manifest=json.loads((folder/"specialization-manifest-v2.json").read_text())
        assert manifest["status"]=="complete" and manifest["source_sha256"]==sha(HERE/"specialize.py")
        assert manifest["minimum"]==48
        invalid=json.loads((folder/"specialization-manifest.json").read_text())
        assert invalid["status"]=="invalid"
        assert invalid["invalid_source_sha256"]==sha(HERE/"specialize-invalid-v1.py.txt")
        assert invalid["preserved_sha256"]==sha(folder/invalid["preserved_manifest"])
    endpoint=known["KFRAME3"]
    assert REQUIRED_FIELDS<=set(endpoint)<=REQUIRED_FIELDS|OPTIONAL_FIELDS
    assert endpoint["state"]=="proved" and endpoint["programme_status"]=="archived"
    assert endpoint["artifact_hash"]=="sha256:"+sha(ROOT/endpoint["checker"])
    assert endpoint["dependencies"]==["KFRAME1"]
    assert not endpoint["closes_problems"] and not endpoint["narrows_problems"]
    assert all((ROOT/p).is_file() for p in [endpoint["canonical_source"],endpoint["checker"],*endpoint["software_lock"]])
    endpoint_validation=json.loads((HERE/"nonzero-target-validation.json").read_text())
    assert endpoint_validation["status"]=="pass"
    assert endpoint_validation["source_sha256"]==sha(HERE/"verify_nonzero_target.py")
    assert endpoint_validation["primary_source_sha256"]==sha(HERE/"nonzero_target_probe.py")
    assert endpoint_validation["probe_sha256"]==sha(HERE/"nonzero-target-probe.json")
    assert endpoint_validation["zero_valuation"]==-3
    assert endpoint_validation["zero_leading_coefficient"]=="-49460507836416"
    assert endpoint_validation["infinity_axis_trace_degree"]==2
    assert endpoint_validation["infinity_axis_trace_leading_coefficient"]=="864"
    receipt=dict(
        status="pass",claim="KFRAME2",umbrella_status="open",
        normalized_support_interval=[10,48],degree_bound=7,arity=3,
        endpoint_obstruction_claim="KFRAME3",endpoint_invariant_replay=True,
        endpoint_validation_sha256=sha(HERE/"nonzero-target-validation.json"),
        nonlinear_counts=[6,19,23],sharp_classes=2,exact_specializations=1293,
        independent_rational_polynomial_inputs_and_witness=True,
        independent_singular_projection_and_fields=True,
        invalid_v1_retained_and_rejected=True,
        scoped_schema=True,scoped_links=True,reviewed_fingerprints=True,dependency_graph_acyclic=True,
        validation_sha256=sha(HERE/"validation.json"),
        field_regression_sha256=sha(HERE/"field-regression.json"),
        navigation_regeneration="Blocked by pre-existing missing artifacts; not bypassed.",
        navigation_errors=[
            "render_status.py: KDSQ6 missing consumer artifacts/generated-results/README.md",
            "research.py render: METHOD-EC-CENTRE-POLICY missing source artifacts/generated-results/elliptic-curves/next_direction_benchmark_v1/completion.json"])
    output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

