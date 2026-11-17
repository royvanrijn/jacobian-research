#!/usr/bin/env python3
"""Verify complete eight-term coverage and every exact coefficient input/result."""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time
import sympy as sp

HERE=Path(__file__).resolve().parent
SEVEN=HERE.parent/"sparse-seven"
sys.path.insert(0,str(SEVEN))
import search
spec=importlib.util.spec_from_file_location("sealed_seven_verifier",SEVEN/"verify.py")
previous=importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)
RESULTS=HERE/"results"
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay-algebra",action="store_true")
    parser.add_argument("--replay-coverage",action="store_true")
    a=parser.parse_args()
    started=time.monotonic()
    path=RESULTS/"manifest.json"
    manifest=json.loads(path.read_text())
    assert manifest["status"]=="complete" and manifest["all_unit"]
    assert (manifest["degree"],manifest["bound"],manifest["shard_count"])==(7,8,90)
    assert manifest["source_sha256"]==sha(SEVEN/"search.py")
    assert manifest["runner_sha256"]==sha(SEVEN/"run_all.py")
    assert manifest["legacy_sha256"]==sha(search.LEGACY/"jcsearch/low_degree_census.py")
    seeds={(p,q,r) for p in range(2,8) for q in range(2,8) for r in range(q+1,8)}
    assert len(manifest["shards"])==90
    assert {tuple(row["seed"]) for row in manifest["shards"]}==seeds
    oriented=set()
    for row in manifest["shards"]:
        shard_path=RESULTS/row["path"]
        assert sha(shard_path)==row["sha256"]
        state=json.loads(shard_path.read_text())
        assert state["parameters"]==dict(degree=7,bound=8,seed=row["seed"])
        assert state["source_sha256"]==manifest["source_sha256"]
        assert state["legacy_sha256"]==manifest["legacy_sha256"]
        assert state["status"]=="complete" and state["pending"]==0 and not state["stack"]
        oriented.update(search.read_support(s) for s in state["supports"])
    if a.replay_coverage:
        subprocess.run([sys.executable,str(HERE/"verify_coverage.py")],check=True)
    coverage=json.loads((HERE/"coverage-validation.json").read_text())
    assert coverage["status"]=="pass" and (coverage["degree"],coverage["bound"])==(7,8)
    assert coverage["manifest_sha256"]==sha(path)
    assert coverage["source_sha256"]==sha(HERE/"verify_coverage.py")
    assert coverage["balanced"]==len(oriented)==23183
    assert coverage["counts_by_size"]=={"6":450,"7":1876,"8":20857}
    assert sha(RESULTS/"supports.json")==manifest["supports_sha256"]
    supports_record=json.loads((RESULTS/"supports.json").read_text())
    supports=[search.read_support(s) for s in supports_record["supports"]]
    assert supports==sorted({s.canonical_under_collision_stabilizer() for s in oriented})
    assert len(supports)==manifest["orbit_count"]==19808
    expected={"6":450,"7":1876,"8":17482}
    assert manifest["orbit_counts_by_size"]==expected
    assert {str(k):v for k,v in Counter(s.nonlinear_size for s in supports).items()}==expected
    ids=[s.identifier for s in supports]
    assert len(ids)==len(set(ids))
    singular=[];sympy=[]
    for row in manifest["algebra_batches"]:
        batch_path=RESULTS/row["path"]
        assert sha(batch_path)==row["sha256"]
        data=json.loads(batch_path.read_text())
        assert data["support_ids"]==[r["support_id"] for r in data["results"]]
        assert data["support_ids"]==[r["support_id"] for r in data["sympy_results"]]
        singular.extend(data["results"]);sympy.extend(data["sympy_results"])
    assert [r["support_id"] for r in singular]==ids==[r["support_id"] for r in sympy]
    assert all(r["field"]=="QQ" and r["unit_ideal"] and r["basis_length"]==1
               for r in singular+sympy)
    for i,support in enumerate(supports):
        variables,equations=search.coefficient_system(support)
        actual=[dict(sp.Poly(eq,*variables).terms()) for eq in equations]
        assert previous.canonical(actual)==previous.canonical(previous.independent_equations(support))
        if (i+1)%2000==0:
            print(json.dumps(dict(stage="independent-inputs",done=i+1,total=len(supports))),flush=True)
    if a.replay_algebra:
        for offset in range(0,len(supports),100):
            batch=supports[offset:offset+100]
            assert all(r["unit_ideal"] for r in search.singular_batch_status(batch))
            assert all(search.groebner_status(s)["unit_ideal"] for s in batch)
            if (offset+100)%2000==0:
                print(json.dumps(dict(stage="rational-replay",done=offset+100,total=len(supports))),flush=True)
    subprocess.run([sys.executable,str(HERE/"verify_structure.py")],check=True)
    structural=json.loads((HERE/"structural-validation.json").read_text())
    assert structural["supports_sha256"]==manifest["supports_sha256"]
    receipt=dict(
        status="pass",claim="GLDC3",source_sha256=sha(Path(__file__)),
        manifest_sha256=sha(path),coverage_sha256=sha(HERE/"coverage-validation.json"),
        independent_input_compiler_sha256=sha(SEVEN/"verify.py"),
        structural_sha256=sha(HERE/"structural-validation.json"),
        verified_coefficient_ideals=len(supports),eight_term_orbits=17482,
        smaller_regression_orbits=2326,independent_coverage=True,
        independent_coefficient_compilation=True,
        algebra_rerun=a.replay_algebra,coverage_rerun=a.replay_coverage,
        structural_replay=True,elapsed_seconds=time.monotonic()-started,
        scope="No normalized three-variable Keller collision in raw degree <=7 with <=8 nonlinear occurrences.",
        versions=dict(python=sys.version.split()[0],sympy=sp.__version__,
                      singular=subprocess.check_output(["Singular","--version"],text=True).splitlines()[0]),
    )
    search.write_json(HERE/"validation.json",receipt)
    print(json.dumps(receipt),flush=True)
if __name__=="__main__":
    main()

