#!/usr/bin/env python3
"""Verify coverage, independently compile all coefficient ideals, audit exact replays.

--replay-algebra also reruns both rational Groebner backends. The default
verifies the retained result records and independent coefficient inputs, and
re-enumerates the entire finite support cover; hashes alone do not prove it.
"""
import argparse
from collections import defaultdict
from itertools import permutations, product
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import search
import sympy as sp

HERE=Path(__file__).resolve().parent
RESULTS=HERE/"results"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def independent_equations(support):
    n=support.nonlinear_size
    coefficient_index={atom:i+1 for i,atom in enumerate(support.atoms)}
    linear=((1,0,0),(0,1,0),(0,0,1))
    # Derivative entries as (source exponent, coefficient-variable exponent,
    # integer multiplier); the zeroth variable is the saturation variable rho.
    jac=[[[] for j in range(3)] for i in range(3)]
    for i in range(3):
        for exponent in (linear[i],*support.rows[i]):
            powers=[0]*(n+1)
            if exponent!=linear[i]:
                powers[coefficient_index[i,exponent]]=1
            for j in range(3):
                if exponent[j]:
                    derivative=list(exponent)
                    derivative[j]-=1
                    jac[i][j].append((tuple(derivative),tuple(powers),exponent[j]))
    determinant=defaultdict(lambda:defaultdict(int))
    for perm in permutations(range(3)):
        inversions=sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        sign=(-1)**inversions
        for triple in product(*(jac[i][perm[i]] for i in range(3))):
            source=tuple(sum(term[0][j] for term in triple) for j in range(3))
            coeff=tuple(sum(term[1][j] for term in triple) for j in range(n+1))
            value=sign
            for term in triple:
                value*=term[2]
            determinant[source][coeff]+=value
    constant={e:c for e,c in determinant.pop((0,0,0)).items() if c}
    assert constant=={(0,)*(n+1):1}
    expected=[{e:c for e,c in coefficients.items() if c}
              for coefficients in determinant.values()]
    expected=[p for p in expected if p]
    for i,row in enumerate(support.rows):
        axis=[e for e in row if e[1:]==(0,0)]
        if i==0 or axis:
            equation={}
            if i==0:
                equation[(0,)*(n+1)]=1
            for e in axis:
                monomial=[0]*(n+1)
                monomial[coefficient_index[i,e]]=1
                equation[tuple(monomial)]=1
            expected.append(equation)
    expected.append({(1,)*(n+1):1,(0,)*(n+1):-1})
    return expected

def canonical(polys):
    return {tuple(sorted(poly.items())) for poly in polys}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay-algebra",action="store_true")
    a=parser.parse_args()
    started=time.monotonic()
    manifest_path=RESULTS/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    assert manifest["status"]=="complete" and manifest["all_unit"]
    assert (manifest["degree"],manifest["bound"],manifest["shard_count"])==(7,7,90)
    assert manifest["source_sha256"]==sha(HERE/"search.py")
    assert manifest["runner_sha256"]==sha(HERE/"run_all.py")
    assert manifest["legacy_sha256"]==sha(search.LEGACY/"jcsearch/low_degree_census.py")
    expected_seeds={(p,q,r) for p in range(2,8) for q in range(2,8) for r in range(q+1,8)}
    assert len(manifest["shards"])==90
    assert {tuple(row["seed"]) for row in manifest["shards"]}==expected_seeds
    for row in manifest["shards"]:
        path=RESULTS/row["path"]
        assert sha(path)==row["sha256"]
        shard=json.loads(path.read_text())
        assert shard["source_sha256"]==manifest["source_sha256"]
        assert shard["legacy_sha256"]==manifest["legacy_sha256"]
        assert shard["parameters"]==dict(degree=7,bound=7,seed=row["seed"])
        assert shard["status"]=="complete" and shard["pending"]==0 and not shard["stack"]
    subprocess.run([sys.executable,str(HERE/"verify_coverage.py")],check=True)
    assert sha(RESULTS/"supports.json")==manifest["supports_sha256"]
    supports_record=json.loads((RESULTS/"supports.json").read_text())
    supports=[search.read_support(s) for s in supports_record["supports"]]
    oriented={search.read_support(s) for row in manifest["shards"]
              for s in json.loads((RESULTS/row["path"]).read_text())["supports"]}
    assert supports==sorted({s.canonical_under_collision_stabilizer() for s in oriented})
    assert len(supports)==manifest["orbit_count"]==2326
    assert manifest["orbit_counts_by_size"]=={"6":450,"7":1876}
    ids=[s.identifier for s in supports]
    assert len(ids)==len(set(ids))
    results=[]
    replay=[]
    for row in manifest["algebra_batches"]:
        path=RESULTS/row["path"]
        assert sha(path)==row["sha256"]
        batch=json.loads(path.read_text())
        assert batch["support_ids"]==[r["support_id"] for r in batch["results"]]
        assert batch["support_ids"]==[r["support_id"] for r in batch["sympy_results"]]
        results.extend(batch["results"])
        replay.extend(batch["sympy_results"])
    assert [r["support_id"] for r in results]==ids
    assert [r["support_id"] for r in replay]==ids
    assert all(r["field"]=="QQ" and r["unit_ideal"] and r["basis_length"]==1
               for r in results+replay)
    for i,support in enumerate(supports):
        variables,equations=search.coefficient_system(support)
        actual=[dict(sp.Poly(eq,*variables).terms()) for eq in equations]
        assert canonical(actual)==canonical(independent_equations(support)),support.identifier
        if (i+1)%500==0:
            print(json.dumps(dict(stage="independent-inputs",done=i+1,total=len(supports))),flush=True)
    if a.replay_algebra:
        for offset in range(0,len(supports),100):
            batch=supports[offset:offset+100]
            assert all(r["unit_ideal"] for r in search.singular_batch_status(batch))
            assert all(search.groebner_status(s)["unit_ideal"] for s in batch)
    receipt=dict(
        status="pass",claim="GLDC2",source_sha256=sha(Path(__file__)),
        manifest_sha256=sha(manifest_path),coverage_sha256=sha(HERE/"coverage-validation.json"),
        verified_coefficient_ideals=len(supports),
        seven_term_orbits=1876,six_term_regression_orbits=450,
        independent_coverage=True,independent_coefficient_compilation=True,
        rational_backends=["Singular","SymPy"],
        algebra_rerun=a.replay_algebra,
        scope="No normalized degree-at-most-seven Keller collision with at most seven nonlinear monomial occurrences.",
        elapsed_seconds=time.monotonic()-started,
        versions=dict(python=sys.version.split()[0],sympy=sp.__version__,
                      singular=subprocess.check_output(["Singular","--version"],text=True).splitlines()[0]),
    )
    search.write_json(HERE/"validation.json",receipt)
    print(json.dumps(receipt),flush=True)
if __name__=="__main__":
    main()

