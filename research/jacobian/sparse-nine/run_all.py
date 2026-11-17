#!/usr/bin/env python3
"""Run/resume all 90 coordinate-gated nine-term shards and solve residual ideals."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
import gated_search as search

HERE=Path(__file__).resolve().parent
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--degree",type=int,default=7)
    p.add_argument("--bound",type=int,default=9)
    p.add_argument("--seconds",type=int,default=120)
    p.add_argument("--nodes",type=int,default=500000)
    p.add_argument("--output",type=Path,default=HERE/"results")
    a=p.parse_args()
    compiler=search.old.SupportEnumerator(a.degree,a.bound)
    seeds=[(p,q,r) for p in range(2,a.degree+1)
           for q in range(2,a.degree+1) for r in range(q+1,a.degree+1)]
    shards=[];supports=set();gate_visits=Counter()
    started=time.monotonic()
    for seed in seeds:
        path=a.output/("p%d-q%d-r%d.json"%seed)
        shard=search.Shard(a.degree,a.bound,*seed,compiler)
        state=shard.run_gated(path,a.seconds,a.nodes)
        shards.append(dict(seed=list(seed),path=path.name,sha256=sha(path),
                           status=state["status"],visited=state["visited"]))
        gate_visits.update(state["gate_visits"])
        supports.update(search.old.read_support(rows) for rows in state["supports"])
    complete=all(s["status"]=="complete" for s in shards)
    manifest=dict(
        schema="coordinate-gated-keller-cover.v1",
        status="enumeration-complete" if complete else "incomplete",
        enumeration_status="complete" if complete else "incomplete",
        algebra_status="not-run",degree=a.degree,bound=a.bound,
        source_sha256=sha(HERE/"gated_search.py"),runner_sha256=sha(Path(__file__)),
        predecessor_sha256=sha(search.SEVEN/"search.py"),
        legacy_sha256=sha(search.old.LEGACY/"jcsearch/low_degree_census.py"),
        gate_dependency="KPAIR1; degree <=10; coordinate and hyperplane proofs",
        shards=shards,shard_count=len(shards),
        visited=sum(s["visited"] for s in shards),gate_visits=dict(gate_visits),
        oriented_retained=len(supports),
        scope="All hypothetical normalized Keller collisions in raw degree <=7 and support <=9, after the proved necessary coordinate gates; not a count of every raw balanced support.",
    )
    search.old.write_json(a.output/"manifest.json",manifest)
    if not complete:
        print(json.dumps(dict(status="incomplete",pending_shards=[s["seed"] for s in shards if s["status"]!="complete"])),flush=True)
        return
    reps=sorted({s.canonical_under_collision_stabilizer() for s in supports})
    support_record=dict(
        status="complete",degree=a.degree,bound=a.bound,oriented_retained=len(supports),
        orbit_count=len(reps),counts_by_size=dict(sorted(Counter(s.nonlinear_size for s in reps).items())),
        supports=[s.to_json() for s in reps])
    search.old.write_json(a.output/"supports.json",support_record)
    results=[];replay=[];batches=[]
    for offset in range(0,len(reps),100):
        batch=reps[offset:offset+100]
        path=a.output/f"exact-{offset:05d}.json"
        ids=[s.identifier for s in batch]
        if path.exists():
            payload=json.loads(path.read_text())
            assert payload["support_ids"]==ids
        else:
            payload=dict(support_ids=ids,results=search.old.singular_batch_status(batch),
                         sympy_results=[search.old.groebner_status(s) for s in batch])
            assert [r["unit_ideal"] for r in payload["results"]]==[r["unit_ideal"] for r in payload["sympy_results"]]
            search.old.write_json(path,payload)
        results.extend(payload["results"]);replay.extend(payload["sympy_results"])
        batches.append(dict(path=path.name,sha256=sha(path)))
        print(json.dumps(dict(stage="algebra",done=len(results),total=len(reps),
                              nonunit=[r["support_id"] for r in results if not r["unit_ideal"]])),flush=True)
    manifest.update(
        status="complete",algebra_status="complete",
        orbit_count=len(reps),counts_by_size=support_record["counts_by_size"],
        supports_sha256=sha(a.output/"supports.json"),algebra_batches=batches,
        all_unit=all(r["unit_ideal"] for r in results),
        nonunit_ids=[r["support_id"] for r in results if not r["unit_ideal"]],
        algebra_backend_count=2 if reps else 0,elapsed_seconds=time.monotonic()-started)
    search.old.write_json(a.output/"manifest.json",manifest)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ("shards","gate_visits","algebra_batches")}),flush=True)
if __name__=="__main__":
    main()

