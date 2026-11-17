#!/usr/bin/env python3
"""Run/resume the finite 90-shard cover, then solve every residual support over Q."""
import argparse
import hashlib
import itertools
import json
import time
from collections import Counter
from pathlib import Path
import search

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--degree",type=int,default=7)
    p.add_argument("--bound",type=int,default=7)
    p.add_argument("--seconds",type=int,default=120)
    p.add_argument("--nodes",type=int,default=500000)
    p.add_argument("--output",type=Path,default=search.HERE/"results")
    p.add_argument("--independent",action="store_true")
    a=p.parse_args()
    compiler=search.SupportEnumerator(a.degree,a.bound)
    seeds=[(x,y,z) for x in range(2,a.degree+1)
           for y in range(2,a.degree+1) for z in range(y+1,a.degree+1)]
    all_supports=set()
    shards=[]
    incomplete=[]
    started=time.monotonic()
    for x,y,z in seeds:
        path=a.output/f"p{x}-q{y}-r{z}.json"
        shard=search.Shard(a.degree,a.bound,x,y,z,compiler)
        state=shard.run(path,a.seconds,a.nodes)
        shards.append(dict(seed=[x,y,z],path=path.name,
                           sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                           status=state["status"],visited=state["visited"]))
        if state["status"]!="complete":
            incomplete.append([x,y,z])
        all_supports.update(search.read_support(rows) for rows in state["supports"])
    if incomplete:
        search.write_json(a.output/"manifest.json",dict(status="incomplete",shards=shards,
                          incomplete=incomplete))
        return
    representatives=sorted({s.canonical_under_collision_stabilizer() for s in all_supports})
    search.write_json(a.output/"supports.json",dict(
        status="complete",degree=a.degree,bound=a.bound,
        oriented_support_count=len(all_supports),
        counts_by_size=dict(sorted(Counter(s.nonlinear_size for s in all_supports).items())),
        orbit_count=len(representatives),
        supports=[s.to_json() for s in representatives],
    ))
    results=[]
    sympy_results=[]
    batches=[]
    for offset in range(0,len(representatives),100):
        batch=representatives[offset:offset+100]
        path=a.output/f"exact-{offset:05d}.json"
        ids=[s.identifier for s in batch]
        if path.exists():
            payload=json.loads(path.read_text())
            assert payload["support_ids"]==ids
        else:
            rows=search.singular_batch_status(batch,timeout_seconds=300)
            payload=dict(support_ids=ids,results=rows)
            search.write_json(path,payload)
        if a.independent and "sympy_results" not in payload:
            payload["sympy_results"]=[search.groebner_status(s) for s in batch]
            assert [r["unit_ideal"] for r in payload["results"]]==[
                r["unit_ideal"] for r in payload["sympy_results"]]
            search.write_json(path,payload)
        results.extend(payload["results"])
        sympy_results.extend(payload.get("sympy_results",[]))
        batches.append(dict(path=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        print(json.dumps(dict(stage="algebra",done=len(results),total=len(representatives),
                              all_unit=all(r["unit_ideal"] for r in results))),flush=True)
    manifest=dict(
        schema="sparse-seven-cover.v1",status="complete",degree=a.degree,bound=a.bound,
        source_sha256=search.digest(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        legacy_sha256=hashlib.sha256((search.LEGACY/"jcsearch/low_degree_census.py").read_bytes()).hexdigest(),
        shards=shards,
        shard_count=len(seeds),
        oriented_support_count=len(all_supports),
        counts_by_size=dict(sorted(Counter(s.nonlinear_size for s in all_supports).items())),
        orbit_count=len(representatives),
        orbit_counts_by_size=dict(sorted(Counter(s.nonlinear_size for s in representatives).items())),
        supports_sha256=hashlib.sha256((a.output/"supports.json").read_bytes()).hexdigest(),
        algebra_batches=batches,
        all_unit=all(r["unit_ideal"] for r in results),
        independent_sympy=len(sympy_results)==len(results),
        nonunit_ids=[r["support_id"] for r in results if not r["unit_ideal"]],
        elapsed_seconds=time.monotonic()-started,
        scope="Normalized three-variable collisions, raw coordinate degree <= degree, nonlinear support <= bound.",
    )
    search.write_json(a.output/"manifest.json",manifest)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ("shards","algebra_batches")}),flush=True)
if __name__=="__main__":
    main()

