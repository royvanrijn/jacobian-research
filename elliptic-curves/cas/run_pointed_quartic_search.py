#!/usr/bin/env python3
"""Bounded, checkpointed PointedQuarticSearch jobs for any rational family.

A job has id, curve, subgroup, centre, coordinate_policy, height, seconds.
Optional denominator_start/end split a projective box into disjoint shards.
The MW18 input adapter selects one certified specialization and caller-supplied
centres; it does not launch the 178-fibre campaign or choose a search budget.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

from pointed_quartic_search import PointedQuarticSearch, ROOT, canonical, sources
from research_runtime.search_state import reduction_cache
from research_runtime.mw_state import MWState
from research_runtime.store import checkpoint
from research_runtime.supervisor import preserve_previous

MW18 = ROOT/"artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json"
SCHEMA = "elliptic-curves.pointed-quartic-jobs.v1"


def mw18_jobs(document, candidate_id, centres, coordinate_policy, height, seconds):
    matches = [r for r in document["candidates"] if r["candidate_id"] == candidate_id]
    if len(matches) != 1:
        raise ValueError("MW18 candidate id is absent or ambiguous")
    row = matches[0]
    points = row["specialized_points"]
    subgroup = points["generic_R17"] + [points["cover_section"]]
    if len(subgroup) != 18 or not points["all_section_identities_verified_exactly"]:
        raise ValueError("expected eighteen exactly specialized sections")
    return {"schema": SCHEMA, "jobs": [
        {"id": f"{candidate_id}:centre:{index}", "curve": row["raw_short_model"],
         "subgroup": subgroup, "centre": centre, "coordinate_policy": coordinate_policy,
         "height": height, "seconds": seconds}
        for index, centre in enumerate(centres)]}


def execute(document, output, checkpoint_dir, max_new=None):
    if document.get("schema") != SCHEMA or not document.get("jobs"):
        raise ValueError("expected a nonempty pointed-quartic-jobs.v1 manifest")
    jobs = document["jobs"]
    if len({j["id"] for j in jobs}) != len(jobs):
        raise ValueError("job ids must be unique")
    # Validate the complete job envelope before starting any search.
    required = {"id", "curve", "subgroup", "centre", "coordinate_policy", "height", "seconds"}
    for job in jobs:
        if required-set(job) or set(job)-required-{"denominator_start", "denominator_end"}:
            raise ValueError("job has missing or unexpected fields")
    identity = {"manifest_sha256": sha256(canonical(document).encode()).hexdigest(), "sources": sources()}
    result = {"schema": "elliptic-curves.pointed-quartic-results.v1", **identity,
              "manifest": document, "status": "PARTIAL_CHECKPOINT", "results": [],
              "claim_boundary": "Exact point discovery and bounded coverage only. Quotient gains require a separate exact independence certificate."}
    output = Path(output)
    if output.exists():
        result = json.loads(output.read_text())
        if any(result[k] != v for k, v in identity.items()) or result["manifest"] != document:
            raise ArithmeticError("job checkpoint belongs to another manifest or implementation")
    records = {row["id"]: row for row in result["results"]}
    if len(records) != len(result["results"]) or not records.keys() <= {j["id"] for j in jobs}:
        raise ArithmeticError("job checkpoint contains an unrequested id")
    cache=reduction_cache()
    if 'arithmetic_facts' in result:cache.store.import_snapshot(result['arithmetic_facts'])
    states={};initial_states={}
    def admission(job,search,outcome):
        key=sha256(canonical({'curve':job['curve'],'subgroup':job['subgroup']}).encode()).hexdigest()
        initial_states.setdefault(key,search.state.record())
        before=states.get(key,search.state);after=before
        for point in outcome.curve_points:
            after=after.adjoin(point,cache=cache,extra_primes=[211,223,227,229,233,239,241,251])
        states[key]=after
        return {'group':key,'state_before':before.key,'state_after':after.key,'certified_rank_gain':after.rank-before.rank}
    new = 0
    for job in jobs:
        old = records.get(job["id"])
        search = PointedQuarticSearch(job["curve"], job["subgroup"], job["centre"], job["coordinate_policy"])
        if old and old["search"]["status"] == "bounded_search_complete":
            bounds = {"height_bound": job["height"], "timeout_seconds": job["seconds"],
                      "denominator_start": job.get("denominator_start", 1),
                      "denominator_end": job.get("denominator_end", job["height"])}
            if any(old["search"][k] != v for k, v in bounds.items()):
                raise ArithmeticError("saved job changed its search bounds")
            outcome=search.verify_record(old["search"])
            actual=admission(job,search,outcome)
            if 'state_admission' in old and old['state_admission']!=actual:raise ArithmeticError('retained incremental point admission differs')
            old['state_admission']=actual
            continue
        if max_new is not None and new >= max_new:
            break
        outcome = search.search(job["height"], job["seconds"],
            denominator_start=job.get("denominator_start", 1), denominator_end=job.get("denominator_end"),
            checkpoint_dir=checkpoint_dir)
        records[job["id"]] = {"id": job["id"], "search": outcome.record, "state_admission":admission(job,search,outcome)}
        new += 1
        result["results"] = [records[j["id"]] for j in jobs if j["id"] in records]
        complete = len(records) == len(jobs) and all(r["search"]["status"] == "bounded_search_complete" for r in records.values())
        result["status"] = "COMPLETE" if complete else "PARTIAL_CHECKPOINT"
        result.update({'initial_mw_states':initial_states,'mw_states':{k:v.record() for k,v in states.items()},
                       'arithmetic_facts':cache.store.snapshot(),
                       'claim_boundary':'Exact independent points are admitted incrementally. Bounded misses are not rank upper bounds.'})
        preserve_previous(output);checkpoint(output,result)
        print(f"POINTED|{job['id']}|{outcome.record['status']}|points={len(outcome.curve_points)}", flush=True)
    if states:
        result.update({'initial_mw_states':initial_states,'mw_states':{k:v.record() for k,v in states.items()},
                       'arithmetic_facts':cache.store.snapshot()})
        preserve_previous(output);checkpoint(output,result)
    return result


def verify_jobs(retained):
    """Replay maps and incremental columns, without a point census."""
    cache=reduction_cache()
    cache.store.import_snapshot(retained['arithmetic_facts'])
    states={k:MWState.from_record(v,cache=cache) for k,v in retained['initial_mw_states'].items()}
    jobs={j['id']:j for j in retained['manifest']['jobs']}
    initialized=set()
    for row in retained['results']:
        job=jobs[row['id']]
        search=PointedQuarticSearch(job['curve'],job['subgroup'],job['centre'],job['coordinate_policy'])
        outcome=search.verify_record(row['search']);admission=row['state_admission'];before=states[admission['group']]
        key=sha256(canonical({'curve':job['curve'],'subgroup':job['subgroup']}).encode()).hexdigest()
        if key!=admission['group']:raise ArithmeticError('point admission belongs to a different input group')
        if key not in initialized:
            if before!=search.state:raise ArithmeticError('initial MWState differs from the requested generators')
            initialized.add(key)
        if before.key!=admission['state_before']:raise ArithmeticError('point admission starts in a different state')
        after=before
        for point in outcome.curve_points:after=after.adjoin(point,cache=cache,extra_primes=[211,223,227,229,233,239,241,251])
        if after.key!=admission['state_after'] or after.rank-before.rank!=admission['certified_rank_gain']:
            raise ArithmeticError('incremental point admission differs')
        states[key]=after
    if {k:v.record() for k,v in states.items()}!=retained['mw_states']:raise ArithmeticError('retained subgroup states differ')
    return {'status':'PASS_RETAINED_POINTS_AND_INCREMENTAL_INDEPENDENCE','chart_count':len(retained['results']),
            'census_regenerated':False}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    choice = p.add_mutually_exclusive_group(required=True)
    choice.add_argument("--input", type=Path, help="generic job manifest")
    choice.add_argument("--verify",type=Path,help="replay retained chart and MWState witnesses")
    choice.add_argument("--mw18-candidate", help="one exact candidate id in the pinned MW18 specialization file")
    p.add_argument("--mw18-input", type=Path, default=MW18)
    p.add_argument("--centres", type=Path, help="JSON list of explicit centre specifications for the MW18 adapter")
    p.add_argument("--coordinate-policy", default="metric:16")
    p.add_argument("--height", type=int)
    p.add_argument("--seconds", type=float)
    p.add_argument("--export-jobs", type=Path, help="write a reviewable manifest without running it")
    p.add_argument("--output", type=Path)
    p.add_argument("--checkpoint-dir", type=Path, default=ROOT/"artifacts/local/elliptic-curves/pointed-quartic-search/charts")
    p.add_argument("--max-new", type=int)
    args = p.parse_args()
    if args.max_new is not None and args.max_new < 1:
        p.error("--max-new must be positive")
    if args.verify:
        result=verify_jobs(json.loads(args.verify.read_text()))
        if args.output:checkpoint(args.output,result)
        print(json.dumps(result,sort_keys=True));return
    if args.input:
        document = json.loads(args.input.read_text())
    else:
        if args.centres is None or args.height is None or args.seconds is None:
            p.error("the MW18 adapter requires explicit --centres, --height and --seconds")
        document = mw18_jobs(json.loads(args.mw18_input.read_text()), args.mw18_candidate,
            json.loads(args.centres.read_text()), args.coordinate_policy, args.height, args.seconds)
        document["source"] = {"path": str(args.mw18_input), "sha256": sha256(args.mw18_input.read_bytes()).hexdigest()}
    if args.export_jobs:
        args.export_jobs.parent.mkdir(parents=True, exist_ok=True)
        args.export_jobs.write_text(json.dumps(document, indent=2, sort_keys=True)+"\n")
    if args.output:
        execute(document, args.output, args.checkpoint_dir, args.max_new)
    elif not args.export_jobs:
        p.error("specify --output or --export-jobs")


if __name__ == "__main__":
    main()
