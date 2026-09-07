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
            search.verify_record(old["search"])
            continue
        if max_new is not None and new >= max_new:
            break
        outcome = search.search(job["height"], job["seconds"],
            denominator_start=job.get("denominator_start", 1), denominator_end=job.get("denominator_end"),
            checkpoint_dir=checkpoint_dir)
        records[job["id"]] = {"id": job["id"], "search": outcome.record}
        new += 1
        result["results"] = [records[j["id"]] for j in jobs if j["id"] in records]
        complete = len(records) == len(jobs) and all(r["search"]["status"] == "bounded_search_complete" for r in records.values())
        result["status"] = "COMPLETE" if complete else "PARTIAL_CHECKPOINT"
        output.parent.mkdir(parents=True, exist_ok=True)
        temporary = output.with_suffix(".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
        temporary.replace(output)
        print(f"POINTED|{job['id']}|{outcome.record['status']}|points={len(outcome.curve_points)}", flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    choice = p.add_mutually_exclusive_group(required=True)
    choice.add_argument("--input", type=Path, help="generic job manifest")
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
