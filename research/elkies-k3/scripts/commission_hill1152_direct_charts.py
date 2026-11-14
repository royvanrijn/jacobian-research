#!/usr/bin/env python3
"""Direct-chart successor to the zero-exposure Hill1152 CVP commissioning.

Uses the identical frozen eight fibres and exact generic section transport.
No outcome-based fibre selection. Two semantic centres: first parent section
and the new cover section. Undo N^4,N^6 polynomial scaling without factoring.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/"elliptic-curves/cas"
sys.path.insert(0, str(CAS))
from research_runtime.search_state import raw_state, reduction_cache
from research_runtime.mw_state import MWState
from research_runtime.supervisor import Limits, capture
from research_runtime.store import checkpoint, digest as exact_digest
from pointed_quartic_search import PointedQuarticSearch

FIRST = ROOT/"artifacts/local/hill1152-specializations-20261114-v1"
GENERIC = ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2/four-generic-covers.json"


def digest(path): return sha256(path.read_bytes()).hexdigest()


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(obj, stream, indent=2, sort_keys=True); stream.write("\n")


def value(cs, t):
    out = Q(0)
    for c in reversed(cs): out = out*t+Q(c)
    return out


def freeze(out):
    old = json.loads((FIRST/"input.json").read_text())
    generic = json.loads(GENERIC.read_text())
    jobs = []
    for i, row in enumerate(old["jobs"]):
        before = FIRST/f"job-{i:02d}/request.json"
        request = json.loads(before.read_text())
        mapping = next(f for f in generic["families"] if f["index"] == row["index"])
        n = value(mapping["N"], Q(row["parameter"]))
        assert n
        curve = [Q(a) for a in request["curve"]]
        curve[3] /= n**8; curve[4] /= n**12
        points = [[Q(x)/n**4, Q(y)/n**6] for x, y in request["points"]]
        assert all(y*y == x*x*x+curve[3]*x+curve[4] for x, y in points)
        job = {"id":i, **row, "source_request_sha256":digest(before),
            "scale_n":str(n),"curve":list(map(str,curve)),
            "points":[list(map(str,P)) for P in points],
            "centres":[[int(k==j) for k in range(18)] for j in [0,17]],
            "coordinate_policy":{"kind":"metric","weight":"16","matrix":[1,0,0,1]},
            "height":10000,"seconds":"0.5","prime_bound":500}
        write(out/f"job-{i:02d}/request.json",job); jobs.append(job)
    write(out/"input.json",{
        "schema":"hill1152-direct-chart-input-v1","producer_sha256":digest(Path(__file__)),
        "first_panel_input_sha256":digest(FIRST/"input.json"),
        "generic_certificate_sha256":digest(GENERIC),"jobs":jobs,
        "worker_seconds":15,"total_wall_seconds":180,"rss_bytes":1073741824,
        "reason":"The first panel has no point-search exposure: CVP budget stops. Preserve all eight fibres and use deterministic direct centres in an exactly smaller coordinate model.",
        "boundary":"This is a separately frozen successor; no first-panel budget is extended or receipt overwritten."})


def worker(request_path, output):
    request = json.loads(request_path.read_text())
    cache = reduction_cache()
    state = raw_state(request["curve"],request["points"],cache=cache,prime_bound=request["prime_bound"])
    if state.rank != 18: raise ArithmeticError("initial specialized subgroup remains unresolved")
    initial = state
    result={"schema":"hill1152-direct-chart-result-v1","request_sha256":digest(request_path),
            "status":"RUNNING","initial_state":initial.record(),"charts":[]}
    checkpoint(output,result)
    for centre in request["centres"]:
        # The original eighteen labelled generators remain a certified prefix.
        centre = centre+[0]*(state.rank-18)
        search = PointedQuarticSearch(state=state,centre={"coefficients":centre},
                                     coordinate_policy=request["coordinate_policy"])
        before = state
        state,outcome=search.search_state(request["height"],float(request["seconds"]),
            extra_primes=[211,223,227,229,233,239,241,251],checkpoint_dir=output.parent/"charts")
        result["charts"].append({"centre":centre,"search":outcome.record,
            "state_before":before.key,"state_after":state.key,"gain":state.rank-before.rank})
        result.update({"final_state":state.record(),"arithmetic_facts":cache.store.snapshot()})
        checkpoint(output,result)
    result.update({"status":"COMPLETE_DECLARED_CHARTS","initial_rank":initial.rank,
        "final_rank":state.rank,"gain":state.rank-initial.rank,
        "boundary":"Independent-point lower bounds only. Bounded misses do not imply upper bounds."})
    checkpoint(output,result)


def replay(path):
    result=json.loads(path.read_text())
    cache=reduction_cache(); cache.store.import_snapshot(result["arithmetic_facts"])
    state=MWState.from_record(result["initial_state"],cache=cache)
    assert state.rank == 18
    request=json.loads(path.with_name("request.json").read_text())
    assert result["request_sha256"]==digest(path.with_name("request.json"))
    for row in result["charts"]:
        search=PointedQuarticSearch(state=state,centre={"coefficients":row["centre"]},
                                   coordinate_policy=request["coordinate_policy"])
        outcome=search.verify_record(row["search"]); before=state
        for point in outcome.curve_points:
            state=state.adjoin(point,cache=cache,extra_primes=[211,223,227,229,233,239,241,251])
        assert (before.key,state.key,state.rank-before.rank)==(row["state_before"],row["state_after"],row["gain"])
    assert state.record()==result["final_state"]
    return {"initial_rank":18,"final_rank":state.rank,"gain":state.rank-18,"charts":len(result["charts"]),"replay":"PASS"}


def run(out):
    started=time.monotonic()
    packet=json.loads((out/"input.json").read_text())
    assert packet["producer_sha256"]==digest(Path(__file__))
    rows=[]
    for job in packet["jobs"]:
        # Reserve both discovery and replay budgets before starting a job.
        if time.monotonic()-started+2*packet["worker_seconds"]+2 >= packet["total_wall_seconds"]: break
        request=out/f"job-{job['id']:02d}/request.json"; output=request.with_name("search.json")
        row={"id":job["id"],"index":job["index"],"parameter":job["parameter"]}
        try:
            capture([sys.executable,str(Path(__file__)),"worker","--request",str(request),"--output",str(output)],
                limits=Limits(packet["worker_seconds"],packet["rss_bytes"]),log_path=request.with_name("worker.log"))
            checks=replay(output)
            row.update({"status":"COMPLETE",**checks})
            write(request.with_name("replay.json"),checks)
        except Exception as error:
            row.update({"status":"UNRESOLVED_WORKER_OR_REPLAY","error":type(error).__name__+": "+str(error)})
            if output.exists():
                retained=json.loads(output.read_text()); row["retained_search_status"]=retained["status"]
                row["completed_charts"]=len(retained["charts"])
        rows.append(row);write(request.with_name("receipt.json"),row)
        print(json.dumps(row),flush=True)
    write(out/"result.json",{"schema":"hill1152-direct-chart-panel-v1","input_sha256":digest(out/"input.json"),
        "jobs":rows,"planned":8,"unrun_jobs":packet["jobs"][len(rows):],
        "charged_wall_seconds":round(time.monotonic()-started,6),
        "boundary":"Complete subgroup and searched-point replay is required for a positive result. Partial workers, misses and unrun jobs remain UNKNOWN; no automatic expansion."})


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode",choices=["freeze","run","worker","replay"])
    p.add_argument("--output",type=Path);p.add_argument("--request",type=Path)
    a=p.parse_args()
    if a.mode=="freeze": freeze(a.output)
    elif a.mode=="run": run(a.output)
    elif a.mode=="worker": worker(a.request,a.output)
    else: print(json.dumps(replay(a.request)))
