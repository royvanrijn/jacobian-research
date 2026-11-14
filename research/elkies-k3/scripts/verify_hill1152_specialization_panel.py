#!/usr/bin/env python3
"""Replay the fixed Hill1152 eight-fibre/16-chart panel from its portable bundle.

No point enumeration, CVP, factorization, descent or specialization sweep.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import itertools
import json
from pathlib import Path
import tempfile
from zipfile import ZipFile

ROOT=Path(__file__).resolve().parents[2]
PRODUCER=Path(__file__).with_name("commission_hill1152_direct_charts.py")
GENERIC=ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2/four-generic-covers.json"
DEFAULT=ROOT/"artifacts/generated-results/elliptic-curves/hill1152_specializations_20261114_v1.zip"


def require(ok, message):
    if not ok: raise ValueError(message)


def digest(path): return sha256(path.read_bytes()).hexdigest()


def verify(bundle):
    with ZipFile(bundle) as archive:
        names=archive.namelist()
        require(len(names)==len(set(names)), "duplicate bundle member")
        require(all(not Path(n).is_absolute() and ".." not in Path(n).parts for n in names), "unsafe bundle path")
        data={n:archive.read(n) for n in names}
    packet=json.loads(data["direct/input.json"])
    panel=json.loads(data["direct/result.json"])
    first=json.loads(data["first/input.json"])
    require(packet["producer_sha256"]==digest(PRODUCER), "producer drift")
    require(packet["generic_certificate_sha256"]==digest(GENERIC), "generic certificate drift")
    require(packet["first_panel_input_sha256"]==sha256(data["first/input.json"]).hexdigest(), "first panel binding")
    require(panel["input_sha256"]==sha256(data["direct/input.json"]).hexdigest(), "panel binding")
    roster=[(i,str(v)) for i in (0,1,2,4) for v in (0,1)]
    require([(j["index"],j["parameter"]) for j in packet["jobs"]]==roster, "changed fibre roster")
    require([(j["index"],j["parameter"]) for j in first["jobs"]]==roster, "first roster changed")
    require(len(panel["jobs"])==8 and not panel["unrun_jobs"], "incomplete panel")
    spec=importlib.util.spec_from_file_location("hill1152_direct_replay",PRODUCER)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    rows=[]
    with tempfile.TemporaryDirectory(prefix="hill1152-panel-replay-") as temp:
        temp=Path(temp)
        for i,job in enumerate(packet["jobs"]):
            prefix=f"job-{i:02d}/"
            oldraw=data["first/"+prefix+"request.json"]
            require(job["source_request_sha256"]==sha256(oldraw).hexdigest(), "source request drift")
            old=json.loads(oldraw);new=json.loads(data["direct/"+prefix+"request.json"])
            require(new==job, "job input differs from frozen packet")
            scale=Q(job["scale_n"])
            oldcurve,newcurve=[list(map(Q,d["curve"])) for d in (old,new)]
            require(newcurve[:3]==oldcurve[:3]==[0,0,0], "short model")
            require(newcurve[3]*scale**8==oldcurve[3] and newcurve[4]*scale**12==oldcurve[4], "exact model transport")
            require(len(old["points"])==len(new["points"])==18, "starting point roster")
            for (ox,oy),(x,y) in zip(old["points"],new["points"]):
                ox,oy,x,y=map(Q,(ox,oy,x,y))
                require(x*scale**4==ox and y*scale**6==oy, "labelled point transport")
                require(y*y==x**3+newcurve[3]*x+newcurve[4], "off-curve starting point")
            require(new["centres"]==[[int(k==j) for k in range(18)] for j in (0,17)], "centre selection")
            search=json.loads(data["direct/"+prefix+"search.json"])
            require(search["request_sha256"]==sha256(data["direct/"+prefix+"request.json"]).hexdigest(), "search binding")
            require(len(search["charts"])==2 and search["status"]=="COMPLETE_DECLARED_CHARTS", "incomplete search")
            require(all(c["search"]["status"]=="bounded_search_complete" for c in search["charts"]), "partial point search")
            path=temp/prefix;path.mkdir(parents=True)
            (path/"request.json").write_bytes(data["direct/"+prefix+"request.json"])
            (path/"search.json").write_bytes(data["direct/"+prefix+"search.json"])
            from research_runtime.finite_reduction import ReductionCache
            from research_runtime.store import FactStore
            from research_runtime.mw_state import MWState
            fresh = ReductionCache(FactStore(temp/("fresh-cache-%02d" % i)))
            initial = MWState.from_record(search["initial_state"], cache=fresh)
            final = MWState.from_record(search["final_state"], cache=fresh)
            require(tuple(map(Q, initial.model.coefficients)) == tuple(newcurve), "initial state model differs from request")
            require(tuple(tuple(map(Q, P)) for P in initial.basis) == tuple(tuple(map(Q, P)) for P in new["points"]), "initial labelled basis differs from request")
            require(initial.rank == 18, "fresh finite arithmetic did not certify starting rank")
            checks=mod.replay(path/"search.json")
            require(final.rank == checks["final_rank"], "final rank did not replay from fresh arithmetic")
            require(checks=={k:panel["jobs"][i][k] for k in checks}, "panel summary not supported by replay")
            rows.append({"index":job["index"],"parameter":job["parameter"],**checks})
    return {"schema":"hill1152-specialization-panel-replay-v1","status":"PASS",
        "bundle_sha256":digest(bundle),"fibres":8,"complete_bounded_charts":16,
        "finite_reductions_recomputed_from_empty_caches":True,
        "rank_lower_bounds":[r["final_rank"] for r in rows],"new_certified_directions":sum(r["gain"] for r in rows),
        "results":rows,"discovery_wall_seconds":panel["charged_wall_seconds"],
        "boundary":"Eight independently admitted specialized subgroups and exact searched-map/point replay. No exact fibre rank, saturation, new global record, novelty, or population productivity conclusion."}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle",type=Path,default=DEFAULT);p.add_argument("--output",type=Path)
    a=p.parse_args();result=verify(a.bundle)
    if a.output:
        with a.output.open("x") as f: json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps({k:result[k] for k in ["status","fibres","complete_bounded_charts","rank_lower_bounds","new_certified_directions"]}))
