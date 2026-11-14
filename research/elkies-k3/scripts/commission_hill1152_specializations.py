#!/usr/bin/env sage -python
"""Frozen eight-fibre commissioning panel for four newly certified MW18 covers.

Authorization: active November EC research goal. No known exceptional fibre
selects the panel. Fixed v in {0,1}, two raw factor-free lazy charts per fibre,
12 seconds per supervised worker, <=180 seconds total charged wall time.
"""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
import time

from sage.all import EllipticCurve, QQ

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/"elliptic-curves/cas"
sys.path.insert(0, str(CAS))
from research_runtime.supervisor import Limits, WorkerFailure, capture

LOADER = Path(__file__).with_name("verify_hill1152_four_covers.py")
PILOT = ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write("\n")


def load():
    spec = importlib.util.spec_from_file_location("hill1152_commissioning_loader", LOADER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def freeze(out):
    write(out/"input.json", {
        "schema": "hill1152-specialization-commissioning-v1",
        "loader_sha256": digest(LOADER), "producer_sha256": digest(Path(__file__)),
        "cover_certificate_sha256": digest(PILOT/"four-generic-covers.json"),
        "jobs": [{"index": i, "parameter": str(v)} for i in [0, 1, 2, 4] for v in [0, 1]],
        "next_holes": 2, "height": 10000, "seconds_per_chart": "0.5",
        "cvp_node_budget": 2000, "prime_bound": 500,
        "worker_wall_seconds": 12, "worker_rss_bytes": 1073741824,
        "total_wall_seconds": 180,
        "selection": "All four certified generic-only cover indices, Cartesian product with fixed rational parameters0,1. No exceptional fibre or specialized rank selects a job.",
        "boundary": "Commissioning panel only; no population productivity, rank upper bound or automatic extension."})


def run(out):
    started = time.monotonic()
    packet = json.loads((out/"input.json").read_text())
    assert packet["loader_sha256"] == digest(LOADER)
    assert packet["producer_sha256"] == digest(Path(__file__))
    assert packet["cover_certificate_sha256"] == digest(PILOT/"four-generic-covers.json")
    loader = load()
    cert = json.loads((PILOT/"four-generic-covers.json").read_text())
    rows = []
    seen = set()
    for number, job in enumerate(packet["jobs"]):
        if time.monotonic()-started+packet["worker_wall_seconds"]+2 >= packet["total_wall_seconds"]:
            break
        jobdir = out/f"job-{number:02d}"
        row = {"number": number, **job}
        try:
            generic, sections, _, mapping = loader.load_family(job["index"])
            v = QQ(job["parameter"])
            values = [a(v) for a in generic.a_invariants()]
            E = EllipticCurve(QQ, values)
            assert E.discriminant()
            points = [E(P[0](v), P[1](v)) for P in sections]
            u = QQ(0)
            # The actual old parameter is recorded only after the roster is frozen.
            # Formula uses the exact map, not a rank/exceptional-point label.
            pilot = loader.load_pilot()
            R = pilot.build(json.loads(pilot.SOURCE.read_text()))[0]
            n, d = R(mapping["N"])(v), R(mapping["D"])(v)
            old_t = QQ(json.loads((PILOT/"result.json").read_text())["common_chart"]["c"])+d/n
            identity = (str(E.j_invariant()), tuple(str(a) for a in E.a_invariants()))
            if identity in seen:
                row["status"] = "EXACT_MODEL_DUPLICATE"
                rows.append(row); write(jobdir/"receipt.json", row); continue
            seen.add(identity)
            H = next(f["height_matrix"] for f in cert["families"] if f["index"] == job["index"])
            request = {
                "schema": "elliptic-curves.lazy-mw-search.v1",
                "curve": [str(a) for a in E.a_invariants()],
                "points": [[str(x), str(y)] for P in points for x, y in [P.xy()]],
                "metric_gram": H,
                "policy": {"enumeration_backend": "gmp-pointed-sieve", "chart_metric_weight": "16", "diversity_window": 2},
                "next_holes": packet["next_holes"], "height": packet["height"],
                "seconds_per_chart": packet["seconds_per_chart"],
                "cvp_node_budget": packet["cvp_node_budget"], "prime_bound": packet["prime_bound"],
            }
            requestpath = jobdir/"request.json"
            write(requestpath, request)
            output = jobdir/"search.json"
            row.update({"old_parent_parameter": str(old_t), "request_sha256": digest(requestpath),
                        "starting_point_count": len(points), "specialized_starting_rank": "UNKNOWN"})
            result = capture([sys.executable, str(CAS/"run_mw_search.py"), "--request", str(requestpath),
                "--output", str(output), "--worker"],
                limits=Limits(packet["worker_wall_seconds"], packet["worker_rss_bytes"]),
                log_path=jobdir/"worker.log")
            search = json.loads(output.read_text())
            assert search["schema"] == "elliptic-curves.lazy-mw-search-result.v1"
            assert search["status"] != "RUNNING"
            row.update({"status": search["status"],
                        "specialized_starting_rank": len(search["initial_state"]["basis"]),
                        "final_certified_rank": len(search["final_state"]["basis"]),
                        "charts_completed": len(search["charts"]),
                        "gain": search["certified_rank_gain"],
                        "cvp_status": search["cvp_status"]})
            # Discovery receipts are immutable; replay uses a separate output.
            replay = jobdir/"replayed.json"
            capture([sys.executable, str(CAS/"run_mw_search.py"), "--verify", str(output),
                "--output", str(replay), "--worker"],
                limits=Limits(packet["worker_wall_seconds"], packet["worker_rss_bytes"]),
                log_path=jobdir/"replay.log")
            row["replay"] = "PASS"
        except WorkerFailure as error:
            row.update({"status": "WORKER_STOPPED", "worker_record": error.record,
                        "rank_upper_bound": "UNKNOWN"})
        except Exception as error:
            row.update({"status": "UNRESOLVED_PREPARATION_OR_REPLAY", "error": type(error).__name__+": "+str(error),
                        "rank_upper_bound": "UNKNOWN"})
        rows.append(row)
        write(jobdir/"receipt.json", row)
        print(json.dumps({k: row.get(k) for k in ["number", "index", "parameter", "status", "specialized_starting_rank", "final_certified_rank", "gain", "replay"]}), flush=True)
    result = {"schema": "hill1152-specialization-commissioning-result-v1",
        "input_sha256": digest(out/"input.json"), "jobs": rows,
        "planned_jobs": len(packet["jobs"]), "completed_receipts": len(rows),
        "unrun_jobs": packet["jobs"][len(rows):],
        "charged_wall_seconds": round(time.monotonic()-started, 6),
        "boundary": "Only exact subgroup admissions and replayed point-search results count. Timeouts, admission failures and CVP node stops do not bound fibre ranks or justify a larger campaign."}
    write(out/"result.json", result)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=["freeze", "run"])
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    (freeze if args.mode == "freeze" else run)(args.output)
