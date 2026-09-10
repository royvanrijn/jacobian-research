#!/usr/bin/env python3
"""Rolling conductor queue for sealed high-rank-foundry-v3 endpoints.

The source foundry is still running, so this queue deliberately discovers only
sealed PASS_CERTIFIED_SEARCH packets and adds them one by one.  Selection is
independent of conductor outcomes: every sealed packet is admitted, converted
to an integral model with an explicit Q-isomorphism, and sent through the
maintained conductor worker followed by an independent replay.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from fractions import Fraction
from pathlib import Path

from research_runtime.supervisor import Limits, run as supervise
from v3_warm_support import atomic, process_info, read, require, sha

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SOURCE = ROOT / "artifacts" / "generated-results" / "elliptic-curves" / "high-rank-foundry-v3"
D = ROOT / "artifacts" / "local" / "elliptic-curves" / "foundry-v3-conductors-v1"
OUT = ROOT / "artifacts" / "generated-results" / "elliptic-curves" / "foundry_v3_conductors_v1"
WORKER = CAS / "inventory_conductor_worker_v2.py"
SAGE = Path.home() / ".local" / "bin" / "sage"
SELF = Path(__file__).resolve()
DEFAULT_WORKERS = 2
BUILD_SECONDS = 1800
REPLAY_SECONDS = 120
RSS_BYTES = 2 * 1024**3


def source_jobs():
    rows = []
    for path in sorted(SOURCE.glob("job-*.json")):
        try:
            value = read(path)
        except (OSError, json.JSONDecodeError):
            continue
        if value.get("status") != "PASS_CERTIFIED_SEARCH":
            continue
        packet = value.get("packet")
        identifier = value.get("id")
        if not isinstance(packet, dict) or not identifier:
            continue
        require(value.get("certificate_replay", {}).get("status") ==
                "PASS_TWO_FINITE_IMPLEMENTATIONS", f"unreplayed source endpoint: {path}")
        rows.append((str(identifier), path, value))
    dedup = {}
    for identifier, path, value in rows:
        old = dedup.get(identifier)
        if old is None or path.name > old[1].name:
            dedup[identifier] = (identifier, path, value)
    return [dedup[k] for k in sorted(dedup)]


def integral_packet(identifier, record, source_path):
    packet = record["packet"]
    model = [Fraction(x) for x in packet["curve"]]
    require(len(model) == 5, f"short model expected: {identifier}")
    weights = (1, 2, 3, 4, 6)
    scale = 1
    for coefficient, weight in zip(model, weights):
        scale = math.lcm(scale, coefficient.denominator ** weight)
    curve = [str(coefficient * scale ** weight)
             for coefficient, weight in zip(model, weights)]
    return {
        "id": identifier,
        "original_curve": packet["curve"],
        "curve": curve,
        "integral_model_scale": str(scale),
        "rank_lower_bound": int(packet["rank_lower_bound"]),
        "known_primes": [],
        "factor_hints": [],
        "source_endpoint": str(source_path.relative_to(ROOT)),
        "source_endpoint_sha256": sha(source_path),
        "packet_sha256": record.get("packet_sha256"),
        "selection": "all sealed PASS_CERTIFIED_SEARCH endpoints; conductor-blind",
    }


def sync():
    D.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    discovered = []
    for identifier, source_path, record in source_jobs():
        case = D / "cases" / identifier
        input_path = case / "input.json"
        packet = integral_packet(identifier, record, source_path)
        atomic(input_path, packet, immutable=True)
        discovered.append({
            "id": identifier,
            "rank_lower_bound": packet["rank_lower_bound"],
            "source": packet["source_endpoint"],
            "source_sha256": packet["source_endpoint_sha256"],
            "input": str(input_path.relative_to(ROOT)),
        })
    manifest = {
        "schema": "elliptic-curves.foundry-v3-conductor-autorun.v1",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_status": read(SOURCE / "LIVE_STATUS.json").get("status")
        if (SOURCE / "LIVE_STATUS.json").exists() else "UNKNOWN",
        "workers": int(os.environ.get("EC_CONDUCTOR_WORKERS", DEFAULT_WORKERS)),
        "build_seconds": BUILD_SECONDS,
        "replay_seconds": REPLAY_SECONDS,
        "rss_bytes": RSS_BYTES,
        "selection": "all sealed PASS_CERTIFIED_SEARCH endpoints; no conductor-conditioned filtering",
        "source_script_sha256": sha(SELF),
        "worker_sha256": sha(WORKER),
        "ids": [r["id"] for r in discovered],
        "records": discovered,
    }
    atomic(D / "manifest.json", manifest)
    return manifest


def case(identifier):
    folder = D / "cases" / identifier
    input_path = folder / "input.json"
    output = folder / "conductor.json"
    result_path = folder / "result.json"
    if result_path.exists():
        return read(result_path)
    folder.mkdir(parents=True, exist_ok=True)
    packet = [str(SAGE), "-python", "-u", str(WORKER),
              "--packet", str(input_path), "--output", str(output)]
    build = folder / "build.supervisor.json"
    if not build.exists():
        build_record = supervise(packet, limits=Limits(BUILD_SECONDS, RSS_BYTES),
                                 log_path=folder / "build.log",
                                 checkpoint_path=build, cwd=ROOT)
    else:
        build_record = read(build)
    if not output.exists():
        value = {"id": identifier, "status": "RESOURCE_UNRESOLVED",
                 "build_outcome": build_record.get("outcome")}
        atomic(result_path, value, immutable=True)
        atomic(OUT / f"{identifier}.json", value, immutable=True)
        return value
    replay = folder / "replay.supervisor.json"
    if not replay.exists():
        replay_record = supervise(packet + ["--check"],
                                  limits=Limits(REPLAY_SECONDS, RSS_BYTES),
                                  log_path=folder / "replay.log",
                                  checkpoint_path=replay, cwd=ROOT)
    else:
        replay_record = read(replay)
    require(replay_record.get("outcome") == "completed",
            f"independent replay failed for {identifier}")
    certificate = read(output)
    value = {
        "id": identifier,
        "status": "PASS_INDEPENDENT_CONDUCTOR_REPLAY",
        "certificate": certificate,
        "build_outcome": build_record.get("outcome"),
        "bindings": {
            str(path.relative_to(ROOT)): sha(path)
            for path in (input_path, output, build, replay)
        },
    }
    atomic(result_path, value, immutable=True)
    atomic(OUT / f"{identifier}.json", value, immutable=True)
    return value


def dispatch(identifier):
    log = D / "cases" / identifier / "controller.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("ab", buffering=0) as stream:
        return subprocess.call([sys.executable, "-u", str(SELF), "case", "--id", identifier],
                               cwd=ROOT, stdin=subprocess.DEVNULL,
                               stdout=stream, stderr=subprocess.STDOUT)


def summary(manifest):
    rows = []
    for identifier in manifest.get("ids", []):
        path = D / "cases" / identifier / "result.json"
        if path.exists():
            rows.append(read(path))
    result = {
        "status": "RUNNING",
        "discovered": len(manifest.get("ids", [])),
        "completed": len(rows),
        "exact": sum(r.get("certificate", {}).get("status") == "EXACT" for r in rows),
        "unknown": sum(r.get("certificate", {}).get("status") == "UNKNOWN" for r in rows),
        "resource_unresolved": sum(r.get("status") == "RESOURCE_UNRESOLVED" for r in rows),
        "source_status": manifest.get("source_status"),
        "updated_at": time.time(),
    }
    if result["source_status"] == "COMPLETE" and result["completed"] == result["discovered"]:
        result["status"] = "COMPLETE"
    atomic(D / "summary.json", result)
    return result


def controller():
    pid = os.getpid()
    state = {"status": "RUNNING", "controller": process_info(pid),
             "started_unix": time.time(), "active": []}
    atomic(D / "state.json", state)
    workers = int(os.environ.get("EC_CONDUCTOR_WORKERS", DEFAULT_WORKERS))
    try:
        while True:
            manifest = sync()
            pending = [identifier for identifier in manifest["ids"]
                       if not (D / "cases" / identifier / "result.json").exists()]
            with ThreadPoolExecutor(max_workers=workers) as pool:
                active = {}
                while pending or active:
                    while pending and len(active) < workers:
                        identifier = pending.pop(0)
                        active[pool.submit(dispatch, identifier)] = identifier
                    state.update(active=list(active.values()), discovered=len(manifest["ids"]))
                    atomic(D / "state.json", state)
                    done, _ = wait(active, timeout=5, return_when=FIRST_COMPLETED)
                    for future in done:
                        identifier = active.pop(future)
                        code = future.result()
                        if code != 0:
                            atomic(D / "state.json", {**state, "status": "STOPPED_REVIEW_REQUIRED",
                                                       "failed_id": identifier, "returncode": code})
                            return
                        summary(manifest)
            current = read(SOURCE / "LIVE_STATUS.json").get("status") \
                if (SOURCE / "LIVE_STATUS.json").exists() else "UNKNOWN"
            manifest = sync()
            if current == "COMPLETE":
                result = summary(manifest)
                if result["completed"] == result["discovered"]:
                    state.update(status="COMPLETE", active=[], finished_unix=time.time())
                    atomic(D / "state.json", state)
                    return
            state.update(active=[], waiting_for_new_source_endpoints=True)
            atomic(D / "state.json", state)
            time.sleep(60)
    finally:
        lock = D / "controller.lock"
        try:
            lock.unlink()
        except FileNotFoundError:
            pass


def launch():
    D.mkdir(parents=True, exist_ok=True)
    sync()
    lock = D / "controller.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        raise SystemExit(f"controller lock exists: {lock}")
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    log = D / "autorun.log"
    child = subprocess.Popen([sys.executable, "-u", str(SELF), "controller"],
                             cwd=ROOT, stdin=subprocess.DEVNULL,
                             stdout=log.open("ab"), stderr=subprocess.STDOUT,
                             start_new_session=True)
    atomic(D / "launch.json", {"pid": child.pid, "launched_unix": time.time(),
                               "script_sha256": sha(SELF)})
    print("FOUNDRY_V3_CONDUCTOR_AUTORUN_LAUNCHED", child.pid, flush=True)


def status():
    manifest = sync()
    value = read(D / "summary.json") if (D / "summary.json").exists() else summary(manifest)
    state = read(D / "state.json") if (D / "state.json").exists() else {"status": "NOT_STARTED"}
    launch_record = read(D / "launch.json") if (D / "launch.json").exists() else {}
    value.update({"controller": state, "launch": launch_record,
                  "source_status": manifest.get("source_status"),
                  "discovered_now": len(manifest.get("ids", []))})
    print(json.dumps(value, indent=2))


def diagnose():
    status()
    for path in (D / "autorun.log", D / "state.json"):
        if path.exists():
            print(f"\n--- {path.relative_to(ROOT)} ---")
            text = path.read_text(errors="replace")
            print(text[-8000:])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["launch", "controller", "case", "status", "diagnose", "sync"])
    parser.add_argument("--id")
    args = parser.parse_args()
    if args.action == "launch":
        launch()
    elif args.action == "controller":
        controller()
    elif args.action == "case":
        require(args.id, "--id required")
        case(args.id)
    elif args.action == "sync":
        print(json.dumps(sync(), indent=2))
    elif args.action == "status":
        status()
    else:
        diagnose()
