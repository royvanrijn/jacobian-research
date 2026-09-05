#!/usr/bin/env python3
"""Run fail-closed residual 2-Selmer descents on the MW18 priority cohort.

The cohort is deterministic: the highest Nagao-score specialization on each
of the nine covers, plus every distinct known extreme-anchor fibre that occurs
among the 178 finalists.  Exact duplicate raw models are descended once.

Each worker supplies all eighteen independently certified points to PARI
``ellrank``.  Only a completed unconditional result becomes a 2-Selmer bound.
Timeouts, memory stops, and backend failures remain incomplete and do not
authorize the later serious point search.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
SPECIALIZATIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json"
)
GATE_HELPER = CAS / "elkies_residual_selmer_gate.py"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-residual-2selmer-priority-v1.json"
)
CHECKPOINT_DIRECTORY = (
    ROOT / "artifacts/local/elkies-k3/r17-extreme-anchored-mw18-residual-2selmer"
)
DEFAULT_SAGE = Path("/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python")
PROTOCOL = "R17ANCHORMW18SELMER"


WORKER = r'''import json
from pathlib import Path
import sys
import time

sys.path.insert(0, CAS_PATH)
from sage.all import EllipticCurve, QQ, pari
from elkies_residual_selmer_gate import pari_ellrank_total_two_selmer_dimension

payload = json.loads(Path(PAYLOAD_PATH).read_text())
pari.allocatemem(PARI_STACK_BYTES)
started = time.monotonic()
curve = EllipticCurve(QQ, [QQ(value) for value in payload["raw_short_model"]])
two_torsion_dimension = int(curve.two_torsion_rank())
points = payload["points"]
known = pari([[point["x"], point["y"]] for point in points])
result = curve.pari_curve().ellrank(0, known)
lower = int(result[0])
upper = int(result[1])
cassels = int(result[2])
total = pari_ellrank_total_two_selmer_dimension(
    rank_lower=lower,
    rank_upper=upper,
    cassels_pairing_rank=cassels,
    two_torsion_dimension=two_torsion_dimension,
)
print("R17ANCHORMW18SELMER_JSON=" + json.dumps({
    "pari_ellrank_lower": lower,
    "pari_ellrank_upper": upper,
    "pari_cassels_pairing_quotient_rank": cassels,
    "returned_independent_point_count": len(result[3]),
    "two_torsion_dimension": two_torsion_dimension,
    "total_two_selmer_dimension": total,
    "worker_seconds": time.monotonic() - started,
}, sort_keys=True), flush=True)
'''


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def read_rss_bytes(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return 0


def stop_owned(process: subprocess.Popen[str], sig: signal.Signals) -> None:
    if process.poll() is None:
        os.killpg(process.pid, sig)


def parse_worker(stdout: str):
    prefix = "R17ANCHORMW18SELMER_JSON="
    rows = [line[len(prefix) :] for line in stdout.splitlines() if line.startswith(prefix)]
    return json.loads(rows[0]) if len(rows) == 1 else None


def raw_model_key(candidate) -> str:
    return sha256(
        json.dumps(candidate["raw_short_model"], separators=(",", ":")).encode()
    ).hexdigest()


def priority_cohort(candidates):
    certified = [
        candidate
        for candidate in candidates
        if candidate["independence"]["status"]
        == "CERTIFIED_INTEGRALLY_INDEPENDENT_RANK_AT_LEAST_18"
    ]
    by_cover = {}
    for candidate in certified:
        label = candidate["cover_label"]
        if label not in by_cover or candidate["nagao"]["total_score"] > by_cover[label]["nagao"]["total_score"]:
            by_cover[label] = candidate
    selected = [(candidate, "top_nagao_score_on_cover") for candidate in by_cover.values()]
    selected.extend(
        (candidate, "known_extreme_anchor_surviving_nagao")
        for candidate in certified
        if candidate["nagao"]["is_certified_anchor"]
    )
    groups = {}
    for candidate, reason in selected:
        key = raw_model_key(candidate)
        group = groups.setdefault(
            key,
            {
                "representative": candidate,
                "candidate_ids": [],
                "selection_reasons": [],
            },
        )
        if candidate["candidate_id"] not in group["candidate_ids"]:
            group["candidate_ids"].append(candidate["candidate_id"])
        if reason not in group["selection_reasons"]:
            group["selection_reasons"].append(reason)
    result = list(groups.values())
    result.sort(
        key=lambda group: (
            -float(group["representative"]["nagao"]["total_score"]),
            group["representative"]["candidate_id"],
        )
    )
    return result


def supervise(group, args):
    candidate = group["representative"]
    points = candidate["specialized_points"]["generic_R17"] + [
        candidate["specialized_points"]["cover_section"]
    ]
    payload = {
        "candidate_id": candidate["candidate_id"],
        "raw_short_model": candidate["raw_short_model"],
        "points": points,
    }
    with tempfile.TemporaryDirectory(prefix="r17-anchored-mw18-selmer-") as directory:
        directory_path = Path(directory)
        payload_path = directory_path / "payload.json"
        worker_path = directory_path / "worker.py"
        payload_path.write_text(json.dumps(payload, separators=(",", ":")) + "\n")
        worker_path.write_text(
            f"CAS_PATH={str(CAS)!r}\n"
            f"PAYLOAD_PATH={str(payload_path)!r}\n"
            f"PARI_STACK_BYTES={args.pari_stack_bytes}\n"
            + WORKER
        )
        started = time.monotonic()
        process = subprocess.Popen(
            [str(args.sage_python), str(worker_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        peak_rss = 0
        outcome = "running"
        while process.poll() is None:
            elapsed = time.monotonic() - started
            if elapsed >= args.timeout_seconds:
                outcome = "strict_wall_timeout"
                stop_owned(process, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    stop_owned(process, signal.SIGKILL)
                    process.wait()
                break
            try:
                peak_rss = max(peak_rss, read_rss_bytes(process.pid))
            except (FileNotFoundError, ProcessLookupError):
                pass
            if peak_rss > args.rss_limit_bytes:
                outcome = "strict_rss_limit"
                stop_owned(process, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    stop_owned(process, signal.SIGKILL)
                    process.wait()
                break
            time.sleep(0.2)
        stdout, stderr = process.communicate()
        wall_seconds = time.monotonic() - started
    worker = parse_worker(stdout)
    if outcome == "running":
        outcome = "completed" if process.returncode == 0 and worker else "backend_failure"
    return {
        "representative_candidate_id": candidate["candidate_id"],
        "candidate_ids_for_exact_raw_model": sorted(group["candidate_ids"]),
        "selection_reasons": sorted(group["selection_reasons"]),
        "cover_label": candidate["cover_label"],
        "anchor_id": candidate["anchor_id"],
        "r": candidate["r"],
        "base_t": candidate["base_t"],
        "nagao_score": candidate["nagao"]["total_score"],
        "raw_model_maximum_coefficient_bits": candidate[
            "raw_model_maximum_coefficient_bits"
        ],
        "backend_result": worker,
        "supervisor": {
            "outcome": outcome,
            "returncode": process.returncode,
            "wall_seconds": wall_seconds,
            "peak_observed_rss_bytes": peak_rss,
            "timeout_seconds": args.timeout_seconds,
            "rss_limit_bytes": args.rss_limit_bytes,
            "pari_stack_bytes": args.pari_stack_bytes,
            "stderr_tail": stderr[-4000:],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specializations", type=Path, default=SPECIALIZATIONS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--checkpoint-directory", type=Path, default=CHECKPOINT_DIRECTORY)
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--rss-limit-bytes", type=int, default=3_000_000_000)
    parser.add_argument("--pari-stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    if args.timeout_seconds <= 0 or args.rss_limit_bytes <= 0:
        parser.error("resource bounds must be positive")
    if not args.sage_python.is_file():
        parser.error(f"Sage Python is unavailable: {args.sage_python}")

    source_bytes = args.specializations.read_bytes()
    source = json.loads(source_bytes)
    if source.get("status") != "COMPLETE_EXACT_MW18_FINALIST_SPECIALIZATION_AUDIT":
        raise ArithmeticError("the exact finalist specialization audit is not complete")
    if source.get("certified_rank_at_least_18_count") != 178:
        raise ArithmeticError("not all finalists have certified rank at least 18")
    cohort = priority_cohort(source["candidates"])
    if len({group["representative"]["cover_label"] for group in cohort}) != 9:
        raise ArithmeticError("the priority cohort lost a cover")
    input_hashes = {
        relative(args.specializations): sha256(source_bytes).hexdigest(),
        relative(GATE_HELPER): digest(GATE_HELPER),
        relative(Path(__file__).resolve()): digest(Path(__file__).resolve()),
    }
    input_key = sha256(
        json.dumps(input_hashes, separators=(",", ":"), sort_keys=True).encode()
    ).hexdigest()
    args.checkpoint_directory.mkdir(parents=True, exist_ok=True)

    import sys

    sys.path.insert(0, str(CAS))
    from elkies_residual_selmer_gate import gate_record

    records = []
    for index, group in enumerate(cohort, 1):
        candidate = group["representative"]
        checkpoint = args.checkpoint_directory / f"{candidate['candidate_id']}.json"
        record = None
        if not args.no_resume and checkpoint.is_file():
            saved = json.loads(checkpoint.read_text())
            if saved.get("input_key") == input_key:
                record = saved["record"]
        if record is None:
            print(
                f"{PROTOCOL}|candidate={index}/{len(cohort)}|"
                f"id={candidate['candidate_id']}|stage=descent|status=START",
                flush=True,
            )
            record = supervise(group, args)
            checkpoint.write_text(
                json.dumps(
                    {"input_key": input_key, "record": record},
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
        worker = record["backend_result"]
        if record["supervisor"]["outcome"] == "completed" and worker is not None:
            gate = gate_record(
                total_two_selmer_dimension=int(worker["total_two_selmer_dimension"]),
                known_generic_rank=18,
                target_rank=32,
                two_torsion_dimension=int(worker["two_torsion_dimension"]),
            )
            record["gate"] = gate
            record["status"] = gate["status"]
        else:
            record["gate"] = {
                "known_generic_rank": 18,
                "target_rank": 32,
                "required_residual_dimension": 14,
                "residual_two_selmer_quotient_dimension": None,
                "expensive_search_authorized": False,
                "decision": "no completed Selmer upper bound; serious point search remains forbidden",
            }
            record["status"] = "INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN"
        records.append(record)
        print(
            f"{PROTOCOL}|candidate={index}/{len(cohort)}|"
            f"id={candidate['candidate_id']}|outcome={record['supervisor']['outcome']}|"
            f"status={record['status']}",
            flush=True,
        )

    survivors = [
        record["representative_candidate_id"]
        for record in records
        if record["gate"].get("expensive_search_authorized") is True
    ]
    rejected = [
        record["representative_candidate_id"]
        for record in records
        if record["status"] == "REJECT_RANK32_BY_RESIDUAL_2_SELMER"
    ]
    incomplete = [
        record["representative_candidate_id"]
        for record in records
        if record["status"] == "INCOMPLETE_NO_SELMER_BOUND_SEARCH_FORBIDDEN"
    ]
    payload = {
        "schema": "elliptic-curves.r17-extreme-anchored-mw18-residual-2selmer-priority.v1",
        "status": "COMPLETE_FAIL_CLOSED_RESIDUAL_2SELMER_PRIORITY_CAMPAIGN",
        "selection": {
            "rule": "top Nagao score on each of nine covers plus each distinct surviving extreme-anchor raw fibre",
            "unique_raw_model_count": len(cohort),
            "source_finalist_count": 178,
        },
        "known_rank_lower_bound": 18,
        "target_rank": 32,
        "required_residual_dimension": 14,
        "records": records,
        "summary": {
            "exact_gate_survivor_count": len(survivors),
            "exact_gate_survivors": survivors,
            "exact_rejection_count": len(rejected),
            "exact_rejections": rejected,
            "incomplete_count": len(incomplete),
            "incomplete": incomplete,
            "serious_point_search_authorized_count": len(survivors),
        },
        "inputs": input_hashes,
        "input_key": input_key,
        "backend": {
            "name": "PARI ellrank through Sage",
            "algorithm": "complete 2-descent over Q when the worker returns",
            "pari_effort": 0,
            "all_eighteen_certified_points_supplied": True,
            "point_search_enabled": False,
            "unconditional_only_on_completed_rows": True,
        },
        "stop_rule": (
            "Only rows with PASS_RANK32_RESIDUAL_2_SELMER_GATE may enter the "
            "serious point-search stage; incomplete rows are not survivors."
        ),
        "claim_boundary": (
            "A worker timeout, memory stop, or failure is no Selmer or rank bound. "
            "Unselected finalists remain exact rank-at-least-18 candidates, not "
            "Selmer-screened candidates."
        ),
        "reproducing_command": (
            "python3 elliptic-curves/cas/"
            "run_r17_extreme_anchored_mw18_residual_selmer.py --no-resume"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|cohort={len(cohort)}|survivors={len(survivors)}|"
        f"rejected={len(rejected)}|incomplete={len(incomplete)}|"
        f"output={relative(args.output)}|status=COMPLETE",
        flush=True,
    )


if __name__ == "__main__":
    main()
