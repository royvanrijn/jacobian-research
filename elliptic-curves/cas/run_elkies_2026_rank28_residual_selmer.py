#!/usr/bin/env python3
"""Run a resource-bounded genuine 2-descent on the published rank-28 fibre.

The worker calls PARI ``ellrank`` through Sage with all 28 certified public
points supplied.  Unlike the BNF-free signature layer, a completed result is
the actual 2-Selmer dimension: PARI performs the cubic-field class-group and
local-solubility parts of Simon's 2-descent.  The supervisor records timeout
or memory stops as incomplete and therefore search-forbidden.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
CONTROL_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_residual_2selmer_gate_v1.json"
)
DEFAULT_SAGE = Path("/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python")

import sys

sys.path.insert(0, str(CAS))
from elkies_residual_selmer_gate import (  # noqa: E402
    INCOMPLETE_STATUS,
    SCHEMA,
    gate_record,
)


WORKER = r'''import json, sys, time
sys.path.insert(0, WORKER_CAS)
from sage.all import EllipticCurve, QQ, pari
from elkies_rank28 import GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS

started = time.monotonic()
curve = EllipticCurve(QQ, list(GENERAL_WEIERSTRASS_COEFFICIENTS))
two_torsion_dimension = int(curve.two_torsion_rank())
known = pari([[str(x), str(y)] for x, y in POINTS])
result = curve.pari_curve().ellrank(0, known)
lower = int(result[0])
upper = int(result[1])
sha_two_dimension = int(result[2])
total_selmer_dimension = upper + two_torsion_dimension + sha_two_dimension
print("ELKIES_R28_SELMER_JSON=" + json.dumps({
    "pari_ellrank_lower": lower,
    "pari_ellrank_upper": upper,
    "pari_sha_two_dimension": sha_two_dimension,
    "returned_independent_point_count": len(result[3]),
    "two_torsion_dimension": two_torsion_dimension,
    "total_two_selmer_dimension": total_selmer_dimension,
    "worker_seconds": time.monotonic() - started,
}, sort_keys=True), flush=True)
'''


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read_rss_bytes(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return 0


def stop_owned(process: subprocess.Popen[str], sig: signal.Signals) -> None:
    if process.poll() is None:
        os.killpg(process.pid, sig)


def parse_worker(stdout: str) -> dict[str, object] | None:
    prefix = "ELKIES_R28_SELMER_JSON="
    rows = [line[len(prefix) :] for line in stdout.splitlines() if line.startswith(prefix)]
    if len(rows) != 1:
        return None
    return json.loads(rows[0])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--rss-limit-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.timeout <= 0 or args.rss_limit_bytes < 64_000_000:
        parser.error("timeout must be positive and RSS limit at least 64MB")

    controls = json.loads(CONTROL_CERTIFICATE.read_text())
    rank28 = controls["fibres"][-1]
    if (
        controls.get("status")
        != "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS"
        or rank28["parameter"] != "-9529/5471"
        or rank28["locally_certified_rank_lower_bound"] != 28
    ):
        raise SystemExit("the exact rank-28 positive control is not available")
    sage_python = shutil.which(str(args.sage_python))
    if sage_python is None:
        raise SystemExit(f"Sage Python is unavailable: {args.sage_python}")

    worker_text = WORKER.replace("WORKER_CAS", repr(str(CAS)))
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as handle:
        handle.write(worker_text)
        worker_path = Path(handle.name)
    try:
        with tempfile.TemporaryFile(mode="w+") as stdout_file, tempfile.TemporaryFile(
            mode="w+"
        ) as stderr_file:
            process = subprocess.Popen(
                [sage_python, str(worker_path)],
                stdout=stdout_file,
                stderr=stderr_file,
                text=True,
                start_new_session=True,
            )
            started = time.monotonic()
            peak_rss = 0
            outcome = "running"
            while process.poll() is None:
                elapsed = time.monotonic() - started
                if elapsed >= args.timeout:
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
                time.sleep(0.25)
            wall_seconds = time.monotonic() - started
            stdout_file.seek(0)
            stderr_file.seek(0)
            stdout = stdout_file.read()
            stderr = stderr_file.read()
    finally:
        worker_path.unlink(missing_ok=True)

    worker = parse_worker(stdout)
    if outcome == "running":
        outcome = "completed" if process.returncode == 0 and worker else "backend_failure"
    if outcome == "completed" and worker is not None:
        gate = gate_record(
            total_two_selmer_dimension=int(worker["total_two_selmer_dimension"]),
            known_generic_rank=17,
            target_rank=32,
            two_torsion_dimension=int(worker["two_torsion_dimension"]),
        )
        status = str(gate["status"])
    else:
        gate = {
            "known_generic_rank": 17,
            "target_rank": 32,
            "required_residual_dimension": 15,
            "residual_two_selmer_quotient_dimension": None,
            "expensive_search_authorized": False,
            "decision": "no completed Selmer upper bound; expensive search remains forbidden",
        }
        status = INCOMPLETE_STATUS

    document = {
        "schema": SCHEMA,
        "status": status,
        "parameter": "-9529/5471",
        "global_minimal_model": rank28["minimal_model"],
        "known_rank_lower_bound": 28,
        "known_additional_directions_beyond_generic_17": 11,
        "directions_still_needed_for_rank_32": 4,
        "positive_control_certificate": {
            "path": str(CONTROL_CERTIFICATE.resolve()),
            "sha256": file_sha256(CONTROL_CERTIFICATE),
        },
        "descent_backend": {
            "name": "PARI ellrank through Sage",
            "algorithm": "Simon's complete 2-descent for curves over Q",
            "unconditional": outcome == "completed",
            "class_group_completeness_completed": outcome == "completed",
            "all_local_solubility_conditions_completed": outcome == "completed",
            "known_points_supplied": 28,
            "proof_boundary": (
                "These completeness flags are true only after PARI ellrank returns. "
                "A timeout is not converted into a Selmer or rank bound."
            ),
        },
        "backend_result": worker,
        "gate": gate,
        "supervisor": {
            "outcome": outcome,
            "returncode": process.returncode,
            "wall_seconds": wall_seconds,
            "peak_observed_rss_bytes": peak_rss,
            "timeout_seconds": args.timeout,
            "rss_limit_bytes": args.rss_limit_bytes,
            "stderr": stderr,
        },
        "stop_rule": (
            "No two-cover solving, ratpoints, slope-box, or other expensive point "
            "search is authorized unless status is PASS_RANK32_RESIDUAL_2_SELMER_GATE."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(
        f"ELKIESR28SELMER|outcome={outcome}|status={status}|"
        f"residual={gate['residual_two_selmer_quotient_dimension']}|output={args.output}"
    )


if __name__ == "__main__":
    main()
