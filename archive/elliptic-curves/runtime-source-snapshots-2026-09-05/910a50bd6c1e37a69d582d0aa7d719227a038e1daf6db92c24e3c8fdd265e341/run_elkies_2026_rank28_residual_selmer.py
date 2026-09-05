#!/usr/bin/env python3
"""Run a resource-bounded genuine 2-descent on the published rank-28 fibre.

The PARI worker calls ``ellrank`` with all 28 certified public points supplied;
the independent eclib worker runs its invariant-quartic 2-descent in
``selmer_only`` mode with both point-search bounds zero.  A separate factored
PARI backend first supplies the proved factorization of the 2-division cubic
discriminant, avoiding a repeated hidden factorization inside ``ellrank``.
Unlike the BNF-free signature layer, a completed result from any backend is
the actual 2-Selmer dimension.  The supervisor records timeout or memory stops
as incomplete and therefore search-forbidden.
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
BAD_PLACE_LEDGER = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_rank28_bad_place_kummer_ledger_v1.json"
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
from build_elkies_2026_rank28_bad_place_ledger import (  # noqa: E402
    DISCRIMINANT_FACTORIZATION,
)


PARI_WORKER = r'''import json, sys, time
sys.path.insert(0, WORKER_CAS)
from sage.all import EllipticCurve, QQ, pari
from elkies_rank28 import GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS
from elkies_residual_selmer_gate import pari_ellrank_total_two_selmer_dimension

started = time.monotonic()
curve = EllipticCurve(QQ, list(GENERAL_WEIERSTRASS_COEFFICIENTS))
two_torsion_dimension = int(curve.two_torsion_rank())
known = pari([[str(x), str(y)] for x, y in POINTS])
result = curve.pari_curve().ellrank(0, known)
lower = int(result[0])
upper = int(result[1])
cassels_pairing_rank = int(result[2])
total_selmer_dimension = pari_ellrank_total_two_selmer_dimension(
    rank_lower=lower,
    rank_upper=upper,
    cassels_pairing_rank=cassels_pairing_rank,
    two_torsion_dimension=two_torsion_dimension,
)
print("ELKIES_R28_SELMER_JSON=" + json.dumps({
    "pari_ellrank_lower": lower,
    "pari_ellrank_upper": upper,
    "pari_cassels_pairing_quotient_rank": cassels_pairing_rank,
    "returned_independent_point_count": len(result[3]),
    "two_torsion_dimension": two_torsion_dimension,
    "total_two_selmer_dimension": total_selmer_dimension,
    "worker_seconds": time.monotonic() - started,
}, sort_keys=True), flush=True)
'''


PARI_FACTORED_WORKER = r'''import json, sys, time
sys.path.insert(0, WORKER_CAS)
from sage.all import EllipticCurve, QQ, pari
from build_elkies_2026_rank28_bad_place_ledger import DISCRIMINANT_FACTORIZATION
from elkies_rank28 import GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS
from elkies_residual_selmer_gate import pari_ellrank_total_two_selmer_dimension

started = time.monotonic()
pari.allocatemem(PARI_STACK_BYTES)
factor_hint_primes = [prime for prime, _exponent in DISCRIMINANT_FACTORIZATION]
pari.addprimes(factor_hint_primes)
curve = EllipticCurve(QQ, list(GENERAL_WEIERSTRASS_COEFFICIENTS))
two_torsion_dimension = int(curve.two_torsion_rank())
known = pari([[str(x), str(y)] for x, y in POINTS])
result = curve.pari_curve().ellrank(0, known)
lower = int(result[0])
upper = int(result[1])
cassels_pairing_rank = int(result[2])
total_selmer_dimension = pari_ellrank_total_two_selmer_dimension(
    rank_lower=lower,
    rank_upper=upper,
    cassels_pairing_rank=cassels_pairing_rank,
    two_torsion_dimension=two_torsion_dimension,
)
print("ELKIES_R28_SELMER_JSON=" + json.dumps({
    "pari_ellrank_lower": lower,
    "pari_ellrank_upper": upper,
    "pari_cassels_pairing_quotient_rank": cassels_pairing_rank,
    "returned_independent_point_count": len(result[3]),
    "two_torsion_dimension": two_torsion_dimension,
    "total_two_selmer_dimension": total_selmer_dimension,
    "factorization_supplied": True,
    "factor_hint_prime_count": len(factor_hint_primes),
    "pari_stack_bytes": PARI_STACK_BYTES,
    "worker_seconds": time.monotonic() - started,
}, sort_keys=True), flush=True)
'''


ECLIB_WORKER = r'''import json, sys, time
sys.path.insert(0, WORKER_CAS)
from sage.all import EllipticCurve, QQ
from elkies_rank28 import GENERAL_WEIERSTRASS_COEFFICIENTS

started = time.monotonic()
curve = EllipticCurve(QQ, list(GENERAL_WEIERSTRASS_COEFFICIENTS))
two_torsion_dimension = int(curve.two_torsion_rank())
backend = curve.mwrank_curve()
backend.two_descent(
    verbose=False,
    selmer_only=True,
    first_limit=0,
    second_limit=0,
    n_aux=40,
    second_descent=True,
)
total_selmer_dimension = int(backend.selmer_rank())
print("ELKIES_R28_SELMER_JSON=" + json.dumps({
    "eclib_rank_lower": int(backend.rank()),
    "eclib_rank_upper": int(backend.rank_bound()),
    "two_torsion_dimension": two_torsion_dimension,
    "total_two_selmer_dimension": total_selmer_dimension,
    "selmer_only": True,
    "first_point_search_limit": 0,
    "second_point_search_limit": 0,
    "auxiliary_prime_count": 40,
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


def validate_factor_hint_certificate(rank28_model: list[str]) -> dict[str, object]:
    """Reject a stale or partial bad-place ledger before starting PARI."""

    ledger = json.loads(BAD_PLACE_LEDGER.read_text())
    expected_factors = [
        {"prime": str(prime), "exponent": exponent}
        for prime, exponent in DISCRIMINANT_FACTORIZATION
    ]
    valid = (
        ledger.get("status")
        == "COMPLETE_ALL_BAD_PLACE_KUMMER_IMAGES_NOT_A_SELMER_BOUND"
        and ledger.get("parameter") == "-9529/5471"
        and ledger.get("specialization", {}).get("global_minimal_model")
        == rank28_model
        and ledger.get("factorization") == expected_factors
        and ledger.get("factorization_product_verified") is True
        and ledger.get("factor_primality_proof_completed") is True
        and ledger.get("all_bad_place_blocks_completed") is True
    )
    if not valid:
        raise ValueError("the exact rank-28 discriminant factor certificate is stale")
    return {
        "path": str(BAD_PLACE_LEDGER.resolve()),
        "sha256": file_sha256(BAD_PLACE_LEDGER),
        "discriminant": ledger["descent_cubic_discriminant"],
        "factor_count": len(expected_factors),
        "all_factors_proved_prime": True,
        "claim_boundary": (
            "This certificate supplies exact factor hints only. Its known-point "
            "Kummer rows are not themselves a Selmer upper bound."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=3600.0)
    parser.add_argument("--rss-limit-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE)
    parser.add_argument(
        "--backend", choices=("pari", "pari-factored", "eclib"), default="pari"
    )
    parser.add_argument("--pari-stack-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if (
        args.timeout <= 0
        or args.rss_limit_bytes < 64_000_000
        or args.pari_stack_bytes < 64_000_000
    ):
        parser.error("timeout must be positive and memory limits at least 64MB")

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

    factor_hint_certificate = None
    if args.backend == "pari-factored":
        factor_hint_certificate = validate_factor_hint_certificate(rank28["minimal_model"])
    worker_template = {
        "pari": PARI_WORKER,
        "pari-factored": PARI_FACTORED_WORKER,
        "eclib": ECLIB_WORKER,
    }[args.backend]
    worker_text = worker_template.replace("WORKER_CAS", repr(str(CAS))).replace(
        "PARI_STACK_BYTES", str(args.pari_stack_bytes)
    )
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
        "factor_hint_certificate": factor_hint_certificate,
        "descent_backend": {
            "name": (
                "eclib/mwrank through Sage"
                if args.backend == "eclib"
                else "PARI ellrank through Sage"
            ),
            "algorithm": (
                "complete invariant-quartic 2-descent in Selmer-only mode"
                if args.backend == "eclib"
                else (
                    "Simon's complete 2-descent over Q with a proved complete "
                    "2-division-discriminant factor table"
                    if args.backend == "pari-factored"
                    else "Simon's complete 2-descent for curves over Q"
                )
            ),
            "unconditional": outcome == "completed",
            "class_group_completeness_completed": outcome == "completed",
            "all_local_solubility_conditions_completed": outcome == "completed",
            "known_points_supplied": 0 if args.backend == "eclib" else 28,
            "point_search_enabled": False,
            "pari_effort": 0 if args.backend.startswith("pari") else None,
            "factorization_supplied": args.backend == "pari-factored",
            "proof_boundary": (
                "These completeness flags are true only after the selected backend "
                "returns its full Selmer result. For eclib the complete invariant-"
                "quartic calculation replaces a separate cubic class-group input. "
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
            "pari_stack_bytes": (
                args.pari_stack_bytes if args.backend == "pari-factored" else None
            ),
            "backend": args.backend,
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
        f"ELKIESR28SELMER|backend={args.backend}|outcome={outcome}|status={status}|"
        f"residual={gate['residual_two_selmer_quotient_dimension']}|output={args.output}"
    )


if __name__ == "__main__":
    main()
