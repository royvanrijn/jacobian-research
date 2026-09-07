#!/usr/bin/env python3
"""Run a bounded, stage-aware PARI S-class computation for Elkies rank 28.

This isolates the class-group layer which precedes local 2-descent.  The
worker consumes the proved factorization of the rank-28 2-division cubic
discriminant, certifies the maximal order, and then runs one of two jobs:

* ``class-quotient`` uses ``bnfinit(..., 0)`` and ``bnfcertify(..., 1)``.
  On completion, the computed class group provably surjects onto the true
  class group.  Its quotient by the bad-prime ideal classes consequently
  gives an unconditional upper bound for ``Cl(O_K[S^-1]) / 2``.
* ``full-units`` uses ``bnfinit(..., 1)`` and full ``bnfcertify``.  This is the
  stronger input needed to materialize global units for later squareclasses.

Neither mode tests local solubility or computes a 2-Selmer group.  Timeout,
memory exhaustion, provisional BNF data, and even a completed S-class bound
all leave expensive point search forbidden.
"""

from __future__ import annotations

import argparse
from collections import deque
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tempfile
import threading
import time
from typing import BinaryIO


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
    / "elkies_2026_rank28_s_class_pari_v1.json"
)
DEFAULT_SAGE = Path("/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python")
SCHEMA = "elliptic-curves.elkies-2026-rank28-s-class-pari.v1"
PROTOCOL = "ELKIESR28SCLASS"
RESULT_PREFIX = f"{PROTOCOL}|result="
TAIL_LIMIT = 64 * 1024


WORKER_TEMPLATE = r'''import json, time
from sage.all import pari

PROTOCOL = "ELKIESR28SCLASS"
COEFFICIENTS = __COEFFICIENTS__
BAD_PRIMES = __BAD_PRIMES__
PARI_STACK_BYTES = __PARI_STACK_BYTES__
BNF_FLAG = __BNF_FLAG__
CERTIFY_FLAG = __CERTIFY_FLAG__
TECH = __TECH__
DEBUG = __DEBUG__
FIELD_MODEL = __FIELD_MODEL__

def stage(name, status, **fields):
    suffix = "".join(f"|{key}={value}" for key, value in sorted(fields.items()))
    print(f"{PROTOCOL}|stage={name}|status={status}{suffix}", flush=True)

def packed_rank(rows):
    pivots = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return len(pivots)

started = time.monotonic()
pari.allocatemem(PARI_STACK_BYTES)
pari.default("debug", DEBUG)
pari.addprimes(BAD_PRIMES)
original_polynomial = pari(
    f"x^3+({COEFFICIENTS[2]})*x^2+({COEFFICIENTS[1]})*x+({COEFFICIENTS[0]})"
)
polynomial = original_polynomial
original_generator_in_field_model = "Mod(x,original_polynomial)"
if FIELD_MODEL == "polredabs":
    stage("polredabs", "start")
    stage_started = time.monotonic()
    reduced = pari.polredabs(original_polynomial, 1)
    polynomial = reduced[0]
    original_generator_in_field_model = str(reduced[1])
    stage(
        "polredabs",
        "complete",
        seconds=f"{time.monotonic()-stage_started:.6f}",
        polynomial=str(polynomial),
        original_generator=original_generator_in_field_model,
    )
elif FIELD_MODEL != "original":
    raise ValueError(f"unsupported field model {FIELD_MODEL}")
stage("nfinit", "start", factorization_supplied="true")
stage_started = time.monotonic()
nf = pari.nfinit([polynomial, BAD_PRIMES])
stage(
    "nfinit",
    "complete",
    seconds=f"{time.monotonic()-stage_started:.6f}",
    polynomial_discriminant=str(pari.poldisc(polynomial)),
    field_discriminant=str(nf[2]),
    defining_order_index=str(nf[3]),
    signature=":".join(str(value) for value in nf.nf_get_sign()),
)

stage("nfcertify", "start")
stage_started = time.monotonic()
obstructions = list(pari.nfcertify(nf))
if obstructions:
    raise ArithmeticError(f"maximal-order certification failed: {obstructions}")
stage("nfcertify", "complete", seconds=f"{time.monotonic()-stage_started:.6f}")

stage("bnfinit", "start", flag=BNF_FLAG, tech=":".join(str(value) for value in TECH))
stage_started = time.monotonic()
bnf = pari.bnfinit(nf, BNF_FLAG, TECH)
stage("bnfinit", "complete", seconds=f"{time.monotonic()-stage_started:.6f}")

stage("bnfcertify", "start", flag=CERTIFY_FLAG)
stage_started = time.monotonic()
certified = bool(pari.bnfcertify(bnf, CERTIFY_FLAG))
if not certified:
    raise ArithmeticError("PARI did not certify the requested BNF claim")
stage("bnfcertify", "complete", seconds=f"{time.monotonic()-stage_started:.6f}")

stage("s_ideal_classes", "start")
cyclics = [int(value) for value in bnf.bnf_get_cyc()]
class_number = int(bnf.bnf_get_no())
even_indices = [index for index, value in enumerate(cyclics) if value % 2 == 0]
s_rows = []
packed_rows = []
for rational_prime in BAD_PRIMES:
    for decomposition_index, prime_ideal in enumerate(pari.idealprimedec(nf, rational_prime)):
        coordinates = [
            int(value) for value in pari.bnfisprincipal(bnf, prime_ideal, 0)
        ]
        if len(coordinates) != len(cyclics):
            raise ArithmeticError("PARI returned a malformed ideal-class coordinate")
        mask = sum(
            (coordinates[index] & 1) << output_index
            for output_index, index in enumerate(even_indices)
        )
        packed_rows.append(mask)
        s_rows.append({
            "rational_prime": str(rational_prime),
            "decomposition_index": decomposition_index,
            "ramification_index": int(prime_ideal[2]),
            "residue_degree": int(prime_ideal[3]),
            "norm": str(pari.idealnorm(nf, prime_ideal)),
            "class_group_coordinates": [str(value) for value in coordinates],
            "mod2_even_cyclic_mask_hex": hex(mask),
        })
s_span_dimension = packed_rank(packed_rows)
class_mod2_dimension = len(even_indices)
result = {
    "pari_version": ".".join(str(value) for value in pari.version()),
    "field_model": FIELD_MODEL,
    "original_polynomial": str(original_polynomial),
    "polynomial": str(polynomial),
    "original_generator_in_field_model": original_generator_in_field_model,
    "polynomial_discriminant": str(pari.poldisc(polynomial)),
    "defining_order_index": str(nf[3]),
    "field_discriminant": str(nf[2]),
    "field_signature": [int(value) for value in nf.nf_get_sign()],
    "bnf_flag": BNF_FLAG,
    "bnfcertify_flag": CERTIFY_FLAG,
    "bnfcertify_completed": certified,
    "class_number": str(class_number),
    "class_group_cyclic_invariants": [str(value) for value in cyclics],
    "class_group_mod2_dimension": class_mod2_dimension,
    "even_cyclic_indices": even_indices,
    "s_prime_ideal_count": len(s_rows),
    "s_ideal_class_rows": s_rows,
    "s_ideal_class_span_dimension_mod2": s_span_dimension,
    "s_class_group_mod2_dimension_upper_bound": (
        class_mod2_dimension - s_span_dimension
    ),
    "worker_seconds": time.monotonic() - started,
}
stage("s_ideal_classes", "complete", count=len(s_rows))
print(PROTOCOL + "|result=" + json.dumps(result, sort_keys=True), flush=True)
'''


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def packed_rank(rows: list[int]) -> int:
    """Dependency-free GF(2) rank used by tests of the worker formula."""

    pivots: dict[int, int] = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = row
                break
            row ^= previous
    return len(pivots)


def s_class_mod2_upper_bound(cyclics: list[int], class_rows: list[list[int]]) -> int:
    """Dimension of the displayed class quotient modulo two."""

    even_indices = [index for index, value in enumerate(cyclics) if value % 2 == 0]
    if any(len(row) != len(cyclics) for row in class_rows):
        raise ValueError("class row length does not match cyclic invariants")
    masks = [
        sum(
            (row[index] & 1) << output_index
            for output_index, index in enumerate(even_indices)
        )
        for row in class_rows
    ]
    return len(even_indices) - packed_rank(masks)


class StreamCapture:
    """Drain one worker stream while retaining a hash and bounded tail."""

    def __init__(self, stream: BinaryIO, limit: int = TAIL_LIMIT) -> None:
        self.stream = stream
        self.limit = limit
        self.byte_count = 0
        self.digest = sha256()
        self.chunks: deque[bytes] = deque()
        self.retained = 0
        self.thread = threading.Thread(target=self._drain, daemon=True)

    def _drain(self) -> None:
        while True:
            chunk = self.stream.read(8192)
            if not chunk:
                return
            self.byte_count += len(chunk)
            self.digest.update(chunk)
            self.chunks.append(chunk)
            self.retained += len(chunk)
            while self.retained > self.limit and self.chunks:
                excess = self.retained - self.limit
                first = self.chunks[0]
                if len(first) <= excess:
                    self.chunks.popleft()
                    self.retained -= len(first)
                else:
                    self.chunks[0] = first[excess:]
                    self.retained -= excess

    def start(self) -> None:
        self.thread.start()

    def finish(self) -> dict[str, object]:
        self.thread.join()
        tail = b"".join(self.chunks).decode("utf-8", errors="replace")
        return {
            "byte_count": self.byte_count,
            "sha256": self.digest.hexdigest(),
            "tail": tail,
            "tail_truncated": self.byte_count > self.limit,
        }


def read_rss_bytes(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return 0


def stop_owned(process: subprocess.Popen[bytes], sig: signal.Signals) -> None:
    if process.poll() is None:
        os.killpg(process.pid, sig)


def parse_protocol(stdout_tail: str) -> tuple[list[dict[str, str]], dict | None]:
    events: list[dict[str, str]] = []
    results = []
    for line in stdout_tail.splitlines():
        if line.startswith(RESULT_PREFIX):
            results.append(json.loads(line[len(RESULT_PREFIX) :]))
        elif line.startswith(f"{PROTOCOL}|stage="):
            fields = {}
            for item in line.split("|")[1:]:
                key, value = item.split("=", 1)
                fields[key] = value
            events.append(fields)
    return events, results[0] if len(results) == 1 else None


def parse_pari_progress(stderr_tail: str) -> dict[str, object]:
    """Extract stable diagnostics from PARI debug text without interpreting it."""

    interesting = []
    relation_requests = []
    needles = (
        "LIMC =",
        "LIMC2 =",
        "factorbase",
        "relations needed",
        "relations remaining",
        "rnd_rel",
    )
    for line in stderr_tail.splitlines():
        if any(needle.lower() in line.lower() for needle in needles):
            interesting.append(line.strip())
        match = re.search(
            r"Look for (\d+) relations in (\d+) ideals \(([^)]+)\)", line
        )
        if match:
            relation_requests.append(
                {
                    "requested_relations": int(match.group(1)),
                    "candidate_ideal_count": int(match.group(2)),
                    "method": match.group(3),
                }
            )
    return {
        "matched_line_count_in_retained_tail": len(interesting),
        "matched_lines": interesting[-80:],
        "relation_request_count_in_retained_tail": len(relation_requests),
        "first_relation_request": relation_requests[0] if relation_requests else None,
        "last_relation_request": relation_requests[-1] if relation_requests else None,
        "minimum_candidate_ideal_count": (
            min(row["candidate_ideal_count"] for row in relation_requests)
            if relation_requests
            else None
        ),
        "interpretation": (
            "These are diagnostic PARI messages, not certified relation counts or "
            "a class-group completeness claim."
        ),
    }


def validate_inputs() -> tuple[dict, dict, list[int], list[int]]:
    controls = json.loads(CONTROL_CERTIFICATE.read_text())
    rank28 = controls["fibres"][-1]
    if (
        controls.get("status")
        != "PASS_EXACT_ELKIES_2026_HIGH_RANK_POSITIVE_CONTROLS"
        or rank28.get("parameter") != "-9529/5471"
        or rank28.get("locally_certified_rank_lower_bound") != 28
    ):
        raise ValueError("the exact rank-28 positive control is stale")

    ledger = json.loads(BAD_PLACE_LEDGER.read_text())
    factors = ledger.get("factorization", [])
    primes = [int(record["prime"]) for record in factors]
    coefficients = [int(value) for value in ledger["descent_cubic_coefficients_ascending"]]
    factor_product = 1
    for record in factors:
        factor_product *= int(record["prime"]) ** int(record["exponent"])
    valid = (
        ledger.get("status")
        == "COMPLETE_ALL_BAD_PLACE_KUMMER_IMAGES_NOT_A_SELMER_BOUND"
        and ledger.get("parameter") == "-9529/5471"
        and ledger.get("specialization", {}).get("global_minimal_model")
        == rank28["minimal_model"]
        and ledger.get("factorization_product_verified") is True
        and ledger.get("factor_primality_proof_completed") is True
        and ledger.get("all_bad_place_blocks_completed") is True
        and len(coefficients) == 4
        and coefficients[-1] == 1
        and len(primes) == 12
        and len(set(primes)) == 12
        and factor_product == int(ledger["descent_cubic_discriminant"])
    )
    if not valid:
        raise ValueError("the exact rank-28 bad-place ledger is stale")
    return controls, ledger, coefficients, primes


def worker_source(
    *,
    coefficients: list[int],
    primes: list[int],
    stack_bytes: int,
    mode: str,
    tech: list[float | int],
    debug: int,
    field_model: str = "original",
) -> str:
    bnf_flag = 0 if mode == "class-quotient" else 1
    certify_flag = 1 if mode == "class-quotient" else 0
    return (
        WORKER_TEMPLATE.replace("__COEFFICIENTS__", repr(coefficients))
        .replace("__BAD_PRIMES__", repr(primes))
        .replace("__PARI_STACK_BYTES__", str(stack_bytes))
        .replace("__BNF_FLAG__", str(bnf_flag))
        .replace("__CERTIFY_FLAG__", str(certify_flag))
        .replace("__TECH__", repr(tech))
        .replace("__DEBUG__", str(debug))
        .replace("__FIELD_MODEL__", repr(field_model))
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--rss-limit-bytes", type=int, default=12_000_000_000)
    parser.add_argument("--pari-stack-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE)
    parser.add_argument(
        "--mode", choices=("class-quotient", "full-units"), default="class-quotient"
    )
    parser.add_argument("--c1", type=float, default=0.1)
    parser.add_argument("--c2", type=float, default=4.0)
    parser.add_argument("--nrpid", type=int, default=20)
    parser.add_argument("--pari-debug", type=int, choices=(0, 1, 2), default=1)
    parser.add_argument(
        "--field-model", choices=("original", "polredabs"), default="original"
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("timeout must be positive")
    if args.rss_limit_bytes < 64_000_000 or args.pari_stack_bytes < 64_000_000:
        parser.error("memory limits must be at least 64 MB")
    if args.c1 <= 0 or args.c2 <= 0 or args.nrpid < 1:
        parser.error("BNF technical parameters must be positive")

    controls, ledger, coefficients, primes = validate_inputs()
    sage_python = shutil.which(str(args.sage_python))
    if sage_python is None:
        raise SystemExit(f"Sage Python is unavailable: {args.sage_python}")
    source = worker_source(
        coefficients=coefficients,
        primes=primes,
        stack_bytes=args.pari_stack_bytes,
        mode=args.mode,
        tech=[args.c1, args.c2, args.nrpid],
        debug=args.pari_debug,
        field_model=args.field_model,
    )

    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as handle:
        handle.write(source)
        worker_path = Path(handle.name)
    try:
        process = subprocess.Popen(
            [sage_python, str(worker_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=True,
        )
        assert process.stdout is not None and process.stderr is not None
        stdout_capture = StreamCapture(process.stdout)
        stderr_capture = StreamCapture(process.stderr)
        stdout_capture.start()
        stderr_capture.start()
        started = time.monotonic()
        peak_rss = 0
        outcome = "running"
        while process.poll() is None:
            elapsed = time.monotonic() - started
            if elapsed >= args.timeout:
                outcome = "strict_wall_timeout"
                stop_owned(process, signal.SIGTERM)
            try:
                peak_rss = max(peak_rss, read_rss_bytes(process.pid))
            except (FileNotFoundError, ProcessLookupError):
                pass
            if peak_rss > args.rss_limit_bytes and process.poll() is None:
                outcome = "strict_rss_limit"
                stop_owned(process, signal.SIGTERM)
            if outcome != "running":
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    stop_owned(process, signal.SIGKILL)
                    process.wait()
                break
            time.sleep(0.25)
        wall_seconds = time.monotonic() - started
        stdout_record = stdout_capture.finish()
        stderr_record = stderr_capture.finish()
    finally:
        worker_path.unlink(missing_ok=True)

    events, result = parse_protocol(str(stdout_record["tail"]))
    if outcome == "running":
        outcome = "completed" if process.returncode == 0 and result else "backend_failure"
    completed = outcome == "completed" and result is not None
    if completed:
        status = "PASS_UNCONDITIONAL_S_CLASS_MOD2_UPPER_BOUND_NOT_A_SELMER_BOUND"
    else:
        status = "INCOMPLETE_S_CLASS_COMPUTATION_SEARCH_FORBIDDEN"
    last_stage = events[-1] if events else None

    document = {
        "schema": SCHEMA,
        "status": status,
        "parameter": "-9529/5471",
        "global_minimal_model": controls["fibres"][-1]["minimal_model"],
        "mode": args.mode,
        "field_model": args.field_model,
        "input_certificates": {
            "positive_controls": {
                "path": str(CONTROL_CERTIFICATE.resolve()),
                "sha256": file_sha256(CONTROL_CERTIFICATE),
            },
            "bad_place_kummer_ledger": {
                "path": str(BAD_PLACE_LEDGER.resolve()),
                "sha256": file_sha256(BAD_PLACE_LEDGER),
                "factor_count": len(primes),
                "factorization_supplied_to_nfinit": True,
                "all_factors_proved_prime": True,
            },
        },
        "worker_source_sha256": sha256(source.encode()).hexdigest(),
        "backend_result": result,
        "class_group_claim": {
            "unconditional": completed,
            "scope": (
                "On completion, bnfcertify flag 1 proves that the true class "
                "group is a quotient of the computed group, so the displayed "
                "S-class mod-2 dimension is an upper bound. Full-units mode "
                "instead requires full bnfcertify."
            ),
        },
        "selmer_claim": {
            "completed": False,
            "residual_two_selmer_quotient_dimension": None,
            "all_local_solubility_conditions_completed": False,
            "expensive_search_authorized": False,
            "reason": (
                "An S-class quotient bound is only the global class-group layer; "
                "the norm condition and every local-solubility condition remain."
            ),
        },
        "supervisor": {
            "outcome": outcome,
            "returncode": process.returncode,
            "wall_seconds": wall_seconds,
            "peak_observed_rss_bytes": peak_rss,
            "timeout_seconds": args.timeout,
            "rss_limit_bytes": args.rss_limit_bytes,
            "pari_stack_bytes": args.pari_stack_bytes,
            "bnf_technical_parameters": {
                "c1": args.c1,
                "c2": args.c2,
                "nrpid": args.nrpid,
            },
            "pari_debug": args.pari_debug,
            "stage_events": events,
            "last_stage_event": last_stage,
            "stdout": stdout_record,
            "stderr": stderr_record,
            "pari_progress": parse_pari_progress(str(stderr_record["tail"])),
        },
        "stop_rule": (
            "No two-cover solving, ratpoints, slope-box, or other expensive point "
            "search is authorized by this artifact."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with args.output.open(mode) as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
    bound = None if result is None else result["s_class_group_mod2_dimension_upper_bound"]
    print(
        f"{PROTOCOL}|outcome={outcome}|status={status}|sclass_bound={bound}|"
        f"last_stage={None if last_stage is None else last_stage.get('stage')}|"
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()
