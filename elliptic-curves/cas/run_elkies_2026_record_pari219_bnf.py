#!/usr/bin/env python3
"""Build a certified PARI 2.19 BNF with proved discriminant-factor hints.

This is the record-fibre counterpart of the pinned generic BNF benchmark.
The prime list must come from an exact factorization certificate; it is passed
to PARI only to avoid rediscovering known factors.  No checkpoint survives
unless both ``bnfinit`` and ``bnfcertify`` complete.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import signal
import subprocess
import time


from run_elkies_2026_pari219_bnf_benchmark import (
    file_sha256,
    gp_quote,
    gp_version,
    parse_progress,
    read_rss_bytes,
    stop_group,
)


SCHEMA = "elliptic-curves.elkies-2026-record-pari219-bnf.v1"
PROTOCOL = "ELKIESR17RECORDPARI219BNF"
FACTOR_CERTIFICATE_STATUS = "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES"


def parse_factor_primes(value: str) -> tuple[int, ...]:
    primes = tuple(int(part) for part in value.split(",") if part.strip())
    if not primes or any(prime <= 1 for prime in primes):
        raise argparse.ArgumentTypeError("a nonempty list of primes exceeding one is required")
    return primes


def certified_factor_primes(path: Path, curve_id: int) -> tuple[int, ...]:
    document = json.loads(path.read_text())
    if document.get("status") != FACTOR_CERTIFICATE_STATUS:
        raise ArithmeticError("the supplied factor certificate is not passing")
    row = next(
        (record for record in document.get("records", []) if int(record["id"]) == curve_id),
        None,
    )
    if row is None:
        raise ArithmeticError("the factor certificate does not contain this curve")
    return tuple(int(prime) for prime in row["bad_primes"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp", type=Path, required=True)
    parser.add_argument("--polynomial", required=True)
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--curve-id", type=int, required=True)
    parser.add_argument("--factor-primes", type=parse_factor_primes, required=True)
    parser.add_argument("--factor-certificate", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=1800.0)
    parser.add_argument("--stack-bytes", type=int, default=5_000_000_000)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--relation-threads",
        type=int,
        help=(
            "PARI bnfinit usethr value; omit to reuse --threads, or pass zero "
            "to retain the serial collector's early-abort strategies"
        ),
    )
    parser.add_argument("--c1", type=float, default=0.3)
    parser.add_argument("--c2", type=float, default=4.0)
    parser.add_argument("--nrpid", type=int, default=20)
    parser.add_argument("--max-factorizations", type=int, default=10_000)
    parser.add_argument("--ideal-power", type=int, default=4)
    parser.add_argument("--pari-debug", type=int, default=1)
    parser.add_argument("--source-commit", default="unknown")
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.timeout_seconds <= 0 or args.threads <= 0 or args.stack_bytes <= 0:
        parser.error("timeouts, thread counts, and stack sizes must be positive")
    if args.relation_threads is not None and args.relation_threads < 0:
        parser.error("--relation-threads must be nonnegative")
    gp = args.gp.resolve()
    certificate = args.factor_certificate.resolve()
    for path in (gp, certificate):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.factor_primes != certified_factor_primes(certificate, args.curve_id):
        raise ArithmeticError("factor primes differ from the exact curve certificate")
    for path in (args.log, args.checkpoint, args.output):
        if path.exists() and not args.overwrite:
            raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
    args.checkpoint.unlink(missing_ok=True)

    relation_threads = (
        args.threads if args.relation_threads is None else args.relation_threads
    )
    tech = [
        args.c1,
        args.c2,
        args.nrpid,
        args.max_factorizations,
        args.ideal_power,
        relation_threads,
    ]
    tech_text = ",".join(str(value) for value in tech)
    factor_text = ",".join(str(prime) for prime in args.factor_primes)
    program = f'''default(nbthreads,{args.threads});
setdebug("bnf",{args.pari_debug});
addprimes([{factor_text}]);
f={args.polynomial};
print("{PROTOCOL}|stage=bnfinit|status=start|tech={tech_text}");
iferr(b=bnfinit(f,1,[{tech_text}]),E,print("{PROTOCOL}|stage=bnfinit|status=error|message=",E);quit(2));
if(type(b)!="t_VEC",print("{PROTOCOL}|stage=bnfinit|status=error|message=non_bnf_type_",type(b));quit(3));
print("{PROTOCOL}|stage=bnfinit|status=done|no=",b.no,"|cyc=",b.cyc);
print("{PROTOCOL}|stage=bnfcertify|status=start");
iferr(c=bnfcertify(b),E,print("{PROTOCOL}|stage=bnfcertify|status=error|message=",E);quit(4));
if(!c,print("{PROTOCOL}|stage=bnfcertify|status=error|message=returned_zero");quit(5));
iferr(writebin("{gp_quote(args.checkpoint.resolve())}",b),E,print("{PROTOCOL}|stage=checkpoint|status=error|message=",E);quit(6));
iferr(bb=read("{gp_quote(args.checkpoint.resolve())}"),E,print("{PROTOCOL}|stage=checkpoint|status=error|message=",E);quit(7));
if(type(bb)!="t_VEC",print("{PROTOCOL}|stage=checkpoint|status=error|message=non_bnf_type_",type(bb));quit(8));
iferr(cc=bnfcertify(bb),E,print("{PROTOCOL}|stage=checkpoint|status=error|message=",E);quit(9));
if(!cc,print("{PROTOCOL}|stage=checkpoint|status=error|message=reload_certification_returned_zero");quit(10));
print("{PROTOCOL}|stage=checkpoint|status=done|reload_certified=1");
print("{PROTOCOL}|stage=bnfcertify|status=done|certified=1");
'''

    started = time.monotonic()
    peak_rss = 0
    outcome = "running"
    with args.log.open("w") as handle:
        process = subprocess.Popen(
            [str(gp), "-q", "-f", "-s", str(args.stack_bytes)],
            stdin=subprocess.PIPE,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )
        assert process.stdin is not None
        process.stdin.write(program)
        process.stdin.close()
        try:
            while process.poll() is None:
                peak_rss = max(peak_rss, read_rss_bytes(process.pid))
                if time.monotonic() - started >= args.timeout_seconds:
                    outcome = "strict_wall_timeout"
                    stop_group(process, signal.SIGTERM)
                    try:
                        process.wait(timeout=15)
                    except subprocess.TimeoutExpired:
                        stop_group(process, signal.SIGKILL)
                        process.wait()
                    break
                time.sleep(0.25)
        except BaseException:
            stop_group(process, signal.SIGTERM)
            process.wait(timeout=15)
            raise

    elapsed = time.monotonic() - started
    log_text = args.log.read_text(errors="replace")
    certified = (
        f"{PROTOCOL}|stage=bnfcertify|status=done|certified=1" in log_text
        and f"{PROTOCOL}|stage=checkpoint|status=done|reload_certified=1" in log_text
        and "  ***" not in log_text
        and process.returncode == 0
        and args.checkpoint.is_file()
    )
    if outcome == "running":
        outcome = "completed_certified_bnf" if certified else "backend_failure"
    if not certified:
        args.checkpoint.unlink(missing_ok=True)

    result = {
        "schema": SCHEMA,
        "status": outcome,
        "case_id": args.case_id,
        "curve_id": args.curve_id,
        "claim_boundary": [
            "Only completed_certified_bnf supplies reusable class and unit data.",
            "A timeout or relation deficit is not a Selmer or rank bound.",
            "Factor hints accelerate exact arithmetic and do not replace bnfcertify.",
        ],
        "backend": {
            "version_vector": gp_version(gp),
            "source_commit": args.source_commit,
            "binary": str(gp),
            "binary_sha256": file_sha256(gp),
        },
        "input": {
            "reduced_cubic": args.polynomial,
            "bnf_tech": tech,
            "factor_hint_primes": [str(prime) for prime in args.factor_primes],
            "factor_certificate": str(certificate),
            "factor_certificate_sha256": file_sha256(certificate),
            "threads": args.threads,
            "relation_threads": relation_threads,
            "stack_bytes": args.stack_bytes,
            "timeout_seconds": args.timeout_seconds,
        },
        "measurement": {
            "wall_seconds": elapsed,
            "peak_observed_rss_bytes": peak_rss,
            "returncode": process.returncode,
            **parse_progress(args.log),
        },
        "log": str(args.log),
        "log_sha256": file_sha256(args.log),
        "checkpoint": str(args.checkpoint) if certified else None,
        "checkpoint_sha256": file_sha256(args.checkpoint) if certified else None,
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|case={args.case_id}|status={outcome}|seconds={elapsed:.3f}"
        f"|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
