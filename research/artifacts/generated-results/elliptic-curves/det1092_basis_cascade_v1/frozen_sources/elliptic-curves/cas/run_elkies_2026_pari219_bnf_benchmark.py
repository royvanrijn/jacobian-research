#!/usr/bin/env python3
"""Benchmark PARI 2.19's threaded BNF relation collector reproducibly.

This is a narrow global-engine probe, not a Selmer calculation.  A successful
run must finish both ``bnfinit`` and ``bnfcertify`` before the binary BNF
checkpoint is retained.  A timeout only records relation-collection progress.
"""

from __future__ import annotations

from research_runtime.supervisor import Limits, run as supervised_run, captured_run, preserve_previous
from research_runtime.store import checkpoint as atomic_checkpoint

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time
from typing import Any


SCHEMA = "elliptic-curves.elkies-2026-pari219-bnf-benchmark.v1"
PROTOCOL = "ELKIESR17PARI219BNF"
RANK21_REDUCED_CUBIC = (
    "x^3-x^2-774250153578278482962797863407542*x+"
    "4105678984643853583390832544029019669185034999158"
)
STAGE_RE = re.compile(
    r"^LIMC = (?P<limc>\d+), LIMC2 = (?P<limc2>\d+)$"
    r"|^KCZ = (?P<kcz>\d+), KC = (?P<kc>\d+), n = (?P<n>\d+)$"
)
DEFICIT_RE = re.compile(
    r"^#### Look for (?P<relations>\d+) relations in "
    r"(?P<ideals>\d+) ideals \(rnd_rel\)$"
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()






def gp_quote(path: Path) -> str:
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def parse_progress(log_path: Path) -> dict[str, Any]:
    stages: list[dict[str, int]] = []
    pending: dict[str, int] | None = None
    latest_deficit = None
    random_rounds = 0
    for line in log_path.read_text(errors="replace").splitlines():
        match = STAGE_RE.match(line)
        if match and match.group("limc") is not None:
            pending = {
                "limc": int(match.group("limc")),
                "limc2": int(match.group("limc2")),
            }
        elif match and pending is not None and match.group("kc") is not None:
            pending.update(
                {
                    "kcz": int(match.group("kcz")),
                    "factorbase_ideals": int(match.group("kc")),
                    "relation_target_with_units": int(match.group("n")),
                }
            )
            stages.append(pending)
            pending = None
        deficit = DEFICIT_RE.match(line)
        if deficit:
            random_rounds += 1
            latest_deficit = {
                "relations_requested": int(deficit.group("relations")),
                "ideals_requested": int(deficit.group("ideals")),
            }
    return {
        "factorbase_stages": stages,
        "random_relation_round_count": random_rounds,
        "latest_random_relation_deficit": latest_deficit,
    }


def gp_version(gp: Path) -> str:
    completed = captured_run(
        [str(gp), "-q", "-f"],
        input="print(version())\n",
        text=True,
        capture_output=True,
        check=True,
    )
    return completed.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp", type=Path, required=True)
    parser.add_argument("--polynomial", default=RANK21_REDUCED_CUBIC)
    parser.add_argument("--case-id", default="control-r21-t3_8")
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--rss-bytes",type=int,default=8_000_000_000)
    parser.add_argument("--stack-bytes", type=int, default=2_000_000_000)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--c1", type=float, default=0.03)
    parser.add_argument("--c2", type=float, default=4.0)
    parser.add_argument("--nrpid", type=int, default=-1)
    parser.add_argument("--max-factorizations", type=int, default=500)
    parser.add_argument("--ideal-power", type=int, default=1)
    parser.add_argument("--pari-debug", type=int, default=1)
    parser.add_argument("--source-commit", default="unknown")
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.timeout_seconds <= 0 or args.threads <= 0 or args.stack_bytes <= 0:
        parser.error("timeouts, thread counts, and stack sizes must be positive")
    for path in (args.log, args.checkpoint, args.output):
        if path.exists() and not args.overwrite:
            raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
    preserve_previous(args.checkpoint)

    gp = args.gp.resolve()
    if not gp.is_file():
        raise FileNotFoundError(gp)
    tech = [
        args.c1,
        args.c2,
        args.nrpid,
        args.max_factorizations,
        args.ideal_power,
        args.threads,
    ]
    tech_text = ",".join(str(value) for value in tech)
    program = f'''default(nbthreads,{args.threads});
setdebug("bnf",{args.pari_debug});
f={args.polynomial};
print("{PROTOCOL}|stage=bnfinit|status=start|tech={tech_text}");
b=bnfinit(f,1,[{tech_text}]);
print("{PROTOCOL}|stage=bnfinit|status=done|no=",b.no,"|cyc=",b.cyc);
print("{PROTOCOL}|stage=bnfcertify|status=start");
c=bnfcertify(b);
if(!c,error("bnfcertify returned zero"));
writebin("{gp_quote(args.checkpoint)}",b);
print("{PROTOCOL}|stage=bnfcertify|status=done|certified=1");
'''

    source_path=args.output.with_suffix('.gp')
    preserve_previous(source_path);source_path.write_text(program)
    supervision=supervised_run([str(gp),'-q','-f','-s',str(args.stack_bytes)],
        input_text=program,log_path=args.log,checkpoint_path=args.output.with_suffix('.supervisor.json'),
        limits=Limits(args.timeout_seconds,args.rss_bytes,pari_stack_bytes=args.stack_bytes))
    peak_rss=supervision['peak_observed_rss_bytes']
    elapsed=wall_seconds=supervision['wall_seconds']
    outcome='running' if supervision['outcome']=='completed' else supervision['outcome']
    log_text = args.log.read_text(errors="replace")
    certified = (
        f"{PROTOCOL}|stage=bnfcertify|status=done|certified=1" in log_text
        and args.checkpoint.is_file()
    )
    if outcome == "running":
        outcome = "completed_certified_bnf" if certified else "backend_failure"
    if not certified:
        preserve_previous(args.checkpoint)

    result = {
        "schema": SCHEMA,
        "status": outcome,
        "case_id": args.case_id,
        "claim_boundary": [
            "Only completed_certified_bnf supplies reusable global class/unit data.",
            "A timeout or relation deficit is not a class-group, Selmer, or rank bound.",
        ],
        "backend": {
            "name": "PARI/GP threaded Buchmann relation collector",
            "license": "GPL-2.0-or-later",
            "version_vector": gp_version(gp),
            "source_commit": args.source_commit,
            "binary": str(gp),
            "binary_sha256": file_sha256(gp),
        },
        "input": {
            "reduced_cubic": args.polynomial,
            "bnf_flag": 1,
            "bnf_tech": tech,
            "threads": args.threads,
            "stack_bytes": args.stack_bytes,
            "timeout_seconds": args.timeout_seconds,
        },
        "measurement": {
            "wall_seconds": wall_seconds,
            "peak_observed_rss_bytes": peak_rss,
            "returncode": supervision["returncode"],
            **parse_progress(args.log),
        },
        "log": str(args.log),
        "log_sha256": file_sha256(args.log),
        "checkpoint": str(args.checkpoint) if certified else None,
        "checkpoint_sha256": file_sha256(args.checkpoint) if certified else None,
    }
    result["supervision"]=supervision
    result["source"]={"path":str(source_path),"sha256":file_sha256(source_path)}
    preserve_previous(args.output)
    atomic_checkpoint(args.output,result)
    print(
        f"{PROTOCOL}|stage=complete|status={outcome}|seconds={wall_seconds:.3f}"
        f"|peak_rss={peak_rss}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
