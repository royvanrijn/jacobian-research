#!/usr/bin/env sage -python
"""Measure known-section naive heights on minimal models for rational newfamily specializations.

Reads records from an existing JSON screen/search result, rebuilds the homogeneous
integral specialization and its 11 known hidden sections, maps to a global minimal
model, and reports the logarithmic naive x-heights of those known sections.

This is a calibration tool: it tells us whether an eclib search height such as
14 or 16 is even large enough to rediscover the known rank-11 subgroup on a
specialization. It does not prove rank.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import traceback

from sage.all import QQ, RR, ZZ

from search_unseeded_extra_points import build_integral_specialization


def naive_log_x(P):
    if P.is_zero():
        return 0.0
    x = QQ(P[0])
    n = ZZ(abs(x.numerator()))
    d = ZZ(x.denominator())
    m = max(n, d)
    if m <= 1:
        return 0.0
    return float(RR(m).log())


def run_single(args):
    started = time.monotonic()
    try:
        print("PHASE build_start", flush=True)
        E, known, _ = build_integral_specialization(
            args.numerator, args.denominator, args.sections_sobj
        )
        print("PHASE build_done", flush=True)

        print("PHASE minimal_start", flush=True)
        ms = time.monotonic()
        Emin = E.global_minimal_model()
        minimal_seconds = time.monotonic() - ms
        iso = E.isomorphism_to(Emin)
        print(f"PHASE minimal_done seconds={minimal_seconds:.6f}", flush=True)

        print("PHASE map_points_start", flush=True)
        points = [iso(P) for P in known]
        print("PHASE map_points_done", flush=True)

        heights = [naive_log_x(P) for P in points]
        heights_sorted = sorted(heights)

        result = {
            "status": "completed",
            "parameter": f"{args.numerator}/{args.denominator}",
            "numerator": args.numerator,
            "denominator": args.denominator,
            "discovery": args.discovery,
            "held": args.held,
            "minimal_seconds": minimal_seconds,
            "known_naive_heights": heights,
            "known_naive_min": heights_sorted[0],
            "known_naive_median": heights_sorted[5],
            "known_naive_max": heights_sorted[-1],
            "known_below_12": sum(h <= 12 for h in heights),
            "known_below_14": sum(h <= 14 for h in heights),
            "known_below_16": sum(h <= 16 for h in heights),
            "known_below_18": sum(h <= 18 for h in heights),
            "wall_seconds": time.monotonic() - started,
        }
        print(json.dumps(result, sort_keys=True), flush=True)
        return 0
    except Exception as exc:
        traceback.print_exc()
        print(json.dumps({
            "status": "error",
            "parameter": f"{args.numerator}/{args.denominator}",
            "numerator": args.numerator,
            "denominator": args.denominator,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "wall_seconds": time.monotonic() - started,
        }, sort_keys=True), flush=True)
        return 1


def run_parent(args):
    source = json.loads(Path(args.input_json).read_text())
    eligible = [r for r in source if r.get("status") == "completed"]
    if args.limit is not None:
        eligible = eligible[:args.limit]

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    records = []

    print(f"eligible={len(eligible)} timeout={args.timeout}", flush=True)
    for pos, r in enumerate(eligible, 1):
        a = int(r["numerator"])
        b = int(r["denominator"])
        print(f"[{pos}/{len(eligible)}] T={a}/{b}", flush=True)
        cmd = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--single",
            "--sections-sobj", args.sections_sobj,
            "--numerator", str(a),
            "--denominator", str(b),
            "--discovery", repr(r.get("discovery", 0.0)),
            "--held", repr(r.get("held", 0.0)),
        ]
        try:
            cp = subprocess.run(
                cmd,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=args.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            partial = exc.stdout or ""
            if isinstance(partial, bytes):
                partial = partial.decode(errors="replace")
            rec = {
                "status": "timeout",
                "parameter": f"{a}/{b}",
                "numerator": a,
                "denominator": b,
                "timeout_seconds": args.timeout,
                "output_tail": "\n".join(partial.splitlines()[-20:]),
            }
            records.append(rec)
            out.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
            print("  TIMEOUT", flush=True)
            for line in partial.splitlines()[-8:]:
                print("   ", line, flush=True)
            continue

        lines = [x for x in cp.stdout.splitlines() if x.strip()]
        try:
            rec = json.loads(lines[-1])
        except Exception:
            rec = {
                "status": "error",
                "parameter": f"{a}/{b}",
                "numerator": a,
                "denominator": b,
                "returncode": cp.returncode,
                "output_tail": "\n".join(lines[-30:]),
            }
        records.append(rec)
        out.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
        print(
            "  status=%s min=%.3f med=%.3f max=%.3f <=14=%s <=16=%s min_s=%.3f" % (
                rec.get("status"),
                rec.get("known_naive_min", float("nan")),
                rec.get("known_naive_median", float("nan")),
                rec.get("known_naive_max", float("nan")),
                rec.get("known_below_14"),
                rec.get("known_below_16"),
                rec.get("minimal_seconds", float("nan")),
            ),
            flush=True,
        )
        if rec.get("status") == "error":
            if rec.get("error"):
                print(
                    f"    ERROR {rec.get('error_type')}: {rec.get('error')}",
                    flush=True,
                )
            elif rec.get("output_tail"):
                for line in rec["output_tail"].splitlines()[-8:]:
                    print("   ", line, flush=True)

    done = [r for r in records if r.get("status") == "completed"]
    minima = sorted(r["known_naive_min"] for r in done)
    print(json.dumps({
        "attempted": len(records),
        "completed": len(done),
        "timeouts": sum(r.get("status") == "timeout" for r in records),
        "output": str(out),
        "minimum_known_naive_height": minima[0] if minima else None,
        "median_of_minima": minima[len(minima)//2] if minima else None,
    }, sort_keys=True), flush=True)
    return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-json")
    p.add_argument("--sections-sobj", default="/tmp/newfamily_hidden_sections_complete.sobj")
    p.add_argument("--timeout", type=int, default=20)
    p.add_argument("--limit", type=int)
    p.add_argument(
        "--output",
        default="artifacts/local/elliptic-curves/newfamily/specialized_section_heights.json",
    )
    p.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--numerator", type=int, help=argparse.SUPPRESS)
    p.add_argument("--denominator", type=int, help=argparse.SUPPRESS)
    p.add_argument("--discovery", type=float, default=0.0, help=argparse.SUPPRESS)
    p.add_argument("--held", type=float, default=0.0, help=argparse.SUPPRESS)
    args = p.parse_args()

    if args.single:
        if args.numerator is None or args.denominator is None:
            p.error("single mode requires numerator and denominator")
        return run_single(args)
    if not args.input_json:
        p.error("parent mode requires --input-json")
    return run_parent(args)


if __name__ == "__main__":
    raise SystemExit(main())
