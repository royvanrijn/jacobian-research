#!/usr/bin/env sage -python
"""Search for rank gain above the explicit rank-11 newfamily subgroup.

This consumes records produced by ``screen_seeded_rational_candidates.py``.
Only specializations that already completed with baseline subgroup rank 11 are
eligible.  Each specialization is rebuilt in an isolated Sage subprocess,
its eleven hidden generic sections are specialized and processed by eclib,
and only then is ``mwrank_MordellWeil.search`` called.

Thus a reported change 11 -> 12+ is a genuine rank gain in the processed
rational subgroup; unlike earlier numerical-height triage it cannot be caused
by floating matrix rank.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from sage.all import EllipticCurve, QQ, ZZ
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

from screen_seeded_rational_candidates import (
    ROOTS,
    load_builder,
    load_sections_from_sobj,
    qbits,
    to_mwrank_triple,
)


def rebuild_seeded_group(args):
    parameter = QQ(args.numerator) / QQ(args.denominator)
    sections = load_sections_from_sobj(args.sections_sobj)
    family = load_builder()(ROOTS)

    A = QQ(family["Amin"](parameter))
    B = QQ(family["Bmin"](parameter))
    E = EllipticCurve(QQ, [0, 0, 0, A, B])
    if E.discriminant() == 0:
        raise RuntimeError("singular specialization")

    specialized = []
    for x_function, y_function in sections:
        P = E([QQ(x_function(parameter)), QQ(y_function(parameter))])
        specialized.append(P)

    if len(set(specialized)) != 11:
        raise RuntimeError("hidden sections collide at specialization")

    minimal_started = time.monotonic()
    Emin = E.global_minimal_model()
    minimal_seconds = time.monotonic() - minimal_started
    iso = E.isomorphism_to(Emin)
    points = [iso(P) for P in specialized]
    bits = [max(qbits(P[0]), qbits(P[1])) for P in points]

    order = sorted(range(11), key=lambda j: bits[j])
    mwcurve = mwrank_EllipticCurve([ZZ(value) for value in Emin.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1, maxr=32)

    seed_steps = []
    for j in order:
        started = time.monotonic()
        before = len(mw.points())
        mw.process([to_mwrank_triple(points[j])], saturation_bound=0)
        after = len(mw.points())
        seed_steps.append(
            {
                "section": j,
                "bits": bits[j],
                "rank_before": before,
                "rank_after": after,
                "seconds": time.monotonic() - started,
            }
        )

    baseline = len(mw.points())
    if baseline != 11:
        raise RuntimeError(f"expected seeded baseline 11, got {baseline}")

    return E, Emin, mw, bits, seed_steps, minimal_seconds


def run_single(args) -> int:
    started = time.monotonic()
    print(
        f"STAGE build T={args.numerator}/{args.denominator}",
        flush=True,
    )

    E, Emin, mw, bits, seed_steps, minimal_seconds = rebuild_seeded_group(args)

    print(
        "STAGE seeded baseline=11 "
        f"max_bits={max(bits)} minimal_seconds={minimal_seconds:.3f}",
        flush=True,
    )

    root_number = int(Emin.root_number())
    before_points = [tuple(map(int, P)) for P in mw.points()]
    before_rank = len(before_points)

    print(
        f"STAGE search height={args.height} root={root_number}",
        flush=True,
    )
    search_started = time.monotonic()
    mw.search(args.height, verbose=args.verbose_search)
    search_seconds = time.monotonic() - search_started

    after_points = [tuple(map(int, P)) for P in mw.points()]
    after_rank = len(after_points)
    gain = after_rank - before_rank

    result = {
        "status": "completed",
        "parameter": f"{args.numerator}/{args.denominator}",
        "numerator": args.numerator,
        "denominator": args.denominator,
        "discovery": args.discovery,
        "held": args.held,
        "root_number": root_number,
        "height_limit": args.height,
        "baseline_rank": before_rank,
        "final_subgroup_rank": after_rank,
        "rank_gain": gain,
        "minimal_generator_bits_min": min(bits),
        "minimal_generator_bits_median": sorted(bits)[5],
        "minimal_generator_bits_max": max(bits),
        "minimal_discriminant_bits": ZZ(abs(Emin.discriminant())).nbits(),
        "minimal_seconds": minimal_seconds,
        "seed_steps": seed_steps,
        "search_seconds": search_seconds,
        "basis_after": after_points if gain else None,
        "wall_seconds": time.monotonic() - started,
    }
    print(json.dumps(result, sort_keys=True), flush=True)
    return 0


def rank_priority(record):
    # Root +1 is interesting for an odd known rank under parity heuristics;
    # this is only a prioritization signal, never a rank proof.
    return (
        record.get("root_number") == 1,
        -record.get("minimal_generator_bits_max", 10**9),
        record.get("held", 0.0),
        record.get("discovery", 0.0),
    )


def run_parent(args) -> int:
    source = json.loads(Path(args.screen_json).read_text())
    eligible = [
        record
        for record in source
        if record.get("status") == "completed"
        and record.get("baseline_rank") == 11
    ]
    eligible.sort(key=rank_priority, reverse=True)
    if args.limit is not None:
        eligible = eligible[: args.limit]

    print(
        f"eligible={len(eligible)} height={args.height} timeout={args.timeout}",
        flush=True,
    )

    records = []
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    for position, record in enumerate(eligible, 1):
        numerator = int(record["numerator"])
        denominator = int(record["denominator"])
        print(
            f"[{position}/{len(eligible)}] T={numerator}/{denominator} "
            f"root={record.get('root_number')} "
            f"bits={record.get('minimal_generator_bits_max')} "
            f"D={record.get('discovery', 0):.4f} H={record.get('held', 0):.4f}",
            flush=True,
        )

        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--single",
            "--sections-sobj",
            args.sections_sobj,
            "--numerator",
            str(numerator),
            "--denominator",
            str(denominator),
            "--discovery",
            repr(record.get("discovery", 0.0)),
            "--held",
            repr(record.get("held", 0.0)),
            "--height",
            str(args.height),
        ]
        if args.verbose_search:
            command.append("--verbose-search")

        try:
            completed = subprocess.run(
                command,
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
            timeout_record = {
                "status": "timeout",
                "parameter": f"{numerator}/{denominator}",
                "numerator": numerator,
                "denominator": denominator,
                "height_limit": args.height,
                "timeout_seconds": args.timeout,
                "partial_output_tail": "\n".join(partial.splitlines()[-30:]),
            }
            records.append(timeout_record)
            print("  TIMEOUT", flush=True)
            if partial:
                for line in partial.splitlines()[-5:]:
                    print(f"    {line}", flush=True)
            continue

        lines = [line for line in completed.stdout.splitlines() if line.strip()]
        try:
            result = json.loads(lines[-1])
        except Exception:
            result = {
                "status": "error",
                "parameter": f"{numerator}/{denominator}",
                "returncode": completed.returncode,
                "output_tail": "\n".join(lines[-30:]),
            }
        records.append(result)
        print(
            "  status=%s rank=%s->%s gain=%s search_s=%s wall=%s"
            % (
                result.get("status"),
                result.get("baseline_rank"),
                result.get("final_subgroup_rank"),
                result.get("rank_gain"),
                result.get("search_seconds"),
                result.get("wall_seconds"),
            ),
            flush=True,
        )
        if result.get("rank_gain", 0) > 0:
            print("  *** EXACT SEEDED RANK GAIN ***", flush=True)

    output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
    completed_records = [r for r in records if r.get("status") == "completed"]
    gains = [r for r in completed_records if r.get("rank_gain", 0) > 0]
    print(
        json.dumps(
            {
                "output": str(output),
                "attempted": len(records),
                "completed": len(completed_records),
                "timeouts": sum(r.get("status") == "timeout" for r in records),
                "rank_gain_hits": len(gains),
                "maximum_final_subgroup_rank": max(
                    (r.get("final_subgroup_rank", -1) for r in completed_records),
                    default=-1,
                ),
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--screen-json",
        default="artifacts/local/elliptic-curves/newfamily/seeded_screen_top25.json",
    )
    parser.add_argument(
        "--sections-sobj",
        default="/tmp/newfamily_hidden_sections_complete.sobj",
    )
    parser.add_argument("--height", type=float, default=8.0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--verbose-search", action="store_true")
    parser.add_argument(
        "--output",
        default="artifacts/local/elliptic-curves/newfamily/seeded_rank_gain_h8.json",
    )

    parser.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--numerator", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--denominator", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--discovery", type=float, default=0.0, help=argparse.SUPPRESS)
    parser.add_argument("--held", type=float, default=0.0, help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args.single:
        if args.numerator is None or args.denominator is None:
            parser.error("single mode requires numerator and denominator")
        return run_single(args)
    return run_parent(args)


if __name__ == "__main__":
    raise SystemExit(main())
