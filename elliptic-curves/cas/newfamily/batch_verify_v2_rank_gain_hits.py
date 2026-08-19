#!/usr/bin/env sage -python
"""Batch exact verification of v2 Schur-complement rank-gain hits.

Consumes one or more JSON outputs from search_unseeded_extra_points_v2.py.
All stable Schur hits are grouped by rational specialization.  Each
specialization is rebuilt/minimalized once, its eleven known hidden sections
are processed once through eclib, and then every distinct flagged candidate
point is processed into the same Mordell--Weil object.

This makes the exact stage efficient and, importantly, can expose more than one
new independent direction at a specialization: 11 -> 12 -> 13 -> ... .

A reported processed subgroup rank r proves rank(E(Q)) >= r.  This script does
not prove an upper bound or exact full Mordell--Weil rank.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import subprocess
import sys
import time

from sage.all import QQ, ZZ
from sage.libs.eclib.interface import mwrank_EllipticCurve, mwrank_MordellWeil

from search_unseeded_extra_points import build_integral_specialization
from screen_seeded_rational_candidates_fast import qbits, to_mwrank_triple


def parse_parameter(text: str):
    if "/" in text:
        a, b = text.split("/", 1)
        return int(a), int(b)
    return int(text), 1


def candidate_key(hit):
    # Minimal-model x/y are exact rational strings and stable within a source
    # specialization.  Normalize sign so Q and -Q are one discovery direction.
    x = QQ(hit["point_minimal"][0])
    y = QQ(hit["point_minimal"][1])
    yp = str(y)
    yn = str(-y)
    return (str(x), min(yp, yn))


def collect_hits(paths):
    grouped = defaultdict(dict)
    provenance = defaultdict(list)
    metadata = {}

    for path_text in paths:
        path = Path(path_text)
        rows = json.loads(path.read_text())
        for row in rows:
            if row.get("status") != "completed":
                continue
            parameter = row.get("parameter")
            if not parameter:
                continue
            metadata.setdefault(parameter, {
                "numerator": row.get("numerator"),
                "denominator": row.get("denominator"),
                "discovery": row.get("discovery"),
                "held": row.get("held"),
            })
            for index, hit in enumerate(row.get("candidate_points", [])):
                if not hit.get("numerical_new_direction"):
                    continue
                key = candidate_key(hit)
                old = grouped[parameter].get(key)
                if old is None:
                    grouped[parameter][key] = hit
                else:
                    # Keep the numerically stronger representative if duplicated
                    # across H levels / scans.
                    try:
                        old_rel = abs(QQ(old.get("high_relative_residual", "0")))
                        new_rel = abs(QQ(hit.get("high_relative_residual", "0")))
                        if new_rel > old_rel:
                            grouped[parameter][key] = hit
                    except Exception:
                        pass
                provenance[(parameter, key)].append({
                    "input_json": str(path),
                    "row_parameter": parameter,
                    "candidate_index": index,
                    "height_limit": row.get("height_limit"),
                    "high_relative_residual": hit.get("high_relative_residual"),
                })

    return grouped, provenance, metadata


def verify_specialization(args):
    payload = json.loads(Path(args.payload).read_text())
    parameter = payload["parameter"]
    hits = payload["hits"]
    a, b = parse_parameter(parameter)

    started = time.monotonic()
    print(f"STAGE build T={parameter} hits={len(hits)}", flush=True)
    E, known, _ = build_integral_specialization(a, b, args.sections_sobj)

    ms = time.monotonic()
    Emin = E.global_minimal_model()
    minimal_seconds = time.monotonic() - ms
    iso = E.isomorphism_to(Emin)
    known_min = [iso(P) for P in known]

    hit_points = []
    for h in hits:
        x = QQ(h["point_minimal"][0])
        y = QQ(h["point_minimal"][1])
        Q = Emin([x, y])
        if Q.is_zero():
            raise RuntimeError("flagged candidate is zero")
        hit_points.append(Q)

    # Exact point de-duplication up to sign after reconstruction.
    unique_hit_points = []
    unique_hits = []
    seen = set()
    for h, Q in zip(hits, hit_points):
        x = QQ(Q[0])
        y = QQ(Q[1])
        key = (str(x), min(str(y), str(-y)))
        if key in seen:
            continue
        seen.add(key)
        unique_hits.append(h)
        unique_hit_points.append(Q)

    labels = [f"U{i}" for i in range(11)] + [f"Q{i}" for i in range(len(unique_hit_points))]
    points = known_min + unique_hit_points
    bits = [max(qbits(P[0]), qbits(P[1])) for P in points]

    # Small points first dramatically helps eclib processing on these curves.
    order = sorted(range(len(points)), key=lambda i: bits[i])

    print(
        f"STAGE minimal seconds={minimal_seconds:.6f} "
        f"disc_bits={ZZ(abs(Emin.discriminant())).nbits()} "
        f"points={len(points)} max_bits={max(bits)}",
        flush=True,
    )

    root_number = None
    try:
        rs = time.monotonic()
        root_number = int(Emin.root_number())
        root_seconds = time.monotonic() - rs
    except Exception:
        root_seconds = None

    mwcurve = mwrank_EllipticCurve([ZZ(v) for v in Emin.ainvs()])
    mw = mwrank_MordellWeil(mwcurve, verbose=False, pp=1, maxr=max(32, len(points) + 8))

    steps = []
    rank_after_known = None
    hit_rank_increases = []
    known_processed = set()

    for position, index in enumerate(order, 1):
        before = len(mw.points())
        ss = time.monotonic()
        label = labels[index]
        print(
            f"STAGE process {position}/{len(points)} label={label} "
            f"bits={bits[index]} rank_before={before}",
            flush=True,
        )
        mw.process([to_mwrank_triple(points[index])], saturation_bound=0)
        after = len(mw.points())
        seconds = time.monotonic() - ss
        is_known = index < 11
        if is_known:
            known_processed.add(index)
        elif after > before:
            hit_index = index - 11
            hit_rank_increases.append({
                "hit_index": hit_index,
                "label": label,
                "rank_before": before,
                "rank_after": after,
                "point_minimal": unique_hits[hit_index].get("point_minimal"),
                "high_relative_residual": unique_hits[hit_index].get("high_relative_residual"),
            })
        print(
            f"STAGE processed label={label} rank={before}->{after} seconds={seconds:.6f}",
            flush=True,
        )
        steps.append({
            "label": label,
            "is_known_section": is_known,
            "bits": bits[index],
            "rank_before": before,
            "rank_after": after,
            "seconds": seconds,
        })

        # Because ordering interleaves hits and known sections, record the rank
        # at the first moment all 11 known sections have been processed.
        if rank_after_known is None and len(known_processed) == 11:
            rank_after_known = after

    final_rank = len(mw.points())
    result = {
        "status": "completed",
        "parameter": parameter,
        "numerator": a,
        "denominator": b,
        "input_hit_count": len(hits),
        "unique_hit_count": len(unique_hit_points),
        "known_subgroup_rank_when_all_known_processed": rank_after_known,
        "processed_subgroup_rank": final_rank,
        "proved_rank_at_least": final_rank,
        "proved_rank_at_least_12": bool(final_rank >= 12),
        "proved_rank_at_least_13": bool(final_rank >= 13),
        "independent_hit_count_observed": len(hit_rank_increases),
        "hit_rank_increases": hit_rank_increases,
        "root_number": root_number,
        "root_number_seconds": root_seconds,
        "minimal_seconds": minimal_seconds,
        "minimal_discriminant_bits": ZZ(abs(Emin.discriminant())).nbits(),
        "point_bits_min": min(bits),
        "point_bits_median": sorted(bits)[len(bits)//2],
        "point_bits_max": max(bits),
        "steps": steps,
        "wall_seconds": time.monotonic() - started,
    }
    print(json.dumps(result, sort_keys=True), flush=True)
    if final_rank >= 12:
        print(f"*** EXACT SUBGROUP RANK >= {final_rank} VERIFIED ***", flush=True)
    return 0


def run_parent(args):
    grouped, provenance, metadata = collect_hits(args.input_json)
    parameters = sorted(
        grouped,
        key=lambda p: (
            -len(grouped[p]),
            -float(metadata.get(p, {}).get("held") or 0.0),
            -float(metadata.get(p, {}).get("discovery") or 0.0),
            p,
        ),
    )

    if args.parameter:
        wanted = set(args.parameter)
        parameters = [p for p in parameters if p in wanted]
    if args.limit is not None:
        parameters = parameters[:args.limit]

    total_hits = sum(len(grouped[p]) for p in parameters)
    print(
        f"specializations={len(parameters)} unique_hits={total_hits} "
        f"timeout={args.timeout}",
        flush=True,
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload_dir = out.parent / (out.stem + "_payloads")
    payload_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for pos, parameter in enumerate(parameters, 1):
        hits = list(grouped[parameter].values())
        print(
            f"[{pos}/{len(parameters)}] T={parameter} hits={len(hits)} "
            f"D={metadata.get(parameter,{}).get('discovery')} "
            f"H={metadata.get(parameter,{}).get('held')}",
            flush=True,
        )

        safe = parameter.replace("/", "_")
        payload_path = payload_dir / f"T_{safe}.json"
        payload_path.write_text(json.dumps({
            "parameter": parameter,
            "metadata": metadata.get(parameter, {}),
            "hits": hits,
            "provenance": [
                {"key": list(key), "sources": provenance[(parameter, key)]}
                for key in grouped[parameter]
            ],
        }, indent=2, sort_keys=True) + "\n")

        cmd = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--single",
            "--payload", str(payload_path),
            "--sections-sobj", args.sections_sobj,
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
                "parameter": parameter,
                "unique_hit_count": len(hits),
                "timeout_seconds": args.timeout,
                "output_tail": "\n".join(partial.splitlines()[-30:]),
            }
            records.append(rec)
            out.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
            print("  TIMEOUT", flush=True)
            for line in partial.splitlines()[-8:]:
                print("   ", line, flush=True)
            continue

        lines = [line for line in cp.stdout.splitlines() if line.strip()]
        json_lines = []
        for line in lines:
            try:
                value = json.loads(line)
                if isinstance(value, dict) and value.get("status"):
                    json_lines.append(value)
            except Exception:
                pass
        if json_lines:
            rec = json_lines[-1]
        else:
            rec = {
                "status": "error",
                "parameter": parameter,
                "returncode": cp.returncode,
                "output_tail": "\n".join(lines[-30:]),
            }
        records.append(rec)
        out.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")

        print(
            "  status=%s hits=%s independent_hits=%s rank=%s root=%s wall=%s" % (
                rec.get("status"),
                rec.get("unique_hit_count"),
                rec.get("independent_hit_count_observed"),
                rec.get("processed_subgroup_rank"),
                rec.get("root_number"),
                rec.get("wall_seconds"),
            ),
            flush=True,
        )
        if (rec.get("processed_subgroup_rank") or 0) >= 12:
            print(
                f"  *** EXACT RANK LOWER BOUND >= {rec['processed_subgroup_rank']} ***",
                flush=True,
            )

    done = [r for r in records if r.get("status") == "completed"]
    verified12 = [r for r in done if (r.get("processed_subgroup_rank") or 0) >= 12]
    verified13 = [r for r in done if (r.get("processed_subgroup_rank") or 0) >= 13]
    maximum = max((r.get("processed_subgroup_rank", -1) for r in done), default=-1)
    print(json.dumps({
        "output": str(out),
        "input_files": args.input_json,
        "attempted_specializations": len(records),
        "completed_specializations": len(done),
        "timeouts": sum(r.get("status") == "timeout" for r in records),
        "verified_rank_at_least_12_specializations": len(verified12),
        "verified_rank_at_least_13_specializations": len(verified13),
        "maximum_processed_subgroup_rank": maximum,
    }, sort_keys=True), flush=True)
    return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-json", action="append", help="v2 search JSON; repeatable")
    p.add_argument("--sections-sobj", default="/tmp/newfamily_hidden_sections_complete.sobj")
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--limit", type=int)
    p.add_argument("--parameter", action="append", help="verify only this T; repeatable")
    p.add_argument(
        "--output",
        default="artifacts/local/elliptic-curves/newfamily/batch_exact_rank_gain_verify.json",
    )
    p.add_argument("--single", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--payload", help=argparse.SUPPRESS)
    args = p.parse_args()

    if args.single:
        if not args.payload:
            p.error("single mode requires --payload")
        return verify_specialization(args)
    if not args.input_json:
        p.error("provide at least one --input-json")
    return run_parent(args)


if __name__ == "__main__":
    raise SystemExit(main())
