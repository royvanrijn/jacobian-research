#!/usr/bin/env sage -python
"""Deep-search the strongest pinned newfamily specializations.

The target set is selected from the committed rank-gain batch.  For each target
and each requested search height, this driver synthesizes the tiny screen JSON
expected by ``search_unseeded_extra_points_v3.py``, runs the git-only discovery
stage, and if stable Schur hits appear, immediately runs the git-only exact
batch verifier on that discovery output.

This is a search driver, not a rank proof by itself.  Exact gains are accepted
only from ``batch_verify_v3_rank_gain_hits.py``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def parse_heights(text: str):
    vals = tuple(float(x.strip()) for x in text.split(",") if x.strip())
    if not vals or any(x <= 0 for x in vals):
        raise ValueError("--heights must be positive comma-separated numbers")
    return tuple(dict.fromkeys(vals))


def select_targets(batch, min_rank, parameters):
    records = batch.get("records", [])
    if parameters:
        wanted = set(parameters)
        rows = [r for r in records if r.get("parameter") in wanted]
        missing = wanted - {r.get("parameter") for r in rows}
        if missing:
            raise RuntimeError(f"parameters absent from pinned batch: {sorted(missing)}")
    else:
        rows = [r for r in records if (r.get("processed_subgroup_rank") or 0) >= min_rank]
    rows.sort(
        key=lambda r: (
            -(r.get("processed_subgroup_rank") or 0),
            -(r.get("exact_rank_gain_over_known") or 0),
            r.get("point_bits_median") or 10**9,
            r.get("parameter") or "",
        )
    )
    return rows


def run(cmd, timeout=None):
    return subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )


def last_json_line(text):
    for line in reversed(text.splitlines()):
        try:
            value = json.loads(line)
        except Exception:
            continue
        if isinstance(value, dict):
            return value
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--batch",
        type=Path,
        default=Path("artifacts/generated-results/elliptic-curves/newfamily_rank_gain_batch_v1.json"),
    )
    ap.add_argument("--parameter", action="append", help="specific T=a/b; repeatable")
    ap.add_argument("--min-rank", type=int, default=13)
    ap.add_argument("--heights", default="16,18,20")
    ap.add_argument("--timeout", type=int, default=900, help="per discovery specialization")
    ap.add_argument("--verify-timeout", type=int, default=900, help="per exact verification specialization")
    ap.add_argument("--precision", type=int, default=180)
    ap.add_argument("--verify-precision", type=int, default=260)
    ap.add_argument("--maxr", type=int, default=64)
    ap.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/local/elliptic-curves/newfamily/deep_rank_gain"),
    )
    args = ap.parse_args()

    batch = json.loads(args.batch.read_text())
    targets = select_targets(batch, args.min_rank, args.parameter)
    heights = parse_heights(args.heights)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    this_dir = Path(__file__).resolve().parent
    discovery_script = this_dir / "search_unseeded_extra_points_v3.py"
    verify_script = this_dir / "batch_verify_v3_rank_gain_hits.py"

    summary = []
    print(
        f"DEEP|targets={len(targets)}|heights={','.join(str(h) for h in heights)}"
        f"|timeout={args.timeout}|verify_timeout={args.verify_timeout}",
        flush=True,
    )

    for target in targets:
        parameter = target["parameter"]
        a, b = int(target["numerator"]), int(target["denominator"])
        safe = parameter.replace("/", "_")
        print(
            f"DEEP|target={parameter}|pinned_rank={target['processed_subgroup_rank']}"
            f"|root={target.get('root_number')}|hits={target.get('unique_hit_count')}",
            flush=True,
        )

        screen_path = args.output_dir / f"screen_T_{safe}.json"
        screen_path.write_text(json.dumps([{
            "status": "completed",
            "parameter": parameter,
            "numerator": a,
            "denominator": b,
            "height_rank": 11,
            "positive_definite": True,
            "discovery": 0.0,
            "held": 0.0,
            "generator_bits_max": target.get("point_bits_max"),
        }], indent=2) + "\n")

        target_summary = {
            "parameter": parameter,
            "pinned_rank": target["processed_subgroup_rank"],
            "root_number": target.get("root_number"),
            "runs": [],
        }

        for height in heights:
            htag = str(height).replace(".", "p")
            discovery_out = args.output_dir / f"T_{safe}_h{htag}_discovery.json"
            discovery_log = args.output_dir / f"T_{safe}_h{htag}_discovery.log"
            cmd = [
                sys.executable,
                str(discovery_script),
                "--screen-json", str(screen_path),
                "--height", str(height),
                "--timeout", str(args.timeout),
                "--precision", str(args.precision),
                "--verify-precision", str(args.verify_precision),
                "--maxr", str(args.maxr),
                "--output", str(discovery_out),
            ]
            print(f"DEEP|stage=discover|T={parameter}|height={height}|status=start", flush=True)
            started = time.monotonic()
            try:
                cp = run(cmd, timeout=args.timeout + 60)
                text = cp.stdout
                discovery_log.write_text(text)
                elapsed = time.monotonic() - started
                rows = json.loads(discovery_out.read_text()) if discovery_out.exists() and discovery_out.read_text().strip() else []
                rec = rows[0] if rows else None
                hits = (rec or {}).get("numerical_new_direction_hits", 0)
                print(
                    f"DEEP|stage=discover|T={parameter}|height={height}|status={(rec or {}).get('status','error')}"
                    f"|found={(rec or {}).get('mwrank_points_found')}|schur_hits={hits}|seconds={elapsed:.3f}",
                    flush=True,
                )
            except subprocess.TimeoutExpired as exc:
                text = exc.stdout or ""
                if isinstance(text, bytes):
                    text = text.decode(errors="replace")
                discovery_log.write_text(text)
                rec = {"status": "driver_timeout"}
                hits = 0
                elapsed = time.monotonic() - started
                print(f"DEEP|stage=discover|T={parameter}|height={height}|status=driver_timeout", flush=True)

            run_row = {
                "height": height,
                "discovery_status": rec.get("status") if rec else "error",
                "mwrank_points_found": rec.get("mwrank_points_found") if rec else None,
                "schur_hits": hits,
                "discovery_seconds": elapsed,
                "discovery_json": str(discovery_out),
                "discovery_log": str(discovery_log),
            }

            if hits:
                verify_out = args.output_dir / f"T_{safe}_h{htag}_exact.json"
                verify_log = args.output_dir / f"T_{safe}_h{htag}_exact.log"
                vcmd = [
                    sys.executable,
                    str(verify_script),
                    "--input-json", str(discovery_out),
                    "--parameter", parameter,
                    "--timeout", str(args.verify_timeout),
                    "--output", str(verify_out),
                ]
                print(f"DEEP|stage=verify|T={parameter}|height={height}|status=start", flush=True)
                started_v = time.monotonic()
                try:
                    vcp = run(vcmd, timeout=args.verify_timeout + 60)
                    verify_log.write_text(vcp.stdout)
                    vrows = json.loads(verify_out.read_text()) if verify_out.exists() and verify_out.read_text().strip() else []
                    vr = vrows[0] if vrows else None
                    run_row.update({
                        "verify_status": vr.get("status") if vr else "error",
                        "known_subgroup_rank": vr.get("known_subgroup_rank") if vr else None,
                        "exact_gain": vr.get("exact_rank_gain_over_known") if vr else None,
                        "processed_rank": vr.get("processed_subgroup_rank") if vr else None,
                        "verify_seconds": time.monotonic() - started_v,
                        "verify_json": str(verify_out),
                        "verify_log": str(verify_log),
                    })
                    print(
                        f"DEEP|stage=verify|T={parameter}|height={height}|status={run_row['verify_status']}"
                        f"|baseline={run_row['known_subgroup_rank']}|gain={run_row['exact_gain']}"
                        f"|rank={run_row['processed_rank']}",
                        flush=True,
                    )
                except subprocess.TimeoutExpired as exc:
                    text = exc.stdout or ""
                    if isinstance(text, bytes):
                        text = text.decode(errors="replace")
                    verify_log.write_text(text)
                    run_row["verify_status"] = "driver_timeout"
                    print(f"DEEP|stage=verify|T={parameter}|height={height}|status=driver_timeout", flush=True)

            target_summary["runs"].append(run_row)
            if (run_row.get("processed_rank") or 0) > target["processed_subgroup_rank"]:
                print(
                    f"DEEP|BREAKTHROUGH|T={parameter}|height={height}"
                    f"|old_rank={target['processed_subgroup_rank']}|new_rank={run_row['processed_rank']}",
                    flush=True,
                )
                # Continue later heights: additional directions may still exist.

        summary.append(target_summary)
        summary_path = args.output_dir / "summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(f"DEEP|done|targets={len(summary)}|summary={args.output_dir / 'summary.json'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
