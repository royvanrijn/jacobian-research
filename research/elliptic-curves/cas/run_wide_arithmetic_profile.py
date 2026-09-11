#!/usr/bin/env python3
"""Frozen-population equation-only arithmetic profiler for the wide R17 search.

The controller never launches a rational-point search.  It consumes a frozen
search census, normalizes one row per curve, and runs bounded Sage workers for:

  BASE  exact curve / 2-division invariants (fast, no class group),
  LOCAL exact factorization, maximal-order and Brumer--Kramer local terms,
  CLASS optional provisional (proof=False) class-group probes on a frozen panel.

Every timeout/error stays UNKNOWN.  Lower-bound ranks are never treated as exact.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import wide_arithmetic_profile_core as core

SELF = Path(__file__).resolve()
CAS = SELF.parent
ROOT = CAS.parents[1]
DEFAULT_OUTPUT = ROOT / "artifacts/local/elliptic-curves/wide-arithmetic-profile-v1"
WORKER = CAS / "wide_arithmetic_worker.sage"


def require(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return core.sha256_file(Path(path))


def save(path, obj):
    core.write_immutable(Path(path), obj)


def rel(path):
    p = Path(path).resolve()
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def source_hashes():
    return {rel(p): sha(p) for p in (SELF, WORKER, CAS / "wide_arithmetic_profile_core.py")}


def parse_source(source: Path, cap_files: int = 20000):
    files = core.discover_source_files(source)
    require(len(files) <= cap_files, f"too many candidate census files ({len(files)} > {cap_files}); point --source at the frozen campaign folder")
    rows = []
    diagnostics = []
    for p in files:
        try:
            got = core.read_rows_from_file(p)
        except Exception as exc:
            diagnostics.append({"file": str(p), "status": "PARSE_ERROR", "error": str(exc)})
            continue
        if got:
            rows.extend(got)
            diagnostics.append({"file": str(p), "status": "ROWS", "rows": len(got), "sha256": sha(p)})
    return files, rows, diagnostics


def command_probe(args):
    files, rows, diagnostics = parse_source(args.source, args.max_source_files)
    curves = core.aggregate_rows(rows)
    payload = {
        "schema": "elliptic-curves.wide-arithmetic-profile-probe.v1",
        "source": str(args.source.resolve()),
        "candidate_files": len(files),
        "raw_rows": len(rows),
        "unique_curves": len(curves),
        "rank_counts": {},
        "observed_followup_improvements": sum(1 for r in curves if r.get("improved_since_initial") is True),
        "curves_with_bound_initial_stage": sum(1 for r in curves if r.get("initial_rank_lower_bound") is not None),
        "sample_curves": curves[:20],
        "diagnostics_with_rows": [d for d in diagnostics if d["status"] == "ROWS"][:50],
        "parse_errors": [d for d in diagnostics if d["status"] == "PARSE_ERROR"][:50],
    }
    counts = {}
    for row in curves:
        r = str(row["final_rank_lower_bound"])
        counts[r] = counts.get(r, 0) + 1
    payload["rank_counts"] = dict(sorted(counts.items(), key=lambda kv: int(kv[0])))
    print(core.stable_json(payload), end="")
    if args.probe_output:
        Path(args.probe_output).write_text(core.stable_json(payload))
    return 0


def command_prepare(args):
    out = args.output.resolve()
    require(not out.exists(), f"output already exists: {out}; use a fresh folder")
    files, rows, diagnostics = parse_source(args.source, args.max_source_files)
    curves = core.aggregate_rows(rows)
    require(curves, "no curve/rank rows found")
    if args.expected_fibres:
        require(len(curves) == args.expected_fibres, f"expected {args.expected_fibres} curves, found {len(curves)}; inspect `probe` before proceeding")
    if args.expected_tail:
        expected = {}
        for token in args.expected_tail.split(","):
            if not token.strip():
                continue
            r, c = token.split("=", 1)
            expected[int(r.strip())] = int(c.strip())
        actual = {}
        for row in curves:
            r = int(row["final_rank_lower_bound"])
            if r >= min(expected):
                actual[r] = actual.get(r, 0) + 1
        require(actual == expected, f"final-rank tail fingerprint mismatch: expected {expected}, found {actual}; likely wrong campaign source")
    if args.expected_improved >= 0:
        improved = sum(1 for row in curves if row.get("improved_since_initial") is True)
        require(improved == args.expected_improved, f"follow-up improvement fingerprint mismatch: expected {args.expected_improved}, found {improved}; inspect source/history binding")
    out.mkdir(parents=True)
    (out / "inputs").mkdir()
    (out / "base").mkdir()
    (out / "local").mkdir()
    (out / "class").mkdir()
    (out / "logs").mkdir()

    # Persist one canonical row per curve; workers consume these exact bytes.
    for row in curves:
        (out / "inputs" / f"{row['curve_key']}.json").write_text(core.stable_json(row))

    source_manifest = []
    for d in diagnostics:
        if d["status"] == "ROWS":
            source_manifest.append(d)
    panel = core.deterministic_controls(curves, args.high_threshold, args.controls_per_bucket, args.control_salt)
    plan = {
        "schema": "elliptic-curves.wide-arithmetic-profile-plan.v1",
        "created_unix": int(time.time()),
        "source": str(args.source.resolve()),
        "expected_fibres": args.expected_fibres,
        "expected_tail": args.expected_tail,
        "expected_improved": args.expected_improved,
        "population_count": len(curves),
        "raw_rows": len(rows),
        "source_files_with_rows": source_manifest,
        "source_parse_error_count": sum(1 for d in diagnostics if d["status"] == "PARSE_ERROR"),
        "source_parse_errors": [d for d in diagnostics if d["status"] == "PARSE_ERROR"][:100],
        "sources": source_hashes(),
        "population_sha256": hashlib.sha256(core.stable_json(curves).encode()).hexdigest(),
        "frozen_class_panel": panel,
        "boundaries": [
            "All rank values are certified lower bounds, never exact ranks unless separately proved elsewhere.",
            "Follow-up gain is an adaptive search outcome and is not an unbiased rank sample.",
            "BASE/LOCAL inputs use equations and frozen search metadata only; no exceptional points enter arithmetic profiling.",
            "Class probes are optional proof=False / provisional computations and cannot produce unconditional upper bounds.",
            "Timeout/error remains UNKNOWN; no missing arithmetic is imputed as zero.",
        ],
    }
    save(out / "plan.json", plan)
    save(out / "population.json", {"schema": core.SCHEMA, "curves": curves})
    print(f"WIDE_ARITH_PREPARE|status=PASS|curves={len(curves)}|output={out}")
    return 0


def check_plan(out: Path):
    plan = read(out / "plan.json")
    for path, digest in plan["sources"].items():
        p = Path(path)
        if not p.is_absolute():
            p = ROOT / p
        require(p.exists() and sha(p) == digest, f"source changed since prepare: {path}")
    pop = read(out / "population.json")["curves"]
    require(hashlib.sha256(core.stable_json(pop).encode()).hexdigest() == plan["population_sha256"], "population changed")
    return plan, pop


def run_one(out: Path, row: dict, mode: str, timeout: int, memory_gb: float, sage: str):
    target = out / mode / f"{row['curve_key']}.json"
    if target.exists():
        old = read(target)
        if old.get("input_sha256") == sha(out / "inputs" / f"{row['curve_key']}.json"):
            return old.get("status", "UNKNOWN_EXISTING")
        raise RuntimeError(f"existing result input mismatch: {target}")
    temp = target.with_suffix(".tmp.json")
    cmd = [sage, "-python", str(WORKER), "--input", str(out / "inputs" / f"{row['curve_key']}.json"), "--output", str(temp), "--mode", mode, "--memory-gb", str(memory)]
    start = time.monotonic()
    try:
        cp = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            timeout=timeout,
        )
        if temp.exists():
            result = read(temp)
        else:
            result = {
                "schema": f"elliptic-curves.wide-arithmetic-{mode}-controller.v1",
                "curve_key": row["curve_key"], "status": "UNKNOWN_WORKER_NO_OUTPUT",
                "returncode": cp.returncode, "stderr_tail": cp.stderr[-4000:],
            }
    except subprocess.TimeoutExpired as exc:
        result = {
            "schema": f"elliptic-curves.wide-arithmetic-{mode}-controller.v1",
            "curve_key": row["curve_key"], "status": "UNKNOWN_TIMEOUT",
            "timeout_seconds": timeout,
            "stderr_tail": (exc.stderr or "")[-4000:] if isinstance(exc.stderr, str) else "",
        }
    except Exception as exc:
        result = {
            "schema": f"elliptic-curves.wide-arithmetic-{mode}-controller.v1",
            "curve_key": row["curve_key"], "status": "UNKNOWN_CONTROLLER_FAILURE",
            "error_type": type(exc).__name__, "error": str(exc),
        }
    result["input_sha256"] = sha(out / "inputs" / f"{row['curve_key']}.json")
    result["wall_seconds"] = round(time.monotonic() - start, 6)
    temp.unlink(missing_ok=True)
    target.write_text(core.stable_json(result))
    return result["status"]


def run_mode(args, mode: str, selected_keys: set[str] | None = None):
    out = args.output.resolve()
    plan, pop = check_plan(out)
    rows = [r for r in pop if selected_keys is None or r["curve_key"] in selected_keys]
    rows = sorted(rows, key=lambda r: r["curve_key"])
    if getattr(args, "limit", 0):
        rows = rows[:args.limit]
    timeout = {"base": args.base_timeout, "local": args.local_timeout, "class": args.class_timeout}[mode]
    memory = {"base": args.base_memory_gb, "local": args.local_memory_gb, "class": args.class_memory_gb}[mode]
    counts = {}
    done = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_one, out, row, mode, timeout, memory, args.sage): row for row in rows}
        for fut in concurrent.futures.as_completed(futures):
            row = futures[fut]
            try:
                status = fut.result()
            except Exception as exc:
                status = "UNKNOWN_CONTROLLER_EXCEPTION"
                (out / "logs" / f"{mode}-{row['curve_key']}.txt").write_text(repr(exc) + "\n")
            counts[status] = counts.get(status, 0) + 1
            done += 1
            if done % 25 == 0 or done == len(rows):
                print(f"WIDE_ARITH_{mode.upper()}|done={done}|total={len(rows)}|status_counts={json.dumps(counts,sort_keys=True)}", flush=True)
    return counts


def command_run(args):
    out = args.output.resolve()
    check_plan(out)
    run_mode(args, "base")
    run_mode(args, "local")
    if args.with_class_probe:
        plan = read(out / "plan.json")
        selected = set(plan["frozen_class_panel"]["high"] + plan["frozen_class_panel"]["controls"])
        run_mode(args, "class", selected)
    summarize(out)
    print("WIDE_ARITH_RUN|status=PASS_CONTROLLER_COMPLETED|note=individual UNKNOWN results remain explicit")
    return 0


def command_class_probe(args):
    out = args.output.resolve()
    plan, _pop = check_plan(out)
    selected = set(plan["frozen_class_panel"]["high"] + plan["frozen_class_panel"]["controls"])
    run_mode(args, "class", selected)
    summarize(out)
    return 0


def load_result(path):
    return read(path) if path.exists() else None


def summarize(out: Path):
    plan, pop = check_plan(out)
    flat = []
    for row in pop:
        base = load_result(out / "base" / f"{row['curve_key']}.json")
        local = load_result(out / "local" / f"{row['curve_key']}.json")
        cls = load_result(out / "class" / f"{row['curve_key']}.json")
        flat.append(core.flatten_profile(row, base, local, cls))
    summary = core.summarize_profiles(flat)
    summary["stage_status_counts"] = {}
    for mode in ("base", "local", "class"):
        counts = {}
        for p in (out / mode).glob("*.json"):
            st = read(p).get("status", "UNKNOWN")
            counts[st] = counts.get(st, 0) + 1
        summary["stage_status_counts"][mode] = counts
    (out / "profiles.json").write_text(core.stable_json({"schema": core.SCHEMA, "rows": flat}))
    fields = sorted({k for r in flat for k in r})
    with (out / "profiles.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(flat)
    (out / "summary.json").write_text(core.stable_json(summary))
    (out / "SUMMARY.md").write_text(core.markdown_summary(summary, flat))
    base_count = sum(1 for _ in (out / "base").glob("*.json"))
    local_count = sum(1 for _ in (out / "local").glob("*.json"))
    full = base_count == len(pop) and local_count == len(pop)
    report = {
        "schema": "elliptic-curves.wide-arithmetic-profile-report.v1",
        "status": "PASS_PROFILE_SUMMARY_FULL" if full else "PASS_PARTIAL_PROFILE_SUMMARY",
        "population_count": len(pop),
        "population_sha256": plan["population_sha256"],
        "profiles_sha256": sha(out / "profiles.json"),
        "summary_sha256": sha(out / "summary.json"),
        "markdown_sha256": sha(out / "SUMMARY.md"),
        "stage_status_counts": summary["stage_status_counts"],
        "base_checkpoint_count": base_count,
        "local_checkpoint_count": local_count,
        "full_population_checkpointed": full,
    }
    (out / "REPORT.json").write_text(core.stable_json(report))
    return report


def command_summarize(args):
    report = summarize(args.output.resolve())
    print(core.stable_json(report), end="")
    return 0


def command_status(args):
    out = args.output.resolve()
    if not (out / "plan.json").exists():
        print("WIDE_ARITH_STATUS|state=NOT_PREPARED")
        return 1
    plan = read(out / "plan.json")
    parts = []
    for mode in ("base", "local", "class"):
        counts = {}
        for p in (out / mode).glob("*.json"):
            st = read(p).get("status", "UNKNOWN")
            counts[st] = counts.get(st, 0) + 1
        parts.append(f"{mode}={sum(counts.values())}/{plan['population_count']}:{json.dumps(counts,sort_keys=True)}")
    print("WIDE_ARITH_STATUS|" + "|".join(parts))
    return 0


def command_check(args):
    out = args.output.resolve()
    plan, pop = check_plan(out)
    report = read(out / "REPORT.json")
    require(report["population_count"] == len(pop), "report population mismatch")
    require(report.get("full_population_checkpointed") is True, "profile is partial; BASE and LOCAL must each checkpoint all population rows before final check")
    require(report["population_sha256"] == plan["population_sha256"], "report population hash mismatch")
    # Recompute summary from immutable per-curve checkpoints and compare bytes.
    before = {name: (out / name).read_bytes() for name in ("profiles.json", "summary.json", "SUMMARY.md", "REPORT.json")}
    summarize(out)
    after = {name: (out / name).read_bytes() for name in before}
    require(before == after, "deterministic summary replay mismatch")
    print("WIDE_ARITH_CHECK|status=PASS|deterministic_summary=PASS")
    return 0


def add_runtime_args(ap):
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--sage", default=os.environ.get("SAGE", "sage"))
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0, help="smoke only: process first N deterministic rows; 0 means all")
    ap.add_argument("--base-timeout", type=int, default=60)
    ap.add_argument("--local-timeout", type=int, default=180)
    ap.add_argument("--class-timeout", type=int, default=120)
    ap.add_argument("--base-memory-gb", type=float, default=0.0, help="optional RLIMIT_AS; 0 disables because Sage virtual-memory maps vary by host")
    ap.add_argument("--local-memory-gb", type=float, default=0.0, help="optional RLIMIT_AS; 0 disables")
    ap.add_argument("--class-memory-gb", type=float, default=0.0, help="optional RLIMIT_AS; 0 disables")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("probe", help="inspect a frozen campaign source without writing a plan")
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--max-source-files", type=int, default=20000)
    p.add_argument("--probe-output", type=Path)
    p.set_defaults(func=command_probe)

    p = sub.add_parser("prepare", help="freeze and normalize the campaign population")
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--expected-fibres", type=int, default=2080)
    p.add_argument("--expected-tail", default="23=25,24=4,25=1", help="exact final lower-bound tail fingerprint; empty string disables")
    p.add_argument("--expected-improved", type=int, default=-1, help="optional campaign-history gate; use 157 only after probe confirms pre-follow-up ranks")
    p.add_argument("--max-source-files", type=int, default=20000)
    p.add_argument("--high-threshold", type=int, default=23)
    p.add_argument("--controls-per-bucket", type=int, default=30)
    p.add_argument("--control-salt", default="wide-arithmetic-profile-v1")
    p.set_defaults(func=command_prepare)

    for name, func in (("run", command_run), ("class-probe", command_class_probe)):
        p = sub.add_parser(name)
        add_runtime_args(p)
        if name == "run":
            p.add_argument("--with-class-probe", action="store_true")
        p.set_defaults(func=func)
    p = sub.add_parser("summarize"); add_runtime_args(p); p.set_defaults(func=command_summarize)
    p = sub.add_parser("status"); add_runtime_args(p); p.set_defaults(func=command_status)
    p = sub.add_parser("check"); add_runtime_args(p); p.set_defaults(func=command_check)
    args = ap.parse_args()
    if hasattr(args, "jobs"):
        require(args.jobs >= 1, "jobs must be positive")
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
