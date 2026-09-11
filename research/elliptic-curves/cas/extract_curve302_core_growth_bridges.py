#!/usr/bin/env python3
"""Extract the 27 decisive Curve302 common-core bridge stages.

Usage:
  python3 research/elliptic-curves/cas/extract_curve302_core_growth_bridges.py run [--replay DIR] [--output DIR] [--raw-root DIR] [--skip-raw-schema]
  python3 research/elliptic-curves/cas/extract_curve302_core_growth_bridges.py check [--replay DIR] [--output DIR] [--raw-root DIR] [--skip-raw-schema]

The arithmetic tables use only the sealed replay ledger.  Optional raw-schema
enrichment reopens only the historical chart definitions for the decisive stages
and records centre/mapping/search descriptors; it does not rerun point search or
point recognition.
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

from curve302_short_core_controls import require
from curve302_positive_exposure_tables import exposure_index
from curve302_core_growth_bridges import analyze_core_growth, object_shape
from extract_curve302_positive_exposure_tables import load

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
DEFAULT_REPLAY = ROOT / "artifacts/local/elliptic-curves/curve302-chart-replay-v1"


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def hash_obj(obj):
    return sha256(canonical(obj).encode()).hexdigest()


def atomic(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, sort_keys=True, indent=2, allow_nan=False)
            f.write("\n")
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        Path(path).write_text("")
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with Path(path).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list, tuple)) else v
                        for k, v in row.items()})


def md_table(rows, columns):
    def text(v):
        if v is None:
            return "—"
        if isinstance(v, bool):
            return "yes" if v else "no"
        return str(v)
    out = ["| " + " | ".join(label for _, label in columns) + " |",
           "|" + "|".join("---" for _ in columns) + "|"]
    for row in rows:
        out.append("| " + " | ".join(text(row.get(k)) for k, _ in columns) + " |")
    return "\n".join(out)


def discover_raw_root(replay: Path, explicit: Path | None):
    if explicit is not None:
        p = Path(explicit).resolve()
        require(p.is_dir(), f"raw root missing: {p}")
        return p
    # The downstream exposure experiment was launched with ``--ledger``, so
    # its plan records the replay output directory rather than the historical
    # transcript root.  Recover the latter from the exact MW-state provenance
    # retained by the replay adapter.
    replay_ledger = replay / "chart-exposure-ledger.json"
    if replay_ledger.is_file():
        chosen = (read(replay_ledger).get("adapter", {}).get("basis", {})
                  .get("chosen_states", {}))
        source_files = [str(row["file"]) for row in chosen.values()
                        if isinstance(row, dict) and row.get("file")]
        if source_files:
            common = Path(os.path.commonpath(source_files)).resolve()
            if common.is_dir():
                return common
    plan = replay / "experiments" / "plan.json"
    if plan.is_file():
        raw = read(plan).get("raw_root")
        if raw and Path(raw).is_dir():
            return Path(raw).resolve()
    return None


def raw_schema_enrichment(raw_root: Path, analysis, ledger):
    """Load raw chart definitions only for the decisive seed/epoch pairs."""
    from curve302_chart_replay_adapter import load_stage_charts

    wanted = sorted({(row["seed"], int(row["epoch"])) for row in analysis["stages"]})
    expected = {}
    for seed, epoch in wanted:
        expected.setdefault(seed, set()).add(epoch)
    raw = load_stage_charts(raw_root, tuple(ledger["direction_ids"]), expected)

    decisive_chart_ids = set()
    for row in analysis["bridge_candidates"]:
        decisive_chart_ids.update(row.get("chart_ids", ()))

    descriptors = []
    index = {}
    for (seed, epoch), charts in raw.items():
        for chart in charts:
            if chart["chart_id"] not in decisive_chart_ids:
                continue
            mapping = chart.get("mapping")
            search = chart.get("search")
            desc = {
                "seed": seed,
                "epoch": int(epoch),
                "chart_id": chart["chart_id"],
                "order": int(chart["order"]),
                "index": int(chart["index"]),
                "centre": chart.get("centre"),
                "mapping": mapping,
                "search": search,
                "mapping_schema": object_shape(mapping),
                "search_schema": object_shape(search),
                "mapping_schema_hash": hash_obj(object_shape(mapping)),
                "search_schema_hash": hash_obj(object_shape(search)),
                "mapping_value_hash": hash_obj(mapping),
                "search_value_hash": hash_obj(search),
                "source_file": Path(chart["source_file"]).name,
            }
            descriptors.append(desc)
            index[(seed, int(epoch), chart["chart_id"])] = desc

    schema_groups = {}
    for desc in descriptors:
        key = (desc["mapping_schema_hash"], desc["search_schema_hash"])
        row = schema_groups.setdefault(key, {
            "mapping_schema_hash": key[0],
            "search_schema_hash": key[1],
            "charts": 0,
            "stages": set(),
            "actual_bridge_exposures": 0,
            "alternative_tied_bridge_exposures": 0,
        })
        row["charts"] += 1
        row["stages"].add((desc["seed"], desc["epoch"]))

    for bridge in analysis["bridge_candidates"]:
        role_actual = bool(bridge["actual_member"])
        for cid in bridge.get("chart_ids", ()):
            desc = index.get((bridge["seed"], int(bridge["epoch"]), cid))
            if desc is None:
                continue
            key = (desc["mapping_schema_hash"], desc["search_schema_hash"])
            if role_actual:
                schema_groups[key]["actual_bridge_exposures"] += 1
            elif bridge["tied_best_positive"]:
                schema_groups[key]["alternative_tied_bridge_exposures"] += 1

    clusters = []
    for row in schema_groups.values():
        row = dict(row)
        row["stages"] = len(row["stages"])
        clusters.append(row)
    clusters.sort(key=lambda r: (-r["stages"], -r["charts"], r["mapping_schema_hash"], r["search_schema_hash"]))
    return {
        "status": "PASS_RAW_DECISIVE_CHART_SCHEMA_ENRICHMENT",
        "raw_root": str(raw_root),
        "decisive_chart_descriptors": descriptors,
        "schema_clusters": clusters,
        "boundary": "Raw centre/mapping/search descriptors are read only for charts positively exposing actual or tied-best bridge directions in the decisive stages.",
    }


def summary_md(analysis, schema):
    s = analysis["summary"]
    stage_rows = analysis["stages"]
    by_class = []
    for name in ("ACTUAL_UNIQUE_BEST", "ACTUAL_TIED_WITH_ALTERNATIVES", "ACTUAL_NOT_BEST"):
        by_class.append({"classification": name, "stages": sum(r["classification"] == name for r in stage_rows)})
    top_axis = analysis["clusters"]["new_axis_signature"][:12]
    top_core = analysis["clusters"]["core_bridge_line"][:12]
    lines = [
        "# Curve302 decisive core-growth bridges",
        "",
        "Positive chart evidence only. Missing historical hits are never interpreted as non-exposure.",
        "",
        f"Decisive one-step common-core growth stages: **{s['decisive_core_growth_stages']}**.",
        f"Positive candidate-stage pairs inside those stages: **{s['candidate_stage_pairs']}**.",
        "",
        "## Actual versus positive alternatives",
        "",
        md_table(by_class, [("classification", "Classification"), ("stages", "Stages")]),
        "",
        "A stage is `ACTUAL_UNIQUE_BEST` only when no explicitly positive alternative matches the actual gain's exact next-core rank increment. `ACTUAL_TIED_WITH_ALTERNATIVES` means at least one positive alternative produces the same best one-step core gain after saturation.",
        "",
        "## The decisive stages",
        "",
        md_table(stage_rows, [
            ("seed", "Seed"), ("epoch", "Epoch"), ("pre_dimension", "Pre dim"),
            ("post_dimension", "Post dim"), ("current_next_core_intersection_rank", "Core before"),
            ("actual_batch_next_core_delta", "Actual batch Δ"), ("best_actual_single_delta", "Best actual Δ"),
            ("best_alternative_single_delta", "Best alt Δ"), ("tied_best_alternative_count", "Tied alts"),
            ("classification", "Class"),
        ]),
        "",
        "## Most common newly enabled named-axis signatures",
        "",
        md_table(top_axis, [("key", "New axes"), ("stages", "Stages"), ("occurrences", "Occurrences"),
                            ("actual_occurrences", "Actual"), ("alternative_tied_best_occurrences", "Alt tied")]),
        "",
        "## Most common new core bridge lines",
        "",
        md_table(top_core, [("key", "Core bridge line"), ("stages", "Stages"), ("occurrences", "Occurrences"),
                            ("actual_occurrences", "Actual"), ("alternative_tied_best_occurrences", "Alt tied")]),
        "",
    ]
    if schema is not None:
        lines += [
            "## Raw chart schema clusters",
            "",
            f"Raw decisive chart descriptors: **{len(schema['decisive_chart_descriptors'])}**.",
            "",
            md_table(schema["schema_clusters"][:12], [
                ("mapping_schema_hash", "Mapping schema"), ("search_schema_hash", "Search schema"),
                ("charts", "Charts"), ("stages", "Stages"),
                ("actual_bridge_exposures", "Actual exposures"),
                ("alternative_tied_bridge_exposures", "Alt tied exposures"),
            ]),
            "",
        ]
    else:
        lines += [
            "## Raw chart schema clusters",
            "",
            "Not generated: no accessible raw-root was supplied/discovered, or `--skip-raw-schema` was used. Arithmetic bridge tables are unaffected.",
            "",
        ]
    return "\n".join(lines)


def compute(args):
    replay, exp, data, ledger, ledger_path, stats = load(args.replay)
    analysis = analyze_core_growth(data, ledger, 14)
    require(analysis["summary"]["decisive_core_growth_stages"] > 0, "no decisive core-growth stages found")
    # The current sealed result is known to have 27; fail if the source changes
    # unless the user explicitly opts out of this calibration gate.
    if not args.allow_stage_count_change:
        require(analysis["summary"]["decisive_core_growth_stages"] == 27,
                f"expected 27 decisive stages, got {analysis['summary']['decisive_core_growth_stages']}; pass --allow-stage-count-change only for an intentional new source")
    schema = None
    raw_root = None if args.skip_raw_schema else discover_raw_root(replay, args.raw_root)
    if raw_root is not None:
        schema = raw_schema_enrichment(raw_root, analysis, ledger)
    return replay, data, ledger_path, stats, analysis, schema


def run(args):
    replay, data, ledger_path, stats, analysis, schema = compute(args)
    out = Path(args.output or replay / "core-growth-bridges-v1").resolve()
    require(not out.exists(), "output exists; use a fresh --output or run check")
    out.mkdir(parents=True)
    atomic(out / "core-growth-analysis.json", analysis)
    if schema is not None:
        atomic(out / "raw-chart-schema.json", schema)
    write_csv(out / "core-growth-stages.csv", analysis["stages"])
    write_csv(out / "core-growth-candidates.csv", analysis["candidates"])
    write_csv(out / "bridge-candidates.csv", analysis["bridge_candidates"])
    for name, rows in analysis["clusters"].items():
        write_csv(out / f"cluster-{name}.csv", rows)
    if schema is not None:
        write_csv(out / "chart-schema-clusters.csv", schema["schema_clusters"])
    summary = summary_md(analysis, schema)
    (out / "SUMMARY.md").write_text(summary)
    report = {
        "status": "PASS_CURVE302_DECISIVE_CORE_GROWTH_BRIDGES",
        "frozen_ledger_sha256": sha(ledger_path),
        "ledger_stats": stats,
        "decisive_core_growth_stages": analysis["summary"]["decisive_core_growth_stages"],
        "analysis_sha256": sha(out / "core-growth-analysis.json"),
        "raw_schema_sha256": sha(out / "raw-chart-schema.json") if schema is not None else None,
        "summary_sha256": sha(out / "SUMMARY.md"),
        "boundary": "Positive-evidence-only exact lattice analysis. No missing chart hit is treated as a negative; raw chart schema enrichment is descriptive only.",
    }
    atomic(out / "REPORT.json", report)
    print(summary, end="")
    print(f"CURVE302_CORE_GROWTH_BRIDGES|status=PASS|output={out}")


def check(args):
    replay, data, ledger_path, stats, analysis, schema = compute(args)
    out = Path(args.output or replay / "core-growth-bridges-v1").resolve()
    report = read(out / "REPORT.json")
    require(report.get("status") == "PASS_CURVE302_DECISIVE_CORE_GROWTH_BRIDGES", "report not passed")
    require(report["frozen_ledger_sha256"] == sha(ledger_path), "frozen replay ledger changed")
    analysis_bytes = json.dumps(analysis, sort_keys=True, indent=2, allow_nan=False) + "\n"
    require((out / "core-growth-analysis.json").read_text() == analysis_bytes,
            "core-growth analysis deterministic mismatch")
    require(report["analysis_sha256"] == sha(out / "core-growth-analysis.json"), "analysis hash mismatch")
    if schema is not None:
        schema_bytes = json.dumps(schema, sort_keys=True, indent=2, allow_nan=False) + "\n"
        require((out / "raw-chart-schema.json").read_text() == schema_bytes,
                "raw chart schema deterministic mismatch")
        require(report["raw_schema_sha256"] == sha(out / "raw-chart-schema.json"), "raw schema hash mismatch")
    else:
        require(report["raw_schema_sha256"] is None, "report expected raw schema but recomputation has none")
    require((out / "SUMMARY.md").read_text() == summary_md(analysis, schema), "SUMMARY deterministic mismatch")
    require(report["summary_sha256"] == sha(out / "SUMMARY.md"), "summary hash mismatch")
    print("CURVE302_CORE_GROWTH_BRIDGES_CHECK|status=PASS")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("command", choices=("run", "check"))
    p.add_argument("--replay", type=Path, default=DEFAULT_REPLAY)
    p.add_argument("--output", type=Path)
    p.add_argument("--raw-root", type=Path)
    p.add_argument("--skip-raw-schema", action="store_true")
    p.add_argument("--allow-stage-count-change", action="store_true")
    args = p.parse_args()
    run(args) if args.command == "run" else check(args)


if __name__ == "__main__":
    main()
