#!/usr/bin/env python3
"""Extract Curve302's exact closure transition graph and rank-29 control.

Usage:
  python3 research/elliptic-curves/cas/extract_curve302_closure_transition_graph.py run [--replay DIR] [--bridges DIR] [--output DIR]
  python3 research/elliptic-curves/cas/extract_curve302_closure_transition_graph.py check [--replay DIR] [--bridges DIR] [--output DIR]
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import tempfile

from curve302_short_core_controls import require
from curve302_closure_transition_graph import transition_graph, rank29_control
from extract_curve302_positive_exposure_tables import load

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
DEFAULT_REPLAY = ROOT / "artifacts/local/elliptic-curves/curve302-chart-replay-v1"


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def atomic(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, sort_keys=True, indent=2, allow_nan=False); f.write("\n")
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        Path(path).write_text(""); return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with Path(path).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list, tuple)) else v for k, v in r.items()})


def discover_bridges(replay: Path, explicit: Path | None):
    if explicit:
        p = Path(explicit).resolve(); require((p / "REPORT.json").is_file(), "bridge REPORT missing"); return p
    candidates = []
    for p in replay.glob("core-growth-bridges-v*"):
        report = p / "REPORT.json"
        if not report.is_file(): continue
        try: obj = read(report)
        except Exception: continue
        if obj.get("status") != "PASS_CURVE302_DECISIVE_CORE_GROWTH_BRIDGES": continue
        m = re.search(r"-v(\d+)$", p.name)
        version = int(m.group(1)) if m else -1
        candidates.append((version, p.name, p))
    require(candidates, "no passed core-growth-bridges-v* bundle found; use --bridges")
    candidates.sort()
    return candidates[-1][2].resolve()


def md_table(rows, cols):
    def t(v):
        if v is None: return "—"
        if isinstance(v, bool): return "yes" if v else "no"
        return str(v)
    out = ["| " + " | ".join(label for _, label in cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for r in rows:
        out.append("| " + " | ".join(t(r.get(k)) for k, _ in cols) + " |")
    return "\n".join(out)


def summary_md(graph, stall):
    s = graph["summary"]
    bottlenecks = sorted(graph["core_transition_types"], key=lambda r: (r.get("observed_bridge_multiplicity", 999), -r["stages"], r["core_edge_id"]))
    ub = graph["unique_best_stages"]
    fam = stall["successful_late_bridge_families_mod_stall"][:12]
    lines = [
        "# Curve302 closure transition graph",
        "",
        "Exact post-processing of the sealed positive bridge analysis. Missing chart hits are never interpreted as negatives.",
        "",
        "## Graph census",
        "",
        f"- decisive stages: **{s['decisive_stages']}**",
        f"- positive actual/tied bridge pairs: **{s['bridge_candidate_pairs']}**",
        f"- literal saturated state edges: **{s['literal_saturated_edges']}**",
        f"- exact common-core transition types: **{s['common_core_transition_types']}**",
        f"- narrow observed gates (one bridge line): **{s['narrow_observed_core_gates']}**",
        f"- broad observed basins (>1 bridge line): **{s['broad_observed_core_basins']}**",
        "",
        "## Observed common-core transitions",
        "",
        md_table(bottlenecks, [
            ("target_post_dimension", "Post dim"), ("core_before_rank", "Core before"),
            ("core_after_rank", "Core after"), ("stages", "Stages"),
            ("observed_bridge_multiplicity", "Bridge lines"), ("actual_occurrences", "Actual"),
            ("alternative_tied_occurrences", "Alt tied"), ("bottleneck_class", "Class"),
        ]),
        "",
        "## The six unique-best historical stages",
        "",
        md_table(ub, [
            ("seed", "Seed"), ("epoch", "Epoch"), ("pre_dimension", "Pre dim"),
            ("post_dimension", "Post dim"), ("core_before_rank", "Core before"), ("core_delta", "Δ core"),
        ]),
        "",
        "## Rank-29 exception",
        "",
        f"Stalled seed: **{stall['stalled_seed']}**; terminal prefix rank **{stall['terminal_prefix_rank']}**.",
        f"Successful dimension-13 common core deficit from the stalled state: **{stall['dimension13_deficit']}**.",
        f"Full dimension-14 deficit from the stalled state: **{stall['dimension14_deficit']}**.",
        f"Recorded terminal rationally-new positive directions: **{stall['terminal_positive_new_directions']}**.",
        f"Recorded terminal C13 improvers: **{stall['terminal_positive_c13_improvers']}**.",
        f"Successful late bridge-line families modulo the stalled prefix: **{stall['successful_late_bridge_line_families_mod_stall']}**.",
        f"Such families positively seen at the stalled terminal stage: **{stall['successful_line_families_seen_positive_at_terminal']}**.",
        f"First late-gate status: **{stall['first_missing_gate_status']}**.",
        "",
        "### Successful late bridge families projected into the stalled two-dimensional quotient",
        "",
        md_table(fam, [
            ("quotient_line_mod_stall", "Line mod stall"), ("successful_stages", "Stages"),
            ("successful_occurrences", "Occurrences"), ("actual_occurrences", "Actual"),
            ("max_c13_delta_if_added_to_stall", "C13 Δ"), ("max_c14_delta_if_added_to_stall", "C14 Δ"),
            ("terminal_positive_match", "Seen terminal +"),
        ]),
        "",
        "The rank-29 comparison is positive-evidence-only: `NO_POSITIVE_TERMINAL_*` means no matching/improving bridge was recorded positively, not that exhaustive chart absence has been proved.",
        "",
    ]
    return "\n".join(lines)


def compute(args):
    replay, exp, data, ledger, ledger_path, stats = load(args.replay)
    bridges = discover_bridges(replay, args.bridges)
    br = read(bridges / "core-growth-analysis.json")
    require(br.get("status") == "PASS_27_CORE_GROWTH_BRIDGE_ANALYSIS", "bridge analysis not passed")
    require(br["summary"]["decisive_core_growth_stages"] == 27, "expected 27 decisive stages")
    require(br["summary"]["actual_tied_with_alternatives_stages"] == 21, "expected 21 tied stages")
    require(br["summary"]["actual_unique_best_stages"] == 6, "expected 6 unique-best stages")
    require(br["summary"]["actual_not_best_stages"] == 0, "expected zero not-best stages")
    graph = transition_graph(data, ledger, br, 14)
    stall = rank29_control(data, ledger, br, graph, 14)
    return replay, bridges, ledger_path, graph, stall


def run(args):
    replay, bridges, ledger_path, graph, stall = compute(args)
    out = Path(args.output or replay / "closure-transition-graph-v1").resolve()
    require(not out.exists(), "output exists; choose a fresh --output or use check")
    out.mkdir(parents=True)
    atomic(out / "transition-graph.json", graph)
    atomic(out / "rank29-comparison.json", stall)
    write_csv(out / "core-transition-types.csv", graph["core_transition_types"])
    write_csv(out / "unique-best-stages.csv", graph["unique_best_stages"])
    write_csv(out / "bridge-edges.csv", graph["bridge_edges"])
    write_csv(out / "rank29-late-bridge-families.csv", stall["successful_late_bridge_families_mod_stall"])
    (out / "SUMMARY.md").write_text(summary_md(graph, stall))
    report = {
        "status": "PASS_CURVE302_CLOSURE_TRANSITION_GRAPH_AND_RANK29_CONTROL",
        "bridge_report_sha256": sha(bridges / "REPORT.json"),
        "bridge_analysis_sha256": sha(bridges / "core-growth-analysis.json"),
        "frozen_ledger_sha256": sha(ledger_path),
        "outputs": {
            "transition-graph.json": sha(out / "transition-graph.json"),
            "rank29-comparison.json": sha(out / "rank29-comparison.json"),
            "SUMMARY.md": sha(out / "SUMMARY.md"),
        },
        "boundary": "Exact positive-evidence transition collapse plus rank-29 control; no chart-completeness assumption and no new point search.",
    }
    atomic(out / "REPORT.json", report)
    print((out / "SUMMARY.md").read_text(), end="")
    print(f"CURVE302_TRANSITION_GRAPH|status=PASS|output={out}")


def check(args):
    replay, bridges, ledger_path, graph, stall = compute(args)
    out = Path(args.output or replay / "closure-transition-graph-v1").resolve()
    report = read(out / "REPORT.json")
    require(report.get("status") == "PASS_CURVE302_CLOSURE_TRANSITION_GRAPH_AND_RANK29_CONTROL", "report not passed")
    require(report["bridge_report_sha256"] == sha(bridges / "REPORT.json"), "bridge report changed")
    require(report["bridge_analysis_sha256"] == sha(bridges / "core-growth-analysis.json"), "bridge analysis changed")
    require(report["frozen_ledger_sha256"] == sha(ledger_path), "ledger changed")
    require(read(out / "transition-graph.json") == graph, "transition graph deterministic mismatch")
    require(read(out / "rank29-comparison.json") == stall, "rank29 comparison deterministic mismatch")
    require((out / "SUMMARY.md").read_text() == summary_md(graph, stall), "summary deterministic mismatch")
    for name in ("transition-graph.json", "rank29-comparison.json", "SUMMARY.md"):
        require(report["outputs"][name] == sha(out / name), f"{name} hash mismatch")
    print("CURVE302_TRANSITION_GRAPH_CHECK|status=PASS")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("command", choices=("run", "check"))
    p.add_argument("--replay", type=Path, default=DEFAULT_REPLAY)
    p.add_argument("--bridges", type=Path)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    run(a) if a.command == "run" else check(a)


if __name__ == "__main__":
    main()
