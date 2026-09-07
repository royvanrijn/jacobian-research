#!/usr/bin/env python3
"""Summarize and independently audit the finite 302 subgroup landscape."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
LANDSCAPE = ART / "curve302_exceptional_subgroup_landscape_v1.json"
OUTPUT = ART / "curve302_exceptional_subgroup_landscape_report_v1.json"
TSV = ART / "curve302_exceptional_subgroup_landscape_v1.tsv"
SVG = ART / "curve302_exceptional_subgroup_landscape_v1.svg"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def qheight(numerator: int) -> float:
    return numerator / 4_000_000.0


def order_metrics(order, start_mask, retained):
    mask = start_mask
    scores = []
    for target in order:
        score = retained[mask][target]
        if score is None:
            raise ArithmeticError("order attempted a direction already in the subgroup")
        scores.append(score)
        mask |= 1 << target
    return {"scores": scores, "maximum_numerator": max(scores, default=0), "sum_numerator": sum(scores), "final_mask": mask}


def optimize_from(start_mask, allowed, retained):
    """Exact DP over a declared remaining-coordinate cube."""

    state = {0: (0, 0, [])}
    for small_mask in range(1 << len(allowed)):
        bottleneck, total, path = state[small_mask]
        full_mask = start_mask | sum((1 << allowed[index]) for index in range(len(allowed)) if (small_mask >> index) & 1)
        for local, target in enumerate(allowed):
            if (small_mask >> local) & 1:
                continue
            score = retained[full_mask][target]
            next_mask = small_mask | (1 << local)
            candidate = (max(bottleneck, score), total + score, path + [target])
            old = state.get(next_mask)
            if old is None or candidate[:2] < old[:2] or (candidate[:2] == old[:2] and candidate[2] < old[2]):
                state[next_mask] = candidate
    best_bottleneck = min(state[(1 << len(allowed)) - 1] for _ in [0])
    # A separate total optimum is needed because the tuple above is lexicographic
    # in bottleneck then total.
    all_paths = []
    for order in itertools.permutations(allowed):
        metrics = order_metrics(order, start_mask, retained)
        all_paths.append((metrics["maximum_numerator"], metrics["sum_numerator"], list(order)))
    return {
        "minimax": {"maximum_numerator": best_bottleneck[0], "sum_numerator": best_bottleneck[1], "order_indices": best_bottleneck[2]},
        "minimum_total": {
            "maximum_numerator": min(all_paths, key=lambda row: (row[1], row[0], row[2]))[0],
            "sum_numerator": min(all_paths, key=lambda row: (row[1], row[0], row[2]))[1],
            "order_indices": min(all_paths, key=lambda row: (row[1], row[0], row[2]))[2],
        },
    }


def monotonicity_audit(fresh, retained, directions):
    raw_increases = 0
    raw_edges = 0
    largest = None
    retained_edges = 0
    for mask in range(len(fresh)):
        for direction in range(directions):
            if (mask >> direction) & 1:
                continue
            for added in range(directions):
                if added == direction or (mask >> added) & 1:
                    continue
                child = mask | (1 << added)
                old, new = fresh[mask][direction], fresh[child][direction]
                if old is None or new is None:
                    raise ArithmeticError("missing stage-local cost")
                raw_edges += 1
                if new > old:
                    raw_increases += 1
                    candidate = (new / old, mask, added, direction, old, new)
                    if largest is None or candidate > largest:
                        largest = candidate
                held_old, held_new = retained[mask][direction], retained[child][direction]
                retained_edges += 1
                if held_new > held_old:
                    raise ArithmeticError("retained finite atlas lost an admissible witness")
    return {
        "raw_stage_local_edges": raw_edges,
        "raw_stage_local_increases": raw_increases,
        "retained_edges_checked": retained_edges,
        "retained_increases": 0,
        "largest_raw_increase": None if largest is None else {
            "ratio": f"{largest[0]:.12g}", "state_mask": largest[1], "adjoin_direction_index": largest[2],
            "measured_direction_index": largest[3], "old_numerator": largest[4], "new_numerator": largest[5],
        },
    }


def plot_svg(names, tail_rows, optimal_rows):
    """Small static figure: raw and retained costs along two actual paths."""

    width, height, left, top, panel_h = 1024, 520, 74, 42, 190
    series = [("Historical M24→M31 tail", tail_rows), ("Global minimax diagnostic order", optimal_rows)]
    values = [qheight(row[key]) for _, rows in series for row in rows for key in ("stage_local_numerator", "retained_numerator")]
    lo, hi = min(values), max(values)
    lo = max(lo * 0.9, 1.0e-10); hi *= 1.1
    def xy(index, value, panel):
        x = left + index * (width - left - 30) / max(1, len(series[panel][1]) - 1)
        y0 = top + panel * (panel_h + 64)
        y = y0 + panel_h * (1 - (math.log10(value) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)))
        return x, y
    fragments = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
                 '<rect width="100%" height="100%" fill="white"/>',
                 '<style>text{font-family:system-ui,sans-serif;fill:#1f2937} .axis{stroke:#64748b;stroke-width:1} .raw{fill:none;stroke:#d97706;stroke-width:2} .held{fill:none;stroke:#0f766e;stroke-width:2} .dotraw{fill:#d97706}.dotheld{fill:#0f766e}</style>',
                 '<text x="74" y="20" font-size="16" font-weight="600">302 half-lattice score: stage-local versus retained finite atlas</text>',
                 '<text x="760" y="20" font-size="12" fill="#0f766e">teal = retained C*, amber = local C</text>']
    for panel, (title, rows) in enumerate(series):
        y0 = top + panel * (panel_h + 64)
        fragments += [f'<rect x="{left}" y="{y0}" width="{width-left-30}" height="{panel_h}" fill="none" stroke="#94a3b8"/>',
                      f'<text x="{left}" y="{y0-10}" font-size="13" font-weight="600">{title}</text>',
                      f'<text x="16" y="{y0+panel_h/2}" transform="rotate(-90 16 {y0+panel_h/2})" font-size="11">log10(C/4·10⁶)</text>']
        for label, key, cls, dot in (("local", "stage_local_numerator", "raw", "dotraw"), ("retained", "retained_numerator", "held", "dotheld")):
            points = " ".join(f"{x:.2f},{y:.2f}" for index, row in enumerate(rows) for x, y in [xy(index, qheight(row[key]), panel)])
            fragments.append(f'<polyline class="{cls}" points="{points}"/>')
            for index, row in enumerate(rows):
                x, y = xy(index, qheight(row[key]), panel)
                fragments.append(f'<circle class="{dot}" cx="{x:.2f}" cy="{y:.2f}" r="3"/><text x="{x:.2f}" y="{y0+panel_h+17}" text-anchor="middle" font-size="9">{row["add_direction"].replace("residual-strict-", "s").replace("recovered-local-", "l").replace("recovered-strict-", "r")}</text>')
        for tick in (lo, math.sqrt(lo * hi), hi):
            _, y = xy(0, tick, panel)
            fragments.append(f'<line class="axis" x1="{left}" x2="{width-30}" y1="{y:.2f}" y2="{y:.2f}" opacity=".25"/><text x="{left-8}" y="{y+4:.2f}" text-anchor="end" font-size="10">{math.log10(tick):.1f}</text>')
    fragments.append('</svg>')
    return "\n".join(fragments) + "\n"


def build() -> tuple[dict, str, str]:
    landscape = read(LANDSCAPE)
    if landscape["status"] != "PASS_RETROSPECTIVE_FINITE_SUBGROUP_LANDSCAPE":
        raise ArithmeticError("the complete 2^14 landscape is unavailable")
    names = landscape["fixed_basis"]["direction_ids"]
    directions = len(names)
    states = landscape["subset_states"]
    if directions != 14 or len(states) != 1 << directions:
        raise ArithmeticError("landscape state count does not equal 2^14")
    fresh = [row["stage_local_numerators"] for row in states]
    retained = [row["retained_numerators"] for row in states]
    name_index = {name: index for index, name in enumerate(names)}
    recovered = [index for index, name in enumerate(names) if name.startswith("recovered-")]
    tail_names = landscape["historical_comparison"]["attested_M24_to_M31_order"]
    tail = [name_index[name] for name in tail_names]
    if len(recovered) != 7 or len(tail) != 7 or set(recovered) & set(tail):
        raise ArithmeticError("not a 7+7 historical decomposition")
    m24_mask = sum(1 << index for index in recovered)
    tail_actual = order_metrics(tail, m24_mask, retained)
    tail_optimum = optimize_from(m24_mask, tail, retained)
    compatible = []
    for order in itertools.permutations(recovered):
        block = order_metrics(order, 0, retained)
        full_scores = block["scores"] + tail_actual["scores"]
        compatible.append((max(full_scores), sum(full_scores), list(order)))
    compatible_minimax = min(compatible, key=lambda row: (row[0], row[1], row[2]))
    compatible_total = min(compatible, key=lambda row: (row[1], row[0], row[2]))
    all_maxima = sorted(row[0] for row in compatible)
    all_totals = sorted(row[1] for row in compatible)
    audit = monotonicity_audit(fresh, retained, directions)
    tail_rows = landscape["historical_comparison"]["attested_tail_rows"]
    minimax_rows = landscape["optimal_path_rows"]["minimax"]
    tsv_rows = ["path\trank_before\tdirection\tstage_local_numerator\tretained_numerator\tstage_local_scaled_height\tretained_scaled_height"]
    for label, rows in (("historical_tail", tail_rows), ("global_minimax", minimax_rows), ("global_minimum_total", landscape["optimal_path_rows"]["minimum_total"])):
        for row in rows:
            tsv_rows.append("\t".join(map(str, [label, row["rank_before"], row["add_direction"], row["stage_local_numerator"], row["retained_numerator"], f"{qheight(row['stage_local_numerator']):.12g}", f"{qheight(row['retained_numerator']):.12g}"])))
    report = {
        "schema": "elliptic-curves.curve302-exceptional-subgroup-landscape-report.v1",
        "status": "PASS_COMPLETE_FINITE_ORDER_DIAGNOSTIC",
        "inputs": {str(LANDSCAPE.relative_to(ROOT)): sha(LANDSCAPE), str(Path(__file__).relative_to(ROOT)): sha(Path(__file__))},
        "scope": "Retrospective order diagnostic on all 16384 displayed M17+S subgroups. It is not a prospective selector and the CVP values remain finite nearest-plane/local-descent upper bounds.",
        "historical_order_provenance": {
            "M17_to_M24": "Only a seven-dimensional recovered block is attested. No total order on its fixed diagnostic coordinate directions is present in the sealed source, so all 5040 compatible internal orders are summarized instead of inventing one.",
            "M24_to_M31_attested_tail": tail_names,
            "tail_actual": tail_actual,
            "tail_optimum_with_fixed_M24_start": tail_optimum,
            "all_5040_M17_to_M24_compatible_then_historical_tail": {
                "minimum_bottleneck_numerator": compatible_minimax[0], "maximum_bottleneck_numerator": all_maxima[-1],
                "median_bottleneck_numerator": all_maxima[len(all_maxima)//2], "minimum_total_numerator": compatible_total[1],
                "maximum_total_numerator": all_totals[-1], "median_total_numerator": all_totals[len(all_totals)//2],
                "best_bottleneck_compatible_prefix": [names[index] for index in compatible_minimax[2]],
                "best_total_compatible_prefix": [names[index] for index in compatible_total[2]],
            },
        },
        "global_14_direction_optima": landscape["optimization"],
        "raw_vs_cumulative_monotonicity": audit,
        "autonomous_V1_comparison": landscape["historical_comparison"]["autonomous_V1"],
        "artifacts": {"table": str(TSV.relative_to(ROOT)), "figure": str(SVG.relative_to(ROOT))},
        "boundary": "A low rounded half-lattice score does not prove a low pointed-quartic coordinate or that an orbit would have been scheduled by a target-blind policy. The comparison with autonomous V1 is by exactly verified subgroup words, not a claimed direction-by-direction recovery order.",
        "reproducing_command": "python3 elliptic-curves/cas/report_curve302_exceptional_subgroup_landscape.py --check",
    }
    return report, "\n".join(tsv_rows) + "\n", plot_svg(names, tail_rows, minimax_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.build == args.check:
        parser.error("choose exactly one of --build or --check")
    report, table, figure = build()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.build:
        if any(path.exists() for path in (OUTPUT, TSV, SVG)):
            raise FileExistsError("preserve immutable subgroup-landscape report artifacts")
        OUTPUT.write_text(rendered)
        TSV.write_text(table)
        SVG.write_text(figure)
    elif not (OUTPUT.read_text() == rendered and TSV.read_text() == table and SVG.read_text() == figure):
        raise ArithmeticError("subgroup-landscape report replay differs")
    print("CURVE302_SUBGROUP_LANDSCAPE_REPORT|states=16384|compatible_orders=5040|status=PASS", flush=True)


if __name__ == "__main__":
    main()
