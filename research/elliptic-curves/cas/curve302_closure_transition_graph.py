"""Exact finite transition graph for Curve302 common-core closure.

Consumes the sealed replay/trajectory data and the decisive core-growth bridge
analysis.  No chart-completeness assumption is used and no new point search or
point recognition is performed.

The key mathematical convention is that every historical prefix was already
certified saturated in Z^14.  Hence a saturated state is uniquely determined by
its rational span, so canonical RREF is an exact identifier for the saturated
subgroup span_Q(S) intersect Z^14.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from typing import Sequence

from sympy import Matrix

from curve302_short_core_controls import integer, primitive, require
from curve302_chart_exposure import stage_prefixes, rationally_unknown
from curve302_positive_exposure_tables import exposure_index, stage_positive_words, span_intersection_rank
from curve302_core_growth_bridges import quotient_line_mod_span


def _q(x) -> F:
    return x if isinstance(x, F) else F(str(x))


def matrix_rank(rows, n: int) -> int:
    if not rows:
        return 0
    return int(Matrix([[_q(x) for x in r] for r in rows]).rank())


def canonical_rref(rows: Sequence[Sequence[int | F]], n: int):
    """Unique exact rational subspace representative."""
    if not rows:
        return tuple()
    M = Matrix([[_q(x) for x in r] for r in rows])
    R, _ = M.rref()
    out = []
    for i in range(R.rows):
        row = tuple(_q(R[i, j]) for j in range(n))
        if any(row):
            out.append(row)
    return tuple(out)


def intersection_rref(rows_a, rows_b, n: int):
    """Canonical RREF basis of span_Q(A) intersect span_Q(B)."""
    # U∩V = (U^perp + V^perp)^perp.
    def annihilator(rows):
        if not rows:
            return tuple(tuple(F(int(i == j)) for i in range(n)) for j in range(n))
        M = Matrix([[_q(x) for x in r] for r in rows])
        return tuple(tuple(_q(v[i]) for i in range(n)) for v in M.nullspace())

    cons = annihilator(rows_a) + annihilator(rows_b)
    if not cons:
        return canonical_rref(tuple(tuple(int(i == j) for i in range(n)) for j in range(n)), n)
    C = Matrix([list(r) for r in cons])
    ns = C.nullspace()
    if not ns:
        return tuple()
    return canonical_rref(tuple(tuple(_q(v[i]) for i in range(n)) for v in ns), n)


def _rat_text(x: F) -> str:
    x = _q(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def encode_basis(rows):
    return [[_rat_text(x) for x in r] for r in rows]


def basis_key(rows):
    return tuple(tuple(_rat_text(x) for x in r) for r in rows)


def stable_id(prefix: str, payload) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return f"{prefix}-{sha256(raw).hexdigest()[:16]}"


def _stage_map(data):
    out = {}
    for run in data["runs"]:
        for st in run["stages"]:
            out[(run["seed"], integer(st["epoch"]))] = st
    return out


def _run_map(data):
    return {r["seed"]: r for r in data["runs"]}


def transition_graph(data, ledger, bridge_analysis, n: int = 14):
    """Collapse decisive stage/candidate pairs into exact transition families."""
    require(bridge_analysis.get("status") == "PASS_27_CORE_GROWTH_BRIDGE_ANALYSIS",
            "bridge analysis not passed")
    prefixes = stage_prefixes(data["runs"], n)
    stages_by_key = {(s["seed"], int(s["epoch"])): s for s in bridge_analysis["stages"]}
    final_dim = max(data["cores"])
    final_core = data["cores"][final_dim]

    nodes = {}
    edge_rows = []
    for row in bridge_analysis["bridge_candidates"]:
        seed, epoch = row["seed"], int(row["epoch"])
        stage = stages_by_key[(seed, epoch)]
        prefix = prefixes[(seed, epoch)]
        word = tuple(integer(x) for x in row["candidate_word"])
        before = canonical_rref(prefix, n)
        after = canonical_rref(tuple(prefix) + (word,), n)
        require(len(after) == len(before) + 1, "bridge candidate does not raise saturated rank by one")
        post_dim = int(stage["post_dimension"])
        target = data["cores"].get(post_dim, ())
        core_before = intersection_rref(prefix, target, n) if target else tuple()
        core_after = intersection_rref(tuple(prefix) + (word,), target, n) if target else tuple()
        final_before = intersection_rref(prefix, final_core, n)
        final_after = intersection_rref(tuple(prefix) + (word,), final_core, n)

        b_payload = encode_basis(before)
        a_payload = encode_basis(after)
        b_id = stable_id("S", b_payload)
        a_id = stable_id("S", a_payload)
        nodes.setdefault(b_id, {"state_id": b_id, "rank": len(before), "basis_rref": b_payload})
        nodes.setdefault(a_id, {"state_id": a_id, "rank": len(after), "basis_rref": a_payload})

        core_edge_key = {
            "target_post_dimension": post_dim,
            "before": encode_basis(core_before),
            "after": encode_basis(core_after),
        }
        edge_rows.append({
            "seed": seed,
            "epoch": epoch,
            "pre_dimension": int(stage["pre_dimension"]),
            "post_dimension": post_dim,
            "classification": stage["classification"],
            "actual_member": bool(row["actual_member"]),
            "tied_best_positive": bool(row["tied_best_positive"]),
            "candidate_word": list(word),
            "before_state_id": b_id,
            "after_state_id": a_id,
            "literal_edge_id": stable_id("E", {"before": b_payload, "after": a_payload}),
            "core_edge_id": stable_id("C", core_edge_key),
            "core_before_rank": len(core_before),
            "core_after_rank": len(core_after),
            "core_delta": len(core_after) - len(core_before),
            "final_common_before_rank": len(final_before),
            "final_common_after_rank": len(final_after),
            "quotient_line_mod_prefix": list(row["quotient_line_mod_prefix"]),
            "core_bridge_line": list(row["core_bridge_line"]) if row.get("core_bridge_line") is not None else None,
            "new_axis_signature": list(row.get("new_axis_signature", ())),
            "chart_count_lower_bound": int(row.get("chart_count_lower_bound", 0)),
        })

    def aggregate(key_name):
        groups = defaultdict(list)
        for r in edge_rows:
            groups[r[key_name]].append(r)
        out = []
        for key, members in groups.items():
            qlines = {tuple(r["quotient_line_mod_prefix"]) for r in members}
            words = {tuple(r["candidate_word"]) for r in members}
            stages = {(r["seed"], r["epoch"]) for r in members}
            classes = sorted({r["classification"] for r in members})
            row = {
                key_name: key,
                "occurrences": len(members),
                "stages": len(stages),
                "distinct_candidate_words": len(words),
                "distinct_bridge_lines_mod_prefix": len(qlines),
                "actual_occurrences": sum(r["actual_member"] for r in members),
                "alternative_tied_occurrences": sum((not r["actual_member"]) and r["tied_best_positive"] for r in members),
                "unique_best_stage_occurrences": len({(r["seed"], r["epoch"]) for r in members if r["classification"] == "ACTUAL_UNIQUE_BEST"}),
                "classifications": classes,
                "examples": [{"seed": r["seed"], "epoch": r["epoch"], "word": r["candidate_word"]} for r in members[:8]],
            }
            if key_name == "core_edge_id":
                first = members[0]
                row.update({
                    "target_post_dimension": first["post_dimension"],
                    "core_before_rank": first["core_before_rank"],
                    "core_after_rank": first["core_after_rank"],
                    "core_delta": first["core_delta"],
                    "observed_bridge_multiplicity": len(qlines),
                    "bottleneck_class": "NARROW_OBSERVED_GATE" if len(qlines) == 1 else "BROAD_OBSERVED_BASIN",
                })
            out.append(row)
        out.sort(key=lambda r: (-r["stages"], -r["occurrences"], r[key_name]))
        return out

    literal = aggregate("literal_edge_id")
    core = aggregate("core_edge_id")
    unique_stages = [s for s in bridge_analysis["stages"] if s["classification"] == "ACTUAL_UNIQUE_BEST"]
    unique_rows = []
    for s in unique_stages:
        members = [r for r in edge_rows if r["seed"] == s["seed"] and r["epoch"] == int(s["epoch"]) and r["actual_member"]]
        unique_rows.append({
            "seed": s["seed"], "epoch": int(s["epoch"]),
            "pre_dimension": int(s["pre_dimension"]), "post_dimension": int(s["post_dimension"]),
            "core_before_rank": int(s["current_next_core_intersection_rank"]),
            "core_delta": int(s["actual_batch_next_core_delta"]),
            "actual_core_edge_ids": sorted({r["core_edge_id"] for r in members}),
            "actual_bridge_lines": [r["quotient_line_mod_prefix"] for r in members],
        })

    return {
        "status": "PASS_EXACT_CLOSURE_TRANSITION_GRAPH",
        "summary": {
            "decisive_stages": len(bridge_analysis["stages"]),
            "bridge_candidate_pairs": len(edge_rows),
            "literal_saturated_edges": len(literal),
            "common_core_transition_types": len(core),
            "narrow_observed_core_gates": sum(r.get("bottleneck_class") == "NARROW_OBSERVED_GATE" for r in core),
            "broad_observed_core_basins": sum(r.get("bottleneck_class") == "BROAD_OBSERVED_BASIN" for r in core),
            "unique_best_stages": len(unique_rows),
        },
        "nodes": sorted(nodes.values(), key=lambda r: (r["rank"], r["state_id"])),
        "bridge_edges": edge_rows,
        "literal_edges": literal,
        "core_transition_types": core,
        "unique_best_stages": unique_rows,
        "boundary": "Positive bridge evidence only. Saturated states are identified exactly by rational RREF because all historical prefix lattices were independently certified saturated.",
    }


def rank29_control(data, ledger, bridge_analysis, graph, n: int = 14, *, stall_dimension: int = 12, next_dimension: int = 13, final_dimension: int = 14):
    """Compare the unique rank-29 terminal state with successful late bridge gates."""
    prefixes = stage_prefixes(data["runs"], n)
    eidx = exposure_index(ledger)
    stalled = [r for r in data["runs"] if int(r["final_dimension"]) == stall_dimension]
    require(len(stalled) == 1, f"expected one quotient-dimension-{stall_dimension} stalled run")
    run = stalled[0]
    seed = run["seed"]

    terminal = None
    for st in run["stages"]:
        if not st.get("gains"):
            terminal = st
    if terminal is not None:
        epoch = integer(terminal["epoch"])
        prefix = prefixes[(seed, epoch)]
        terminal_charts = eidx.get((seed, epoch), tuple())
    else:
        # Compatibility with old trajectory snapshots lacking the explicit
        # terminal no-gain stage: reconstruct the final prefix exactly.
        rows = [tuple(int(i == run["seed_index"]) for i in range(n))]
        for st in run["stages"]:
            for gain in st.get("gains", ()):
                rows.append(tuple(integer(x) for x in gain["word"]))
        prefix = tuple(rows)
        epoch = None
        terminal_charts = tuple()
    require(matrix_rank(prefix, n) == stall_dimension, f"stalled terminal prefix is not rank {stall_dimension}")

    terminal_pos = stage_positive_words(terminal_charts, prefix, n) if terminal_charts else {}
    c13 = data["cores"].get(next_dimension, ())
    c14 = data["cores"].get(final_dimension, ())
    require(c13 and c14, f"successful dimension-{next_dimension}/{final_dimension} common cores unavailable")
    before13 = span_intersection_rank(prefix, c13, n)
    before14 = span_intersection_rank(prefix, c14, n)

    terminal_rows = []
    for w, meta in sorted(terminal_pos.items()):
        d13 = span_intersection_rank(tuple(prefix) + (w,), c13, n) - before13
        d14 = span_intersection_rank(tuple(prefix) + (w,), c14, n) - before14
        q = quotient_line_mod_span(prefix, w, n)["ambient_line"]
        terminal_rows.append({
            "word": list(w), "quotient_line_mod_stall": list(q),
            "c13_delta": d13, "c14_delta": d14,
            "chart_count_lower_bound": int(meta["chart_count"]),
            "chart_ids": list(meta.get("chart_ids", ())),
        })

    # Successful bridge representatives occurring strictly after quotient dim12.
    late = [r for r in graph["bridge_edges"] if int(r["pre_dimension"]) >= stall_dimension]
    projected = []
    for r in late:
        w = tuple(integer(x) for x in r["candidate_word"])
        if not rationally_unknown(prefix, w, n):
            continue
        q = quotient_line_mod_span(prefix, w, n)["ambient_line"]
        d13 = span_intersection_rank(tuple(prefix) + (w,), c13, n) - before13
        d14 = span_intersection_rank(tuple(prefix) + (w,), c14, n) - before14
        projected.append({
            "source_seed": r["seed"], "source_epoch": r["epoch"],
            "source_pre_dimension": r["pre_dimension"], "actual_member": r["actual_member"],
            "quotient_line_mod_stall": list(q), "candidate_word": r["candidate_word"],
            "c13_delta_if_added_to_stall": d13, "c14_delta_if_added_to_stall": d14,
        })

    def line_groups(rows):
        groups = defaultdict(list)
        for r in rows:
            groups[tuple(r["quotient_line_mod_stall"])].append(r)
        out = []
        terminal_lines = {tuple(r["quotient_line_mod_stall"]) for r in terminal_rows}
        for line, members in groups.items():
            out.append({
                "quotient_line_mod_stall": list(line),
                "successful_occurrences": len(members),
                "successful_stages": len({(r["source_seed"], r["source_epoch"]) for r in members}),
                "actual_occurrences": sum(r["actual_member"] for r in members),
                "max_c13_delta_if_added_to_stall": max(r["c13_delta_if_added_to_stall"] for r in members),
                "max_c14_delta_if_added_to_stall": max(r["c14_delta_if_added_to_stall"] for r in members),
                "terminal_positive_match": line in terminal_lines,
            })
        out.sort(key=lambda r: (-r["successful_stages"], -r["successful_occurrences"], r["quotient_line_mod_stall"]))
        return out

    families = line_groups(projected)
    terminal_c13 = [r for r in terminal_rows if r["c13_delta"] > 0]
    terminal_c14 = [r for r in terminal_rows if r["c14_delta"] > 0]
    c13_deficit = len(c13) - before13
    c14_deficit = len(c14) - before14
    if c13_deficit > 0:
        first_gate_status = "POSITIVE_TERMINAL_C13_BRIDGE_EXISTS" if terminal_c13 else "NO_POSITIVE_TERMINAL_C13_BRIDGE"
    elif c14_deficit > 0:
        first_gate_status = "POSITIVE_TERMINAL_C14_BRIDGE_EXISTS" if terminal_c14 else "NO_POSITIVE_TERMINAL_C14_BRIDGE"
    else:
        first_gate_status = "NO_LATE_COMMON_CORE_DEFICIT"

    return {
        "status": "PASS_RANK29_LATE_GATE_CONTROL",
        "stalled_seed": seed,
        "terminal_epoch": epoch,
        "terminal_prefix_rank": stall_dimension,
        "dimension13_common_core_rank": len(c13),
        "dimension13_intersection_rank": before13,
        "dimension13_deficit": c13_deficit,
        "dimension14_common_core_rank": len(c14),
        "dimension14_intersection_rank": before14,
        "dimension14_deficit": c14_deficit,
        "terminal_positive_new_directions": len(terminal_rows),
        "terminal_positive_c13_improvers": len(terminal_c13),
        "terminal_positive_c14_improvers": len(terminal_c14),
        "successful_late_bridge_pairs": len(late),
        "successful_late_bridge_pairs_new_mod_stall": len(projected),
        "successful_late_bridge_line_families_mod_stall": len(families),
        "successful_line_families_seen_positive_at_terminal": sum(r["terminal_positive_match"] for r in families),
        "first_missing_gate_status": first_gate_status,
        "terminal_positive_candidates": terminal_rows,
        "successful_late_bridge_families_mod_stall": families,
        "successful_late_bridge_projections": projected,
        "boundary": "Positive terminal evidence only. Absence of a recorded terminal bridge is not chart-completeness evidence; it is a lower-bound comparison against positively recorded successful bridge families.",
    }
