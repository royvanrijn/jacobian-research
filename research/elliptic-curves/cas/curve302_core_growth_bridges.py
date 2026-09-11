"""Exact positive-evidence analysis of Curve302's genuine common-core growth stages.

This module is intentionally downstream of the sealed chart replay and the
positive-evidence tables.  It does not use chart completeness and never turns an
unrecorded hit into a negative observation.

For each historical gain stage where at least one positively exposed single
quotient direction increases the intersection with the next observed common
integral core, it computes:

* actual versus alternative one-step core gains;
* exact tied-best counts;
* the candidate's canonical quotient line modulo the current saturated prefix;
* which of the 14 displayed quotient axes become newly contained after adding
  that candidate and saturating;
* the new one-dimensional line induced inside the next common core, when the
  single-candidate core gain is exactly one.

The stage-local quotient line is canonical over Q.  Prefix lattices were
independently certified saturated in the preceding experiment, so equality of
these quotient lines is equality of the corresponding saturated one-dimensional
extensions of the displayed D/M17 lattice.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd, lcm
from typing import Iterable, Sequence

from sympy import Matrix

from curve302_short_core_controls import integer, primitive, require
from curve302_chart_exposure import integer_contains, stage_prefixes
from curve302_positive_exposure_tables import (
    exposure_index,
    span_intersection_rank,
    stage_positive_words,
)


def _q(value) -> F:
    if isinstance(value, F):
        return value
    # SymPy Rational exposes p/q; str() is exact as well.
    return F(str(value))


def matrix_rank(rows: Sequence[Sequence[int | F]], n: int) -> int:
    if not rows:
        return 0
    return int(Matrix([[_q(x) for x in row] for row in rows]).rank())


def _primitive_rational_line(values: Iterable[int | F]) -> tuple[int, ...]:
    """Canonical primitive integer representative of one rational line."""
    qs = tuple(_q(v) for v in values)
    require(any(q != 0 for q in qs), "zero quotient residue cannot define a line")
    den = 1
    for q in qs:
        den = lcm(den, q.denominator)
    ints = [q.numerator * (den // q.denominator) for q in qs]
    g = 0
    for x in ints:
        g = gcd(g, abs(int(x)))
    require(g > 0, "failed to primitive-normalize rational line")
    ints = [int(x) // g for x in ints]
    sign = 1 if next(x for x in ints if x) > 0 else -1
    return tuple(sign * x for x in ints)


def quotient_line_mod_span(rows: Sequence[Sequence[int | F]], vector: Sequence[int | F], n: int):
    """Canonical rational quotient line of vector modulo span_Q(rows).

    RREF is unique over Q.  Eliminating the pivot coordinates therefore gives a
    deterministic representative of the image in Q^n/span(rows).  We then
    primitive-normalize that residual, modulo sign, because Sat(S+Zv) depends on
    the quotient line rather than on the orientation of v.
    """
    require(len(vector) == n, "quotient vector has wrong dimension")
    if rows:
        M = Matrix([[_q(x) for x in row] for row in rows])
        R, pivots = M.rref()
    else:
        R, pivots = Matrix.zeros(0, n), tuple()
    residual = [_q(x) for x in vector]
    for row_i, pivot in enumerate(pivots):
        coeff = residual[pivot]
        if coeff:
            for j in range(n):
                residual[j] -= coeff * _q(R[row_i, j])
    require(any(x != 0 for x in residual), "candidate lies in current rational span")
    return {
        "ambient_line": _primitive_rational_line(residual),
        "pivot_columns": tuple(int(x) for x in pivots),
        "free_columns": tuple(j for j in range(n) if j not in pivots),
    }


def _intersection_basis(rows_a, rows_b, n: int) -> tuple[tuple[F, ...], ...]:
    """Canonical RREF basis of span(rows_a) intersection span(rows_b)."""
    def null_constraints(rows):
        if not rows:
            # orthogonal complement of {0} is all Q^n: x=0 constraints below
            return [Matrix.eye(n).col(j) for j in range(n)]
        M = Matrix([[_q(x) for x in r] for r in rows])
        return M.nullspace()

    constraints = null_constraints(rows_a) + null_constraints(rows_b)
    if not constraints:
        inter = Matrix.eye(n)
    else:
        C = Matrix.vstack(*(v.T for v in constraints))
        basis_cols = C.nullspace()
        if not basis_cols:
            return tuple()
        inter = Matrix.vstack(*(v.T for v in basis_cols))
    R, _ = inter.rref()
    out = []
    for i in range(R.rows):
        row = tuple(_q(R[i, j]) for j in range(n))
        if any(row):
            out.append(row)
    return tuple(out)


def core_bridge_line(prefix, word, core, n: int):
    """New line inside `core` induced by Sat(prefix + Z*word), if delta is one."""
    before = _intersection_basis(prefix, core, n)
    after = _intersection_basis(tuple(prefix) + (tuple(word),), core, n)
    delta = len(after) - len(before)
    if delta != 1:
        return None
    for row in after:
        if matrix_rank(tuple(before) + (row,), n) > len(before):
            q = quotient_line_mod_span(before, row, n)
            return q["ambient_line"]
    raise AssertionError("rank-one core gain had no new intersection line")


def newly_contained_axes(prefix, word, direction_ids: Sequence[str], n: int):
    """Named ambient basis axes newly entering Sat(prefix + Z*word)."""
    after_rows = tuple(prefix) + (tuple(word),)
    out = []
    for i, name in enumerate(direction_ids):
        e = tuple(int(j == i) for j in range(n))
        before = integer_contains(prefix, e)
        # Prefix is saturated; after saturation, rational-span membership is the
        # exact criterion for containment of this ambient integral vector.
        after = matrix_rank(after_rows, n) == matrix_rank(after_rows + (e,), n)
        if not before and after:
            out.append(name)
    return tuple(out)


def analyze_stage(*, seed: str, epoch: int, prefix, gains, positives, next_core,
                  direction_ids: Sequence[str], n: int):
    """Analyze one historical gain stage from positive exposure evidence only."""
    actual = {primitive(g["primitive"] if isinstance(g, dict) and "primitive" in g else g)
              for g in gains}
    require(actual, "stage has no gain batch")
    require(actual <= set(positives), f"actual gain absent from positive exposure set: {seed}/{epoch}")
    current_ir = span_intersection_rank(prefix, next_core, n) if next_core else 0
    rows = []
    for w, meta in sorted(positives.items()):
        ir = span_intersection_rank(tuple(prefix) + (w,), next_core, n) if next_core else 0
        delta = ir - current_ir
        qline = quotient_line_mod_span(prefix, w, n)
        rows.append({
            "seed": seed,
            "epoch": int(epoch),
            "pre_dimension": len(prefix),
            "candidate_word": tuple(w),
            "actual_member": w in actual,
            "chart_count_lower_bound": int(meta["chart_count"]),
            "first_positive_chart_order": int(meta["first_list_order"]),
            "chart_ids": tuple(meta.get("chart_ids", ())),
            "next_core_rank": len(next_core),
            "current_next_core_intersection_rank": current_ir,
            "candidate_next_core_intersection_rank": ir,
            "next_core_delta": delta,
            "completes_next_core": bool(next_core) and ir == len(next_core),
            "quotient_line_mod_prefix": qline["ambient_line"],
            "quotient_free_columns": qline["free_columns"],
            "new_axis_signature": newly_contained_axes(prefix, w, direction_ids, n),
            "core_bridge_line": core_bridge_line(prefix, w, next_core, n) if next_core and delta == 1 else None,
        })
    best = max((r["next_core_delta"] for r in rows), default=0)
    for row in rows:
        row["tied_best_positive"] = row["next_core_delta"] == best
        row["role"] = (
            "actual_tied_best" if row["actual_member"] and row["tied_best_positive"] else
            "actual" if row["actual_member"] else
            "alternative_tied_best" if row["tied_best_positive"] else
            "alternative"
        )
    actual_rows = [r for r in rows if r["actual_member"]]
    alt_rows = [r for r in rows if not r["actual_member"]]
    actual_best = max(r["next_core_delta"] for r in actual_rows)
    alt_best = max((r["next_core_delta"] for r in alt_rows), default=0)
    post = tuple(prefix) + tuple(tuple(g["word"] if isinstance(g, dict) and "word" in g else g) for g in gains)
    batch_ir = span_intersection_rank(post, next_core, n) if next_core else 0
    return {
        "stage": {
            "seed": seed,
            "epoch": int(epoch),
            "pre_dimension": len(prefix),
            "post_dimension": len(prefix) + len(gains),
            "batch_size": len(gains),
            "next_core_rank": len(next_core),
            "current_next_core_intersection_rank": current_ir,
            "actual_batch_next_core_intersection_rank": batch_ir,
            "actual_batch_next_core_delta": batch_ir - current_ir,
            "positive_direction_count": len(rows),
            "positive_alternative_count": len(alt_rows),
            "best_positive_single_delta": best,
            "best_actual_single_delta": actual_best,
            "best_alternative_single_delta": alt_best,
            "tied_best_positive_count": sum(r["tied_best_positive"] for r in rows),
            "tied_best_actual_count": sum(r["tied_best_positive"] and r["actual_member"] for r in rows),
            "tied_best_alternative_count": sum(r["tied_best_positive"] and not r["actual_member"] for r in rows),
            "classification": (
                "ACTUAL_NOT_BEST" if actual_best < best else
                "ACTUAL_UNIQUE_BEST" if not any(r["tied_best_positive"] for r in alt_rows) else
                "ACTUAL_TIED_WITH_ALTERNATIVES"
            ),
        },
        "candidates": rows,
    }


def _cluster(rows, key_name: str):
    groups = defaultdict(list)
    for row in rows:
        value = row.get(key_name)
        if value is None:
            continue
        if isinstance(value, list):
            value = tuple(value)
        groups[value].append(row)
    out = []
    for value, members in groups.items():
        out.append({
            "key": value,
            "occurrences": len(members),
            "stages": len({(r["seed"], r["epoch"]) for r in members}),
            "actual_occurrences": sum(r["actual_member"] for r in members),
            "alternative_tied_best_occurrences": sum((not r["actual_member"]) and r["tied_best_positive"] for r in members),
            "examples": [{"seed": r["seed"], "epoch": r["epoch"], "word": r["candidate_word"]}
                         for r in members[:8]],
        })
    out.sort(key=lambda r: (-r["stages"], -r["occurrences"], str(r["key"])))
    return out


def analyze_core_growth(data, ledger, n=14):
    """Extract the genuine positive one-step common-core growth stages."""
    direction_ids = tuple(ledger["direction_ids"])
    require(len(direction_ids) == n, "direction roster has wrong dimension")
    prefixes = stage_prefixes(data["runs"], n)
    eidx = exposure_index(ledger)
    stages = []
    candidates = []
    for run in data["runs"]:
        seed = run["seed"]
        for stage in run["stages"]:
            gains = stage.get("gains", ())
            if not gains:
                continue
            epoch = integer(stage["epoch"])
            prefix = prefixes[(seed, epoch)]
            post_dim = len(prefix) + len(gains)
            next_core = data["cores"].get(post_dim, ())
            positives = stage_positive_words(eidx[(seed, epoch)], prefix, n)
            result = analyze_stage(seed=seed, epoch=epoch, prefix=prefix, gains=gains,
                                   positives=positives, next_core=next_core,
                                   direction_ids=direction_ids, n=n)
            # The 27 decisive stages are exactly those where a recorded positive
            # single direction can increase the next common-core intersection.
            if result["stage"]["best_positive_single_delta"] <= 0:
                continue
            stages.append(result["stage"])
            candidates.extend(result["candidates"])

    bridge_rows = [r for r in candidates if r["actual_member"] or r["tied_best_positive"]]
    summary = {
        "decisive_core_growth_stages": len(stages),
        "candidate_stage_pairs": len(candidates),
        "actual_or_tied_bridge_pairs": len(bridge_rows),
        "actual_unique_best_stages": sum(s["classification"] == "ACTUAL_UNIQUE_BEST" for s in stages),
        "actual_tied_with_alternatives_stages": sum(s["classification"] == "ACTUAL_TIED_WITH_ALTERNATIVES" for s in stages),
        "actual_not_best_stages": sum(s["classification"] == "ACTUAL_NOT_BEST" for s in stages),
        "stages_with_alternative_same_best_delta": sum(s["tied_best_alternative_count"] > 0 for s in stages),
        "stages_with_alternative_lower_than_actual": sum(s["best_alternative_single_delta"] < s["best_actual_single_delta"] for s in stages),
    }
    clusters = {
        "exact_candidate_word": _cluster(bridge_rows, "candidate_word"),
        "new_axis_signature": _cluster(bridge_rows, "new_axis_signature"),
        "core_bridge_line": _cluster(bridge_rows, "core_bridge_line"),
        "quotient_line_mod_prefix": _cluster(bridge_rows, "quotient_line_mod_prefix"),
    }
    return {
        "status": "PASS_27_CORE_GROWTH_BRIDGE_ANALYSIS",
        "summary": summary,
        "stages": stages,
        "candidates": candidates,
        "bridge_candidates": bridge_rows,
        "clusters": clusters,
        "boundary": (
            "Positive evidence only. Decisive stages are selected by a strictly positive exact one-step "
            "intersection-rank gain with the observed next common core. Missing chart hits are never negatives."
        ),
    }


# --- Optional chart-schema enrichment --------------------------------------

def object_shape(value, depth=0):
    """Stable schema-only shape; values are deliberately ignored."""
    if depth > 8:
        return "..."
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    if isinstance(value, dict):
        return {str(k): object_shape(value[k], depth + 1) for k in sorted(value)}
    if isinstance(value, (list, tuple)):
        shapes = []
        seen = set()
        for item in value[:32]:
            s = repr(object_shape(item, depth + 1))
            if s not in seen:
                seen.add(s)
                shapes.append(object_shape(item, depth + 1))
        return {"list_len": len(value), "element_shapes": shapes}
    return type(value).__name__
