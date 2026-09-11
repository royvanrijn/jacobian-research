#!/usr/bin/env python3
"""Three read-only, retrospective Curve302 follow-ups. No point-search backend.

Python 3.10+, NumPy and psutil (also usable with `sage -python`).
Use `run` for a new folder, `resume` only between sealed stages, and `check`
for deterministic recomputation. The previous experiment's committed REPORT
and its three output files are required. See the companion runbook.
"""
from __future__ import annotations

import argparse
import ast
from bisect import bisect_left, bisect_right
from collections import Counter
from contextlib import contextmanager
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
import math
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import tempfile
import time

# Set before importing NumPy, including when launched using Sage's Python.
for _env in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_env] = "1"

NAMES = tuple([f"recovered-local-{i:02d}" for i in range(1, 5)] +
              [f"recovered-strict-{i:02d}" for i in range(1, 4)] +
              [f"residual-strict-{i:02d}" for i in range(1, 8)])
STAGES = ("persistence", "strict-filtration", "next-moves")
SCHEMA = "curve302-closure-followup.v1"
BOUNDARY = (
    "Retrospective analysis of a supplied rounded-height form, finite retained atlas, "
    "and already reconstructed trajectories on one known curve. No new point search, "
    "rank claim, exact canonical-height assertion, prospective predictor or propagation theorem. "
    "A successful check recomputes these analyses; it does not repeat the source's EC coordinate proofs."
)
MAX_JSON_BYTES = 256 * 1024**2


def require(test, message):
    if not test:
        raise ValueError(message)


def integer(value):
    require(type(value) is int, f"expected a JSON integer, not {value!r}")
    return value


def read(path):
    path = Path(path)
    require(path.stat().st_size <= MAX_JSON_BYTES, f"oversize input: {path}")
    def pairs(items):
        obj = {}
        for key, value in items:
            require(key not in obj, f"duplicate JSON key {key}")
            obj[key] = value
        return obj
    def bad(value):
        raise ValueError(f"non-finite JSON number {value}")
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=pairs, parse_constant=bad)


def digest(path):
    h = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024**2), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def atomic(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix="." + path.name + ".", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(obj, stream, sort_keys=True, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def parity(row):
    return sum((integer(x) & 1) << i for i, x in enumerate(row))


def primitive(row):
    row = tuple(integer(v) for v in row)
    g = 0
    for v in row:
        g = math.gcd(g, abs(v))
    if not g:
        return row
    sign = 1 if next(v for v in row if v) > 0 else -1
    return tuple(sign * v // g for v in row)


def rref(rows, n):
    pivots = {}
    for raw in rows:
        v = integer(raw)
        require(0 <= v < 1 << n, "binary row exceeds ambient dimension")
        for p in sorted(pivots):
            if v >> p & 1:
                v ^= pivots[p]
        if not v:
            continue
        p = (v & -v).bit_length() - 1
        for q in list(pivots):
            if pivots[q] >> p & 1:
                pivots[q] ^= v
        pivots[p] = v
    return tuple(pivots[p] for p in sorted(pivots))


def in_span(v, rows, n):
    for b in rref(rows, n):
        p = (b & -b).bit_length() - 1
        if v >> p & 1:
            v ^= b
    return v == 0


def orthogonal(rows, n):
    reduced = rref(rows, n)
    pivots = {(b & -b).bit_length()-1: b for b in reduced}
    answer = []
    for f in range(n):
        if f in pivots:
            continue
        v = 1 << f
        for p, row in pivots.items():
            if row >> f & 1:
                v |= 1 << p
        answer.append(v)
    return rref(answer, n)


def intersection(a, b, n):
    return orthogonal((*orthogonal(a, n), *orthogonal(b, n)), n)


def span_word(target, generators, n):
    """Return an exact F2 combination of original generators, not RREF labels."""
    pivots = {}
    for i, v in enumerate(generators):
        word = 1 << i
        for p in sorted(pivots):
            if v >> p & 1:
                v ^= pivots[p][0]
                word ^= pivots[p][1]
        if v:
            p = (v & -v).bit_length()-1
            pivots[p] = (v, word)
    word, v = 0, target
    for p in sorted(pivots):
        if v >> p & 1:
            v ^= pivots[p][0]
            word ^= pivots[p][1]
    require(v == 0, "intersection vector has no generating word")
    check = 0
    for i in bits(word):
        check ^= generators[i]
    require(check == target, "intersection witness failed")
    return list(bits(word))


def rational_spd(matrix):
    """Exact LDL positivity check; no tolerance and no Sage dependency."""
    q = [[F(v) for v in row] for row in matrix]
    n = len(q)
    require(n > 0 and all(len(row) == n for row in q), "non-square quotient form")
    require(all(q[i][j] == q[j][i] for i in range(n) for j in range(n)), "asymmetric quotient form")
    lower = [[F(i == j) for j in range(n)] for i in range(n)]
    diagonal = []
    for i in range(n):
        d = q[i][i] - sum(lower[i][k]**2 * diagonal[k] for k in range(i))
        require(d > 0, "quotient form is not exactly positive definite")
        diagonal.append(d)
        for j in range(i+1, n):
            lower[j][i] = (q[j][i] - sum(lower[j][k]*lower[i][k]*diagonal[k] for k in range(i))) / d
    return q


def validate_landscape(source, expected_names=NAMES):
    require(source.get("schema") == "elliptic-curves.curve302-exceptional-subgroup-landscape.v1", "wrong landscape schema")
    require(source.get("status") == "PASS_RETROSPECTIVE_FINITE_SUBGROUP_LANDSCAPE", "landscape is not passed")
    names = tuple(source["fixed_basis"]["direction_ids"])
    require(names == tuple(expected_names), "landscape directions changed")
    n = len(names)
    require(integer(source["fixed_basis"]["generic_rank"]) == 17, "wrong generic rank")
    require(integer(source["fixed_basis"]["rounded_metric_scale"]) > 0, "invalid metric scale")
    rows = source["subset_states"]
    require(len(rows) == 1 << n, "incomplete subset cube")
    costs = []
    for mask, row in enumerate(rows):
        require(integer(row["state_mask"]) == mask, "reordered subset cube")
        values = row["retained_numerators"]
        require(len(values) == n, "wrong edge width")
        for j, value in enumerate(values):
            if mask >> j & 1:
                require(value is None, "active target should not have an edge")
            else:
                require(integer(value) >= 0, "negative cost")
        costs.append(values)
    for mask, values in enumerate(costs):
        for i in bits(mask):
            parent = costs[mask ^ (1 << i)]
            require(all(values[j] <= parent[j] for j in range(n) if not mask >> j & 1), "retained atlas is not monotone")
    return names, costs


def parsed_rref(value, n):
    if isinstance(value, str):
        value = ast.literal_eval(value)  # numeric tuples only; never eval
    require(isinstance(value, (list, tuple)), "invalid RREF encoding")
    original = tuple(integer(v) for v in value)
    require(rref(original, n) == original, "noncanonical source RREF")
    return original


def validate_results(inputs, expected_names=NAMES):
    """Bind the new analysis to the previous sealed output bundle."""
    inputs = Path(inputs)
    source = read(inputs / "landscape.json")
    laws = read(inputs / "closure-laws.json")
    rel = read(inputs / "quotient-relations.json")
    traj = read(inputs / "trajectories.json")
    report = read(inputs / "REPORT.json")
    require(report.get("status") == "PASS_THREE_CLOSURE_EXPERIMENTS", "source experiment not sealed/passed")
    for stage in ("closure-laws", "quotient-relations", "trajectories"):
        require(report["outputs"][stage] == digest(inputs / f"{stage}.json"), f"source report hash mismatch: {stage}")
    require(laws.get("schema") == "curve302-closure-laws.v1" and laws.get("status") == "PASS_EXHAUSTIVE_16384_CLOSURE_ANALYSIS", "source closure report is invalid")
    require(laws["input_sha256"] == digest(inputs / "landscape.json"), "source landscape binding mismatch")
    names, costs = validate_landscape(source, expected_names)
    n = len(names)
    require(rel.get("schema") == "curve302-quotient-near-relations.v1" and rel.get("status") == "PASS_QUOTIENT_RELATION_ANALYSIS", "source quotient report is invalid")
    require(traj.get("schema") == "curve302-actual-v3-closure-trajectories.v1" and traj.get("status") == "PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED", "source trajectory report is invalid")
    for obj in (laws, rel, traj):
        require(tuple(obj["direction_ids"]) == names, "cross-artifact direction mismatch")
    require(rel["rounded_metric_scale"] == source["fixed_basis"]["rounded_metric_scale"], "metric units differ")
    q = rational_spd(rel["schur_quotient"])
    require(len(q) == n, "quotient form dimension differs")
    strict = parsed_rref(rel["strict_quotient_rref_bitmasks"], n)
    require(len(strict) == n - 4, "strict quotient is not codimension four")
    require(len(traj["runs"]) == n and {r["seed"] for r in traj["runs"]} == set(names), "incomplete or duplicated trajectory roster")
    runs = []
    for raw in sorted(traj["runs"], key=lambda row: names.index(row["seed"])):
        r = {"seed": raw["seed"], "final_rank": integer(raw["final_rank"]), "charts": integer(raw["charts"]), "events": []}
        seed = names.index(r["seed"])
        rows = [1 << seed]
        rank = 18
        for epoch, stage in enumerate(raw["stages"]):
            require(integer(stage["epoch"]) == epoch, "missing/reordered epoch")
            require(integer(stage["before"]) == rank, "noncontiguous trajectory")
            after = integer(stage["after"])
            require(after == rank + len(stage["new"]), "gain count differs from acquired rows")
            for obj in stage["new"]:
                require(integer(obj["denominator"]) == 1 and obj.get("integral") is True, "nonintegral trajectory coordinates")
                v = tuple(integer(x) for x in obj["quotient_word"])
                require(len(v) == n and any(v), "zero or wrong-width transition")
                p = primitive(v)
                require(tuple(obj["primitive_quotient_word"]) == p, "primitive transition differs")
                mask = parity(v)
                require(integer(obj["mod2_bitmask"]) == mask, "parity differs from original word")
                require(not in_span(mask, rows, n), "acquisition is mod-2 dependent")
                require(obj["mod2_strict"] == in_span(mask, strict, n), "strict label differs")
                rows.append(mask)
                r["events"].append({"word": v, "primitive": p, "epoch": epoch})
            rank = after
            require(integer(stage["quotient_mod2_rank"]) == rank - 17 == len(rref(rows, n)), "stage rank differs")
            require(parsed_rref(stage["quotient_subspace_rref"], n) == rref(rows, n), "stage subspace differs")
        require(rank == r["final_rank"] and 18 <= rank <= 17+n, "terminal rank differs")
        r["seed_index"] = seed
        runs.append(r)
    return {"names": names, "costs": costs, "form": q, "strict": strict, "runs": runs,
            "scale": rel["rounded_metric_scale"], "source": source,
            "dataset_kind": "SYNTHETIC_TEST_FIXTURE" if source.get("synthetic_fixture") else "SUPPLIED_SEALED_RESULTS"}


# Experiment 1 -------------------------------------------------------------
def closure(costs, start, threshold):
    n = len(costs[0])
    mask = start
    while True:
        add = sum(1 << j for j in range(n) if not mask >> j & 1 and costs[mask][j] <= threshold)
        if not add:
            return mask
        mask |= add


def minimax_arrivals(costs, start):
    """Exact DAG minimax. Monotone costs imply reachable unions are reachable."""
    n, size = len(costs[0]), len(costs)
    best = [None] * size
    previous = [None] * size
    best[start] = 0
    arrivals, masks = [None] * n, [None] * n
    for j in bits(start):
        arrivals[j], masks[j] = 0, start
    for mask in range(size):
        if best[mask] is None:
            continue
        for j in range(n):
            if mask >> j & 1:
                continue
            child = mask | (1 << j)
            candidate = max(best[mask], costs[mask][j])
            if best[child] is None or candidate < best[child]:
                best[child], previous[child] = candidate, (mask, j)
            if arrivals[j] is None or (candidate, child) < (arrivals[j], masks[j]):
                arrivals[j], masks[j] = candidate, child
    require(all(v is not None for v in arrivals), "unreachable target")
    require(best[-1] == max(arrivals), "reachable-union minimax mismatch")
    witnesses = []
    for j, target_mask in enumerate(masks):
        path, mask = [], target_mask
        while mask != start:
            parent, target = previous[mask]
            path.append(target)
            mask = parent
        path.reverse()
        state, peak = start, 0
        for target in path:
            peak = max(peak, costs[state][target])
            state |= 1 << target
        require(state >> j & 1 and peak <= arrivals[j], "minimax path witness failed")
        witnesses.append(path)
    return arrivals, best[-1], witnesses


def growth_steps(arrivals):
    return [{"threshold": t, "mask": sum(1 << i for i, a in enumerate(arrivals) if a <= t),
             "size": sum(a <= t for a in arrivals)} for t in sorted({0, *arrivals})]


def persistence(data, policy):
    names, costs = data["names"], data["costs"]
    curves = []
    cutoffs = sorted({value for row in costs for value in row if value is not None})
    n = len(names)
    for start in [0, *(1 << i for i in range(n))]:
        arrival, full, witnesses = minimax_arrivals(costs, start)
        # This checks BOTH sides of every breakpoint using a separate algorithm.
        probes = sorted({0, *arrival, *(t-1 for t in arrival if t > 0)})
        for t in probes:
            predicted = sum(1 << j for j, a in enumerate(arrival) if a <= t)
            require(closure(costs, start, t) == predicted, "parametric curve disagrees with direct fixed-point closure")
        steps = growth_steps(arrival)
        curves.append({"seed": None if not start else names[start.bit_length()-1], "start_mask": start,
                       "arrivals": arrival, "full_threshold": full, "growth": steps,
                       "arrival_path_witnesses": witnesses, "direct_boundary_checks": len(probes)})
    t0 = curves[0]["full_threshold"]
    empty = curves[0]["arrivals"]
    for row in curves[1:]:
        require(row["full_threshold"] <= t0, "seed worsens monotone closure")
        row["full_threshold_advantage"] = t0-row["full_threshold"]
        row["advantage_fraction"] = None if t0 == 0 else (t0-row["full_threshold"])/t0
        # Exact continuous-threshold area, with all breakpoints integers.
        row["closure_size_advantage_area_to_empty_full"] = sum(max(0, t0-a) for a in row["arrivals"]) - sum(max(0, t0-a) for a in empty)
        t = max(0, t0-1)
        row["size_just_below_empty_full"] = sum(a <= t for a in row["arrivals"]) if t0 else None
    return {"schema": "curve302-followup-persistence.v1", "status": "COMPLETE",
            "distinct_edge_cutoff_count": len(cutoffs), "edge_cutoffs_sha256": sha256(canonical(cutoffs)).hexdigest(),
            "curve_encoding": "Right-continuous step functions: these breakpoints encode closures at EVERY retained edge cutoff, without repeated full scans.",
            "curves": curves, "empty_full_threshold": t0,
            "strictly_helpful_singletons": sum(row["full_threshold_advantage"] > 0 for row in curves[1:]),
            "boundary": "Weighted monotone hypergraph accessibility on named axes, not an ordinary graph component, intrinsic lattice invariant or runtime bound."}


# Experiment 2 -------------------------------------------------------------
def filtration(data, policy):
    names, strict = data["names"], data["strict"]
    n = len(names)
    by_dim, runs = {}, []
    for run in data["runs"]:
        generators = [1 << run["seed_index"]]
        states = []
        for step in range(len(run["events"])+1):
            if step:
                generators.append(parity(run["events"][step-1]["word"]))
            subspace = rref(generators, n)
            meet = intersection(subspace, strict, n)
            r, a = len(subspace), len(meet)
            b = r-a
            require(0 <= a <= len(strict) and 0 <= b <= n-len(strict), "invalid strict/local dimensions")
            witnesses = [{"strict_vector": v, "generator_indices": span_word(v, generators, n)} for v in meet]
            record = {"step": step, "quotient_dimension": r, "strict_dimension": a, "local_image_dimension": b,
                      "subspace": subspace, "strict_intersection": meet, "strict_combination_witnesses": witnesses}
            if states:
                delta = (a-states[-1]["strict_dimension"], b-states[-1]["local_image_dimension"])
                require(delta in ((1, 0), (0, 1)), "one acquisition did not add one filtered dimension")
                record["strict_increment"], record["local_increment"] = delta
                record["acquired_vector_itself_strict"] = in_span(generators[-1], strict, n)
            states.append(record)
            by_dim.setdefault(r, []).append((run["seed"], subspace))
        terminal = states[-1]
        runs.append({"seed": run["seed"], "states": states,
                     "first_full_local_quotient_dimension": next((s["quotient_dimension"] for s in states if s["local_image_dimension"] == n-len(strict)), None),
                     "endpoint": {k: terminal[k] for k in ("quotient_dimension", "strict_dimension", "local_image_dimension")},
                     "remaining_strict_dimensions": len(strict)-terminal["strict_dimension"],
                     "remaining_local_dimensions": n-len(strict)-terminal["local_image_dimension"]})
    comparisons = []
    for r, entries in sorted(by_dim.items()):
        common = entries[0][1]
        for _, sub in entries[1:]:
            common = intersection(common, sub, n)
        pairs = []
        lower = max(0, 2*r-n)
        for (name_a, sub_a), (name_b, sub_b) in combinations(entries, 2):
            meet = intersection(sub_a, sub_b, n)
            dim = len(meet)
            require(dim >= lower, "intersection violates ambient dimension bound")
            pairs.append({"a": name_a, "b": name_b, "intersection_dimension": dim,
                          "grassmann_distance": r-dim, "intersection_excess_over_dimension_bound": dim-lower})
        comparisons.append({"quotient_dimension": r, "runs_present": len(entries),
                            "distinct_subspaces": len({sub for _, sub in entries}),
                            "all_run_common_dimension": len(common), "all_run_common_subspace": common,
                            "unavoidable_pairwise_intersection_lower_bound": lower, "pairs": pairs,
                            "terminal_convergence_is_tautological": r == n})
    return {"schema": "curve302-followup-strict-filtration.v1", "status": "COMPLETE", "runs": runs,
            "subspace_comparisons": comparisons,
            "boundary": "Exact F2 combinations in the specified quotient. Strict membership is modulo the generic image, not a new local/class-group computation. Large late intersections have a dimension-forced component."}


# Experiment 3 -------------------------------------------------------------
def fixed_vocabulary(n, profile="full"):
    """Frozen independently of transition outcomes; no target insertion."""
    candidates = {tuple(int(i == j) for i in range(n)) for j in range(n)}
    sizes = (2,) if profile == "smoke" else (2, 3, 4)
    for size in sizes:
        if size > n:
            continue
        choices = (-1, 1) if size == 4 else (-2, -1, 1, 2)
        for indices in combinations(range(n), size):
            for coefficients in product(choices, repeat=size):
                values = [0]*n
                for i, c in zip(indices, coefficients):
                    values[i] = c
                candidates.add(primitive(values))
    return sorted(candidates)


def rank_interval(scores, target, ascending, atol=1e-10, rtol=1e-9):
    """Conservative tie interval; no tie-break by the observed target's index."""
    import numpy as np
    s = np.asarray(scores, dtype=float)
    value = float(s[target])
    tol = atol + rtol*abs(value)
    delta = (s-value) if ascending else (value-s)
    better = int(np.count_nonzero(delta < -tol))
    equal = int(np.count_nonzero(np.abs(delta) <= tol))
    require(equal > 0, "target was lost to numeric tie logic")
    return better+1, better+equal


def all_percentiles(scores, ascending):
    import numpy as np
    values = np.asarray(scores, dtype=float)
    if not ascending:
        values = -values
    if len(values) == 1:
        return np.array([0.5])
    ordered = np.sort(values)
    tol = 1e-10 + 1e-9*np.abs(values)
    lower = np.searchsorted(ordered, values-tol, side="left")
    upper = np.searchsorted(ordered, values+tol, side="right")
    return (lower+(upper-lower-1)/2)/(len(values)-1)


def projection_scores(q, prior, candidates):
    """Whitened QR is more stable than subtracting two large Schur matrices."""
    import numpy as np
    q = np.asarray(q, dtype=float)
    require(np.all(np.isfinite(q)), "non-finite quotient form")
    lower = np.linalg.cholesky(q)
    prior = np.asarray(prior, dtype=float)
    require(np.all(np.isfinite(prior)) and np.max(np.abs(prior)) < 2**40, "active words too large for the frozen numerical policy; UNKNOWN")
    basis, upper = np.linalg.qr((prior @ lower).T, mode="reduced")
    require(np.min(np.abs(np.diag(upper))) > 1e-12*np.linalg.norm(prior @ lower), "ill-conditioned active subgroup; numerical result UNKNOWN")
    c = np.asarray(candidates, dtype=float)
    require(np.all(np.isfinite(c)) and np.max(np.abs(c)) < 2**40, "candidate words too large for the frozen numerical policy; UNKNOWN")
    z = c @ lower
    residual = z - (z @ basis) @ basis.T
    residual_norm = np.einsum("ij,ij->i", residual, residual)
    full_norm = np.einsum("ij,ij->i", z, z)
    # Here candidates are already outside the mod-2 span; hence outside the
    # rational span as the admitted prefix has full mod-2 rank.
    require(np.all(residual_norm > 1e-12 * np.maximum(full_norm, 1)), "numerically zero residual for eligible candidate; UNKNOWN")
    target = lower - (lower @ basis) @ basis.T
    target_norm = np.einsum("ij,ij->i", target, target)
    present = target_norm > 1e-12*np.maximum(np.diag(q), 1)
    require(np.any(present), "no remaining target directions")
    unit_targets = target[present] / np.sqrt(target_norm[present])[:, None]
    correlations = (residual / np.sqrt(residual_norm)[:, None]) @ unit_targets.T
    fractional = np.sum(correlations**2, axis=1)
    require(np.all(np.isfinite(fractional)), "non-finite unlock score")
    return {"residual_norm": residual_norm, "axis_fractional_unlock": fractional,
            "static_full_norm": full_norm}


def next_moves(data, policy):
    import numpy as np
    names, runs = data["names"], data["runs"]
    n = len(names)
    q = np.array([[float(v) for v in row] for row in data["form"]])
    require(np.all(np.isfinite(q)), "quotient form overflow")
    normalization = float(np.median(np.diag(q)))
    require(normalization > 0, "invalid normalization")
    q /= normalization
    eigenvalues = np.linalg.eigvalsh(q)
    require(eigenvalues[0] > 0 and eigenvalues[-1]/eigenvalues[0] < 1e10, "ill-conditioned quotient form; no numeric rankings claimed")
    fixed = fixed_vocabulary(n, policy["profile"])
    axes = {tuple(int(i == j) for i in range(n)) for j in range(n)}
    rng = np.random.default_rng(policy["random_seed"])
    scores_direction = {"residual_norm": True, "axis_fractional_unlock": False, "static_full_norm": True}
    summary, arms = {}, {}
    for arm in ("fixed_bounded", "leave_one_run_out"):
        records = []
        per_run = []
        for run in runs:
            vocabulary = fixed if arm == "fixed_bounded" else sorted(axes | {
                event["primitive"] for other in runs if other["seed"] != run["seed"] for event in other["events"]})
            lookup = {v: i for i, v in enumerate(vocabulary)}
            masks = [parity(v) for v in vocabulary]
            prior = [[int(i == run["seed_index"]) for i in range(n)]]
            binary = [1 << run["seed_index"]]
            run_records = []
            reference_draws = {s: np.zeros(policy["random_draws"]) for s in scores_direction}
            covered_count = 0
            for step, event in enumerate(run["events"]):
                span = {0}
                for b in rref(binary, n):
                    span |= {v ^ b for v in list(span)}
                indices = [i for i, v in enumerate(masks) if v not in span]
                candidates = [vocabulary[i] for i in indices]
                index_map = {index: j for j, index in enumerate(indices)}
                wanted = lookup.get(event["primitive"])
                actual = index_map.get(wanted)
                record = {"seed": run["seed"], "step": step, "epoch": event["epoch"],
                          "quotient_rank_before": len(prior), "remaining_real_dimension": n-len(prior),
                          "original_quotient_word": list(event["word"]), "primitive_target": list(event["primitive"]),
                          "vocabulary_size": len(vocabulary), "eligible_candidates": len(candidates),
                          "target_in_vocabulary": wanted is not None, "covered": actual is not None, "scores": {}}
                if candidates:
                    values = projection_scores(q, prior, candidates)
                else:
                    values = {}
                if actual is not None:
                    covered_count += 1
                draws = rng.integers(len(candidates), size=policy["random_draws"]) if actual is not None else None
                for score, ascending in scores_direction.items():
                    if actual is None:
                        record["scores"][score] = {"rank_best": None, "rank_worst": None, "percentile_loss": 1.0,
                                                    "conservative_top1": False, "conservative_top5": False, "conservative_top20": False,
                                                    "status": "OUT_OF_VOCABULARY_OR_NO_ELIGIBLE_TARGET"}
                        reference_draws[score] += 1
                        continue
                    best, worst = rank_interval(values[score], actual, ascending)
                    percentiles = all_percentiles(values[score], ascending)
                    record["scores"][score] = {
                        "rank_best": best, "rank_worst": worst, "score": float(values[score][actual]),
                        "percentile_loss": float(percentiles[actual]),
                        "conservative_top1": worst <= 1, "conservative_top5": worst <= 5, "conservative_top20": worst <= 20,
                        "all_candidates_tied": best == 1 and worst == len(candidates), "status": "RANKED"}
                    reference_draws[score] += percentiles[draws]
                run_records.append(record)
                prior.append(list(event["word"]))
                binary.append(parity(event["word"]))
            require(len(run_records) > 0, "empty trajectory has no prediction exposure")
            stats = {"seed": run["seed"], "events": len(run_records), "covered": covered_count,
                     "scores": {s: {"mean_percentile_loss_oov_penalized": sum(r["scores"][s]["percentile_loss"] for r in run_records)/len(run_records),
                                       "conditional_uniform_reference_interval": [float(x) for x in np.quantile(reference_draws[s]/len(run_records), [0.025, 0.5, 0.975])]}
                                for s in scores_direction}}
            per_run.append(stats)
            records.extend(run_records)
        score_summary = {}
        for score in scores_direction:
            rows = [record["scores"][score] for record in records]
            covered = [row for row in rows if row["status"] == "RANKED"]
            score_summary[score] = {
                "run_balanced_mean_percentile_loss_oov_penalized": sum(r["scores"][score]["mean_percentile_loss_oov_penalized"] for r in per_run)/len(per_run),
                "covered_mean_percentile_loss": None if not covered else sum(r["percentile_loss"] for r in covered)/len(covered),
                "conservative_top1_all_events": sum(r["conservative_top1"] for r in rows),
                "conservative_top5_all_events": sum(r["conservative_top5"] for r in rows),
                "conservative_top20_all_events": sum(r["conservative_top20"] for r in rows),
                "all_tied_covered_events": sum(r.get("all_candidates_tied", False) for r in covered)}
        arm_summary = {"events": len(records), "covered": sum(r["covered"] for r in records), "scores": score_summary,
                       "adaptive_improvement_over_static": {
                           s: score_summary["static_full_norm"]["run_balanced_mean_percentile_loss_oov_penalized"] - score_summary[s]["run_balanced_mean_percentile_loss_oov_penalized"]
                           for s in ("residual_norm", "axis_fractional_unlock")}}
        arms[arm] = {"summary": arm_summary, "per_run": per_run, "events": records}
        summary[arm] = arm_summary
    return {"schema": "curve302-followup-next-moves.v1", "status": "COMPLETE", "summary": summary, "arms": arms,
            "form_normalization_median_diagonal": normalization,
            "fixed_vocabulary_size": len(fixed), "fixed_vocabulary_sha256": sha256(canonical(fixed)).hexdigest(),
            "numeric_policy": {"rank_tie_atol": 1e-10, "rank_tie_rtol": 1e-9, "max_condition": 1e10, "random_seed": policy["random_seed"], "random_draws": policy["random_draws"]},
            "boundary": (
                "Vocabulary is never enlarged with a held-out target. Fixed bounded vectors are independent of acquisition outcomes; "
                "leave-one-run-out uses other runs plus axes. BOTH retain the known M31 form and are retrospective. "
                "Uniform eligible-candidate intervals are conditional reference distributions, NOT independent-trial p-values. "
                "These candidates were not necessarily offered by V3's actual chart schedule. Scores compare literal primitive vectors, "
                "not all equivalent rational lines modulo the prefix. Directional unlock necessarily ties when residual dimension is one. "
                "OOV events remain in the denominator, with loss 1 and no top-k hit."
            )}


# Runner and input provenance ---------------------------------------------
def runtime():
    import numpy as np
    import psutil
    return {"python": platform.python_version(), "numpy": np.__version__, "psutil": psutil.__version__,
            "platform": platform.platform(), "executable": sys.executable}


def find_results(repo, explicit=None):
    if explicit is not None:
        paths = [Path(explicit).resolve()]
    else:
        root = Path(repo)/"research/artifacts/generated-results/elliptic-curves"
        require(root.is_dir(), "no generated-results root; supply --results and --landscape explicitly")
        paths = []
        for directory, dirs, files in os.walk(root):
            dirs[:] = [d for d in dirs if d not in {"frozen_sources", "runtime", ".git", "partial_source_copy"}]
            if "trajectories.json" in files and "quotient-relations.json" in files:
                paths.append(Path(directory))
    valid = [p for p in paths if all((p/f).is_file() for f in ("closure-laws.json", "quotient-relations.json", "trajectories.json", "REPORT.json"))]
    require(bool(valid), "No complete source bundle found. --results must name a directory containing closure-laws.json, quotient-relations.json, trajectories.json and REPORT.json; no synthetic fallback.")
    by_hash = {}
    for p in valid:
        key = tuple(digest(p/f) for f in ("closure-laws.json", "quotient-relations.json", "trajectories.json", "REPORT.json"))
        by_hash.setdefault(key, []).append(p)
    require(len(by_hash) == 1, "Multiple different source bundles found; select --results explicitly: " + ", ".join(map(str, valid)))
    return sorted(valid, key=str)[0]


@contextmanager
def lock(folder):
    import fcntl
    with (Path(folder)/"LOCK").open("a+") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("another controller/check owns this output folder") from exc
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def prepare(folder, repo, results, landscape, policy):
    require(not folder.exists(), "output folder exists; choose resume/status/check or a new folder")
    source = find_results(repo, results)
    landscape = Path(landscape) if landscape else Path(repo)/"research/artifacts/generated-results/elliptic-curves/curve302_exceptional_subgroup_landscape_v1.json"
    require(landscape.is_file(), "landscape missing; supply --landscape")
    software = runtime()
    folder.mkdir(parents=True, exist_ok=False)
    with lock(folder):
        # Copy only to the newly owned folder. Frozen inputs isolate concurrent
        # commits elsewhere and prevent a TOCTOU change in research results.
        (folder/"inputs").mkdir()
        paths = {"landscape.json": landscape, **{f: source/f for f in ("closure-laws.json", "quotient-relations.json", "trajectories.json", "REPORT.json")}}
        origins = {}
        for name, path in paths.items():
            before = digest(path)
            shutil.copyfile(path, folder/"inputs"/name)
            require(before == digest(path) == digest(folder/"inputs"/name), "source changed during snapshot")
            origins[name] = {"path": str(path.resolve()), "sha256": before}
        data = validate_results(folder/"inputs")
        source_file = folder/"source/runner.py"
        source_file.parent.mkdir()
        shutil.copyfile(Path(__file__), source_file)
        try:
            commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True, stderr=subprocess.DEVNULL, timeout=10).strip()
        except (OSError, subprocess.SubprocessError):
            commit = None  # file hashes, not an inferred HEAD, bind the run
        plan = {"schema": SCHEMA, "policy": policy, "software": software, "source_commit_observed": commit,
                "dataset_kind": data["dataset_kind"],
                "source_sha256": digest(source_file), "inputs": origins, "source_acquisitions": sum(len(r["events"]) for r in data["runs"]),
                "stages": STAGES, "boundary": BOUNDARY}
        atomic(folder/"plan.json", plan)
        atomic(folder/"manifest.json", {"plan_sha256": digest(folder/"plan.json")})
    return plan


def guard(folder):
    plan = read(folder/"plan.json")
    require(plan.get("schema") == SCHEMA, "wrong plan schema")
    require(digest(folder/"plan.json") == read(folder/"manifest.json")["plan_sha256"], "plan changed")
    require(tuple(plan["stages"]) == STAGES, "stage list changed")
    require(digest(folder/"source/runner.py") == plan["source_sha256"], "frozen source changed")
    require(runtime() == plan["software"], "runtime changed; use the same Python environment")
    for name, row in plan["inputs"].items():
        require(name in {"landscape.json", "closure-laws.json", "quotient-relations.json", "trajectories.json", "REPORT.json"}, "unexpected input path")
        require(digest(folder/"inputs"/name) == row["sha256"], "frozen input changed: "+name)
    return plan


def kill_owned_group(proc):
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=2)
    except subprocess.TimeoutExpired:
        pass
    # Descendants may outlive their process-group leader.
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    proc.wait()


def supervise(command, folder, phase, seconds, memory):
    import psutil
    started = time.monotonic()
    outcome, peak, code = "UNKNOWN", 0, None
    supervisor_error = None
    proc = None
    with (phase/"worker.log").open("wb") as stream:
        try:
            proc = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                    start_new_session=True, env={**os.environ, "PYTHONUNBUFFERED": "1"})
            while True:
                code = proc.poll()
                if code is not None:
                    outcome = "COMPLETED" if code == 0 else "FAILED"
                    break
                if (folder/"STOP").exists():
                    outcome = "STOPPED_CENSORED"; break
                if time.monotonic()-started > seconds:
                    outcome = "TIMEOUT_CENSORED"; break
                try:
                    root = psutil.Process(proc.pid)
                    rss = 0
                    for p in [root, *root.children(recursive=True)]:
                        try:
                            rss += p.memory_info().rss
                        except psutil.NoSuchProcess:
                            pass
                    peak = max(peak, rss)
                    if rss > memory:
                        outcome = "MEMORY_CENSORED"; break
                except psutil.NoSuchProcess:
                    pass
                time.sleep(0.1)
        except KeyboardInterrupt:
            outcome = "INTERRUPTED_CENSORED"
        except Exception as error:
            outcome = "SUPERVISOR_ERROR"
            supervisor_error = repr(error)
        finally:
            if proc is not None:
                kill_owned_group(proc)
                code = proc.returncode
    row = {"outcome": outcome, "returncode": code, "seconds": time.monotonic()-started, "peak_rss_bytes": peak, "supervisor_error": supervisor_error}
    atomic(phase/"receipt.json", row)
    return row


def sealed_phase(folder, stage):
    phase = folder/"phases"/stage
    if not (phase/"seal.json").exists():
        return False
    seal = read(phase/"seal.json")
    require(seal["receipt_sha256"] == digest(phase/"receipt.json"), "receipt changed")
    receipt = read(phase/"receipt.json")
    require(receipt["outcome"] == "COMPLETED" and receipt["returncode"] == 0, "sealed failed phase")
    require(seal["output_sha256"] == digest(folder/f"{stage}.json"), "sealed output changed")
    require(read(folder/f"{stage}.json").get("status") == "COMPLETE", "sealed result incomplete")
    return True


def execute(folder):
    with lock(folder):
        plan = guard(folder)
        for stage in STAGES:
            if sealed_phase(folder, stage):
                continue
            phase = folder/"phases"/stage
            require(not phase.exists(), "failed, censored or unreceipted phase exists; preserve it and use a new folder")
            require(not (folder/"STOP").exists(), "STOP marker present; no phase launched")
            phase.mkdir(parents=True)
            atomic(phase/"started.json", {"stage": stage, "plan_sha256": digest(folder/"plan.json")})
            result = supervise([sys.executable, str(folder/"source/runner.py"), "_stage", "--folder", str(folder), "--stage", stage],
                               folder, phase, plan["policy"]["phase_seconds"], plan["policy"]["memory_mib"]*1024**2)
            guard(folder)
            if result["outcome"] != "COMPLETED":
                atomic(folder/"REPORT.json", {"status": "UNKNOWN_FAILED_OR_CENSORED", "stage": stage, "receipt": result})
                raise RuntimeError(f"{stage}: {result['outcome']}; inspect {phase/'worker.log'}")
            require(read(folder/f"{stage}.json").get("status") == "COMPLETE", "stage lacks complete output")
            atomic(phase/"seal.json", {"receipt_sha256": digest(phase/"receipt.json"), "output_sha256": digest(folder/f"{stage}.json")})
            print(f"SEALED {stage} ({result['seconds']:.2f}s)", flush=True)
        report = {"status": "COMPLETE_THREE_RETROSPECTIVE_FOLLOWUPS", "outputs": {s: digest(folder/f"{s}.json") for s in STAGES},
                  "source_acquisitions": plan["source_acquisitions"], "dataset_kind": plan["dataset_kind"], "boundary": BOUNDARY}
        atomic(folder/"REPORT.json", report)
        render_summary(folder)
        print(json.dumps(report, indent=2))


def render_summary(folder):
    p, f, m = (read(folder/f"{s}.json") for s in STAGES)
    text = ["# Curve302 closure follow-up\n", BOUNDARY+"\n",
            f"Empty-seed full threshold: {p['empty_full_threshold']}; strictly helpful singletons: {p['strictly_helpful_singletons']}.\n",
            "## Endpoint strict/local dimensions\n"]
    for run in f["runs"]:
        e = run["endpoint"]
        text.append(f"- {run['seed']}: quotient {e['quotient_dimension']}, strict {e['strict_dimension']}, local {e['local_image_dimension']}.")
    text.append("\n## Next-move tests\n")
    for name, arm in m["summary"].items():
        text.append(f"- {name}: coverage {arm['covered']}/{arm['events']}; adaptive-minus-static improvements {arm['adaptive_improvement_over_static']}.")
    text.append("\nSee JSON for tie intervals, OOV penalties, conditional reference intervals, and individual exact F2 witnesses.\n")
    (folder/"SUMMARY.md").write_text("\n".join(text), encoding="utf-8")


def check(folder):
    with lock(folder):
        plan = guard(folder)
        report = read(folder/"REPORT.json")
        require(report.get("status") == "COMPLETE_THREE_RETROSPECTIVE_FOLLOWUPS", "run not complete")
        for stage in STAGES:
            require(sealed_phase(folder, stage), "stage not sealed")
            require(report["outputs"][stage] == digest(folder/f"{stage}.json"), "report binding differs")
        # Keep every check's logs/receipts, including failed recomputations.
        checks = folder/"checks"
        checks.mkdir(exist_ok=True)
        target = Path(tempfile.mkdtemp(prefix="attempt-", dir=checks))
        for stage in STAGES:
            phase = target/stage
            phase.mkdir()
            result = supervise([sys.executable, str(folder/"source/runner.py"), "_stage", "--folder", str(folder),
                                "--stage", stage, "--check-output", str(target/f"{stage}.json")],
                               folder, phase, plan["policy"]["phase_seconds"], plan["policy"]["memory_mib"]*1024**2)
            require(result["outcome"] == "COMPLETED", f"recomputation failed: {stage}; evidence: {phase}")
            require(read(target/f"{stage}.json") == read(folder/f"{stage}.json"), f"recomputation mismatch: {stage}; evidence: {target}")
        atomic(target/"CHECK.json", {"status": "PASS_DETERMINISTIC_RECOMPUTATION", "source_report_sha256": digest(folder/"REPORT.json")})

    print("CHECK PASSED: deterministic recomputation of all three analyses (not a new EC proof)")


def status(folder):
    if not folder.exists():
        print("NOT_PREPARED"); return
    active = False
    try:
        with lock(folder):
            pass
    except RuntimeError:
        active = True
    stages = {}
    for stage in STAGES:
        p = folder/"phases"/stage
        if (p/"seal.json").exists():
            stages[stage] = "SEALED_UNCHECKED"
        elif (p/"receipt.json").exists():
            stages[stage] = read(p/"receipt.json")["outcome"]
        elif p.exists():
            stages[stage] = "RUNNING" if active else "INTERRUPTED_UNRECEIPTED"
        else:
            stages[stage] = "PENDING"
    print(json.dumps({"controller_lock_held": active, "stages": stages,
                      "report": read(folder/"REPORT.json") if (folder/"REPORT.json").exists() else None}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "run", "resume", "status", "check", "_stage"))
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--results", type=Path)
    parser.add_argument("--landscape", type=Path)
    parser.add_argument("--folder", type=Path)
    parser.add_argument("--phase-seconds", type=int, default=1800)
    parser.add_argument("--memory-mib", type=int, default=3072)
    parser.add_argument("--random-draws", type=int, default=1000)
    parser.add_argument("--profile", choices=("full", "smoke"), default="full")
    parser.add_argument("--stage", choices=STAGES, help=argparse.SUPPRESS)
    parser.add_argument("--check-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    repo = args.repo.resolve()
    folder = (args.folder or repo/"research/artifacts/local/elliptic-curves/curve302-closure-followup-v1").resolve()
    require(args.phase_seconds > 0 and 128 <= args.memory_mib <= 16384 and 10 <= args.random_draws <= 10000, "invalid finite limits")
    if args.command == "status":
        status(folder); return
    if args.command == "_stage":
        require(args.stage is not None, "missing stage")
        plan = guard(folder)
        if args.check_output is not None:
            destination = args.check_output.resolve()
            require(destination.is_relative_to(folder/"checks") and destination.name == f"{args.stage}.json",
                    "check output must be inside this run's checks directory")
        else:
            phase = folder/"phases"/args.stage
            require((phase/"started.json").is_file() and not (phase/"seal.json").exists(), "stage not started or already sealed")
            require(not (folder/f"{args.stage}.json").exists(), "unsealed stage output already exists; preserve it")
        data = validate_results(folder/"inputs")
        fn = {"persistence": persistence, "strict-filtration": filtration, "next-moves": next_moves}[args.stage]
        result = fn(data, plan["policy"])
        guard(folder)
        atomic(args.check_output or folder/f"{args.stage}.json", result)
        return
    if args.command == "check":
        check(folder); return
    if args.command in ("prepare", "run"):
        policy = {"phase_seconds": args.phase_seconds, "memory_mib": args.memory_mib,
                  "random_draws": args.random_draws, "random_seed": 3020911, "profile": args.profile}
        prepare(folder, repo, args.results, args.landscape, policy)
        print(f"PREPARED {folder}", flush=True)
        if args.command == "prepare":
            return
    execute(folder)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError) as error:
        print(f"FAIL CLOSED: {error}", file=sys.stderr)
        sys.exit(2)
