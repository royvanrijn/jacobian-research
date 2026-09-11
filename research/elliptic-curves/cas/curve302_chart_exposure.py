"""Exact offline helpers for Curve302 historical chart-exposure controls.

This module deliberately does *not* infer exposure from quotient height.  A chart
exposes a direction only when the historical audit/replay contains an explicit
14-dimensional quotient word for that chart result (or a normalized exposure
entry produced from such evidence).  Unsupported point-only schemas are kept as
schema diagnostics and never promoted to mathematical evidence.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
import math
import random
from pathlib import Path
from typing import Iterable, Sequence

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from curve302_short_core_controls import QuotientMap, integer, primitive, qnorm, require


def saturation_index(rows, n: int) -> int:
    """Index of the row lattice in its saturation inside Z^n (exact)."""
    rows = [tuple(integer(x) for x in row) for row in rows]
    require(all(len(row) == n for row in rows), "saturation row width mismatch")
    if not rows:
        return 1
    A = Matrix(rows)
    rank = A.rank()
    if rank == 0:
        return 1
    D = smith_normal_form(A, domain=ZZ)
    diagonal = []
    for i in range(min(D.rows, D.cols)):
        d = abs(integer(D[i, i]))
        if d:
            diagonal.append(d)
    require(len(diagonal) == rank, "Smith rank mismatch")
    out = 1
    for d in diagonal:
        out *= d
    return out


WORD_KEYS = (
    "primitive_quotient_word", "quotient_word", "candidate_word",
    "direction_word", "represented_word", "extension_word",
)
HEIGHT_KEYS = (
    "parameter_height", "slope_height", "rational_height", "search_height",
)
BOUND_KEYS = ("search_bound", "height_bound", "bound", "height")
EXPOSURE_CONTAINER_KEYS = (
    "exposures", "points", "cloud", "hits", "found", "found_points",
    "rational_points", "replay_points", "results",
)
CHART_CONTAINER_KEYS = ("charts", "chart_results", "attempts", "calls")
SEED_KEYS = ("seed", "seed_id", "seed_name", "direction", "direction_id")
EPOCH_KEYS = ("epoch", "stage", "stage_index", "round")
CHART_ID_KEYS = ("chart_id", "id", "chart", "name")
ORDER_KEYS = ("order", "chart_order", "call", "call_index", "index")
SCORE_BAND_KEYS = ("score_band", "band", "priority_band", "bucket")
COEFFICIENT_KEYS = (
    "quartic", "quartic_coefficients", "coefficients", "model_coefficients",
    "binary_quartic", "pointed_quartic",
)
COMPLETE_KEYS = ("exhaustive", "complete", "search_complete", "replay_complete")


class UnsupportedChartSchema(ValueError):
    pass


def read_json(path: Path):
    def unique(pairs):
        out = {}
        for key, value in pairs:
            require(key not in out, f"duplicate JSON key {key} in {path}")
            out[key] = value
        return out
    def reject(value):
        raise ValueError(f"nonfinite JSON value {value} in {path}")
    return json.loads(path.read_text(), object_pairs_hook=unique, parse_constant=reject)


def as_fraction(value) -> F:
    if isinstance(value, F):
        return value
    require(not isinstance(value, bool), "boolean is not a rational")
    if isinstance(value, int):
        return F(value)
    return F(str(value))


def canonical_word(value, n: int):
    if not isinstance(value, (list, tuple)) or len(value) != n:
        return None
    try:
        word = tuple(integer(v) for v in value)
    except Exception:
        return None
    if not any(word):
        return None
    return primitive(word)


def explicit_word(record, n: int):
    """Return one explicit quotient direction, rejecting contradictory aliases."""
    if not isinstance(record, dict):
        return None
    found = []
    for key in WORD_KEYS:
        if key in record:
            word = canonical_word(record[key], n)
            if word is not None:
                found.append((key, word))
    if not found:
        return None
    unique = {word for _, word in found}
    require(len(unique) == 1, f"contradictory explicit quotient words: {[k for k,_ in found]}")
    return found[0][1]


def first_scalar(record, keys):
    if not isinstance(record, dict):
        return None
    for key in keys:
        if key in record and isinstance(record[key], (str, int)) and not isinstance(record[key], bool):
            return record[key]
    return None


def parameter_height(record):
    raw = first_scalar(record, HEIGHT_KEYS)
    if raw is None:
        # Some replays store the actual rational parameter, from which the
        # standard projective height is exact and unambiguous.
        for key in ("parameter", "slope", "u", "t"):
            if isinstance(record, dict) and key in record and isinstance(record[key], (str, int)):
                try:
                    q = as_fraction(record[key])
                except Exception:
                    continue
                return max(abs(q.numerator), q.denominator)
        return None
    q = as_fraction(raw)
    require(q.denominator == 1 and q >= 0, "parameter height must be a nonnegative integer")
    return q.numerator


def search_bound(record):
    raw = first_scalar(record, BOUND_KEYS)
    if raw is None:
        return None
    q = as_fraction(raw)
    require(q.denominator == 1 and q >= 0, "search bound must be a nonnegative integer")
    return q.numerator


def bool_complete(record):
    if not isinstance(record, dict):
        return None
    values = [record[k] for k in COMPLETE_KEYS if k in record and isinstance(record[k], bool)]
    if not values:
        return None
    require(len(set(values)) == 1, "contradictory complete/exhaustive flags")
    return values[0]


def max_integer_bits(value):
    """Maximum numerator/denominator bit length in a coefficient payload."""
    best = 0
    if isinstance(value, bool) or value is None:
        return 0
    if isinstance(value, int):
        return abs(value).bit_length()
    if isinstance(value, str):
        try:
            q = F(value)
        except Exception:
            return 0
        return max(abs(q.numerator).bit_length(), q.denominator.bit_length())
    if isinstance(value, (list, tuple)):
        return max((max_integer_bits(v) for v in value), default=0)
    if isinstance(value, dict):
        return max((max_integer_bits(v) for v in value.values()), default=0)
    return 0


def quartic_bits(record):
    if not isinstance(record, dict):
        return None
    values = [max_integer_bits(record[k]) for k in COEFFICIENT_KEYS if k in record]
    values = [v for v in values if v]
    return max(values) if values else None


def score_band(record):
    value = first_scalar(record, SCORE_BAND_KEYS)
    return None if value is None else str(value)


def _walk_exposure_entries(value, n: int, path: str, depth: int = 0):
    """Yield explicit quotient-word records only from recognized result containers."""
    if depth > 8:
        return
    if isinstance(value, dict):
        word = explicit_word(value, n)
        if word is not None:
            yield path, value, word
        for key in EXPOSURE_CONTAINER_KEYS:
            if key in value:
                child = value[key]
                if isinstance(child, list):
                    for i, item in enumerate(child):
                        yield from _walk_exposure_entries(item, n, f"{path}.{key}[{i}]", depth + 1)
                elif isinstance(child, dict):
                    yield from _walk_exposure_entries(child, n, f"{path}.{key}", depth + 1)
    elif isinstance(value, list):
        for i, item in enumerate(value):
            yield from _walk_exposure_entries(item, n, f"{path}[{i}]", depth + 1)


def looks_like_point_only(record):
    if not isinstance(record, dict):
        return False
    keys = set(record)
    return ({"x", "y"} <= keys or "point" in keys or "elliptic_point" in keys) and not any(k in keys for k in WORD_KEYS)


def normalize_chart(record, order: int, n: int, *, default_bound=None, source_path="?"):
    require(isinstance(record, dict), "chart record is not an object")
    chart_id = next((str(record[k]) for k in CHART_ID_KEYS if k in record and isinstance(record[k], (str, int))),
                    f"chart-{order:06d}")
    explicit_order = first_scalar(record, ORDER_KEYS)
    if explicit_order is not None:
        try:
            explicit_order = integer(explicit_order)
        except Exception:
            explicit_order = None
    chart_order = order if explicit_order is None else explicit_order
    bound = search_bound(record)
    if bound is None:
        bound = default_bound
    complete = bool_complete(record)
    bits = quartic_bits(record)
    band = score_band(record)

    rows = []
    seen = {}
    point_only = 0
    for exp_order, (path, entry, word) in enumerate(_walk_exposure_entries(record, n, source_path)):
        if path == source_path and entry is record:
            # A chart-level nominal word is not automatically a discovered/exposed
            # point. Only accept it when an explicit result semantic is present.
            semantic = str(record.get("representation_kind", record.get("kind", ""))).lower()
            if semantic not in {"exposure", "exposed_direction", "search_hit", "found_point", "replay_point"}:
                continue
        h = parameter_height(entry)
        within = entry.get("within_bound") if isinstance(entry.get("within_bound"), bool) else None
        if within is None and h is not None and bound is not None:
            within = h <= bound
        if within is False:
            continue
        prev = seen.get(word)
        row = {"word": word, "parameter_height": h, "exposure_order": exp_order,
               "source_path": path, "within_bound": within}
        if prev is None or (h is not None and (prev["parameter_height"] is None or h < prev["parameter_height"])):
            seen[word] = row
    # Count unsupported point-only result records for diagnostics. This is kept
    # separate from evidence and prevents silent under-interpretation.
    def count_points(obj, depth=0):
        nonlocal point_only
        if depth > 8:
            return
        if isinstance(obj, dict):
            if looks_like_point_only(obj):
                point_only += 1
            for key in EXPOSURE_CONTAINER_KEYS:
                if key in obj:
                    count_points(obj[key], depth + 1)
        elif isinstance(obj, list):
            for item in obj:
                count_points(item, depth + 1)
    count_points(record)
    rows = sorted(seen.values(), key=lambda row: (row["exposure_order"], row["word"]))
    return {"chart_id": chart_id, "order": chart_order, "score_band": band,
            "quartic_coefficient_bits": bits, "search_bound": bound,
            "complete_flag": complete, "exposures": rows,
            "unsupported_point_only_records": point_only}


def _context_seed(record, names, fallback=None):
    if isinstance(record, dict):
        for key in SEED_KEYS:
            value = record.get(key)
            if value in names:
                return value
    return fallback


def _context_epoch(record, fallback=None):
    if isinstance(record, dict):
        for key in EPOCH_KEYS:
            value = record.get(key)
            if isinstance(value, int) and not isinstance(value, bool) and value >= 0:
                return value
            if isinstance(value, str) and value.isdigit():
                return int(value)
    return fallback


def _seed_from_path(path: Path, names):
    text = str(path).lower().replace("_", "-")
    matches = [name for name in names if name.lower() in text]
    return matches[0] if len(matches) == 1 else None


def _epoch_from_path(path: Path):
    import re
    text = path.name.lower()
    patterns = (r"epoch[-_]?([0-9]+)", r"stage[-_]?([0-9]+)", r"round[-_]?([0-9]+)")
    found = []
    for pattern in patterns:
        found.extend(int(m) for m in re.findall(pattern, text))
    return found[0] if len(set(found)) == 1 else None


def discover_chart_blocks(value, names, *, seed=None, epoch=None, path="$", depth=0):
    """Find candidate historical stage chart arrays with inherited context."""
    if depth > 12:
        return
    if isinstance(value, dict):
        seed = _context_seed(value, names, seed)
        epoch = _context_epoch(value, epoch)
        default_bound = search_bound(value)
        for key in CHART_CONTAINER_KEYS:
            charts = value.get(key)
            if isinstance(charts, list) and charts and all(isinstance(c, dict) for c in charts):
                yield {"seed": seed, "epoch": epoch, "path": f"{path}.{key}",
                       "default_bound": default_bound, "charts": charts}
        for key, child in value.items():
            if key in CHART_CONTAINER_KEYS:
                continue
            if isinstance(child, (dict, list)):
                yield from discover_chart_blocks(child, names, seed=seed, epoch=epoch,
                                                 path=f"{path}.{key}", depth=depth + 1)
    elif isinstance(value, list):
        for i, child in enumerate(value):
            if isinstance(child, (dict, list)):
                yield from discover_chart_blocks(child, names, seed=seed, epoch=epoch,
                                                 path=f"{path}[{i}]", depth=depth + 1)


def normalize_explicit_ledger(payload, names, n):
    require(payload.get("schema") == "curve302-chart-exposure-ledger.v1", "wrong normalized ledger schema")
    require(tuple(payload.get("direction_ids", ())) == tuple(names), "normalized ledger direction roster changed")
    runs = []
    for raw_run in payload.get("runs", []):
        seed = raw_run.get("seed")
        require(seed in names, "unknown normalized ledger seed")
        stages = []
        for raw_stage in raw_run.get("stages", []):
            epoch = integer(raw_stage["epoch"])
            charts = [normalize_chart(c, i, n, default_bound=raw_stage.get("search_bound"),
                                      source_path=f"normalized:{seed}:{epoch}:{i}")
                      for i, c in enumerate(raw_stage.get("charts", []))]
            stages.append({"epoch": epoch, "charts": charts, "source_file": "normalized-ledger",
                           "source_path": f"runs[{seed}].stages[{epoch}]"})
        runs.append({"seed": seed, "stages": stages})
    return {"schema": "curve302-chart-exposure-ledger.v1", "direction_ids": list(names), "runs": runs,
            "source_files": [], "normalization": "explicit-ledger"}


def normalize_raw_tree(raw_root: Path, names: Sequence[str], expected_epochs: dict[str, set[int]], n: int):
    """Best-effort strict adapter for historical JSON audits.

    It accepts only explicit quotient-word result records.  Competing candidate
    chart arrays for the same seed/epoch are resolved by *strict superset* chart
    count, then exact canonical equality; otherwise the source is ambiguous.
    """
    raw_root = Path(raw_root).resolve()
    require(raw_root.is_dir(), f"raw root is not a directory: {raw_root}")
    json_files = sorted(p for p in raw_root.rglob("*.json") if p.is_file())
    require(json_files, "raw root contains no JSON files")
    candidates = defaultdict(list)
    probe = {"files_scanned": len(json_files), "candidate_blocks": 0, "resolved_blocks": 0,
             "unresolved_context_blocks": [], "unsupported_point_only_records": 0, "sample_key_paths": []}
    used_files = set()
    for file in json_files:
        try:
            payload = read_json(file)
        except Exception:
            continue
        if isinstance(payload, dict) and payload.get("schema") == "curve302-chart-exposure-ledger.v1":
            ledger = normalize_explicit_ledger(payload, names, n)
            ledger["source_files"] = [str(file)]
            ledger["probe"] = probe
            return ledger
        file_seed = _seed_from_path(file.relative_to(raw_root), names)
        file_epoch = _epoch_from_path(file)
        for block in discover_chart_blocks(payload, names, seed=file_seed, epoch=file_epoch):
            probe["candidate_blocks"] += 1
            seed, epoch = block["seed"], block["epoch"]
            if seed not in names or epoch is None or epoch not in expected_epochs.get(seed, set()):
                if len(probe["unresolved_context_blocks"]) < 50:
                    probe["unresolved_context_blocks"].append({"file": str(file), "path": block["path"],
                                                                 "seed": seed, "epoch": epoch,
                                                                 "charts": len(block["charts"])})
                continue
            normalized = [normalize_chart(c, i, n, default_bound=block["default_bound"],
                                          source_path=f"{file}:{block['path']}[{i}]")
                          for i, c in enumerate(block["charts"])]
            unsupported = sum(c["unsupported_point_only_records"] for c in normalized)
            explicit = sum(len(c["exposures"]) for c in normalized)
            canonical = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
            candidates[(seed, epoch)].append({"file": file, "path": block["path"], "charts": normalized,
                                               "chart_count": len(normalized), "explicit": explicit,
                                               "unsupported": unsupported, "canonical": canonical})
    resolved = {}
    for key, rows in candidates.items():
        # Prefer the most information-rich explicit replay, then largest chart
        # list.  Ties must be byte-for-byte equivalent after normalization.
        best_score = max((row["explicit"], row["chart_count"]) for row in rows)
        top = [row for row in rows if (row["explicit"], row["chart_count"]) == best_score]
        canon = {row["canonical"] for row in top}
        require(len(canon) == 1, f"ambiguous historical chart blocks for {key}; pass a normalized --ledger")
        chosen = top[0]
        resolved[key] = chosen
        used_files.add(chosen["file"])
        probe["unsupported_point_only_records"] += chosen["unsupported"]
        probe["resolved_blocks"] += 1
    missing = [(seed, epoch) for seed, epochs in expected_epochs.items() for epoch in sorted(epochs)
               if (seed, epoch) not in resolved]
    require(not missing, f"missing historical chart blocks for {missing[:10]}{'...' if len(missing)>10 else ''}")
    # Exposure-free stages are not evidence for the proposed experiments.  Keep
    # the schema probe precise rather than treating zero as a mathematical zero.
    no_explicit = [(seed, epoch) for (seed, epoch), row in resolved.items() if row["explicit"] == 0]
    if no_explicit:
        raise UnsupportedChartSchema(
            "historical chart blocks were found but contain no explicit 14D quotient-word exposure records "
            f"for {len(no_explicit)} stages (first {no_explicit[:5]}). Use the schema probe / a narrow adapter; "
            "point-only records are not silently interpreted.")
    runs = []
    for seed in names:
        stages = []
        for epoch in sorted(expected_epochs[seed]):
            row = resolved[(seed, epoch)]
            stages.append({"epoch": epoch, "charts": row["charts"], "source_file": str(row["file"]),
                           "source_path": row["path"]})
        runs.append({"seed": seed, "stages": stages})
    return {"schema": "curve302-chart-exposure-ledger.v1", "direction_ids": list(names), "runs": runs,
            "source_files": [str(p) for p in sorted(used_files)], "normalization": "raw-audit-explicit-word-v1",
            "probe": probe}


def validate_ledger(ledger, names, runs, *, expected_total_charts=None):
    require(ledger.get("schema") == "curve302-chart-exposure-ledger.v1", "wrong ledger schema")
    require(tuple(ledger.get("direction_ids", ())) == tuple(names), "ledger direction roster mismatch")
    by_seed = {row["seed"]: row for row in ledger.get("runs", [])}
    require(set(by_seed) == set(names), "ledger seed roster incomplete")
    total_charts = total_exposures = unsupported = 0
    complete_true = complete_false = complete_unknown = 0
    for run in runs:
        seed = run["seed"]
        stages = by_seed[seed]["stages"]
        expected_epochs = [integer(e) for e in run.get("epochs", range(len(run["events"])))]
        require([integer(s["epoch"]) for s in stages] == expected_epochs,
                f"ledger epochs differ for {seed}")
        for stage in stages:
            charts = stage["charts"]
            require(charts, f"empty chart list for {seed} epoch {stage['epoch']}")
            orders = [integer(c["order"]) for c in charts]
            # Historical list order is authoritative; explicit order fields may
            # have gaps but must be unique.
            require(len(set(orders)) == len(orders), "duplicate chart order")
            total_charts += len(charts)
            for chart in charts:
                flag = chart.get("complete_flag")
                if flag is True: complete_true += 1
                elif flag is False: complete_false += 1
                else: complete_unknown += 1
                unsupported += integer(chart.get("unsupported_point_only_records", 0))
                seen = set()
                for exp in chart["exposures"]:
                    word = tuple(integer(v) for v in exp["word"])
                    require(len(word) == len(names) and primitive(word) == word, "noncanonical exposure word")
                    require(word not in seen, "duplicate chart exposure direction")
                    seen.add(word); total_exposures += 1
                    h = exp.get("parameter_height")
                    if h is not None:
                        require(integer(h) >= 0, "negative parameter height")
    if expected_total_charts is not None:
        require(total_charts == expected_total_charts,
                f"historical chart count mismatch: ledger {total_charts}, trajectory audit {expected_total_charts}")
    return {"charts": total_charts, "exposures": total_exposures,
            "complete_charts": complete_true, "incomplete_charts": complete_false,
            "unknown_completeness_charts": complete_unknown,
            "unsupported_point_only_records": unsupported}


def vector_key(v):
    return ",".join(str(int(x)) for x in v)


def parse_vector_key(text):
    return tuple(int(x) for x in text.split(","))


def rank_band(index_zero: int):
    rank = index_zero + 1
    lo = 1
    while lo * 10 <= rank:
        lo *= 10
    return lo, lo * 10


def build_candidate_population(vocab, runs, *, static_limit=1000):
    require(1 <= static_limit <= len(vocab), "invalid static candidate limit")
    vectors = [vocab.vector(i) for i in range(static_limit)]
    positions = {v: i for i, v in enumerate(vectors)}
    for run in runs:
        for event in run["events"]:
            v = event["primitive"]
            if v not in positions:
                positions[v] = vocab.positions[v]
                vectors.append(v)
    vectors = sorted(set(vectors), key=lambda v: (positions.get(v, vocab.positions.get(v, 10**30)), v))
    final_positions = {}
    for v in vectors:
        if v in positions:
            final_positions[v] = positions[v]
        else:
            require(v in vocab.positions, "candidate position unavailable")
            final_positions[v] = vocab.positions[v]
    return vectors, final_positions


def stage_prefixes(runs, n):
    result = {}
    for run in runs:
        rows = [tuple(int(i == run["seed_index"]) for i in range(n))]
        for epoch, event in enumerate(run["events"]):
            result[(run["seed"], epoch)] = tuple(rows)
            rows.append(tuple(event["word"]))
    return result


def rationally_unknown(prefix, candidate, n):
    state = QuotientMap(n)
    for row in prefix:
        kind, _ = state.admission(row)
        require(kind == "PRIMITIVE", "historical prefix extension is not primitive")
        state = state.extend(row)
    return any(state.image(candidate))


def integer_contains(rows: Sequence[Sequence[int]], vector: Sequence[int]):
    if not rows:
        return not any(vector)
    A = Matrix([[int(x) for x in row] for row in rows]).T
    b = Matrix([int(x) for x in vector])
    try:
        sol, params = A.gauss_jordan_solve(b)
    except ValueError:
        return False
    if params.rows:
        # With independent prefix rows, the coefficient solution is unique.
        # Keep this conservative for malformed synthetic inputs.
        return False
    return all(getattr(x, "q", 1) == 1 for x in sol)


def contains_core_integrally(rows, core):
    return all(integer_contains(rows, v) for v in core)


def candidate_metrics_for_stage(charts, candidates, positions, prefix, actual, vocab, form, n):
    state = QuotientMap(n)
    for row in prefix:
        kind, _ = state.admission(row)
        require(kind == "PRIMITIVE", "historical prefix extension is not primitive")
        state = state.extend(row)
    candidate_set = {v for v in candidates if any(state.image(v))}
    exposure = {v: {"count": 0, "first_order": None, "min_parameter_height": None,
                    "min_quartic_bits": None, "chart_ids": []} for v in candidate_set}
    chart_rows = []
    for list_order, chart in enumerate(charts):
        exposed = []
        for item in chart["exposures"]:
            word = tuple(item["word"])
            if word not in candidate_set:
                continue
            row = exposure[word]
            row["count"] += 1
            if row["first_order"] is None:
                row["first_order"] = list_order
            h = item.get("parameter_height")
            if h is not None:
                h = integer(h)
                row["min_parameter_height"] = h if row["min_parameter_height"] is None else min(row["min_parameter_height"], h)
            bits = chart.get("quartic_coefficient_bits")
            if bits is not None:
                bits = integer(bits)
                row["min_quartic_bits"] = bits if row["min_quartic_bits"] is None else min(row["min_quartic_bits"], bits)
            row["chart_ids"].append(chart["chart_id"])
            exposed.append({"word": word, "exposure_order": integer(item.get("exposure_order", 0))})
        exposed.sort(key=lambda x: (x["exposure_order"], positions[x["word"]], x["word"]))
        chart_rows.append({"chart_id": chart["chart_id"], "original_order": list_order,
                           "score_band": chart.get("score_band"),
                           "quartic_coefficient_bits": chart.get("quartic_coefficient_bits"),
                           "complete_flag": chart.get("complete_flag"),
                           "exposed": [x["word"] for x in exposed]})
    actual = primitive(actual)
    require(actual in candidate_set, "actual acquisition absent from candidate population or already known")
    ai = positions[actual]
    band = vocab.band(ai)
    lo_i, hi_i = band["index_start"], band["index_stop"]
    controls = [v for v in candidate_set if v != actual and lo_i <= positions[v] < hi_i]
    # If decade-matched population is empty, retain an explicit UNKNOWN rather
    # than widening after looking at the outcome.
    def summary(v):
        row = exposure[v]
        return {"word": v, "static_rank": positions[v] + 1, "static_norm": str(qnorm(form, v)), **row}
    actual_row = summary(actual)
    control_rows = [summary(v) for v in controls]
    coverage_complete = all(chart.get("complete_flag") is True for chart in charts)
    return {"actual": actual_row, "controls": control_rows,
            "rank_band": [band["nominal_rank_low"], band["nominal_rank_high"]],
            "rank_band_index_range": [lo_i, hi_i],
            "candidate_count": len(candidate_set), "chart_rows": chart_rows,
            "coverage_complete": coverage_complete,
            "complete_chart_count": sum(chart.get("complete_flag") is True for chart in charts),
            "incomplete_or_unknown_chart_count": sum(chart.get("complete_flag") is not True for chart in charts),
            "exposed_candidate_count": sum(row["count"] > 0 for row in exposure.values())}


def tied_percentile(actual_value, controls, *, higher_better: bool):
    values = [v for v in controls if v is not None]
    if actual_value is None or not values:
        return None
    if higher_better:
        better = sum(v > actual_value for v in values); equal = sum(v == actual_value for v in values)
    else:
        better = sum(v < actual_value for v in values); equal = sum(v == actual_value for v in values)
    # Midrank percentile: 1.0 means best, 0.0 worst.
    return str(F(len(values) - better - F(equal, 2), len(values)))


def multiplicity_comparison(stage_rows):
    records = []
    for stage in stage_rows:
        a = stage["actual"]; controls = stage["controls"]
        complete = stage.get("coverage_complete") is True
        row = {"seed": stage["seed"], "epoch": stage["epoch"], "actual_word": a["word"],
               "actual_static_rank": a["static_rank"], "rank_band": stage["rank_band"],
               "control_count": len(controls), "coverage_complete": complete,
               "actual_exposure_count": a["count"],
               "actual_min_parameter_height": a["min_parameter_height"],
               "actual_min_quartic_bits": a["min_quartic_bits"]}
        if complete:
            row["status"] = "COMPLETE_CHART_COVERAGE"
            row["exposure_count_percentile"] = tied_percentile(a["count"], [c["count"] for c in controls], higher_better=True)
            row["parameter_height_percentile"] = tied_percentile(a["min_parameter_height"], [c["min_parameter_height"] for c in controls], higher_better=False)
            row["quartic_bits_percentile"] = tied_percentile(a["min_quartic_bits"], [c["min_quartic_bits"] for c in controls], higher_better=False)
        else:
            row["status"] = "UNKNOWN_INCOMPLETE_CHART_COVERAGE"
            row["exposure_count_percentile"] = None
            row["parameter_height_percentile"] = None
            row["quartic_bits_percentile"] = None
        row["actual_is_exposed"] = a["count"] > 0
        row["control_exposed"] = sum(c["count"] > 0 for c in controls)
        records.append(row)
    usable = [r for r in records if r["control_count"] and r["coverage_complete"]]
    exposed_actual = sum(r["actual_is_exposed"] for r in records)
    return {"status": "PASS_RECORDED_EXPOSURE_MULTIPLICITY_WITH_COMPLETENESS_GATE",
            "stages": records,
            "summary": {"stage_count": len(records), "complete_coverage_stages": sum(r["coverage_complete"] for r in records),
                        "incomplete_or_unknown_coverage_stages": sum(not r["coverage_complete"] for r in records),
                        "usable_rank_band_stages": len(usable),
                        "actual_exposed_stages": exposed_actual,
                        "actual_not_exposed_stages": len(records) - exposed_actual,
                        "median_exposure_percentile": median_fraction([r["exposure_count_percentile"] for r in usable]),
                        "median_parameter_height_percentile": median_fraction([r["parameter_height_percentile"] for r in usable]),
                        "median_quartic_bits_percentile": median_fraction([r["quartic_bits_percentile"] for r in usable])},
            "boundary": "Recorded exposure counts are always retained. Norm-matched absence/percentile claims are promoted only for stages whose chart searches/replays are explicitly complete. Percentiles are retrospective descriptive ranks, not p-values or prospective success probabilities."}


def median_fraction(values):
    vals = sorted(F(v) for v in values if v is not None)
    if not vals:
        return None
    m = len(vals) // 2
    answer = vals[m] if len(vals) % 2 else (vals[m-1] + vals[m]) / 2
    return str(answer)


def first_exposed(chart_rows, order):
    # chart_rows were already intersected with the stage's rationally unknown
    # candidate population, so no repeated span test is needed here.
    for chart_index in order:
        chart = chart_rows[chart_index]
        if chart["exposed"]:
            return {"chart_index": chart_index, "chart_id": chart["chart_id"],
                    "word": chart["exposed"][0]}
    return None


def _rng(label):
    return random.Random(int.from_bytes(sha256(label.encode()).digest(), "big"))


def counterfactual_stage(stage, prefix_rows, cores_by_dim, names, positions, *, random_orders=256, master_seed="302-chart-exposure-v1"):
    charts = stage["chart_rows"]
    n = len(names); indices = list(range(len(charts)))
    if stage.get("coverage_complete") is not True:
        policies = []
        for name in ("original", "reverse", "quartic_bits", "random", "random_within_score_band"):
            policies.append({"policy": name, "status": "UNKNOWN_INCOMPLETE_CHART_COVERAGE"})
        return {"seed": stage["seed"], "epoch": stage["epoch"],
                "next_dimension": len(prefix_rows) + 1, "chart_count": len(charts),
                "score_bands_available": all(c.get("score_band") is not None for c in charts),
                "coverage_complete": False, "policies": policies}
    original = indices
    reverse = list(reversed(indices))
    coeff = sorted(indices, key=lambda i: (charts[i]["quartic_coefficient_bits"] is None,
                                           charts[i]["quartic_coefficient_bits"] if charts[i]["quartic_coefficient_bits"] is not None else 10**100,
                                           i))
    orders = {"original": [original], "reverse": [reverse], "quartic_bits": [coeff]}
    rng = _rng(f"{master_seed}|{stage['seed']}|{stage['epoch']}|random")
    random_set = []
    for _ in range(random_orders):
        o = indices[:]; rng.shuffle(o); random_set.append(o)
    orders["random"] = random_set
    bands_available = all(charts[i].get("score_band") is not None for i in indices)
    if bands_available:
        grouped = []
        for band in []:
            pass
        first_seen = []
        groups = defaultdict(list)
        for i in indices:
            b = charts[i]["score_band"]
            if b not in groups:
                first_seen.append(b)
            groups[b].append(i)
        rngb = _rng(f"{master_seed}|{stage['seed']}|{stage['epoch']}|within-band")
        band_orders = []
        for _ in range(random_orders):
            o = []
            for b in first_seen:
                g = groups[b][:]; rngb.shuffle(g); o.extend(g)
            band_orders.append(o)
        orders["random_within_score_band"] = band_orders
    else:
        orders["random_within_score_band"] = None

    next_dim = len(prefix_rows) + 1
    core = cores_by_dim.get(next_dim, ())
    axes = {"L2": tuple(int(i == 1) for i in range(n)),
            "L1": tuple(int(i == 0) for i in range(n)),
            "L4": tuple(int(i == 3) for i in range(n))}
    actual = tuple(stage["actual"]["word"])

    # Exact lattice properties depend only on the stage prefix and exposed word,
    # not on the chart ordering.  Cache them once per distinct exposed word so
    # randomized counterfactuals do not recompute Smith/Hermite-style tests.
    properties = {}
    for word in {tuple(v) for chart in charts for v in chart["exposed"]}:
        rows = (*prefix_rows, word)
        properties[word] = {
            "primitive_extension": saturation_index(rows, len(word)) == 1,
            "contains_observed_next_core_integrally": contains_core_integrally(rows, core),
            "contains_axes_integrally": {key: integer_contains(rows, axis) for key, axis in axes.items()},
        }

    def evaluate(name, block):
        if block is None:
            return {"policy": name, "status": "UNKNOWN_SCORE_BANDS_UNAVAILABLE"}
        outcomes = []
        for o in block:
            hit = first_exposed(charts, o)
            if hit is None:
                outcomes.append({"status": "NO_EXPOSED_UNKNOWN"})
                continue
            word = tuple(hit["word"]); props = properties[word]
            outcomes.append({"status": "GAIN", "word": word, "chart_id": hit["chart_id"],
                             "matches_actual": primitive(word) == primitive(actual),
                             **props})
        gains = [r for r in outcomes if r["status"] == "GAIN"]
        return {"policy": name, "draws": len(outcomes), "gains": len(gains),
                "no_exposed_unknown": len(outcomes)-len(gains),
                "matches_actual": sum(r["matches_actual"] for r in gains),
                "primitive_extensions": sum(r["primitive_extension"] for r in gains),
                "contains_observed_next_core_integrally": sum(r["contains_observed_next_core_integrally"] for r in gains),
                "contains_axes_integrally": {key: sum(r["contains_axes_integrally"][key] for r in gains) for key in axes},
                "outcomes": outcomes if len(outcomes) <= 2 else None}

    return {"seed": stage["seed"], "epoch": stage["epoch"], "next_dimension": next_dim,
            "chart_count": len(charts), "score_bands_available": bands_available,
            "policies": [evaluate(name, block) for name, block in orders.items()]}


def summarize_counterfactuals(rows):
    policies = defaultdict(lambda: Counter())
    unknown = defaultdict(int)
    for row in rows:
        for policy in row["policies"]:
            name = policy["policy"]
            if policy.get("status", "").startswith("UNKNOWN"):
                unknown[name] += 1; continue
            for key in ("draws", "gains", "no_exposed_unknown", "matches_actual", "primitive_extensions",
                        "contains_observed_next_core_integrally"):
                policies[name][key] += integer(policy[key])
            for axis, value in policy["contains_axes_integrally"].items():
                policies[name]["axis_"+axis] += integer(value)
    summary = []
    for name in sorted(set(policies) | set(unknown)):
        row = {"policy": name, **dict(policies[name]), "unknown_stages": unknown[name]}
        summary.append(row)
    return {"status": "PASS_STAGE_LOCAL_CHART_ORDER_CONTROLS", "stages": rows, "summary": summary,
            "boundary": "Each historical stage is reordered independently and only stages with explicitly complete chart coverage are promoted. A counterfactual gain is never propagated into later historical stages; these controls separate recorded atlas exposure from within-stage scheduling only."}


def schema_probe(raw_root: Path, names: Sequence[str], *, max_files=10000, max_samples=80):
    """Non-promoting inventory to make a schema compatibility failure actionable."""
    raw_root = Path(raw_root).resolve()
    files = sorted(p for p in raw_root.rglob("*.json") if p.is_file())
    out = {"raw_root": str(raw_root), "json_files": len(files), "files_examined": 0,
           "chart_blocks": 0, "chart_records": 0, "explicit_word_records": 0,
           "point_only_records": 0, "samples": [], "top_level_schemas": Counter(),
           "chart_container_keys": Counter(), "exposure_container_keys": Counter(),
           "record_keys": Counter()}
    for file in files[:max_files]:
        try:
            payload = read_json(file)
        except Exception as exc:
            if len(out["samples"]) < max_samples:
                out["samples"].append({"file": str(file), "kind": "parse_error", "detail": str(exc)})
            continue
        out["files_examined"] += 1
        if isinstance(payload, dict) and "schema" in payload:
            out["top_level_schemas"][str(payload["schema"])] += 1
        file_seed = _seed_from_path(file.relative_to(raw_root), names)
        file_epoch = _epoch_from_path(file)
        for block in discover_chart_blocks(payload, names, seed=file_seed, epoch=file_epoch):
            out["chart_blocks"] += 1
            out["chart_records"] += len(block["charts"])
            key = block["path"].split(".")[-1]
            out["chart_container_keys"][key] += 1
            explicit = 0; point_only = 0
            for chart in block["charts"][:2000]:
                for k in chart:
                    out["record_keys"][k] += 1
                for container in EXPOSURE_CONTAINER_KEYS:
                    if container in chart:
                        out["exposure_container_keys"][container] += 1
                for _, entry, word in _walk_exposure_entries(chart, len(names), "$"):
                    explicit += 1
                def count(obj, depth=0):
                    nonlocal point_only
                    if depth > 6:
                        return
                    if isinstance(obj, dict):
                        if looks_like_point_only(obj):
                            point_only += 1
                        for c in EXPOSURE_CONTAINER_KEYS:
                            if c in obj:
                                count(obj[c], depth+1)
                    elif isinstance(obj, list):
                        for item in obj:
                            count(item, depth+1)
                count(chart)
            out["explicit_word_records"] += explicit
            out["point_only_records"] += point_only
            if len(out["samples"]) < max_samples:
                out["samples"].append({"file": str(file), "kind": "chart_block", "path": block["path"],
                                       "seed": block["seed"], "epoch": block["epoch"],
                                       "charts": len(block["charts"]), "explicit_word_records": explicit,
                                       "point_only_records": point_only,
                                       "sample_chart_keys": sorted(block["charts"][0].keys())[:80]})
    for key in ("top_level_schemas", "chart_container_keys", "exposure_container_keys", "record_keys"):
        out[key] = dict(out[key].most_common(100))
    out["boundary"] = "Probe is schema inventory only. Point-only records and heuristic chart blocks are not exposure evidence."
    return out
