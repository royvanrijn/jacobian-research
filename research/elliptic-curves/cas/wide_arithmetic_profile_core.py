#!/usr/bin/env python3
"""Pure-Python helpers for the wide elliptic-curve arithmetic profile.

This module deliberately contains no Sage imports so the census normalizer,
stratification, deterministic control selection and reporting can be tested in
ordinary Python.  Arithmetic claims are produced only by the Sage worker.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import csv
import ast
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
from typing import Any, Iterable

SCHEMA = "elliptic-curves.wide-arithmetic-profile.v1"

RANK_KEYS_FINAL = (
    "final_rank_lower_bound", "certified_rank_lower_bound", "rank_lower_bound",
    "best_rank_lower_bound", "certified_lower_bound", "final_rank", "best_rank",
    "current_rank", "rank",
)
RANK_KEYS_INITIAL_STAGE = (
    "initial_stage_rank_lower_bound", "initial_search_rank_lower_bound",
    "pre_followup_rank_lower_bound", "pre_follow_up_rank_lower_bound",
    "first_pass_rank_lower_bound", "first_round_rank_lower_bound",
)
GENERIC_RANK_KEYS = (
    "generic_rank", "generic_rank_lower_bound", "starting_generic_rank",
    "seed_rank", "starting_rank", "start_rank",
)
STAGE_KEYS = ("stage", "round", "phase", "pass", "search_stage", "continuation_round")
ID_KEYS = (
    "live_case_id", "case_id", "fibre_id", "fiber_id", "curve_id", "id", "name",
)
T_KEYS = ("t", "parameter", "specialization_parameter", "specialisation_parameter")
AINV_KEYS = ("ainvs", "a_invariants", "a-invariants", "weierstrass_coefficients")
NESTED_CURVE_KEYS = ("curve", "elliptic_curve", "model", "weierstrass", "equation")


def stable_json(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_immutable(path: Path, obj: Any) -> None:
    path = Path(path)
    raw = stable_json(obj).encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != raw:
            raise RuntimeError(f"refusing to overwrite different checkpoint: {path}")
        return
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_bytes(raw)
    tmp.replace(path)


def parse_int(value: Any) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        s = value.strip().replace("≥", "").replace(">=", "").strip()
        m = re.fullmatch(r"[-+]?\d+", s)
        if m:
            return int(s)
    return None


def rational_string(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, (int, float, str)):
        try:
            q = Fraction(str(value))
            return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
        except Exception:
            return str(value)
    if isinstance(value, dict):
        n = value.get("numerator", value.get("num"))
        d = value.get("denominator", value.get("den"))
        if n is not None and d is not None:
            try:
                q = Fraction(int(n), int(d))
                return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
            except Exception:
                pass
    return None


def _as_int_list(value: Any) -> list[int] | None:
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("[") and text.endswith("]"):
            try:
                value = ast.literal_eval(text)
            except Exception:
                return None
    if not isinstance(value, (list, tuple)) or len(value) != 5:
        return None
    out = []
    for x in value:
        v = parse_int(x)
        if v is None:
            return None
        out.append(v)
    return out


def extract_ainvs(obj: Any, depth: int = 0) -> list[int] | None:
    """Extract [a1,a2,a3,a4,a6] from common repository packet shapes.

    This intentionally does not parse free-form equation strings.  Unsupported
    rows remain unbound rather than risking a silently wrong curve.
    """
    if depth > 4:
        return None
    direct = _as_int_list(obj)
    if direct is not None:
        return direct
    if not isinstance(obj, dict):
        return None
    for key in AINV_KEYS:
        if key in obj:
            got = _as_int_list(obj[key])
            if got is not None:
                return got
    if all(k in obj for k in ("a1", "a2", "a3", "a4", "a6")):
        got = [parse_int(obj[k]) for k in ("a1", "a2", "a3", "a4", "a6")]
        if all(v is not None for v in got):
            return [int(v) for v in got]
    # Short Weierstrass packets.
    if "A" in obj and "B" in obj:
        A, B = parse_int(obj["A"]), parse_int(obj["B"])
        if A is not None and B is not None:
            return [0, 0, 0, A, B]
    for key in NESTED_CURVE_KEYS + ("result", "summary", "metadata", "certificate", "state"):
        if key in obj:
            got = extract_ainvs(obj[key], depth + 1)
            if got is not None:
                return got
    return None


def recursive_values(obj: Any, keys: Iterable[str], depth: int = 0) -> list[Any]:
    if depth > 5:
        return []
    out = []
    if isinstance(obj, dict):
        for k in keys:
            if k in obj:
                out.append(obj[k])
        for k, v in obj.items():
            if k in NESTED_CURVE_KEYS or k in ("result", "summary", "metadata", "certificate", "search", "state"):
                out.extend(recursive_values(v, keys, depth + 1))
    return out


def first_int(obj: dict, keys: Iterable[str]) -> int | None:
    for value in recursive_values(obj, keys):
        got = parse_int(value)
        if got is not None:
            return got
    return None


def first_string(obj: dict, keys: Iterable[str]) -> str | None:
    for value in recursive_values(obj, keys):
        if isinstance(value, (str, int)):
            return str(value)
    return None


def first_rational(obj: dict, keys: Iterable[str]) -> str | None:
    for value in recursive_values(obj, keys):
        got = rational_string(value)
        if got is not None:
            return got
    return None


def normalize_row(obj: dict, source: str, ordinal: int) -> dict | None:
    ainvs = extract_ainvs(obj)
    final_rank = first_int(obj, RANK_KEYS_FINAL)
    if ainvs is None or final_rank is None:
        return None
    initial_stage_rank = first_int(obj, RANK_KEYS_INITIAL_STAGE)
    generic_rank = first_int(obj, GENERIC_RANK_KEYS)
    stage = first_string(obj, STAGE_KEYS)
    row_id = first_string(obj, ID_KEYS)
    t = first_rational(obj, T_KEYS)
    if row_id is None:
        row_id = f"anon-{hashlib.sha256((source + ':' + str(ordinal) + ':' + repr(ainvs)).encode()).hexdigest()[:16]}"
    return {
        "id": row_id,
        "t": t,
        "ainvs": ainvs,
        "declared_initial_stage_rank_lower_bound": initial_stage_rank,
        "generic_rank_lower_bound": generic_rank,
        "stage": stage,
        "final_rank_lower_bound": final_rank,
        "source": source,
        "source_ordinal": ordinal,
    }


def walk_candidate_dicts(obj: Any, max_nodes: int = 2_000_000) -> Iterable[dict]:
    stack = [obj]
    seen = 0
    while stack:
        cur = stack.pop()
        seen += 1
        if seen > max_nodes:
            raise RuntimeError("candidate JSON exceeds node safety cap")
        if isinstance(cur, dict):
            yield cur
            stack.extend(cur.values())
        elif isinstance(cur, list):
            stack.extend(cur)


def read_rows_from_file(path: Path) -> list[dict]:
    path = Path(path)
    out: list[dict] = []
    suffix = path.suffix.lower()
    if suffix == ".csv":
        with path.open(newline="") as f:
            for i, row in enumerate(csv.DictReader(f)):
                got = normalize_row(row, str(path), i)
                if got:
                    out.append(got)
        return out
    if suffix in (".jsonl", ".ndjson"):
        with path.open() as f:
            for i, line in enumerate(f):
                if not line.strip():
                    continue
                obj = json.loads(line)
                if isinstance(obj, dict):
                    got = normalize_row(obj, str(path), i)
                    if got:
                        out.append(got)
        return out
    if suffix != ".json":
        return out
    data = json.loads(path.read_text())
    for i, obj in enumerate(walk_candidate_dicts(data)):
        got = normalize_row(obj, str(path), i)
        if got:
            out.append(got)
    return out


def curve_key(row: dict) -> str:
    # Equation is canonical for the census. t is intentionally not required:
    # some final packets omit the original specialization parameter.
    return sha256_bytes(stable_json(row["ainvs"]).encode())[:24]


def aggregate_rows(rows: Iterable[dict]) -> list[dict]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[curve_key(row)].append(row)
    out = []
    for key, items in grouped.items():
        finals = [int(x["final_rank_lower_bound"]) for x in items]
        declared_initials = [x.get("declared_initial_stage_rank_lower_bound") for x in items if x.get("declared_initial_stage_rank_lower_bound") is not None]
        generic_ranks = [x.get("generic_rank_lower_bound") for x in items if x.get("generic_rank_lower_bound") is not None]
        observed_ranks = sorted(set(finals))
        stage_initials = []
        for x in items:
            label = str(x.get("stage") or "").lower().replace("_", "-")
            if label in ("initial", "initial-stage", "first", "first-pass", "round-0", "0"):
                stage_initials.append(int(x["final_rank_lower_bound"]))
        distinct_sources = {str(x.get("source")) for x in items}
        if declared_initials:
            initial = min(map(int, declared_initials)); initial_source = "declared_initial_stage"
        elif stage_initials:
            initial = min(stage_initials); initial_source = "stage_label"
        elif len(observed_ranks) > 1:
            initial = min(observed_ranks); initial_source = "earliest_observed_rank"
        elif len(distinct_sources) > 1:
            # Same bound persisted into at least one later source file.
            initial = observed_ranks[0]; initial_source = "multi_source_same_bound"
        else:
            initial = None; initial_source = None
        chosen = sorted(items, key=lambda r: (str(r.get("id")), str(r.get("source")), int(r.get("source_ordinal", 0))))[0]
        ids = sorted({str(x.get("id")) for x in items if x.get("id")})
        ts = sorted({str(x.get("t")) for x in items if x.get("t") is not None})
        out.append({
            "curve_key": key,
            "id": ids[0] if ids else chosen["id"],
            "aliases": ids,
            "t": ts[0] if len(ts) == 1 else None,
            "t_candidates": ts,
            "ainvs": chosen["ainvs"],
            "generic_rank_lower_bound": min(map(int, generic_ranks)) if generic_ranks else None,
            "initial_rank_lower_bound": initial,
            "initial_rank_source": initial_source,
            "final_rank_lower_bound": max(finals),
            "observed_rank_lower_bounds": observed_ranks,
            "record_count": len(items),
            "improved_since_initial": (max(finals) > initial) if initial is not None else None,
            "sources": sorted({str(x["source"]) for x in items}),
        })
    return sorted(out, key=lambda r: (r["curve_key"], r["id"]))


def discover_source_files(source: Path) -> list[Path]:
    source = Path(source)
    if source.is_file():
        return [source]
    if not source.is_dir():
        raise FileNotFoundError(source)
    files = []
    for p in source.rglob("*"):
        if p.is_file() and p.suffix.lower() in (".json", ".jsonl", ".ndjson", ".csv"):
            # Skip huge point clouds / per-chart artifacts by filename when obvious.
            low = p.name.lower()
            if any(token in low for token in ("cloud-", "points-", "charts-", "transcript")):
                continue
            files.append(p)
    return sorted(files)


def brumer_kramer_local_term(delta_sign: int, phi_m_count: int, additive_split_counts: Iterable[int]) -> dict:
    if delta_sign == 0:
        raise ValueError("discriminant sign must be nonzero")
    u = 1 if delta_sign < 0 else 2
    additive = [int(n) for n in additive_split_counts]
    if any(n < 1 or n > 3 for n in additive):
        raise ValueError("cubic splitting count must be in 1..3")
    n_term = int(phi_m_count) + sum(n - 1 for n in additive)
    return {"u": u, "n": n_term, "local_term": u + n_term}


def rank_bucket(rank: int) -> str:
    if rank >= 25:
        return "25+"
    if rank == 24:
        return "24"
    if rank == 23:
        return "23"
    if rank >= 21:
        return "21-22"
    if rank >= 19:
        return "19-20"
    return "17-18"


def median(values: list[float]) -> float | None:
    vals = [float(x) for x in values if x is not None and math.isfinite(float(x))]
    return statistics.median(vals) if vals else None


def _average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        r = (i + 1 + j) / 2.0  # 1-based average rank
        for k in range(i, j):
            ranks[order[k]] = r
        i = j
    return ranks


def spearman(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) != len(ys) or len(xs) < 3:
        return None
    rx, ry = _average_ranks(list(map(float, xs))), _average_ranks(list(map(float, ys)))
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def deterministic_controls(rows: list[dict], high_threshold: int = 23, per_bucket: int = 30, salt: str = "wide-arithmetic-profile-v1") -> dict:
    high = [r for r in rows if int(r["final_rank_lower_bound"]) >= high_threshold]
    buckets: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if int(row["final_rank_lower_bound"]) >= high_threshold:
            continue
        buckets[rank_bucket(int(row["final_rank_lower_bound"]))].append(row)
    controls = []
    for bucket in sorted(buckets):
        ordered = sorted(
            buckets[bucket],
            key=lambda r: (hashlib.sha256((salt + "/" + r["curve_key"]).encode()).hexdigest(), r["curve_key"]),
        )
        controls.extend(ordered[:per_bucket])
    return {
        "high": [r["curve_key"] for r in sorted(high, key=lambda x: x["curve_key"])],
        "controls": [r["curve_key"] for r in controls],
        "control_buckets": {b: min(per_bucket, len(v)) for b, v in sorted(buckets.items())},
        "salt": salt,
    }


def flatten_profile(row: dict, base: dict | None, local: dict | None, class_probe: dict | None) -> dict:
    out = {
        "curve_key": row["curve_key"],
        "id": row["id"],
        "t": row.get("t"),
        "initial_rank_lower_bound": row.get("initial_rank_lower_bound"),
        "final_rank_lower_bound": row["final_rank_lower_bound"],
        "improved_since_initial": row.get("improved_since_initial"),
        "rank_bucket": rank_bucket(int(row["final_rank_lower_bound"])),
    }
    if base and base.get("status") == "PASS":
        for k in (
            "delta_sign", "log2_abs_minimal_discriminant", "log2_abs_polynomial_discriminant",
            "rational_2torsion_rank", "cubic_irreducible",
        ):
            out[k] = base.get(k)
    if local and local.get("status") == "PASS":
        for k in (
            "field_signature", "field_discriminant", "log2_abs_field_discriminant",
            "field_ramified_prime_count", "phi_m_count", "phi_a_count", "bk_u", "bk_n",
            "bk_local_term", "root_number", "conductor", "log2_conductor",
        ):
            out[k] = local.get(k)
        if local.get("bk_local_term") is not None:
            out["forced_class_2rank_lower_from_known_rank"] = max(
                0, int(row["final_rank_lower_bound"]) - int(local["bk_local_term"])
            )
            out["g_upper_needed_to_close_known_rank"] = int(row["final_rank_lower_bound"]) - int(local["bk_local_term"])
    if class_probe and class_probe.get("status") == "PASS_PROVISIONAL_GRH":
        out["provisional_class_2rank"] = class_probe.get("class_2rank")
        if out.get("bk_local_term") is not None:
            out["provisional_grh_selmer_upper"] = int(class_probe["class_2rank"]) + int(out["bk_local_term"])
            out["provisional_grh_closes_known_rank"] = out["provisional_grh_selmer_upper"] <= int(row["final_rank_lower_bound"])
    return out


def summarize_profiles(flat: list[dict]) -> dict:
    numeric_metrics = [
        "bk_local_term", "bk_n", "phi_m_count", "phi_a_count",
        "log2_abs_minimal_discriminant", "log2_abs_field_discriminant",
        "field_ramified_prime_count", "log2_conductor", "forced_class_2rank_lower_from_known_rank",
    ]
    groups: dict[str, list[dict]] = defaultdict(list)
    for row in flat:
        groups[row["rank_bucket"]].append(row)
    by_bucket = {}
    for bucket, rows in sorted(groups.items()):
        entry = {
            "count": len(rows),
            "improved_known_count": sum(1 for r in rows if r.get("improved_since_initial") is True),
            "root_number_minus_one_count": sum(1 for r in rows if r.get("root_number") == -1),
        }
        for metric in numeric_metrics:
            vals = [r[metric] for r in rows if r.get(metric) is not None]
            entry[f"median_{metric}"] = median(vals)
            entry[f"n_{metric}"] = len(vals)
        by_bucket[bucket] = entry
    correlations = {}
    for metric in numeric_metrics:
        pairs = [(float(r[metric]), float(r["final_rank_lower_bound"])) for r in flat if r.get(metric) is not None]
        correlations[f"spearman_{metric}_vs_final_lower_bound"] = spearman(
            [p[0] for p in pairs], [p[1] for p in pairs]
        )
        pairs2 = [
            (float(r[metric]), float(int(r["final_rank_lower_bound"]) - int(r["initial_rank_lower_bound"])))
            for r in flat if r.get(metric) is not None and r.get("initial_rank_lower_bound") is not None
        ]
        correlations[f"spearman_{metric}_vs_followup_gain"] = spearman(
            [p[0] for p in pairs2], [p[1] for p in pairs2]
        )
    return {
        "schema": SCHEMA,
        "population_count": len(flat),
        "rank_lower_bound_counts": dict(sorted((str(k), v) for k, v in _counts(int(r["final_rank_lower_bound"]) for r in flat).items())),
        "bucket_summary": by_bucket,
        "correlations": correlations,
        "interpretation_boundary": (
            "Rank values are certified lower bounds, not exact ranks. Follow-up gain is search-outcome data and is adaptively sampled. "
            "Correlations are descriptive only; no p-values or causal interpretation are claimed."
        ),
    }


def _counts(values: Iterable[int]) -> dict[int, int]:
    d: dict[int, int] = defaultdict(int)
    for v in values:
        d[int(v)] += 1
    return dict(d)


def markdown_summary(summary: dict, flat: list[dict]) -> str:
    lines = [
        "# Wide-search equation-only arithmetic profile",
        "",
        f"Population: **{summary['population_count']}** frozen fibres.",
        "",
        "> Rank values below are certified lower bounds, not exact ranks. Follow-up improvements are search outcomes from an adaptive campaign. No p-values or causal claims are made.",
        "",
        "## Final certified lower-bound distribution",
        "",
        "| Lower bound | Fibres |",
        "|---:|---:|",
    ]
    for rank, count in sorted(((int(k), v) for k, v in summary["rank_lower_bound_counts"].items()), reverse=True):
        lines.append(f"| ≥{rank} | {count} |")
    lines += ["", "## Arithmetic by lower-bound stratum", ""]
    headers = ["Stratum", "n", "median BK local", "median log2|D_K|", "median ramified primes", "median forced g lower"]
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for bucket, e in summary["bucket_summary"].items():
        def fmt(key):
            v = e.get(key)
            return "—" if v is None else (f"{v:.3f}" if isinstance(v, float) else str(v))
        lines.append(
            f"| {bucket} | {e['count']} | {fmt('median_bk_local_term')} | {fmt('median_log2_abs_field_discriminant')} | "
            f"{fmt('median_field_ramified_prime_count')} | {fmt('median_forced_class_2rank_lower_from_known_rank')} |"
        )
    lines += [
        "",
        "## Descriptive correlations",
        "",
        "These are diagnostics for prioritising theory work. They are confounded by the search design and lower-bound censoring.",
        "",
        "| Metric | Spearman vs final lower bound | Spearman vs follow-up gain |",
        "|---|---:|---:|",
    ]
    for key in sorted(k for k in summary["correlations"] if k.endswith("_vs_final_lower_bound")):
        metric = key[len("spearman_"):-len("_vs_final_lower_bound")]
        a = summary["correlations"].get(key)
        b = summary["correlations"].get(f"spearman_{metric}_vs_followup_gain")
        fa = "—" if a is None else f"{a:.4f}"
        fb = "—" if b is None else f"{b:.4f}"
        lines.append(f"| `{metric}` | {fa} | {fb} |")
    provisional = [r for r in flat if r.get("provisional_class_2rank") is not None]
    if provisional:
        lines += [
            "", "## Bounded provisional class-group probes", "",
            "These values are **conditional/provisional** (`proof=False` / GRH-dependent class-group computation). They are not unconditional rank bounds.", "",
            "| ID | rank LB | BK local | provisional g | provisional Selmer upper | closes LB? |",
            "|---|---:|---:|---:|---:|---|",
        ]
        for r in sorted(provisional, key=lambda x: (-int(x["final_rank_lower_bound"]), x["id"])):
            lines.append(
                f"| {r['id']} | {r['final_rank_lower_bound']} | {r.get('bk_local_term','—')} | {r['provisional_class_2rank']} | "
                f"{r.get('provisional_grh_selmer_upper','—')} | {r.get('provisional_grh_closes_known_rank','—')} |"
            )
    lines += [
        "", "## Boundaries", "",
        "- `bk_local_term = u(E)+n(E)` is equation-derived after exact local reduction and cubic-field splitting.",
        "- `forced_class_2rank_lower_from_known_rank = max(0, rank_lower_bound - bk_local_term)` uses the known lower bound and is therefore not an independent predictor.",
        "- Timed-out factoring, maximal-order or class-group work remains `UNKNOWN`; it is never recorded as zero.",
        "- Curves with rational 2-torsion are reported but the Brumer–Kramer bound used here is marked not applicable.",
        "",
    ]
    return "\n".join(lines)
