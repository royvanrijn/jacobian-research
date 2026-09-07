#!/usr/bin/env python3
"""Attach exact known-bisection outcomes to the frozen R17 training cohort.

This is an exhaustive evaluation of an existing certified atlas, not a new
cover or rational-point search.  Nonzero split covers are specialized exactly
and their gain over the seventeen generic sections is certified in finite
reduction quotients.  Failure to find a split cover is a zero for this
mechanism only; it is not a rank label or a Mordell--Weil dependence claim.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from math import log
from multiprocessing import get_context
import os
from pathlib import Path
from resource import RUSAGE_SELF, getrusage
import sys
from time import perf_counter, process_time
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
ECSEARCH = ELLIPTIC_ROOT / "ecsearch"
CAS = ELLIPTIC_ROOT / "cas"
SCRIPTS = ELLIPTIC_ROOT / "scripts"
sys.path[:0] = [str(ELLIPTIC_ROOT), str(ECSEARCH), str(CAS), str(SCRIPTS)]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    homogeneous_value,
    load_q12o5867_data,
)
from r17_bisection_ranker import semantic_label_sha256  # noqa: E402
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    matrix_rank_and_pivots_mod_prime,
    verify_finite_quotient_certificate,
)
from evaluate_elkies_2026_bisections_at_controls import (  # noqa: E402
    construct_hits,
    fraction_text,
    point_record,
    rows_sha256,
    short_add,
    short_linear_combination,
    signature_rows,
    solve_columns_mod_two,
)
from r17_training_data import (  # noqa: E402
    EMBARGOED_PARAMETERS,
    normalize_quadratic,
    normalized_parameter,
    split_quadratic_indices,
)


MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
DEFAULT_INPUT = ROOT / "artifacts/local/elliptic-curves/r17-training-selected.jsonl"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elliptic-curves/r17-training-bisection-labels.jsonl"
DEFAULT_SUMMARY = ROOT / "artifacts/local/elliptic-curves/r17-training-bisection-labels-summary.json"
SCHEMA = "elliptic-curves.r17-training-bisection-label.v1"


Quadratic = tuple[tuple[int, int, int], int]
_WORKER_QUADRATICS: tuple[Quadratic, ...] = ()


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def _scan_parameter(task: tuple[int, tuple[int, int]]) -> tuple[int, list[int], list[int], float]:
    row_index, parameter = task
    started = process_time()
    split, ramified = split_quadratic_indices(parameter, _WORKER_QUADRATICS)
    return row_index, split, ramified, process_time() - started


def scan_all(
    parameters: Sequence[tuple[int, int]], quadratics: tuple[Quadratic, ...], workers: int
) -> list[tuple[list[int], list[int], float]]:
    global _WORKER_QUADRATICS
    _WORKER_QUADRATICS = quadratics
    tasks = list(enumerate(parameters))
    results: list[tuple[list[int], list[int], float] | None] = [None] * len(tasks)
    if workers == 1:
        iterator = map(_scan_parameter, tasks)
        pool = None
    else:
        # The atlas is read-only and inherited copy-on-write on the supported
        # Linux research host.  No 39,120-row pickle is sent with each task.
        pool = get_context("fork").Pool(processes=workers)
        iterator = pool.imap(_scan_parameter, tasks, chunksize=4)
    try:
        for completed, (index, split, ramified, cpu_seconds) in enumerate(iterator, 1):
            results[index] = (split, ramified, cpu_seconds)
            if completed % 250 == 0 or completed == len(tasks):
                print(f"R17BISECTIONLABEL|phase=square-census|rows={completed}/{len(tasks)}", flush=True)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    if any(result is None for result in results):
        raise AssertionError("the parallel square census omitted a row")
    return [result for result in results if result is not None]


def exact_j_invariant(data, parameter: tuple[int, int]) -> Fraction:
    a, b = parameter
    coefficient_a = homogeneous_value(data.a_coefficients, a, b, 8)
    coefficient_b = homogeneous_value(data.b_coefficients, a, b, 12)
    denominator = 4 * coefficient_a**3 + 27 * coefficient_b**2
    if denominator == 0:
        raise ValueError(f"singular R17 fibre at {a}/{b}")
    return Fraction(6912 * coefficient_a**3, denominator)


def x_naive_height(point: tuple[Fraction, Fraction]) -> int:
    x_coordinate = Fraction(point[0])
    return max(abs(x_coordinate.numerator), x_coordinate.denominator)


def compact_hit_record(hit: dict[str, Any], *, escapes_generic_span: bool) -> dict[str, Any]:
    source = hit["source_record"]
    return {
        "label": source["label"],
        "lattice_orbit_mask": int(source["lattice_orbit_mask"]),
        "q_value": fraction_text(hit["q_value"]),
        "canonical_positive_square_root": fraction_text(hit["square_root"]),
        "positive_source_point": point_record(hit["positive_source_point"]),
        "negative_source_point": point_record(hit["negative_source_point"]),
        "trace_source_point": point_record(hit["trace_source_point"]),
        "positive_x_naive_height": str(x_naive_height(hit["positive_source_point"])),
        "negative_x_naive_height": str(x_naive_height(hit["negative_source_point"])),
        "escapes_generic_span_in_stored_finite_quotients": escapes_generic_span,
        "exact_verification": {
            "both_branches_on_source_fibre": True,
            "sum_of_branches_equals_stored_trace": True,
            "stored_trace_equals_published_basis_word": True,
        },
    }


def certify_split_hits(
    *,
    data,
    atlas_records: Sequence[dict[str, Any]],
    hit_indices: Sequence[int],
    parameter: tuple[int, int],
    prime_bound: int,
) -> dict[str, Any]:
    selected_records = [atlas_records[index] for index in hit_indices]
    hits = construct_hits(selected_records, *parameter)
    if [hit["source_record"]["label"] for hit in hits] != [
        record["label"] for record in selected_records
    ]:
        raise AssertionError("integer and rational split-cover tests disagree")

    specialization = evaluate_projective_specialization(data, *parameter)
    coefficient_a = specialization.coefficient_a
    for hit in hits:
        for name in ("positive_source_point", "negative_source_point", "trace_source_point"):
            if not is_on_weierstrass_curve(specialization.model, hit[name]):
                raise ArithmeticError(f"{hit['source_record']['label']} missed the source fibre")
        if short_add(
            hit["positive_source_point"], hit["negative_source_point"], coefficient_a
        ) != hit["trace_source_point"]:
            raise ArithmeticError("the two bisection branches missed their stored trace")
        trace_word = [int(value) for value in hit["source_record"]["published_basis_w"]]
        if short_linear_combination(specialization.points, trace_word, coefficient_a) != hit[
            "trace_source_point"
        ]:
            raise ArithmeticError("a bisection trace missed its generic-basis word")

    all_points = specialization.points + tuple(hit["positive_source_point"] for hit in hits)
    certificate = build_finite_quotient_certificate(
        specialization.model,
        all_points,
        relation_prime=2,
        prime_bound=prime_bound,
    )
    verify_finite_quotient_certificate(specialization.model, all_points, certificate)
    rows = signature_rows(certificate)
    baseline_rank, _baseline_pivots = matrix_rank_and_pivots_mod_prime(
        [row[:17] for row in rows], 17, 2
    )
    joint_rank = int(certificate["combined_rank_over_relation_field"])
    torsion_witness = certificate["torsion_witness"]
    if baseline_rank != 17 or torsion_witness is None:
        return {
            "status": "CENSORED_FINITE_QUOTIENT_BASELINE_NOT_CERTIFIED",
            "blocker": (
                f"generic finite-quotient rank={baseline_rank}; "
                f"torsion_witness_present={torsion_witness is not None}"
            ),
            "baseline_rank": baseline_rank,
            "joint_rank": joint_rank,
            "gain": None,
            "minimum_height": None,
            "hit_records": [compact_hit_record(hit, escapes_generic_span=False) for hit in hits],
            "certificate": certificate,
        }

    gain = joint_rank - 17
    if gain < 0:
        raise AssertionError("the joint finite quotient lost a generic direction")
    escaping_indices = []
    for hit_index in range(len(hits)):
        if solve_columns_mod_two(rows, list(range(17)), 17 + hit_index) is None:
            escaping_indices.append(hit_index)
    if bool(gain) != bool(escaping_indices):
        raise AssertionError("joint gain and individual escape tests disagree")
    minimum_height = None
    if escaping_indices:
        minimum_height = min(
            x_naive_height(hits[index][branch])
            for index in escaping_indices
            for branch in ("positive_source_point", "negative_source_point")
        )
    pivot_hit_indices = [
        index - 17
        for index in certificate["pivot_columns_zero_based"]
        if int(index) >= 17
    ]
    if len(pivot_hit_indices) != gain:
        raise AssertionError("finite-quotient pivots disagree with the certified gain")
    hit_records = [
        compact_hit_record(hit, escapes_generic_span=index in escaping_indices)
        for index, hit in enumerate(hits)
    ]
    return {
        "status": "PASS_EXACT_KNOWN_BISECTION_GAIN_LOWER_BOUND",
        "blocker": None,
        "baseline_rank": baseline_rank,
        "joint_rank": joint_rank,
        "gain": gain,
        "minimum_height": minimum_height,
        "independent_hit_labels": [hits[index]["source_record"]["label"] for index in pivot_hit_indices],
        "hit_records": hit_records,
        "certificate": certificate,
        "stacked_rows_sha256": rows_sha256(rows),
    }


def make_label(
    *,
    row: dict[str, Any],
    data,
    atlas_records: Sequence[dict[str, Any]],
    scan: tuple[list[int], list[int], float],
    j_invariant: Fraction,
    j_group_size: int,
    j_group_splits: Sequence[str],
    prime_bound: int,
) -> dict[str, Any]:
    hit_indices, ramified_indices, scan_cpu_seconds = scan
    started = process_time()
    certification = None
    if hit_indices:
        certification = certify_split_hits(
            data=data,
            atlas_records=atlas_records,
            hit_indices=hit_indices,
            parameter=tuple(row["projective_pair"]),
            prime_bound=prime_bound,
        )
        gain = certification["gain"]
        minimum_height = certification["minimum_height"]
        status = certification["status"]
        blocker = certification["blocker"]
        hit_records = certification["hit_records"]
    else:
        gain = 0
        minimum_height = None
        status = "PASS_EXACT_KNOWN_BISECTION_CENSUS_NO_NONZERO_SPLIT"
        blocker = None
        hit_records = []
    cpu_seconds = scan_cpu_seconds + (process_time() - started)
    cost_per_direction = None if not gain else cpu_seconds / gain
    j_text = fraction_text(j_invariant)
    return {
        "schema": SCHEMA,
        "parameter": row["parameter"],
        "projective_pair": row["projective_pair"],
        "split": row["split"],
        "selection_lanes": row["selection_lanes"],
        "label_status": status,
        "arithmetic_grouping": {
            "exact_j_invariant": j_text,
            "exact_j_sha256": sha256(j_text.encode()).hexdigest(),
            "same_j_selected_row_count": j_group_size,
            "same_j_splits": list(j_group_splits),
            "global_minimal_q_isomorphism_key": None,
            "boundary": (
                "Exact j groups all rational twists conservatively. Global-minimal-model "
                "keys remain required before model fitting if a same-j group occurs."
            ),
        },
        "outcomes": {
            "new_independent_directions_lower_bound": gain,
            "finite_quotient_gain_lower_bound": gain,
            "minimum_exceptional_naive_height": (
                None if minimum_height is None else str(minimum_height)
            ),
            "minimum_exceptional_log_naive_height": (
                None if minimum_height is None else log(minimum_height)
            ),
            "cpu_seconds": cpu_seconds,
            "cpu_seconds_per_new_direction": cost_per_direction,
            "peak_rss_bytes": None,
            "completed_search_bound": {
                "mechanism": "complete_preexisting_certified_r17_bisection_atlas",
                "nonzero_split_square_tests": len(atlas_records),
                "finite_reduction_prime_bound": prime_bound if hit_indices else None,
            },
            "residual_two_selmer_dimension": None,
            "unrestricted_point_search_authorized": False,
            "censoring_or_blocker": blocker,
        },
        "bisection_census": {
            "atlas_size": len(atlas_records),
            "nonzero_split_bisection_count": len(hit_indices),
            "ramified_bisection_count": len(ramified_indices),
            "ramified_labels": [atlas_records[index]["label"] for index in ramified_indices],
            "hits": hit_records,
        },
        "finite_quotient_certificate": (
            None if certification is None else certification["certificate"]
        ),
        "finite_quotient_audit": (
            None
            if certification is None
            else {
                key: value
                for key, value in certification.items()
                if key not in {"certificate", "hit_records", "minimum_height", "gain"}
            }
        ),
        "proof_boundary": (
            "The outcome is an exact lower bound for points exposed by the already-certified "
            "bisection atlas. Zero means no gain from its nonzero split specializations, not "
            "rank 17, dependence, or failure of any other point-search mechanism."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--prime-bound", type=int, default=199)
    parser.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--limit", type=int, help="diagnostic prefix only; never a production corpus")
    args = parser.parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be positive")
    if args.prime_bound < 3:
        raise SystemExit("--prime-bound must be at least 3")

    started = perf_counter()
    sys.set_int_max_str_digits(0)
    rows = [json.loads(line) for line in args.input.read_text().splitlines() if line.strip()]
    if args.limit is not None:
        if args.limit < 1:
            raise SystemExit("--limit must be positive")
        rows = rows[: args.limit]
    parameters = [normalized_parameter(*row["projective_pair"]) for row in rows]
    if EMBARGOED_PARAMETERS.intersection(parameters):
        raise AssertionError("an embargoed positive control entered the label cohort")
    if len(parameters) != len(set(parameters)):
        raise AssertionError("the selected cohort contains duplicate parameters")

    data = load_q12o5867_data(MODEL, SECTIONS)
    atlas = json.loads(BISECTIONS.read_text())
    atlas_records = atlas["bisections"]
    if len(atlas_records) != 39_120:
        raise AssertionError("the certified bisection atlas changed size")
    quadratics = tuple(
        normalize_quadratic(record["residual_chord"]["q_coefficients"])
        for record in atlas_records
    )
    scans = scan_all(parameters, quadratics, args.workers)

    j_values = [exact_j_invariant(data, parameter) for parameter in parameters]
    j_groups: dict[Fraction, list[int]] = defaultdict(list)
    for index, value in enumerate(j_values):
        j_groups[value].append(index)
    cross_split_groups = []
    for value, indices in j_groups.items():
        splits = sorted({rows[index]["split"] for index in indices})
        if len(splits) > 1:
            cross_split_groups.append(
                {
                    "exact_j_sha256": sha256(fraction_text(value).encode()).hexdigest(),
                    "parameters": [rows[index]["parameter"] for index in indices],
                    "splits": splits,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    counters = Counter()
    lane_counters: dict[str, Counter] = defaultdict(Counter)
    with args.output.open("w") as handle:
        for index, row in enumerate(rows):
            group_indices = j_groups[j_values[index]]
            label = make_label(
                row=row,
                data=data,
                atlas_records=atlas_records,
                scan=scans[index],
                j_invariant=j_values[index],
                j_group_size=len(group_indices),
                j_group_splits=sorted({rows[group]["split"] for group in group_indices}),
                prime_bound=args.prime_bound,
            )
            handle.write(json.dumps(label, sort_keys=True, separators=(",", ":")) + "\n")
            outcome = label["outcomes"]
            gain = outcome["finite_quotient_gain_lower_bound"]
            split_count = label["bisection_census"]["nonzero_split_bisection_count"]
            counters["rows"] += 1
            counters["nonzero_split_covers"] += split_count
            counters["ramified_covers"] += label["bisection_census"]["ramified_bisection_count"]
            counters["rows_with_split"] += bool(split_count)
            counters["censored_rows"] += gain is None
            counters["rows_with_certified_gain"] += gain is not None and gain > 0
            if gain is not None:
                counters["certified_gain_sum"] += gain
            for lane in row["selection_lanes"]:
                lane_counters[lane]["rows"] += 1
                lane_counters[lane]["rows_with_split"] += bool(split_count)
                lane_counters[lane]["rows_with_certified_gain"] += gain is not None and gain > 0
                if gain is not None:
                    lane_counters[lane]["certified_gain_sum"] += gain
            if (index + 1) % 100 == 0 or index + 1 == len(rows):
                print(f"R17BISECTIONLABEL|phase=exact-certification|rows={index + 1}/{len(rows)}", flush=True)

    elapsed = perf_counter() - started
    summary = {
        "schema": "elliptic-curves.r17-training-bisection-label-summary.v1",
        "status": (
            "EXPERIMENTAL_PREFIX_EXACT_KNOWN_BISECTION_LABELS"
            if args.limit is not None
            else "EXACT_KNOWN_BISECTION_LABELS_NO_FULL_RANK_CLAIM"
        ),
        "row_count": len(rows),
        "atlas_size": len(atlas_records),
        "total_exact_square_tests": len(rows) * len(atlas_records),
        "outcomes": dict(counters),
        "lane_outcomes": {lane: dict(values) for lane, values in sorted(lane_counters.items())},
        "arithmetic_duplicate_audit": {
            "exact_j_group_count": len(j_groups),
            "same_j_group_count": sum(len(indices) > 1 for indices in j_groups.values()),
            "cross_split_same_j_group_count": len(cross_split_groups),
            "cross_split_same_j_groups": cross_split_groups,
            "model_fitting_gate": (
                "PASS_NO_SAME_J_DUPLICATES"
                if len(j_groups) == len(rows)
                else "REQUIRES_GLOBAL_MINIMAL_KEY_GROUPING_BEFORE_FIT"
            ),
        },
        "embargo": {
            "parameters": [f"{a}/{b}" for a, b in sorted(EMBARGOED_PARAMETERS)],
            "all_absent": True,
        },
        "finite_quotient": {
            "relation_prime": 2,
            "reduction_prime_bound": args.prime_bound,
        },
        "generation": {
            "command": " ".join(
                (
                    "python3 elliptic-curves/scripts/label_r17_training_bisections.py",
                    f"--input {relative(args.input)}",
                    f"--output {relative(args.output)}",
                    f"--summary {relative(args.summary)}",
                    f"--workers {args.workers}",
                    f"--prime-bound {args.prime_bound}",
                )
            ),
            "workers": args.workers,
            "runtime_seconds": elapsed,
            "peak_parent_rss_bytes": int(getrusage(RUSAGE_SELF).ru_maxrss) * 1024,
            "script_sha256": file_sha256(Path(__file__)),
            "inputs": {
                relative(path): file_sha256(path)
                for path in (args.input, MODEL, SECTIONS, BISECTIONS)
            },
            "output": str(args.output.resolve()),
            "output_sha256": file_sha256(args.output),
            "semantic_label_sha256_excluding_timings": semantic_label_sha256(
                json.loads(line) for line in args.output.read_text().splitlines() if line.strip()
            ),
        },
        "proof_boundary": [
            "The complete preexisting bisection atlas was evaluated; no new cover or point search was run.",
            "Certified gains are finite-reduction lower bounds over the generic seventeen.",
            "Zero is mechanism-specific and is not a rank, dependence, or Selmer claim.",
            "Ramified cover values are inventoried but are not counted as nonzero two-branch splits.",
            "No residual two-Selmer dimension or unrestricted deep-search outcome was manufactured.",
        ],
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
