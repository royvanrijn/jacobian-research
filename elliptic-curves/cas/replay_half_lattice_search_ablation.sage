#!/usr/bin/env sage -python
"""Fixture-blind equal-budget ablation of half-lattice chart selection.

The development phase contains only the four published-R17 specialization
parameters used for the +8 through +11 controls.  The holdout phase reads the
frozen generic-only inputs for alternate-Q80 curve 12 and ICARM curves 356 and
385.  Neither phase imports exceptional-point fixtures.

Every selected mask uses the same specialized shortest representative, PARI
minimization/reduction, and one ``hyperellratpoints`` search.  The arms are the
generic deepest 43, specialized deepest 43, their union, five deterministic
hash-random 43-sets, a generic shallowest-depth 43-set, and a generic
median-depth 43-set.  Shared covers are run once; their measured CPU and wall
cost is charged to every arm containing them.
"""

from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import resource
import sys
import time
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
RANK29_INPUT = ELLIPTIC / "data/half_lattice_rank29_control_inputs_v1.json"
ENGINE = CAS / "half_lattice_fake_descent_replay.sage"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
DEVELOPMENT_OUTPUT = ART / "half_lattice_search_ablation_r17_development_blind_v1.json"
HOLDOUT_OUTPUT = ART / "half_lattice_search_ablation_rank29_holdout_blind_v1.json"

DIMENSION = 17
ARM_SIZE = 43
RANDOM_ARM_COUNT = 5
RANDOM_DOMAIN_SEED = "half-lattice-search-ablation-v1/random-domain"
SHALLOW_TIE_SEED = "half-lattice-search-ablation-v1/shallow"
MEDIAN_TIE_SEED = "half-lattice-search-ablation-v1/median"
SPECIALIZED_SCALE = 1_000_000
SPECIALIZED_CHECK_SCALE = 100_000
DEVELOPMENT_FROZEN_SOURCE_SHA256 = "22b6530c9edd1891ab63604f06409e21d8de2f262f4c0feb5c80bd14922f9b6e"
GENERIC_RANKING_CACHE: dict[tuple[tuple[int, ...], ...], tuple[list[tuple[int, int]], dict[str, Any]]] = {}

# Parameters are development labels only.  No rank, point, or quotient data is
# present in this executable.
R17_DEVELOPMENT = (
    ("r17-development-a", -2, 377),
    ("r17-development-b", -308, 251),
    ("r17-development-c", 2456, 135),
    ("r17-development-d", -9529, 5471),
)

sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from elliptic_candidate_record import source_point_to_target  # noqa: E402
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from search_nagao_u135_alternate_covers import relation_proposals  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def cpu_clock() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def bit_height(value: Fraction) -> int:
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def point_key(point: tuple[Fraction, Fraction]):
    return (bit_height(point[0]), bit_height(point[1]), point)


def hash_order(seed: str, mask: int) -> bytes:
    return sha256(seed.encode("ascii") + b"\0" + int(mask).to_bytes(3, "big")).digest()


def fixed_random_arms() -> dict[str, tuple[int, ...]]:
    masks = list(range(1, 1 << DIMENSION))
    masks.sort(key=lambda mask: (hash_order(RANDOM_DOMAIN_SEED, mask), mask))
    return {
        f"random43-{index + 1}": tuple(
            sorted(masks[index * ARM_SIZE : (index + 1) * ARM_SIZE])
        )
        for index in range(RANDOM_ARM_COUNT)
    }


def choose_hash_sample(masks: Iterable[int], seed: str) -> tuple[int, ...]:
    ordered = sorted(masks, key=lambda mask: (hash_order(seed, mask), mask))
    if len(ordered) < ARM_SIZE:
        raise ArithmeticError(f"only {len(ordered)} masks available for {seed}")
    return tuple(sorted(ordered[:ARM_SIZE]))


def r17_cases() -> list[dict[str, Any]]:
    family = load_q12o5867_data(MODEL, SECTIONS)
    cases = []
    for label, numerator, denominator in R17_DEVELOPMENT:
        specialization = evaluate_projective_specialization(family, numerator, denominator)
        minimal_model, minimal_change, unused = global_minimal_model_with_change(
            specialization.model
        )
        minimal_points = tuple(
            source_point_to_target(point, minimal_change) for point in specialization.points
        )
        short_model, short_change = short_certificate_model(minimal_model)
        generic_points = tuple(
            source_point_to_target(point, short_change) for point in minimal_points
        )
        cases.append(
            {
                "label": label,
                "parameter": f"{numerator}/{denominator}",
                "projective_parameter": [numerator, denominator],
                "short_model": tuple(Fraction(value) for value in short_model),
                "generic_points": generic_points,
                "generic_gram": tuple(tuple(row) for row in engine.GENERIC_GRAM),
                "lineage": "published-r17",
            }
        )
    return cases


def rank29_cases() -> list[dict[str, Any]]:
    frozen = json.loads(RANK29_INPUT.read_text())
    if frozen.get("status") != "FROZEN_GENERIC_SUBGROUP_ONLY_NO_EXCEPTIONAL_COORDINATES":
        raise ValueError("rank29 ablation input is not frozen generic-only data")
    if frozen["boundary"]["output_contains_exceptional_point_coordinates"] is not False:
        raise ValueError("rank29 ablation input contains exceptional coordinates")
    return [
        {
            "label": row["label"],
            "parameter": None,
            "projective_parameter": None,
            "short_model": tuple(Fraction(value) for value in row["short_model"]),
            "generic_points": tuple(
                (Fraction(point[0]), Fraction(point[1])) for point in row["generic_points"]
            ),
            "generic_gram": tuple(
                tuple(int(value) for value in gram_row)
                for gram_row in row["generic_height_gram"]
            ),
            "lineage": row["lineage"],
        }
        for row in frozen["cases"]
    ]


def generic_ranking(gram: Sequence[Sequence[int]]) -> tuple[list[tuple[int, int]], dict[str, Any]]:
    key = tuple(tuple(int(value) for value in row) for row in gram)
    cached = GENERIC_RANKING_CACHE.get(key)
    if cached is not None:
        rows, original_timing = cached
        return rows, {
            **original_timing,
            "reused_from_process_cache": True,
            "incremental_cpu_seconds": 0.0,
            "incremental_wall_seconds": 0.0,
        }
    started_wall = time.monotonic()
    started_cpu = cpu_clock()
    oracle = engine.CosetOracle(key)
    rows = []
    histogram: Counter[int] = Counter()
    maximum_error = 0.0
    for mask in range(1 << DIMENSION):
        norm, unused_representative, error = oracle.solve(mask)
        rows.append((norm, mask))
        histogram[norm] += 1
        maximum_error = max(maximum_error, error)
    timing = {
        "cpu_seconds": cpu_clock() - started_cpu,
        "wall_seconds": time.monotonic() - started_wall,
        "reused_from_process_cache": False,
        "incremental_cpu_seconds": cpu_clock() - started_cpu,
        "incremental_wall_seconds": time.monotonic() - started_wall,
        "minimum_norm_histogram": {str(key): value for key, value in sorted(histogram.items())},
        "maximum_enumeration_error": maximum_error,
    }
    GENERIC_RANKING_CACHE[key] = (rows, timing)
    return rows, timing


def specialized_ranking(
    model: Sequence[Fraction], generic_points: Sequence[tuple[Fraction, Fraction]]
) -> tuple[Any, dict[int, tuple[int, ...]], dict[int, Decimal], tuple[int, ...], dict[str, Any]]:
    started_wall = time.monotonic()
    started_cpu = cpu_clock()
    height_gram = engine.canonical_height_gram(model, generic_points)
    runs = {}
    representative_by_mask = {}
    depth_by_mask = {}
    for scale in (SPECIALIZED_CHECK_SCALE, SPECIALIZED_SCALE):
        rounded = tuple(
            tuple(int((value * Decimal(scale)).to_integral_value()) for value in row)
            for row in height_gram
        )
        oracle = engine.CosetOracle(rounded)
        top_rows = []
        all_rows = [] if scale == SPECIALIZED_SCALE else None
        for mask in range(1 << DIMENSION):
            unused_norm, representative, unused_error = oracle.solve(mask)
            depth = engine.quadratic_decimal(height_gram, representative) / 4
            if all_rows is not None:
                representative_by_mask[mask] = representative
                depth_by_mask[mask] = depth
                all_rows.append((depth, mask))
            else:
                top_rows.append((depth, mask))
        rows = all_rows if all_rows is not None else top_rows
        rows.sort(key=lambda row: (-row[0], row[1]))
        runs[scale] = tuple(mask for unused_depth, mask in rows[:ARM_SIZE])
    return (
        height_gram,
        representative_by_mask,
        depth_by_mask,
        runs[SPECIALIZED_SCALE],
        {
            "cpu_seconds": cpu_clock() - started_cpu,
            "wall_seconds": time.monotonic() - started_wall,
            "scales": [SPECIALIZED_CHECK_SCALE, SPECIALIZED_SCALE],
            "top43_order_stable": runs[SPECIALIZED_CHECK_SCALE] == runs[SPECIALIZED_SCALE],
            "top43_set_stable": set(runs[SPECIALIZED_CHECK_SCALE])
            == set(runs[SPECIALIZED_SCALE]),
        },
    )


def construct_arms(
    generic_rows: Sequence[tuple[int, int]], specialized_top: Sequence[int]
) -> tuple[dict[str, tuple[int, ...]], dict[str, Any]]:
    descending = sorted(generic_rows, key=lambda row: (-row[0], row[1]))
    generic_top = tuple(mask for unused_norm, mask in descending[:ARM_SIZE])
    maximum_norm = descending[0][0]
    maximum_norm_stratum_count = sum(
        norm == maximum_norm for norm, unused_mask in generic_rows
    )
    if len(generic_top) != ARM_SIZE:
        raise ArithmeticError("generic deepest arm does not contain 43 classes")

    positive = sorted((norm, mask) for norm, mask in generic_rows if mask)
    minimum_positive_norm = positive[0][0]
    median_norm = positive[len(positive) // 2][0]
    shallow = choose_hash_sample(
        (mask for norm, mask in positive if norm == minimum_positive_norm),
        SHALLOW_TIE_SEED,
    )
    median = choose_hash_sample(
        (mask for norm, mask in positive if norm == median_norm), MEDIAN_TIE_SEED
    )
    arms = {
        "generic-deepest43": tuple(sorted(generic_top)),
        "specialized-deepest43": tuple(sorted(specialized_top)),
        "deep-union": tuple(sorted(set(generic_top) | set(specialized_top))),
        **fixed_random_arms(),
        "generic-shallowest43": shallow,
        "generic-median43": median,
    }
    return arms, {
        "maximum_generic_norm": maximum_norm,
        "maximum_generic_norm_stratum_count": maximum_norm_stratum_count,
        "generic_deepest_tie_rule": "largest norm, then smallest mask",
        "minimum_positive_generic_norm": minimum_positive_norm,
        "median_nonzero_generic_norm": median_norm,
        "generic_specialized_intersection_count": len(set(generic_top) & set(specialized_top)),
    }


def generic_relation_filter(
    model,
    generic_points,
    candidates,
    *,
    chunk_size: int,
    timeout_seconds: float,
    stack_bytes: int,
):
    retained = []
    failed_chunks = []
    exact_count = 0
    for start in range(0, len(candidates), chunk_size):
        chunk = candidates[start : start + chunk_size]
        try:
            proposals = relation_proposals(
                model,
                generic_points,
                chunk,
                timeout=timeout_seconds,
                stack_bytes=stack_bytes,
            )
        except Exception as error:
            failed_chunks.append(
                {
                    "start": start,
                    "count": len(chunk),
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            )
            retained.extend((point, "filter_failed") for point in chunk)
            continue
        for point, (unused_relation, exact) in zip(chunk, proposals):
            if exact:
                exact_count += 1
            else:
                retained.append((point, "no_exact_generic_relation"))
    return retained, exact_count, failed_chunks


def checkpoint_payload(
    *, phase: str, args, results: Sequence[dict[str, Any]], status: str
) -> dict[str, Any]:
    inputs = {str(ENGINE.relative_to(ROOT)): digest(ENGINE)}
    if phase == "development":
        inputs.update(
            {
                str(MODEL.relative_to(ROOT)): digest(MODEL),
                str(SECTIONS.relative_to(ROOT)): digest(SECTIONS),
            }
        )
    else:
        inputs[str(RANK29_INPUT.relative_to(ROOT))] = digest(RANK29_INPUT)
    inputs[str(Path(__file__).resolve().relative_to(ROOT))] = digest(Path(__file__).resolve())
    return {
        "schema": "elliptic-curves.half-lattice-search-ablation-blind.v1",
        "status": status,
        "phase": phase,
        "blindness_boundary": {
            "exceptional_point_fixture_loaded": False,
            "development_parameters_used_without_public_points": phase == "development",
            "holdout_input_contains_only_curve_and_generic_subgroup": phase == "holdout",
            "arm_rules_and_hash_seeds_are_source_frozen": True,
        },
        "source_history": {
            "development_completed_source_sha256": DEVELOPMENT_FROZEN_SOURCE_SHA256,
            "holdout_presearch_fix": (
                "Removed only the assertion that the maximum generic-depth stratum has "
                "exactly 43 elements. The already-active norm-then-mask top-43 rule, all "
                "hash seeds, representatives, and search budgets are unchanged. The "
                "alternate-Q80 holdout had 49 maximum-depth classes and failed before "
                "any holdout quartic search."
            ),
        },
        "arm_definition": {
            "arm_size": ARM_SIZE,
            "random_arm_count": RANDOM_ARM_COUNT,
            "random_domain_seed": RANDOM_DOMAIN_SEED,
            "random_sets_are_disjoint_sha256_permutation_chunks": True,
            "shallow_tie_seed": SHALLOW_TIE_SEED,
            "median_tie_seed": MEDIAN_TIE_SEED,
            "representative_policy_for_every_arm": (
                "shortest representative at specialized height scale 10^6"
            ),
        },
        "declared_budget": {
            "height_bound_each_cover": args.height_bound,
            "timeout_seconds_each_cover": args.timeout_seconds,
            "stack_bytes_each_cover": args.stack_bytes,
            "generic_relation_filter_chunk_size": args.relation_chunk_size,
            "generic_relation_filter_timeout_each_chunk": args.relation_timeout_seconds,
            "single_search_per_distinct_cover": True,
            "shared_cover_outcome_charged_to_each_containing_arm": True,
        },
        "input_hashes": inputs,
        "completed_case_count": len(results),
        "results": list(results),
        "claim_boundary": [
            "All quartic-to-curve maps and returned curve points are checked exactly.",
            "Generic CVP depths are exact; specialized canonical-height CVPs are numerical.",
            "The blind artifact contains no target quotient coordinates or exceptional points.",
            "Point-search misses are bounded by the declared budget.",
            "CPU seconds include parent and completed child user-plus-system time per cover.",
            "Exact target-relative quotient ranks are supplied only by the separate verifier.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            f"elliptic-curves/cas/replay_half_lattice_search_ablation.sage --phase {phase}"
        ),
    }


def run_case(case: dict[str, Any], args) -> dict[str, Any]:
    label = case["label"]
    model = case["short_model"]
    generic_points = case["generic_points"]
    signatures = find_mod2_reduction_certificate(model, generic_points, prime_bound=1000)
    if combined_mod2_rank(signatures, DIMENSION) != DIMENSION:
        raise ArithmeticError(f"{label}: generic subgroup lost exact independence")

    generic_rows, generic_timing = generic_ranking(case["generic_gram"])
    (
        unused_height_gram,
        representatives,
        specialized_depths,
        specialized_top,
        specialized_timing,
    ) = specialized_ranking(model, generic_points)
    arms, selection = construct_arms(generic_rows, specialized_top)
    generic_norms = {mask: norm for norm, mask in generic_rows}
    selected_masks = tuple(sorted({mask for masks in arms.values() for mask in masks}))

    discoveries: dict[tuple[Fraction, Fraction], set[int]] = {}
    cover_records = []
    for position, mask in enumerate(selected_masks, 1):
        started_cpu = cpu_clock()
        outcome = engine.run_quartic_search(
            mask=mask,
            representative=representatives[mask],
            short_model=model,
            generic_points=generic_points,
            height_bound=args.height_bound,
            timeout_seconds=args.timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        cpu_seconds = cpu_clock() - started_cpu
        for point in outcome.curve_points:
            discoveries.setdefault(point, set()).add(mask)
        record = outcome.record
        cover_records.append(
            {
                "mask": mask,
                "hex": f"0x{mask:05x}",
                "generic_depth": str(Fraction(generic_norms[mask], 4)),
                "specialized_depth": str(specialized_depths[mask]),
                "specialized_representative": list(representatives[mask]),
                "status": record["status"],
                "cpu_seconds": cpu_seconds,
                "wall_seconds": record["wall_seconds"],
                "finite_curve_point_count": len(outcome.curve_points),
                "integral_coefficient_bits": record[
                    "integral_model_maximum_coefficient_bits"
                ],
                "reduced_coefficient_bits": (
                    record["reduced_model"]["maximum_coefficient_bits"]
                    if record["status"] == "bounded_search_complete"
                    else None
                ),
                "search_milliseconds": record.get("search_milliseconds"),
            }
        )
        print(
            f"HALFABLATE|case={label}|cover={position}/{len(selected_masks)}|"
            f"mask={mask:#07x}|status={record['status']}|points={len(outcome.curve_points)}|"
            f"cpu={cpu_seconds:.3f}",
            flush=True,
        )

    basis_signs = {
        signed for point in generic_points for signed in (point, (point[0], -point[1]))
    }
    nonbasis = tuple(sorted((point for point in discoveries if point not in basis_signs), key=point_key))
    retained, exact_generic_count, failed_chunks = generic_relation_filter(
        model,
        generic_points,
        nonbasis,
        chunk_size=args.relation_chunk_size,
        timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    candidate_rows = [
        {
            "point": point_record(point),
            "source_masks": sorted(discoveries[point]),
            "source_hex": [f"0x{mask:05x}" for mask in sorted(discoveries[point])],
            "generic_relation_filter": filter_status,
        }
        for point, filter_status in retained
    ]
    mask_to_record = {row["mask"]: row for row in cover_records}
    arm_rows = []
    for arm_id, masks in arms.items():
        mask_set = set(masks)
        candidate_indices = [
            index
            for index, row in enumerate(candidate_rows)
            if mask_set.intersection(row["source_masks"])
        ]
        nonbasis_count = sum(
            any(mask in mask_set for mask in source_masks)
            for point, source_masks in discoveries.items()
            if point not in basis_signs
        )
        arm_rows.append(
            {
                "id": arm_id,
                "class_count": len(masks),
                "masks": list(masks),
                "hex": [f"0x{mask:05x}" for mask in masks],
                "cover_cpu_seconds": sum(mask_to_record[mask]["cpu_seconds"] for mask in masks),
                "cover_wall_seconds": sum(mask_to_record[mask]["wall_seconds"] for mask in masks),
                "distinct_nonbasis_candidate_count": nonbasis_count,
                "retained_candidate_count": len(candidate_indices),
                "candidate_point_indices": candidate_indices,
                "generic_deepest_overlap": len(
                    set(masks) & set(arms["generic-deepest43"])
                ),
                "specialized_deepest_overlap": len(
                    set(masks) & set(arms["specialized-deepest43"])
                ),
            }
        )
    return {
        "label": label,
        "parameter": case["parameter"],
        "projective_parameter": case["projective_parameter"],
        "lineage": case["lineage"],
        "short_model": [str(value) for value in model],
        "generic_points": [point_record(point) for point in generic_points],
        "generic_subgroup_mod2_rank": DIMENSION,
        "generic_subgroup_certificate_primes": [row.prime for row in signatures],
        "ranking": {
            "generic": generic_timing,
            "specialized": specialized_timing,
            **selection,
        },
        "searched_distinct_class_count": len(selected_masks),
        "cover_records": cover_records,
        "candidate_filter": {
            "distinct_nonbasis_candidate_count": len(nonbasis),
            "exact_generic_relations_removed": exact_generic_count,
            "retained_candidate_count": len(candidate_rows),
            "failed_chunk_count": len(failed_chunks),
            "failed_chunks": failed_chunks,
        },
        "candidate_points": candidate_rows,
        "arms": arm_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("development", "holdout"), required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=128)
    parser.add_argument("--relation-timeout-seconds", type=float, default=60.0)
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid quartic search budget")
    if args.relation_chunk_size <= 0 or not 0 < args.relation_timeout_seconds <= 300:
        raise SystemExit("invalid relation-filter budget")
    output = args.output or (
        DEVELOPMENT_OUTPUT if args.phase == "development" else HOLDOUT_OUTPUT
    )
    cases = r17_cases() if args.phase == "development" else rank29_cases()
    results = []
    for index, case in enumerate(cases, 1):
        print(
            f"HALFABLATE|phase={args.phase}|case={index}/{len(cases)}|"
            f"label={case['label']}|status=START",
            flush=True,
        )
        results.append(run_case(case, args))
        payload = checkpoint_payload(
            phase=args.phase, args=args, results=results, status="PARTIAL_CHECKPOINT"
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(
            f"HALFABLATE|phase={args.phase}|case={case['label']}|status=CHECKPOINT",
            flush=True,
        )
    payload = checkpoint_payload(
        phase=args.phase, args=args, results=results, status="PASS_BLIND_ABLATION_SEARCH"
    )
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"HALFABLATE|phase={args.phase}|status=PASS|cases={len(results)}|"
        f"output={output.relative_to(ROOT)}",
        flush=True,
    )


engine = SourceFileLoader("half_lattice_ablation_engine", str(ENGINE)).load_module()


if __name__ == "__main__":
    main()
