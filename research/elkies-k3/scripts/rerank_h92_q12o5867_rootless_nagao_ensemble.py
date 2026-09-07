#!/usr/bin/env python3
"""Fresh-prime ensemble rerank for q12/orbit5867 Nagao candidates.

The two stored CRT/Gauss populations are reconstructed from their exact beam
states (and checked against the pinned population hashes), then united with
the complete H=10000 finalist population.  Candidate selection uses only a
new prime interval, disjoint from every construction and validation prime in
the two CRT artifacts.

At each prime the local Nagao contribution is centered and standardized over
the good fibres of P^1(F_p).  Missing values at singular reductions are
imputed by the prime mean (standardized contribution zero), rather than being
rewarded or punished.  Round-robin prime blocks provide an ensemble stability
score.  This remains a heuristic rerank, not rank evidence.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from math import sqrt
from pathlib import Path
import shlex
import statistics
import sys
from time import perf_counter
from typing import Iterable, Sequence


SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parents[1]
sys.path.insert(0, str(SCRIPTS))

from construct_h92_q12o5867_rootless_nagao_crt import (  # noqa: E402
    short_representatives,
)
from search_h92_q12o5867_rootless_nagao import (  # noqa: E402
    LocalSymbol,
    build_residue_tables,
    is_prime,
    load_family_model,
)


DEFAULT_BOX = ROOT / "artifacts/local/elkies-k3/q12o5867-rootless-nagao-cpp-h10000.json"
DEFAULT_CRT = (
    ROOT / "artifacts/local/elkies-k3/q12o5867-rootless-nagao-crt-gauss.json",
    ROOT
    / "artifacts/local/elkies-k3/q12o5867-rootless-nagao-crt-gauss-c199-313-v503-997.json",
)
DEFAULT_SPECIALIZATIONS = (
    ROOT / "artifacts/local/elliptic-curves/q12o5867-specializations"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/q12o5867-fresh-prime-ensemble-shortlist-v1.json"
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def population_sha256(rows: Iterable[tuple[int, int, int, int]]) -> str:
    text = "".join(f"{a}/{b}\n" for a, b, _, _ in rows).encode()
    return sha256(text).hexdigest()


def reconstruct_crt_population(
    path: Path, old_pairs: set[tuple[int, int]]
) -> tuple[dict[str, object], tuple[tuple[int, int, int, int], ...]]:
    document = json.loads(path.read_text())
    if document.get("status") != "PASS_BOUNDED_HEURISTIC_PROJECTIVE_CRT_GAUSS_CONSTRUCTOR":
        raise ValueError(f"not a q12/orbit5867 CRT/Gauss artifact: {path}")
    bounds = document["bounds"]
    answers: dict[tuple[int, int], tuple[int, int, int, int]] = {}
    for state_index, state in enumerate(document["construction"]["beam_states"]):
        for numerator, denominator, height in short_representatives(
            int(state["finite_residue"]),
            int(state["finite_modulus"]),
            int(state["infinity_modulus"]),
            coefficient_radius=int(bounds["final_coefficient_radius"]),
        ):
            if height < int(bounds["minimum_height"]):
                continue
            key = numerator, denominator
            candidate = numerator, denominator, height, state_index
            previous = answers.get(key)
            if previous is None or state_index < previous[3]:
                answers[key] = candidate
    complete = tuple(
        sorted(
            answers.values(),
            key=lambda row: (row[2], abs(row[0]), row[1], row[0]),
        )
    )
    novel = tuple(row for row in complete if row[:2] not in old_pairs)
    expected = document["construction"]
    checks = {
        "enumerated_parameter_count": len(complete),
        "novel_parameter_count": len(novel),
        "novel_parameter_population_sha256": population_sha256(novel),
    }
    expected_checks = {
        "enumerated_parameter_count": int(
            expected["enumerated_parameter_count_before_prior_dedup"]
        ),
        "novel_parameter_count": int(expected["novel_parameter_count"]),
        "novel_parameter_population_sha256": expected[
            "novel_parameter_population_sha256"
        ],
    }
    if checks != expected_checks:
        raise ValueError(
            f"CRT population reconstruction mismatch for {path}: "
            f"observed={checks} expected={expected_checks}"
        )
    return document, novel


def load_box_population(path: Path) -> tuple[dict[str, object], tuple[tuple[int, int, int, int], ...]]:
    document = json.loads(path.read_text())
    finalists = document.get("finalists")
    if not isinstance(finalists, list) or len(finalists) != document.get(
        "final_survivor_count"
    ):
        raise ValueError("H=10000 input does not contain its complete finalist population")
    rows = tuple(
        (
            int(record["projective_pair"][0]),
            int(record["projective_pair"][1]),
            int(record["projective_height"]),
            index,
        )
        for index, record in enumerate(finalists)
    )
    if len({row[:2] for row in rows}) != len(rows):
        raise ValueError("H=10000 finalist population contains duplicate parameters")
    return document, rows


def previous_promotions(directory: Path) -> dict[str, list[str]]:
    answer: dict[str, list[str]] = defaultdict(list)
    if not directory.exists():
        return {}
    for path in sorted(directory.glob("*.json")):
        document = json.loads(path.read_text())
        parameter = document.get("parameter")
        if isinstance(parameter, dict):
            value = parameter.get("affine_value")
            if isinstance(value, str):
                answer[value].append(str(path.resolve()))
    return dict(answer)


def source_prime_sets(document: dict[str, object]) -> tuple[set[int], set[int]]:
    bounds = document["bounds"]
    construction_values = bounds.get(
        "all_construction_primes", bounds.get("discovery_primes")
    )
    validation_values = bounds.get(
        "validation_usable_primes", bounds.get("heldout_usable_primes")
    )
    if construction_values is None or validation_values is None:
        raise ValueError("CRT artifact does not identify its scoring primes")
    construction = set(map(int, construction_values))
    validation = set(map(int, validation_values))
    return construction, validation


def round_robin_blocks(primes: Sequence[int], count: int) -> tuple[tuple[int, ...], ...]:
    if count < 2 or count > len(primes):
        raise ValueError("ensemble block count must lie between 2 and the prime count")
    blocks = [[] for _ in range(count)]
    for index, prime in enumerate(primes):
        blocks[index % count].append(prime)
    return tuple(tuple(block) for block in blocks)


def prime_standardization(table: Sequence[object]) -> tuple[float, float]:
    values = [float(symbol.contribution_units) for symbol in table if symbol.good_reduction]
    if len(values) < 2:
        raise ValueError("a fresh-prime table has fewer than two good fibres")
    mean = statistics.fmean(values)
    variance = statistics.fmean((value - mean) ** 2 for value in values)
    if variance <= 0:
        raise ValueError("a fresh-prime table has zero Nagao variance")
    return mean, sqrt(variance)


def assign_competition_ranks(
    records: list[dict[str, object]], key: str, rank_key: str
) -> None:
    ordered = sorted(
        records,
        key=lambda row: (
            -float(row[key]),
            int(row["fresh_bad_reduction_prime_count"]),
            int(row["projective_height"]),
            int(row["projective_pair"][1]),
            int(row["projective_pair"][0]),
        ),
    )
    prior_value: float | None = None
    rank = 0
    for position, row in enumerate(ordered, 1):
        value = float(row[key])
        if prior_value is None or value != prior_value:
            rank = position
            prior_value = value
        row[rank_key] = rank


def score_population(
    population: Sequence[dict[str, object]],
    tables: dict[int, tuple[object, ...]],
    blocks: Sequence[Sequence[int]],
) -> list[dict[str, object]]:
    prime_to_block = {
        prime: block_index
        for block_index, block in enumerate(blocks)
        for prime in block
    }
    standards = {prime: prime_standardization(table) for prime, table in tables.items()}
    inverse_cache: dict[tuple[int, int], int] = {}
    records: list[dict[str, object]] = []
    for source in population:
        numerator, denominator = map(int, source["projective_pair"])
        block_z_sums = [0.0] * len(blocks)
        block_raw_units = [0] * len(blocks)
        block_bad = [0] * len(blocks)
        for prime, table in tables.items():
            if denominator % prime == 0:
                index = prime
            else:
                residue = denominator % prime
                cache_key = prime, residue
                inverse = inverse_cache.get(cache_key)
                if inverse is None:
                    inverse = pow(residue, -1, prime)
                    inverse_cache[cache_key] = inverse
                index = numerator % prime * inverse % prime
            symbol = table[index]
            block_index = prime_to_block[prime]
            if symbol.good_reduction:
                mean, standard_deviation = standards[prime]
                block_z_sums[block_index] += (
                    float(symbol.contribution_units) - mean
                ) / standard_deviation
                block_raw_units[block_index] += int(symbol.contribution_units)
            else:
                # Mean imputation is standardized value zero.  Keep the count
                # visible because this is a heuristic choice, not local data.
                block_bad[block_index] += 1
        block_signals = [
            block_z_sums[index] / sqrt(len(blocks[index]))
            for index in range(len(blocks))
        ]
        mean_signal = statistics.fmean(block_signals)
        signal_sd = statistics.stdev(block_signals)
        robust_lcb = mean_signal - signal_sd / sqrt(len(blocks))
        total_z = sum(block_z_sums) / sqrt(len(tables))
        record = dict(source)
        record.update(
            {
                "fresh_standardized_block_signals": block_signals,
                "fresh_standardized_total_z": total_z,
                "fresh_ensemble_mean": mean_signal,
                "fresh_ensemble_sample_sd": signal_sd,
                "fresh_ensemble_one_se_lcb": robust_lcb,
                "fresh_raw_nagao_block_score_units_1e12": block_raw_units,
                "fresh_raw_nagao_total_score_units_1e12": sum(block_raw_units),
                "fresh_bad_reduction_count_by_block": block_bad,
                "fresh_bad_reduction_prime_count": sum(block_bad),
                "fresh_good_reduction_prime_count": len(tables) - sum(block_bad),
            }
        )
        records.append(record)

    assign_competition_ranks(
        records, "fresh_ensemble_one_se_lcb", "fresh_ensemble_lcb_rank"
    )
    assign_competition_ranks(
        records, "fresh_standardized_total_z", "fresh_standardized_total_rank"
    )
    assign_competition_ranks(
        records,
        "fresh_raw_nagao_total_score_units_1e12",
        "fresh_raw_nagao_total_rank",
    )
    for omitted in range(len(blocks)):
        key = f"leave_block_{omitted + 1}_out_z"
        denominator = sqrt(len(tables) - len(blocks[omitted]))
        for row in records:
            z_sums = [
                signal * sqrt(len(blocks[index]))
                for index, signal in enumerate(
                    row["fresh_standardized_block_signals"]
                )
            ]
            row[key] = (sum(z_sums) - z_sums[omitted]) / denominator
        assign_competition_ranks(records, key, f"{key}_rank")
    for row in records:
        loo_ranks = [
            int(row[f"leave_block_{index + 1}_out_z_rank"])
            for index in range(len(blocks))
        ]
        row["fresh_leave_one_block_out_worst_rank"] = max(loo_ranks)
        row["fresh_leave_one_block_out_rank_sum"] = sum(loo_ranks)
    return records


def select_shortlist(
    records: list[dict[str, object]],
    previously_promoted: dict[str, list[str]],
    *,
    stable_count: int,
    total_count: int,
    source_count: int,
) -> list[dict[str, object]]:
    eligible = [
        row for row in records if row["parameter"] not in previously_promoted
    ]
    selected: dict[str, dict[str, object]] = {}

    def add(rows: Iterable[dict[str, object]], tier: str, limit: int) -> None:
        additions = 0
        for row in rows:
            parameter = str(row["parameter"])
            if parameter in selected:
                selected[parameter]["selection_tiers"].append(tier)
                continue
            if additions >= limit:
                break
            copied = dict(row)
            copied["selection_tiers"] = [tier]
            selected[parameter] = copied
            additions += 1

    stable_order = sorted(
        eligible,
        key=lambda row: (
            -float(row["fresh_ensemble_one_se_lcb"]),
            -float(row["fresh_standardized_total_z"]),
            int(row["fresh_leave_one_block_out_worst_rank"]),
            int(row["projective_height"]),
        ),
    )
    add(stable_order, "global_stable_lcb", stable_count)
    total_order = sorted(
        eligible,
        key=lambda row: (
            -float(row["fresh_standardized_total_z"]),
            -float(row["fresh_ensemble_one_se_lcb"]),
            int(row["fresh_leave_one_block_out_worst_rank"]),
            int(row["projective_height"]),
        ),
    )
    add(total_order, "global_fresh_total", total_count)
    sources = sorted({str(row["source_population"]) for row in eligible})
    for source in sources:
        source_order = [row for row in stable_order if row["source_population"] == source]
        add(source_order, f"source_stable_lcb:{source}", source_count)

    result = sorted(
        selected.values(),
        key=lambda row: (
            -float(row["fresh_ensemble_one_se_lcb"]),
            -float(row["fresh_standardized_total_z"]),
            int(row["projective_height"]),
        ),
    )
    for priority, row in enumerate(result, 1):
        row["promotion_priority"] = priority
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fresh-prime ensemble rerank of q12/orbit5867 candidates."
    )
    parser.add_argument("--box", type=Path, default=DEFAULT_BOX)
    parser.add_argument(
        "--additional-box",
        type=Path,
        action="append",
        default=[],
        help=(
            "additional complete bounded-sieve finalist population; may be "
            "repeated (the primary --box remains the CRT prior-dedup anchor)"
        ),
    )
    parser.add_argument("--crt", type=Path, action="append")
    parser.add_argument("--specializations", type=Path, default=DEFAULT_SPECIALIZATIONS)
    parser.add_argument(
        "--include-specialized",
        action="store_true",
        help=(
            "do not exclude parameters merely because an exact specialization "
            "artifact already exists; useful when reranking a newly enlarged population"
        ),
    )
    parser.add_argument(
        "--exclude-parameter",
        action="append",
        default=[],
        help=(
            "explicit previously promoted affine parameter; when supplied at "
            "least once, these values replace the live specialization-directory scan"
        ),
    )
    parser.add_argument("--prime-min", type=int, default=1009)
    parser.add_argument("--prime-max", type=int, default=1499)
    parser.add_argument("--blocks", type=int, default=6)
    parser.add_argument("--confirmation-prime-min", type=int, default=1511)
    parser.add_argument("--confirmation-prime-max", type=int, default=1999)
    parser.add_argument("--confirmation-blocks", type=int, default=6)
    parser.add_argument("--confirmed-count", type=int, default=16)
    parser.add_argument("--stable-count", type=int, default=24)
    parser.add_argument("--total-count", type=int, default=12)
    parser.add_argument("--source-count", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    crt_paths = tuple(args.crt) if args.crt else DEFAULT_CRT
    if args.prime_min < 5 or args.prime_max < args.prime_min:
        raise SystemExit("invalid fresh prime interval")
    if (
        args.confirmation_prime_min < 5
        or args.confirmation_prime_max < args.confirmation_prime_min
    ):
        raise SystemExit("invalid confirmation prime interval")
    if min(args.stable_count, args.total_count, args.source_count) < 0:
        raise SystemExit("shortlist tier counts must be nonnegative")
    if args.confirmed_count < 1:
        raise SystemExit("confirmed shortlist count must be positive")

    started = perf_counter()
    box_document, box_rows = load_box_population(args.box)
    old_pairs = {row[:2] for row in box_rows}
    crt_documents = []
    source_rows: list[tuple[str, tuple[tuple[int, int, int, int], ...]]] = [
        ("h10000_finalists", box_rows)
    ]
    additional_box_documents = []
    for index, path in enumerate(args.additional_box, 1):
        document, rows = load_box_population(path)
        additional_box_documents.append((path, document, rows))
        source_rows.append((f"additional_box_population_{index}", rows))
    occupied_primes: set[int] = set()
    for index, path in enumerate(crt_paths, 1):
        document, rows = reconstruct_crt_population(path, old_pairs)
        crt_documents.append((path, document))
        construction, validation = source_prime_sets(document)
        occupied_primes.update(construction)
        occupied_primes.update(validation)
        source_rows.append((f"crt_population_{index}", rows))

    population_map: dict[tuple[int, int], dict[str, object]] = {}
    source_counts = {}
    for source, rows in source_rows:
        source_counts[source] = len(rows)
        for numerator, denominator, height, source_index in rows:
            key = numerator, denominator
            parameter = f"{numerator}/{denominator}"
            previous = population_map.get(key)
            if previous is None:
                population_map[key] = {
                    "parameter": parameter,
                    "projective_pair": [numerator, denominator],
                    "projective_height": height,
                    "source_population": source,
                    "source_index": source_index,
                    "also_in_source_populations": [],
                }
            else:
                previous["also_in_source_populations"].append(source)
    population = sorted(
        population_map.values(),
        key=lambda row: (
            int(row["projective_height"]),
            abs(int(row["projective_pair"][0])),
            int(row["projective_pair"][1]),
            int(row["projective_pair"][0]),
        ),
    )
    population_text = "".join(f"{row['parameter']}\n" for row in population).encode()

    requested_primes = tuple(
        prime
        for prime in range(args.prime_min, args.prime_max + 1)
        if is_prime(prime)
    )
    overlap = occupied_primes.intersection(requested_primes)
    if overlap:
        raise ValueError(f"fresh prime interval overlaps prior CRT evidence: {sorted(overlap)}")
    model = load_family_model()
    table_started = perf_counter()
    table_blocks, rejected = build_residue_tables(model, (requested_primes,))
    table_seconds = perf_counter() - table_started
    tables = table_blocks[0]
    usable_primes = tuple(tables)
    blocks = round_robin_blocks(usable_primes, args.blocks)
    score_started = perf_counter()
    records = score_population(population, tables, blocks)
    score_seconds = perf_counter() - score_started
    promotions = (
        {}
        if args.include_specialized
        else (
            {
                parameter: ["explicit_command_line_exclusion"]
                for parameter in args.exclude_parameter
            }
            if args.exclude_parameter
            else previous_promotions(args.specializations)
        )
    )
    shortlist = select_shortlist(
        records,
        promotions,
        stable_count=args.stable_count,
        total_count=args.total_count,
        source_count=args.source_count,
    )

    confirmation_requested = tuple(
        prime
        for prime in range(
            args.confirmation_prime_min, args.confirmation_prime_max + 1
        )
        if is_prime(prime)
    )
    confirmation_overlap = (
        occupied_primes.union(usable_primes).intersection(confirmation_requested)
    )
    if confirmation_overlap:
        raise ValueError(
            "confirmation interval overlaps earlier evidence: "
            f"{sorted(confirmation_overlap)}"
        )
    confirmation_table_started = perf_counter()
    confirmation_table_blocks, confirmation_rejected = build_residue_tables(
        model, (confirmation_requested,)
    )
    confirmation_table_seconds = perf_counter() - confirmation_table_started
    confirmation_tables = confirmation_table_blocks[0]
    confirmation_usable = tuple(confirmation_tables)
    confirmation_blocks = round_robin_blocks(
        confirmation_usable, args.confirmation_blocks
    )
    confirmation_score_started = perf_counter()
    confirmation_rows = score_population(
        [dict(row) for row in shortlist],
        confirmation_tables,
        confirmation_blocks,
    )
    confirmation_score_seconds = perf_counter() - confirmation_score_started
    confirmation_by_parameter = {
        str(row["parameter"]): row for row in confirmation_rows
    }
    for row in shortlist:
        confirmation = confirmation_by_parameter[str(row["parameter"])]
        for key, value in confirmation.items():
            if key.startswith("fresh_"):
                row[f"confirmation_{key[6:]}"] = value
        row["combined_selection_confirmation_standardized_z"] = (
            float(row["fresh_standardized_total_z"]) * sqrt(len(tables))
            + float(confirmation["fresh_standardized_total_z"])
            * sqrt(len(confirmation_tables))
        ) / sqrt(len(tables) + len(confirmation_tables))
        row["passes_positive_confirmation_gate"] = (
            float(confirmation["fresh_standardized_total_z"]) > 0
            and float(confirmation["fresh_ensemble_one_se_lcb"]) > 0
        )
    confirmed = sorted(
        (row for row in shortlist if row["passes_positive_confirmation_gate"]),
        key=lambda row: (
            -float(row["confirmation_ensemble_one_se_lcb"]),
            -float(row["confirmation_standardized_total_z"]),
            -float(row["fresh_ensemble_one_se_lcb"]),
            int(row["projective_height"]),
        ),
    )[: args.confirmed_count]
    for priority, row in enumerate(confirmed, 1):
        row["confirmed_promotion_priority"] = priority

    output = {
        "schema": "h92-q12o5867-rootless-fresh-prime-ensemble-rerank-v1",
        "status": "PASS_BOUNDED_HEURISTIC_FRESH_PRIME_ENSEMBLE_SHORTLIST",
        "proof_boundary": (
            "Fresh-prime centered/standardized Nagao scores and ensemble stability "
            "are heuristic parameter-ranking evidence only. They do not prove a "
            "rank jump, a point, an upper bound, or parity."
        ),
        "model_sha256": model.source_sha256,
        "inputs": {
            "h10000": {
                "path": str(args.box.resolve()),
                "sha256": file_sha256(args.box),
                "population_count": len(box_rows),
            },
            "additional_boxes": [
                {
                    "path": str(path.resolve()),
                    "sha256": file_sha256(path),
                    "population_count": len(rows),
                }
                for path, _document, rows in additional_box_documents
            ],
            "crt": [
                {
                    "path": str(path.resolve()),
                    "sha256": file_sha256(path),
                    "reconstructed_novel_population_count": int(
                        document["construction"]["novel_parameter_count"]
                    ),
                    "reconstructed_population_sha256": document["construction"][
                        "novel_parameter_population_sha256"
                    ],
                }
                for path, document in crt_documents
            ],
            "specialization_directory": str(args.specializations.resolve()),
            "previously_promoted_parameters": promotions,
        },
        "population": {
            "source_counts_before_union_dedup": source_counts,
            "unique_count": len(population),
            "unique_ordered_sha256": sha256(population_text).hexdigest(),
            "cross_source_duplicate_count": sum(
                len(row["also_in_source_populations"]) for row in population
            ),
        },
        "fresh_primes": {
            "requested_interval": [args.prime_min, args.prime_max],
            "requested_primes": list(requested_primes),
            "usable_primes": list(usable_primes),
            "rejected_primes": list(rejected),
            "disjoint_from_all_crt_construction_and_validation_primes": True,
            "round_robin_blocks": [list(block) for block in blocks],
        },
        "confirmation_primes": {
            "requested_interval": [
                args.confirmation_prime_min,
                args.confirmation_prime_max,
            ],
            "requested_primes": list(confirmation_requested),
            "usable_primes": list(confirmation_usable),
            "rejected_primes": list(confirmation_rejected),
            "disjoint_from_all_prior_and_selection_primes": True,
            "round_robin_blocks": [list(block) for block in confirmation_blocks],
            "fixed_input_shortlist_count": len(shortlist),
        },
        "method": {
            "local_statistic": (
                "(Nagao contribution - P1 good-fibre mean) / P1 good-fibre "
                "population standard deviation"
            ),
            "singular_reduction_imputation": (
                "prime mean, hence standardized contribution zero; bad count retained"
            ),
            "block_signal": "sum standardized local statistics / sqrt(block prime count)",
            "stable_score": "mean(block signals) - sample_sd(block signals)/sqrt(block count)",
            "stability_audit": "competition ranks for every leave-one-block-out total",
            "selection": {
                "global_stable_lcb": args.stable_count,
                "global_fresh_total": args.total_count,
                "per_source_stable_lcb": args.source_count,
                "previously_promoted_excluded": not args.include_specialized,
            },
            "confirmation_gate": (
                "fixed selection shortlist rescored on a second disjoint interval; "
                "retain positive confirmation total z and positive confirmation "
                "one-SE ensemble lower score, then order by confirmation lower score"
            ),
        },
        "shortlist_count": len(shortlist),
        "shortlist": shortlist,
        "positive_confirmation_gate_count": sum(
            bool(row["passes_positive_confirmation_gate"]) for row in shortlist
        ),
        "confirmed_shortlist_count": len(confirmed),
        "confirmed_shortlist": confirmed,
        "top_previously_promoted_for_comparison": [
            row
            for row in sorted(
                records,
                key=lambda row: (
                    -float(row["fresh_ensemble_one_se_lcb"]),
                    -float(row["fresh_standardized_total_z"]),
                ),
            )
            if row["parameter"] in promotions
        ][:10],
        "runtime": {
            "table_seconds": table_seconds,
            "score_seconds": score_seconds,
            "confirmation_table_seconds": confirmation_table_seconds,
            "confirmation_score_seconds": confirmation_score_seconds,
            "total_seconds": perf_counter() - started,
        },
        "reproducing_command": shlex.join(sys.argv),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS population={len(population)} fresh_primes={len(tables)} "
        f"shortlist={len(shortlist)} confirmed={len(confirmed)} "
        f"top={confirmed[0]['parameter'] if confirmed else 'none'} "
        f"seconds={output['runtime']['total_seconds']:.3f} output={args.output}"
    )


if __name__ == "__main__":
    main()
