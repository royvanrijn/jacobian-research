#!/usr/bin/env python3
"""Held-out-prime rerank of bounded q12/orbit5867 Nagao survivors.

The input population is fixed by a prior staged scan.  This script constructs
fresh complete P^1(F_p) tables on a disjoint prime interval, scores the fixed
population by lookup, and reports discovery/holdout rank overlap.  It does not
evaluate sections, search for points, or provide rank evidence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import sys
from time import perf_counter

from search_h92_q12o5867_rootless_nagao import (
    Candidate,
    SCORE_SCALE,
    build_residue_tables,
    is_prime,
    load_family_model,
    local_symbol_record,
    score_block,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/local/elkies-k3/q12o5867-rootless-nagao-cpp-h10000.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/local/elkies-k3/q12o5867-rootless-nagao-h10000-holdout-p199-499.json"
)


def load_population(path: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    document = json.loads(path.read_text())
    if document.get("status") not in {
        "PASS_BOUNDED_HEURISTIC_PROJECTIVE_NAGAO_CPP_SIEVE",
        "PASS_BOUNDED_HEURISTIC_PROJECTIVE_NAGAO_SIEVE",
    }:
        raise ValueError("input is not a bounded q12/orbit5867 Nagao sieve output")
    records = document.get("finalists")
    if not isinstance(records, list) or not records:
        raise ValueError("input has no finalist population")
    if len(records) != document.get("final_survivor_count"):
        raise ValueError("input finalists are truncated; reranking requires every survivor")
    seen: set[tuple[int, int]] = set()
    for record in records:
        pair = tuple(record["projective_pair"])
        if len(pair) != 2 or pair in seen:
            raise ValueError("input projective pairs are malformed or duplicated")
        seen.add(pair)
    return document, records


def rerank_population(
    discovery_records: list[dict[str, object]],
    tables: dict[int, tuple[object, ...]],
) -> tuple[list[dict[str, object]], list[str], list[str]]:
    inverse_cache: dict[tuple[int, int], int | None] = {}
    records: list[dict[str, object]] = []
    for discovery_rank, source in enumerate(discovery_records, start=1):
        numerator, denominator = map(int, source["projective_pair"])
        height = int(source["projective_height"])
        held = score_block(
            Candidate(numerator=numerator, denominator=denominator, height=height),
            tables,
            inverse_cache,
        )
        records.append(
            {
                "parameter": source["parameter"],
                "projective_pair": [numerator, denominator],
                "projective_height": height,
                "discovery_rank": discovery_rank,
                "discovery_total_score_units_1e12": int(
                    source["total_score_units_1e12"]
                ),
                "discovery_total_score": int(source["total_score_units_1e12"])
                / SCORE_SCALE,
                "discovery_block_score_units_1e12": list(
                    source["block_score_units_1e12"]
                ),
                "holdout_score_units_1e12": held.total_score_units,
                "holdout_score": held.total_score_units / SCORE_SCALE,
                "holdout_good_prime_count": held.good_primes,
                "holdout_bad_reduction_prime_count": held.bad_primes,
            }
        )

    holdout_order = sorted(
        records,
        key=lambda record: (
            -record["holdout_score_units_1e12"],
            -record["holdout_good_prime_count"],
            record["holdout_bad_reduction_prime_count"],
            record["projective_height"],
            record["projective_pair"][1],
            record["projective_pair"][0],
        ),
    )
    for holdout_rank, record in enumerate(holdout_order, start=1):
        record["holdout_rank"] = holdout_rank
        record["worst_of_discovery_holdout_rank"] = max(
            record["discovery_rank"], holdout_rank
        )
        record["discovery_holdout_rank_sum"] = record["discovery_rank"] + holdout_rank

    robust_order = sorted(
        records,
        key=lambda record: (
            record["worst_of_discovery_holdout_rank"],
            record["discovery_holdout_rank_sum"],
            record["projective_height"],
            record["projective_pair"][1],
            record["projective_pair"][0],
        ),
    )
    return (
        records,
        [record["parameter"] for record in holdout_order],
        [record["parameter"] for record in robust_order],
    )


def overlap_summary(
    discovery_order: list[str], holdout_order: list[str], cutoffs: tuple[int, ...]
) -> list[dict[str, object]]:
    summary = []
    for cutoff in cutoffs:
        bounded = min(cutoff, len(discovery_order))
        discovery_top = set(discovery_order[:bounded])
        shared = [parameter for parameter in holdout_order[:bounded] if parameter in discovery_top]
        summary.append(
            {
                "top_cutoff": bounded,
                "overlap_count": len(shared),
                "overlap_parameters_in_holdout_order": shared,
            }
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rerank q12/orbit5867 Nagao survivors on disjoint held-out primes."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--prime-min", type=int, default=199)
    parser.add_argument("--prime-max", type=int, default=499)
    args = parser.parse_args()
    if args.prime_min < 5 or args.prime_max < args.prime_min:
        raise SystemExit("invalid held-out prime interval")

    started = perf_counter()
    discovery, discovery_records = load_population(args.input)
    discovery_primes = {
        int(prime) for stage in discovery["stages"] for prime in stage["primes"]
    }
    requested_primes = tuple(
        prime
        for prime in range(args.prime_min, args.prime_max + 1)
        if is_prime(prime)
    )
    overlap = discovery_primes.intersection(requested_primes)
    if overlap:
        raise ValueError(f"held-out interval overlaps discovery primes: {sorted(overlap)}")

    model = load_family_model()
    table_started = perf_counter()
    blocks, rejected = build_residue_tables(model, (requested_primes,))
    table_seconds = perf_counter() - table_started
    tables = blocks[0]
    records, holdout_order, robust_order = rerank_population(discovery_records, tables)
    discovery_order = [record["parameter"] for record in records]
    runtime = perf_counter() - started

    document = {
        "schema": "h92-q12o5867-rootless-nagao-heldout-rerank-v1",
        "status": "PASS_BOUNDED_HEURISTIC_DISJOINT_PRIME_RERANK",
        "proof_boundary": (
            "Discovery and held-out Nagao scores are heuristics, not rank evidence. "
            "No section or specialized rational point was evaluated."
        ),
        "input": {
            "path": str(args.input.resolve()),
            "schema": discovery["schema"],
            "population": len(records),
            "discovery_primes": sorted(discovery_primes),
        },
        "model_sha256": model.source_sha256,
        "holdout": {
            "requested_interval": [args.prime_min, args.prime_max],
            "requested_primes": list(requested_primes),
            "usable_primes": list(tables),
            "rejected_primes": list(rejected),
            "disjoint_from_discovery": True,
            "nagao_formula": "((2-a_p)/(p+1-a_p))*log(p)",
            "integer_scale": SCORE_SCALE,
        },
        "local_tables": {
            str(prime): [local_symbol_record(symbol) for symbol in table]
            for prime, table in tables.items()
        },
        "candidates_in_discovery_order": records,
        "holdout_order": holdout_order,
        "robust_order_by_worst_then_sum_rank": robust_order,
        "top_overlap": overlap_summary(
            discovery_order, holdout_order, (10, 25, 50, 100)
        ),
        "runtime": {
            "table_seconds": table_seconds,
            "total_seconds": runtime,
        },
        "reproducing_command": shlex.join(sys.argv),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS population={len(records)} usable_primes={len(tables)} "
        f"top_holdout={holdout_order[0]} top_robust={robust_order[0]} "
        f"seconds={runtime:.3f} output={args.output}"
    )


if __name__ == "__main__":
    main()
