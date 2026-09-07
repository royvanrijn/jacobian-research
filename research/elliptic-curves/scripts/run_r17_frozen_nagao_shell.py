#!/usr/bin/env python3
"""Rank a prospective compact-R17 height shell with the frozen Nagao rule.

The target ordering is copied verbatim from the complete H<=10,000 positive-
control artifact: descending weakest standardized block, followed by its
existing tie-break.  This command may change only the population and the
amount of downstream work allocated by that ordering.  It also materializes
matched pooled-Nagao and deterministic-random control lanes.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
K3_SCRIPTS = ROOT / "elkies-k3/scripts"
sys.path.insert(0, str(K3_SCRIPTS))

from search_h92_q12o5867_rootless_nagao import (  # noqa: E402
    build_residue_tables,
    export_cpp_tables,
    load_family_model,
)


SCANNER = K3_SCRIPTS / "scan_h92_q12o5867_rootless_nagao.cpp"
FROZEN_REFERENCE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_compact_t_nagao_positive_control_h10000_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "r17_frozen_nagao_shell_h10001_30000_v1.json"
)
DEFAULT_COHORT = (
    ROOT / "artifacts/local/elliptic-curves/r17-frozen-shell-h10001-30000-cohort.jsonl"
)
DEFAULT_LOCAL = ROOT / "artifacts/local/elliptic-curves/r17-frozen-shell-h10001-30000"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def scoring_contract(reference: dict[str, object]) -> dict[str, object]:
    scoring = reference["scoring"]
    assert isinstance(scoring, dict)
    expected = {
        "primary_ranking_key": "minimum block signal",
        "tie_breaker": (
            "mean block signal, good primes, bad primes, height, denominator, numerator"
        ),
        "per_prime_standardization": (
            "center and population-standardize over good fibres of P1(F_p)"
        ),
        "singular_fibre_policy": "mean imputation (standardized contribution zero)",
        "block_normalization": "sum(z_p)/sqrt(number of primes in block)",
    }
    for key, value in expected.items():
        if scoring.get(key) != value:
            raise ValueError(f"frozen scoring field changed: {key}")
    ensembles = scoring.get("prime_ensembles")
    if not isinstance(ensembles, list) or len(ensembles) != 3:
        raise ValueError("the frozen rule no longer has exactly three blocks")
    flattened = [int(prime) for block in ensembles for prime in block]
    if len(flattened) != len(set(flattened)):
        raise ValueError("the frozen prime blocks are no longer disjoint")
    return {key: scoring[key] for key in (*expected, "prime_ensembles")}


def selected_rows(
    ranking: dict[str, object], *, target_count: int, control_count: int
) -> list[dict[str, object]]:
    target_pool = ranking["ranked_prefix"]
    ordinary_pool = ranking["ordinary_nagao_control_prefix"]
    random_pool = ranking["random_control_lane"]
    assert isinstance(target_pool, list)
    assert isinstance(ordinary_pool, list)
    assert isinstance(random_pool, list)
    chosen: list[dict[str, object]] = []
    seen: set[tuple[int, int]] = set()

    def add_lane(pool, lane: str, count: int) -> None:
        added = 0
        for source in pool:
            score = source["score"]
            pair = tuple(map(int, score["projective_pair"]))
            if pair in seen:
                continue
            seen.add(pair)
            added += 1
            frozen_rank = source.get("population_rank")
            chosen.append(
                {
                    "parameter": score["parameter"],
                    "projective_pair": list(pair),
                    "split": "prospective_shell_h10001_30000",
                    "selection_lanes": [lane],
                    "lane_rank": added,
                    "frozen_population_rank": frozen_rank,
                    "ranking_scores": {
                        "worst_block_signal": score["worst_block_signal"],
                        "mean_block_signal": score["mean_block_signal"],
                        "standardized_block_signals": score[
                            "standardized_block_signals"
                        ],
                    },
                    "search_depth_tier": (
                        "deepest" if lane == "frozen_weakest_block" and added <= 16
                        else "medium" if lane == "frozen_weakest_block" and added <= 64
                        else "base"
                    ),
                }
            )
            if added == count:
                return
        raise ValueError(f"the {lane} pool did not provide {count} disjoint rows")

    add_lane(target_pool, "frozen_weakest_block", target_count)
    add_lane(ordinary_pool, "ordinary_nagao_control", control_count)
    add_lane(random_pool, "random_control", control_count)
    return chosen


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--minimum-height", type=int, default=10_001)
    parser.add_argument("--maximum-height", type=int, default=30_000)
    parser.add_argument("--ranked-prefix-count", type=int, default=256)
    parser.add_argument("--control-pool-count", type=int, default=256)
    parser.add_argument("--target-search-count", type=int, default=128)
    parser.add_argument("--control-search-count", type=int, default=128)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--local-directory", type=Path, default=DEFAULT_LOCAL)
    args = parser.parse_args()
    if not (1 < args.minimum_height <= args.maximum_height):
        raise SystemExit("the height shell must be nonempty and exclude infinity")
    if min(
        args.ranked_prefix_count,
        args.control_pool_count,
        args.target_search_count,
        args.control_search_count,
    ) < 1:
        raise SystemExit("all prefix and cohort counts must be positive")
    if args.target_search_count > args.ranked_prefix_count:
        raise SystemExit("target search count exceeds the ranked prefix")
    if args.control_search_count * 2 > args.control_pool_count:
        raise SystemExit("control pool must be at least twice the searched control count")

    reference = json.loads(FROZEN_REFERENCE.read_text())
    if reference.get("status") != "PASS_POSITIVE_CONTROL_SCORING_GATE":
        raise SystemExit("the frozen H<=10,000 scoring gate is unavailable")
    reference_search = reference.get("search", {})
    if reference_search.get("height_limit") != 10_000:
        raise SystemExit("the frozen development population changed")
    if args.minimum_height <= int(reference_search["height_limit"]):
        raise SystemExit("the prospective shell overlaps score development")
    frozen = scoring_contract(reference)
    prime_blocks = tuple(tuple(map(int, block)) for block in frozen["prime_ensembles"])

    compiler = shutil.which("g++")
    if compiler is None:
        raise SystemExit("g++ is required")
    args.local_directory.mkdir(parents=True, exist_ok=True)
    table_path = args.local_directory / "frozen-tables.txt"
    binary_path = args.local_directory / "scan-frozen-region"
    raw_path = args.local_directory / "ranking-raw.json"
    model = load_family_model()
    tables, rejected = build_residue_tables(model, prime_blocks)
    if rejected:
        raise SystemExit(f"a frozen score prime became unusable: {rejected}")
    export_cpp_tables(table_path, model, tables)
    compile_command = [
        compiler,
        "-O3",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-pedantic",
        str(SCANNER),
        "-o",
        str(binary_path),
    ]
    scan_command = [
        str(binary_path),
        str(table_path),
        str(args.maximum_height),
        str(args.maximum_height),
        "1",
        "1,1,1",
        str(args.ranked_prefix_count),
        str(raw_path),
        "1",
        "--rank-region",
        str(args.minimum_height),
        str(args.control_pool_count),
    ]
    started = perf_counter()
    subprocess.run(compile_command, check=True)
    subprocess.run(scan_command, check=True)
    result = json.loads(raw_path.read_text())
    if result.get("status") != "PASS_COMPLETE_FROZEN_RULE_REGION_RANKING":
        raise SystemExit("the complete region scan did not pass")
    if result["frozen_ranking"]["prime_ensembles"] != frozen["prime_ensembles"]:
        raise AssertionError("the region scan changed the frozen prime windows")
    if result["frozen_ranking"]["primary"] != "minimum standardized block signal":
        raise AssertionError("the region scan changed the primary comparator")

    cohort = selected_rows(
        result,
        target_count=args.target_search_count,
        control_count=args.control_search_count,
    )
    args.cohort.parent.mkdir(parents=True, exist_ok=True)
    args.cohort.write_text(
        "".join(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n" for row in cohort)
    )
    result["prospective_experiment"] = {
        "score_development_population": "complete compact-t H<=10000 box",
        "prospective_population": (
            f"complete compact-t shell {args.minimum_height}<=H<={args.maximum_height}"
        ),
        "sets_are_disjoint_by_height": True,
        "ranking_rule_frozen_before_shell_opened": True,
        "permitted_optimization": "downstream search-depth allocation only",
        "depth_allocation": {
            "frozen ranks 1-16": "deepest",
            "frozen ranks 17-64": "medium",
            "frozen ranks 65-128": "base",
            "ordinary and random controls": "base",
        },
        "selected_cohort": {
            "path": relative(args.cohort),
            "sha256": digest(args.cohort),
            "target_count": args.target_search_count,
            "ordinary_control_count": args.control_search_count,
            "random_control_count": args.control_search_count,
            "lanes_are_parameter_disjoint": True,
        },
    }
    result["frozen_reference"] = {
        "path": relative(FROZEN_REFERENCE),
        "sha256": digest(FROZEN_REFERENCE),
        "scoring_contract": frozen,
    }
    result["orchestration"] = {
        "script": relative(Path(__file__)),
        "script_sha256": digest(Path(__file__)),
        "scanner": relative(SCANNER),
        "scanner_sha256": digest(SCANNER),
        "table_sha256": digest(table_path),
        "compile_command": shlex.join(compile_command),
        "scan_command": shlex.join(scan_command),
        "wall_seconds": perf_counter() - started,
        "raw_local_result": relative(raw_path),
        "raw_local_result_sha256": digest(raw_path),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        f"R17FROZENSHELL|population={result['population_count']}|"
        f"ranked={len(result['ranked_prefix'])}|cohort={len(cohort)}|"
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()
