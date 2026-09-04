#!/usr/bin/env python3
"""Merge target-free A1/MW16 parameter-search checkpoints fail-closed."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elliptic-curves/a1_mw16_target_free_parameter_candidates_h300_v1.json"
SHARD_DIR = ROOT / "artifacts/local/elliptic-curves/a1-mw16-target-free-parameter-search-direct"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/a1_mw16_target_free_parameter_search_h300_v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_digest(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def chart_status(record) -> str:
    search = record.get("search")
    if not isinstance(search, dict) or not isinstance(search.get("status"), str):
        raise ArithmeticError("cover record has no search status")
    return search["status"]


def compact_result(result):
    quotient = result.get("exact_quotient_rank_recovered")
    if quotient is None or int(quotient) >= 1:
        return {"retention": "full_positive_or_fail_closed", "result": result}
    if result.get("status") != "PASS_COMPLETE_INITIAL_HALF_LATTICE_WAVE":
        raise ArithmeticError("zero result does not have the expected exact gate status")
    covers = result.get("cover_records", [])
    counts = Counter(chart_status(row) for row in covers)
    deepest = result.get("generic_half_lattice", {})
    if len(covers) != int(deepest.get("deepest_class_count", -1)):
        raise ArithmeticError("chart count does not equal exact deepest-stratum size")
    if result.get("discoveries", []):
        raise ArithmeticError("zero quotient result unexpectedly retains discoveries")
    saturation = result.get("discovered_group_saturation", {})
    if saturation.get("status") != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        raise ArithmeticError("zero quotient result lacks exact group classification")
    return {
        "retention": "compact_zero_discovery",
        "candidate_id": result["candidate_id"],
        "presentation_id": result["presentation_id"],
        "fibration_id": result["fibration_id"],
        "parameter": result["parameter"],
        "q_isomorphism_class_id": result["q_isomorphism_class_id"],
        "nagao": result["nagao"],
        "status": result["status"],
        "generic_mod2_independence_rank": result["generic_mod2_independence_rank"],
        "deepest_class_count": deepest["deepest_class_count"],
        "deepest_masks": deepest["deepest_masks"],
        "specialized_ranking": result["specialized_ranking"],
        "chart_status_counts": dict(sorted(counts.items())),
        "exact_quotient_rank_recovered": 0,
        "basis_rank_after": result["basis_rank_after"],
        "discovered_group_status": saturation["status"],
    }


def build(args):
    source = json.loads(args.source.read_text())
    if source.get("status") != "PASS_TARGET_FREE_A1_MW16_PARAMETER_CANDIDATES":
        raise ArithmeticError("source parameter ledger is not passing")
    candidates = source.get("candidates", [])
    expected_ids = [row["candidate_id"] for row in candidates]
    if len(expected_ids) != 104 or len(expected_ids) != len(set(expected_ids)):
        raise ArithmeticError("source is not the frozen 104-candidate ledger")
    candidate_by_id = {row["candidate_id"]: row for row in candidates}

    shard_paths = sorted(args.shard_dir.glob("shard-*.json"))
    if len(shard_paths) != args.expected_shards:
        raise ArithmeticError(f"expected {args.expected_shards} shards, found {len(shard_paths)}")
    results = []
    shard_source_hashes = set()
    budgets = []
    for path in shard_paths:
        shard = json.loads(path.read_text())
        if shard.get("status") != "PASS_TARGET_FREE_A1_MW16_PARAMETER_SEARCH_SHARD":
            raise ArithmeticError(f"shard is not complete: {relative(path)}")
        if shard.get("completed_candidate_count") != len(shard.get("candidate_ids", [])):
            raise ArithmeticError(f"shard completion count changed: {relative(path)}")
        if [row["candidate_id"] for row in shard["results"]] != shard["candidate_ids"]:
            raise ArithmeticError(f"shard result order changed: {relative(path)}")
        source_key = relative(args.source)
        old_hash = shard.get("inputs", {}).get(source_key)
        if not old_hash:
            raise ArithmeticError(f"shard lacks source provenance: {relative(path)}")
        shard_source_hashes.add(old_hash)
        budgets.append(shard["declared_budget"])
        results.extend(shard["results"])
    if len(shard_source_hashes) != 1:
        raise ArithmeticError("shards do not share one source-ledger hash")
    if any(row != budgets[0] for row in budgets[1:]):
        raise ArithmeticError("shard budgets differ")
    if [row["candidate_id"] for row in results] != expected_ids:
        raise ArithmeticError("shards do not form the exact candidate partition")

    positives = []
    failed_closed = []
    chart_counts = Counter()
    deepest_counts = Counter()
    signed_affine_points_reported = 0
    finite_curve_points_reported = 0
    search_milliseconds = []
    per_presentation = defaultdict(lambda: Counter(candidates=0, positives=0, failed_closed=0))
    per_fibration = defaultdict(lambda: Counter(candidates=0, positives=0, failed_closed=0))
    compact = []
    for result in results:
        candidate = candidate_by_id[result["candidate_id"]]
        for key in (
            "presentation_id", "fibration_id", "parameter",
            "q_isomorphism_class_id", "nagao",
        ):
            if result.get(key) != candidate.get(key):
                raise ArithmeticError(f"candidate/result mismatch for {key}")
        quotient = result.get("exact_quotient_rank_recovered")
        counts = (
            per_presentation[result["presentation_id"]],
            per_fibration[result["fibration_id"]],
        )
        for counter in counts:
            counter["candidates"] += 1
        if quotient is None:
            failed_closed.append(result["candidate_id"])
            for counter in counts:
                counter["failed_closed"] += 1
        elif int(quotient) >= 1:
            positives.append(result["candidate_id"])
            for counter in counts:
                counter["positives"] += 1
        deepest_counts[result.get("generic_half_lattice", {}).get("deepest_class_count")] += 1
        for cover in result.get("cover_records", []):
            search = cover["search"]
            chart_counts[chart_status(cover)] += 1
            signed_affine_points_reported += int(
                search.get("signed_affine_points_reported", 0)
            )
            finite_curve_points_reported += len(
                search.get("finite_curve_points", [])
            )
            if "search_milliseconds" in search:
                search_milliseconds.append(int(search["search_milliseconds"]))
        compact.append(compact_result(result))

    return {
        "schema": "elliptic-curves.a1-mw16-target-free-parameter-search-merged.v1",
        "status": "PASS_COMPLETE_TARGET_FREE_A1_MW16_PARAMETER_SEARCH",
        "source_candidate_count": len(candidates),
        "completed_candidate_count": len(results),
        "positive_candidate_ids": positives,
        "positive_candidate_count": len(positives),
        "failed_closed_candidate_ids": failed_closed,
        "failed_closed_candidate_count": len(failed_closed),
        "chart_status_counts": dict(sorted(chart_counts.items())),
        "deepest_class_count_histogram": {
            str(key): value for key, value in sorted(deepest_counts.items())
        },
        "signed_affine_points_reported": signed_affine_points_reported,
        "finite_curve_points_reported": finite_curve_points_reported,
        "quartic_search_milliseconds": {
            "total": sum(search_milliseconds),
            "maximum": max(search_milliseconds),
        },
        "declared_budget": budgets[0],
        "candidate_ledger_sha256": canonical_digest(candidates),
        "source_ledger_sha256_at_merge": digest(args.source),
        "source_ledger_sha256_recorded_by_shards": next(iter(shard_source_hashes)),
        "source_candidate_content_matches_shards": True,
        "per_presentation": {
            key: dict(sorted(value.items()))
            for key, value in sorted(per_presentation.items())
        },
        "per_fibration": {
            key: dict(sorted(value.items()))
            for key, value in sorted(per_fibration.items())
        },
        "results": compact,
        "selmer_gate": {
            "positive_definition": "exact_quotient_rank_recovered >= 1",
            "authorized_candidate_ids": positives,
            "complete_residual_2_selmer_calls_required": len(positives),
            "expensive_continuation_authorized_before_complete_selmer": False,
        },
        "inputs": {
            relative(args.source): digest(args.source),
            relative(Path(__file__)): digest(Path(__file__)),
            **{relative(path): digest(path) for path in shard_paths},
        },
        "claim_boundary": [
            "The campaign selected and searched 104 anonymous A1/MW16 parameter fibres without known-record controls.",
            "An exact positive means a rational point survived exact group-law and finite-reduction independence checks.",
            "A zero is only the rank recovered inside the declared bounded attempts.",
            "A bounded_search_timeout chart did not exhaust its height box and is censored, not negative evidence.",
            "No adaptive quotient wave, unrestricted point search, or Selmer computation occurs here.",
        ],
        "reproducing_command": (
            "python3 elliptic-curves/cas/merge_a1_mw16_target_free_parameter_search_shards.py --check"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--shard-dir", type=Path, default=SHARD_DIR)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--expected-shards", type=int, default=8)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build(args)
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != encoded:
            raise SystemExit("FAIL: merged target-free certificate changed")
        print("PASS: merged target-free certificate is unchanged")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    print(
        "A1MW16MERGE|"
        f"completed={payload['completed_candidate_count']}|"
        f"positive={payload['positive_candidate_count']}|"
        f"failed_closed={payload['failed_closed_candidate_count']}|"
        f"output={relative(args.output)}|status={payload['status']}"
    )


if __name__ == "__main__":
    main()
