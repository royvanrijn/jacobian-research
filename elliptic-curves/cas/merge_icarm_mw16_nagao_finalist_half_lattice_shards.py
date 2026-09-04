#!/usr/bin/env python3
"""Merge checkpointed MW16 finalist half-lattice shards fail-closed.

The raw checkpoints contain enormous rational base points and belong under
``artifacts/local``.  This script verifies the frozen candidate partition and
emits a compact generated certificate.  Full records are retained for every
positive or structurally rejected candidate; zero-discovery records retain
the exact gate result and chart completion/timeout accounting.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_specializations_h300_v1.json"
SHARD_DIR = ROOT / "artifacts/local/elliptic-curves/mw16-finalist-half-lattice"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_half_lattice_h300_v1.json"
EXPECTED_CANDIDATE_LEDGER_SHA256 = "1b408f67e44bd390d6643d9e21be706677e085a90238360fbef136529d2bc610"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_digest(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def strip_runtime(value):
    """Remove non-mathematical wall/CPU measurements before hashing."""

    if isinstance(value, dict):
        return {
            key: strip_runtime(item)
            for key, item in value.items()
            if key not in {"wall_seconds", "cpu_seconds", "runtime_seconds"}
        }
    if isinstance(value, list):
        return [strip_runtime(item) for item in value]
    return value


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def chart_status(record) -> str:
    search = record.get("search")
    if not isinstance(search, dict) or not isinstance(search.get("status"), str):
        raise ArithmeticError("cover record has no search status")
    return search["status"]


def rational_bits(value) -> int:
    value = Fraction(value)
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def compact_result(result, candidate):
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
    discoveries = result.get("discoveries", [])
    if discoveries:
        raise ArithmeticError("zero quotient result unexpectedly retains discoveries")
    saturation = result.get("discovered_group_saturation", {})
    if saturation.get("status") != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        raise ArithmeticError("zero quotient result lacks exact group classification")
    quartic_bits = [
        int(row["search"]["raw_rational_coefficient_maximum_bits"])
        for row in covers
    ]
    return {
        "retention": "compact_zero_discovery",
        "candidate_id": result["candidate_id"],
        "parent_id": result["parent_id"],
        "parameter": result["parameter"],
        "q_isomorphism_class_id": result["q_isomorphism_class_id"],
        "nagao": result["nagao"],
        "status": result["status"],
        "generic_mod2_independence_rank": result["generic_mod2_independence_rank"],
        "deepest_class_count": deepest["deepest_class_count"],
        "deepest_masks": deepest["deepest_masks"],
        "specialized_ranking": result["specialized_ranking"],
        "raw_short_model_maximum_bits": max(
            rational_bits(value) for value in candidate["raw_short_model"]
        ),
        "raw_quartic_coefficient_maximum_bits": max(quartic_bits),
        "chart_status_counts": dict(sorted(counts.items())),
        "exact_quotient_rank_recovered": 0,
        "basis_rank_after": result["basis_rank_after"],
        "discovered_group_status": saturation["status"],
    }


def build(args):
    source = json.loads(args.source.read_text())
    if source.get("status") != "PASS_EXACT_MW16_NAGAO_FINALIST_SPECIALIZATIONS":
        raise ArithmeticError("source specialization ledger is not passing")
    candidates = source.get("candidates", [])
    expected_ids = [row["candidate_id"] for row in candidates]
    if len(expected_ids) != 104 or len(expected_ids) != len(set(expected_ids)):
        raise ArithmeticError("source is not the frozen 104-candidate ledger")
    candidate_by_id = {row["candidate_id"]: row for row in candidates}
    candidate_ledger_digest = canonical_digest(candidates)
    if candidate_ledger_digest != EXPECTED_CANDIDATE_LEDGER_SHA256:
        raise ArithmeticError("frozen candidate content changed")

    shard_paths = sorted(args.shard_dir.glob("shard-*.json"))
    if len(shard_paths) != args.expected_shards:
        raise ArithmeticError(
            f"expected {args.expected_shards} shards, found {len(shard_paths)}"
        )
    results = []
    shard_result_groups = []
    shard_source_hashes = set()
    budgets = []
    for path in shard_paths:
        shard = json.loads(path.read_text())
        if shard.get("status") != "PASS_BOUNDED_PREFIX_NAGAO_FINALIST_HALF_LATTICE_GATE":
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
        shard_result_groups.append((path, shard["results"]))
    if len(shard_source_hashes) != 1:
        raise ArithmeticError("shards do not share one source-ledger hash")
    if any(row != budgets[0] for row in budgets[1:]):
        raise ArithmeticError("shard budgets differ")
    if [row["candidate_id"] for row in results] != expected_ids:
        raise ArithmeticError("shards do not form the exact frozen candidate partition")

    positives = []
    failed_closed = []
    chart_counts = Counter()
    raw_model_bits = []
    raw_quartic_bits = []
    per_parent = defaultdict(lambda: Counter(candidates=0, positives=0, failed_closed=0))
    compact = []
    for result in results:
        candidate = candidate_by_id[result["candidate_id"]]
        raw_model_bits.append(
            max(rational_bits(value) for value in candidate["raw_short_model"])
        )
        for key in ("parent_id", "parameter", "q_isomorphism_class_id", "nagao"):
            if result.get(key) != candidate.get(key):
                raise ArithmeticError(f"candidate/result mismatch for {key}")
        quotient = result.get("exact_quotient_rank_recovered")
        parent_counts = per_parent[result["parent_id"]]
        parent_counts["candidates"] += 1
        if quotient is None:
            failed_closed.append(result["candidate_id"])
            parent_counts["failed_closed"] += 1
        elif int(quotient) >= 1:
            positives.append(result["candidate_id"])
            parent_counts["positives"] += 1
        for cover in result.get("cover_records", []):
            chart_counts[chart_status(cover)] += 1
            raw_quartic_bits.append(
                int(cover["search"]["raw_rational_coefficient_maximum_bits"])
            )
        compact.append(compact_result(result, candidate))

    selmer_calls = len(positives)
    local_shard_witnesses = []
    for path, shard_results in shard_result_groups:
        normalized = [
            strip_runtime(compact_result(row, candidate_by_id[row["candidate_id"]]))
            for row in shard_results
        ]
        local_shard_witnesses.append(
            {
                "path": relative(path),
                "candidate_count": len(shard_results),
                "first_candidate_id": shard_results[0]["candidate_id"],
                "last_candidate_id": shard_results[-1]["candidate_id"],
                "normalized_result_sha256": canonical_digest(normalized),
            }
        )

    return {
        "schema": "elliptic-curves.icarm-mw16-nagao-finalist-half-lattice-merged.v1",
        "status": "PASS_COMPLETE_FROZEN_NAGAO_FINALIST_HALF_LATTICE_GATE",
        "source_finalist_count": len(candidates),
        "completed_candidate_count": len(results),
        "positive_candidate_ids": positives,
        "positive_candidate_count": len(positives),
        "failed_closed_candidate_ids": failed_closed,
        "failed_closed_candidate_count": len(failed_closed),
        "chart_status_counts": dict(sorted(chart_counts.items())),
        "arithmetic_size_diagnostic": {
            "raw_short_model_maximum_bits_range": [
                min(raw_model_bits),
                max(raw_model_bits),
            ],
            "raw_quartic_coefficient_maximum_bits_range": [
                min(raw_quartic_bits),
                max(raw_quartic_bits),
            ],
            "all_chart_attempts_timed_out": (
                chart_counts == Counter({"bounded_search_timeout": len(raw_quartic_bits)})
            ),
            "next_engineering_gate": (
                "exact Q-isomorphic model and section transport that materially "
                "reduces arithmetic size before any larger half-lattice budget"
            ),
        },
        "declared_budget": budgets[0],
        "candidate_ledger_sha256": candidate_ledger_digest,
        "source_ledger_sha256_at_merge": digest(args.source),
        "shards_share_one_source_ledger": len(shard_source_hashes) == 1,
        "source_candidate_content_matches_shards": (
            candidate_ledger_digest == EXPECTED_CANDIDATE_LEDGER_SHA256
        ),
        "local_shard_witnesses": local_shard_witnesses,
        "per_parent": {
            key: dict(sorted(value.items())) for key, value in sorted(per_parent.items())
        },
        "results": compact,
        "selmer_gate": {
            "positive_definition": "exact_quotient_rank_recovered >= 1",
            "authorized_candidate_ids": positives,
            "complete_residual_2_selmer_calls_required": selmer_calls,
            "expensive_continuation_authorized_before_complete_selmer": False,
        },
        "inputs": {
            relative(args.source): digest(args.source),
            relative(Path(__file__)): digest(Path(__file__)),
        },
        "claim_boundary": [
            "Nagao fixed candidate order only and contributes no rank evidence.",
            "An exact positive means a rational point survived exact group-law and finite-reduction independence checks.",
            "A zero is only the rank of the group recovered inside the declared bounded attempts.",
            "A bounded_search_timeout chart did not exhaust its height box and is censored, not a negative search result.",
            "No adaptive quotient wave, unrestricted point search, or Selmer computation occurs in this artifact.",
            "Only exact positives are passed to the complete same-minimal-curve residual 2-Selmer gate.",
        ],
        "reproducing_command": (
            "python3 elliptic-curves/cas/merge_icarm_mw16_nagao_finalist_half_lattice_shards.py"
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
            raise SystemExit("FAIL: compact merged certificate changed")
        print("PASS: compact merged certificate is unchanged")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    print(
        "MW16FINALISTMERGE|"
        f"completed={payload['completed_candidate_count']}|"
        f"positive={payload['positive_candidate_count']}|"
        f"failed_closed={payload['failed_closed_candidate_count']}|"
        f"output={relative(args.output)}|status={payload['status']}"
    )


if __name__ == "__main__":
    main()
