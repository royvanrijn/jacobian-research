#!/usr/bin/env python3
"""Audit whether foundry routes have enough data for inverse-ADE planning."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
FOUNDRY = GENERATED / "elkies-k3-lattice-foundry-v1.json"
SOURCE_ATTEMPTS = ROOT / "elkies-k3/data/lattice-foundry/source-equation-attempts-v1.json"
BENCHMARK = GENERATED / "elkies-k3-inverse-ade-target-planner-benchmark-v1.json"
OUTPUT = GENERATED / "elkies-k3-inverse-ade-foundry-readiness-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    foundry = json.loads(FOUNDRY.read_text())
    attempts = json.loads(SOURCE_ATTEMPTS.read_text())
    benchmark = json.loads(BENCHMARK.read_text())
    companions = {row["frame_id"]: row for row in foundry["companion_fibrations"]}
    targets = {row["frame_id"]: row for row in foundry["rootless_targets"]}
    frames = {
        row["frame_id"]: row
        for ns_row in foundry["ns_classes"]
        for row in ns_row["frames"]
    }

    route_rows = []
    missing_counts: Counter[str] = Counter()
    for route in foundry["route_ledger"]:
        source = companions[route["source_frame_id"]]
        target = targets[route["target_frame_id"]]
        missing = [
            "common_marked_NS_basis_and_source_U",
            "rank15_core_plus_rank2_bridge_decomposition",
            "graph_glue_generator",
            "good_neighbor_prime_plan",
            "prescribed_survival_birth_templates",
            "elliptic_neighbor_transport_or_relative_U_lift",
        ]
        for item in missing:
            missing_counts[item] += 1
        route_rows.append(
            {
                "route_id": route["route_id"],
                "ns_id": route["ns_id"],
                "source_frame_id": source["frame_id"],
                "source_root_rank": source["root_rank"],
                "source_root_type": source["root_type"],
                "target_frame_id": target["frame_id"],
                "source_and_target_positive_gram_present": bool(
                    frames[source["frame_id"]].get("gram")
                    and frames[target["frame_id"]].get("gram")
                ),
                "planner_ready": False,
                "missing_inputs": missing,
            }
        )

    exact_equation_sources = [
        row["source_id"] for row in attempts["attempts"] if row["equation_success"]
    ]
    benchmark_gate = {
        row["corridor"]: row["summary"] for row in benchmark["corridors"]
    }
    payload = {
        "schema": "elkies-k3.inverse-ade-foundry-readiness.v1",
        "status": "BLOCKED_FAIL_CLOSED_NO_FOUNDRY_PLANNER_RUN",
        "inputs": {
            relative(FOUNDRY): digest(FOUNDRY),
            relative(SOURCE_ATTEMPTS): digest(SOURCE_ATTEMPTS),
            relative(BENCHMARK): digest(BENCHMARK),
        },
        "benchmark_gate": benchmark_gate,
        "accounting": {
            "rootless_target_frames": len(targets),
            "target_ns_classes": len({row["ns_id"] for row in targets.values()}),
            "source_target_route_pairs": len(route_rows),
            "routes_with_positive_source_and_target_grams": sum(
                row["source_and_target_positive_gram_present"] for row in route_rows
            ),
            "routes_with_complete_inverse_ade_planner_inputs": sum(
                row["planner_ready"] for row in route_rows
            ),
            "equation_shortlist_entries": len(foundry["equation_shortlist"]),
            "curated_equation_source_attempts": len(attempts["attempts"]),
            "certified_characteristic_zero_equation_sources": len(exact_equation_sources),
            "source_root_rank_histogram": dict(
                sorted(Counter(row["source_root_rank"] for row in route_rows).items())
            ),
            "missing_input_counts": dict(sorted(missing_counts.items())),
        },
        "exact_equation_source_ids": exact_equation_sources,
        "routes": route_rows,
        "decision": (
            "Do not launch 936 unconstrained searches. The withheld benchmark does "
            "not pass its all-corridor recovery/10x gate, and no foundry route stores "
            "the bridge/glue plus marked survival-birth templates required by this "
            "planner. Moreover none of the six curated source attempts is yet a "
            "characteristic-zero equation source."
        ),
        "proof_boundary": (
            "This is an exact data-readiness audit, not a nonexistence theorem. "
            "Positive frame Gram matrices alone do not supply a same-surface marked "
            "elliptic-neighbour route, an equation, or arithmetic descent."
        ),
        "reproduce": "python3 elkies-k3/scripts/audit_inverse_ade_foundry_readiness.py",
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print(payload["status"])
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
