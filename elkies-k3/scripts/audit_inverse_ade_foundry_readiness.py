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
BENCHMARK = GENERATED / "elkies-k3-inverse-ade-target-planner-benchmark-v2.json"
CURATED_ROUTE = (
    ROOT / "elkies-k3/data/lattice-foundry/planner-ready-h3-a1-r17-v1.json"
)
CURATED_CERTIFICATE = (
    GENERATED / "elkies-k3-single-planner-ready-foundry-route-v1.json"
)
OUTPUT = GENERATED / "elkies-k3-inverse-ade-foundry-readiness-v2.json"

CURATED_REQUIRED_FIELDS = [
    "common_marked_NS_basis_and_source_U",
    "rank15_core_plus_rank2_bridge_decomposition",
    "graph_glue_generator",
    "good_neighbor_prime_plan",
    "prescribed_survival_birth_templates",
    "elliptic_neighbor_transport_or_relative_U_lift",
]


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
    curated_manifest = json.loads(CURATED_ROUTE.read_text())
    curated_certificate = json.loads(CURATED_CERTIFICATE.read_text())
    assert curated_manifest["route_count"] == 1
    assert len(curated_manifest["routes"]) == 1
    curated = curated_manifest["routes"][0]
    curated_missing = [
        field for field in CURATED_REQUIRED_FIELDS if field not in curated
    ]
    curated_ready = (
        not curated_missing
        and curated_certificate["route_count"] == 1
        and curated_certificate["route_id"] == curated["route_id"]
        and curated_certificate["status"]
        == "PASS_ONE_BLIND_PLANNER_READY_ROUTE_TO_MW17"
    )
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
            "marked_target_core_in_source_rational_space",
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
    if curated_ready:
        exact_equation_sources.append("H3-fixed-final-A1-reverse-QQ")
    benchmark_gate = {
        row["corridor"]: row["summary"] for row in benchmark["corridors"]
    }
    payload = {
        "schema": "elkies-k3.inverse-ade-foundry-readiness.v2",
        "status": (
            "PASS_EXACTLY_ONE_CURATED_PLANNER_READY_ROUTE"
            if curated_ready
            else "BLOCKED_CURATED_ROUTE_NOT_CERTIFIED"
        ),
        "inputs": {
            relative(FOUNDRY): digest(FOUNDRY),
            relative(SOURCE_ATTEMPTS): digest(SOURCE_ATTEMPTS),
            relative(BENCHMARK): digest(BENCHMARK),
            relative(CURATED_ROUTE): digest(CURATED_ROUTE),
            relative(CURATED_CERTIFICATE): digest(CURATED_CERTIFICATE),
        },
        "benchmark_gate": benchmark_gate,
        "accounting": {
            "rootless_target_frames": len(targets),
            "target_ns_classes": len({row["ns_id"] for row in targets.values()}),
            "source_target_route_pairs": len(route_rows),
            "bulk_source_target_route_pairs": len(route_rows),
            "all_route_records_after_curated_addition": len(route_rows) + 1,
            "routes_with_positive_source_and_target_grams": sum(
                row["source_and_target_positive_gram_present"] for row in route_rows
            ),
            "routes_with_complete_inverse_ade_planner_inputs": sum(
                row["planner_ready"] for row in route_rows
            )
            + int(curated_ready),
            "bulk_routes_with_complete_inverse_ade_planner_inputs": sum(
                row["planner_ready"] for row in route_rows
            ),
            "curated_planner_ready_routes": int(curated_ready),
            "equation_shortlist_entries": len(foundry["equation_shortlist"]),
            "curated_equation_source_attempts": len(attempts["attempts"]),
            "certified_characteristic_zero_equation_sources": len(exact_equation_sources),
            "source_root_rank_histogram": dict(
                sorted(Counter(row["source_root_rank"] for row in route_rows).items())
            ),
            "missing_input_counts": dict(sorted(missing_counts.items())),
        },
        "exact_equation_source_ids": exact_equation_sources,
        "curated_routes": [
            {
                "route_id": curated["route_id"],
                "source_mw_rank_if_rho_19": curated["selection"][
                    "source_mw_rank_if_rho_19"
                ],
                "target_frame_id": curated["selection"]["target_frame_id"],
                "target_ade_type": curated["selection"]["target_ade_type"],
                "planner_ready": curated_ready,
                "missing_inputs": curated_missing,
                "bulk_route_ledger_member": False,
            }
        ],
        "routes": route_rows,
        "decision": (
            "Run the one certified H3 A1/MW16 to NS0001-F001 positive control; do "
            "not launch the 936 unconstrained bulk searches. The 936 original rows "
            "remain unchanged and unready. The single added route has an exact QQ "
            "equation source and all six manually curated planner/transport fields."
        ),
        "proof_boundary": (
            "This is an exact data-readiness audit, not a nonexistence theorem. "
            "The one positive control is separate from the 936 original route-ledger "
            "rows. Positive frame Gram matrices alone do not supply a same-surface "
            "marked elliptic-neighbour route, an equation, or arithmetic descent."
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
