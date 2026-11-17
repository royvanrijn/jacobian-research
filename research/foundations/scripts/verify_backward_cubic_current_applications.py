#!/usr/bin/env python3
"""Apply the backward cubic ledger to active restricted-minima archives."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
ARTIFACTS = ROOT / "artifacts" / "generated-results"
OUTPUT = ARTIFACTS / "backward_cubic_current_applications.json"

from jcsearch.backward_cubic import BackwardTerminalProfile  # noqa: E402
from search_restricted_bcw_circuits import (  # noqa: E402
    replay_encoded_plan,
    terminal_profile,
)


ARCHIVES = (
    "restricted_bcw_circuit_search_xxs_w32_d28.json",
    "restricted_bcw_circuit_search_yvyb_structural_w22.json",
    "restricted_bcw_circuit_search_ayb_yvyb_w36.json",
    "restricted_bcw_circuit_search_xvvz_v2h_structural_w25.json",
    "restricted_bcw_circuit_search_xvvz_v2h_mixed_w25.json",
    "restricted_bcw_circuit_search_aspert_m12.json",
    "restricted_bcw_circuit_search_xxs_rank_hybrid_w24.json",
    "restricted_bcw_circuit_search_xxs_v2r_w16.json",
    "restricted_bcw_circuit_search_xxs_y2vb_w16.json",
    "restricted_bcw_circuit_search_all_w64.json",
)


def retained_records() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for name in ARCHIVES:
        artifact = json.loads((ARTIFACTS / name).read_text())
        assert artifact["format"] == "restricted-bcw-circuit-search-v2"
        for index, terminal in enumerate(artifact["pareto_terminals"]):
            stored = terminal["profile"]
            rank = int(stored["cubic_output_rank"])
            homogeneous_dimension = int(stored["homogeneous_dimension"])
            base_dimension = homogeneous_dimension - rank - 1
            profile = BackwardTerminalProfile(base_dimension, rank)
            assert profile.homogeneous_dimension == homogeneous_dimension
            records.append(
                {
                    "archive": name,
                    "pareto_index": index,
                    "base_dimension": base_dimension,
                    "cubic_output_rank": rank,
                    "raw_homogeneous_dimension": homogeneous_dimension,
                    "constant_kernel_quotient_dimension": int(
                        stored["quotient_dimension"]
                    ),
                    "direct_cubic_key": list(profile.direct_cubic_key),
                    "homogeneous_key": list(profile.homogeneous_key),
                    "plan": terminal["plan"],
                }
            )
    return records


def replay_record(record: dict[str, object]) -> dict[str, object]:
    state = replay_encoded_plan(record["plan"])
    exact = terminal_profile(state, hessian_power=False)
    assert exact["direct_cubic_key"] == record["direct_cubic_key"]
    assert exact["homogeneous_key"] == record["homogeneous_key"]
    assert exact["quotient_dimension"] == record[
        "constant_kernel_quotient_dimension"
    ]
    source_pairs = exact["source_collision_pairs"]
    projected_pairs = exact["projected_collision_pairs"]
    assert source_pairs
    assert projected_pairs
    return {
        "archive": record["archive"],
        "pareto_index": record["pareto_index"],
        "direct_cubic_key": exact["direct_cubic_key"],
        "homogeneous_key": exact["homogeneous_key"],
        "constant_kernel_quotient_dimension": exact["quotient_dimension"],
        "source_collision_pairs": source_pairs,
        "projected_collision_pairs": projected_pairs,
    }


def main() -> None:
    records = retained_records()
    assert records
    best_direct = min(records, key=lambda row: tuple(row["direct_cubic_key"]))
    best_raw_homogeneous = min(
        records, key=lambda row: tuple(row["homogeneous_key"])
    )
    best_quotient = min(
        records,
        key=lambda row: (
            row["constant_kernel_quotient_dimension"],
            tuple(row["homogeneous_key"]),
        ),
    )

    assert best_direct["direct_cubic_key"] == [18, 26, 7]
    assert best_raw_homogeneous["homogeneous_key"] == [26, 18, 7]
    assert best_quotient["constant_kernel_quotient_dimension"] == 22

    # Replay representatives of the best direct source and the best final
    # homogeneous quotient through the newly integrated pair-aware pipeline.
    replayed = [
        replay_record(best_direct),
        replay_record(best_quotient),
    ]

    payload = {
        "format": "backward-cubic-current-applications-v1",
        "status": (
            "exact audit of retained bounded-search representatives; "
            "not a lower bound and no new counterexample"
        ),
        "current_fields": {
            "arbitrary_degree_three_dimension": {
                "objective": "(base dimension, homogeneous dimension, cubic rank)",
                "best_retained_key": best_direct["direct_cubic_key"],
                "incumbent_macfarlane_key": [13, 20, 6],
            },
            "cubic_homogeneous_dimension": {
                "objective": "(raw homogeneous dimension, base dimension, cubic rank)",
                "best_retained_key": best_raw_homogeneous["homogeneous_key"],
                "incumbent_macfarlane_dimension": 20,
                "best_retained_constant_kernel_quotient_dimension": (
                    best_quotient["constant_kernel_quotient_dimension"]
                ),
            },
            "homogeneous_quartic_hessian_nilpotent_dimension": {
                "upper_bound": 40,
                "construction": (
                    "standard homogeneous cotangent lift of the pinned "
                    "MacFarlane G20 certificate"
                ),
                "internal_independent_replay_upper_bound": 42,
            },
            "restricted_rank_minima": {
                "policy": (
                    "retain its existing rank/index/Hessian Pareto objective "
                    "and record the backward keys as independent metadata"
                ),
                "collision_policy": (
                    "at least one exact collision pair must remain distinct"
                ),
                "homogenizing_level_policy": (
                    "normalize to t=1; all nonzero parent fibers scale to "
                    "the base map and the t=0 fiber is injective"
                ),
            },
        },
        "archive_scope": {
            "files": list(ARCHIVES),
            "retained_representative_count": len(records),
            "warning": (
                "the historical archives were selected by restricted-rank "
                "objectives, so their absence of a smaller direct key cannot "
                "be interpreted as a failed direct-dimension search"
            ),
        },
        "exact_replays": replayed,
        "conclusion": (
            "none of the retained restricted-minima representatives improves "
            "MacFarlane in either backward objective; future searches must "
            "maintain a separate direct-dimension archive"
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    digest = hashlib.sha256(OUTPUT.read_bytes()).hexdigest()
    print(
        "PASS backward applications: audited "
        f"{len(records)} retained restricted-minima representatives"
    )
    print(
        "PASS backward applications: best retained direct key is "
        f"{tuple(best_direct['direct_cubic_key'])}"
    )
    print(
        "PASS backward applications: exact replays preserve at least one "
        "collision pair"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")
    print(f"SHA256 {digest}")


if __name__ == "__main__":
    main()
