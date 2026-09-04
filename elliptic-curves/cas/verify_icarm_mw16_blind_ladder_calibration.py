#!/usr/bin/env python3
"""Compare complement-blind MW16 ladder returns with the held-out A1 atlas."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-icarm-11952-norm8-a1-mw16-atlas-v1.json"
INITIAL = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_ladder_blind_v1.json"
PRESENTATIONS = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_parent_presentation_audit_v1.json"
CURVE398 = ROOT / "artifacts/generated-results/elliptic-curves/curve398_mw16_adaptive_half_lattice_blind_v1.json"
CURVE398_VERIFY = ROOT / "artifacts/generated-results/elliptic-curves/curve398_mw16_adaptive_half_lattice_verification_v1.json"
CURVE400 = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_curve400_adaptive_calibration_v1.json"
PROSPECTIVE = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_half_lattice_h300_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_mw16_blind_ladder_calibration_v1.json"
GENERIC_RANK = 16
CURVES = (398, 400, 401, 542, 548)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def build():
    atlas = json.loads(ATLAS.read_text())
    initial = json.loads(INITIAL.read_text())
    presentations = json.loads(PRESENTATIONS.read_text())
    curve398 = json.loads(CURVE398.read_text())
    curve398_verify = json.loads(CURVE398_VERIFY.read_text())
    curve400 = json.loads(CURVE400.read_text())
    prospective = json.loads(PROSPECTIVE.read_text())
    if atlas.get("status") != "PASS_EXACT_COMPLETE_PRIORITY_ICARM_A1_MW16_ATLAS":
        raise ArithmeticError("complete A1 atlas is not passing")
    if initial.get("status") != "PASS_COMPLETE_NINE_PARENT_INITIAL_HALF_LATTICE_LADDER":
        raise ArithmeticError("nine-presentation blind initial ladder is not passing")
    if presentations.get("status") != "PASS_EXACT_NINE_PRESENTATIONS_FIVE_FIBRATIONS":
        raise ArithmeticError("presentation audit is not passing")
    if presentations.get("exact_fibration_class_count") != 5:
        raise ArithmeticError("presentation audit no longer gives five fibrations")
    if curve398.get("status") != "STOPPED_AT_DECLARED_LIFT_LIMIT":
        raise ArithmeticError("curve-398 adaptive calibration status changed")
    if curve398_verify.get("status") != "PASS_EXACT_CROSS_FIBRATION_RANK30_REDISCOVERY":
        raise ArithmeticError("curve-398 held-out verification is not passing")
    if curve400.get("status") != "PASS_COMPLETE_CURVE400_ADAPTIVE_CALIBRATION":
        raise ArithmeticError("curve-400 adaptive calibration is not passing")
    if (
        prospective.get("status")
        != "PASS_COMPLETE_FROZEN_NAGAO_FINALIST_HALF_LATTICE_GATE"
        or prospective.get("completed_candidate_count") != 104
        or prospective.get("positive_candidate_count") != 0
        or prospective.get("failed_closed_candidate_count") != 0
        or prospective.get("chart_status_counts")
        != {"bounded_search_timeout": 856}
    ):
        raise ArithmeticError("prospective MW16 half-lattice gate changed")

    atlas_hits = {
        int(row["curve_id"]): row
        for row in atlas["targets"]
        if row["outcome"] == "HIT" and int(row["curve_id"]) in CURVES
    }
    if set(atlas_hits) != set(CURVES):
        raise ArithmeticError("held-out atlas hit set changed")
    parents_by_curve = {curve_id: [] for curve_id in CURVES}
    for row in initial["parents"]:
        curve_id = int(row["curve_id"])
        if curve_id in parents_by_curve:
            parents_by_curve[curve_id].append(row)
    if [len(parents_by_curve[c]) for c in CURVES] != [2, 2, 1, 1, 3]:
        raise ArithmeticError("nine-presentation nesting changed")

    curve398_total = len(curve398["current_basis"]) - GENERIC_RANK
    if (
        curve398_total != 14
        or curve398_verify["blind_transition"][
            "after_complete_five_bit_adaptive_wave"
        ]
        != 30
    ):
        raise ArithmeticError("curve-398 adaptive total changed")
    curve400_total = int(curve400["exact_quotient_rank_recovered_total"])
    if curve400_total != 12:
        raise ArithmeticError("curve-400 adaptive total changed")
    curve400_initial = curve400["initial"]
    curve400_adaptive = curve400["adaptive"]
    if (
        curve400_initial["basis_rank_after"] != 21
        or curve400_initial["exact_quotient_rank_recovered"] != 5
        or len(curve400_initial["cover_records"]) != 4
        or Counter(
            row["search"]["status"] for row in curve400_initial["cover_records"]
        )
        != Counter({"bounded_search_complete": 4})
        or curve400_adaptive["basis_rank_before"] != 21
        or curve400_adaptive["basis_rank_after"] != 28
        or curve400_adaptive["incremental_exact_quotient_rank_recovered"] != 7
        or len(curve400_adaptive["cover_records"]) != 124
        or Counter(
            row["search"]["status"] for row in curve400_adaptive["cover_records"]
        )
        != Counter({"bounded_search_complete": 124})
        or Counter(
            row["type"]
            for row in curve400_adaptive["discovered_group_saturation"]["events"]
        )
        != Counter({"NEW_Q_INDEPENDENT_DIRECTION": 7})
    ):
        raise ArithmeticError("curve-400 complete adaptive certificate changed")

    rows = []
    for curve_id in CURVES:
        parents = parents_by_curve[curve_id]
        initial_values = {
            int(row["exact_quotient_rank_recovered"]) for row in parents
        }
        if len(initial_values) != 1:
            raise ArithmeticError("pseudoreplicate initial responses disagree")
        initial_recovered = initial_values.pop()
        demonstrated = int(atlas_hits[curve_id]["snapshot_rank_lower_bound"]) - GENERIC_RANK
        adaptive_recovered = {
            398: curve398_total,
            400: curve400_total,
        }.get(curve_id, initial_recovered)
        adaptive_status = {
            398: "complete_five_bit_wave_then_declared_stop",
            400: "complete_five_bit_wave",
            401: "not_run_initial_ten_bit_quotient",
            542: "not_needed_initial_full_recovery",
            548: "not_needed_initial_full_recovery",
        }[curve_id]
        rows.append(
            {
                "curve_id": curve_id,
                "presentation_count": len(parents),
                "parent_ids": [row["parent_id"] for row in parents],
                "atlas_rank_lower_bound": int(atlas_hits[curve_id]["snapshot_rank_lower_bound"]),
                "atlas_demonstrated_jump_over_mw16": demonstrated,
                "initial_exact_quotient_rank_recovered": initial_recovered,
                "initial_shortfall_from_demonstrated_jump": demonstrated - initial_recovered,
                "best_blind_ladder_exact_quotient_rank_recovered": adaptive_recovered,
                "remaining_shortfall_from_demonstrated_jump": demonstrated - adaptive_recovered,
                "adaptive_status": adaptive_status,
            }
        )

    initial_full = sum(row["initial_shortfall_from_demonstrated_jump"] == 0 for row in rows)
    initial_within_one = sum(row["initial_shortfall_from_demonstrated_jump"] <= 1 for row in rows)
    demonstrated_total = sum(row["atlas_demonstrated_jump_over_mw16"] for row in rows)
    initial_total = sum(row["initial_exact_quotient_rank_recovered"] for row in rows)
    best_total = sum(row["best_blind_ladder_exact_quotient_rank_recovered"] for row in rows)
    if (initial_full, initial_within_one, demonstrated_total, initial_total, best_total) != (
        2,
        3,
        55,
        38,
        54,
    ):
        raise ArithmeticError("curve-level ladder summary changed")

    return {
        "schema": "elliptic-curves.icarm-mw16-blind-ladder-calibration.v1",
        "status": "PASS_EXACT_FIVE_CURVE_BLIND_LADDER_CALIBRATION",
        "observation_unit": {
            "primary": "target_curve",
            "curve_count": 5,
            "presentation_count": 9,
            "exact_fibration_class_count": 5,
            "pseudoreplicate_rule": (
                "repeated coordinate presentations are never counted as independent outcomes"
            ),
        },
        "curve_results": rows,
        "curve_level_summary": {
            "demonstrated_jump_directions_total": demonstrated_total,
            "initial_exact_directions_recovered_total": initial_total,
            "best_blind_ladder_exact_directions_recovered_total": best_total,
            "initial_full_recovery_curve_count": initial_full,
            "initial_within_one_direction_curve_count": initial_within_one,
            "best_blind_ladder_full_recovery_curve_count": 4,
            "only_remaining_demonstrated_direction": {
                "curve_id": 401,
                "count": 1,
                "complete_next_wave_chart_count": 8 * ((1 << 10) - 1),
                "complete_next_wave_run": False,
            },
        },
        "prospective_gate": {
            "nagao_finalist_count": prospective["source_finalist_count"],
            "exact_specialized_fibre_count": prospective[
                "completed_candidate_count"
            ],
            "exact_positive_candidate_count": prospective[
                "positive_candidate_count"
            ],
            "fail_closed_candidate_count": prospective[
                "failed_closed_candidate_count"
            ],
            "chart_status_counts": prospective["chart_status_counts"],
            "wholly_timeout_censored": prospective[
                "arithmetic_size_diagnostic"
            ]["all_chart_attempts_timed_out"],
            "residual_selmer_calls_authorized": prospective["selmer_gate"][
                "complete_residual_2_selmer_calls_required"
            ],
            "next_engineering_gate": prospective[
                "arithmetic_size_diagnostic"
            ]["next_engineering_gate"],
        },
        "inputs": {
            relative(path): digest(path)
            for path in (
                ATLAS,
                INITIAL,
                PRESENTATIONS,
                CURVE398,
                CURVE398_VERIFY,
                CURVE400,
                PROSPECTIVE,
                Path(__file__),
            )
        },
        "claim_boundary": [
            "Atlas rank lower bounds and jump labels are loaded only here, after the blind searches completed.",
            "The five target curves, not nine coordinate presentations, are the statistical observations.",
            "This purposive five-curve calibration supports detector engineering, not a population success-rate estimate.",
            "Matching a demonstrated jump lower bound does not prove the target curve has no further rational directions.",
            "The unrun 8184-chart curve-401 adaptive wave is recorded as unrun, not inferred from the other curves.",
            "The prospective zero is wholly timeout-censored and supplies no negative rank or Nagao-ordering evidence.",
            "No prospective candidate reaches residual Selmer or expensive continuation.",
        ],
        "reproducing_command": (
            "python3 elliptic-curves/cas/verify_icarm_mw16_blind_ladder_calibration.py --check"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != encoded:
            raise SystemExit("FAIL: MW16 blind-ladder calibration changed")
        print("PASS: MW16 blind-ladder calibration is unchanged")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    print(
        "MW16BLINDCAL|curves=5|presentations=9|"
        f"initial={payload['curve_level_summary']['initial_exact_directions_recovered_total']}/55|"
        f"best={payload['curve_level_summary']['best_blind_ladder_exact_directions_recovered_total']}/55|"
        f"output={relative(args.output)}|status={payload['status']}"
    )


if __name__ == "__main__":
    main()
