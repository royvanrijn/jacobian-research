#!/usr/bin/env python3
"""Audit and promote the curve-385 primary sparse rank-32 search ledger.

The search ledger is deliberately verbose.  This script checks its frozen
protocol, stage accounting, bounded-search outcomes, and exact discovered-
group classifications before writing a compact manifest and a deterministic
gzip copy of the complete ledger.  The result is a bounded negative search
outcome, never a rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
LOCAL_LEDGER = ROOT / "artifacts/local/curve385_sparse_quotient_rank32_search_v1.json"
BLIND = ART / "curve385_iterated_half_lattice_blind_v1.json"
PROTOCOL = ART / "curve385_sparse_quotient_rank32_protocol_v1.json"
OUTPUT = ART / "curve385_sparse_quotient_rank32_primary_v1.json"
ARCHIVE = ART / "curve385_sparse_quotient_rank32_primary_ledger_v1.json.gz"

EXPECTED_LEDGER_SHA256 = "17600ab552c8c4c5184d8ec02c6743c475424998e2672a8eeacc3ee75df5b77d"
EXPECTED_BLIND_SHA256 = "356001898f738f607d984e081663a015825e11de0c606d35055af156eb2d7502"
EXPECTED_PROTOCOL_SHA256 = "2c9150f50f305b8aa3763590cd5e81c4d7e121f9373177827780789ce472834f"
EXPECTED_PROTOCOL_DEFINITION_SHA256 = "5723679da2907e036095f90376cdabde457a4f7ba5bc284ad4a4ca3edea1aa37"
EXPECTED_BASIS_SHA256 = "5a6d5ff70b409de413722e35b71e430025852639a71a166f2b7cfc288f27a230"

EXPECTED_CONFIGURATION = {
    "checkpoint_every_completed_searches": 10,
    "gp_stack_bytes_each_quartic": 1_000_000_000,
    "height_bound_each_quartic": 100_000,
    "maximum_lattice_states": 4,
    "maximum_stage_each_lattice_state": 2,
    "relation_chunk_size": 64,
    "relation_timeout_seconds_each_chunk": 180.0,
    "retries": 0,
    "wall_timeout_seconds_each_quartic": 15.0,
}


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def point_key(record: dict[str, str]) -> tuple[Fraction, Fraction]:
    point = Fraction(record["x"]), Fraction(record["y"])
    opposite = point[0], -point[1]
    return min(point, opposite)


def discovery_map(rows: list[dict[str, Any]]) -> dict[tuple[Fraction, Fraction], set[str]]:
    answer: dict[tuple[Fraction, Fraction], set[str]] = {}
    for row in rows:
        key = point_key(row["point"])
        if key in answer:
            raise ArithmeticError("duplicate canonical point in discovery ledger")
        answer[key] = set(map(str, row["sources"]))
    return answer


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def validate_and_summarize(raw: bytes) -> dict[str, Any]:
    require(sha256(raw).hexdigest() == EXPECTED_LEDGER_SHA256, "full ledger hash changed")
    require(file_sha256(BLIND) == EXPECTED_BLIND_SHA256, "blind M29 input changed")
    require(file_sha256(PROTOCOL) == EXPECTED_PROTOCOL_SHA256, "sparse protocol changed")
    ledger = json.loads(raw)
    blind = json.loads(BLIND.read_text())
    protocol = json.loads(PROTOCOL.read_text())

    require(
        ledger.get("schema") == "elliptic-curves.curve385-sparse-quotient-rank32-search.v1",
        "wrong search-ledger schema",
    )
    require(
        ledger.get("status") == "STOPPED_AFTER_DECLARED_SPARSE_STAGE_LIMIT",
        "primary sparse campaign did not reach its declared terminal state",
    )
    require(ledger.get("configuration") == EXPECTED_CONFIGURATION, "search budget changed")
    require(
        ledger.get("stop")
        == {"basis_rank": 29, "maximum_stage": 2, "no_rank_upper_bound_claimed": True},
        "search stop boundary changed",
    )
    require(
        ledger.get("protocol")
        == {
            "definition_sha256": EXPECTED_PROTOCOL_DEFINITION_SHA256,
            "path": relative(PROTOCOL),
            "whole_file_sha256": EXPECTED_PROTOCOL_SHA256,
        },
        "search ledger does not identify the frozen protocol",
    )
    require(
        protocol.get("protocol_definition_hash") == EXPECTED_PROTOCOL_DEFINITION_SHA256,
        "protocol definition changed",
    )
    require(
        ledger.get("source_blind_ledger_sha256") == EXPECTED_BLIND_SHA256,
        "search ledger identifies the wrong M29 source",
    )
    require(ledger.get("curve") == blind.get("curve"), "curve data changed after M29")
    require(ledger.get("old_deep43") == blind.get("old_deep43"), "deep-43 classes changed")
    require(len(ledger["old_deep43"]["masks"]) == 43, "deep-class count is not 43")

    states = ledger.get("lattice_states", [])
    require(len(states) == 1, "primary no-growth run must contain one lattice state")
    state = states[0]
    require(
        (state.get("index"), state.get("basis_rank"), state.get("quotient_bit_count"))
        == (1, 29, 12),
        "unexpected primary lattice state",
    )
    require(
        state.get("status") == "DECLARED_STAGE_LIMIT_REACHED_WITHOUT_GROWTH",
        "lattice state did not stop at the declared stage limit",
    )
    require(state.get("basis_sha256") == EXPECTED_BASIS_SHA256, "starting basis changed")
    require(canonical_hash(state["basis"]) == EXPECTED_BASIS_SHA256, "starting basis hash invalid")
    require(ledger.get("current_basis") == state["basis"], "final basis differs from start")
    require(len(ledger["current_basis"]) == 29, "final discovered rank is not 29")

    starting_discoveries = discovery_map(blind["discoveries"])
    final_discoveries = discovery_map(ledger["discoveries"])
    for point, sources in starting_discoveries.items():
        require(
            point in final_discoveries and sources <= final_discoveries[point],
            "a frozen M29 discovery or provenance tag was lost",
        )
    starting_keys = set(map(str, blind["searched_base_point_keys"]))
    final_keys = set(map(str, ledger["searched_base_point_keys"]))
    require(starting_keys <= final_keys, "a previously searched chart key was lost")
    require(len(final_keys) == len(ledger["searched_base_point_keys"]), "duplicate final chart key")

    stage_plans = protocol["stage_plans_by_quotient_bit_count"]["12"][:2]
    stages = state.get("stages", [])
    require(len(stages) == 2, "primary campaign must contain exactly two stages")
    old_masks = set(map(int, ledger["old_deep43"]["masks"]))
    all_new_base_keys: set[str] = set()
    stage_summaries = []
    total_occurrences = 0
    total_charts_with_points = 0
    for stage, plan in zip(stages, stage_plans, strict=True):
        stage_index = int(stage["index"])
        require(stage.get("id") == plan["id"], "stage id differs from protocol")
        require(stage.get("policy") == plan, "stage policy differs from protocol")
        require(stage.get("status") == "CLASSIFIED", "sparse stage is not classified")
        require(stage.get("group_changed") is False, "primary stage changed the group")
        require(stage.get("basis_rank_after") == 29, "primary stage changed the rank")
        require(stage.get("basis_after") == state["basis"], "primary stage changed the basis")
        require(
            stage.get("basis_after_sha256") == EXPECTED_BASIS_SHA256,
            "post-stage basis hash changed",
        )
        require(stage.get("new_independent_direction_count") == 0, "rank event present")
        require(stage.get("finite_index_saturation_event_count") == 0, "index event present")

        ranking = stage["ranking"]
        require(
            ranking.get("quotient_bit_count") == 12
            and ranking.get("physical_word_count") == len(plan["new_physical_words"])
            and ranking.get("ranked_lift_count") == plan["new_chart_count"],
            "ranked sparse-lift accounting changed",
        )
        require(
            ranking.get("priority_order_identical_between_scales") is True,
            "CVP audit and operative priority orders disagree",
        )

        records = stage["cover_records"]
        skipped = list(map(str, stage["unchanged_previously_searched_chart_keys"]))
        require(len(records) == stage["searched_new_chart_count"], "searched count mismatch")
        require(len(skipped) == stage["unchanged_previously_searched_chart_count"], "skip count mismatch")
        require(len(set(skipped)) == len(skipped), "duplicate skipped chart key")
        require(
            len(records) + len(skipped) == plan["new_chart_count"],
            "planned charts are not fully accounted for",
        )
        require(stage["bounded_complete_count"] == len(records), "bounded completion mismatch")
        require(stage["timeout_count"] == 0, "timeout present in completed primary stage")
        require(stage["pari_failure_count"] == 0, "PARI failure present in completed primary stage")

        words = set(map(int, plan["new_physical_words"]))
        pairs: set[tuple[int, int]] = set()
        priorities: set[int] = set()
        stage_base_keys: set[str] = set()
        finite_occurrences = 0
        charts_with_points = 0
        for record in records:
            pair = int(record["old_mask"]), int(record["physical_quotient_word"])
            require(pair[0] in old_masks and pair[1] in words, "cover lies outside stage policy")
            require(pair not in pairs, "duplicate searched quotient chart")
            pairs.add(pair)
            priority = int(record["priority"])
            require(1 <= priority <= plan["new_chart_count"], "invalid sparse priority")
            require(priority not in priorities, "duplicate sparse priority")
            priorities.add(priority)
            require(
                record["search"].get("status") == "bounded_search_complete",
                "non-complete search record in complete stage",
            )
            base_key = str(record["base_point_key"])
            require(base_key not in stage_base_keys, "duplicate searched base point in stage")
            require(base_key not in starting_keys, "previously searched base point was rerun")
            require(base_key not in all_new_base_keys, "base point was rerun across stages")
            stage_base_keys.add(base_key)
            all_new_base_keys.add(base_key)
            points = record["search"]["finite_curve_points"]
            finite_occurrences += len(points)
            charts_with_points += bool(points)
            source = (
                f"sparse:state:1:stage:{stage_index}:old:{pair[0]}:q:{pair[1]}"
            )
            for point_record in points:
                point = point_key(point_record)
                require(
                    point in final_discoveries and source in final_discoveries[point],
                    "returned rational point is missing from discovery provenance",
                )

        saturation = stage["discovered_group_saturation"]
        require(
            saturation.get("status") == "PASS_BASIS_EQUALS_DISCOVERED_GROUP",
            "discovered points were not exactly classified",
        )
        require(saturation.get("basis_rank") == 29, "classified basis rank changed")
        require(saturation.get("basis_sha256") == EXPECTED_BASIS_SHA256, "classification basis changed")
        require(saturation.get("events") == [], "unexpected discovered-group event")
        require(saturation.get("unresolved") == [], "unresolved discovery present")
        require(
            saturation.get("exact_integral_relation_count")
            == saturation.get("discovered_nonbasis_point_count")
            == len(saturation.get("exact_integral_relations", [])),
            "exact relation accounting is incomplete",
        )

        source_prefix = f"sparse:state:1:stage:{stage_index}:"
        stage_discovery_count = sum(
            any(source.startswith(source_prefix) for source in sources)
            for sources in final_discoveries.values()
        )
        total_occurrences += finite_occurrences
        total_charts_with_points += charts_with_points
        stage_summaries.append(
            {
                "index": stage_index,
                "id": stage["id"],
                "basis_id": plan["basis_id"],
                "basis_weight_shell": plan["basis_weight_shell"],
                "planned_chart_count": plan["new_chart_count"],
                "fresh_search_count": len(records),
                "previously_searched_exact_chart_count": len(skipped),
                "bounded_complete_count": stage["bounded_complete_count"],
                "timeout_count": stage["timeout_count"],
                "pari_failure_count": stage["pari_failure_count"],
                "charts_returning_finite_points": charts_with_points,
                "finite_point_occurrence_count": finite_occurrences,
                "distinct_discoveries_with_stage_source": stage_discovery_count,
                "rank_before": 29,
                "rank_after": 29,
                "new_independent_direction_count": 0,
                "finite_index_saturation_event_count": 0,
                "exact_integral_relation_count": saturation["exact_integral_relation_count"],
                "ranked_lifts_sha256": ranking["ranked_lifts_sha256"],
                "cover_records_sha256": canonical_hash(records),
                "skipped_chart_keys_sha256": canonical_hash(skipped),
                "exact_integral_relations_sha256": canonical_hash(
                    saturation["exact_integral_relations"]
                ),
                "wall_seconds": stage["wall_seconds"],
                "cpu_seconds_completed_searches": stage["cpu_seconds_completed_searches"],
            }
        )

    require(
        final_keys == starting_keys | all_new_base_keys,
        "final searched-key set is not the exact old-plus-fresh union",
    )
    require(
        len(all_new_base_keys) == sum(row["fresh_search_count"] for row in stage_summaries),
        "fresh chart keys are not unique",
    )
    primary_sources = {
        point
        for point, sources in final_discoveries.items()
        if any(source.startswith("sparse:state:1:stage:") for source in sources)
    }

    return {
        "schema": "elliptic-curves.curve385-sparse-quotient-rank32-primary-result.v1",
        "status": "PASS_COMPLETE_PRIMARY_SPARSE_CAMPAIGN_BOUNDED_NO_GROWTH",
        "program": {
            "path": relative(Path(__file__)),
            "sha256": file_sha256(Path(__file__)),
        },
        "full_ledger": {
            "uncompressed_sha256": EXPECTED_LEDGER_SHA256,
            "uncompressed_bytes": len(raw),
        },
        "frozen_inputs": {
            "blind_M29": {"path": relative(BLIND), "sha256": EXPECTED_BLIND_SHA256},
            "sparse_protocol": {
                "path": relative(PROTOCOL),
                "sha256": EXPECTED_PROTOCOL_SHA256,
                "definition_sha256": EXPECTED_PROTOCOL_DEFINITION_SHA256,
            },
        },
        "configuration": ledger["configuration"],
        "campaign": {
            "curve_id": "curve385-rank29",
            "starting_certified_rank_lower_bound": 29,
            "target_rank_lower_bound": 32,
            "final_certified_rank_lower_bound": 29,
            "quotient_bit_count": 12,
            "distinguished_generic_deep_class_count": 43,
            "basis_sha256": EXPECTED_BASIS_SHA256,
            "stage_summaries": stage_summaries,
            "total_planned_chart_count": sum(
                row["planned_chart_count"] for row in stage_summaries
            ),
            "total_fresh_search_count": len(all_new_base_keys),
            "total_previously_searched_exact_chart_count": sum(
                row["previously_searched_exact_chart_count"] for row in stage_summaries
            ),
            "total_charts_returning_finite_points": total_charts_with_points,
            "total_finite_point_occurrence_count": total_occurrences,
            "distinct_discoveries_with_primary_campaign_source": len(primary_sources),
            "new_distinct_discoveries_beyond_source_ledger": len(
                set(final_discoveries) - set(starting_discoveries)
            ),
            "starting_discovery_count": len(starting_discoveries),
            "final_discovery_count": len(final_discoveries),
            "exact_group_growth": False,
            "target_reached": False,
        },
        "interpretation": (
            "Under the frozen per-chart bounds, natural quotient weights one and two "
            "produce no exact rank or finite-index enlargement of the discovered M29 subgroup."
        ),
        "next_precommitted_stage": {
            "index": 3,
            "id": protocol["stage_plans_by_quotient_bit_count"]["12"][2]["id"],
            "requires_explicit_stage_limit_escalation": True,
            "new_chart_count": protocol["stage_plans_by_quotient_bit_count"]["12"][2][
                "new_chart_count"
            ],
        },
        "claim_boundary": [
            "This is complete accounting only for the two predeclared primary sparse stages at the recorded per-chart bounds.",
            "The bounded miss is not a rank upper bound and does not prove that the displayed M29 subgroup is saturated in E(Q).",
            "The campaign found many rational points, and exact integral relations place every discovered point in the same M29 subgroup.",
            "Rank at least 32 remains open; alternate quotient bases and higher natural weights were not searched in this campaign.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=LOCAL_LEDGER)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.check:
        if not args.archive.exists() or not args.output.exists():
            raise SystemExit("missing promoted primary-search artifact")
        raw = gzip.decompress(args.archive.read_bytes())
        if args.source.exists() and args.source.read_bytes() != raw:
            raise ArithmeticError("local ledger differs from the promoted full ledger")
        output = validate_and_summarize(raw)
        output["full_ledger"].update(
            {
                "path": relative(args.archive),
                "gzip_sha256": file_sha256(args.archive),
                "gzip_bytes": args.archive.stat().st_size,
            }
        )
        if args.output.read_bytes() != canonical_bytes(output):
            raise SystemExit(f"stale or invalid compact artifact: {args.output}")
        print(
            "C385SPARSEPRIMARY|status=PASS|rank=29|planned=3354|fresh=3116|"
            f"ledger={EXPECTED_LEDGER_SHA256}|manifest={file_sha256(args.output)}"
        )
        return

    if not args.source.exists():
        raise SystemExit(f"source search ledger does not exist: {args.source}")
    raw = args.source.read_bytes()
    output = validate_and_summarize(raw)
    args.archive.parent.mkdir(parents=True, exist_ok=True)
    args.archive.write_bytes(gzip.compress(raw, compresslevel=9, mtime=0))
    output["full_ledger"].update(
        {
            "path": relative(args.archive),
            "gzip_sha256": file_sha256(args.archive),
            "gzip_bytes": args.archive.stat().st_size,
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical_bytes(output))
    print(
        "C385SPARSEPRIMARY|status=WROTE|rank=29|planned=3354|fresh=3116|"
        f"ledger={EXPECTED_LEDGER_SHA256}|manifest={file_sha256(args.output)}"
    )


if __name__ == "__main__":
    main()
