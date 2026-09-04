#!/usr/bin/env python3
"""Freeze the immutable inputs and execution order for MW17-jump-v2.

This builder only selects from already frozen populations.  It does not
specialize a curve, run a point search, or overwrite any source population.
The resulting manifest is the sole candidate input to the exact detector.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
H10000 = GENERATED / "elliptic-curves/elkies_2026_compact_t_nagao_positive_control_h10000_v1.json"
CRT = GENERATED / "elkies-k3-r17-prospective-crt-frozen-cohorts-v1.json"
CRT_FEATURES = GENERATED / "elkies-k3-r17-prospective-crt-arithmetic-features-v1.json.gz"
ATLAS_HOLDOUT = GENERATED / "elkies-k3-r17-norm12-prospective-ordinary-family-holdout-v1.json"
LINEAGE = GENERATED / "elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
DIRECT = {
    family: GENERATED / f"elkies-k3-r17-norm12-orbit{family}-direct-fibration-v1.json"
    for family in ("07ca9", "08234", "08f72", "11952", "103b2")
}
ATLAS = GENERATED / "elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
RUNNER = ROOT / "elliptic-curves/cas/run_mw17_jump_v2.sage"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
POLICY = ROOT / "elliptic-curves/cas/half_lattice_chart_policy.py"
LADDER_V2 = ROOT / "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind_v2.sage"
LADDER_V1 = ROOT / "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind.sage"
OUTPUT = GENERATED / "elkies-k3-mw17-jump-v2-campaign-v1.json"

SCHEMA = "elkies-k3.mw17-jump-v2-campaign.v1"
STATUS = "FROZEN_SELECTED_IMMUTABLE_POPULATIONS_BEFORE_V2_EVALUATION"
H10000_TAKE = 191
CRT_PER_ANCHOR = 256
CRT_LANE_QUOTAS = {"full": 52, "C": 51, "D": 51, "E": 51, "F": 51}
FAMILY_ORDER = ("07ca9", "08234", "08f72", "11952", "074d9", "0e80b")
DIRECT_FAMILY = {**{family: family for family in FAMILY_ORDER}, "0e80b": "103b2"}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def load(path: Path):
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as handle:
        return json.load(handle)


def crt_lane(cohort: str) -> str:
    if cohort in ("A_356_full", "B_385_full"):
        return "full"
    lane = cohort.split("_", 1)[0]
    if lane not in {"C", "D", "E", "F"}:
        raise ArithmeticError(f"unknown CRT cohort {cohort}")
    return lane


def frozen_inputs() -> list[Path]:
    return [
        H10000,
        CRT,
        CRT_FEATURES,
        ATLAS_HOLDOUT,
        LINEAGE,
        ATLAS,
        *DIRECT.values(),
    ]


def implementation_inputs() -> list[Path]:
    return [Path(__file__), RUNNER, LEGACY, ENGINE, POLICY, LADDER_V1, LADDER_V2]


def build() -> dict[str, Any]:
    h10000 = load(H10000)
    crt = load(CRT)
    crt_features = load(CRT_FEATURES)
    atlas_holdout = load(ATLAS_HOLDOUT)
    if h10000.get("status") != "PASS_POSITIVE_CONTROL_SCORING_GATE":
        raise ArithmeticError("the H=10000 population is not the frozen scoring artifact")
    if h10000.get("population_count") != 121_589_944 or len(h10000["finalists"]) != 1_000:
        raise ArithmeticError("the H=10000 source inventory changed")
    if crt.get("status") != "FROZEN_UNOPENED_MATCHED_CRT_AND_ABLATION_COHORTS":
        raise ArithmeticError("the CRT source population is not frozen and unopened")
    if len(crt["rows"]) != 2_560:
        raise ArithmeticError("the CRT source inventory changed")
    if crt_features.get("status") != "PASS_COMPLETE_EXACT_PRESEARCH_ARITHMETIC_PANEL":
        raise ArithmeticError("the CRT feature panel is not complete")
    if atlas_holdout.get("status") != "FROZEN_UNOPENED_PROSPECTIVE_WHOLE_FAMILY_HOLDOUT":
        raise ArithmeticError("the atlas-family holdout is not frozen and unopened")
    if len(atlas_holdout["rows"]) != 1_536:
        raise ArithmeticError("the atlas-family holdout inventory changed")

    rows: list[dict[str, Any]] = []

    # Highest-priority source: all 256 ordinary holdout rows in each of the
    # six representative classes, in the declared class order.  Within a
    # class the original selection_counter is the only order.
    by_family: dict[str, list[tuple[int, dict[str, Any]]]] = defaultdict(list)
    for source_index, row in enumerate(atlas_holdout["rows"]):
        by_family[row["family"].rsplit("-", 1)[1]].append((source_index, row))
    if set(by_family) != set(FAMILY_ORDER):
        raise ArithmeticError("the prospective representative-family inventory changed")
    for priority, family in enumerate(FAMILY_ORDER):
        family_rows = sorted(
            by_family[family], key=lambda item: (int(item[1]["selection_counter"]), item[1]["sample_id"])
        )
        if len(family_rows) != 256:
            raise ArithmeticError(f"family {family} no longer has 256 frozen holdout rows")
        for family_position, (source_index, row) in enumerate(family_rows):
            rows.append(
                {
                    "sample_id": f"atlas:{row['sample_id']}",
                    "source_sample_id": row["sample_id"],
                    "source_population": "atlas_ordinary_family_holdout",
                    "source_index": source_index,
                    "source_position_within_population": family_position,
                    "priority_tranche": priority,
                    "family": family,
                    "direct_model_family": DIRECT_FAMILY[family],
                    "frame_class": row["frame_class"],
                    "parameter": row["parameter"],
                    "projective_pair": row["projective_pair"],
                    "selection_diagnostic": {
                        "kind": "frozen_ordinary_counter_hash_draw_no_score",
                        "selection_counter": row["selection_counter"],
                        "used_as_detector_filter": False,
                    },
                }
            )

    # The repository stores 1,000 ranked H=10000 finalists, not a distinct
    # 191-row file.  Freeze the requested tranche as the first 191 entries of
    # that immutable ranking and record the discrepancy explicitly below.
    for source_index, row in enumerate(h10000["finalists"][:H10000_TAKE]):
        rows.append(
            {
                "sample_id": f"h10000:{source_index:04d}",
                "source_sample_id": None,
                "source_population": "bounded_box_h10000_requested_finalist_tranche",
                "source_index": source_index,
                "source_position_within_population": source_index,
                "priority_tranche": len(FAMILY_ORDER),
                "family": "074d9",
                "direct_model_family": "074d9",
                "frame_class": "published-R17",
                "parameter": row["parameter"],
                "projective_pair": row["projective_pair"],
                "selection_diagnostic": {
                    "kind": "legacy_local_nagao_ranking_only",
                    "rank_in_stored_top1000": source_index + 1,
                    "raw_block_score_units_1e12": row["raw_block_score_units_1e12"],
                    "used_as_detector_filter": False,
                },
            }
        )

    # Treat the two anchors as the two CRT populations.  Allocate 256 rows to
    # each with lane sizes differing by at most one.  Local/Nagao features
    # choose order inside a lane only; they never affect the reported result.
    feature_by_id = {row["sample_id"]: row for row in crt_features["rows"]}
    if len(feature_by_id) != len(crt["rows"]):
        raise ArithmeticError("the CRT feature join is not one-to-one")
    crt_groups: dict[tuple[int, str], list[tuple[int, dict[str, Any], dict[str, Any]]]] = defaultdict(list)
    for source_index, row in enumerate(crt["rows"]):
        feature = feature_by_id.get(row["sample_id"])
        if feature is None or feature["parameter"] != row["parameter"]:
            raise ArithmeticError("the CRT feature join changed a parameter")
        crt_groups[(int(row["anchor_curve_id"]), crt_lane(row["cohort"]))].append(
            (source_index, row, feature)
        )
    for anchor_offset, anchor in enumerate((356, 385)):
        selected: list[tuple[str, int, int, dict[str, Any], dict[str, Any]]] = []
        for lane in ("full", "C", "D", "E", "F"):
            candidates = crt_groups[(anchor, lane)]
            candidates.sort(
                key=lambda item: (
                    -min(item[2]["nagao_comparison"]["block_score_units_1e12"]),
                    -item[2]["nagao_comparison"]["total_score_units_1e12"],
                    item[1]["sample_id"],
                )
            )
            take = CRT_LANE_QUOTAS[lane]
            if len(candidates) < take:
                raise ArithmeticError(f"CRT anchor {anchor} lane {lane} is too small")
            selected.extend((lane, lane_rank, *item) for lane_rank, item in enumerate(candidates[:take], 1))
        if len(selected) != CRT_PER_ANCHOR:
            raise ArithmeticError("the balanced CRT quota does not sum to 256")
        selected.sort(
            key=lambda item: (
                -min(item[4]["nagao_comparison"]["block_score_units_1e12"]),
                -item[4]["nagao_comparison"]["total_score_units_1e12"],
                item[3]["sample_id"],
            )
        )
        for anchor_position, (lane, lane_rank, source_index, row, feature) in enumerate(selected):
            nagao = feature["nagao_comparison"]
            rows.append(
                {
                    "sample_id": f"crt:{row['sample_id']}",
                    "source_sample_id": row["sample_id"],
                    "source_population": f"crt_anchor_{anchor}_balanced256",
                    "source_index": source_index,
                    "source_position_within_population": anchor_position,
                    "priority_tranche": len(FAMILY_ORDER) + 1 + anchor_offset,
                    "family": "074d9",
                    "direct_model_family": "074d9",
                    "frame_class": "published-R17",
                    "parameter": row["parameter"],
                    "projective_pair": row["projective_pair"],
                    "selection_diagnostic": {
                        "kind": "balanced_crt_lane_then_weakest_block_nagao",
                        "anchor_curve_id": anchor,
                        "cohort": row["cohort"],
                        "lane": lane,
                        "lane_rank": lane_rank,
                        "weakest_block_score_units_1e12": min(nagao["block_score_units_1e12"]),
                        "total_score_units_1e12": nagao["total_score_units_1e12"],
                        "used_as_detector_filter": False,
                    },
                }
            )

    for campaign_index, row in enumerate(rows):
        row["campaign_index"] = campaign_index
    if len(rows) != 2_239 or len({row["sample_id"] for row in rows}) != len(rows):
        raise ArithmeticError("the MW17-jump-v2 candidate inventory is not 2,239 unique addresses")

    source_counts = Counter(row["source_population"] for row in rows)
    family_counts = Counter(row["family"] for row in rows)
    document: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "campaign_name": "MW17-jump-v2",
        "candidate_count": len(rows),
        "candidate_list_sha256": canonical_hash(rows),
        "rows": rows,
        "selection": {
            "execution_order": [
                "atlas 07ca9",
                "atlas 08234",
                "alternate-Q80 08f72",
                "alternate-Q80 11952",
                "other published-R17 074d9",
                "other published-R17 0e80b (evaluated on exactly PGL2-equivalent 103b2 chart)",
                "legacy H=10000 requested finalist tranche",
                "balanced CRT anchor 356",
                "balanced CRT anchor 385",
            ],
            "counts_by_source_population": dict(sorted(source_counts.items())),
            "counts_by_representative_family": dict(sorted(family_counts.items())),
            "h10000_source_fact": {
                "scanned_parameter_count": 121_589_944,
                "stored_ranked_finalist_count": 1_000,
                "selected_count": H10000_TAKE,
                "selection_rule": "first 191 rows of the immutable stored ranking",
                "warning": (
                    "No distinct 191-row source object is present. This manifest does not claim that "
                    "the remaining 809 stored finalists were previously tested or excluded."
                ),
            },
            "crt_rule": {
                "population_definition": "one population per anchor curve id",
                "selected_per_anchor": CRT_PER_ANCHOR,
                "lane_quotas": CRT_LANE_QUOTAS,
                "within_lane_order": "decreasing minimum Nagao block, decreasing total Nagao score, sample id",
            },
        },
        "detector": {
            "name": "exact generic-depth to adaptive half-lattice detector",
            "generic_rank": 17,
            "complete_generic_half_class_census": 1 << 17,
            "initial_chart_count": 43,
            "adaptive_chart_count_after_nonzero_certified_initial_gain": 301,
            "height_bound_each_chart": 100_000,
            "wall_timeout_seconds_each_chart": 15.0,
            "gp_stack_bytes_each_chart": 1_000_000_000,
            "relation_chunk_size": 64,
            "relation_timeout_seconds": 180.0,
            "finite_reduction_prime_bound": 1_000,
            "normalization_policy": (
                "certify and classify on the p=2-minimal-to-short model; 074d9 also searches there, "
                "while other exact direct-family charts search on the direct short integral model "
                "with exact point transport between the two isomorphic models"
            ),
            "ranking_field": "actual_certified_quotient_rank_gain",
            "rank_lower_bound_field": "certified_rank_lower_bound",
            "rank_lower_bound_formula": "17 + actual_certified_quotient_rank_gain",
            "initial_gain_policy": (
                "Initial 43-chart gain is a detector phase and diagnostic, never a candidate-selection "
                "or leaderboard filter. Gains 1 through 17 enter the adaptive phase; initial zero remains "
                "a bounded detector result and does not reject the immutable candidate."
            ),
            "global_termination": "write stop sentinel immediately after any certified gain at least 15",
        },
        "immutability": {
            "source_files_are_read_only": True,
            "source_file_sha256": {relative(path): digest(path) for path in frozen_inputs()},
            "new_outputs": [
                "artifacts/local/elkies-k3/mw17-jump-v2/chunk-*.json",
                "artifacts/generated-results/elkies-k3-mw17-jump-v2-ledger-v1.json",
            ],
        },
        "implementation_hashes": {relative(path): digest(path) for path in implementation_inputs()},
        "resource_gate": {
            "explicit_large_campaign_scope": "the 2,239 selected fibres in this manifest only",
            "forbidden_rescan": "the 121,589,944-parameter H=10000 scan is not run by this campaign",
            "checkpoint_unit": "one completed fibre per chunk file rewrite",
            "worker_wall_timeout_seconds": 7_200,
            "worker_address_space_bytes": 16_000_000_000,
        },
        "claim_boundary": [
            "A score changes only after exact equation checks and an exact finite-reduction independence certificate.",
            "The detector certifies a rank lower bound for its discovered subgroup, not the full Mordell-Weil rank.",
            "Every zero or miss is bounded by the frozen chart, height, timeout, and backend budgets.",
            "Local and Nagao values only freeze selection order; they are absent from result ranking.",
        ],
        "generation": {
            "script": relative(Path(__file__)),
            "command": "python3 elkies-k3/scripts/build_mw17_jump_v2_campaign.py",
        },
    }
    definition = {key: value for key, value in document.items() if key != "protocol_definition_sha256"}
    document["protocol_definition_sha256"] = canonical_hash(definition)
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = build()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != text:
            raise SystemExit("MW17-jump-v2 campaign artifact is stale")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(
        f"MW17JUMPV2CAMPAIGN|status=FROZEN|candidates={payload['candidate_count']}|"
        f"sha256={digest(args.output)}|output={relative(args.output)}"
    )


if __name__ == "__main__":
    main()
