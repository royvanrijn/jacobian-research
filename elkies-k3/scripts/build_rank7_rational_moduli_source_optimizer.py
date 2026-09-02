#!/usr/bin/env python3
"""Build the live source-first queue for rational-moduli MW17 surfaces."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
PARETO = GEN / "elkies-k3-rank7-surface-pareto-v1.json"
K304_PROMOTION = GEN / "elkies-k3-k3-04b86146cc6b284b-equation-first-promotion-v1.json"
K304_CORRIDOR = (
    GEN
    / "elkies-k3-k3-04b86146cc6b284b-same-ns-compiler-routes-rankfirst-cap2000-v1.json"
)
K310_REJECTION = GEN / "elkies-k3-k3-10a14a46c14b3150-semistable-source-rejection-v1.json"
K342_ADAPTERS = (
    GEN / "elkies-k3-k3-3425921cd7db891f-source-search-target-partner1-v1.json",
    GEN / "elkies-k3-k3-3425921cd7db891f-source-search-target-partner2-v1.json",
)
K342_SEARCHES = (
    GEN / "elkies-k3-k3-3425921cd7db891f-semistable-mw0-2-sources-large-a-partner1-v1.json",
    GEN / "elkies-k3-k3-3425921cd7db891f-semistable-mw0-2-sources-large-a-partner2-v1.json",
)
K349_ADAPTERS = (
    GEN / "elkies-k3-k3-49b947f9626a0481-source-search-target-partner1-v1.json",
    GEN / "elkies-k3-k3-49b947f9626a0481-source-search-target-partner2-v1.json",
)
K349_SEARCHES = (
    GEN / "elkies-k3-k3-49b947f9626a0481-semistable-mw0-2-sources-large-a-partner1-v1.json",
    GEN / "elkies-k3-k3-49b947f9626a0481-semistable-mw0-2-sources-large-a-partner2-v1.json",
)
DEFAULT_OUTPUT = GEN / "elkies-k3-rank7-rational-moduli-source-optimizer-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_negative(search: dict) -> None:
    if search["status"] != "PASS_EXACT_SEARCH_NO_SUCCESS_HIT":
        raise ArithmeticError("expected a scoped negative ideal-source search")
    scope = search["search_scope"]
    if not scope["complete"] or not scope["all_a_only"]:
        raise ArithmeticError("ideal-source search is not complete in its declared cut")
    if scope["source_root_rank"] != [15, 17] or scope["source_support_count"] != [1, 3]:
        raise ArithmeticError("ideal source window changed")
    if search["accounting"]["distinct_reduced_gram_sources"]:
        raise ArithmeticError("negative source inventory is inconsistent")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    pareto = json.loads(PARETO.read_text())
    k304 = json.loads(K304_PROMOTION.read_text())
    k304_corridor = json.loads(K304_CORRIDOR.read_text())
    k310 = json.loads(K310_REJECTION.read_text())
    k342_adapters = [json.loads(path.read_text()) for path in K342_ADAPTERS]
    k342_searches = [json.loads(path.read_text()) for path in K342_SEARCHES]
    k349_adapters = [json.loads(path.read_text()) for path in K349_ADAPTERS]
    k349_searches = [json.loads(path.read_text()) for path in K349_SEARCHES]

    rational_mw17 = sorted(
        (
            row
            for row in pareto["surfaces"]
            if row["maximum_generic_mw_rank"] == 17
            and row["moduli"]["rationality"] is True
        ),
        key=lambda row: (row["determinant"], row["surface_id"]),
    )
    expected_ids = [
        "K3-04b86146cc6b284b",
        "K3-10a14a46c14b3150",
        "K3-3425921cd7db891f",
        "K3-49b947f9626a0481",
        "K3-99a0b9b18de6e19b",
        "K3-dc0e324e4ac40dbc",
    ]
    if [row["surface_id"] for row in rational_mw17] != expected_ids:
        raise ArithmeticError("rational-moduli MW17 inventory changed")
    if k304["status"] != "PASS_FORMALLY_SMOOTH_Z7_PROMOTION_FIRST_QQ_POINT_REJECTED_DET20":
        raise ArithmeticError("determinant-500 promotion status changed")
    if k304_corridor["status"] != "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_EMPTY":
        raise ArithmeticError("determinant-500 bounded corridor status changed")
    corridor_result = k304_corridor["results"]
    if len(corridor_result) != 1 or corridor_result[0]["case"] != "k304b":
        raise ArithmeticError("determinant-500 bounded corridor case changed")
    if corridor_result[0]["best_routes_by_target"]:
        raise ArithmeticError("determinant-500 bounded corridor unexpectedly has a hit")
    if k310["status"] != "PASS_COMPLETE_TWO_AUXILIARY_LARGE_AMBIENT_IDEAL_SOURCE_REJECTION":
        raise ArithmeticError("determinant-750 rejection status changed")
    if [row["auxiliary"]["partner_index_one_based"] for row in k342_adapters] != [1, 2]:
        raise ArithmeticError("determinant-864 auxiliary coverage changed")
    for search in k342_searches:
        exact_negative(search)
    if [row["auxiliary"]["partner_index_one_based"] for row in k349_adapters] != [1, 2]:
        raise ArithmeticError("determinant-1296 auxiliary coverage changed")
    for search in k349_searches:
        exact_negative(search)

    catalogue_rows = {}
    for row in rational_mw17:
        catalogue_rows[row["surface_id"]] = {
            "surface_id": row["surface_id"],
            "determinant": row["determinant"],
            "target_mw_rank": row["maximum_generic_mw_rank"],
            "target_frame_ids": row["maximum_mw_frame_ids"],
            "modular_curve": row["moduli"]["t_arithmetic"]["curve_label"],
            "moduli_genus": row["moduli"]["genus"],
            "moduli_rational": row["moduli"]["rationality"],
            "auxiliary_classes": row["partner_auxiliary_count"],
        }

    active = {
        **catalogue_rows[expected_ids[0]],
        "source_status": "FORMALLY_SMOOTH_MARKED_MW1_BRANCH_RATIONAL_POINT_STILL_OPEN",
        "best_source": {
            "root_type": "A3+A4+A9",
            "mw_rank": 1,
            "support_count": 3,
            "semistable": True,
            "basis_pole_profile": [1],
        },
        "bounded_rational_scan": k304["equation_gate"]["bounded_integral_parameter_scan"],
        "same_surface_mw2_fallback": k304["same_surface_mw2_fallback"]["decision"],
        "same_ns_corridor": {
            "status": k304_corridor["status"],
            "search": corridor_result[0]["search"],
            "best_root_rank_by_depth": [
                [row["depth"], row["best_root_rank"]]
                for row in corridor_result[0]["accounting"]
            ],
        },
        "priority": 1,
    }
    rejected = [
        {
            **catalogue_rows[expected_ids[1]],
            "source_status": "NO_SEMISTABLE_MW0_2_SOURCE_IN_COMPLETE_TWO_AUXILIARY_SIX_LARGE_AMBIENT_CUT",
            "priority": None,
        },
        {
            **catalogue_rows[expected_ids[2]],
            "source_status": "NO_SEMISTABLE_MW0_2_SOURCE_IN_COMPLETE_TWO_AUXILIARY_SIX_LARGE_AMBIENT_CUT",
            "priority": None,
        },
        {
            **catalogue_rows[expected_ids[3]],
            "source_status": "NO_SEMISTABLE_MW0_2_SOURCE_IN_COMPLETE_TWO_AUXILIARY_SIX_LARGE_AMBIENT_CUT",
            "priority": None,
        },
    ]
    queued = [
        {
            **catalogue_rows[surface_id],
            "source_status": "IDEAL_SOURCE_SEARCH_NOT_YET_RUN",
            "priority": index,
        }
        for index, surface_id in enumerate(expected_ids[4:], start=2)
    ]

    payload = {
        "schema": "elkies-k3.rank7-rational-moduli-source-optimizer.v1",
        "status": "PASS_LIVE_SOURCE_FIRST_QUEUE",
        "objective_order": [
            "rational moduli and target MW rank",
            "source MW rank",
            "semistability and reducible support count",
            "rational marking and section pole profile",
            "coefficient conditions",
            "neighbor corridor cost",
            "target multisection richness",
        ],
        "accounting": {
            "rational_moduli_rootless_mw17_surfaces": len(rational_mw17),
            "active_promotions": 1,
            "scoped_ideal_source_rejections": 3,
            "queued_ideal_source_searches": len(queued),
        },
        "active_candidate": active,
        "next_search_queue": queued,
        "scoped_rejections": rejected,
        "decision": (
            "Continue rational algebraization on determinant 500 and widen its corridor "
            "search beyond the completed q=4 degree-2 pole-1 beam. Test determinant 1500 "
            "next. Determinants 750, 864, and 1296 should not receive equation work unless "
            "the source cut is widened."
        ),
        "proof_boundary": {
            "proved": (
                "The six rational-moduli rootless MW17 surfaces and their ordering are copied "
                "from the exact Pareto ledger. The determinant-500 source gates and the "
                "bounded determinant-500 corridor miss, and the declared two-auxiliary "
                "ideal-source misses at determinants 750, 864, and 1296 are "
                "hash-pinned exact computations."
            ),
            "not_proved": (
                "The queue is an optimization policy, not a completeness theorem. The two "
                "negative rows exclude only the declared six-large-ambient semistable MW0-2 "
                "cut, and the queued surfaces have no new low-MW source search yet."
            ),
        },
        "inputs": {
            relative(PARETO): digest(PARETO),
            relative(K304_PROMOTION): digest(K304_PROMOTION),
            relative(K304_CORRIDOR): digest(K304_CORRIDOR),
            relative(K310_REJECTION): digest(K310_REJECTION),
            **{relative(path): digest(path) for path in K342_ADAPTERS + K342_SEARCHES},
            **{relative(path): digest(path) for path in K349_ADAPTERS + K349_SEARCHES},
        },
        "reproduce": "python3 elkies-k3/scripts/build_rank7_rational_moduli_source_optimizer.py",
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print("RANK7RATIONALOPT|active=1|rejected=3|queued=2|status=PASS")


if __name__ == "__main__":
    main()
