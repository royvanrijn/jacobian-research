#!/usr/bin/env python3
"""Consolidate the equation-first audit for the determinant-654 MW16 K3."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
ADAPTER = GEN / "elkies-k3-k3-14ad03cd7c1848b2-source-search-target-partner1-lattice-only-v1.json"
SOURCES = GEN / "elkies-k3-k3-14ad03cd7c1848b2-semistable-mw0-2-sources-large-a-partner1-v1.json"
FIBRE_I8_I9 = GEN / "elkies-k3-k3-14ad03cd7c1848b2-a7-a8-mw2-fibre-ansatz-mod5-v1.json"
CORRIDOR = GEN / "elkies-k3-k3-14ad03cd7c1848b2-same-ns-compiler-routes-rankfirst-cap2000-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-14ad03cd7c1848b2-equation-first-candidate-v1.json"

MARKING_STEMS = {
    "K3-14ad03cd7c1848b2-S0050": (
        "a2-a5-a8-mw2-s0050",
        ("mod5-square", "mod5-nonsquare", "mod7-square", "mod7-nonsquare"),
    ),
    "K3-14ad03cd7c1848b2-S0093": (
        "a2-a5-a8-mw2-s0093",
        ("mod5-square", "mod5-nonsquare", "mod7-square", "mod7-nonsquare"),
    ),
    "K3-14ad03cd7c1848b2-S0360": (
        "3a5-mw2-s0360",
        ("mod5-square", "mod5-nonsquare", "mod7-square", "mod7-nonsquare"),
    ),
    "K3-14ad03cd7c1848b2-S0071": (
        "a5-a10-mw2-s0071",
        ("mod5-square", "mod5-nonsquare"),
    ),
    "K3-14ad03cd7c1848b2-S0197": (
        "a7-a8-mw2-s0197",
        ("mod5-square", "mod5-nonsquare"),
    ),
}

EXPECTED_SOURCE_PROFILES = {
    "K3-14ad03cd7c1848b2-S0050": ("A2+A5+A8", [0, 0]),
    "K3-14ad03cd7c1848b2-S0093": ("A2+A5+A8", [0, 0]),
    "K3-14ad03cd7c1848b2-S0360": ("3A5", [0, 0]),
    "K3-14ad03cd7c1848b2-S0071": ("A10+A5", [0, 1]),
    "K3-14ad03cd7c1848b2-S0197": ("A7+A8", [0, 1]),
}


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def marking_paths() -> dict[str, list[Path]]:
    result = {}
    for source_id, (stem, tags) in MARKING_STEMS.items():
        result[source_id] = [
            GEN / f"elkies-k3-k3-14ad03cd7c1848b2-{stem}-marking-{tag}-v1.json"
            for tag in tags
        ]
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    adapter = load(ADAPTER)
    sources_payload = load(SOURCES)
    fibre = load(FIBRE_I8_I9)
    corridor = load(CORRIDOR)
    if adapter["surface_id"] != "K3-14ad03cd7c1848b2":
        raise ValueError("surface adapter changed")
    if adapter["determinant"] != 654 or adapter["frame"]["root_type"] != "A1":
        raise ValueError("target is no longer determinant 654 with root type A1")
    if adapter["equation_work_authorized"]:
        raise ValueError("lattice-only target unexpectedly authorizes equation work")

    sources = sources_payload["sources"]
    rank_histogram = Counter(row["source"]["mw_rank_for_rho_19"] for row in sources)
    if len(sources) != 420 or rank_histogram != Counter({1: 9, 2: 411}):
        raise ValueError("determinant-654 source census changed")
    complete_basis_rows = sum(
        int(
            row["source"]["mw_rank_for_rho_19"] == 2
            and bool(row["source"]["pole_audit"]["basis_with_all_poles_at_most_two"])
        )
        for row in sources
    )
    if complete_basis_rows != 127:
        raise ValueError("complete low-pole MW2 basis count changed")

    source_by_id = {row["source_id"]: row["source"] for row in sources}
    selected_sources = []
    for source_id, (root_type, poles) in EXPECTED_SOURCE_PROFILES.items():
        source = source_by_id[source_id]
        actual_poles = [row["pole_order"] for row in source["pole_audit"]["basis"]]
        if (
            source["root_type"] != root_type
            or source["mw_rank_for_rho_19"] != 2
            or source["torsion"] != 1
            or actual_poles != poles
        ):
            raise ValueError(f"selected source profile changed: {source_id}")
        selected_sources.append(
            {
                "source_id": source_id,
                "source_gram_sha256": source["gram_sha256"],
                "root_type": root_type,
                "support_count": source["support_count"],
                "basis_pole_profile": poles,
                "mw_height_gram": source["mw_height_gram"],
            }
        )

    if (
        fibre["status"] != "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ"
        or not fibre["scan"]["exhausted"]
        or fibre["ansatz"]["normalized_reducible_supports"]
        != ["0:I8", "infinity:I9"]
    ):
        raise ValueError("I8+I9 fibre census is not exhaustive")

    marking_summary = []
    all_marking_paths = marking_paths()
    for source_id, paths in all_marking_paths.items():
        scans = []
        for path in paths:
            marking = load(path)
            if (
                marking["source"]["source_id"] != source_id
                or marking["status"]
                != "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW2_BASIS_LOCUS"
                or not marking["scope"]["fibre_census_exhaustive"]
                or marking["accounting"]["marked_ordered_basis_pairs"] != 0
            ):
                raise ValueError(f"marking gate is not an exhaustive empty result: {path}")
            scans.append(
                {
                    "prime": marking["prime"],
                    "twist_square_class": marking["quadratic_twist_square_class"],
                    "marked_generator_sections": marking["accounting"]["marked_generator_sections"],
                    "models_with_both_generator_classes": marking["accounting"][
                        "models_with_both_generator_section_classes"
                    ],
                    "component_matched_pairs": marking["accounting"][
                        "component_matched_pair_candidates"
                    ],
                    "pairs_meeting_singular_fibres": marking["accounting"][
                        "pairs_meeting_singular_fibres"
                    ],
                    "pairs_with_wrong_smooth_intersection": marking["accounting"][
                        "pairs_with_wrong_smooth_intersection"
                    ],
                }
            )
        marking_summary.append({"source_id": source_id, "scans": scans})

    if corridor["status"] != "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_EMPTY":
        raise ValueError("same-NS corridor status changed")
    route = corridor["results"][0]
    depths = [
        {"depth": row["depth"], "best_root_rank": row["best_root_rank"]}
        for row in route["accounting"]
    ]
    if route["case"] != "k314ad" or route["best_routes_by_target"] or depths != [
        {"depth": 1, "best_root_rank": 13},
        {"depth": 2, "best_root_rank": 9},
        {"depth": 3, "best_root_rank": 4},
        {"depth": 4, "best_root_rank": 3},
        {"depth": 5, "best_root_rank": 3},
        {"depth": 6, "best_root_rank": 3},
        {"depth": 7, "best_root_rank": 3},
        {"depth": 8, "best_root_rank": 3},
    ]:
        raise ValueError("bounded route accounting changed")

    inputs = [ADAPTER, SOURCES, FIBRE_I8_I9, CORRIDOR]
    inputs.extend(path for paths in all_marking_paths.values() for path in paths)
    payload = {
        "schema": "elkies-k3.k3-14ad-equation-first-candidate.v1",
        "status": "PASS_LATTICE_SOURCE_SUPPLY_WITH_CHEAP_MARKING_GATES_EMPTY",
        "surface": {
            "surface_id": "K3-14ad03cd7c1848b2",
            "determinant": 654,
            "target_frame_id": adapter["frame"]["frame_id"],
            "target_root_type": "A1",
            "target_mw_rank_for_rho_19": 16,
            "t_arithmetic_gate": adapter["t_arithmetic_pre_solver_gate"],
            "equation_work_authorized": False,
        },
        "source_census": {
            "semistable_mw0_to_2_rows": len(sources),
            "mw_rank_histogram": {str(key): value for key, value in sorted(rank_histogram.items())},
            "mw2_rows_with_complete_basis_through_pole_two": complete_basis_rows,
        },
        "selected_cheapest_source_profiles": selected_sources,
        "normalized_marking_gates": marking_summary,
        "same_ns_corridor": {
            "source_id": route["source"]["source_id"],
            "search": route["search"],
            "best_root_rank_by_depth": depths,
            "target_hit": False,
        },
        "optimizer_decision": {
            "classification": "DEMOTE_BEHIND_DET500_FORMALLY_SMOOTH_MW1_SOURCE",
            "reason": (
                "The surface has abundant abstract MW1/MW2 sources, but none of the five "
                "cheapest complete marking profiles survives its exhaustive normalized "
                "finite-field gate. The selected pole-zero source also misses the MW16 "
                "target in the declared bounded degree-two beam."
            ),
            "remaining_open": (
                "Higher-pole source profiles, other normalizations or reduction behavior, "
                "wider/mixed-degree corridors, characteristic-zero equations, and the "
                "rootful-target D.F=2,3,4 multisection spectrum remain untested."
            ),
        },
        "inputs": {relative(path): digest(path) for path in inputs},
        "proof_boundary": {
            "proved": (
                "The exact lattice source counts, five selected low-pole profiles, every "
                "displayed exhaustive normalized finite-field marking census, and the "
                "declared bounded same-NS route search are replayed from pinned artifacts."
            ),
            "not_proved": (
                "The empty modular charts are not an absolute nonexistence theorem over Q, "
                "and the beam-pruned capped route miss is not a graph obstruction."
            ),
        },
        "reproduce": "python3 elkies-k3/scripts/certify_k3_14ad_equation_first_candidate.py",
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "K314ADEQUATIONFIRST|sources=420|profiles=5|marking_pairs=0|"
        "corridor_hit=0|status=PASS"
    )


if __name__ == "__main__":
    main()
