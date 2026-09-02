#!/usr/bin/env python3
"""Build an equation-first audit for the exact degree-three top-five surfaces.

The source-ranking ledger is intentionally rank first.  That can hide a
slightly higher-rank source whose complete Mordell--Weil basis has much smaller
section poles.  This report keeps the original rank-first leader, separate
minimum-pole and minimum-basis-pole leaders, a semistable minimum-basis-pole
leader, and the full nondominated source-metric frontier for each surface whose
degree-three spectrum has now been enumerated completely.

No rational marking or neighbour route is inferred from lattice data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RANKING = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-source-ranking-v2.json"
)
DEFAULT_DEGREE3 = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-current-source-top5-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-equation-first-degree3-top5-v1.json"
)
NS0028_MARKING_ARTIFACTS = [
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-pole0-section-pairs-mod5.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-pole0-section-pairs-mod5-nonsquare-twist.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-pole0-section-pairs-mod7.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0028-pole0-section-pairs-mod7-nonsquare-twist.json",
]
NS0005_FIBRE_ARTIFACTS = [
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod5-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod7-v1.json",
]
NS0005_POLE_ZERO_ARTIFACTS = [
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-infinity-pole0-sections-mod5-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-infinity-pole0-sections-mod5-nonsquare-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-infinity-pole0-sections-mod7-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-infinity-pole0-sections-mod7-nonsquare-v1.json",
]
NS0005_PAIR_ARTIFACTS = [
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-pole1-pole0-pairs-mod5-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-pole1-pole0-pairs-mod7-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-pole1-pole0-pairs-mod7-nonsquare-v1.json",
]
NS0005_MARKING_ARTIFACTS = [
    *NS0005_FIBRE_ARTIFACTS,
    *NS0005_POLE_ZERO_ARTIFACTS,
    *NS0005_PAIR_ARTIFACTS,
]
NS0031_MARKING_ARTIFACTS = [
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod5-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod5-nonsquare-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod7-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod7-nonsquare-v1.json",
]
NS0031_TANGENT_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-marked-gf7-hensel-v1.json"
)
NS0031_MULTISECTION_PILOT_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-rootless-mw17-multisection-pilot-v1.json"
)
NS0031_DEGREE3_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-ns0031-pilot-top5-v1.json"
)
UNKNOWN = 10**9


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def pole_value(record: dict, field: str) -> int:
    value = record.get(field)
    return UNKNOWN if value is None else int(value)


def source_metrics(row: dict) -> tuple[int, ...]:
    return (
        int(row["source_mw_rank"]),
        int(row["reducible_fibre_support_count"]),
        0 if row["semistable_configuration_compatible"] else 1,
        pole_value(row["minimum_nonzero_section_pole"], "pole_order"),
        pole_value(row["minimum_complete_mw_basis_pole"], "maximum_pole_order"),
        int(row["expected_additional_coefficient_conditions"]),
    )


def dominates(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    return all(a <= b for a, b in zip(left, right)) and any(
        a < b for a, b in zip(left, right)
    )


def compact_source(row: dict) -> dict:
    return {
        "source_id": row["source_id"],
        "source_artifact": row["source_artifact"],
        "source_mw_rank": int(row["source_mw_rank"]),
        "source_root_type": row["source_root_type"],
        "reducible_fibre_support_count": int(
            row["reducible_fibre_support_count"]
        ),
        "semistable_configuration_compatible": bool(
            row["semistable_configuration_compatible"]
        ),
        "expected_additional_coefficient_conditions": int(
            row["expected_additional_coefficient_conditions"]
        ),
        "minimum_nonzero_section_pole": row["minimum_nonzero_section_pole"],
        "minimum_complete_mw_basis_pole": row["minimum_complete_mw_basis_pole"],
        "arithmetic_marking": row["arithmetic_marking"],
        "certified_neighbor_route": row["certified_neighbor_route"],
        "metric_vector": [
            None if value == UNKNOWN else value for value in source_metrics(row)
        ],
    }


def select_canonical(rows: list[dict], key) -> dict:
    return min(
        rows,
        key=lambda row: (
            key(row),
            row["source_root_type"],
            row["source_artifact"],
            row["source_id"],
        ),
    )


def build(ranking_path: Path, degree3_path: Path) -> dict:
    ranking = json.loads(ranking_path.read_text())
    degree3 = json.loads(degree3_path.read_text())
    if ranking.get("schema") != "elkies-k3.lattice-foundry-source-ranking.v2":
        raise ValueError("unexpected source-ranking schema")
    if degree3.get("schema") != "elkies-k3.lattice-foundry-degree3-complete-spectrum.v1":
        raise ValueError("unexpected degree-three schema")

    selected_frames = degree3["selection"][
        "current_mw2_source_ranked_top_five_frame_ids"
    ]
    ns0028_marking_rows = []
    for path in NS0028_MARKING_ARTIFACTS:
        payload = json.loads(path.read_text())
        if (
            payload.get("schema")
            != "elkies-k3.lattice-foundry-ns0028-pole0-section-pairs-modp.v1"
        ):
            raise ValueError(f"unexpected NS0028 marking schema: {path}")
        ns0028_marking_rows.append(
            {
                "artifact": relative(path),
                "prime": int(payload["prime"]),
                "quadratic_twist_square_class": payload[
                    "quadratic_twist_square_class"
                ],
                "status": payload["status"],
                "models_with_P": int(payload["accounting"]["models_with_P"]),
                "models_with_Q": int(payload["accounting"]["models_with_Q"]),
                "models_with_marked_pairs": int(
                    payload["accounting"]["models_with_marked_pairs"]
                ),
                "total_P_sections": int(
                    payload["accounting"]["total_P_sections"]
                ),
                "total_Q_sections": int(
                    payload["accounting"]["total_Q_sections"]
                ),
                "total_marked_pairs": int(
                    payload["accounting"]["total_marked_pairs"]
                ),
            }
        )
    ns0005_fibre_rows = []
    for path in NS0005_FIBRE_ARTIFACTS:
        payload = json.loads(path.read_text())
        if (
            payload.get("schema")
            != "elkies-k3.lattice-foundry-three-support-semistable-source-ansatz-modp.v1"
        ):
            raise ValueError(f"unexpected NS0005 fibre schema: {path}")
        ns0005_fibre_rows.append(
            {
                "artifact": relative(path),
                "prime": int(payload["prime"]),
                "status": payload["status"],
                "squarefree_fibre_models": int(
                    payload["accounting"]["squarefree_examples_with_signs"]
                ),
            }
        )
    ns0005_pole_zero_rows = []
    for path in NS0005_POLE_ZERO_ARTIFACTS:
        payload = json.loads(path.read_text())
        if (
            payload.get("schema")
            != "elkies-k3.lattice-foundry-three-support-infinity-pole0-sections-modp.v1"
        ):
            raise ValueError(f"unexpected NS0005 pole-zero schema: {path}")
        ns0005_pole_zero_rows.append(
            {
                "artifact": relative(path),
                "prime": int(payload["prime"]),
                "quadratic_twist_square_class": payload[
                    "quadratic_twist_square_class"
                ],
                "status": payload["status"],
                "models_with_pole_zero_section": int(
                    payload["accounting"]["models_with_Q"]
                ),
                "pole_zero_sections": int(
                    payload["accounting"]["total_Q_sections"]
                ),
            }
        )
    ns0005_pair_rows = []
    for path in NS0005_PAIR_ARTIFACTS:
        payload = json.loads(path.read_text())
        if (
            payload.get("schema")
            != "elkies-k3.lattice-foundry-ns0005-pole1-pole0-pairs-modp.v1"
        ):
            raise ValueError(f"unexpected NS0005 pair schema: {path}")
        ns0005_pair_rows.append(
            {
                "artifact": relative(path),
                "prime": int(payload["prime"]),
                "quadratic_twist_square_class": payload[
                    "quadratic_twist_square_class"
                ],
                "status": payload["status"],
                "pole_one_sections_with_required_components": int(
                    payload["accounting"][
                        "pole_one_sections_with_required_components"
                    ]
                ),
                "smooth_pair_intersection_histogram": payload["accounting"][
                    "smooth_pair_intersection_histogram"
                ],
                "marked_mw2_pairs": int(
                    payload["accounting"]["marked_mw2_pairs"]
                ),
            }
        )
    ns0031_marking_rows = []
    for path in NS0031_MARKING_ARTIFACTS:
        payload = json.loads(path.read_text())
        if (
            payload.get("schema")
            != "elkies-k3.lattice-foundry-ns0031-a1-2a7-marking-modp.v1"
        ):
            raise ValueError(f"unexpected NS0031 marking schema: {path}")
        ns0031_marking_rows.append(
            {
                "artifact": relative(path),
                "prime": int(payload["prime"]),
                "quadratic_twist_square_class": payload[
                    "quadratic_twist_square_class"
                ],
                "status": payload["status"],
                "models_with_pole_zero_section": int(
                    payload["accounting"]["models_with_pole_zero_section"]
                ),
                "pole_zero_sections": int(
                    payload["accounting"]["pole_zero_sections"]
                ),
                "pole_one_sections": int(
                    payload["accounting"]["pole_one_sections"]
                ),
                "marked_mw2_pairs": int(
                    payload["accounting"]["marked_mw2_pairs"]
                ),
                "smooth_pair_intersection_histogram": payload["accounting"][
                    "smooth_pair_intersection_histogram"
                ],
            }
        )
    ns0031_sources = [
        row
        for row in ranking["candidates"]
        if row["ns_id"] == "NS0031"
        and row["source_id"] == "NS0031-S001"
        and row["source_root_type"] == "A1+2A7"
        and int(row["source_mw_rank"]) == 2
    ]
    if len(ns0031_sources) != 1:
        raise ArithmeticError(
            f"expected one NS0031 A1+2A7/MW2 source, found {len(ns0031_sources)}"
        )
    ns0031_source = ns0031_sources[0]
    ns0031_mw17_targets = sorted(
        target["frame_id"]
        for target in ns0031_source["admissible_high_rank_targets"]
        if int(target["mw_rank"]) == 17 and target["root_type"] == "0"
    )
    ns0031_tangent = json.loads(NS0031_TANGENT_ARTIFACT.read_text())
    if ns0031_tangent.get("schema") != (
        "elkies-k3.lattice-foundry-ns0031-marked-gf7-hensel.v1"
    ):
        raise ValueError("unexpected NS0031 marked-tangent schema")
    if ns0031_tangent.get("status") != (
        "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_LIFT_TO_REQUESTED_PRECISION"
    ):
        raise ArithmeticError("NS0031 marked-tangent gate did not pass")
    ns0031_pilot = json.loads(NS0031_MULTISECTION_PILOT_ARTIFACT.read_text())
    if ns0031_pilot.get("schema") != (
        "elkies-k3.lattice-foundry-multisection-spectrum.v1"
    ):
        raise ValueError("unexpected NS0031 multisection-pilot schema")
    ns0031_degree3 = json.loads(NS0031_DEGREE3_ARTIFACT.read_text())
    if ns0031_degree3.get("schema") != (
        "elkies-k3.lattice-foundry-degree3-complete-spectrum.v1"
    ):
        raise ValueError("unexpected NS0031 complete degree-three schema")
    ns0031_degree3_rows = sorted(
        ns0031_degree3["spectra"],
        key=lambda row: (
            -int(row["rational_trisection_translation_cosets"]),
            -int(row["genus_one_trisection_translation_cosets"]),
            row["frame_id"],
        ),
    )
    spectra = {row["frame_id"]: row for row in degree3["spectra"]}
    surfaces = []
    for frame_id in selected_frames:
        spectrum = spectra[frame_id]
        ns_id = spectrum["ns_id"]
        candidates = [
            row
            for row in ranking["candidates"]
            if row["ns_id"] == ns_id
            and int(row["source_mw_rank"]) <= 2
            and any(
                target["frame_id"] == frame_id
                and int(target["mw_rank"]) == 17
                and target["root_type"] == "0"
                for target in row["admissible_high_rank_targets"]
            )
        ]
        if not candidates:
            raise ArithmeticError(f"no MW<=2 source candidates for {frame_id}")

        rank_first = select_canonical(candidates, lambda row: tuple(row["score_tuple"]))
        exact_nonzero = [
            row
            for row in candidates
            if row["minimum_nonzero_section_pole"].get("pole_order") is not None
        ]
        exact_basis = [
            row
            for row in exact_nonzero
            if row["minimum_complete_mw_basis_pole"].get("maximum_pole_order")
            is not None
        ]
        if not exact_basis:
            raise ArithmeticError(f"no exact complete-basis pole audit for {frame_id}")
        exact_semistable_basis = [
            row for row in exact_basis if row["semistable_configuration_compatible"]
        ]
        if not exact_semistable_basis:
            raise ArithmeticError(
                f"no exact semistable complete-basis pole audit for {frame_id}"
            )
        ideal_cut = [
            row
            for row in exact_semistable_basis
            if int(row["reducible_fibre_support_count"]) <= 3
        ]
        if not ideal_cut:
            raise ArithmeticError(f"no source passes the ideal cut for {frame_id}")
        minimum_nonzero = select_canonical(
            exact_basis,
            lambda row: (
                pole_value(row["minimum_nonzero_section_pole"], "pole_order"),
                pole_value(
                    row["minimum_complete_mw_basis_pole"], "maximum_pole_order"
                ),
                int(row["source_mw_rank"]),
                int(row["reducible_fibre_support_count"]),
                0 if row["semistable_configuration_compatible"] else 1,
            ),
        )
        minimum_basis = select_canonical(
            exact_basis,
            lambda row: (
                pole_value(
                    row["minimum_complete_mw_basis_pole"], "maximum_pole_order"
                ),
                pole_value(row["minimum_nonzero_section_pole"], "pole_order"),
                int(row["source_mw_rank"]),
                int(row["reducible_fibre_support_count"]),
                0 if row["semistable_configuration_compatible"] else 1,
            ),
        )
        minimum_semistable_basis = select_canonical(
            exact_semistable_basis,
            lambda row: (
                pole_value(
                    row["minimum_complete_mw_basis_pole"], "maximum_pole_order"
                ),
                pole_value(row["minimum_nonzero_section_pole"], "pole_order"),
                int(row["source_mw_rank"]),
                int(row["reducible_fibre_support_count"]),
            ),
        )
        ideal_cut_leader = select_canonical(
            ideal_cut,
            lambda row: (
                pole_value(
                    row["minimum_complete_mw_basis_pole"], "maximum_pole_order"
                ),
                pole_value(row["minimum_nonzero_section_pole"], "pole_order"),
                int(row["source_mw_rank"]),
                int(row["reducible_fibre_support_count"]),
                int(row["expected_additional_coefficient_conditions"]),
            ),
        )

        metric_representatives = {}
        for row in candidates:
            metric = source_metrics(row)
            key = (metric, row["source_root_type"])
            current = metric_representatives.get(key)
            if current is None or (
                row["source_artifact"], row["source_id"]
            ) < (current["source_artifact"], current["source_id"]):
                metric_representatives[key] = row
        representatives = list(metric_representatives.values())
        frontier = [
            row
            for row in representatives
            if not any(
                dominates(source_metrics(other), source_metrics(row))
                for other in representatives
            )
        ]
        frontier.sort(
            key=lambda row: (
                source_metrics(row),
                row["source_root_type"],
                row["source_artifact"],
                row["source_id"],
            )
        )

        surfaces.append(
            {
                "ns_id": ns_id,
                "target_frame_id": frame_id,
                "target_mw_rank": 17,
                "target_root_type": "0",
                "determinant": int(spectrum["determinant"]),
                "complete_degree_three": {
                    "translation_cosets": int(spectrum["translation_cosets"]),
                    "rational_trisection_translation_cosets": int(
                        spectrum["rational_trisection_translation_cosets"]
                    ),
                    "genus_one_trisection_translation_cosets": int(
                        spectrum["genus_one_trisection_translation_cosets"]
                    ),
                    "maximum_coset_minimum_norm": int(
                        spectrum["maximum_coset_minimum_norm"]
                    ),
                },
                "mw_at_most_two_candidate_rows": len(candidates),
                "rank_first_leader": compact_source(rank_first),
                "minimum_nonzero_pole_leader": compact_source(minimum_nonzero),
                "minimum_complete_basis_pole_leader": compact_source(minimum_basis),
                "minimum_semistable_complete_basis_pole_leader": compact_source(
                    minimum_semistable_basis
                ),
                "ideal_mw_at_most_two_semistable_support_at_most_three_leader": compact_source(
                    ideal_cut_leader
                ),
                "modular_marking_evidence": (
                    {
                        "status": "PASS_EXACT_TWO_PRIME_NORMALIZED_CHARTS_WITH_NO_MARKED_PAIR",
                        "charts": ns0028_marking_rows,
                        "boundary": (
                            "This rejects only the displayed normalized GF(5) and "
                            "GF(7) charts; it is not a characteristic-zero obstruction."
                        ),
                    }
                    if ns_id == "NS0028"
                    else (
                        {
                            "status": (
                                "PASS_EXACT_TWO_PRIME_NORMALIZED_CHARTS_WITH_"
                                "INDIVIDUAL_GENERATORS_BUT_NO_MARKED_MW2_PAIR"
                            ),
                            "fibre_charts": ns0005_fibre_rows,
                            "pole_zero_charts": ns0005_pole_zero_rows,
                            "pair_charts": ns0005_pair_rows,
                            "boundary": (
                                "This rejects only the displayed normalized GF(5) "
                                "and GF(7) charts; it is not a characteristic-zero "
                                "obstruction."
                            ),
                        }
                        if ns_id == "NS0005"
                        else {
                            "status": "UNKNOWN_NOT_AUDITED_IN_THIS_REPORT",
                            "charts": [],
                        }
                    )
                ),
                "source_metric_frontier": [compact_source(row) for row in frontier],
            }
        )

    surfaces.sort(
        key=lambda row: (
            -row["complete_degree_three"][
                "rational_trisection_translation_cosets"
            ],
            -row["complete_degree_three"][
                "genus_one_trisection_translation_cosets"
            ],
            row["target_frame_id"],
        )
    )
    return {
        "schema": "elkies-k3.lattice-foundry-equation-first-degree3-top5.v1",
        "status": "PASS_EXACT_AGGREGATION_WITH_ARITHMETIC_MARKING_AND_ROUTE_OPEN",
        "inputs": {
            relative(path): digest(path)
            for path in [
                ranking_path,
                degree3_path,
                *NS0028_MARKING_ARTIFACTS,
                *NS0005_MARKING_ARTIFACTS,
                *NS0031_MARKING_ARTIFACTS,
                NS0031_TANGENT_ARTIFACT,
                NS0031_MULTISECTION_PILOT_ARTIFACT,
                NS0031_DEGREE3_ARTIFACT,
            ]
        },
        "selection": {
            "target_scope": (
                "the five rootless MW17 frames with complete degree-three spectra "
                "selected by the current MW2 source-ranked ledger"
            ),
            "source_scope": "all attached MW0--2 candidates in source-ranking v2",
            "metric_vector_minimize": [
                "source_mw_rank",
                "reducible_fibre_support_count",
                "nonsemistable_penalty",
                "minimum_nonzero_section_pole_order_or_infinity",
                "minimum_complete_mw_basis_pole_order_or_infinity",
                "expected_additional_coefficient_conditions",
            ],
            "policy": (
                "Retain the rank-first leader and the pole leaders separately; "
                "do not collapse incomparable equation-complexity tradeoffs to one weight."
            ),
        },
        "additional_equation_level_leads": [
            {
                "ns_id": "NS0031",
                "source": compact_source(ns0031_source),
                "rootless_mw17_target_frame_ids": ns0031_mw17_targets,
                "complete_degree_three_spectrum": {
                    "status": ns0031_degree3["status"],
                    "artifact": relative(NS0031_DEGREE3_ARTIFACT),
                    "pilot_artifact": relative(
                        NS0031_MULTISECTION_PILOT_ARTIFACT
                    ),
                    "selection": (
                        "The five pilot leaders by sampled rational trisections, "
                        "breaking ties by sampled genus-one trisections; final "
                        "ordering uses only complete counts."
                    ),
                    "audited_frame_count": len(ns0031_degree3_rows),
                    "total_rootless_mw17_frame_count": len(ns0031_mw17_targets),
                    "spectra_by_exact_rational_trisection_count": [
                        {
                            "frame_id": row["frame_id"],
                            "translation_cosets": int(row["translation_cosets"]),
                            "rational_trisection_translation_cosets": int(
                                row["rational_trisection_translation_cosets"]
                            ),
                            "genus_one_trisection_translation_cosets": int(
                                row["genus_one_trisection_translation_cosets"]
                            ),
                            "maximum_coset_minimum_norm": int(
                                row["maximum_coset_minimum_norm"]
                            ),
                        }
                        for row in ns0031_degree3_rows
                    ],
                },
                "modular_marking_evidence": {
                    "status": (
                        "PASS_EXACT_GF7_SQUARE_TWIST_WITH_MARKED_MW2_PAIRS;_"
                        "GF5_BOTH_TWISTS_AND_GF7_NONSQUARE_TWIST_EMPTY"
                    ),
                    "charts": ns0031_marking_rows,
                    "boundary": (
                        "The GF(7) square-twist chart proves existence of two "
                        "complete marked pairs only in that finite-field normalized "
                        "chart. It does not construct a characteristic-zero family, "
                        "a Q-rational marking, or a neighbour corridor."
                    ),
                },
                "marked_tangent_evidence": {
                    "artifact": relative(NS0031_TANGENT_ARTIFACT),
                    "status": ns0031_tangent["status"],
                    "variables": int(ns0031_tangent["system"]["variable_count"]),
                    "equations": int(ns0031_tangent["system"]["equation_count"]),
                    "jacobian_rank_mod_7": int(
                        ns0031_tangent["jacobian_certificate"]["rank_mod_7"]
                    ),
                    "tangent_dimension": int(
                        ns0031_tangent["jacobian_certificate"]["tangent_dimension"]
                    ),
                    "maximal_minor_mod_7": int(
                        ns0031_tangent["jacobian_certificate"][
                            "pivot_minor_determinant_mod_7"
                        ]
                    ),
                    "explicit_lift_precision_exponent": int(
                        ns0031_tangent["finite_precision_lift"][
                            "achieved_precision_exponent"
                        ]
                    ),
                    "explicit_lift_modulus": int(
                        ns0031_tangent["finite_precision_lift"]["modulus"]
                    ),
                    "boundary": ns0031_tangent["proof_boundary"]["not_proved"],
                },
                "priority": (
                    "FIRST_CURRENT_EQUATION_LEVEL_POSITIVE_MW2_TO_ROOTLESS_MW17_"
                    "SOURCE_AND_STRONGEST_COMPLETE_RATIONAL_TRISECTION_COUNT;_"
                    "PROVE_THE_FORMAL_LIFT_AND_ENUMERATE_A_NEIGHBOUR_CORRIDOR"
                ),
            }
        ],
        "surfaces": surfaces,
        "proof_boundary": {
            "proved": (
                "All displayed lattice-source, section-pole, and degree-three counts "
                "are copied from hash-pinned exact artifacts. Frontier dominance is "
                "an exact comparison of the declared integer/Boolean metrics."
            ),
            "not_proved": (
                "No source marking is descended to Q, no Weierstrass equation or "
                "rational parameterization is constructed, and no neighbour corridor "
                "or specialization rank jump is certified."
            ),
        },
        "reproduce": (
            "python3 elkies-k3/scripts/build_lattice_foundry_equation_first_shortlist.py"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ranking", type=Path, default=DEFAULT_RANKING)
    parser.add_argument("--degree3", type=Path, default=DEFAULT_DEGREE3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    ranking_path = arguments.ranking.resolve()
    degree3_path = arguments.degree3.resolve()
    output_path = arguments.output.resolve()
    result = build(ranking_path, degree3_path)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("equation-first degree-three shortlist is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "EQUATIONFIRSTD3|"
        f"surfaces={len(result['surfaces'])}|"
        f"frontier_rows={sum(len(row['source_metric_frontier']) for row in result['surfaces'])}|"
        f"status={result['status']}"
    )


if __name__ == "__main__":
    main()
