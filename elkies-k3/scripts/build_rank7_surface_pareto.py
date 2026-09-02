#!/usr/bin/env python3
"""Build a typed Pareto ledger for the rank-seven auxiliary catalogue.

The universal frontier uses only metrics available for every imported
surface: maximum catalogued MW rank, easiest *known* exact source MW rank,
source reducible-fibre support count, and determinant.  A known external
source certificate outranks the surface's own catalogued frames; otherwise
the easiest catalogued frame is used and its limitation is explicit.

Equation, field, corridor, short-vector, symmetry, and multisection data are
never imputed.  Separate enriched frontiers compare only rows sharing the
required evidence.  This is a discovery ordering, not a completeness or
arithmetic-rank theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CATALOGUE = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json"
)
SOURCE_RANKING = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-source-ranking-v2.json"
)
MULTISECTIONS = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-multisection-spectrum-v1.json"
)
UMBRAL = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-umbral-orbits-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-surface-pareto-v1.json"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def support_count(root_type: str) -> int:
    if root_type == "0":
        return 0
    result = 0
    for component in root_type.split("+"):
        match = re.fullmatch(r"(?:(\d+))?([ADE]\d+)", component)
        if match is None:
            raise ValueError(f"unrecognized root component {component!r}")
        result += int(match.group(1) or 1)
    return result


def legacy_frame_ids(frame: dict) -> list[str]:
    return sorted(
        {
            row["legacy_frame_id"]
            for row in frame["provenance"]
            if "legacy_frame_id" in row
        }
    )


def dominates(left: dict, right: dict, fields: tuple[str, ...]) -> bool:
    values_left = tuple(left[field] for field in fields)
    values_right = tuple(right[field] for field in fields)
    return all(a <= b for a, b in zip(values_left, values_right)) and any(
        a < b for a, b in zip(values_left, values_right)
    )


def pareto_layers(records: list[dict], fields: tuple[str, ...]) -> list[list[str]]:
    remaining = list(records)
    layers = []
    while remaining:
        frontier = [
            row
            for row in remaining
            if not any(
                other is not row and dominates(other, row, fields)
                for other in remaining
            )
        ]
        assert frontier
        frontier.sort(
            key=lambda row: (
                tuple(row[field] for field in fields), row["surface_id"]
            )
        )
        layers.append([row["surface_id"] for row in frontier])
        frontier_ids = {id(row) for row in frontier}
        remaining = [row for row in remaining if id(row) not in frontier_ids]
    return layers


def route_metrics(route: dict) -> dict:
    if route.get("status") != "PASS_EXACT_CERTIFIED_ROUTE":
        return {
            "status": route.get("status", "UNKNOWN_NOT_YET_ENUMERATED"),
            "cost": None,
        }
    return {"status": route["status"], "cost": route["cost"]}


def symmetry_evidence(
    surface: dict, umbral_by_legacy_frame: dict[str, dict]
) -> dict:
    evidence = []
    for frame in surface["frames"]:
        for legacy_id in legacy_frame_ids(frame):
            if legacy_id not in umbral_by_legacy_frame:
                continue
            target = umbral_by_legacy_frame[legacy_id]
            full = target["full_ambient_stabilizer"]
            evidence.append(
                {
                    "frame_id": frame["frame_id"],
                    "legacy_frame_id": legacy_id,
                    "status": "PASS_EXACT_FULL_AMBIENT_STABILIZER",
                    "ambient_stabilizer_order": int(full["order"]),
                    "umbral_image_order": int(full["umbral_image_order"]),
                    "umbral_image_classes": full["umbral_image_classes"],
                }
            )
        for provenance in frame["provenance"]:
            stabilizer = provenance.get("literal_section_stabilizer")
            if stabilizer is None:
                stabilizer = provenance.get("literal_residual_stabilizer")
            if stabilizer is None:
                continue
            evidence.append(
                {
                    "frame_id": frame["frame_id"],
                    "source_frame_id": provenance.get("source_frame_id"),
                    "status": "PASS_EXACT_LITERAL_SECTION_STABILIZER_LOWER_BOUND",
                    "ambient_stabilizer_order_lower_bound": len(stabilizer),
                    "classes": [row["class"] for row in stabilizer],
                    "moved_dimensions_mod_2": [
                        row["moved_dimension_mod_2"] for row in stabilizer
                    ],
                }
            )
    if not evidence:
        return {
            "status": "UNKNOWN_NOT_COMPUTED",
            "evidence": [],
            "certified_order_lower_bound": 1,
            "lower_bound_is_identity_only": True,
        }
    lower_bounds = [
        row.get(
            "ambient_stabilizer_order",
            row.get("ambient_stabilizer_order_lower_bound", 1),
        )
        for row in evidence
    ]
    return {
        "status": "PASS_EXACT_AVAILABLE_STABILIZER_EVIDENCE",
        "evidence": evidence,
        "certified_order_lower_bound": max(lower_bounds),
        "lower_bound_is_identity_only": False,
    }


def target_short_vector_evidence(surface: dict) -> dict:
    maximum_mw = max(frame["mw_rank_for_rho_19"] for frame in surface["frames"])
    rows = []
    for frame in surface["frames"]:
        if frame["mw_rank_for_rho_19"] != maximum_mw:
            continue
        intrinsics = frame.get("rootless_intrinsics") or {}
        if not intrinsics:
            continue
        rows.append(
            {
                "frame_id": frame["frame_id"],
                "status": "PASS_EXACT_ROOTLESS_FRAME_INTRINSICS",
                "minimum_squared_norm": intrinsics.get("minimum_squared_norm"),
                "norm_four_vectors": intrinsics.get("norm_four_vectors"),
                "norm_four_unoriented_pairs": intrinsics.get(
                    "norm_four_unoriented_pairs"
                ),
                "automorphism_group_order": intrinsics.get(
                    "automorphism_group_order"
                ),
                "theta_coefficients_by_squared_norm_through_bound": intrinsics.get(
                    "theta_coefficients_by_squared_norm_through_bound"
                ),
            }
        )
    return {
        "status": (
            "PASS_EXACT_AVAILABLE_TARGET_SHORT_VECTOR_EVIDENCE"
            if rows
            else "UNKNOWN_MW_SHORT_VECTOR_SPECTRUM_NOT_COMPUTED"
        ),
        "frames": rows,
    }


def best_multisection_evidence(
    surface: dict, multisection_by_legacy_frame: dict[str, dict]
) -> dict:
    rows = []
    for frame in surface["frames"]:
        for legacy_id in legacy_frame_ids(frame):
            if legacy_id not in multisection_by_legacy_frame:
                continue
            target = multisection_by_legacy_frame[legacy_id]
            rows.append(
                {
                    "frame_id": frame["frame_id"],
                    "legacy_frame_id": legacy_id,
                    "mw_rank": frame["mw_rank_for_rho_19"],
                    "richness_coordinates": target["richness_coordinates"],
                    "degree_two": target["degree_two"],
                    "degree_three_sample": target["degree_three_sample"],
                    "degree_four_sample": target["degree_four_sample"],
                }
            )
    return {
        "status": (
            "PASS_EXACT_D2_BOUNDED_D3_D4_EVIDENCE"
            if rows
            else "UNKNOWN_LOW_DEGREE_MULTISECTION_SPECTRUM_NOT_COMPUTED"
        ),
        "targets": rows,
        "sampling_boundary": (
            "Degree two is exact in the source artifact; degree three and degree "
            "four retain their declared bounded/sample status."
            if rows
            else None
        ),
    }


def build(
    catalogue: dict,
    source_ranking: dict,
    multisections: dict,
    umbral: dict,
) -> dict:
    assert catalogue["schema"] == "elkies-k3.rank7-auxiliary-catalogue.v1"
    assert catalogue["status"] == (
        "PASS_EXACT_SURFACE_FIRST_IMPORTED_CATALOGUE_FULL_ORBIT_CENSUS_OPEN"
    )
    assert source_ranking["schema"] == "elkies-k3.lattice-foundry-source-ranking.v2"
    assert source_ranking["status"] == (
        "PASS_EXACT_SOURCE_METRICS_WITH_TYPED_OPEN_ARITHMETIC_AND_ROUTE_GATES"
    )
    assert multisections["schema"] == (
        "elkies-k3.lattice-foundry-multisection-spectrum.v1"
    )
    assert umbral["schema"] == "elkies-k3.lattice-foundry-umbral-orbits.v1"

    source_by_ns = {row["ns_id"]: row for row in source_ranking["surface_leaders"]}
    multisection_by_frame = {
        row["frame_id"]: row for row in multisections["targets"]
    }
    umbral_by_frame = {row["frame_id"]: row for row in umbral["targets"]}

    records = []
    for surface in catalogue["surfaces"]:
        frames = surface["frames"]
        maximum_mw = max(frame["mw_rank_for_rho_19"] for frame in frames)
        maximum_frames = [
            frame for frame in frames if frame["mw_rank_for_rho_19"] == maximum_mw
        ]
        external_sources = [
            source_by_ns[ns_id]
            for ns_id in surface["legacy_ns_ids"]
            if ns_id in source_by_ns
        ]
        if external_sources:
            source = min(
                external_sources,
                key=lambda row: (
                    row["source_mw_rank"],
                    row["reducible_fibre_support_count"],
                    row["source_root_type"],
                ),
            )
            source_record = {
                "status": "PASS_EXACT_EXTERNAL_SOURCE_CERTIFICATE",
                "source_kind": "external_rootful_source_inventory",
                "legacy_ns_id": source["ns_id"],
                "source_id": source["source_id"],
                "source_artifact": source["source_artifact"],
                "certificate_scope": source["certificate_scope"],
                "root_type": source["source_root_type"],
                "root_rank": int(source["source_root_rank"]),
                "mw_rank": int(source["source_mw_rank"]),
                "reducible_fibre_support_count": int(
                    source["reducible_fibre_support_count"]
                ),
                "semistable_configuration_compatible": source[
                    "semistable_configuration_compatible"
                ],
                "minimum_nonzero_section_pole": source[
                    "minimum_nonzero_section_pole"
                ],
                "arithmetic_marking": source["arithmetic_marking"],
                "expected_additional_coefficient_conditions": source[
                    "expected_additional_coefficient_conditions"
                ],
                "certified_neighbor_route": route_metrics(
                    source["certified_neighbor_route"]
                ),
            }
        else:
            source_frame = min(
                frames,
                key=lambda frame: (
                    frame["mw_rank_for_rho_19"],
                    support_count(frame["root_type"]),
                    frame["frame_id"],
                ),
            )
            source_record = {
                "status": "PASS_EXACT_CATALOGUE_FRAME_ONLY_SOURCE",
                "source_kind": "easiest_currently_catalogued_frame",
                "frame_id": source_frame["frame_id"],
                "root_type": source_frame["root_type"],
                "root_rank": int(source_frame["root_rank"]),
                "mw_rank": int(source_frame["mw_rank_for_rho_19"]),
                "reducible_fibre_support_count": support_count(
                    source_frame["root_type"]
                ),
                "semistable_configuration_compatible": all(
                    component.startswith("A")
                    or re.fullmatch(r"\d+A\d+", component)
                    for component in source_frame["root_type"].split("+")
                ),
                "minimum_nonzero_section_pole": {
                    "status": "UNKNOWN_NOT_COMPUTED",
                    "pole_order": None,
                },
                "arithmetic_marking": {
                    "status": "UNKNOWN_NOT_INFERRED_FROM_LATTICE",
                    "rational_source_marking": None,
                    "characteristic_zero_galois_orbit_size": None,
                    "rational_parameterization": None,
                },
                "expected_additional_coefficient_conditions": source_frame[
                    "mw_rank_for_rho_19"
                ],
                "certified_neighbor_route": {
                    "status": "UNKNOWN_NOT_YET_ENUMERATED",
                    "cost": None,
                },
                "scope_warning": (
                    "No external low-MW source inventory is attached to this surface; "
                    "the displayed source is only the easiest imported frame."
                ),
            }

        symmetry = symmetry_evidence(surface, umbral_by_frame)
        short_vectors = target_short_vector_evidence(surface)
        multisection = best_multisection_evidence(
            surface, multisection_by_frame
        )
        rank_gain = maximum_mw - source_record["mw_rank"]
        record = {
            "surface_id": surface["surface_id"],
            "legacy_ns_ids": surface["legacy_ns_ids"],
            "determinant": int(surface["determinant"]),
            "determinant_band": surface["determinant_band"],
            "maximum_generic_mw_rank": maximum_mw,
            "maximum_mw_frame_ids": [frame["frame_id"] for frame in maximum_frames],
            "easiest_known_source": source_record,
            "known_source_to_target_rank_gain": rank_gain,
            "catalogued_frame_count": len(frames),
            "partner_auxiliary_count": len(surface["partner_auxiliaries"]),
            "source_equation_complexity": {
                "status": (
                    "PARTIAL_EXACT_LATTICE_METRICS_EQUATION_NOT_CONSTRUCTED"
                ),
                "source_mw_rank": source_record["mw_rank"],
                "reducible_fibre_support_count": source_record[
                    "reducible_fibre_support_count"
                ],
                "minimum_nonzero_section_pole": source_record[
                    "minimum_nonzero_section_pole"
                ],
                "expected_additional_coefficient_conditions": source_record[
                    "expected_additional_coefficient_conditions"
                ],
            },
            "source_field_of_definition": source_record["arithmetic_marking"],
            "shortest_physical_neighbor_corridor": source_record[
                "certified_neighbor_route"
            ],
            "target_mw_short_vector_spectrum": short_vectors,
            "low_degree_multisection_spectrum": multisection,
            "determinant_conductor_prospects": {
                "determinant": int(surface["determinant"]),
                "conductor_status": "UNKNOWN_NO_ARITHMETIC_SPECIALIZATION_SELECTED",
            },
            "automorphism_symmetry": symmetry,
            "moduli": {
                "complex_lattice_polarized_dimension": 1,
                "dimension_status": "EXACT_FROM_PICARD_RANK_19",
                "genus": None,
                "rationality": None,
                "genus_rationality_status": "UNKNOWN_NOT_COMPUTED",
            },
            "core_objectives_minimize": {
                "negative_maximum_generic_mw_rank": -maximum_mw,
                "easiest_known_source_mw_rank": source_record["mw_rank"],
                "source_support_count": source_record[
                    "reducible_fibre_support_count"
                ],
                "determinant": int(surface["determinant"]),
            },
        }
        record.update(record["core_objectives_minimize"])
        records.append(record)

    records.sort(key=lambda row: row["surface_id"])
    core_fields = (
        "negative_maximum_generic_mw_rank",
        "easiest_known_source_mw_rank",
        "source_support_count",
        "determinant",
    )
    core_layers = pareto_layers(records, core_fields)
    core_rank = {
        surface_id: rank
        for rank, layer in enumerate(core_layers, start=1)
        for surface_id in layer
    }
    for record in records:
        record["core_pareto_layer"] = core_rank[record["surface_id"]]

    pole_records = [
        row
        for row in records
        if row["easiest_known_source"]["minimum_nonzero_section_pole"].get(
            "pole_order"
        )
        is not None
    ]
    for row in pole_records:
        row["exact_minimum_section_pole_order"] = row["easiest_known_source"][
            "minimum_nonzero_section_pole"
        ]["pole_order"]
    pole_fields = core_fields + ("exact_minimum_section_pole_order",)
    pole_layers = pareto_layers(pole_records, pole_fields) if pole_records else []

    symmetry_records = [
        row
        for row in records
        if not row["automorphism_symmetry"]["lower_bound_is_identity_only"]
    ]
    for row in symmetry_records:
        row["negative_certified_symmetry_order_lower_bound"] = -row[
            "automorphism_symmetry"
        ]["certified_order_lower_bound"]
    symmetry_fields = core_fields + (
        "negative_certified_symmetry_order_lower_bound",
    )
    symmetry_layers = (
        pareto_layers(symmetry_records, symmetry_fields) if symmetry_records else []
    )

    route_records = [
        row
        for row in records
        if row["shortest_physical_neighbor_corridor"]["status"]
        == "PASS_EXACT_CERTIFIED_ROUTE"
    ]
    for row in route_records:
        cost = row["shortest_physical_neighbor_corridor"]["cost"]
        row["certified_corridor_maximum_old_fibre_degree"] = int(
            cost["maximum_old_fibre_degree"]
        )
        row["certified_corridor_edge_count"] = int(cost["edge_count"])
    route_fields = core_fields + (
        "certified_corridor_maximum_old_fibre_degree",
        "certified_corridor_edge_count",
    )
    route_layers = pareto_layers(route_records, route_fields) if route_records else []

    multisection_records = [
        row
        for row in records
        if row["low_degree_multisection_spectrum"]["status"]
        == "PASS_EXACT_D2_BOUNDED_D3_D4_EVIDENCE"
    ]
    multi_frame_surfaces = [
        row for row in records if row["catalogued_frame_count"] > 1
    ]

    assert len(records) == catalogue["accounting"]["surface_classes_by_T_NS"]
    assert len(records) == 161
    assert sum(
        row["easiest_known_source"]["source_kind"]
        == "external_rootful_source_inventory"
        for row in records
    ) == 48
    assert sum(
        row["easiest_known_source"]["source_kind"]
        == "easiest_currently_catalogued_frame"
        for row in records
    ) == 113
    assert len(multi_frame_surfaces) == 81
    assert len(pole_records) == 39
    assert len(route_records) == 0
    assert len(symmetry_records) == 117
    assert len(multisection_records) == 9
    assert len(core_layers) == 74
    assert core_layers[0] == [
        "K3-ebaf00b3723751ba",
        "K3-8188cdcda8c57b2d",
        "K3-f43753fb154e3406",
        "K3-14ad03cd7c1848b2",
    ]
    assert len(pole_layers[0]) == 8
    assert symmetry_layers[0] == [
        "K3-ebaf00b3723751ba",
        "K3-8188cdcda8c57b2d",
    ]

    return {
        "schema": "elkies-k3.rank7-surface-pareto.v1",
        "status": "PASS_EXACT_TYPED_SURFACE_PARETO_LEDGER_INCOMPLETE_METRICS_NOT_IMPUTED",
        "proof_scope": {
            "proved": (
                "All imported catalogue surfaces are ranked by universally available "
                "exact lattice metrics. External source certificates, exact minimum "
                "poles, certified routes, stabilizers, rootless short-vector data, and "
                "multisection spectra are attached only when their source artifacts "
                "provide them. Pareto dominance is exact for each declared metric set."
            ),
            "not_proved": (
                "The imported catalogue and source inventories remain bounded. Missing "
                "equations, fields, corridors, MW spectra, conductor prospects, and "
                "moduli genus/rationality are not inferred. A Pareto leader is a search "
                "priority, not an arithmetic-rank or optimality theorem."
            ),
        },
        "objective_policy": {
            "core_minimize": list(core_fields),
            "missing_data": (
                "Never impute. Use the core frontier for all surfaces and separate "
                "coverage-restricted enriched frontiers."
            ),
            "requested_metrics": [
                "maximum generic MW rank",
                "easiest source MW rank",
                "source support count",
                "source equation complexity",
                "source field of definition",
                "shortest physical neighbour corridor",
                "MW short-vector spectrum",
                "d=2,3,4 multisection spectrum",
                "determinant/conductor prospects",
                "automorphism/symmetry",
                "moduli genus/rationality",
            ],
        },
        "accounting": {
            "surfaces": len(records),
            "surfaces_with_external_low_mw_source_evidence": sum(
                row["easiest_known_source"]["source_kind"]
                == "external_rootful_source_inventory"
                for row in records
            ),
            "surfaces_using_catalogued_frame_as_source_fallback": sum(
                row["easiest_known_source"]["source_kind"]
                == "easiest_currently_catalogued_frame"
                for row in records
            ),
            "surfaces_with_multiple_catalogued_frames": len(multi_frame_surfaces),
            "surfaces_with_exact_minimum_pole": len(pole_records),
            "surfaces_with_certified_route": len(route_records),
            "surfaces_with_nontrivial_symmetry_evidence": len(symmetry_records),
            "surfaces_with_multisection_evidence": len(multisection_records),
            "surfaces_with_target_short_vector_evidence": sum(
                row["target_mw_short_vector_spectrum"]["status"]
                == "PASS_EXACT_AVAILABLE_TARGET_SHORT_VECTOR_EVIDENCE"
                for row in records
            ),
            "core_pareto_layers": len(core_layers),
            "core_pareto_frontier_size": len(core_layers[0]),
            "exact_pole_pareto_frontier_size": len(pole_layers[0]) if pole_layers else 0,
            "symmetry_evidence_pareto_frontier_size": (
                len(symmetry_layers[0]) if symmetry_layers else 0
            ),
            "certified_route_pareto_frontier_size": (
                len(route_layers[0]) if route_layers else 0
            ),
        },
        "frontiers": {
            "core_all_surfaces": {
                "fields": list(core_fields),
                "layers": core_layers,
            },
            "exact_minimum_pole_coverage": {
                "fields": list(pole_fields),
                "layers": pole_layers,
            },
            "symmetry_evidence_coverage": {
                "fields": list(symmetry_fields),
                "layers": symmetry_layers,
            },
            "certified_route_coverage": {
                "fields": list(route_fields),
                "layers": route_layers,
            },
        },
        "multi_frame_surface_ids": [
            row["surface_id"] for row in multi_frame_surfaces
        ],
        "surfaces": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=CATALOGUE)
    parser.add_argument("--source-ranking", type=Path, default=SOURCE_RANKING)
    parser.add_argument("--multisections", type=Path, default=MULTISECTIONS)
    parser.add_argument("--umbral", type=Path, default=UMBRAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    paths = {
        "catalogue": arguments.catalogue.resolve(),
        "source_ranking": arguments.source_ranking.resolve(),
        "multisections": arguments.multisections.resolve(),
        "umbral": arguments.umbral.resolve(),
    }
    result = build(
        json.loads(paths["catalogue"].read_text()),
        json.loads(paths["source_ranking"].read_text()),
        json.loads(paths["multisections"].read_text()),
        json.loads(paths["umbral"].read_text()),
    )
    result["inputs"] = {
        relative(path): digest(path) for path in paths.values()
    }
    result["reproduce"] = (
        "python3 elkies-k3/scripts/build_rank7_surface_pareto.py"
    )
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("rank-seven surface Pareto artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    accounting = result["accounting"]
    print(
        "RANK7PARETO|surfaces={}|core_front={}|pole={}|route={}|symmetry={}|status=PASS".format(
            accounting["surfaces"],
            accounting["core_pareto_frontier_size"],
            accounting["surfaces_with_exact_minimum_pole"],
            accounting["surfaces_with_certified_route"],
            accounting["surfaces_with_nontrivial_symmetry_evidence"],
        )
    )


if __name__ == "__main__":
    main()
