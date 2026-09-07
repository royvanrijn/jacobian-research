#!/usr/bin/env sage-python
"""Canonicalize the exact N(3D8) fixed-coordinate shell scan.

Quotient every accepted coordinate seed by the complete glue-preserving S3
residual group. Every row-module equality, stabilizer element, and canonical
representative is computed over Z by Hermite normal form. Compute literal
complement stabilizers, exact ternary gates, and (T,NS)-first
surface/auxiliary/frame deduplication.
"""

from __future__ import annotations

import argparse
import json
import runpy
from collections import Counter, defaultdict
from pathlib import Path

from sage.all import GF, ZZ, identity_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
RESIDUAL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-3d8-glue-residual-group-v1.json"
)
PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-3d8-fixed-coordinate-shell-probe-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-3d8-fixed-coordinate-shells-v1.json"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
SURFACE_COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "canonicalize_8a3_fixed_coordinate_shells.sage"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_fixed_seed_common")
SURFACE_COMMON = runpy.run_path(
    str(SURFACE_COMMON_SOURCE), run_name="_fixed_shell_surface_common"
)

rows = COMMON["rows"]
digest = COMMON["digest"]
matrix_key = COMMON["matrix_key"]
row_module_basis = COMMON["row_module_basis"]
row_module_key = COMMON["row_module_key"]
discriminant_invariants = COMMON["discriminant_invariants"]
induced_action = COMMON["induced_action"]
root_type = COMMON["root_type"]
ternary_gate = COMMON["ternary_gate"]
finite_order = SURFACE_COMMON["finite_order"]
rootless_intrinsics = SURFACE_COMMON["rootless_intrinsics"]
surface_rows = SURFACE_COMMON["surface_rows"]


def canonicalize(catalog, residual, probe, config=None):
    config = config or {
        "ambient_label": "3D8",
        "backend_id": "ROOTED-3D8",
        "residual_schema": "elkies-k3.3d8-glue-residual-group.v1",
        "residual_status": "PASS_EXACT_3D8_GLUE_AND_RESIDUAL_GROUP",
        "probe_schema": "elkies-k3.3d8-fixed-coordinate-shell-probe.v1",
        "probe_status": (
            "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_3D8_COORDINATE_SHELL_SCAN"
        ),
        "expected_group_order": 6,
        "expected_seed_count": 40,
        "orbit_id_prefix": "3D8-CF",
    }
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert residual["schema"] == config["residual_schema"]
    assert residual["status"] == config["residual_status"]
    assert probe["schema"] == config["probe_schema"]
    assert probe["status"] == config["probe_status"]
    assert probe["parameters"]["determinant_bound"] == 5000
    assert probe["parameters"]["minimum_mw_rank"] == 12

    ambient = matrix(
        ZZ,
        next(
            row["gram"]
            for row in catalog["rooted_niemeier_lattices"]
            if row["label"] == config["ambient_label"]
        ),
    )
    element_rows = residual["residual_group"]["elements"]
    assert residual["residual_group"]["order"] == len(element_rows) == config[
        "expected_group_order"
    ]
    actions = [matrix(ZZ, row["matrix"]) for row in element_rows]
    for action in actions:
        assert action * ambient * action.transpose() == ambient

    class_by_element_index = {}
    class_metadata = {}
    for class_row in residual["residual_group"]["conjugacy_classes"]:
        class_id = class_row["class_id"]
        class_metadata[class_id] = {
            "action_order": class_row["action_order"],
            "fixed_rank": class_row["fixed_rank"],
            "fixed_determinant": class_row["fixed_determinant"],
            "class_size": class_row["class_size"],
        }
        for element_index in class_row["element_indices_zero_based"]:
            assert element_index not in class_by_element_index
            class_by_element_index[element_index] = class_id
    assert len(class_by_element_index) == len(actions)

    group = []
    group_index_by_matrix_key = {}
    for index, (element, action) in enumerate(zip(element_rows, actions)):
        group.append(
            {
                "index_zero_based": index,
                "matrix": action,
                "class_id": class_by_element_index[index],
                "component_permutation_zero_based": element[
                    "component_permutation_zero_based"
                ],
                "component_diagram_variant_indices": element[
                    "component_diagram_variant_indices"
                ],
                "component_discriminant_maps": element.get(
                    "component_discriminant_maps"
                ),
                "order": element["action_order"],
                "fixed_rank": element["fixed_rank"],
            }
        )
        key = matrix_key(action)
        assert key not in group_index_by_matrix_key
        group_index_by_matrix_key[key] = index

    seeds = []
    for class_scan in probe["class_scans"]:
        class_id = class_scan["residual_conjugacy_class_id"]
        for seed in class_scan["accepted_seeds"]:
            basis = matrix(ZZ, seed["auxiliary_basis_in_ambient"])
            seeds.append(
                {
                    "source_class_id": class_id,
                    "coordinate_subset_zero_based": seed[
                        "coordinate_subset_zero_based"
                    ],
                    "auxiliary_basis": basis,
                    "exact_key": row_module_key(basis),
                }
            )
    assert len(seeds) == config["expected_seed_count"]

    exact_key_to_seed_indices = defaultdict(list)
    for index, seed in enumerate(seeds):
        exact_key_to_seed_indices[seed["exact_key"]].append(index)

    unassigned = set(range(len(seeds)))
    raw_orbits = []
    exact_images_checked = 0
    while unassigned:
        seed_index = min(unassigned)
        seed = seeds[seed_index]
        basis = seed["auxiliary_basis"]
        images = []
        for group_index, action in enumerate(actions):
            image_basis = row_module_basis(basis * action)
            images.append((row_module_key(image_basis), group_index, image_basis))
        exact_images_checked += len(images)
        canonical_key, transporter_index, canonical_basis = min(
            images, key=lambda row: (row[0], row[1])
        )
        member_seed_indices = set()
        stabilizer_indices = []
        for image_key, group_index, _image_basis in images:
            if image_key == seed["exact_key"]:
                stabilizer_indices.append(group_index)
            member_seed_indices.update(exact_key_to_seed_indices.get(image_key, []))
        assert seed_index in member_seed_indices
        assert stabilizer_indices
        assert len(actions) % len(stabilizer_indices) == 0
        assert member_seed_indices <= unassigned
        unassigned.difference_update(member_seed_indices)
        raw_orbits.append(
            {
                "canonical_key": canonical_key,
                "representative_basis": canonical_basis,
                "seed_indices": sorted(member_seed_indices),
                "source_stabilizer_indices": stabilizer_indices,
                "transporter_index": transporter_index,
                "orbit_size": len(actions) // len(stabilizer_indices),
            }
        )

    canonical_keys = [orbit["canonical_key"] for orbit in raw_orbits]
    assert len(set(canonical_keys)) == len(canonical_keys)
    raw_orbits.sort(key=lambda orbit: orbit["canonical_key"])

    genera_cache = {}
    form_cache = {}
    orbit_records = []
    for orbit_index, orbit in enumerate(raw_orbits):
        auxiliary_basis = orbit["representative_basis"]
        canonical_key = orbit["canonical_key"]
        assert row_module_key(auxiliary_basis) == canonical_key
        auxiliary = auxiliary_basis * ambient * auxiliary_basis.transpose()
        complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
        frame = complement_basis * ambient * complement_basis.transpose()
        assert frame.det() == auxiliary.det()
        roots = root_type(frame)
        origin_class_ids = sorted(
            {seeds[index]["source_class_id"] for index in orbit["seed_indices"]}
        )

        transporter = actions[orbit["transporter_index"]]
        inverse_transporter = transporter.inverse().change_ring(ZZ)
        canonical_stabilizer_indices = []
        for source_stabilizer_index in orbit["source_stabilizer_indices"]:
            conjugate = (
                inverse_transporter
                * actions[source_stabilizer_index]
                * transporter
            )
            canonical_stabilizer_indices.append(
                group_index_by_matrix_key[matrix_key(conjugate)]
            )
        canonical_stabilizer_indices = sorted(set(canonical_stabilizer_indices))
        assert len(canonical_stabilizer_indices) == len(
            orbit["source_stabilizer_indices"]
        )

        stabilizer = []
        for group_index in canonical_stabilizer_indices:
            group_row = group[group_index]
            assert row_module_key(auxiliary_basis * group_row["matrix"]) == canonical_key
            action = induced_action(complement_basis, group_row["matrix"])
            assert action * frame * action.transpose() == frame
            stabilizer.append(
                {
                    "residual_group_index_zero_based": group_index,
                    "residual_conjugacy_class_id": group_row["class_id"],
                    "component_permutation_zero_based": group_row[
                        "component_permutation_zero_based"
                    ],
                    "component_diagram_variant_indices": group_row[
                        "component_diagram_variant_indices"
                    ],
                    "component_discriminant_maps": group_row[
                        "component_discriminant_maps"
                    ],
                    "order": finite_order(action),
                    "moved_dimension_mod_2": int(
                        matrix(GF(2), action - identity_matrix(ZZ, 17)).rank()
                    ),
                    "fixed_dimension_over_Q": int(
                        17 - (action - identity_matrix(ZZ, 17)).rank()
                    ),
                }
            )
        assert any(
            row["residual_conjugacy_class_id"] in origin_class_ids
            and row["moved_dimension_mod_2"] > 0
            for row in stabilizer
        )
        gate = ternary_gate(frame, genera_cache, form_cache)
        orbit_records.append(
            {
                "orbit_id": f"{config['orbit_id_prefix']}-O{orbit_index + 1:04d}",
                "residual_group_orbit_size": orbit["orbit_size"],
                "source_seed_indices_zero_based": orbit["seed_indices"],
                "source_coordinate_subsets": [
                    {
                        "residual_conjugacy_class_id": seeds[index][
                            "source_class_id"
                        ],
                        "coordinate_subset_zero_based": seeds[index][
                            "coordinate_subset_zero_based"
                        ],
                    }
                    for index in orbit["seed_indices"]
                ],
                "originating_residual_conjugacy_class_ids": origin_class_ids,
                "determinant": int(auxiliary.det()),
                "discriminant_invariants_greater_than_one": (
                    discriminant_invariants(auxiliary)
                ),
                "auxiliary_basis_in_ambient": rows(auxiliary_basis),
                "auxiliary_gram": rows(auxiliary),
                "complement_basis_in_ambient": rows(complement_basis),
                "frame_gram": rows(frame),
                "root_data": roots,
                "rootless_intrinsics": rootless_intrinsics(frame, roots),
                "literal_residual_stabilizer": stabilizer,
                "k3_discriminant_gate": gate,
            }
        )

    surfaces = surface_rows(config["ambient_label"], orbit_records)
    compatible_orbit_ids = {
        orbit_id
        for surface in surfaces
        for orbit_id in surface["residual_group_orbit_ids"]
    }
    orbit_mw = Counter(
        row["root_data"]["mw_rank_for_rho_19"] for row in orbit_records
    )
    frame_mw = Counter(
        frame["root_data"]["mw_rank_for_rho_19"]
        for surface in surfaces
        for frame in surface["frames"]
    )
    orbit_sizes = Counter(
        row["residual_group_orbit_size"] for row in orbit_records
    )
    stabilizer_class_coverage = Counter(
        class_id
        for orbit in orbit_records
        for class_id in {
            action["residual_conjugacy_class_id"]
            for action in orbit["literal_residual_stabilizer"]
            if action["moved_dimension_mod_2"] > 0
        }
    )
    return {
        "ambient_label": config["ambient_label"],
        "backend_id": config["backend_id"],
        "residual_group": {
            "order": len(group),
            "component_permutation_image_order": residual["residual_group"][
                "component_permutation_image_order"
            ],
            "component_kernel_order": residual["residual_group"][
                "component_kernel_order"
            ],
            "conjugacy_classes": class_metadata,
        },
        "canonicalization": {
            "exact_group_images_checked": exact_images_checked,
            "equality_certificate": (
                "every residual image is reduced to its exact primitive "
                "integral row-module HNF; equality is literal HNF equality"
            ),
        },
        "source_probe_accounting": probe["accounting"],
        "accounting": {
            "high_mw_mod2_accepted_seeds_before_residual_dedup": len(seeds),
            "residual_group_embedding_orbits": len(orbit_records),
            "k3_compatible_residual_group_embedding_orbits": len(
                compatible_orbit_ids
            ),
            "residual_group_orbit_size_distribution": {
                str(key): value for key, value in sorted(orbit_sizes.items())
            },
            "orbit_mw_rank_distribution": {
                str(key): value for key, value in sorted(orbit_mw.items())
            },
            "distinct_frame_discriminant_forms": len(form_cache),
            "forms_with_matching_ternary_genus": sum(
                gate["matching_even_ternary_genera"] > 0
                for gate in form_cache.values()
            ),
            "surface_classes_after_T_NS_first_dedup": len(surfaces),
            "partner_auxiliary_isometry_classes_after_surface_dedup": sum(
                len(surface["partner_auxiliaries"]) for surface in surfaces
            ),
            "frame_isometry_classes_after_surface_dedup": sum(
                len(surface["frames"]) for surface in surfaces
            ),
            "post_dedup_frame_mw_rank_distribution": {
                str(key): value for key, value in sorted(frame_mw.items())
            },
            "nontrivial_mod2_stabilizer_class_coverage": {
                key: value
                for key, value in sorted(stabilizer_class_coverage.items())
            },
        },
        "embedding_orbits": orbit_records,
        "surfaces_T_NS_first": surfaces,
    }


def build(catalog, residual, probe):
    backend = canonicalize(catalog, residual, probe)
    assert (
        backend["accounting"]["residual_group_embedding_orbits"],
        backend["accounting"][
            "k3_compatible_residual_group_embedding_orbits"
        ],
        backend["accounting"]["surface_classes_after_T_NS_first_dedup"],
        backend["accounting"][
            "partner_auxiliary_isometry_classes_after_surface_dedup"
        ],
        backend["accounting"][
            "frame_isometry_classes_after_surface_dedup"
        ],
    ) == (40, 25, 7, 7, 7)
    assert backend["accounting"]["post_dedup_frame_mw_rank_distribution"] == {
        "12": 6,
        "13": 1,
    }
    return {
        "schema": "elkies-k3.3d8-fixed-coordinate-shells.v1",
        "status": "PASS_EXACT_DECLARED_3D8_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        "proof_scope": {
            "proved": (
                "complete order-six residual lift group; every 7-of-r coordinate "
                "summand in both nonidentity fixed-class LLL bases; exact "
                "determinant, length, MW12--17, mod-two, residual quotient, "
                "ternary, and T/NS-first gates"
            ),
            "not_proved": (
                "all primitive rank-seven fixed-lattice or ambient sublattices, "
                "full Weyl embedding orbits, ternary class enumeration beyond "
                "genera, or determinant-band completeness"
            ),
        },
        "parameters": probe["parameters"],
        "backends": [backend],
        "accounting": {
            "backends": 1,
            "coordinate_subsets_tested": backend["source_probe_accounting"][
                "coordinate_subsets_tested"
            ],
            "accepted_seeds_before_residual_dedup": backend["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
            "residual_group_embedding_orbits": backend["accounting"][
                "residual_group_embedding_orbits"
            ],
            "k3_compatible_residual_group_embedding_orbits": backend[
                "accounting"
            ]["k3_compatible_residual_group_embedding_orbits"],
            "surface_classes_before_global_cross_backend_dedup": backend[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
            "partner_auxiliary_classes_before_global_cross_backend_dedup": (
                backend["accounting"][
                    "partner_auxiliary_isometry_classes_after_surface_dedup"
                ]
            ),
            "frame_classes_before_global_cross_backend_dedup": backend[
                "accounting"
            ]["frame_isometry_classes_after_surface_dedup"],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--residual", type=Path, default=RESIDUAL)
    parser.add_argument("--probe", type=Path, default=PROBE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build(
        json.loads(arguments.catalog.read_text()),
        json.loads(arguments.residual.read_text()),
        json.loads(arguments.probe.read_text()),
    )
    payload["inputs"] = {
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(arguments.catalog),
        str(arguments.residual.resolve().relative_to(ROOT)): digest(arguments.residual),
        str(arguments.probe.resolve().relative_to(ROOT)): digest(arguments.probe),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
        str(SURFACE_COMMON_SOURCE.resolve().relative_to(ROOT)): digest(
            SURFACE_COMMON_SOURCE
        ),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/canonicalize_3d8_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("3D8 canonical shell artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "3D8CANON|seeds={}|orbits={}|surfaces={}|partners={}|frames={}|status=PASS_EXACT".format(
            payload["accounting"]["accepted_seeds_before_residual_dedup"],
            payload["accounting"]["residual_group_embedding_orbits"],
            payload["accounting"][
                "surface_classes_before_global_cross_backend_dedup"
            ],
            payload["accounting"][
                "partner_auxiliary_classes_before_global_cross_backend_dedup"
            ],
            payload["accounting"][
                "frame_classes_before_global_cross_backend_dedup"
            ],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
