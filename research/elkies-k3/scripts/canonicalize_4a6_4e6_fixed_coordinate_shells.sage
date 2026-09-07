#!/usr/bin/env sage-python
"""Canonicalize the exact 4A6/4E6 component-fixed coordinate shells.

Consume the complete residual sections and exhaustive pre-quotient coordinate
scan.  Close every accepted rank-seven auxiliary under the complete chamber
residual group, compute literal stabilizers and their induced mod-two
complement actions, apply the exact ternary discriminant-form gate, and
deduplicate first by (T,NS), then by auxiliary and frame isometry.

The result is exact for the declared fixed-lattice coordinate languages and
residual quotients.  It is not a full fixed-lattice, Weyl-orbit, or
determinant-band enumeration.
"""

from __future__ import annotations

import argparse
import json
import runpy
from collections import Counter
from pathlib import Path

from sage.all import GF, ZZ, identity_matrix, matrix, pari


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
SECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-4a6-4e6-residual-sections-v1.json"
)
PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-4a6-4e6-fixed-coordinate-shell-probe-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-4a6-4e6-fixed-coordinate-shells-v1.json"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_fixed_seed_common")

rows = COMMON["rows"]
digest = COMMON["digest"]
gram_digest = COMMON["gram_digest"]
compact_digest = COMMON["compact_digest"]
matrix_key = COMMON["matrix_key"]
row_module_basis = COMMON["row_module_basis"]
row_module_key = COMMON["row_module_key"]
negative_form_key = COMMON["negative_form_key"]
discriminant_invariants = COMMON["discriminant_invariants"]
induced_action = COMMON["induced_action"]
root_type = COMMON["root_type"]
ternary_gate = COMMON["ternary_gate"]
find_or_add_isometry_class = COMMON["find_or_add_isometry_class"]


def finite_order(action):
    identity = identity_matrix(ZZ, action.nrows())
    power = identity
    for order in range(1, 25):
        power *= action
        if power == identity:
            return order
    raise AssertionError("unexpected finite action order")


def rootless_intrinsics(frame, roots):
    if roots["root_rank"] != 0:
        return None
    data = pari(frame).qfminim(4)
    assert int(data[1]) == 4
    return {
        "minimum_squared_norm": 4,
        "norm_four_vectors": int(data[0]),
        "norm_four_unoriented_pairs": int(data[0]) // 2,
    }


def surface_rows(label, orbit_records):
    surfaces_by_key = {}
    for orbit_index, orbit in enumerate(orbit_records):
        ns_key = negative_form_key(
            orbit["k3_discriminant_gate"][
                "frame_discriminant_form_normal_key"
            ]
        )
        for ternary in orbit["k3_discriminant_gate"][
            "ternary_genus_representatives"
        ]:
            key = {
                "ns_discriminant_form_key": ns_key,
                "transcendental_gram": ternary["gram"],
            }
            encoded = json.dumps(key, sort_keys=True, separators=(",", ":"))
            if encoded not in surfaces_by_key:
                surfaces_by_key[encoded] = {
                    "surface_id": f"K3-{label}CF-{compact_digest(key)}",
                    "surface_key": key,
                    "determinant": orbit["determinant"],
                    "orbit_indices": [],
                    "partner_classes": [],
                    "frame_classes": [],
                }
            surface = surfaces_by_key[encoded]
            surface["orbit_indices"].append(orbit_index)
            partner = find_or_add_isometry_class(
                surface["partner_classes"],
                matrix(ZZ, orbit["auxiliary_gram"]),
                orbit_index,
            )
            frame = find_or_add_isometry_class(
                surface["frame_classes"],
                matrix(ZZ, orbit["frame_gram"]),
                orbit_index,
            )
            frame.setdefault("partner_class_objects", [])
            if partner not in frame["partner_class_objects"]:
                frame["partner_class_objects"].append(partner)

    result = []
    for encoded in sorted(surfaces_by_key):
        surface = surfaces_by_key[encoded]
        partner_id_by_object = {}
        partners = []
        for index, partner in enumerate(surface["partner_classes"], start=1):
            partner_id = f"{surface['surface_id']}-K{index:03d}"
            partner_id_by_object[id(partner)] = partner_id
            representative = orbit_records[partner["record_indices"][0]]
            partners.append(
                {
                    "partner_id": partner_id,
                    "gram": rows(partner["gram"]),
                    "gram_sha256": gram_digest(partner["gram"]),
                    "determinant": int(partner["gram"].det()),
                    "residual_group_orbit_ids": [
                        orbit_records[value]["orbit_id"]
                        for value in partner["record_indices"]
                    ],
                    "representative_auxiliary_basis_in_ambient": representative[
                        "auxiliary_basis_in_ambient"
                    ],
                }
            )
        frames = []
        for index, frame_class in enumerate(
            surface["frame_classes"], start=1
        ):
            representative = orbit_records[frame_class["record_indices"][0]]
            frames.append(
                {
                    "frame_id": f"{surface['surface_id']}-F{index:03d}",
                    "gram": rows(frame_class["gram"]),
                    "gram_sha256": gram_digest(frame_class["gram"]),
                    "determinant": int(frame_class["gram"].det()),
                    "root_data": representative["root_data"],
                    "rootless_intrinsics": representative[
                        "rootless_intrinsics"
                    ],
                    "partner_ids": sorted(
                        partner_id_by_object[id(partner)]
                        for partner in frame_class[
                            "partner_class_objects"
                        ]
                    ),
                    "residual_group_orbit_ids": [
                        orbit_records[value]["orbit_id"]
                        for value in frame_class["record_indices"]
                    ],
                    "representative_complement_basis_in_ambient": representative[
                        "complement_basis_in_ambient"
                    ],
                    "representative_literal_residual_stabilizer": representative[
                        "literal_residual_stabilizer"
                    ],
                }
            )
        result.append(
            {
                "surface_id": surface["surface_id"],
                "surface_key": surface["surface_key"],
                "determinant": surface["determinant"],
                "partner_auxiliaries": partners,
                "frames": frames,
                "residual_group_orbit_ids": [
                    orbit_records[value]["orbit_id"]
                    for value in surface["orbit_indices"]
                ],
            }
        )
    result.sort(key=lambda row: (row["determinant"], row["surface_id"]))
    return result


def canonicalize_backend(ambient, section_backend, probe_backend):
    label = section_backend["ambient_label"]
    assert probe_backend["ambient_label"] == label
    assert section_backend["residual_group_order"] == len(
        section_backend["elements"]
    )
    group = []
    for index, element in enumerate(section_backend["elements"]):
        action = matrix(ZZ, element["matrix"])
        assert action * ambient * action.transpose() == ambient
        group.append(
            {
                "index_zero_based": index,
                "matrix": action,
                "component_permutation_zero_based": element[
                    "component_permutation_zero_based"
                ],
                "order": element["action_order"],
                "fixed_rank": element["fixed_rank"],
            }
        )
    class_by_matrix_key = {}
    class_metadata = {}
    for class_row in section_backend["conjugacy_classes"]:
        class_id = class_row["class_id"]
        class_metadata[class_id] = {
            "component_cycle_type": class_row["component_cycle_type"],
            "action_order": class_row["action_order"],
            "fixed_rank": class_row["fixed_rank"],
            "fixed_determinant": class_row["fixed_determinant"],
            "class_size": class_row["class_size"],
        }
        for key in class_row["element_matrix_keys"]:
            matrix_tuple = tuple(key)
            assert matrix_tuple not in class_by_matrix_key
            class_by_matrix_key[matrix_tuple] = class_id
    assert len(class_by_matrix_key) == len(group)
    for group_row in group:
        group_row["class_id"] = class_by_matrix_key[
            matrix_key(group_row["matrix"])
        ]

    seeds = []
    for class_scan in probe_backend["class_scans"]:
        class_id = class_scan["residual_conjugacy_class_id"]
        for seed in class_scan["accepted_seeds"]:
            seeds.append(
                {
                    "source_class_id": class_id,
                    "coordinate_subset_zero_based": seed[
                        "coordinate_subset_zero_based"
                    ],
                    "auxiliary_basis": matrix(
                        ZZ, seed["auxiliary_basis_in_ambient"]
                    ),
                }
            )
    assert len(seeds) == probe_backend["accounting"][
        "high_mw_mod2_accepted_seeds_before_residual_dedup"
    ]

    orbit_by_key = {}
    for seed_index, seed in enumerate(seeds):
        images = {}
        for group_row in group:
            image = row_module_basis(
                seed["auxiliary_basis"] * group_row["matrix"]
            )
            images[row_module_key(image)] = image
        canonical_key = min(images)
        if canonical_key not in orbit_by_key:
            orbit_by_key[canonical_key] = {
                "representative_basis": images[canonical_key],
                "image_bases": images,
                "seed_indices": [],
            }
        orbit = orbit_by_key[canonical_key]
        assert set(orbit["image_bases"]) == set(images)
        orbit["seed_indices"].append(seed_index)

    genera_cache = {}
    form_cache = {}
    orbit_records = []
    for orbit_index, canonical_key in enumerate(sorted(orbit_by_key)):
        orbit = orbit_by_key[canonical_key]
        auxiliary_basis = orbit["representative_basis"]
        auxiliary = auxiliary_basis * ambient * auxiliary_basis.transpose()
        complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
        frame = complement_basis * ambient * complement_basis.transpose()
        assert frame.det() == auxiliary.det()
        roots = root_type(frame)
        origin_class_ids = sorted(
            {
                seeds[index]["source_class_id"]
                for index in orbit["seed_indices"]
            }
        )
        stabilizer = []
        for group_row in group:
            if (
                row_module_key(auxiliary_basis * group_row["matrix"])
                != canonical_key
            ):
                continue
            action = induced_action(complement_basis, group_row["matrix"])
            assert action * frame * action.transpose() == frame
            stabilizer.append(
                {
                    "section_index_zero_based": group_row[
                        "index_zero_based"
                    ],
                    "residual_conjugacy_class_id": group_row["class_id"],
                    "component_permutation_zero_based": group_row[
                        "component_permutation_zero_based"
                    ],
                    "order": finite_order(action),
                    "moved_dimension_mod_2": int(
                        matrix(
                            GF(2),
                            action - identity_matrix(ZZ, 17),
                        ).rank()
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
                "orbit_id": f"{label}-CF-O{orbit_index + 1:04d}",
                "residual_group_orbit_size": len(orbit["image_bases"]),
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

    surfaces = surface_rows(label, orbit_records)
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
    k3_compatible_orbit_ids = {
        orbit_id
        for surface in surfaces
        for orbit_id in surface["residual_group_orbit_ids"]
    }
    return {
        "ambient_label": label,
        "backend_id": f"ROOTED-{label}",
        "residual_group": {
            "order": len(group),
            "component_permutation_image_order": section_backend[
                "component_permutation_image_order"
            ],
            "conjugacy_classes": class_metadata,
        },
        "source_probe_accounting": probe_backend["accounting"],
        "accounting": {
            "high_mw_mod2_accepted_seeds_before_residual_dedup": len(seeds),
            "residual_group_embedding_orbits": len(orbit_records),
            "k3_compatible_residual_group_embedding_orbits": len(
                k3_compatible_orbit_ids
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
                len(surface["partner_auxiliaries"])
                for surface in surfaces
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


def build(catalog, sections, probe):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert sections["schema"] == "elkies-k3.4a6-4e6-residual-sections.v1"
    assert sections["status"] == "PASS_EXACT_RESIDUAL_SECTIONS"
    assert probe["schema"] == (
        "elkies-k3.4a6-4e6-fixed-coordinate-shell-probe.v1"
    )
    assert probe["status"] == (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_COORDINATE_SHELL_SCAN"
    )
    assert probe["parameters"]["determinant_bound"] == 5000
    assert probe["parameters"]["minimum_mw_rank"] == 12
    ambient_by_label = {
        row["label"]: matrix(ZZ, row["gram"])
        for row in catalog["rooted_niemeier_lattices"]
    }
    section_by_label = {
        row["ambient_label"]: row for row in sections["backends"]
    }
    probe_by_label = {
        row["ambient_label"]: row for row in probe["backends"]
    }
    backends = [
        canonicalize_backend(
            ambient_by_label[label],
            section_by_label[label],
            probe_by_label[label],
        )
        for label in ("4A6", "4E6")
    ]
    assert [
        (
            row["ambient_label"],
            row["accounting"]["residual_group_embedding_orbits"],
            row["accounting"][
                "k3_compatible_residual_group_embedding_orbits"
            ],
            row["accounting"]["surface_classes_after_T_NS_first_dedup"],
            row["accounting"][
                "partner_auxiliary_isometry_classes_after_surface_dedup"
            ],
            row["accounting"][
                "frame_isometry_classes_after_surface_dedup"
            ],
        )
        for row in backends
    ] == [("4A6", 86, 71, 9, 10, 10), ("4E6", 45, 42, 1, 1, 1)]
    return {
        "schema": "elkies-k3.4a6-4e6-fixed-coordinate-shells.v1",
        "status": "PASS_EXACT_DECLARED_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        "proof_scope": {
            "proved": (
                "complete residual lift groups; every 7-of-r coordinate "
                "summand for every nonidentity residual matrix conjugacy "
                "class of fixed rank at least seven; exact determinant, "
                "length, MW12--17, mod-two, full residual quotient, ternary, "
                "and T/NS-first gates"
            ),
            "not_proved": (
                "all primitive rank-seven fixed-lattice or ambient "
                "sublattices, full Weyl embedding orbits, ternary class "
                "enumeration beyond genera, or determinant-band completeness"
            ),
        },
        "parameters": probe["parameters"],
        "backends": backends,
        "accounting": {
            "backends": len(backends),
            "coordinate_subsets_tested": sum(
                row["source_probe_accounting"]["coordinate_subsets_tested"]
                for row in backends
            ),
            "accepted_seeds_before_residual_dedup": sum(
                row["accounting"][
                    "high_mw_mod2_accepted_seeds_before_residual_dedup"
                ]
                for row in backends
            ),
            "residual_group_embedding_orbits": sum(
                row["accounting"]["residual_group_embedding_orbits"]
                for row in backends
            ),
            "k3_compatible_residual_group_embedding_orbits": sum(
                row["accounting"][
                    "k3_compatible_residual_group_embedding_orbits"
                ]
                for row in backends
            ),
            "surface_classes_before_global_cross_backend_dedup": sum(
                row["accounting"]["surface_classes_after_T_NS_first_dedup"]
                for row in backends
            ),
            "partner_auxiliary_classes_before_global_cross_backend_dedup": sum(
                row["accounting"][
                    "partner_auxiliary_isometry_classes_after_surface_dedup"
                ]
                for row in backends
            ),
            "frame_classes_before_global_cross_backend_dedup": sum(
                row["accounting"][
                    "frame_isometry_classes_after_surface_dedup"
                ]
                for row in backends
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--sections", type=Path, default=SECTIONS)
    parser.add_argument("--probe", type=Path, default=PROBE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build(
        json.loads(arguments.catalog.read_text()),
        json.loads(arguments.sections.read_text()),
        json.loads(arguments.probe.read_text()),
    )
    payload["inputs"] = {
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(
            arguments.catalog
        ),
        str(arguments.sections.resolve().relative_to(ROOT)): digest(
            arguments.sections
        ),
        str(arguments.probe.resolve().relative_to(ROOT)): digest(
            arguments.probe
        ),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/canonicalize_4a6_4e6_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("4A6/4E6 canonical shell artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "A6E6CANON|seeds={}|orbits={}|surfaces={}|partners={}|frames={}|status=PASS_EXACT".format(
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
