#!/usr/bin/env sage-python
"""Canonicalize the exact N(12A2) fixed-coordinate shell scan.

Quotient every accepted coordinate seed by the complete order-190080 ternary
Golay monomial group.  The action becomes a permutation of the 24 simple-root
coordinates.  Batched row reduction over F_2 is used only as a rejection
filter; every possible seed equality, stabilizer element, and canonical
representative in the selected modular bucket is verified by integral HNF.

Then compute the full literal residual stabilizer on each complement, apply
the exact ternary gate, and deduplicate first by (T,NS), then by auxiliary and
frame isometry.
"""

from __future__ import annotations

import argparse
import json
import runpy
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sage.all import GF, QQ, ZZ, Permutation, PermutationGroup, identity_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
RESIDUAL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-12a2-ternary-golay-residual-group-v1.json"
)
PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-12a2-fixed-coordinate-shell-probe-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-12a2-fixed-coordinate-shells-v1.json"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
SURFACE_COMMON_SOURCE = (
    Path(__file__).resolve().parent / "canonicalize_8a3_fixed_coordinate_shells.sage"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_fixed_seed_common")
SURFACE_COMMON = runpy.run_path(
    str(SURFACE_COMMON_SOURCE), run_name="_12a2_fixed_shell_surface_common"
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


def root_coordinate_permutation(root_basis, ambient_action):
    root_action = root_basis * ambient_action * root_basis.inverse()
    assert all(value in (0, 1) for value in root_action.list())
    images = []
    for source in range(24):
        targets = [target for target in range(24) if root_action[source, target]]
        assert len(targets) == 1
        images.append(targets[0] + 1)
    assert sorted(images) == list(range(1, 25))
    return Permutation(images)


def batched_rref_mod_2_packed(images):
    """Return exact F_2 RREF rows packed as seven 24-bit integers."""
    values = np.remainder(images, 2).astype(np.uint32)
    weights = np.uint32(1) << np.arange(24, dtype=np.uint32)
    packed = np.sum(values * weights, axis=2, dtype=np.uint32)
    count, rank = packed.shape
    assert rank == 7
    for pivot_row in range(rank):
        waiting = np.ones(count, dtype=bool)
        for column in range(24):
            indices = np.flatnonzero(waiting)
            block = (
                packed[indices, pivot_row:] >> np.uint32(column)
            ) & np.uint32(1)
            has_pivot = np.any(block != 0, axis=1)
            selected = indices[has_pivot]
            if not len(selected):
                continue
            source_rows = np.argmax(block[has_pivot] != 0, axis=1) + pivot_row
            temporary = packed[selected, pivot_row].copy()
            packed[selected, pivot_row] = packed[selected, source_rows]
            packed[selected, source_rows] = temporary
            for other_row in range(rank):
                if other_row == pivot_row:
                    continue
                nonzero = (
                    (packed[selected, other_row] >> np.uint32(column))
                    & np.uint32(1)
                ).astype(bool)
                targets = selected[nonzero]
                packed[targets, other_row] ^= packed[targets, pivot_row]
            waiting[selected] = False
            if not np.any(waiting):
                break
        assert not np.any(waiting), "primitive auxiliary lost rank modulo 2"
    return packed


def packed_key(value):
    return tuple(map(int, value))


def validate_modular_engine(scaled_root_basis, inverse_images):
    sample_indices = np.linspace(
        0, len(inverse_images) - 1, num=32, dtype=np.int64
    )
    source = np.asarray(scaled_root_basis, dtype=np.int16)
    images = np.transpose(source[:, inverse_images[sample_indices]], (1, 0, 2))
    packed = batched_rref_mod_2_packed(images)
    field = GF(2)
    for row, image in zip(packed, images):
        reference = matrix(field, image).echelon_form()
        reference_packed = tuple(
            sum(int(value) << column for column, value in enumerate(values))
            for values in reference.rows()
        )
        assert packed_key(row) == reference_packed


def minimum_packed_bucket(packed):
    selected = np.ones(len(packed), dtype=bool)
    for column in range(packed.shape[1]):
        minimum = packed[selected, column].min()
        selected &= packed[:, column] == minimum
    return np.flatnonzero(selected)


def permutation_key(element):
    return tuple(int(element(index)) for index in range(1, 25))


def component_data(element):
    permutation = []
    variants = []
    for source in range(12):
        targets = {
            (int(element(2 * source + local + 1)) - 1) // 2
            for local in range(2)
        }
        assert len(targets) == 1
        target = targets.pop()
        permutation.append(target)
        variants.append(
            [
                (int(element(2 * source + local + 1)) - 1) % 2
                for local in range(2)
            ]
        )
    assert sorted(permutation) == list(range(12))
    assert all(sorted(variant) == [0, 1] for variant in variants)
    return permutation, variants


def canonicalize(catalog, residual, probe):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert residual["schema"] == "elkies-k3.12a2-ternary-golay-residual-group.v1"
    assert residual["status"] == (
        "PASS_EXACT_12A2_TERNARY_GOLAY_AND_FULL_RESIDUAL_GROUP"
    )
    assert probe["schema"] == "elkies-k3.12a2-fixed-coordinate-shell-probe.v1"
    assert probe["status"] == (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_12A2_COORDINATE_SHELL_SCAN"
    )
    assert probe["parameters"]["determinant_bound"] == 5000
    assert probe["parameters"]["minimum_mw_rank"] == 12

    ambient = matrix(
        ZZ,
        next(
            row["gram"]
            for row in catalog["rooted_niemeier_lattices"]
            if row["label"] == "12A2"
        ),
    )
    root_basis = matrix(ZZ, residual["root_lattice"]["combined_basis_in_ambient"])
    root_basis_inverse = root_basis.inverse()
    assert abs(root_basis.det()) == 729
    generators = [
        root_coordinate_permutation(root_basis, matrix(ZZ, row["matrix"]))
        for row in residual["residual_group"]["generators"]
    ]
    permutation_group = PermutationGroup(generators)
    assert int(permutation_group.order()) == 190080
    elements = list(permutation_group)
    assert len(elements) == 190080
    element_key_to_index = {
        permutation_key(element): index for index, element in enumerate(elements)
    }
    assert len(element_key_to_index) == len(elements)
    inverse_images = np.asarray(
        [
            [int(element.inverse()(index + 1)) - 1 for index in range(24)]
            for element in elements
        ],
        dtype=np.uint8,
    )

    class_metadata = {}
    class_by_element_key = {}
    for class_row in residual["residual_group"]["conjugacy_classes"]:
        class_id = class_row["class_id"]
        class_metadata[class_id] = {
            "action_order": class_row["action_order"],
            "fixed_rank": class_row["fixed_rank"],
            "fixed_determinant": class_row["fixed_determinant"],
            "class_size": class_row["class_size"],
        }
        representative = permutation_group(
            root_coordinate_permutation(
                root_basis, matrix(ZZ, class_row["representative_matrix"])
            )
        )
        conjugacy_class = permutation_group.conjugacy_class(representative)
        assert int(conjugacy_class.cardinality()) == class_row["class_size"]
        for element in conjugacy_class:
            key = permutation_key(element)
            assert key not in class_by_element_key
            class_by_element_key[key] = class_id
    assert len(class_by_element_key) == len(elements)

    seeds = []
    for class_scan in probe["class_scans"]:
        class_id = class_scan["residual_conjugacy_class_id"]
        for seed in class_scan["accepted_seeds"]:
            basis = matrix(ZZ, seed["auxiliary_basis_in_ambient"])
            scaled_root = 3 * basis * root_basis_inverse
            assert all(value.denominator() == 1 for value in scaled_root.list())
            scaled_root = matrix(ZZ, scaled_root)
            modular = batched_rref_mod_2_packed(
                np.asarray(scaled_root, dtype=np.int16)[None, :, :]
            )[0]
            seeds.append(
                {
                    "source_class_id": class_id,
                    "coordinate_subset_zero_based": seed[
                        "coordinate_subset_zero_based"
                    ],
                    "auxiliary_basis": basis,
                    "scaled_root_basis": scaled_root,
                    "exact_key": row_module_key(scaled_root),
                    "modular_key": packed_key(modular),
                }
            )
    assert len(seeds) == 237
    validate_modular_engine(seeds[0]["scaled_root_basis"], inverse_images)

    exact_key_to_seed_indices = defaultdict(list)
    modular_key_to_exact_keys = defaultdict(set)
    for index, seed in enumerate(seeds):
        exact_key_to_seed_indices[seed["exact_key"]].append(index)
        modular_key_to_exact_keys[seed["modular_key"]].add(seed["exact_key"])
    modular_target_array = np.asarray(
        sorted(modular_key_to_exact_keys), dtype=np.uint32
    )
    modular_target_void = modular_target_array.view(
        np.dtype((np.void, modular_target_array.dtype.itemsize * 7))
    ).reshape(-1)

    action_cache = {}

    def ambient_action(group_index):
        if group_index not in action_cache:
            element = elements[group_index]
            root_action = matrix(ZZ, 24)
            for source in range(24):
                root_action[source, int(element(source + 1)) - 1] = 1
            rational_action = root_basis_inverse * root_action * root_basis
            assert all(value.denominator() == 1 for value in rational_action.list())
            action = matrix(ZZ, rational_action)
            assert abs(action.det()) == 1
            assert action * ambient * action.transpose() == ambient
            action_cache[group_index] = action
        return action_cache[group_index]

    unassigned = set(range(len(seeds)))
    raw_orbits = []
    modular_images_checked = 0
    exact_candidate_images = 0
    while unassigned:
        seed_index = min(unassigned)
        seed = seeds[seed_index]
        scaled_root = np.asarray(seed["scaled_root_basis"], dtype=np.int16)
        images = np.transpose(scaled_root[:, inverse_images], (1, 0, 2))
        modular_images = batched_rref_mod_2_packed(images)
        modular_images_checked += len(modular_images)
        modular_void = modular_images.view(
            np.dtype((np.void, modular_images.dtype.itemsize * 7))
        ).reshape(-1)

        exact_image_cache = {}

        def exact_image(group_index):
            nonlocal exact_candidate_images
            if group_index not in exact_image_cache:
                value = row_module_basis(matrix(ZZ, images[group_index]))
                exact_image_cache[group_index] = (row_module_key(value), value)
                exact_candidate_images += 1
            return exact_image_cache[group_index]

        minimum_indices = minimum_packed_bucket(modular_images)
        minimum_exact = [
            (*exact_image(int(group_index)), int(group_index))
            for group_index in minimum_indices
        ]
        canonical_key, canonical_scaled_basis, transporter_index = min(
            minimum_exact, key=lambda row: (row[0], row[2])
        )

        possible_indices = np.flatnonzero(
            np.isin(modular_void, modular_target_void)
        )
        member_seed_indices = set()
        stabilizer_indices = []
        for group_index_value in possible_indices:
            group_index = int(group_index_value)
            modular_key = packed_key(modular_images[group_index])
            possible_exact_keys = modular_key_to_exact_keys.get(modular_key, set())
            image_key, _image_basis = exact_image(group_index)
            if image_key == seed["exact_key"]:
                stabilizer_indices.append(group_index)
            if image_key in possible_exact_keys:
                member_seed_indices.update(exact_key_to_seed_indices[image_key])
        assert seed_index in member_seed_indices
        assert stabilizer_indices
        assert len(elements) % len(stabilizer_indices) == 0
        assert member_seed_indices <= unassigned
        unassigned.difference_update(member_seed_indices)

        canonical_ambient = QQ(1) / 3 * canonical_scaled_basis * root_basis
        assert all(value.denominator() == 1 for value in canonical_ambient.list())
        canonical_ambient = matrix(ZZ, canonical_ambient)
        assert row_module_key(
            canonical_ambient
        ) == row_module_key(seed["auxiliary_basis"] * ambient_action(transporter_index))
        raw_orbits.append(
            {
                "canonical_key": canonical_key,
                "representative_basis": canonical_ambient,
                "seed_indices": sorted(member_seed_indices),
                "source_seed_index": seed_index,
                "source_stabilizer_indices": stabilizer_indices,
                "transporter_index": transporter_index,
                "orbit_size": len(elements) // len(stabilizer_indices),
            }
        )
        if len(raw_orbits) % 25 == 0:
            print(
                "12A2CANON_PROGRESS|orbits={}|seeds_remaining={}".format(
                    len(raw_orbits), len(unassigned)
                ),
                flush=True,
            )

    canonical_keys = [orbit["canonical_key"] for orbit in raw_orbits]
    assert len(set(canonical_keys)) == len(canonical_keys)
    raw_orbits.sort(key=lambda orbit: orbit["canonical_key"])

    genera_cache = {}
    form_cache = {}
    orbit_records = []
    for orbit_index, orbit in enumerate(raw_orbits):
        auxiliary_basis = orbit["representative_basis"]
        auxiliary = auxiliary_basis * ambient * auxiliary_basis.transpose()
        complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
        frame = complement_basis * ambient * complement_basis.transpose()
        assert frame.det() == auxiliary.det()
        roots = root_type(frame)
        origin_class_ids = sorted(
            {seeds[index]["source_class_id"] for index in orbit["seed_indices"]}
        )

        transporter = ambient_action(orbit["transporter_index"])
        inverse_transporter = transporter.inverse().change_ring(ZZ)
        stabilizer = []
        for source_group_index in orbit["source_stabilizer_indices"]:
            conjugate = (
                inverse_transporter
                * ambient_action(source_group_index)
                * transporter
            )
            assert row_module_key(auxiliary_basis * conjugate) == row_module_key(
                auxiliary_basis
            )
            conjugate_element = root_coordinate_permutation(root_basis, conjugate)
            conjugate_key = permutation_key(conjugate_element)
            group_index = element_key_to_index[conjugate_key]
            complement_action = induced_action(complement_basis, conjugate)
            assert complement_action * frame * complement_action.transpose() == frame
            component_permutation, diagram_variants = component_data(
                conjugate_element
            )
            stabilizer.append(
                {
                    "residual_group_index_zero_based": group_index,
                    "residual_conjugacy_class_id": class_by_element_key[
                        conjugate_key
                    ],
                    "component_permutation_zero_based": component_permutation,
                    "component_diagram_variants": diagram_variants,
                    "root_coordinate_permutation_one_based": list(conjugate_key),
                    "order": finite_order(complement_action),
                    "moved_dimension_mod_2": int(
                        matrix(
                            GF(2),
                            complement_action - identity_matrix(ZZ, 17),
                        ).rank()
                    ),
                    "fixed_dimension_over_Q": int(
                        17
                        - (
                            complement_action - identity_matrix(ZZ, 17)
                        ).rank()
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
                "orbit_id": f"12A2-CF-O{orbit_index + 1:04d}",
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

    surfaces = surface_rows("12A2", orbit_records)
    compatible_orbit_ids = {
        orbit_id
        for surface in surfaces
        for orbit_id in surface["residual_group_orbit_ids"]
    }
    orbit_mw = Counter(
        row["root_data"]["mw_rank_for_rho_19"] for row in orbit_records
    )
    frame_mw = Counter(
        frame_row["root_data"]["mw_rank_for_rho_19"]
        for surface in surfaces
        for frame_row in surface["frames"]
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
        "ambient_label": "12A2",
        "backend_id": "ROOTED-12A2",
        "residual_group": {
            "order": len(elements),
            "component_permutation_image_order": residual["residual_group"][
                "component_permutation_image_order"
            ],
            "central_diagram_kernel_order": residual["residual_group"][
                "central_diagram_kernel_order"
            ],
            "conjugacy_classes": class_metadata,
        },
        "canonicalization": {
            "mod_2_group_images_checked": modular_images_checked,
            "exact_candidate_images_checked": exact_candidate_images,
            "equality_certificate": (
                "F_2 RREF is only a rejection filter; every possible accepted-"
                "seed equality, stabilizer, and minimum modular-bucket image is "
                "verified by exact integral row-module HNF in scaled A2 root "
                "coordinates"
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
        backend["accounting"]["k3_compatible_residual_group_embedding_orbits"],
        backend["accounting"]["surface_classes_after_T_NS_first_dedup"],
        backend["accounting"][
            "partner_auxiliary_isometry_classes_after_surface_dedup"
        ],
        backend["accounting"]["frame_isometry_classes_after_surface_dedup"],
    ) == (214, 210, 99, 108, 151)
    assert backend["accounting"]["residual_group_orbit_size_distribution"] == {
        "4752": 2,
        "9504": 4,
        "11880": 4,
        "15840": 7,
        "23760": 12,
        "31680": 15,
        "47520": 3,
        "63360": 61,
        "95040": 106,
    }
    assert backend["accounting"]["post_dedup_frame_mw_rank_distribution"] == {
        "12": 13,
        "13": 101,
        "14": 25,
        "15": 2,
        "17": 10,
    }
    assert backend["canonicalization"]["mod_2_group_images_checked"] == (
        214 * 190080
    )
    assert backend["canonicalization"]["exact_candidate_images_checked"] == 1914
    return {
        "schema": "elkies-k3.12a2-fixed-coordinate-shells.v1",
        "status": "PASS_EXACT_DECLARED_12A2_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        "proof_scope": {
            "proved": (
                "complete order-190080 residual lift group; every 7-of-r "
                "coordinate summand in all nine eligible fixed-class LLL bases; "
                "exact determinant, length, MW12--17, mod-two, residual quotient, "
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
        "elkies-k3/scripts/canonicalize_12a2_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("12A2 canonical shell artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "12A2CANON|seeds={}|orbits={}|surfaces={}|partners={}|frames={}|status=PASS_EXACT".format(
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
