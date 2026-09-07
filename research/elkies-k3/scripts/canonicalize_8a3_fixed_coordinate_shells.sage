#!/usr/bin/env sage-python
"""Canonicalize the exact N(8A3) fixed-coordinate shell scan.

Quotient every accepted coordinate seed by the complete order-2688 residual
group.  A batched row-reduced form over F_251 is used only as a rejection
filter and to select a full-orbit canonical bucket.  Every possible equality,
stabilizer element, and canonical representative in that bucket is then
verified over Z by Hermite-normalized primitive row modules.

After the residual quotient, compute literal stabilizer actions on the rank-17
complements, apply the exact ternary discriminant-form gate, and deduplicate
first by (T,NS), then by auxiliary and frame isometry.
"""

from __future__ import annotations

import argparse
import json
import runpy
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sage.all import GF, ZZ, identity_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
RESIDUAL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-8a3-glue-code-residual-group-v1.json"
)
PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-8a3-fixed-coordinate-shell-probe-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-8a3-fixed-coordinate-shells-v1.json"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
CANONICAL_COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "canonicalize_4a6_4e6_fixed_coordinate_shells.sage"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_fixed_seed_common")
CANONICAL_COMMON = runpy.run_path(
    str(CANONICAL_COMMON_SOURCE), run_name="_fixed_shell_canonical_common"
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
finite_order = CANONICAL_COMMON["finite_order"]
rootless_intrinsics = CANONICAL_COMMON["rootless_intrinsics"]
surface_rows = CANONICAL_COMMON["surface_rows"]

MODULAR_PRIME = 251


def batched_rref_mod_prime(images, prime=MODULAR_PRIME):
    """Return exact batched RREFs for full-row-rank integer matrices mod p."""
    values = np.remainder(images, prime).astype(np.int64, copy=True)
    count, rank, columns = values.shape
    inverses = np.zeros(prime, dtype=np.int64)
    for value in range(1, prime):
        inverses[value] = pow(value, prime - 2, prime)
    for pivot_row in range(rank):
        waiting = np.ones(count, dtype=bool)
        for column in range(columns):
            indices = np.flatnonzero(waiting)
            if not len(indices):
                break
            block = values[indices, pivot_row:, column]
            has_pivot = np.any(block != 0, axis=1)
            selected = indices[has_pivot]
            if not len(selected):
                continue
            source_rows = (
                np.argmax(block[has_pivot] != 0, axis=1) + pivot_row
            )
            temporary = values[selected, pivot_row, :].copy()
            values[selected, pivot_row, :] = values[
                selected, source_rows, :
            ]
            values[selected, source_rows, :] = temporary
            factors = inverses[values[selected, pivot_row, column]]
            values[selected, pivot_row, :] = (
                values[selected, pivot_row, :] * factors[:, None]
            ) % prime
            for other_row in range(rank):
                if other_row == pivot_row:
                    continue
                factors = values[selected, other_row, column].copy()
                values[selected, other_row, :] = (
                    values[selected, other_row, :]
                    - factors[:, None] * values[selected, pivot_row, :]
                ) % prime
            waiting[selected] = False
        assert not np.any(waiting), "primitive row module lost rank modulo p"
    return values


def modular_keys(rrefs):
    return [row.astype(np.uint8).tobytes() for row in rrefs]


def validate_modular_engine(seed_basis, integer_actions, modular_actions):
    images = np.einsum(
        "ij,njk->nik",
        np.asarray(seed_basis, dtype=np.int64),
        integer_actions,
        optimize=True,
    )
    batched = batched_rref_mod_prime(images)
    field = GF(MODULAR_PRIME)
    seed = matrix(field, seed_basis)
    for index, action in enumerate(modular_actions):
        reference = np.asarray(
            (seed * action).echelon_form(), dtype=np.int64
        )
        assert np.array_equal(batched[index], reference)


def canonicalize(catalog, residual, probe):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert residual["schema"] == "elkies-k3.8a3-glue-code-residual-group.v1"
    assert residual["status"] == "PASS_EXACT_8A3_GLUE_CODE_AND_RESIDUAL_GROUP"
    assert probe["schema"] == "elkies-k3.8a3-fixed-coordinate-shell-probe.v1"
    assert probe["status"] == (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_8A3_COORDINATE_SHELL_SCAN"
    )
    assert probe["parameters"]["determinant_bound"] == 5000
    assert probe["parameters"]["minimum_mw_rank"] == 12

    ambient = matrix(
        ZZ,
        next(
            row["gram"]
            for row in catalog["rooted_niemeier_lattices"]
            if row["label"] == "8A3"
        ),
    )
    element_rows = residual["residual_group"]["elements"]
    assert residual["residual_group"]["order"] == len(element_rows) == 2688
    actions = [matrix(ZZ, row["matrix"]) for row in element_rows]
    for action in actions:
        assert action * ambient * action.transpose() == ambient
    integer_actions = np.asarray(
        [np.asarray(action, dtype=np.int64) for action in actions],
        dtype=np.int64,
    )
    modular_actions = [matrix(GF(MODULAR_PRIME), action) for action in actions]

    class_by_signed_key = {}
    class_metadata = {}
    for class_row in residual["residual_group"]["conjugacy_classes"]:
        class_id = class_row["class_id"]
        class_metadata[class_id] = {
            "action_order": class_row["action_order"],
            "fixed_rank": class_row["fixed_rank"],
            "fixed_determinant": class_row["fixed_determinant"],
            "class_size": class_row["class_size"],
        }
        for permutation, signs in class_row["element_signed_keys"]:
            key = (tuple(permutation), tuple(signs))
            assert key not in class_by_signed_key
            class_by_signed_key[key] = class_id
    assert len(class_by_signed_key) == len(actions)

    group = []
    group_index_by_matrix_key = {}
    for index, (element, action) in enumerate(zip(element_rows, actions)):
        signed_key = (
            tuple(element["component_permutation_zero_based"]),
            tuple(element["component_diagram_signs"]),
        )
        group.append(
            {
                "index_zero_based": index,
                "matrix": action,
                "class_id": class_by_signed_key[signed_key],
                "component_permutation_zero_based": element[
                    "component_permutation_zero_based"
                ],
                "component_diagram_signs": element["component_diagram_signs"],
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
    assert len(seeds) == 1166

    exact_key_to_seed_indices = defaultdict(list)
    modular_key_to_exact_keys = defaultdict(set)
    for index, seed in enumerate(seeds):
        exact_key_to_seed_indices[seed["exact_key"]].append(index)
        reduced = batched_rref_mod_prime(
            np.asarray(seed["auxiliary_basis"], dtype=np.int64)[None, :, :]
        )
        seed["modular_key"] = modular_keys(reduced)[0]
        modular_key_to_exact_keys[seed["modular_key"]].add(seed["exact_key"])

    validate_modular_engine(
        seeds[0]["auxiliary_basis"], integer_actions, modular_actions
    )

    unassigned = set(range(len(seeds)))
    raw_orbits = []
    modular_candidate_images = 0
    exact_candidate_images = 0
    while unassigned:
        seed_index = min(unassigned)
        seed = seeds[seed_index]
        basis = seed["auxiliary_basis"]
        image_arrays = np.einsum(
            "ij,njk->nik",
            np.asarray(basis, dtype=np.int64),
            integer_actions,
            optimize=True,
        )
        reduced = batched_rref_mod_prime(image_arrays)
        image_modular_keys = modular_keys(reduced)

        minimum_modular_key = min(image_modular_keys)
        minimum_indices = [
            index
            for index, key in enumerate(image_modular_keys)
            if key == minimum_modular_key
        ]
        minimum_exact = []
        exact_image_key_cache = {}
        for group_index in minimum_indices:
            image_basis = row_module_basis(basis * actions[group_index])
            key = row_module_key(image_basis)
            exact_image_key_cache[group_index] = key
            minimum_exact.append((key, group_index, image_basis))
        canonical_key, transporter_index, canonical_basis = min(
            minimum_exact, key=lambda row: (row[0], row[1])
        )

        member_seed_indices = set()
        stabilizer_indices = []
        for group_index, modular_key in enumerate(image_modular_keys):
            possible_exact_keys = modular_key_to_exact_keys.get(modular_key)
            if not possible_exact_keys and modular_key != seed["modular_key"]:
                continue
            modular_candidate_images += 1
            if group_index not in exact_image_key_cache:
                image_basis = row_module_basis(basis * actions[group_index])
                exact_image_key_cache[group_index] = row_module_key(image_basis)
            exact_candidate_images += 1
            image_key = exact_image_key_cache[group_index]
            if image_key == seed["exact_key"]:
                stabilizer_indices.append(group_index)
            if image_key in possible_exact_keys:
                member_seed_indices.update(exact_key_to_seed_indices[image_key])

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
                "source_seed_index": seed_index,
                "source_stabilizer_indices": stabilizer_indices,
                "transporter_index": transporter_index,
                "orbit_size": len(actions) // len(stabilizer_indices),
            }
        )
        if len(raw_orbits) % 100 == 0:
            print(
                "8A3CANON_PROGRESS|orbits={}|seeds_remaining={}".format(
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
            assert row_module_key(
                auxiliary_basis * group_row["matrix"]
            ) == canonical_key
            action = induced_action(complement_basis, group_row["matrix"])
            assert action * frame * action.transpose() == frame
            stabilizer.append(
                {
                    "residual_group_index_zero_based": group_index,
                    "residual_conjugacy_class_id": group_row["class_id"],
                    "component_permutation_zero_based": group_row[
                        "component_permutation_zero_based"
                    ],
                    "component_diagram_signs": group_row[
                        "component_diagram_signs"
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
                "orbit_id": f"8A3-CF-O{orbit_index + 1:04d}",
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

    surfaces = surface_rows("8A3", orbit_records)
    k3_compatible_orbit_ids = {
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
        "ambient_label": "8A3",
        "backend_id": "ROOTED-8A3",
        "residual_group": {
            "order": len(group),
            "component_permutation_image_order": residual["residual_group"][
                "component_permutation_image_order"
            ],
            "conjugacy_classes": class_metadata,
        },
        "canonicalization": {
            "modular_rejection_prime": MODULAR_PRIME,
            "modular_engine_validated_against_sage_images": len(group),
            "modular_candidate_images_exactly_checked": modular_candidate_images,
            "exact_candidate_images_checked": exact_candidate_images,
            "equality_certificate": (
                "primitive rank-seven row modules have full rank modulo 251; "
                "unequal modular RREFs cannot be equal integral row modules, "
                "and every equal modular candidate is checked by integral HNF"
            ),
        },
        "source_probe_accounting": probe["accounting"],
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
    ) == (1162, 1160, 435, 523, 574)
    assert backend["accounting"]["post_dedup_frame_mw_rank_distribution"] == {
        "12": 257,
        "13": 245,
        "14": 54,
        "15": 6,
        "17": 12,
    }
    return {
        "schema": "elkies-k3.8a3-fixed-coordinate-shells.v1",
        "status": "PASS_EXACT_DECLARED_8A3_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        "proof_scope": {
            "proved": (
                "complete order-2688 residual lift group; every 7-of-r "
                "coordinate summand in every eligible fixed-class LLL basis; "
                "exact determinant, length, MW12--17, mod-two, full residual "
                "quotient, ternary, and T/NS-first gates"
            ),
            "not_proved": (
                "all primitive rank-seven fixed-lattice or ambient "
                "sublattices, full Weyl embedding orbits, ternary class "
                "enumeration beyond genera, or determinant-band completeness"
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
        str(CANONICAL_COMMON_SOURCE.resolve().relative_to(ROOT)): digest(
            CANONICAL_COMMON_SOURCE
        ),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/canonicalize_8a3_fixed_coordinate_shells.sage"
    )
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("8A3 canonical shell artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "8A3CANON|seeds={}|orbits={}|surfaces={}|partners={}|frames={}|status=PASS_EXACT".format(
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
