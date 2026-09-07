#!/usr/bin/env sage-python
"""Enumerate an order-four symmetry-first shell inside N(4A5+D4).

Recover the four A5 and one D4 root components intrinsically and enumerate
every chamber-preserving component/diagram automorphism that preserves the
Niemeier overlattice.  The resulting exact residual group has order 48 and
component-permutation image S4.  Select the first order-four lift of a double
A5-component swap with fixed rank ten, enumerate all 7-of-10 coordinate
summands of its pinned LLL fixed-lattice basis, and impose determinant,
discriminant-length, MW12--17, and nontrivial mod-2 complement-action gates.

Survivors are canonicalized under the full residual group before the exact
ternary discriminant-form gate and (T,NS)-first deduplication.  The result is
exact in this coordinate language and residual quotient, not a complete
fixed-lattice, Weyl-orbit, or determinant-band census.
"""

from __future__ import annotations

import argparse
import itertools
import json
import runpy
from collections import Counter
from pathlib import Path

from sage.all import GF, Graph, QQ, ZZ, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-4a5-d4-order4-fixed-high-mw-seed-v1.json"
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
primitive_closure_index = COMMON["primitive_closure_index"]
discriminant_invariants = COMMON["discriminant_invariants"]
induced_action = COMMON["induced_action"]
signed_roots = COMMON["signed_roots"]
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
    raise AssertionError("unexpected residual action order")


def cycle_type(permutation):
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        current = min(unseen)
        length = 0
        while current in unseen:
            unseen.remove(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths))


def component_data(gram):
    roots = signed_roots(gram)
    unseen = set(range(len(roots)))
    components = []
    while unseen:
        component = {min(unseen)}
        pending = list(component)
        unseen.difference_update(component)
        while pending:
            current = pending.pop()
            adjacent = {
                index
                for index in unseen
                if roots[current] * gram * roots[index] != 0
            }
            component.update(adjacent)
            unseen.difference_update(adjacent)
            pending.extend(sorted(adjacent))
        component_roots = [roots[index] for index in sorted(component)]
        components.append(
            {
                "invariant": (
                    int(matrix(ZZ, component_roots).rank()),
                    len(component_roots),
                ),
                "roots": component_roots,
            }
        )
    assert [row["invariant"] for row in components] == [
        (4, 24),
        (5, 30),
        (5, 30),
        (5, 30),
        (5, 30),
    ]
    for component in components:
        roots = component["roots"]
        for trial in range(1, 100):
            chamber = vector(
                ZZ,
                [
                    (index + 1) ** 2 + trial * (index + 1) + trial**2
                    for index in range(24)
                ],
            )
            values = [root * gram * chamber for root in roots]
            if all(value != 0 for value in values):
                break
        positive = [root for root, value in zip(roots, values) if value > 0]
        positive_set = {tuple(map(int, root)) for root in positive}
        basis = matrix(
            ZZ,
            [
                root
                for root in positive
                if not any(
                    tuple(map(int, root - other)) in positive_set
                    for other in positive
                )
            ],
        )
        rank = component["invariant"][0]
        cartan = basis * gram * basis.transpose()
        assert basis.nrows() == basis.rank() == rank
        assert all(cartan[index, index] == 2 for index in range(rank))
        assert sum(
            cartan[row, column] == -1
            for row in range(rank)
            for column in range(row)
        ) == rank - 1
        component["simple_basis"] = basis
    return components


def dynkin_graph(basis, gram):
    cartan = basis * gram * basis.transpose()
    graph = Graph()
    graph.add_vertices(range(basis.nrows()))
    graph.add_edges(
        (row, column)
        for row in range(basis.nrows())
        for column in range(row)
        if cartan[row, column] == -1
    )
    return graph


def graph_variants(left, right, gram):
    left_graph = dynkin_graph(left, gram)
    right_graph = dynkin_graph(right, gram)
    isomorphic, initial = left_graph.is_isomorphic(
        right_graph, certificate=True
    )
    assert isomorphic
    return [
        [automorphism(initial[index]) for index in range(left.nrows())]
        for automorphism in right_graph.automorphism_group()
    ]


def residual_group(gram):
    components = component_data(gram)
    bases = [row["simple_basis"] for row in components]
    invariants = [row["invariant"] for row in components]
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    root_basis_inverse = root_basis.inverse()
    groups = {}
    for index, invariant in enumerate(invariants):
        groups.setdefault(invariant, []).append(index)
    variants = {
        (source, target): graph_variants(
            bases[source], bases[target], gram
        )
        for source in range(5)
        for target in range(5)
        if invariants[source] == invariants[target]
    }

    elements = []
    seen = set()
    type_permutations = itertools.product(
        *[list(itertools.permutations(indices)) for indices in groups.values()]
    )
    for choices in type_permutations:
        permutation = list(range(5))
        for indices, images in zip(groups.values(), choices):
            for source, target in zip(indices, images):
                permutation[source] = target
        permutation = tuple(permutation)
        for selections in itertools.product(
            *[
                variants[source, permutation[source]]
                for source in range(5)
            ]
        ):
            target_basis = matrix(
                QQ,
                [
                    bases[permutation[source]][selections[source][row]]
                    for source in range(5)
                    for row in range(bases[source].nrows())
                ],
            )
            action = root_basis_inverse * target_basis
            if not all(entry.denominator() == 1 for entry in action.list()):
                continue
            action = matrix(ZZ, action)
            if abs(action.det()) != 1:
                continue
            assert action * gram * action.transpose() == gram
            key = matrix_key(action)
            assert key not in seen
            seen.add(key)
            order = finite_order(action)
            fixed_rank = int(24 - (action - identity_matrix(ZZ, 24)).rank())
            elements.append(
                {
                    "matrix": action,
                    "component_permutation": permutation,
                    "cycle_type": cycle_type(permutation),
                    "order": order,
                    "fixed_rank": fixed_rank,
                }
            )
    elements.sort(
        key=lambda row: (
            row["component_permutation"],
            row["order"],
            -row["fixed_rank"],
            matrix_key(row["matrix"]),
        )
    )
    assert len(elements) == 48
    assert len({row["component_permutation"] for row in elements}) == 24
    profile = Counter(
        (row["cycle_type"], row["order"], row["fixed_rank"])
        for row in elements
    )
    assert profile == Counter(
        {
            ((1, 1, 1, 1, 1), 1, 24): 1,
            ((1, 1, 1, 1, 1), 2, 16): 1,
            ((1, 1, 1, 2), 2, 16): 12,
            ((1, 2, 2), 4, 10): 6,
            ((1, 1, 3), 3, 12): 8,
            ((1, 1, 3), 6, 8): 8,
            ((1, 4), 8, 6): 12,
        }
    )

    # Exact closure proof by greedily enlarging a generating set inside the
    # enumerated matrix set.
    by_key = {matrix_key(row["matrix"]): row for row in elements}
    identity = identity_matrix(ZZ, 24)
    generated = {matrix_key(identity): identity}
    generators = []

    def close():
        pending = list(generated.values())
        while pending:
            current = pending.pop()
            for generator in generators:
                product = current * generator
                key = matrix_key(product)
                assert key in by_key
                if key not in generated:
                    generated[key] = product
                    pending.append(product)

    for row in elements:
        if matrix_key(row["matrix"]) in generated:
            continue
        generators.append(row["matrix"])
        close()
        if len(generated) == 48:
            break
    assert set(generated) == set(by_key)

    for index, row in enumerate(elements):
        row["section_index_zero_based"] = index
        row["class"] = "CYCLE_{}_ORDER_{}_FIXED_{}".format(
            "_".join(map(str, row["cycle_type"])),
            row["order"],
            row["fixed_rank"],
        )
    selected = next(
        row
        for row in elements
        if row["component_permutation"] == (0, 2, 1, 4, 3)
        and row["order"] == 4
        and row["fixed_rank"] == 10
    )
    return components, elements, selected, generators, profile


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


def build(catalog, determinant_bound, minimum_mw_rank):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row
        for row in catalog["rooted_niemeier_lattices"]
        if row["label"] == "4A5_D4"
    )
    ambient = matrix(ZZ, ambient_row["gram"])
    components, section, selected, generators, profile = residual_group(ambient)
    identity24 = identity_matrix(ZZ, 24)
    fixed = row_module_basis(
        (selected["matrix"] - identity24)
        .transpose()
        .right_kernel_matrix()
        .change_ring(ZZ)
    )
    assert fixed.nrows() == 10 and primitive_closure_index(fixed) == 1
    fixed_gram = fixed * ambient * fixed.transpose()
    assert fixed_gram.det() == 1024
    lll_change = fixed_gram.LLL_gram().transpose()
    assert abs(lll_change.det()) == 1
    reduced_basis = lll_change * fixed
    reduced_gram = reduced_basis * ambient * reduced_basis.transpose()

    counters = Counter()
    seed_records = []
    for combination in itertools.combinations(range(10), 7):
        counters["coordinate_subsets_tested"] += 1
        auxiliary_basis = matrix(
            ZZ, [reduced_basis.row(index) for index in combination]
        )
        assert primitive_closure_index(auxiliary_basis) == 1
        auxiliary = auxiliary_basis * ambient * auxiliary_basis.transpose()
        if auxiliary.det() > determinant_bound:
            counters["determinant_rejected"] += 1
            continue
        invariants = discriminant_invariants(auxiliary)
        if len(invariants) > 3:
            counters["discriminant_length_rejected"] += 1
            continue
        complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
        frame = complement_basis * ambient * complement_basis.transpose()
        assert frame.det() == auxiliary.det()
        roots = root_type(frame)
        counters[f"mw_rank_{roots['mw_rank_for_rho_19']}"] += 1
        if roots["mw_rank_for_rho_19"] < minimum_mw_rank:
            counters["mw_rank_below_factory_floor_rejected"] += 1
            continue
        action = induced_action(complement_basis, selected["matrix"])
        moved_dimension = int(
            matrix(GF(2), action - identity_matrix(ZZ, 17)).rank()
        )
        if moved_dimension == 0:
            counters["mod2_trivial_rejected"] += 1
            continue
        seed_records.append(
            {
                "coordinate_subset_zero_based": list(combination),
                "auxiliary_basis": auxiliary_basis,
                "root_data": roots,
                "selected_action_moved_dimension_mod_2": moved_dimension,
            }
        )

    orbit_by_key = {}
    for seed_index, seed in enumerate(seed_records):
        images = {}
        for group_row in section:
            image = row_module_basis(seed["auxiliary_basis"] * group_row["matrix"])
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
        roots = root_type(frame)
        stabilizer = []
        for group_row in section:
            if row_module_key(auxiliary_basis * group_row["matrix"]) != canonical_key:
                continue
            action = induced_action(complement_basis, group_row["matrix"])
            assert action * frame * action.transpose() == frame
            stabilizer.append(
                {
                    "section_index_zero_based": group_row[
                        "section_index_zero_based"
                    ],
                    "class": group_row["class"],
                    "component_permutation_zero_based": list(
                        group_row["component_permutation"]
                    ),
                    "order": finite_order(action),
                    "moved_dimension_mod_2": int(
                        matrix(
                            GF(2), action - identity_matrix(ZZ, 17)
                        ).rank()
                    ),
                    "fixed_dimension_over_Q": int(
                        17 - (action - identity_matrix(ZZ, 17)).rank()
                    ),
                }
            )
        assert any(
            row["class"] == selected["class"]
            and row["moved_dimension_mod_2"] > 0
            for row in stabilizer
        )
        gate = ternary_gate(frame, genera_cache, form_cache)
        orbit_records.append(
            {
                "orbit_id": f"4A5D4-O4-O{orbit_index + 1:04d}",
                "residual_group_orbit_size": len(orbit["image_bases"]),
                "seed_indices_zero_based": orbit["seed_indices"],
                "seed_coordinate_subsets_zero_based": [
                    seed_records[index]["coordinate_subset_zero_based"]
                    for index in orbit["seed_indices"]
                ],
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
                    "surface_id": f"K3-4A5D4O4-{compact_digest(key)}",
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

    surface_rows = []
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
        for index, frame_class in enumerate(surface["frame_classes"], start=1):
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
                        for partner in frame_class["partner_class_objects"]
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
        surface_rows.append(
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
    surface_rows.sort(key=lambda row: (row["determinant"], row["surface_id"]))

    all_mw = Counter(
        {
            int(key.removeprefix("mw_rank_")): value
            for key, value in counters.items()
            if key.startswith("mw_rank_")
            and key.removeprefix("mw_rank_").isdigit()
        }
    )
    seed_mw = Counter(
        seed["root_data"]["mw_rank_for_rho_19"] for seed in seed_records
    )
    orbit_mw = Counter(
        orbit["root_data"]["mw_rank_for_rho_19"] for orbit in orbit_records
    )
    frame_mw = Counter(
        frame["root_data"]["mw_rank_for_rho_19"]
        for surface in surface_rows
        for frame in surface["frames"]
    )
    moved = Counter(
        max(
            action["moved_dimension_mod_2"]
            for action in orbit["literal_residual_stabilizer"]
            if action["class"] == selected["class"]
        )
        for orbit in orbit_records
    )
    orbit_sizes = Counter(
        orbit["residual_group_orbit_size"] for orbit in orbit_records
    )

    assert counters["coordinate_subsets_tested"] == 120
    assert counters["determinant_rejected"] == 0
    assert counters["discriminant_length_rejected"] == 2
    assert all_mw == Counter({5: 18, 7: 2, 9: 59, 13: 36, 17: 3})
    assert seed_mw == Counter({13: 36, 17: 3})
    assert counters["mod2_trivial_rejected"] == 0
    assert len(orbit_records) == 39
    assert orbit_sizes == Counter({6: 20, 12: 19})
    assert orbit_mw == Counter({13: 36, 17: 3})
    assert moved == Counter({11: 39})
    assert len(form_cache) == 9
    assert all(
        gate["matching_even_ternary_genera"] == 1
        for gate in form_cache.values()
    )
    assert len(surface_rows) == 9
    assert sum(len(row["partner_auxiliaries"]) for row in surface_rows) == 9
    assert sum(len(row["frames"]) for row in surface_rows) == 9
    assert frame_mw == Counter({13: 7, 17: 2})

    return {
        "schema": "elkies-k3.4a5-d4-order4-fixed-high-mw-seed.v1",
        "status": "PASS_EXACT_DECLARED_4A5_D4_ORDER4_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST",
        "proof_scope": {
            "proved": (
                "All 48 chamber-preserving component/diagram lifts are enumerated "
                "and form an exact residual group with S4 component image. All "
                "7-of-10 coordinate summands of the selected order-four fixed "
                "lattice basis are tested. Retained primitive auxiliaries pass the "
                "determinant, length, MW, nontrivial mod-2, full residual quotient, "
                "ternary-genus, and (T,NS)-first deduplication gates."
            ),
            "not_proved": (
                "This coordinate shell is not all primitive rank-seven sublattices "
                "of Fix(g) or N(4A5+D4). The chamber residual quotient is not the "
                "full Weyl quotient, and ternary genera are not class enumerations."
            ),
        },
        "parameters": {
            "ambient": "4A5_D4",
            "selected_component_permutation_zero_based": list(
                selected["component_permutation"]
            ),
            "selected_class": selected["class"],
            "determinant_bound": determinant_bound,
            "discriminant_length_bound": 3,
            "minimum_mw_rank": minimum_mw_rank,
            "maximum_mw_rank": 17,
            "seed_language": "7-of-10 coordinate direct summands of pinned LLL fixed-lattice basis",
        },
        "residual_group": {
            "order": len(section),
            "component_permutation_image_order": len(
                {row["component_permutation"] for row in section}
            ),
            "generator_count": len(generators),
            "class_distribution": {
                "{}|order={}|fixed_rank={}".format(
                    ",".join(map(str, cycles)), order, fixed_rank
                ): count
                for (cycles, order, fixed_rank), count in sorted(profile.items())
            },
            "simple_root_bases_in_ambient": [
                rows(component["simple_basis"]) for component in components
            ],
            "elements": [
                {
                    "section_index_zero_based": row["section_index_zero_based"],
                    "class": row["class"],
                    "order": row["order"],
                    "fixed_rank": row["fixed_rank"],
                    "component_permutation_zero_based": list(
                        row["component_permutation"]
                    ),
                    "ambient_matrix": rows(row["matrix"]),
                }
                for row in section
            ],
        },
        "fixed_lattice": {
            "rank": 10,
            "determinant": int(fixed_gram.det()),
            "basis_in_ambient": rows(fixed),
            "gram": rows(fixed_gram),
            "lll_change": rows(lll_change),
            "lll_basis_in_ambient": rows(reduced_basis),
            "lll_gram": rows(reduced_gram),
            "primitive_in_ambient": True,
        },
        "accounting": {
            "coordinate_subsets_tested": counters["coordinate_subsets_tested"],
            "determinant_rejected": counters["determinant_rejected"],
            "discriminant_length_rejected": counters[
                "discriminant_length_rejected"
            ],
            "all_length_admissible_mw_rank_distribution": {
                str(key): value for key, value in sorted(all_mw.items())
            },
            "mw_rank_below_factory_floor_rejected": counters[
                "mw_rank_below_factory_floor_rejected"
            ],
            "mod2_trivial_rejected": counters["mod2_trivial_rejected"],
            "high_mw_mod2_accepted_seeds": len(seed_records),
            "accepted_seed_mw_rank_distribution": {
                str(key): value for key, value in sorted(seed_mw.items())
            },
            "residual_group_embedding_orbits": len(orbit_records),
            "residual_group_orbit_size_distribution": {
                str(key): value for key, value in sorted(orbit_sizes.items())
            },
            "orbit_mw_rank_distribution": {
                str(key): value for key, value in sorted(orbit_mw.items())
            },
            "post_dedup_frame_mw_rank_distribution": {
                str(key): value for key, value in sorted(frame_mw.items())
            },
            "selected_class_moved_dimension_mod_2_distribution": {
                str(key): value for key, value in sorted(moved.items())
            },
            "distinct_frame_discriminant_forms": len(form_cache),
            "forms_with_matching_ternary_genus": sum(
                gate["matching_even_ternary_genera"] > 0
                for gate in form_cache.values()
            ),
            "surface_classes_after_T_NS_first_dedup": len(surface_rows),
            "partner_auxiliary_isometry_classes_after_surface_dedup": sum(
                len(surface["partner_auxiliaries"]) for surface in surface_rows
            ),
            "frame_isometry_classes_after_surface_dedup": sum(
                len(surface["frames"]) for surface in surface_rows
            ),
        },
        "embedding_orbits": orbit_records,
        "surfaces_T_NS_first": surface_rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--determinant-bound", type=int, default=5000)
    parser.add_argument("--minimum-mw-rank", type=int, default=12)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    assert arguments.determinant_bound == 5000
    assert arguments.minimum_mw_rank == 12
    result = build(
        json.loads(arguments.catalog.read_text()),
        arguments.determinant_bound,
        arguments.minimum_mw_rank,
    )
    result["inputs"] = {
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(
            arguments.catalog
        ),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
    }
    result["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/enumerate_4a5_d4_order4_fixed_high_mw_seed.sage"
    )
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("4A5+D4 order-four fixed seed artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    accounting = result["accounting"]
    print(
        "4A5D4O4|seeds={}|orbits={}|surfaces={}|partners={}|frames={}|status=PASS_EXACT".format(
            accounting["high_mw_mod2_accepted_seeds"],
            accounting["residual_group_embedding_orbits"],
            accounting["surface_classes_after_T_NS_first_dedup"],
            accounting[
                "partner_auxiliary_isometry_classes_after_surface_dedup"
            ],
            accounting["frame_isometry_classes_after_surface_dedup"],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
