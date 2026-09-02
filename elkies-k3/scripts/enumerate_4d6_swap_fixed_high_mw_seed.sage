#!/usr/bin/env sage-python
"""Enumerate a symmetry-first high-MW seed shell inside N(4D6).

Recover the four D6 root components intrinsically, lift their component
permutations through the Niemeier glue, and certify an exact S4 section.  The
selected literal component transposition has a primitive rank-sixteen fixed
lattice.  Enumerate all 7-of-16 coordinate direct summands of a pinned LLL
basis of that fixed lattice, then impose determinant, discriminant-length,
MW12--17, and nontrivial mod-2 complement-action gates.

Survivors are closed and canonicalized under the exact S4 section before the
ternary discriminant-form gate and required (T,NS)-first deduplication.  This
is a deterministic bounded discovery shell, not a complete enumeration of
rank-seven sublattices of N(4D6) or of the transposition fixed lattice.

status: EXACT_DECLARED_4D6_SWAP_FIXED_COORDINATE_SEED_SHELL
"""

from __future__ import annotations

import argparse
import hashlib
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
    / "artifacts/generated-results/elkies-k3-4d6-swap-fixed-high-mw-seed-v1.json"
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
multiplicative_order = COMMON["multiplicative_order"]
signed_roots = COMMON["signed_roots"]
root_type = COMMON["root_type"]
ternary_gate = COMMON["ternary_gate"]
section_is_group = COMMON["section_is_group"]
find_or_add_isometry_class = COMMON["find_or_add_isometry_class"]


def component_roots(gram):
    roots = signed_roots(gram)
    unseen = set(range(len(roots)))
    result = []
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
        vectors = [roots[index] for index in sorted(component)]
        assert matrix(ZZ, vectors).rank() == 6
        assert len(vectors) == 60
        result.append(vectors)
    assert len(result) == 4
    return result


def simple_roots(gram, roots):
    for trial in range(1, 100):
        chamber = vector(
            ZZ,
            [
                (index + 1) ** 2 + trial * (index + 1) + trial**2
                for index in range(gram.nrows())
            ],
        )
        values = [root * gram * chamber for root in roots]
        if all(value != 0 for value in values):
            break
    else:
        raise AssertionError("failed to choose a regular chamber")
    positive = [root for root, value in zip(roots, values) if value > 0]
    positive_set = {tuple(map(int, root)) for root in positive}
    result = matrix(
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
    cartan = result * gram * result.transpose()
    assert result.nrows() == result.rank() == 6
    assert all(cartan[index, index] == 2 for index in range(6))
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(6)
        for column in range(6)
        if row != column
    )
    assert sum(cartan[row, column] == -1 for row in range(6) for column in range(row)) == 5
    return result


def dynkin_graph(cartan):
    graph = Graph()
    graph.add_vertices(range(cartan.nrows()))
    graph.add_edges(
        (row, column)
        for row in range(cartan.nrows())
        for column in range(row)
        if cartan[row, column] == -1
    )
    return graph


def graph_isometry_variants(left_basis, right_basis, gram):
    left = dynkin_graph(left_basis * gram * left_basis.transpose())
    right = dynkin_graph(right_basis * gram * right_basis.transpose())
    isomorphic, initial = left.is_isomorphic(right, certificate=True)
    assert isomorphic
    return [
        [automorphism(initial[index]) for index in range(left_basis.nrows())]
        for automorphism in right.automorphism_group()
    ]


def first_lifted_component_permutation(gram, bases, root_basis, permutation):
    variants = [
        graph_isometry_variants(bases[source], bases[target], gram)
        for source, target in enumerate(permutation)
    ]
    for selections in itertools.product(*variants):
        target_basis = matrix(
            QQ,
            [
                bases[permutation[source]][selections[source][row]]
                for source in range(4)
                for row in range(6)
            ],
        )
        action = root_basis.inverse() * target_basis
        if not all(entry.denominator() == 1 for entry in action.list()):
            continue
        action = matrix(ZZ, action)
        if abs(action.det()) != 1:
            continue
        assert action * gram * action.transpose() == gram
        return action
    return None


def compose_component_permutations(left, right):
    # Row vectors act first by left, then by right.
    return tuple(right[left[index]] for index in range(len(left)))


def component_cycle_type(permutation):
    unseen = set(range(len(permutation)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths))


def section_class(permutation):
    cycle_type = component_cycle_type(permutation)
    return {
        (1, 1, 1, 1): "1A",
        (1, 1, 2): "2A_COMPONENT_SWAP",
        (2, 2): "2B_DOUBLE_COMPONENT_SWAP",
        (1, 3): "3A_COMPONENT_CYCLE",
        (4,): "4A_COMPONENT_CYCLE",
    }[cycle_type]


def generated_section(generators):
    identity = identity_matrix(ZZ, 24)
    identity_permutation = tuple(range(4))
    by_matrix = {
        matrix_key(identity): {
            "matrix": identity,
            "component_permutation": identity_permutation,
        }
    }
    pending = [by_matrix[matrix_key(identity)]]
    while pending:
        current = pending.pop()
        for generator in generators:
            action = current["matrix"] * generator["matrix"]
            permutation = compose_component_permutations(
                current["component_permutation"],
                generator["component_permutation"],
            )
            key = matrix_key(action)
            if key in by_matrix:
                assert by_matrix[key]["component_permutation"] == permutation
                continue
            by_matrix[key] = {
                "matrix": action,
                "component_permutation": permutation,
            }
            pending.append(by_matrix[key])
            assert len(by_matrix) <= 24
    result = sorted(
        by_matrix.values(), key=lambda row: row["component_permutation"]
    )
    for index, row in enumerate(result):
        row["section_index_zero_based"] = index
        row["class"] = section_class(row["component_permutation"])
        row["order"] = multiplicative_order(row["matrix"])
    assert len(result) == 24
    assert {row["component_permutation"] for row in result} == set(
        itertools.permutations(range(4))
    )
    assert section_is_group(result)
    assert Counter(row["class"] for row in result) == Counter(
        {
            "1A": 1,
            "2A_COMPONENT_SWAP": 6,
            "2B_DOUBLE_COMPONENT_SWAP": 3,
            "3A_COMPONENT_CYCLE": 8,
            "4A_COMPONENT_CYCLE": 6,
        }
    )
    return result


def exact_s4_section(ambient):
    components = component_roots(ambient)
    bases = [simple_roots(ambient, roots) for roots in components]
    root_basis = matrix(ZZ, [root for basis in bases for root in basis.rows()])
    assert root_basis.nrows() == root_basis.ncols() == root_basis.rank() == 24
    lifts = {}
    for permutation in itertools.permutations(range(4)):
        lift = first_lifted_component_permutation(
            ambient, bases, root_basis, permutation
        )
        assert lift is not None
        lifts[permutation] = lift
    selected_permutation = (0, 1, 3, 2)
    triple_permutation = (1, 2, 0, 3)
    section = generated_section(
        [
            {
                "matrix": lifts[selected_permutation],
                "component_permutation": selected_permutation,
            },
            {
                "matrix": lifts[triple_permutation],
                "component_permutation": triple_permutation,
            },
        ]
    )
    selected = next(
        row
        for row in section
        if row["component_permutation"] == selected_permutation
    )
    assert selected["order"] == 2
    return bases, section, selected


def build(catalog, determinant_bound, minimum_mw_rank):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    ambient_row = next(
        row
        for row in catalog["rooted_niemeier_lattices"]
        if row["label"] == "4D6"
    )
    ambient = matrix(ZZ, ambient_row["gram"])
    identity24 = identity_matrix(ZZ, 24)
    simple_bases, section, selected = exact_s4_section(ambient)
    fixed = row_module_basis(
        (selected["matrix"] - identity24)
        .transpose()
        .right_kernel_matrix()
        .change_ring(ZZ)
    )
    assert fixed.nrows() == 16 and primitive_closure_index(fixed) == 1
    fixed_gram = fixed * ambient * fixed.transpose()
    assert fixed_gram.det() == 256
    lll_change = fixed_gram.LLL_gram().transpose()
    assert abs(lll_change.det()) == 1
    reduced_basis = lll_change * fixed
    reduced_gram = reduced_basis * ambient * reduced_basis.transpose()

    counters = Counter()
    seed_records = []
    for combination in itertools.combinations(range(16), 7):
        counters["coordinate_subsets_tested"] += 1
        auxiliary_basis = matrix(
            ZZ, [reduced_basis.row(index) for index in combination]
        )
        assert primitive_closure_index(auxiliary_basis) == 1
        auxiliary = auxiliary_basis * ambient * auxiliary_basis.transpose()
        determinant = int(auxiliary.det())
        if determinant > determinant_bound:
            counters["determinant_rejected"] += 1
            continue
        invariants = discriminant_invariants(auxiliary)
        if len(invariants) > 3:
            counters["discriminant_length_rejected"] += 1
            continue
        complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
        complement = complement_basis * ambient * complement_basis.transpose()
        assert complement.det() == auxiliary.det()
        roots = root_type(complement)
        counters[f"mw_rank_{roots['mw_rank_for_rho_19']}"] += 1
        if roots["mw_rank_for_rho_19"] < minimum_mw_rank:
            counters["mw_rank_below_factory_floor_rejected"] += 1
            continue
        selected_action = induced_action(complement_basis, selected["matrix"])
        moved_dimension = int(
            matrix(
                GF(2), selected_action - identity_matrix(ZZ, 17)
            ).rank()
        )
        if moved_dimension == 0:
            counters["mod2_trivial_rejected"] += 1
            continue
        counters["high_mw_mod2_accepted_seeds"] += 1
        seed_records.append(
            {
                "coordinate_subset_zero_based": list(combination),
                "auxiliary_basis": auxiliary_basis,
                "auxiliary_gram": auxiliary,
                "complement_basis": complement_basis,
                "frame_gram": complement,
                "discriminant_invariants": invariants,
                "root_data": roots,
                "selected_swap_moved_dimension_mod_2": moved_dimension,
            }
        )

    orbit_by_key = {}
    for seed_index, seed in enumerate(seed_records):
        images = {}
        for section_row in section:
            image = row_module_basis(
                seed["auxiliary_basis"] * section_row["matrix"]
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
        roots = root_type(frame)
        stabilizer = []
        for section_row in section:
            if (
                row_module_key(auxiliary_basis * section_row["matrix"])
                != canonical_key
            ):
                continue
            action = induced_action(complement_basis, section_row["matrix"])
            assert action * frame * action.transpose() == frame
            stabilizer.append(
                {
                    "section_index_zero_based": section_row[
                        "section_index_zero_based"
                    ],
                    "class": section_row["class"],
                    "component_permutation_zero_based": list(
                        section_row["component_permutation"]
                    ),
                    "order": multiplicative_order(action),
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
            row["class"] == "2A_COMPONENT_SWAP"
            and row["moved_dimension_mod_2"] > 0
            for row in stabilizer
        )
        gate = ternary_gate(frame, genera_cache, form_cache)
        orbit_records.append(
            {
                "orbit_id": f"4D6SWAP-O{orbit_index + 1:04d}",
                "section_orbit_size": len(orbit["image_bases"]),
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
                "literal_section_stabilizer": stabilizer,
                "k3_discriminant_gate": gate,
            }
        )

    surfaces_by_key = {}
    for orbit_index, orbit in enumerate(orbit_records):
        frame_key = orbit["k3_discriminant_gate"][
            "frame_discriminant_form_normal_key"
        ]
        ns_key = negative_form_key(frame_key)
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
                    "surface_id": f"K3-4D6SWAP-{compact_digest(key)}",
                    "surface_key": key,
                    "determinant": orbit["determinant"],
                    "orbit_indices": [],
                    "partner_classes": [],
                    "frame_classes": [],
                }
            surface = surfaces_by_key[encoded]
            surface["orbit_indices"].append(orbit_index)
            auxiliary_class = find_or_add_isometry_class(
                surface["partner_classes"],
                matrix(ZZ, orbit["auxiliary_gram"]),
                orbit_index,
            )
            frame_class = find_or_add_isometry_class(
                surface["frame_classes"],
                matrix(ZZ, orbit["frame_gram"]),
                orbit_index,
            )
            frame_class.setdefault("partner_class_objects", [])
            if auxiliary_class not in frame_class["partner_class_objects"]:
                frame_class["partner_class_objects"].append(auxiliary_class)

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
                    "section_orbit_ids": [
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
            frame_id = f"{surface['surface_id']}-F{index:03d}"
            frames.append(
                {
                    "frame_id": frame_id,
                    "gram": rows(frame_class["gram"]),
                    "gram_sha256": gram_digest(frame_class["gram"]),
                    "determinant": int(frame_class["gram"].det()),
                    "root_data": representative["root_data"],
                    "partner_ids": sorted(
                        partner_id_by_object[id(partner)]
                        for partner in frame_class["partner_class_objects"]
                    ),
                    "section_orbit_ids": [
                        orbit_records[value]["orbit_id"]
                        for value in frame_class["record_indices"]
                    ],
                    "representative_complement_basis_in_ambient": representative[
                        "complement_basis_in_ambient"
                    ],
                    "representative_literal_section_stabilizer": representative[
                        "literal_section_stabilizer"
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
                "section_orbit_ids": [
                    orbit_records[value]["orbit_id"]
                    for value in surface["orbit_indices"]
                ],
            }
        )
    surface_rows.sort(key=lambda row: (row["determinant"], row["surface_id"]))

    seed_mw_distribution = Counter(
        seed["root_data"]["mw_rank_for_rho_19"] for seed in seed_records
    )
    orbit_mw_distribution = Counter(
        orbit["root_data"]["mw_rank_for_rho_19"] for orbit in orbit_records
    )
    post_dedup_frame_mw_distribution = Counter(
        frame["root_data"]["mw_rank_for_rho_19"]
        for surface in surface_rows
        for frame in surface["frames"]
    )
    stabilizer_moved_distribution = Counter(
        max(
            action["moved_dimension_mod_2"]
            for action in orbit["literal_section_stabilizer"]
            if action["class"] == "2A_COMPONENT_SWAP"
        )
        for orbit in orbit_records
    )
    assert counters["coordinate_subsets_tested"] == 11440
    assert len(orbit_records) <= len(seed_records)
    assert all(
        orbit["section_orbit_size"] in (1, 2, 3, 4, 6, 8, 12, 24)
        for orbit in orbit_records
    )

    return {
        "schema": "elkies-k3.4d6-swap-fixed-high-mw-seed.v1",
        "status": "PASS_EXACT_DECLARED_4D6_SWAP_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST",
        "proof_scope": {
            "proved": (
                "The intrinsic D6 root components give an exact lifted S4 section. "
                "All 7-of-16 coordinate direct summands of the pinned LLL basis of "
                "the selected component-swap fixed lattice are tested. Every retained "
                "K is primitive, has determinant at most 5000 and discriminant length "
                "at most three, and has an exactly computed MW12--17 complement. The "
                "literal component-swap action is nontrivial modulo two. Survivors are "
                "canonicalized under the S4 section and deduplicated first by (T,NS)."
            ),
            "not_proved": (
                "This coordinate language does not exhaust primitive rank-seven "
                "sublattices of the fixed lattice or N(4D6). The S4 section is not the "
                "full Weyl group. Ternary genus representatives prove exact K3 "
                "realizability examples but do not enumerate every class in a genus."
            ),
        },
        "parameters": {
            "ambient": "4D6",
            "selected_component_permutation_zero_based": list(
                selected["component_permutation"]
            ),
            "selected_class": selected["class"],
            "determinant_bound": determinant_bound,
            "discriminant_length_bound": 3,
            "minimum_mw_rank": minimum_mw_rank,
            "maximum_mw_rank": 17,
            "seed_language": "7-of-16 coordinate direct summands of pinned LLL fixed-lattice basis",
        },
        "component_permutation_section": {
            "name": "exact lifted S4 section in Aut(N(4D6))",
            "order": len(section),
            "class_distribution": dict(
                sorted(Counter(row["class"] for row in section).items())
            ),
            "simple_root_bases_in_ambient": [rows(basis) for basis in simple_bases],
            "elements": [
                {
                    "section_index_zero_based": row["section_index_zero_based"],
                    "class": row["class"],
                    "order": row["order"],
                    "component_permutation_zero_based": list(
                        row["component_permutation"]
                    ),
                    "ambient_matrix": rows(row["matrix"]),
                }
                for row in section
            ],
        },
        "fixed_lattice": {
            "rank": 16,
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
            "mw_rank_below_factory_floor_rejected": counters[
                "mw_rank_below_factory_floor_rejected"
            ],
            "mod2_trivial_rejected": counters["mod2_trivial_rejected"],
            "high_mw_mod2_accepted_seeds": len(seed_records),
            "accepted_seed_mw_rank_distribution": {
                str(key): value for key, value in sorted(seed_mw_distribution.items())
            },
            "S4_section_embedding_orbits": len(orbit_records),
            "section_orbit_mw_rank_distribution": {
                str(key): value for key, value in sorted(orbit_mw_distribution.items())
            },
            "post_dedup_frame_mw_rank_distribution": {
                str(key): value
                for key, value in sorted(post_dedup_frame_mw_distribution.items())
            },
            "swap_moved_dimension_mod_2_distribution": {
                str(key): value
                for key, value in sorted(stabilizer_moved_distribution.items())
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
        "elkies-k3/scripts/enumerate_4d6_swap_fixed_high_mw_seed.sage"
    )
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("4D6 swap-fixed high-MW seed artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    accounting = result["accounting"]
    print(
        "4D6SWAP|seeds={}|section_orbits={}|surfaces={}|partners={}|frames={}|status=PASS_EXACT".format(
            accounting["high_mw_mod2_accepted_seeds"],
            accounting["S4_section_embedding_orbits"],
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
