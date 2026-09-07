#!/usr/bin/env sage-python
"""Enumerate a symmetry-first high-MW seed shell inside Fix(2C).

This is a deterministic bounded discovery shell, not the complete set of
rank-seven sublattices of the rank-sixteen fixed lattice.  Select the first
2C element in the exact Dih_4 chamber section of N(A7^2 D5^2), compute an
integral LLL basis of its primitive fixed lattice, and enumerate all 7 of 16
coordinate direct summands.  Retain determinant at most 5000, discriminant
length at most three, and MW rank 12 through 17.  Close and canonicalize the
survivors under the full Dih_4 section before applying the ternary gate and
deduplicating first by (T,NS), then by auxiliary and frame isometry.

Every seed K is fixed pointwise by the selected 2C involution.  Its induced
action on M=K^perp is computed exactly and must remain nontrivial modulo two.

status: EXACT_DECLARED_2C_FIXED_COORDINATE_SEED_SHELL
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from sage.all import (
    GF,
    Genus,
    QQ,
    QuadraticForm,
    ZZ,
    identity_matrix,
    matrix,
    pari,
    vector,
)
from sage.quadratic_forms.genera.genus import genera


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
UMBRAL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-umbral-orbits-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-2a7-2d5-2c-fixed-high-mw-seed-v1.json"
)


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def rows_as_strings(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gram_digest(value):
    encoded = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(encoded.encode()).hexdigest()


def compact_digest(value, length=16):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()[:length]


def matrix_key(value):
    return tuple(map(int, value.list()))


def row_module_basis(value):
    return value.row_module(ZZ).basis_matrix()


def row_module_key(value):
    return matrix_key(row_module_basis(value))


def normal_form_key(discriminant_form):
    normal = discriminant_form.normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rows_as_strings(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def negative_form_key(normal):
    return {
        "invariants": normal["invariants"],
        "quadratic_gram": [
            [str(-QQ(entry)) for entry in row]
            for row in normal["quadratic_gram"]
        ],
        "value_module": normal["value_module"],
        "basis": "negative_of_source_frame_discriminant_normal_basis",
    }


def exact_isometric(left, right):
    return (
        left.nrows() == right.nrows()
        and left.det() == right.det()
        and pari(left).qfisom(pari(right)) != 0
    )


def primitive_closure_index(value):
    saturated = value.row_module(ZZ).saturation().basis_matrix()
    coordinates = value * saturated.pseudoinverse()
    assert all(entry.denominator() == 1 for entry in coordinates.list())
    return abs(int(matrix(ZZ, coordinates).det()))


def discriminant_invariants(gram):
    diagonal = gram.smith_form()[0].diagonal()
    return [abs(int(entry)) for entry in diagonal if abs(int(entry)) > 1]


def induced_action(complement, ambient_action):
    result = complement.transpose().solve_right(
        (complement * ambient_action).transpose()
    ).transpose()
    assert all(entry.denominator() == 1 for entry in result.list())
    result = matrix(ZZ, result)
    assert result * complement == complement * ambient_action
    return result


def multiplicative_order(value):
    identity = identity_matrix(ZZ, value.nrows())
    for order in range(1, 9):
        if value**order == identity:
            return order
    raise AssertionError("unexpected finite action order")


def signed_roots(gram):
    data = pari(gram).qfminim(2)
    half = [
        vector(ZZ, column)
        for column in matrix(ZZ, data[2].sage()).columns()
    ]
    result = half + [-root for root in half]
    assert len(result) == int(data[0])
    return result


def root_type(gram):
    roots = signed_roots(gram)
    if not roots:
        return {
            "root_type": "0",
            "root_rank": 0,
            "signed_root_count": 0,
            "root_components": [],
            "root_determinant": 1,
            "mw_rank_for_rho_19": 17,
        }
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
        component_roots = matrix(ZZ, [roots[index] for index in sorted(component)])
        rank = int(component_roots.rank())
        count = len(component)
        if count == rank * (rank + 1):
            label = f"A{rank}"
            determinant = rank + 1
        elif rank >= 4 and count == 2 * rank * (rank - 1):
            label = f"D{rank}"
            determinant = 4
        elif (rank, count) in {(6, 72), (7, 126), (8, 240)}:
            label = f"E{rank}"
            determinant = {6: 3, 7: 2, 8: 1}[rank]
        else:
            raise AssertionError((rank, count))
        components.append(
            {
                "type": label,
                "rank": rank,
                "signed_root_count": count,
                "determinant": determinant,
            }
        )
    components.sort(key=lambda row: (row["type"], row["rank"]))
    multiplicities = Counter(row["type"] for row in components)
    label = "+".join(
        f"{count if count > 1 else ''}{name}"
        for name, count in sorted(multiplicities.items())
    )
    rank = sum(row["rank"] for row in components)
    root_determinant = 1
    for component in components:
        root_determinant *= component["determinant"]
    return {
        "root_type": label,
        "root_rank": rank,
        "signed_root_count": len(roots),
        "root_components": components,
        "root_determinant": root_determinant,
        "mw_rank_for_rho_19": 17 - rank,
    }


def ternary_gate(frame, genera_cache, form_cache):
    determinant = int(frame.det())
    frame_key = normal_form_key(Genus(frame).discriminant_form())
    encoded_key = json.dumps(frame_key, sort_keys=True)
    if determinant not in genera_cache:
        genera_cache[determinant] = genera((2, 1), determinant, even=True)
    if encoded_key not in form_cache:
        matches = []
        for genus_index, genus in enumerate(genera_cache[determinant]):
            if normal_form_key(genus.discriminant_form()) != frame_key:
                continue
            representative = matrix(ZZ, genus.representative())
            assert representative.det() == -determinant
            assert QuadraticForm(QQ, representative).signature() == 1
            matches.append(
                {
                    "genus_index": genus_index,
                    "gram": rows(representative),
                    "local_symbols": [str(symbol) for symbol in genus.local_symbols()],
                }
            )
        form_cache[encoded_key] = {
            "frame_discriminant_form_normal_key": frame_key,
            "all_even_ternary_genera_at_determinant": len(
                genera_cache[determinant]
            ),
            "matching_even_ternary_genera": len(matches),
            "ternary_genus_representatives": matches,
        }
    return form_cache[encoded_key]


def section_is_group(section):
    keys = {matrix_key(row["matrix"]) for row in section}
    return all(
        matrix_key(left["matrix"] * right["matrix"]) in keys
        for left in section
        for right in section
    )


def find_or_add_isometry_class(classes, gram, record_index):
    for row in classes:
        if exact_isometric(gram, row["gram"]):
            row["record_indices"].append(record_index)
            return row
    row = {"gram": gram, "record_indices": [record_index]}
    classes.append(row)
    return row


def build(catalog, umbral, determinant_bound, minimum_mw_rank):
    assert catalog["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert umbral["schema"] == "elkies-k3.lattice-foundry-umbral-orbits.v1"
    ambient_row = next(
        row
        for row in catalog["rooted_niemeier_lattices"]
        if row["label"] == "2A7_2D5"
    )
    ambient = matrix(ZZ, ambient_row["gram"])
    identity24 = identity_matrix(ZZ, 24)
    section = [
        {
            "section_index_zero_based": index,
            "class": row["class"],
            "order": row["order"],
            "matrix": matrix(ZZ, row["ambient_matrix"]),
        }
        for index, row in enumerate(umbral["group_section"])
    ]
    assert len(section) == 8 and section_is_group(section)
    assert Counter(row["class"] for row in section) == Counter(
        {"1A": 1, "2A": 1, "2B": 2, "2C": 2, "4A": 2}
    )
    selected = next(row for row in section if row["class"] == "2C")
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
        assert moved_dimension > 0
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
                "selected_2C_moved_dimension_mod_2": moved_dimension,
            }
        )

    orbit_by_key = {}
    for seed_index, seed in enumerate(seed_records):
        images = {}
        for section_row in section:
            image = row_module_basis(seed["auxiliary_basis"] * section_row["matrix"])
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
            if row_module_key(auxiliary_basis * section_row["matrix"]) != canonical_key:
                continue
            action = induced_action(complement_basis, section_row["matrix"])
            assert action * frame * action.transpose() == frame
            stabilizer.append(
                {
                    "section_index_zero_based": section_row[
                        "section_index_zero_based"
                    ],
                    "class": section_row["class"],
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
            row["class"] == "2C" and row["moved_dimension_mod_2"] > 0
            for row in stabilizer
        )
        gate = ternary_gate(frame, genera_cache, form_cache)
        assert gate["matching_even_ternary_genera"] > 0
        orbit_records.append(
            {
                "orbit_id": f"2CFIX-O{orbit_index + 1:04d}",
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
                    "surface_id": f"K3-2CFIX-{compact_digest(key)}",
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
                        orbit_records[index]["orbit_id"]
                        for index in partner["record_indices"]
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
                        orbit_records[record_index]["orbit_id"]
                        for record_index in frame_class["record_indices"]
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
                    orbit_records[index]["orbit_id"]
                    for index in surface["orbit_indices"]
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
    selected_2c_moved_distribution = Counter(
        next(
            action["moved_dimension_mod_2"]
            for action in orbit["literal_section_stabilizer"]
            if action["class"] == "2C"
        )
        for orbit in orbit_records
    )
    assert counters["coordinate_subsets_tested"] == 11440
    assert counters["determinant_rejected"] == 1
    assert counters["discriminant_length_rejected"] == 106
    assert len(seed_records) == 97
    assert seed_mw_distribution == Counter({13: 77, 12: 16, 15: 3, 14: 1})
    assert len(orbit_records) == 97
    assert all(orbit["section_orbit_size"] == 4 for orbit in orbit_records)
    assert all(
        [action["class"] for action in orbit["literal_section_stabilizer"]]
        == ["1A", "2C"]
        for orbit in orbit_records
    )
    assert selected_2c_moved_distribution == Counter({7: 63, 8: 32, 6: 2})
    assert len(form_cache) == 73
    assert all(
        gate["matching_even_ternary_genera"] == 1
        for gate in form_cache.values()
    )
    assert len(surface_rows) == 73
    assert sum(len(row["partner_auxiliaries"]) for row in surface_rows) == 76
    assert sum(len(row["frames"]) for row in surface_rows) == 86
    assert sum(len(row["frames"]) > 1 for row in surface_rows) == 13
    assert post_dedup_frame_mw_distribution == Counter(
        {13: 66, 12: 16, 15: 3, 14: 1}
    )

    return {
        "schema": "elkies-k3.2a7-2d5-2c-fixed-high-mw-seed.v1",
        "status": "PASS_EXACT_DECLARED_2C_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST",
        "proof_scope": {
            "proved": (
                "All 7-of-16 coordinate direct summands of the pinned integral "
                "LLL basis of Fix(2C) are tested. Every retained K is primitive, "
                "has determinant at most 5000 and discriminant length at most three, "
                "and has an exactly computed complement of MW rank 12 through 17. "
                "The selected 2C action on every complement is nontrivial modulo two. "
                "Survivors are closed and canonicalized under the exact Dih_4 section, "
                "then deduplicated first by (T,NS), and only then by K and frame isometry."
            ),
            "not_proved": (
                "Coordinate direct summands of one deterministic LLL basis do not "
                "exhaust primitive rank-seven sublattices of Fix(2C), and Dih_4-section "
                "canonicalization is not the missing full Weyl quotient. Ternary genus "
                "representatives prove exact K3 realizability examples but do not "
                "enumerate every class inside each ternary genus."
            ),
        },
        "parameters": {
            "ambient": "2A7_2D5",
            "selected_section_index_zero_based": selected[
                "section_index_zero_based"
            ],
            "selected_class": "2C",
            "determinant_bound": determinant_bound,
            "discriminant_length_bound": 3,
            "minimum_mw_rank": minimum_mw_rank,
            "maximum_mw_rank": 17,
            "seed_language": "7-of-16 coordinate direct summands of pinned LLL fixed-lattice basis",
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
            "high_mw_mod2_accepted_seeds": len(seed_records),
            "accepted_seed_mw_rank_distribution": {
                str(key): value for key, value in sorted(seed_mw_distribution.items())
            },
            "Dih_4_section_embedding_orbits": len(orbit_records),
            "section_orbit_mw_rank_distribution": {
                str(key): value for key, value in sorted(orbit_mw_distribution.items())
            },
            "post_dedup_frame_mw_rank_distribution": {
                str(key): value
                for key, value in sorted(post_dedup_frame_mw_distribution.items())
            },
            "selected_2C_moved_dimension_mod_2_distribution": {
                str(key): value
                for key, value in sorted(selected_2c_moved_distribution.items())
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
    parser.add_argument("--umbral", type=Path, default=UMBRAL)
    parser.add_argument("--determinant-bound", type=int, default=5000)
    parser.add_argument("--minimum-mw-rank", type=int, default=12)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    assert arguments.determinant_bound == 5000
    assert arguments.minimum_mw_rank == 12
    result = build(
        json.loads(arguments.catalog.read_text()),
        json.loads(arguments.umbral.read_text()),
        arguments.determinant_bound,
        arguments.minimum_mw_rank,
    )
    result["inputs"] = {
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(arguments.catalog),
        str(arguments.umbral.resolve().relative_to(ROOT)): digest(arguments.umbral),
    }
    result["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
    )
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("2C-fixed high-MW seed artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    accounting = result["accounting"]
    print(
        "2CFIXED|seeds={}|section_orbits={}|surfaces={}|partners={}|frames={}|status=PASS_EXACT".format(
            accounting["high_mw_mod2_accepted_seeds"],
            accounting["Dih_4_section_embedding_orbits"],
            accounting["surface_classes_after_T_NS_first_dedup"],
            accounting["partner_auxiliary_isometry_classes_after_surface_dedup"],
            accounting["frame_isometry_classes_after_surface_dedup"],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
