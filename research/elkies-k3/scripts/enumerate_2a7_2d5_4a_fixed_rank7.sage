#!/usr/bin/env sage-python
"""Enumerate the determinant-5000 rank-seven family fixed by umbral 4A.

Let N be the Niemeier lattice with roots A7^2 D5^2 and let g be either of
the two 4A elements in the exact Dih_4 chamber section.  Their common fixed
lattice F has rank eight.  This script enumerates *all* primitive corank-one
sublattices K of F with det(K) <= 5000, computes M=K^perp in N, and applies
the ternary discriminant-form gate.

If H is the Gram matrix of F and a is a primitive integral covector, then

    K = ker(a:F -> Z),    det(K) = det(H) * a H^(-1) a^t.

Thus the declared family is a complete finite dual-lattice enumeration.  It
is not a census of rank-seven sublattices outside Fix(4A).  The whole Dih_4
section fixes F pointwise, so every retained K has literal 2B, 2C, and 4A
stabilizer elements before any expensive ambient-orbit computation.

status: EXACT_4A_FIXED_CORANK_ONE_FAMILY_DET_LE_5000
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from sage.all import (
    GF,
    Genus,
    QQ,
    QuadraticForm,
    ZZ,
    identity_matrix,
    lcm,
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
    / "artifacts/generated-results/elkies-k3-2a7-2d5-4a-fixed-rank7-v1.json"
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


def row_module_key(value):
    return matrix_key(value.row_module(ZZ).basis_matrix())


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


def discriminant_invariants(gram):
    diagonal = gram.smith_form()[0].diagonal()
    return [abs(int(entry)) for entry in diagonal if abs(int(entry)) > 1]


def primitive_closure_index(value):
    saturated = value.row_module(ZZ).saturation().basis_matrix()
    coordinates = value * saturated.pseudoinverse()
    assert all(entry.denominator() == 1 for entry in coordinates.list())
    return abs(int(matrix(ZZ, coordinates).det()))


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
    raise AssertionError("unexpected action order")


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
        elif rank >= 4 and count == 2 * rank * (rank - 1):
            label = f"D{rank}"
        elif (rank, count) in {(6, 72), (7, 126), (8, 240)}:
            label = f"E{rank}"
        else:
            raise AssertionError((rank, count))
        components.append({"type": label, "rank": rank, "signed_root_count": count})
    components.sort(key=lambda row: (row["type"], row["rank"]))
    multiplicities = Counter(row["type"] for row in components)
    label = "+".join(
        f"{count if count > 1 else ''}{name}"
        for name, count in sorted(multiplicities.items())
    )
    rank = sum(row["rank"] for row in components)
    return {
        "root_type": label,
        "root_rank": rank,
        "signed_root_count": len(roots),
        "root_components": components,
        "mw_rank_for_rho_19": 17 - rank,
    }


def ternary_gate(frame, cache):
    determinant = int(frame.det())
    frame_key = normal_form_key(Genus(frame).discriminant_form())
    cache_key = json.dumps(frame_key, sort_keys=True)
    if cache_key not in cache:
        all_genera = genera((2, 1), determinant, even=True)
        matches = []
        for genus_index, genus in enumerate(all_genera):
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
        cache[cache_key] = {
            "frame_discriminant_form_normal_key": frame_key,
            "all_even_ternary_genera_at_determinant": len(all_genera),
            "matching_even_ternary_genera": len(matches),
            "ternary_genus_representatives": matches,
        }
    return cache[cache_key]


def build(catalog, umbral, determinant_bound):
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
            "class": row["class"],
            "order": row["order"],
            "matrix": matrix(ZZ, row["ambient_matrix"]),
        }
        for row in umbral["group_section"]
    ]
    assert Counter(row["class"] for row in section) == Counter(
        {"1A": 1, "2A": 1, "2B": 2, "2C": 2, "4A": 2}
    )
    four_a = [row["matrix"] for row in section if row["class"] == "4A"]
    assert len(four_a) == 2 and four_a[0] * four_a[1] == identity24
    fixed_bases = [
        (action - identity24).transpose().right_kernel_matrix().change_ring(ZZ)
        for action in four_a
    ]
    assert fixed_bases[0].row_module(ZZ) == fixed_bases[1].row_module(ZZ)
    fixed = fixed_bases[0].row_module(ZZ).basis_matrix()
    assert fixed.nrows() == 8 and primitive_closure_index(fixed) == 1
    fixed_gram = fixed * ambient * fixed.transpose()
    assert fixed_gram.det() == 4096

    section_restrictions = []
    for row in section:
        restriction = fixed.transpose().solve_right(
            (fixed * row["matrix"]).transpose()
        ).transpose()
        assert restriction == identity_matrix(ZZ, 8)
        section_restrictions.append(restriction)

    dual = fixed_gram.inverse()
    dual_scale = lcm(entry.denominator() for entry in dual.list())
    integral_dual_gram = matrix(ZZ, dual_scale * dual)
    dual_bound = int((determinant_bound * dual_scale) // fixed_gram.det())
    short = pari(integral_dual_gram).qfminim(dual_bound)
    covectors = set()
    signed_primitive = 0
    for column in matrix(ZZ, short[2].sage()).columns():
        for signed in (vector(ZZ, column), -vector(ZZ, column)):
            if math.gcd(*[abs(int(entry)) for entry in signed]) != 1:
                continue
            norm = int(signed * integral_dual_gram * signed)
            determinant = int(fixed_gram.det() * QQ(norm) / dual_scale)
            if determinant > determinant_bound:
                continue
            signed_primitive += 1
            value = tuple(map(int, signed))
            covectors.add(min(value, tuple(-entry for entry in value)))
    covectors = sorted(covectors)
    assert signed_primitive == 2 * len(covectors)

    candidates = []
    action_profiles = Counter()
    for covector in covectors:
        coefficient_kernel = matrix(ZZ, 1, 8, covector).right_kernel_matrix()
        auxiliary_basis = coefficient_kernel * fixed
        assert auxiliary_basis.nrows() == 7
        assert primitive_closure_index(auxiliary_basis) == 1
        auxiliary = auxiliary_basis * ambient * auxiliary_basis.transpose()
        dual_norm = int(
            vector(ZZ, covector)
            * integral_dual_gram
            * vector(ZZ, covector)
        )
        predicted_determinant = int(
            fixed_gram.det() * QQ(dual_norm) / dual_scale
        )
        assert auxiliary.det() == predicted_determinant <= determinant_bound
        complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
        complement = complement_basis * ambient * complement_basis.transpose()
        assert complement.nrows() == 17 and complement.det() == auxiliary.det()
        roots = root_type(complement)
        profile = []
        distinct_actions = set()
        for row in section:
            action = induced_action(complement_basis, row["matrix"])
            assert action * complement * action.transpose() == complement
            distinct_actions.add(matrix_key(action))
            entry = (
                row["class"],
                multiplicative_order(action),
                int(matrix(GF(2), action - identity_matrix(ZZ, 17)).rank()),
                int(17 - (action - identity_matrix(ZZ, 17)).rank()),
            )
            profile.append(entry)
        assert len(distinct_actions) == 8
        profile = tuple(profile)
        action_profiles[profile] += 1
        candidates.append(
            {
                "primitive_covector_up_to_sign": list(covector),
                "dual_integral_norm": dual_norm,
                "determinant": int(auxiliary.det()),
                "auxiliary_basis": auxiliary_basis,
                "auxiliary_gram": auxiliary,
                "complement_basis": complement_basis,
                "frame_gram": complement,
                "root_data": roots,
                "action_profile": profile,
            }
        )
    assert len(action_profiles) == 5
    assert sorted(action_profiles.values()) == [4, 12, 12, 12, 296]
    sorted_profiles = sorted(action_profiles)
    profile_id_by_key = {
        profile: f"A4FIX-MOD2-{index:02d}"
        for index, profile in enumerate(sorted_profiles, start=1)
    }
    for profile in sorted_profiles:
        assert [entry[0] for entry in profile] == [
            "1A", "2A", "2B", "2B", "2C", "2C", "4A", "4A"
        ]
        assert [entry[1] for entry in profile] == [1, 2, 2, 2, 2, 2, 4, 4]
        assert [entry[3] for entry in profile] == [17, 9, 9, 9, 9, 9, 1, 1]
        assert all(
            moved > 0
            for label, unused_order, moved, unused_fixed in profile
            if label in {"2B", "2C", "4A"}
        )

    auxiliary_classes = []
    for candidate_index, candidate in enumerate(candidates):
        auxiliary_class = None
        for existing in auxiliary_classes:
            if exact_isometric(candidate["auxiliary_gram"], existing["gram"]):
                auxiliary_class = existing
                break
        if auxiliary_class is None:
            auxiliary_class = {
                "gram": candidate["auxiliary_gram"],
                "candidate_indices": [],
                "frame_classes": [],
            }
            auxiliary_classes.append(auxiliary_class)
        auxiliary_class["candidate_indices"].append(candidate_index)
        frame_class = None
        for existing in auxiliary_class["frame_classes"]:
            if exact_isometric(candidate["frame_gram"], existing["gram"]):
                frame_class = existing
                break
        if frame_class is None:
            frame_class = {
                "gram": candidate["frame_gram"],
                "candidate_indices": [],
                "representative_index": candidate_index,
            }
            auxiliary_class["frame_classes"].append(frame_class)
        frame_class["candidate_indices"].append(candidate_index)

    ternary_cache = {}
    auxiliary_rows = []
    frame_rows = []
    for auxiliary_index, auxiliary_class in enumerate(auxiliary_classes, start=1):
        auxiliary_id = f"A4FIX-K{auxiliary_index:03d}"
        local_frame_ids = []
        for frame_class in auxiliary_class["frame_classes"]:
            frame_id = f"A4FIX-F{len(frame_rows) + 1:03d}"
            local_frame_ids.append(frame_id)
            representative = candidates[frame_class["representative_index"]]
            gate = ternary_gate(frame_class["gram"], ternary_cache)
            frame_rows.append(
                {
                    "frame_id": frame_id,
                    "auxiliary_id": auxiliary_id,
                    "determinant": int(frame_class["gram"].det()),
                    "gram": rows(frame_class["gram"]),
                    "gram_sha256": gram_digest(frame_class["gram"]),
                    "discriminant_invariants_greater_than_one": (
                        discriminant_invariants(frame_class["gram"])
                    ),
                    "discriminant_length": len(
                        discriminant_invariants(frame_class["gram"])
                    ),
                    "root_data": representative["root_data"],
                    "embedding_count_in_declared_family": len(
                        frame_class["candidate_indices"]
                    ),
                    "representative_auxiliary_basis_in_ambient": rows(
                        representative["auxiliary_basis"]
                    ),
                    "representative_complement_basis_in_ambient": rows(
                        representative["complement_basis"]
                    ),
                    "representative_covector_up_to_sign": representative[
                        "primitive_covector_up_to_sign"
                    ],
                    "k3_discriminant_gate": gate,
                }
            )
            for index in frame_class["candidate_indices"]:
                candidates[index]["frame_id"] = frame_id
        representative = candidates[auxiliary_class["candidate_indices"][0]]
        auxiliary_rows.append(
            {
                "auxiliary_id": auxiliary_id,
                "determinant": int(auxiliary_class["gram"].det()),
                "gram": rows(auxiliary_class["gram"]),
                "gram_sha256": gram_digest(auxiliary_class["gram"]),
                "discriminant_invariants_greater_than_one": (
                    discriminant_invariants(auxiliary_class["gram"])
                ),
                "discriminant_length": len(
                    discriminant_invariants(auxiliary_class["gram"])
                ),
                "embedding_count_in_declared_family": len(
                    auxiliary_class["candidate_indices"]
                ),
                "frame_ids": local_frame_ids,
                "representative_auxiliary_basis_in_ambient": rows(
                    representative["auxiliary_basis"]
                ),
            }
        )
        for index in auxiliary_class["candidate_indices"]:
            candidates[index]["auxiliary_id"] = auxiliary_id

    surface_rows_by_key = {}
    for frame in frame_rows:
        ns_key = negative_form_key(
            frame["k3_discriminant_gate"]["frame_discriminant_form_normal_key"]
        )
        for ternary in frame["k3_discriminant_gate"][
            "ternary_genus_representatives"
        ]:
            key = {
                "ns_discriminant_form_key": ns_key,
                "transcendental_gram": ternary["gram"],
            }
            encoded = json.dumps(key, sort_keys=True, separators=(",", ":"))
            if encoded not in surface_rows_by_key:
                surface_rows_by_key[encoded] = {
                    "surface_gate_id": f"K3-A4FIX-{compact_digest(key)}",
                    "surface_key": key,
                    "determinant": frame["determinant"],
                    "auxiliary_ids": [],
                    "frame_ids": [],
                }
            surface = surface_rows_by_key[encoded]
            if frame["auxiliary_id"] not in surface["auxiliary_ids"]:
                surface["auxiliary_ids"].append(frame["auxiliary_id"])
            if frame["frame_id"] not in surface["frame_ids"]:
                surface["frame_ids"].append(frame["frame_id"])
    surface_rows = sorted(
        surface_rows_by_key.values(),
        key=lambda row: (row["determinant"], row["surface_gate_id"]),
    )

    candidate_rows = [
        {
            "candidate_index_zero_based": index,
            "primitive_covector_up_to_sign": candidate[
                "primitive_covector_up_to_sign"
            ],
            "dual_integral_norm": candidate["dual_integral_norm"],
            "determinant": candidate["determinant"],
            "auxiliary_id": candidate["auxiliary_id"],
            "frame_id": candidate["frame_id"],
            "root_type": candidate["root_data"]["root_type"],
            "root_rank": candidate["root_data"]["root_rank"],
            "mw_rank_for_rho_19": candidate["root_data"][
                "mw_rank_for_rho_19"
            ],
            "mod2_action_profile_id": profile_id_by_key[
                candidate["action_profile"]
            ],
        }
        for index, candidate in enumerate(candidates)
    ]
    determinant_distribution = Counter(row["determinant"] for row in candidate_rows)
    mw_distribution = Counter(row["mw_rank_for_rho_19"] for row in candidate_rows)
    root_distribution = Counter(row["root_type"] for row in candidate_rows)
    compatible_frames = sum(
        row["k3_discriminant_gate"]["matching_even_ternary_genera"] > 0
        for row in frame_rows
    )
    assert len(candidate_rows) == 336
    assert determinant_distribution == Counter({4096: 312, 2048: 24})
    assert mw_distribution == Counter({17: 304, 15: 16, 13: 16})
    assert root_distribution == Counter({"0": 304, "2A1": 16, "4A1": 16})
    assert len(auxiliary_rows) == 3 and len(frame_rows) == 5
    assert all(row["discriminant_length"] == 7 for row in auxiliary_rows)
    assert all(row["discriminant_length"] == 7 for row in frame_rows)
    assert compatible_frames == 0 and not surface_rows

    return {
        "schema": "elkies-k3.2a7-2d5-4a-fixed-rank7.v1",
        "status": "PASS_EXACT_4A_FIXED_CORANK_ONE_RANK7_FAMILY_DET_LE_5000",
        "proof_scope": {
            "proved": (
                "All primitive rank-seven corank-one sublattices of the common "
                "4A fixed lattice with determinant at most the declared bound are "
                "enumerated by the exact dual-lattice identity. Every K is primitive "
                "in N(A7^2 D5^2), its complement and roots are exact, and the full "
                "Dih_4 chamber section fixes K pointwise. The induced 2B, 2C, and "
                "4A complement actions are nontrivial modulo two. The exact "
                "ternary gate rejects the whole family: every discriminant group "
                "has length seven, whereas a rank-three T has length at most three."
            ),
            "not_proved": (
                "This is not a census of rank-seven auxiliaries outside Fix(4A), "
                "nor are the 336 embeddings deduplicated under the full Weyl group. "
                "Ternary representatives prove existence in matching genera but do "
                "not enumerate every ternary lattice class in those genera."
            ),
        },
        "parameters": {
            "ambient": "2A7_2D5",
            "determinant_bound": determinant_bound,
            "family": "primitive corank-one sublattices of Fix(4A)",
        },
        "fixed_lattice": {
            "rank": 8,
            "determinant": int(fixed_gram.det()),
            "basis_in_ambient": rows(fixed),
            "gram": rows(fixed_gram),
            "primitive_in_ambient": True,
            "two_4A_elements_are_inverse_with_common_fixed_lattice": True,
            "whole_Dih_4_section_fixes_pointwise": True,
            "integral_dual_scale": int(dual_scale),
            "integral_dual_gram": rows(integral_dual_gram),
            "integral_dual_norm_bound": dual_bound,
            "determinant_identity": "det(K)=det(F)*a*H^(-1)*a^t",
        },
        "stabilizer_action": {
            "literal_section_order": 8,
            "class_distribution": dict(
                sorted(Counter(row["class"] for row in section).items())
            ),
            "complement_action_profiles": [
                {
                    "profile_id": profile_id_by_key[profile],
                    "embedding_count": action_profiles[profile],
                    "actions": [
                        {
                            "class": label,
                            "order": order,
                            "moved_dimension_mod_2": moved,
                            "fixed_dimension_over_Q": fixed_dimension,
                        }
                        for label, order, moved, fixed_dimension in profile
                    ],
                }
                for profile in sorted_profiles
            ],
            "requested_classes_pass_rank_g_minus_identity_mod_2_gate": [
                "2B",
                "2C",
                "4A",
            ],
        },
        "accounting": {
            "signed_primitive_dual_covectors": signed_primitive,
            "primitive_covectors_up_to_sign": len(covectors),
            "primitive_rank7_embeddings_in_declared_family": len(candidate_rows),
            "determinant_distribution": {
                str(key): value for key, value in sorted(determinant_distribution.items())
            },
            "root_type_distribution": dict(sorted(root_distribution.items())),
            "mw_rank_distribution": {
                str(key): value for key, value in sorted(mw_distribution.items())
            },
            "auxiliary_isometry_classes": len(auxiliary_rows),
            "frame_isometry_classes_within_auxiliary_classes": len(frame_rows),
            "distinct_mod2_action_profiles": len(action_profiles),
            "frames_passing_ternary_genus_gate": compatible_frames,
            "surface_rows_after_T_NS_first_dedup": len(surface_rows),
            "family_rejected_by_discriminant_length": True,
            "common_discriminant_length": 7,
        },
        "surface_gate_rows_T_NS_first": surface_rows,
        "auxiliary_classes": auxiliary_rows,
        "frame_classes": frame_rows,
        "candidates": candidate_rows,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--umbral", type=Path, default=UMBRAL)
    parser.add_argument("--determinant-bound", type=int, default=5000)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    assert arguments.determinant_bound == 5000
    result = build(
        json.loads(arguments.catalog.read_text()),
        json.loads(arguments.umbral.read_text()),
        arguments.determinant_bound,
    )
    result["inputs"] = {
        str(arguments.catalog.resolve().relative_to(ROOT)): digest(arguments.catalog),
        str(arguments.umbral.resolve().relative_to(ROOT)): digest(arguments.umbral),
    }
    result["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/enumerate_2a7_2d5_4a_fixed_rank7.sage"
    )
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("4A-fixed rank-seven artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    accounting = result["accounting"]
    print(
        "A4FIXED|embeddings={}|K={}|frames={}|surfaces={}|mw={}|status=PASS_EXACT".format(
            accounting["primitive_rank7_embeddings_in_declared_family"],
            accounting["auxiliary_isometry_classes"],
            accounting["frame_isometry_classes_within_auxiliary_classes"],
            accounting["surface_rows_after_T_NS_first_dedup"],
            accounting["mw_rank_distribution"],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
