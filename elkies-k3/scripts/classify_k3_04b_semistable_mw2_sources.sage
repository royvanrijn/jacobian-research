#!/usr/bin/env sage-python
"""Classify the determinant-500 semistable pole-[0,0] MW2 source cut."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCES = ROOT / "artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-semistable-mw2-source-isometries-v1.json"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def connected_components(gram):
    unseen = set(range(gram.nrows()))
    result = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        todo = [first]
        component = []
        while todo:
            node = todo.pop()
            component.append(node)
            for other in tuple(unseen):
                if gram[node, other]:
                    unseen.remove(other)
                    todo.append(other)
        result.append(sorted(component))
    return result


def depth_from_correction(rank, correction):
    order = rank + 1
    matches = [
        depth
        for depth in range(order // 2 + 1)
        if QQ(depth * (order - depth)) / order == correction
    ]
    if len(matches) != 1:
        raise ArithmeticError("component correction has ambiguous local depth")
    return matches[0]


def marking_profile(source):
    root_rank = int(source["root_rank"])
    root = matrix(QQ, source["root_adapted_gram"])[:root_rank, :root_rank]
    components = connected_components(root)
    if sorted(map(len, components)) != [3, 4, 8]:
        raise ArithmeticError("unexpected A3+A4+A8 components")
    components.sort(key=len)
    labels = [
        vector(QQ, section["simple_root_pairings"])
        for section in source["pole_audit"]["basis"]
    ]
    depth_profiles = [[], []]
    cross_correction = QQ(0)
    for component in components:
        block = root.matrix_from_rows_and_columns(component, component)
        inverse = block.inverse()
        local = [
            vector(QQ, [label[index] for index in component])
            for label in labels
        ]
        corrections = [value * inverse * value for value in local]
        for basis_index in range(2):
            depth_profiles[basis_index].append(
                depth_from_correction(len(component), corrections[basis_index])
            )
        cross_correction += local[0] * inverse * local[1]

    height = matrix(QQ, source["mw_height_gram"])
    coordinates = matrix(
        QQ,
        [
            section["mw_quotient_coordinates"]
            for section in source["pole_audit"]["basis"]
        ],
    )
    basis_height = coordinates * height * coordinates.transpose()
    smooth_intersection = QQ(2) - cross_correction - basis_height[0, 1]
    if smooth_intersection.denominator() != 1 or smooth_intersection < 0:
        raise ArithmeticError("invalid section intersection")

    def candidate(swapped):
        order = (1, 0) if swapped else (0, 1)
        supports = [
            [
                len(component) + 1,
                depth_profiles[order[0]][index],
                depth_profiles[order[1]][index],
            ]
            for index, component in enumerate(components)
        ]
        ordered_height = basis_height.matrix_from_rows_and_columns(order, order)
        key = (tuple(map(tuple, supports)), tuple(ordered_height.list()))
        return key, supports, ordered_height

    unused_key, supports, ordered_height = min(
        candidate(False), candidate(True), key=lambda row: row[0]
    )
    return {
        "support_profile_I_order_left_depth_right_depth": supports,
        "basis_height_gram": rational_rows(ordered_height),
        "component_cross_correction": str(cross_correction),
        "required_smooth_pair_intersection": int(smooth_intersection),
    }


def integral_isometry(left, right):
    raw = pari(left).qfisom(pari(right))
    if raw == 0:
        return None
    candidate = matrix(ZZ, raw)
    for value in (candidate, candidate.transpose()):
        if value * left * value.transpose() == right:
            return value
        if value * right * value.transpose() == left:
            inverse = value.inverse()
            if inverse.denominator() == 1:
                return inverse.change_ring(ZZ)
    raise ArithmeticError("PARI returned an unrecognized isometry orientation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    source_payload = json.loads(SOURCES.read_text())
    candidates = []
    for row in source_payload["sources"]:
        source = row["source"]
        if (
            source["mw_rank_for_rho_19"] == 2
            and source["root_type"] == "A3+A4+A8"
            and source["support_count"] == 3
            and source["torsion"] == 1
            and [
                section["pole_order"] for section in source["pole_audit"]["basis"]
            ]
            == [0, 0]
        ):
            candidates.append(row)
    candidates.sort(key=lambda row: row["source_id"])
    if len(candidates) != 9:
        raise ArithmeticError("semistable pole-[0,0] MW2 cut changed")

    representative = matrix(ZZ, candidates[0]["source"]["gram"])
    members = []
    profile_members = {}
    for row in candidates:
        source_gram = matrix(ZZ, row["source"]["gram"])
        isometry = (
            identity_matrix(ZZ, 17)
            if row is candidates[0]
            else integral_isometry(source_gram, representative)
        )
        if isometry is None or abs(isometry.det()) != 1:
            raise ArithmeticError("ideal MW2 cut split into multiple isometry classes")
        profile = marking_profile(row["source"])
        profile_key = json.dumps(profile, sort_keys=True)
        profile_members.setdefault(profile_key, []).append(row["source_id"])
        members.append(
            {
                "source_id": row["source_id"],
                "source_gram_sha256": row["source"]["gram_sha256"],
                "source_to_representative_isometry": rows(isometry),
                "isometry_determinant": int(isometry.det()),
                "marking_profile": profile,
            }
        )

    marked_profiles = [
        {"marking_profile": json.loads(key), "source_ids": source_ids}
        for key, source_ids in sorted(profile_members.items())
    ]
    if len(marked_profiles) != 2 or sorted(
        len(row["source_ids"]) for row in marked_profiles
    ) != [3, 6]:
        raise ArithmeticError("selected physical-basis profile partition changed")
    equation_conditions = {
        (
            tuple(
                map(
                    tuple,
                    row["marking_profile"][
                        "support_profile_I_order_left_depth_right_depth"
                    ],
                )
            ),
            row["marking_profile"]["required_smooth_pair_intersection"],
        )
        for row in marked_profiles
    }
    if len(equation_conditions) != 1:
        raise ArithmeticError("marked profiles require different equation conditions")

    payload = {
        "schema": "elkies-k3.k3-04b-semistable-mw2-source-isometries.v1",
        "status": "PASS_EXACT_SINGLE_INTEGRAL_ISOMETRY_CLASS_TWO_BASIS_PROFILES",
        "ideal_cut": {
            "surface_id": "K3-04b86146cc6b284b",
            "root_type": "A3+A4+A8",
            "mw_rank": 2,
            "semistable": True,
            "support_count": 3,
            "torsion_order": 1,
            "complete_basis_pole_profile": [0, 0],
        },
        "accounting": {
            "reduced_gram_rows": len(candidates),
            "integral_isometry_classes": 1,
            "class_sizes": [len(candidates)],
            "selected_physical_basis_profiles": len(marked_profiles),
            "selected_physical_basis_profile_sizes": sorted(
                len(row["source_ids"]) for row in marked_profiles
            ),
        },
        "class": {
            "representative_source_id": candidates[0]["source_id"],
            "tested_equation_source_id": "K3-04b86146cc6b284b-S2021",
            "marked_basis_profiles": marked_profiles,
            "common_normalized_equation_conditions": {
                "support_profile_I_order_left_depth_right_depth": [
                    [4, 1, 1], [5, 0, 0], [9, 1, 3]
                ],
                "required_smooth_pair_intersection": 1,
            },
            "members": members,
        },
        "proof_boundary": {
            "proved": (
                "All nine reduced-Gram rows in the declared ideal MW2 cut are "
                "carried to one representative by displayed determinant-one integral "
                "isometries. Their selected physical pole-zero bases split into two "
                "exact profiles of sizes three and six; both profiles impose the same "
                "support depths and smooth section-pair intersection on an equation."
            ),
            "not_proved": (
                "The six-large-ambient source inventory is not the full fibration "
                "classification. Integral marking equivalence does not supply a "
                "Q-rational equation, rational sections, or a neighbour corridor."
            ),
        },
        "inputs": {relative(SOURCES): digest(SOURCES)},
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/classify_k3_04b_semistable_mw2_sources.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print("K304BMW2ISO|rows=9|classes=1|basis_profiles=2|status=PASS", flush=True)


if __name__ == "__main__":
    main()
