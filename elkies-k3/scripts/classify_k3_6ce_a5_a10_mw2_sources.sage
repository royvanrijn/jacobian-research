#!/usr/bin/env sage-python
"""Classify the determinant-384 A5+A10/MW2 pole-[0,1] source cut."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCES = ROOT / "artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-semistable-mw0-2-sources-large-a-partner1-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-source-isometries-v1.json"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
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
    return sorted(result, key=len)


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


def marking_profile(source, expected_component_ranks=(5, 10)):
    root_rank = int(source["root_rank"])
    root = matrix(QQ, source["root_adapted_gram"])[:root_rank, :root_rank]
    components = connected_components(root)
    if tuple(map(len, components)) != tuple(expected_component_ranks):
        raise ArithmeticError("unexpected A-type components")
    basis = source["pole_audit"]["basis"]
    labels = [vector(QQ, section["simple_root_pairings"]) for section in basis]
    depths = [[], []]
    cross_correction = QQ(0)
    for component in components:
        block = root.matrix_from_rows_and_columns(component, component)
        inverse = block.inverse()
        local = [
            vector(QQ, [label[index] for index in component]) for label in labels
        ]
        for basis_index in range(2):
            depths[basis_index].append(
                depth_from_correction(
                    len(component), local[basis_index] * inverse * local[basis_index]
                )
            )
        cross_correction += local[0] * inverse * local[1]

    height = matrix(QQ, source["mw_height_gram"])
    coordinates = matrix(
        QQ, [section["mw_quotient_coordinates"] for section in basis]
    )
    basis_height = coordinates * height * coordinates.transpose()
    poles = [int(section["pole_order"]) for section in basis]
    smooth_intersection = (
        QQ(2 + sum(poles)) - cross_correction - basis_height[0, 1]
    )
    if smooth_intersection.denominator() != 1 or smooth_intersection < 0:
        raise ArithmeticError("invalid section intersection")
    return {
        "basis_pole_profile": poles,
        "support_profile_I_order_left_depth_right_depth": [
            [len(component) + 1, depths[0][index], depths[1][index]]
            for index, component in enumerate(components)
        ],
        "basis_height_gram": rational_rows(basis_height),
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
    candidates = [
        row
        for row in source_payload["sources"]
        if row["source"]["mw_rank_for_rho_19"] == 2
        and row["source"]["root_type"] == "A10+A5"
        and row["source"]["support_count"] == 2
        and row["source"]["torsion"] == 1
        and row["source"]["pole_audit"]["basis"] is not None
        and [
            section["pole_order"]
            for section in row["source"]["pole_audit"]["basis"]
        ]
        == [0, 1]
    ]
    candidates.sort(key=lambda row: row["source_id"])
    if len(candidates) != 24:
        raise ArithmeticError("A5+A10 pole-[0,1] MW2 cut changed")

    classes = []
    for row in candidates:
        gram = matrix(ZZ, row["source"]["gram"])
        match = None
        for class_row in classes:
            isometry = integral_isometry(gram, class_row["representative_gram"])
            if isometry is not None:
                match = (class_row, isometry)
                break
        if match is None:
            class_row = {
                "representative_source_id": row["source_id"],
                "representative_gram": gram,
                "members": [],
            }
            classes.append(class_row)
            isometry = identity_matrix(ZZ, 17)
        else:
            class_row, isometry = match
        profile = marking_profile(row["source"])
        class_row["members"].append(
            {
                "source_id": row["source_id"],
                "source_gram_sha256": row["source"]["gram_sha256"],
                "source_to_representative_isometry": rows(isometry),
                "isometry_determinant": int(isometry.det()),
                "marking_profile": profile,
            }
        )

    output_classes = []
    total_profiles = 0
    for class_index, class_row in enumerate(classes, start=1):
        profile_members = {}
        for member in class_row["members"]:
            key = json.dumps(member["marking_profile"], sort_keys=True)
            profile_members.setdefault(key, []).append(member["source_id"])
        profiles = [
            {"marking_profile": json.loads(key), "source_ids": source_ids}
            for key, source_ids in sorted(profile_members.items())
        ]
        total_profiles += len(profiles)
        output_classes.append(
            {
                "class_index": class_index,
                "representative_source_id": class_row["representative_source_id"],
                "member_count": len(class_row["members"]),
                "marked_basis_profiles": profiles,
                "members": class_row["members"],
            }
        )

    payload = {
        "schema": "elkies-k3.k3-6ce-a5-a10-mw2-source-isometries.v1",
        "status": "PASS_EXACT_INTEGRAL_ISOMETRY_AND_MARKING_PROFILE_CLASSIFICATION",
        "ideal_cut": {
            "surface_id": "K3-6ce16abb9de3c7c5",
            "target_mw_rank": 15,
            "root_type": "A10+A5",
            "mw_rank": 2,
            "semistable": True,
            "support_count": 2,
            "torsion_order": 1,
            "complete_basis_pole_profile": [0, 1],
        },
        "accounting": {
            "reduced_gram_rows": len(candidates),
            "integral_isometry_classes": len(classes),
            "class_sizes": sorted(len(row["members"]) for row in classes),
            "selected_physical_basis_profiles": total_profiles,
        },
        "classes": output_classes,
        "proof_boundary": {
            "proved": (
                "Every selected reduced Gram is placed in an exact integral-isometry "
                "class, and each displayed physical basis profile has exact A-component "
                "depths, height Gram, pole orders, and section intersection."
            ),
            "not_proved": (
                "No rational marking, equation, moduli-curve identification, neighbour "
                "corridor, or equivalence of distinct marked profiles is asserted."
            ),
        },
        "inputs": {relative(SOURCES): digest(SOURCES)},
        "reproduce": "sage -python elkies-k3/scripts/classify_k3_6ce_a5_a10_mw2_sources.sage",
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if not output_path.exists() or output_path.read_text() != serialized:
            raise SystemExit(f"stale artifact: {output_path}")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "K36CEA5A10|"
        f"rows={len(candidates)}|classes={len(classes)}|profiles={total_profiles}|status=PASS"
    )


if __name__ == "__main__":
    main()
