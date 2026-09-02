#!/usr/bin/env sage-python
"""Classify the determinant-720 ideal-cut sources up to integral isometry.

The large prescribed-root inventory merges equal deterministic reduced Grams,
not full integral-isometry orbits.  For the equation-first ideal cut used in
the search—semistable MW2, at most three supports, and a complete pole-[0,0]
basis—this script performs the missing exact pairwise isometry test.  It also
reconstructs the selected physical basis's component depths, height Gram,
cross correction, and required smooth intersection, so equation charts are
deduplicated only when their actual marking profile agrees.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import re
from pathlib import Path

from sage.all import QQ, ZZ, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
SOURCES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
)
POLES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-source-poles-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-ideal-source-isometries-v1.json"
)

engine_path = ROOT / "elkies-k3/scripts/exact_neighbor_engine.sage"
engine = {"__file__": str(engine_path), "__name__": "golay_ideal_source_engine"}
exec(compile(engine_path.read_text(), str(engine_path), "exec"), engine)
deterministic_simple_roots = engine["deterministic_simple_roots"]


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value) -> list[list[str]]:
    return [[str(entry) for entry in row] for row in value.rows()]


def component_ranks(root_type: str) -> list[int]:
    answer = []
    for term in root_type.split("+"):
        match = re.fullmatch(r"(?:(\d+))?A(\d+)", term)
        if match is None:
            raise ValueError("ideal-cut source is not semistable")
        answer.extend([int(match.group(2))] * int(match.group(1) or 1))
    return answer


def connected_components(cartan) -> list[list[int]]:
    unseen = set(range(cartan.nrows()))
    answer = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        stack = [first]
        component = []
        while stack:
            left = stack.pop()
            component.append(left)
            adjacent = [
                right for right in sorted(unseen) if cartan[left, right] != 0
            ]
            for right in adjacent:
                unseen.remove(right)
                stack.append(right)
        answer.append(sorted(component))
    return answer


def ordered_components(cartan, ranks: list[int]) -> list[list[int]]:
    available = connected_components(cartan)
    answer = []
    for rank in ranks:
        match = next(
            (component for component in available if len(component) == rank), None
        )
        if match is None:
            raise ArithmeticError("Cartan components disagree with root type")
        answer.append(match)
        available.remove(match)
    if available:
        raise ArithmeticError("unused Cartan component")
    return answer


def depth_from_correction(rank: int, correction) -> int:
    order = rank + 1
    matches = [
        depth
        for depth in range(order // 2 + 1)
        if QQ(depth * (order - depth)) / order == correction
    ]
    if len(matches) != 1:
        raise ArithmeticError("component correction has ambiguous depth")
    return matches[0]


def marking_profile(source: dict, audit: dict) -> dict:
    frame = matrix(QQ, source["gram"])
    simple, unused_positive, cartan_rows = deterministic_simple_roots(frame)
    cartan = matrix(QQ, cartan_rows)
    ranks = component_ranks(source["root_type"])
    components = ordered_components(cartan, ranks)
    pairings = [
        vector(QQ, section["simple_root_pairings"])
        for section in audit["basis"]
    ]
    depth_profiles = [[], []]
    cross = QQ(0)
    for rank, component in zip(ranks, components):
        block = cartan.matrix_from_rows_and_columns(component, component)
        inverse = block.inverse()
        local = [
            vector(QQ, [pairing[index] for index in component])
            for pairing in pairings
        ]
        corrections = [value * inverse * value for value in local]
        for basis_index in range(2):
            depth_profiles[basis_index].append(
                depth_from_correction(rank, corrections[basis_index])
            )
        cross += local[0] * inverse * local[1]

    height = matrix(QQ, audit["height_gram"])
    coordinates = matrix(
        QQ, [section["free_mw_coordinates"] for section in audit["basis"]]
    )
    basis_height = coordinates * height * coordinates.transpose()
    required_intersection = QQ(2) - cross - basis_height[0, 1]
    if required_intersection.denominator() != 1 or required_intersection < 0:
        raise ArithmeticError("invalid physical basis intersection")
    def candidate(swapped):
        order = (1, 0) if swapped else (0, 1)
        support_profile = sorted(
            [
                (
                    rank + 1,
                    depth_profiles[order[0]][index],
                    depth_profiles[order[1]][index],
                )
                for index, rank in enumerate(ranks)
            ]
        )
        ordered_height = basis_height.matrix_from_rows_and_columns(order, order)
        key = (
            tuple(support_profile),
            tuple(ordered_height.list()),
            cross,
            required_intersection,
        )
        return key, support_profile, ordered_height, swapped

    unused_key, support_profile, ordered_height, swapped = min(
        candidate(False), candidate(True), key=lambda row: row[0]
    )
    return {
        "support_profile_I_order_left_depth_right_depth": [
            list(row) for row in support_profile
        ],
        "basis_height_gram": rational_rows(ordered_height),
        "component_cross_correction": str(cross),
        "required_smooth_pair_intersection": int(required_intersection),
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
            return value.inverse().change_ring(ZZ)
    raise ArithmeticError("PARI returned an unrecognized isometry orientation")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    source_payload = json.loads(SOURCES.read_text())
    pole_payload = json.loads(POLES.read_text())
    audits = {row["source_id"]: row["audit"] for row in pole_payload["audits"]}

    ideal_rows = []
    for row in source_payload["sources"]:
        source = row["source"]
        audit = audits[row["source_id"]]
        if (
            int(source["mw_rank_for_rho_19"]) == 2
            and int(source["support_count"]) <= 3
            and all(component["type"].startswith("A") for component in source["root_components"])
            and audit.get("basis_sorted_pole_profile") == [0, 0]
        ):
            ideal_rows.append(row)
    if len(ideal_rows) != 48:
        raise ArithmeticError("ideal-cut row count changed")

    by_root_type = defaultdict(list)
    for row in ideal_rows:
        by_root_type[row["source"]["root_type"]].append(row)

    classes = []
    for root_type, source_rows in sorted(by_root_type.items()):
        representatives = []
        members = []
        for source_row in sorted(source_rows, key=lambda row: row["source_id"]):
            source_gram = matrix(ZZ, source_row["source"]["gram"])
            for class_index, representative in enumerate(representatives):
                isometry = integral_isometry(source_gram, representative)
                if isometry is not None:
                    members[class_index].append((source_row, isometry))
                    break
            else:
                representatives.append(source_gram)
                members.append(
                    [(source_row, identity_matrix(ZZ, source_gram.nrows()))]
                )

        for representative, class_members in zip(representatives, members):
            representative_row = class_members[0][0]
            profiles = []
            serialized_members = []
            for source_row, isometry in class_members:
                profile = marking_profile(
                    source_row["source"], audits[source_row["source_id"]]
                )
                profiles.append(profile)
                serialized_members.append(
                    {
                        "source_id": source_row["source_id"],
                        "source_gram_sha256": source_row["source"]["gram_sha256"],
                        "source_to_representative_isometry": rows(isometry),
                        "isometry_determinant": int(isometry.det()),
                        "marking_profile": profile,
                    }
                )
            canonical_profile = profiles[0]
            if any(profile != canonical_profile for profile in profiles):
                raise ArithmeticError(
                    "integrally isometric sources have different marking profiles"
                )
            classes.append(
                {
                    "class_id": f"G720-I{len(classes)+1:03d}",
                    "representative_source_id": representative_row["source_id"],
                    "root_type": root_type,
                    "mw_rank": 2,
                    "basis_pole_profile": [0, 0],
                    "reduced_gram_row_count": len(class_members),
                    "marking_profile": canonical_profile,
                    "members": serialized_members,
                }
            )

    if len(classes) != 3 or any(row["reduced_gram_row_count"] not in (4, 9, 35) for row in classes):
        raise ArithmeticError("ideal-cut isometry classification changed")
    output = {
        "schema": "elkies-k3.golay-det720-ideal-source-isometries.v1",
        "status": "PASS_EXACT_IDEAL_CUT_INTEGRAL_ISOMETRY_AND_MARKING_CLASSES",
        "ideal_cut": {
            "mw_rank": 2,
            "semistable": True,
            "maximum_support_count": 3,
            "complete_basis_pole_profile": [0, 0],
        },
        "accounting": {
            "reduced_gram_rows": len(ideal_rows),
            "integral_isometry_classes": len(classes),
            "class_sizes": [row["reduced_gram_row_count"] for row in classes],
        },
        "classes": classes,
        "proof_boundary": {
            "proved": (
                "All 48 reduced-Gram rows satisfying the ideal cut are partitioned by "
                "explicit integral isometries. Every displayed matrix has determinant "
                "plus or minus one and carries its source Gram to the class representative. "
                "The selected physical basis profile is recomputed exactly for every row."
            ),
            "not_proved": (
                "This classifies only the ideal pole-[0,0] semistable cut, not all 4,823 "
                "source rows. Equal marking profiles identify the same normalized equation "
                "conditions; no Q-rational equation or neighbour corridor is inferred."
            ),
        },
        "inputs": {relative(SOURCES): digest(SOURCES), relative(POLES): digest(POLES)},
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/classify_golay_det720_ideal_source_isometries.sage"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("Golay determinant-720 ideal source classes are stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "GOLAY720IDEALISO|rows=48|classes=3|sizes={}|status=PASS".format(
            [row["reduced_gram_row_count"] for row in classes]
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
