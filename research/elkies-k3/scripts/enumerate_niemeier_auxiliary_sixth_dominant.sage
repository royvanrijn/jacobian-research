#!/usr/bin/env sage
"""Enumerate Weyl-dominant norm-12 sixth-generator candidates.

status: ACTIVE_SEARCH
claim: For every D5 anchor, enumerate exactly the norm-12 vectors modulo the
  residual Weyl group by nonnegative Dynkin labels and, where necessary, one
  additional integral coordinate in the root-free fixed line.
inputs: artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json,
  artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json
outputs: artifacts/generated-results/elkies-k3-niemeier-auxiliary-sixth-dominant.json
supersedes/superseded-by: none
"""

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

from sage.all import QQ, RR, ZZ, ceil, floor, matrix, pari, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
ANCHORS = ROOT / "artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-niemeier-auxiliary-sixth-dominant.json"


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def signed_roots(gram):
    minimum_data = pari(gram).qfminim(2)
    representatives = matrix(ZZ, minimum_data[2].sage()).transpose()
    result = [vector(ZZ, row) for row in representatives.rows()]
    return result + [-root for root in result]


def root_components(gram, roots):
    unseen = set(range(len(roots)))
    result = []
    while unseen:
        component = {min(unseen)}
        frontier = list(component)
        unseen.difference_update(component)
        while frontier:
            current = frontier.pop()
            neighbours = {
                index
                for index in unseen
                if roots[current] * gram * roots[index] != 0
            }
            component.update(neighbours)
            unseen.difference_update(neighbours)
            frontier.extend(sorted(neighbours))
        result.append([roots[index] for index in sorted(component)])
    return result


def simple_roots(gram, roots):
    for trial in range(1, 100):
        chamber_vector = vector(
            ZZ,
            [
                (index + 1) ** 2 + trial * (index + 1) + trial**2
                for index in range(gram.nrows())
            ],
        )
        values = [root * gram * chamber_vector for root in roots]
        if all(value != 0 for value in values):
            break
    else:
        raise AssertionError("failed to choose a regular chamber vector")
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
    assert result.rank() == result.nrows()
    return result


def dominant_labels_up_to(cartan_inverse, bound):
    rank = cartan_inverse.nrows()
    current = [ZZ(0)] * rank
    result = []

    def extend(index, norm):
        if index == rank:
            result.append((tuple(map(int, current)), norm))
            return
        cross = sum(
            current[previous] * cartan_inverse[previous, index]
            for previous in range(index)
        )
        coefficient = 0
        while True:
            new_norm = (
                norm
                + 2 * coefficient * cross
                + coefficient * coefficient * cartan_inverse[index, index]
            )
            if new_norm > bound:
                break
            current[index] = coefficient
            extend(index + 1, new_norm)
            coefficient += 1
        current[index] = 0

    extend(0, QQ(0))
    return result


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

catalog_payload = json.loads(CATALOG.read_text())
anchor_payload = json.loads(ANCHORS.read_text())
gram_by_label = {
    entry["label"]: matrix(ZZ, entry["gram"])
    for entry in catalog_payload["rooted_niemeier_lattices"]
}

completed = []
total_candidates = 0
for anchor in anchor_payload["anchors"]:
    label = anchor["niemeier"]
    gram = gram_by_label[label]
    d5_basis = matrix(ZZ, anchor["D5_basis_in_ambient"])
    assert d5_basis * gram * d5_basis.transpose() == matrix(
        ZZ,
        [[2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, -1, -1], [0, 0, -1, 2, 0], [0, 0, -1, 0, 2]],
    )
    ambient_roots = signed_roots(gram)
    orthogonal_roots = [
        root
        for root in ambient_roots
        if root * gram * d5_basis.transpose() == 0
    ]
    components = root_components(gram, orthogonal_roots)
    simple_components = [simple_roots(gram, roots) for roots in components]
    residual_rank = sum(component.nrows() for component in simple_components)
    anchor_key = f"{label}:{anchor['anchor_index']}"
    assert residual_rank in (18, 19)

    residual_simple = matrix(
        ZZ,
        [row for component in simple_components for row in component.rows()],
    )
    complement_basis = (d5_basis * gram).right_kernel_matrix()
    complement_gram = complement_basis * gram * complement_basis.transpose()
    pairing = complement_basis * gram * residual_simple.transpose()
    assert pairing.nrows() == 19 and pairing.ncols() == residual_rank
    assert pairing.rank() == residual_rank
    if residual_rank == 19:
        coordinate_map = pairing
    else:
        for coordinate in range(19):
            extra = matrix(ZZ, 19, 1, [int(index == coordinate) for index in range(19)])
            candidate_map = pairing.augment(extra)
            if candidate_map.rank() == 19:
                coordinate_map = candidate_map
                break
        else:
            raise AssertionError("failed to complement the Dynkin-label map")
    coordinate_inverse = coordinate_map.inverse()
    coordinate_norm = (
        coordinate_inverse
        * complement_gram
        * coordinate_inverse.transpose()
    )
    residual_cartan_inverse = (
        residual_simple * gram * residual_simple.transpose()
    ).inverse()
    if residual_rank == 19:
        assert coordinate_norm == residual_cartan_inverse
    else:
        # Minimizing over the final fixed-line coordinate gives exactly the
        # norm of the projection to the residual root span.
        root_block = coordinate_norm[:18, :18]
        cross = coordinate_norm[:18, 18]
        fixed_norm = coordinate_norm[18, 18]
        assert root_block - cross * cross.transpose() / fixed_norm == residual_cartan_inverse

    component_candidates = []
    component_summaries = []
    for roots, component in zip(components, simple_components):
        cartan = component * gram * component.transpose()
        candidates = dominant_labels_up_to(cartan.inverse(), QQ(12))
        component_candidates.append(candidates)
        norm_histogram = Counter(str(norm) for unused_labels, norm in candidates)
        component_summaries.append(
            {
                "rank": component.nrows(),
                "signed_root_count": len(roots),
                "dominant_labels_norm_at_most_12": len(candidates),
                "norm_histogram": dict(sorted(norm_histogram.items())),
            }
        )

    candidate_rows = []
    pre_integrality_count = 0
    nonprimitive_rejected = 0
    for choice in itertools.product(*component_candidates):
        projected_norm = sum(norm for unused_labels, norm in choice)
        if projected_norm > 12:
            continue
        labels = vector(
            ZZ,
            [entry for component_labels, unused_norm in choice for entry in component_labels],
        )
        if residual_rank == 19:
            if projected_norm != 12:
                continue
            fixed_coordinates = [None]
        else:
            root_block = coordinate_norm[:18, :18]
            cross = coordinate_norm[:18, 18]
            fixed_norm = coordinate_norm[18, 18]
            linear = (labels * cross)[0]
            constant = labels * root_block * labels
            minimum = constant - linear * linear / fixed_norm
            assert minimum == projected_norm
            radius = RR((QQ(12) - minimum) / fixed_norm).sqrt()
            centre = -linear / fixed_norm
            lower = floor(RR(centre) - radius) - 2
            upper = ceil(RR(centre) + radius) + 2
            fixed_coordinates = [
                value
                for value in range(int(lower), int(upper) + 1)
                if fixed_norm * value * value + 2 * linear * value + constant == 12
            ]
        for fixed_coordinate in fixed_coordinates:
            pre_integrality_count += 1
            extended = labels if fixed_coordinate is None else vector(ZZ, list(labels) + [fixed_coordinate])
            complement_coordinates = extended * coordinate_inverse
            if not all(entry.denominator() == 1 for entry in complement_coordinates):
                continue
            complement_coordinates = vector(ZZ, complement_coordinates)
            ambient_vector = complement_coordinates * complement_basis
            assert ambient_vector * gram * ambient_vector == 12
            assert ambient_vector * gram * d5_basis.transpose() == 0
            first_six = d5_basis.stack(matrix(ZZ, [ambient_vector]))
            smith = first_six.smith_form()[0]
            smith_invariants = tuple(
                abs(int(smith[index, index])) for index in range(6)
            )
            if smith_invariants != (1, 1, 1, 1, 1, 1):
                nonprimitive_rejected += 1
                continue
            remaining_roots = sum(
                root * gram * ambient_vector == 0 for root in orthogonal_roots
            )
            row = {
                "dynkin_labels": list(map(int, labels)),
                "ambient_vector": list(map(int, ambient_vector)),
                "first_six_primitive": True,
                "remaining_signed_roots_after_D5_and_sixth": remaining_roots,
            }
            if fixed_coordinate is not None:
                row["fixed_line_coordinate"] = fixed_coordinate
            candidate_rows.append(row)

    total_candidates += len(candidate_rows)
    remaining_histogram = Counter(
        row["remaining_signed_roots_after_D5_and_sixth"] for row in candidate_rows
    )
    completed.append(
        {
            "anchor": anchor_key,
            "niemeier": label,
            "component_type": anchor["component_type"],
            "residual_root_rank": residual_rank,
            "fixed_space_rank": 19 - residual_rank,
            "residual_root_components": component_summaries,
            "dominant_norm12_label_combinations_before_integrality": pre_integrality_count,
            "weyl_dominant_norm12_vectors": len(candidate_rows),
            "nonprimitive_norm12_vectors_rejected": nonprimitive_rejected,
            "remaining_signed_root_histogram": {
                str(key): value for key, value in sorted(remaining_histogram.items())
            },
            "candidates": candidate_rows,
        }
    )

assert len(completed) == 16

result = {
    "schema": "elkies-k3.niemeier-auxiliary-sixth-dominant.v2",
    "status": "PASS_EXACT_ALL_D5_ANCHORS_SIXTH_DOMINANT_ENUMERATION",
    "classification_scope": {
        "proved": (
            "For every D5 anchor, the listed "
            "norm-12 vectors are an exhaustive set modulo the residual Weyl "
            "group, obtained from every nonnegative Dynkin-label tuple of "
            "projected norm at most 12 and every exact fixed-line solution "
            "that belongs to the Niemeier lattice. Candidates whose first "
            "six auxiliary vectors are not primitive are rejected exactly."
        ),
        "not_proved": (
            "The remaining anchor-stabilizer quotient, the seventh auxiliary "
            "generator, primitivity, and "
            "rootless complement classification remain open."
        ),
    },
    "accounting": {
        "D5_anchor_orbits": 16,
        "sixth_generator_anchor_enumerations_completed": len(completed),
        "rank18_residual_root_anchors": sum(
            row["residual_root_rank"] == 18 for row in completed
        ),
        "weyl_dominant_norm12_candidates": total_candidates,
        "nonprimitive_sixth_candidates_rejected": sum(
            row["nonprimitive_norm12_vectors_rejected"] for row in completed
        ),
    },
    "completed_anchors": completed,
}

payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    assert arguments.output.read_text() == payload
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(payload)

print(
    "NIEMEIERAUX6|anchors=16|completed={}|rank18={}|candidates={}|"
    "seventh_complete=0|status=PASS".format(
        len(completed),
        sum(row["residual_root_rank"] == 18 for row in completed),
        total_candidates,
    )
)
