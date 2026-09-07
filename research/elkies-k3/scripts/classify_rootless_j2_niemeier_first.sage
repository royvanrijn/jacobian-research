#!/usr/bin/env sage
"""Classify rootless J2 frames by exhaustive Niemeier-first enumeration.

status: ACTIVE_PROOF
claim: Starting from the pinned rank-7 auxiliary lattice, enumerate every
  primitive embedding into every Niemeier lattice up to enough Weyl symmetry
  to obtain every rootless orthogonal complement, then deduplicate those
  complements by exact integral isometry.
inputs: elkies-k3/data/lattice/rootless_j2_auxiliary_rank7_gram.txt,
  elkies-k3/data/lattice/rank17_gram.txt,
  artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json,
  artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json,
  artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json,
  artifacts/generated-results/elkies-k3-niemeier-auxiliary-sixth-dominant.json
outputs: artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json
supersedes/superseded-by: none
"""

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from sage.all import (
    CartanMatrix,
    Genus,
    QQ,
    ZZ,
    ceil,
    diagonal_matrix,
    floor,
    matrix,
    pari,
    vector,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUXILIARY = ROOT / "elkies-k3/data/lattice/rootless_j2_auxiliary_rank7_gram.txt"
PUBLISHED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
ALTERNATE = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
ANCHORS = ROOT / "artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json"
SIXTH = ROOT / "artifacts/generated-results/elkies-k3-niemeier-auxiliary-sixth-dominant.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json"

D5_GRAM = CartanMatrix(["D", 5])
SEVENTH_PAIRINGS = vector(ZZ, [-1, 0, 0, 0, 1, -6])
SEVENTH_PROJECTED_NORM_BUDGET = QQ(79) / 4


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def gram_sha256(value):
    payload = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


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


def positive_dominant_labels_up_to(cartan_inverse, bound):
    """All strictly positive integral Dynkin labels within an exact bound."""
    rank = cartan_inverse.nrows()
    current = [ZZ(1)] * rank
    result = []

    def extend(index, norm):
        if index == rank:
            result.append((tuple(map(int, current)), norm))
            return
        cross = sum(
            current[previous] * cartan_inverse[previous, index]
            for previous in range(index)
        )
        coefficient = 1
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
        current[index] = 1

    extend(0, QQ(0))
    return result


def rational_ldl(value):
    """Return exact L,D with value = L*diag(D)*L^t and unit lower L."""
    size = value.nrows()
    lower = matrix(QQ, size, size)
    diagonal = [QQ(0)] * size
    for row in range(size):
        lower[row, row] = 1
        diagonal[row] = value[row, row] - sum(
            lower[row, index] ** 2 * diagonal[index]
            for index in range(row)
        )
        assert diagonal[row] > 0
        for later in range(row + 1, size):
            lower[later, row] = (
                value[later, row]
                - sum(
                    lower[later, index]
                    * diagonal[index]
                    * lower[row, index]
                    for index in range(row)
                )
            ) / diagonal[row]
    assert lower * diagonal_matrix(diagonal) * lower.transpose() == value
    return lower, diagonal


def exact_ceil_sqrt(value):
    """A floating-point-free ceiling of sqrt(value) for nonnegative QQ."""
    assert value >= 0
    numerator = int(value.numerator())
    denominator = int(value.denominator())
    result = math.isqrt(numerator // denominator)
    while result * result * denominator < numerator:
        result += 1
    return result


def shifted_ellipsoid_shell(quadratic, centre, norm):
    """Enumerate integral b with (b-centre) Q (b-centre)^t = norm.

    Exact LDL branch bounds are deliberately rounded outwards using integer
    arithmetic.  The final equality and every partial inequality are exact.
    """
    dimension = quadratic.nrows()
    if dimension == 0:
        return [tuple()] if norm == 0 else []
    lower, diagonal = rational_ldl(quadratic)
    current = [ZZ(0)] * dimension
    shifted = [QQ(0)] * dimension
    result = []

    def extend(index, remaining):
        if index < 0:
            if remaining == 0:
                result.append(tuple(map(int, current)))
            return
        tail = sum(
            lower[later, index] * shifted[later]
            for later in range(index + 1, dimension)
        )
        local_centre = centre[index] - tail
        radius_squared = remaining / diagonal[index]
        radius = exact_ceil_sqrt(radius_squared)
        first = int(floor(local_centre)) - radius - 1
        last = int(ceil(local_centre)) + radius + 1
        for coordinate in range(first, last + 1):
            local_shift = QQ(coordinate) - local_centre
            term = diagonal[index] * local_shift**2
            if term > remaining:
                continue
            current[index] = coordinate
            shifted[index] = QQ(coordinate) - centre[index]
            extend(index - 1, remaining - term)

    extend(dimension - 1, norm)
    for candidate in result:
        difference = vector(QQ, candidate) - centre
        assert difference * quadratic * difference == norm
    return result


def standard_auxiliary_gram():
    result = matrix(ZZ, 7, 7)
    result[:5, :5] = D5_GRAM
    result[5, 5] = 12
    result[6, 6] = 24
    for index, pairing in enumerate(SEVENTH_PAIRINGS[:5]):
        result[index, 6] = pairing
        result[6, index] = pairing
    result[5, 6] = result[6, 5] = SEVENTH_PAIRINGS[5]
    return result


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

catalog_payload = json.loads(CATALOG.read_text())
anchor_payload = json.loads(ANCHORS.read_text())
sixth_payload = json.loads(SIXTH.read_text())
alternate_payload = json.loads(ALTERNATE.read_text())

pinned_auxiliary = load_matrix(AUXILIARY)
standard_auxiliary = standard_auxiliary_gram()
standard_isometry = matrix(
    ZZ, pari(standard_auxiliary).qfisom(pari(pinned_auxiliary)).sage()
)
assert abs(standard_isometry.det()) == 1
assert (
    standard_isometry.transpose()
    * pinned_auxiliary
    * standard_isometry
    == standard_auxiliary
)

published = load_matrix(PUBLISHED)
alternate = matrix(ZZ, alternate_payload["rootless_frame"])
target_genus = Genus(published)
assert Genus(alternate) == target_genus

gram_by_label = {
    entry["label"]: matrix(ZZ, entry["gram"])
    for entry in catalog_payload["rooted_niemeier_lattices"]
}
d5_by_anchor = {
    f"{entry['niemeier']}:{entry['anchor_index']}": matrix(
        ZZ, entry["D5_basis_in_ambient"]
    )
    for entry in anchor_payload["anchors"]
}
roots_by_label = {
    label: signed_roots(gram) for label, gram in gram_by_label.items()
}

accounting_by_anchor = []
rootless_embeddings = []
total_sixth = 0
total_feasible_sixth = 0
total_positive_label_combinations = 0
total_ellipsoid_solutions = 0
total_integral_sevenths = 0
total_nonprimitive_rejected = 0

for anchor_row in sixth_payload["completed_anchors"]:
    label = anchor_row["niemeier"]
    anchor_key = anchor_row["anchor"]
    gram = gram_by_label[label]
    d5_basis = d5_by_anchor[anchor_key]
    ambient_roots = roots_by_label[label]
    anchor_counts = Counter()
    fixed_dimension_histogram = Counter()

    for sixth_row in anchor_row["candidates"]:
        total_sixth += 1
        anchor_counts["sixth_candidates"] += 1
        sixth_vector = vector(ZZ, sixth_row["ambient_vector"])
        first_six = d5_basis.stack(matrix(ZZ, [sixth_vector]))
        orthogonal_roots = [
            root
            for root in ambient_roots
            if root * gram * first_six.transpose() == 0
        ]
        components = root_components(gram, orthogonal_roots)
        simple_components = [
            simple_roots(gram, component) for component in components
        ]
        simple_basis = matrix(
            ZZ,
            [row for component in simple_components for row in component.rows()],
        )
        residual_rank = simple_basis.nrows()
        fixed_dimension = 18 - residual_rank
        assert 0 <= fixed_dimension <= 6
        fixed_dimension_histogram[fixed_dimension] += 1

        component_labels = [
            positive_dominant_labels_up_to(
                (component * gram * component.transpose()).inverse(),
                SEVENTH_PROJECTED_NORM_BUDGET,
            )
            for component in simple_components
        ]
        positive_choices = []
        for choice in itertools.product(*component_labels):
            projected_norm = sum(norm for unused_labels, norm in choice)
            if projected_norm <= SEVENTH_PROJECTED_NORM_BUDGET:
                positive_choices.append((choice, projected_norm))
        if not positive_choices:
            continue

        total_feasible_sixth += 1
        anchor_counts["positive_label_feasible_sixth"] += 1
        anchor_counts["positive_label_combinations"] += len(positive_choices)
        total_positive_label_combinations += len(positive_choices)

        constraint_map = (gram * first_six.transpose()).augment(
            gram * simple_basis.transpose()
        )
        extra_coordinates = []
        for coordinate in range(24):
            if constraint_map.ncols() == 24:
                break
            unit = matrix(
                ZZ,
                24,
                1,
                [int(index == coordinate) for index in range(24)],
            )
            candidate = constraint_map.augment(unit)
            if candidate.rank() > constraint_map.rank():
                constraint_map = candidate
                extra_coordinates.append(coordinate)
        assert constraint_map.nrows() == constraint_map.ncols() == 24
        assert len(extra_coordinates) == fixed_dimension
        constraint_inverse = constraint_map.inverse()
        coordinate_norm = (
            constraint_inverse * gram * constraint_inverse.transpose()
        )

        prefix_size = 6 + residual_rank
        fixed_quadratic = coordinate_norm[prefix_size:, prefix_size:]
        fixed_inverse = (
            fixed_quadratic.inverse()
            if fixed_dimension
            else matrix(QQ, 0, 0)
        )
        for choice, projected_norm in positive_choices:
            labels = vector(
                ZZ,
                [entry for component, unused_norm in choice for entry in component],
            )
            assert all(label_value > 0 for label_value in labels)
            prefix = vector(QQ, list(SEVENTH_PAIRINGS) + list(labels))
            prefix_norm = prefix * coordinate_norm[:prefix_size, :prefix_size] * prefix
            if fixed_dimension:
                linear = prefix * coordinate_norm[:prefix_size, prefix_size:]
                centre = -linear * fixed_inverse
                minimum_norm = prefix_norm - linear * fixed_inverse * linear
            else:
                centre = vector(QQ, [])
                minimum_norm = prefix_norm
            assert minimum_norm == QQ(17) / 4 + projected_norm
            remaining_norm = QQ(24) - minimum_norm
            fixed_solutions = shifted_ellipsoid_shell(
                fixed_quadratic, centre, remaining_norm
            )
            total_ellipsoid_solutions += len(fixed_solutions)
            anchor_counts["exact_ellipsoid_solutions"] += len(fixed_solutions)

            for fixed_coordinates in fixed_solutions:
                coordinates = vector(
                    QQ,
                    list(SEVENTH_PAIRINGS)
                    + list(labels)
                    + list(fixed_coordinates),
                )
                seventh_rational = coordinates * constraint_inverse
                if not all(entry.denominator() == 1 for entry in seventh_rational):
                    anchor_counts["nonintegral_sevenths_rejected"] += 1
                    continue
                total_integral_sevenths += 1
                anchor_counts["integral_sevenths"] += 1
                seventh = vector(ZZ, seventh_rational)
                assert seventh * gram * seventh == 24
                assert seventh * gram * first_six.transpose() == SEVENTH_PAIRINGS
                assert seventh * gram * simple_basis.transpose() == labels

                auxiliary_basis = first_six.stack(matrix(ZZ, [seventh]))
                smith = auxiliary_basis.smith_form()[0]
                invariants = tuple(
                    abs(int(smith[index, index])) for index in range(7)
                )
                if invariants != (1, 1, 1, 1, 1, 1, 1):
                    total_nonprimitive_rejected += 1
                    anchor_counts["nonprimitive_auxiliaries_rejected"] += 1
                    continue
                assert auxiliary_basis * gram * auxiliary_basis.transpose() == standard_auxiliary

                complement_basis = (auxiliary_basis * gram).right_kernel_matrix()
                complement_gram = complement_basis * gram * complement_basis.transpose()
                assert complement_basis.nrows() == 17
                assert complement_gram.det() == 948
                assert int(pari(complement_gram).qfminim(2)[0]) == 0
                assert Genus(complement_gram) == target_genus
                anchor_counts["primitive_rootless_embeddings"] += 1
                rootless_embeddings.append(
                    {
                        "anchor": anchor_key,
                        "niemeier": label,
                        "auxiliary_basis": auxiliary_basis,
                        "complement_basis": complement_basis,
                        "complement_gram": complement_gram,
                        "fixed_dimension": fixed_dimension,
                        "positive_dynkin_labels": list(map(int, labels)),
                        "sixth_vector": sixth_vector,
                        "seventh_vector": seventh,
                    }
                )

    accounting_by_anchor.append(
        {
            "anchor": anchor_key,
            "niemeier": label,
            "fixed_dimension_histogram": {
                str(key): value for key, value in sorted(fixed_dimension_histogram.items())
            },
            **{key: int(value) for key, value in sorted(anchor_counts.items())},
        }
    )
    print(
        "NIEMEIERJ2ANCHOR|anchor={}|sixth={}|feasible={}|labels={}|rootless={}".format(
            anchor_key,
            anchor_counts["sixth_candidates"],
            anchor_counts["positive_label_feasible_sixth"],
            anchor_counts["positive_label_combinations"],
            anchor_counts["primitive_rootless_embeddings"],
        ),
        flush=True,
    )

classes = []
for embedding in rootless_embeddings:
    gram = embedding["complement_gram"]
    matched = None
    for class_row in classes:
        if pari(class_row["gram"]).qfisom(pari(gram)) != 0:
            matched = class_row
            break
    if matched is None:
        matched = {
            "gram": gram,
            "embedding_count_in_enumerated_cover": 0,
            "provenance_counts": Counter(),
            "representative": embedding,
        }
        classes.append(matched)
    matched["embedding_count_in_enumerated_cover"] += 1
    matched["provenance_counts"][embedding["anchor"]] += 1

class_rows = []
for index, class_row in enumerate(classes, start=1):
    gram = class_row["gram"]
    representative = class_row["representative"]
    minimum_four = pari(gram).qfminim(4)
    published_match = pari(published).qfisom(pari(gram)) != 0
    alternate_match = pari(alternate).qfisom(pari(gram)) != 0
    assert not (published_match and alternate_match)
    automorphism_order = int(pari(gram).qfauto()[0])
    class_rows.append(
        {
            "class_index": index,
            "gram": rows(gram),
            "gram_sha256": gram_sha256(gram),
            "determinant": int(gram.det()),
            "minimum": 4,
            "norm4_unoriented_pairs": int(minimum_four[0]) // 2,
            "automorphism_group_order": automorphism_order,
            "matches_published_R17": bool(published_match),
            "matches_alternate_Q80": bool(alternate_match),
            "embedding_count_in_enumerated_cover": class_row[
                "embedding_count_in_enumerated_cover"
            ],
            "niemeier_provenance_counts": dict(
                sorted(class_row["provenance_counts"].items())
            ),
            "representative_embedding": {
                "anchor": representative["anchor"],
                "niemeier": representative["niemeier"],
                "fixed_dimension": representative["fixed_dimension"],
                "positive_dynkin_labels": representative["positive_dynkin_labels"],
                "auxiliary_basis_in_ambient": rows(representative["auxiliary_basis"]),
                "complement_basis_in_ambient": rows(representative["complement_basis"]),
                "sixth_vector": list(map(int, representative["sixth_vector"])),
                "seventh_vector": list(map(int, representative["seventh_vector"])),
            },
        }
    )

assert sum(row["matches_published_R17"] for row in class_rows) == 1
assert sum(row["matches_alternate_Q80"] for row in class_rows) == 1

result = {
    "schema": "elkies-k3.rootless-j2-niemeier-first.v1",
    "status": "PASS_COMPLETE_ROOTLESS_J2_CLASSIFICATION",
    "classification_scope": {
        "proved": (
            "Every primitive embedding of the pinned auxiliary rank-7 lattice "
            "in a Niemeier lattice is covered up to the sequential D5-anchor, "
            "residual-Weyl sixth-vector, and residual-Weyl seventh-vector "
            "reductions. A dominant seventh vector has rootless complement "
            "exactly when all residual Dynkin labels are positive. Every such "
            "positive-label tuple and every integral point on its remaining "
            "exact rational ellipsoid is enumerated. The resulting saturated "
            "rank-17 complements are deduplicated by exact integral isometry."
        ),
        "embedding_orbit_boundary": (
            "The enumerated embedding representatives are a complete cover, "
            "but duplicate representatives under anchor stabilizers are not "
            "deduplicated; embedding_count_in_enumerated_cover is therefore "
            "not asserted to be a full-automorphism embedding-orbit count."
        ),
    },
    "standard_auxiliary": {
        "gram": rows(standard_auxiliary),
        "pinned_gram": rows(pinned_auxiliary),
        "standard_isometry_columns_in_pinned_basis": rows(standard_isometry),
        "determinant": int(standard_auxiliary.det()),
        "seventh_pairings_with_D5_and_sixth": list(map(int, SEVENTH_PAIRINGS)),
        "seventh_residual_norm_budget": str(SEVENTH_PROJECTED_NORM_BUDGET),
    },
    "accounting": {
        "rooted_niemeier_classes": 23,
        "leech_excluded_by_roots": 1,
        "D5_admissible_niemeier_classes": 13,
        "D5_anchor_orbits": 16,
        "weyl_dominant_primitive_sixth_candidates": total_sixth,
        "positive_label_feasible_sixth_candidates": total_feasible_sixth,
        "positive_label_combinations": total_positive_label_combinations,
        "exact_fixed_ellipsoid_solutions": total_ellipsoid_solutions,
        "integral_seventh_vectors": total_integral_sevenths,
        "nonprimitive_full_auxiliaries_rejected": total_nonprimitive_rejected,
        "primitive_rootless_embeddings_in_cover": len(rootless_embeddings),
        "rootless_complement_isometry_classes": len(class_rows),
    },
    "anchor_accounting": accounting_by_anchor,
    "rootless_classes": class_rows,
}

payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    assert arguments.output.read_text() == payload
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(payload)

print(
    "NIEMEIERJ2|sixth={}|feasible={}|labels={}|ellipsoid={}|integral={}|"
    "primitive_rootless={}|classes={}|published={}|alternate={}|status=PASS".format(
        total_sixth,
        total_feasible_sixth,
        total_positive_label_combinations,
        total_ellipsoid_solutions,
        total_integral_sevenths,
        len(rootless_embeddings),
        len(class_rows),
        sum(row["matches_published_R17"] for row in class_rows),
        sum(row["matches_alternate_Q80"] for row in class_rows),
    )
)
