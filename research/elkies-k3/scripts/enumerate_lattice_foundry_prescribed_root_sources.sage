#!/usr/bin/env sage-python
"""Enumerate low-MW foundry sources directly in Niemeier lattices.

For each selected foundry Neron--Severi class this script fixes one exact
rank-seven Nishiyama auxiliary whose root system is D5, standardizes it as

    D5 simple roots + one primitive D5-orthogonal vector + one final vector,

and enumerates primitive embeddings of that ordered auxiliary in selected
rooted Niemeier lattices.  The sixth vector is enumerated modulo the residual
Weyl group by dominant Dynkin labels.  For the seventh vector, zero Dynkin
labels are prescribed first so that their induced Dynkin subgraph has source
root rank 15--17 (MW rank 2--0); only then are positive labels and the exact
fixed-space ellipsoid enumerated.

Every emitted source is an exact saturated rank-17 orthogonal complement,
has exact ADE/MW data, and lies in the target frame genus. Exact repetitions
with the same deterministic reduced Gram are merged. Distinct reduced Grams
are deliberately not asserted to be distinct integral-isometry classes;
PARI's general rank-17 isometry search is not used as a discovery gate.
Sharing the auxiliary identifies the same NS/T class, but does not construct
a marked elliptic-neighbour corridor, a rational marking, or an equation.

The practical default is a quick determinant-948 ``NS0001`` smoke window in
``N(D16+E8)``.  The H3 height-form control requires ``--all-ambients``; use
``--all-ns --all-ambients`` only for the full database window.  Any explicit
candidate limits make the output a deterministic prefix search and are
recorded in the proof boundary.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from sage.all import (
    CartanMatrix,
    Genus,
    QQ,
    RR,
    ZZ,
    block_diagonal_matrix,
    ceil,
    diagonal_matrix,
    floor,
    identity_matrix,
    lcm,
    matrix,
    pari,
    vector,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATABASE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
ANCHORS = ROOT / "artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json"
H3_FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-prescribed-root-sources-v1.json"
)
D5_GRAM = CartanMatrix(["D", 5])
DESCRIPTION = __doc__


_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())


def rows(value):
    return [list(map(int, row)) for row in matrix(ZZ, value).rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in matrix(QQ, value).rows()]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def gram_digest(value):
    payload = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def reduced_gram(value):
    """LLL-reduce a positive Gram and retain the integral row change."""
    value = matrix(ZZ, value)
    change = matrix(ZZ, pari(value).qflllgram()).transpose()
    result = change * value * change.transpose()
    assert abs(change.det()) == 1 and result.det() == value.det()
    return result, change


def rational_form_isometric(left, right):
    """Test a small rational positive form by one common integral scaling."""
    left = matrix(QQ, left)
    right = matrix(QQ, right)
    if left.nrows() != right.nrows():
        return False
    scale = lcm(entry.denominator() for entry in list(left) + list(right))
    return (
        pari((scale * left).change_ring(ZZ)).qfisom(
            pari((scale * right).change_ring(ZZ))
        )
        != 0
    )


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(entry) for entry in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def primitive_rows(value):
    value = matrix(ZZ, value)
    smith = value.smith_form()[0]
    return tuple(
        abs(int(smith[index, index])) for index in range(value.nrows())
    ) == (1,) * value.nrows()


def complete_primitive_rows(value, ambient_rank):
    """Complete primitive rows to a deterministic unimodular square matrix."""
    value = matrix(ZZ, value)
    rank = value.nrows()
    smith, left, right = value.smith_form()
    assert smith == left * value * right
    assert smith[:rank, :rank] == identity_matrix(ZZ, rank)
    assert not any(smith[:, rank:].list())
    completion = (
        block_diagonal_matrix(
            left.inverse().change_ring(ZZ),
            identity_matrix(ZZ, ambient_rank - rank),
        )
        * right.inverse().change_ring(ZZ)
    )
    assert abs(completion.det()) == 1
    assert completion[:rank] == value
    return completion


def canonical_sign(value):
    value = vector(ZZ, value)
    first = next((entry for entry in value if entry), ZZ(0))
    return -value if first < 0 else value


def standard_d5_basis(auxiliary):
    simple, unused_positive, unused_cartan = deterministic_simple_roots(auxiliary)
    if simple.nrows() != 5:
        raise ValueError("foundry auxiliary does not have root rank five")
    for permutation in itertools.permutations(range(5)):
        candidate = simple[list(permutation)]
        if candidate * auxiliary * candidate.transpose() == D5_GRAM:
            if not primitive_rows(candidate):
                raise ValueError("auxiliary D5 root lattice is not primitive")
            return candidate
    raise ValueError("foundry auxiliary roots are not D5")


def standardize_auxiliary(auxiliary, sixth_norm_max):
    """Return a canonical D5/v6/v7 basis, or None above the declared v6 bound."""
    auxiliary = matrix(ZZ, auxiliary)
    d5 = standard_d5_basis(auxiliary)
    orthogonal = (d5 * auxiliary).right_kernel_matrix()
    assert orthogonal.nrows() == 2
    orthogonal_gram = orthogonal * auxiliary * orthogonal.transpose()
    minimum = pari(orthogonal_gram).qfminim(sixth_norm_max)
    candidates = []
    for column in matrix(ZZ, minimum[2].sage()).columns():
        for sign in (1, -1):
            sixth = canonical_sign(sign * vector(ZZ, column) * orthogonal)
            first_six = d5.stack(matrix(ZZ, [sixth]))
            if not primitive_rows(first_six):
                continue
            candidates.append(
                (
                    int(sixth * auxiliary * sixth),
                    tuple(map(abs, sixth)),
                    tuple(map(int, sixth)),
                )
            )
    if not candidates:
        return None
    unused_norm, unused_absolute, sixth_tuple = min(set(candidates))
    sixth = vector(ZZ, sixth_tuple)
    first_six = d5.stack(matrix(ZZ, [sixth]))
    basis = complete_primitive_rows(first_six, 7)
    standardized = basis * auxiliary * basis.transpose()
    assert standardized[:5, :5] == D5_GRAM
    assert not any(standardized[5, :5])
    return {
        "basis": basis,
        "gram": standardized,
        "sixth_norm": int(standardized[5, 5]),
        "seventh_norm": int(standardized[6, 6]),
        "seventh_pairings": vector(ZZ, list(standardized.row(6)[:6])),
    }


def signed_roots(gram):
    minimum = pari(gram).qfminim(2)
    half = [
        vector(ZZ, column)
        for column in matrix(ZZ, minimum[2].sage()).columns()
    ]
    result = half + [-root for root in half]
    assert len(result) == int(minimum[0])
    return result


def connected_root_components(gram, roots):
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
        components.append([roots[index] for index in sorted(component)])
    return components


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
        raise AssertionError("failed to choose a regular residual chamber")
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


def residual_simple_data(gram, fixed, ambient_roots):
    orthogonal_roots = [
        root
        for root in ambient_roots
        if root * gram * fixed.transpose() == 0
    ]
    components = connected_root_components(gram, orthogonal_roots)
    simple_components = [simple_roots(gram, component) for component in components]
    simple = matrix(
        ZZ,
        [row for component in simple_components for row in component.rows()],
    )
    cartan = simple * gram * simple.transpose()
    return orthogonal_roots, simple_components, simple, cartan


def dominant_labels_up_to(cartan_inverse, bound, positive_indices=None):
    """Enumerate nonnegative labels, positive exactly on selected indices."""
    rank = cartan_inverse.nrows()
    if positive_indices is None:
        positive_indices = tuple(range(rank))
        allow_zero = True
    else:
        positive_indices = tuple(sorted(positive_indices))
        allow_zero = False
    current = [ZZ(0)] * rank
    result = []

    def extend(position, norm):
        if position == len(positive_indices):
            result.append((tuple(map(int, current)), norm))
            return
        index = positive_indices[position]
        cross = sum(
            current[previous] * cartan_inverse[previous, index]
            for previous in range(index)
        )
        coefficient = 0 if allow_zero else 1
        while True:
            new_norm = (
                norm
                + 2 * coefficient * cross
                + coefficient * coefficient * cartan_inverse[index, index]
            )
            if new_norm > bound:
                break
            current[index] = coefficient
            extend(position + 1, new_norm)
            coefficient += 1
        current[index] = 0

    extend(0, QQ(0))
    return result


def rational_ldl(value):
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
    assert value >= 0
    numerator = int(value.numerator())
    denominator = int(value.denominator())
    result = math.isqrt(numerator // denominator)
    while result * result * denominator < numerator:
        result += 1
    return result


def shifted_ellipsoid_shell(quadratic, centre, norm):
    dimension = quadratic.nrows()
    if dimension == 0:
        return [tuple()] if norm == 0 else []
    if norm < 0:
        return []
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
        radius = exact_ceil_sqrt(remaining / diagonal[index])
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
    return result


def coordinate_model(gram, fixed, simple):
    constraint = (gram * fixed.transpose()).augment(gram * simple.transpose())
    extras = []
    for coordinate in range(gram.nrows()):
        if constraint.ncols() == gram.nrows():
            break
        unit = matrix(
            ZZ,
            gram.nrows(),
            1,
            [int(index == coordinate) for index in range(gram.nrows())],
        )
        candidate = constraint.augment(unit)
        if candidate.rank() > constraint.rank():
            constraint = candidate
            extras.append(coordinate)
    assert constraint.nrows() == constraint.ncols() == gram.nrows()
    inverse = constraint.inverse()
    coordinate_norm = inverse * gram * inverse.transpose()
    prefix_size = fixed.nrows() + simple.nrows()
    fixed_quadratic = coordinate_norm[prefix_size:, prefix_size:]
    return inverse, coordinate_norm, fixed_quadratic, tuple(extras)


def enumerate_sixths(gram, d5, ambient_roots, norm, candidate_limit):
    unused_roots, components, simple, unused_cartan = residual_simple_data(
        gram, d5, ambient_roots
    )
    residual_rank = simple.nrows()
    assert residual_rank in (18, 19)
    component_choices = [
        dominant_labels_up_to(
            (component * gram * component.transpose()).inverse(), QQ(norm)
        )
        for component in components
    ]
    complement = (d5 * gram).right_kernel_matrix()
    complement_gram = complement * gram * complement.transpose()
    pairing = complement * gram * simple.transpose()
    if residual_rank == 19:
        coordinate_map = pairing
    else:
        for coordinate in range(19):
            extra = matrix(
                ZZ, 19, 1, [int(index == coordinate) for index in range(19)]
            )
            candidate = pairing.augment(extra)
            if candidate.rank() == 19:
                coordinate_map = candidate
                break
        else:
            raise AssertionError("failed to complete sixth Dynkin coordinates")
    coordinate_inverse = coordinate_map.inverse()
    coordinate_norm = coordinate_inverse * complement_gram * coordinate_inverse.transpose()
    candidates = []
    for choice in itertools.product(*component_choices):
        projected_norm = sum(value for unused_labels, value in choice)
        if projected_norm > norm:
            continue
        labels = vector(
            ZZ,
            [entry for component_labels, unused in choice for entry in component_labels],
        )
        if residual_rank == 19:
            fixed_values = [None] if projected_norm == norm else []
        else:
            root_block = coordinate_norm[:18, :18]
            cross = coordinate_norm[:18, 18]
            fixed_norm = coordinate_norm[18, 18]
            linear = (labels * cross)[0]
            constant = labels * root_block * labels
            radius = RR((QQ(norm) - projected_norm) / fixed_norm).sqrt()
            centre = -linear / fixed_norm
            fixed_values = [
                value
                for value in range(
                    int(floor(RR(centre) - radius)) - 2,
                    int(ceil(RR(centre) + radius)) + 3,
                )
                if fixed_norm * value * value + 2 * linear * value + constant == norm
            ]
        for fixed_value in fixed_values:
            extended = labels if fixed_value is None else vector(ZZ, list(labels) + [fixed_value])
            coordinates = extended * coordinate_inverse
            if not all(entry.denominator() == 1 for entry in coordinates):
                continue
            sixth = vector(ZZ, coordinates) * complement
            first_six = d5.stack(matrix(ZZ, [sixth]))
            if not primitive_rows(first_six):
                continue
            assert sixth * gram * sixth == norm
            candidates.append(
                {
                    "vector": sixth,
                    "dominant_labels": tuple(map(int, labels)),
                    "fixed_coordinate": None if fixed_value is None else int(fixed_value),
                }
            )
    candidates.sort(key=lambda row: (row["dominant_labels"], tuple(row["vector"])))
    total = len(candidates)
    if candidate_limit:
        candidates = candidates[:candidate_limit]
    return candidates, total


def cartan_components(cartan, indices=None):
    indices = set(range(cartan.nrows())) if indices is None else set(indices)
    result = []
    while indices:
        component = {min(indices)}
        pending = list(component)
        indices.difference_update(component)
        while pending:
            current = pending.pop()
            adjacent = {
                index for index in indices if cartan[current, index] == -1
            }
            component.update(adjacent)
            indices.difference_update(adjacent)
            pending.extend(sorted(adjacent))
        result.append(tuple(sorted(component)))
    return result


def component_name(cartan):
    rank = cartan.nrows()
    determinant = abs(int(cartan.det()))
    signed_roots_count = int(pari(cartan).qfminim(2)[0])
    if determinant == rank + 1 and signed_roots_count == rank * (rank + 1):
        return f"A{rank}"
    if rank >= 4 and determinant == 4 and signed_roots_count == 2 * rank * (rank - 1):
        return f"D{rank}"
    return {(6, 3, 72): "E6", (7, 2, 126): "E7", (8, 1, 240): "E8"}[
        (rank, determinant, signed_roots_count)
    ]


def prescribed_label_rows(
    cartan,
    norm_bound,
    root_rank_min,
    root_rank_max,
    support_min,
    support_max,
    all_a_only,
    label_limit,
):
    rank = cartan.nrows()
    inverse = cartan.inverse()
    result = []
    for root_rank in range(root_rank_min, root_rank_max + 1):
        positive_count = rank - root_rank
        if positive_count < 0:
            continue
        for positive_indices in itertools.combinations(range(rank), positive_count):
            zero_indices = sorted(set(range(rank)) - set(positive_indices))
            components = cartan_components(cartan, zero_indices)
            if not (support_min <= len(components) <= support_max):
                continue
            names = [
                component_name(cartan.matrix_from_rows_and_columns(component, component))
                for component in components
            ]
            if all_a_only and any(not name.startswith("A") for name in names):
                continue
            for labels, label_norm in dominant_labels_up_to(
                inverse, norm_bound, positive_indices=positive_indices
            ):
                result.append(
                    {
                        "labels": labels,
                        "label_norm": label_norm,
                        "prescribed_root_type": "+".join(sorted(names)) if names else "0",
                        "prescribed_root_rank": root_rank,
                        "prescribed_support_count": len(components),
                    }
                )
    result.sort(
        key=lambda row: (
            row["prescribed_root_rank"],
            row["prescribed_root_type"],
            row["label_norm"],
            row["labels"],
        )
    )
    total = len(result)
    if label_limit:
        result = result[:label_limit]
    return result, total


def frame_root_record(frame):
    roots, unused_basis, data = roots_and_data(frame)
    if not roots:
        return "0", 0, 0, 1, []
    simple, unused_positive, cartan_rows = deterministic_simple_roots(frame)
    cartan = matrix(ZZ, cartan_rows)
    components = []
    for indices in cartan_components(cartan):
        block = cartan.matrix_from_rows_and_columns(indices, indices)
        name = component_name(block)
        components.append(
            {
                "type": name,
                "rank": block.nrows(),
                "determinant": abs(int(block.det())),
                "signed_root_count": int(pari(block).qfminim(2)[0]),
            }
        )
    multiplicities = Counter(row["type"] for row in components)
    label = "+".join(
        f"{count if count > 1 else ''}{name}"
        for name, count in sorted(multiplicities.items())
    )
    root_determinant = math.prod(row["determinant"] for row in components)
    return label, int(data[0]), int(data[1]), int(root_determinant), components


def enumerate_sevenths(
    gram,
    first_six,
    ambient_roots,
    standard,
    arguments,
):
    unused_roots, unused_components, simple, cartan = residual_simple_data(
        gram, first_six, ambient_roots
    )
    residual_rank = simple.nrows()
    if residual_rank < arguments.source_root_rank_min:
        return [], {"residual_rank": residual_rank, "label_rows": 0, "ellipsoid_solutions": 0}
    pairings = vector(ZZ, standard["seventh_pairings"])
    norm = ZZ(standard["seventh_norm"])
    first_six_gram = first_six * gram * first_six.transpose()
    base_norm = pairings * first_six_gram.inverse() * pairings
    label_budget = QQ(norm) - base_norm
    if label_budget < 0:
        return [], {"residual_rank": residual_rank, "label_rows": 0, "ellipsoid_solutions": 0}
    label_rows, total_label_rows = prescribed_label_rows(
        cartan,
        label_budget,
        arguments.source_root_rank_min,
        arguments.source_root_rank_max,
        arguments.source_support_min,
        arguments.source_support_max,
        arguments.all_a_only,
        arguments.seventh_label_limit,
    )
    inverse, coordinate_norm, fixed_quadratic, extra_coordinates = coordinate_model(
        gram, first_six, simple
    )
    fixed_dimension = 18 - residual_rank
    assert fixed_quadratic.nrows() == fixed_dimension
    if fixed_dimension > arguments.fixed_dimension_max:
        return [], {
            "residual_rank": residual_rank,
            "fixed_dimension": fixed_dimension,
            "label_rows": total_label_rows,
            "skipped_fixed_dimension": True,
            "ellipsoid_solutions": 0,
        }
    fixed_inverse = fixed_quadratic.inverse() if fixed_dimension else matrix(QQ, 0, 0)
    prefix_size = 6 + residual_rank
    candidates = []
    ellipsoid_count = 0
    for label_row in label_rows:
        labels = vector(ZZ, label_row["labels"])
        prefix = vector(QQ, list(pairings) + list(labels))
        prefix_norm = prefix * coordinate_norm[:prefix_size, :prefix_size] * prefix
        if fixed_dimension:
            linear = prefix * coordinate_norm[:prefix_size, prefix_size:]
            centre = -linear * fixed_inverse
            minimum_norm = prefix_norm - linear * fixed_inverse * linear
        else:
            centre = vector(QQ, [])
            minimum_norm = prefix_norm
        assert minimum_norm == base_norm + label_row["label_norm"]
        fixed_solutions = shifted_ellipsoid_shell(
            fixed_quadratic, centre, QQ(norm) - minimum_norm
        )
        ellipsoid_count += len(fixed_solutions)
        for fixed_coordinates in fixed_solutions:
            coordinates = vector(
                QQ, list(pairings) + list(labels) + list(fixed_coordinates)
            )
            seventh = coordinates * inverse
            if not all(entry.denominator() == 1 for entry in seventh):
                continue
            seventh = vector(ZZ, seventh)
            auxiliary_basis = first_six.stack(matrix(ZZ, [seventh]))
            if not primitive_rows(auxiliary_basis):
                continue
            if auxiliary_basis * gram * auxiliary_basis.transpose() != standard["gram"]:
                continue
            candidates.append(
                {
                    "vector": seventh,
                    "labels": label_row["labels"],
                    "fixed_coordinates": tuple(map(int, fixed_coordinates)),
                    "fixed_dimension": fixed_dimension,
                    "extra_coordinates": extra_coordinates,
                    "prescribed_root_type": label_row["prescribed_root_type"],
                    "prescribed_root_rank": label_row["prescribed_root_rank"],
                    "prescribed_support_count": label_row["prescribed_support_count"],
                }
            )
    return candidates, {
        "residual_rank": residual_rank,
        "fixed_dimension": fixed_dimension,
        "label_rows": total_label_rows,
        "retained_label_rows": len(label_rows),
        "ellipsoid_solutions": ellipsoid_count,
        "integral_primitive_sevenths": len(candidates),
    }


def source_payload(frame, minimized):
    label, root_rank, root_count, root_det, components = frame_root_record(frame)
    mw_height = minimized["mw_height"]
    return {
        "gram": rows(frame),
        "gram_sha256": gram_digest(frame),
        "root_type": label,
        "root_components": components,
        "root_rank": root_rank,
        "signed_root_count": root_count,
        "root_determinant": root_det,
        "mw_rank_for_rho_19": 17 - root_rank,
        "root_lattice_primitive": minimized["root_lattice_primitive"],
        "root_smith_invariants": list(map(int, minimized["root_smith_invariants"])),
        "root_adapted_gram": rows(minimized["frame"]),
        "root_adapted_basis_rows_in_source_basis": rows(minimized["basis"]),
        "mw_height_gram": None if mw_height is None else rational_rows(mw_height),
        "mw_regulator": None if mw_height is None else str(abs(mw_height.det())),
        "torsion": 1 if minimized["root_lattice_primitive"] else "REQUIRES_GLUE_ANALYSIS",
        "determinant": int(frame.det()),
    }


def build(arguments):
    database_path = arguments.database.resolve()
    catalog_path = arguments.catalog.resolve()
    anchors_path = arguments.anchors.resolve()
    database = json.loads(database_path.read_text())
    catalog = json.loads(catalog_path.read_text())
    anchors_payload = json.loads(anchors_path.read_text())
    h3_frame, unused_h3_change = reduced_gram(load_matrix(H3_FRAME))
    h3_height = minimize_child_frame(h3_frame)["mw_height"]
    gram_by_label = {
        row["label"]: matrix(ZZ, row["gram"])
        for row in catalog["rooted_niemeier_lattices"]
    }
    anchors = list(anchors_payload["anchors"])

    requested_ns = set(arguments.ns_id)
    if arguments.all_ns:
        requested_ns = set()
    elif not requested_ns:
        requested_ns = {"NS0001"}
    requested_ambients = set(arguments.ambient_label)
    if arguments.all_ambients:
        requested_ambients = set()
    elif not requested_ambients:
        requested_ambients = {"D16_E8"}
    anchors = [
        row for row in anchors
        if not requested_ambients or row["niemeier"] in requested_ambients
    ]
    if arguments.anchor_limit:
        anchors = anchors[:arguments.anchor_limit]

    selected_ns = []
    skipped_ns = []
    for ns in database["ns_classes"]:
        if requested_ns and ns["ns_id"] not in requested_ns:
            continue
        targets = [
            frame for frame in ns["frames"]
            if arguments.target_mw_min
            <= int(frame["mw_rank_for_rho_19"])
            <= arguments.target_mw_max
        ]
        if not targets:
            skipped_ns.append({"ns_id": ns["ns_id"], "reason": "NO_TARGET_IN_MW_BAND"})
            continue
        selected_ns.append((ns, targets))
    if arguments.ns_limit:
        selected_ns = selected_ns[:arguments.ns_limit]
    known_ns = {row["ns_id"] for row in database["ns_classes"]}
    if requested_ns - known_ns:
        raise ValueError(f"unknown NS ids: {sorted(requested_ns - known_ns)}")
    if not selected_ns:
        raise ValueError("no eligible NS classes selected")
    if not anchors:
        raise ValueError("no D5 anchors selected")

    source_classes = []
    accounting = Counter()
    per_ns = []
    truncated = bool(
        arguments.ns_limit
        or arguments.anchor_limit
        or arguments.sixth_candidate_limit
        or arguments.seventh_label_limit
    )
    for ns, targets in selected_ns:
        auxiliary = matrix(ZZ, ns["auxiliary_gram"])
        standardized = standardize_auxiliary(auxiliary, arguments.sixth_norm_max)
        if standardized is None:
            skipped_ns.append(
                {"ns_id": ns["ns_id"], "reason": "NO_PRIMITIVE_D5_ORTHOGONAL_SIXTH_WITHIN_BOUND"}
            )
            continue
        target_genus = Genus(matrix(ZZ, targets[0]["gram"]))
        ns_classes = []
        ns_classes_by_reduced_digest = {}
        ns_accounting = Counter()
        anchor_rows = []
        for anchor in anchors:
            ambient_label = anchor["niemeier"]
            ambient = gram_by_label[ambient_label]
            ambient_roots = signed_roots(ambient)
            d5 = matrix(ZZ, anchor["D5_basis_in_ambient"])
            assert d5 * ambient * d5.transpose() == D5_GRAM
            sixths, total_sixths = enumerate_sixths(
                ambient,
                d5,
                ambient_roots,
                standardized["sixth_norm"],
                arguments.sixth_candidate_limit,
            )
            ns_accounting["sixth_candidates_total"] += total_sixths
            ns_accounting["sixth_candidates_visited"] += len(sixths)
            anchor_accounting = Counter()
            for sixth_index, sixth_row in enumerate(sixths, 1):
                first_six = d5.stack(matrix(ZZ, [sixth_row["vector"]]))
                sevenths, seventh_accounting = enumerate_sevenths(
                    ambient,
                    first_six,
                    ambient_roots,
                    standardized,
                    arguments,
                )
                anchor_accounting.update(
                    {
                        "seventh_label_rows": seventh_accounting.get("label_rows", 0),
                        "seventh_ellipsoid_solutions": seventh_accounting.get("ellipsoid_solutions", 0),
                        "integral_primitive_sevenths": len(sevenths),
                    }
                )
                for seventh_row in sevenths:
                    auxiliary_basis = first_six.stack(matrix(ZZ, [seventh_row["vector"]]))
                    complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
                    if complement_basis.nrows() != 17 or not primitive_rows(complement_basis):
                        raise AssertionError("Niemeier orthogonal complement is not saturated")
                    raw_frame = complement_basis * ambient * complement_basis.transpose()
                    if raw_frame.det() != auxiliary.det() or Genus(raw_frame) != target_genus:
                        raise AssertionError("source complement left the target genus")
                    frame, frame_change = reduced_gram(raw_frame)
                    label, root_rank, unused_count, unused_det, unused_components = frame_root_record(frame)
                    if not (
                        arguments.source_root_rank_min
                        <= root_rank
                        <= arguments.source_root_rank_max
                    ):
                        raise AssertionError("prescribed zero labels gave the wrong root rank")
                    if root_rank != seventh_row["prescribed_root_rank"]:
                        raise AssertionError("prescribed and actual root ranks disagree")
                    minimized = minimize_child_frame(frame)
                    reduced_digest = gram_digest(frame)
                    match = ns_classes_by_reduced_digest.get(reduced_digest)
                    # Sequential Weyl faces often reach the same frame
                    # repeatedly. General rank-17 qfisom calls can require
                    # gigabytes even for these root-rich forms, so this
                    # discovery ledger merges only exact equality after a
                    # deterministic unimodular LLL reduction.
                    provenance = {
                        "niemeier": ambient_label,
                        "anchor_index": int(anchor["anchor_index"]),
                        "sixth_index": sixth_index,
                        "sixth_vector": list(map(int, sixth_row["vector"])),
                        "sixth_dominant_labels": list(sixth_row["dominant_labels"]),
                        "seventh_vector": list(map(int, seventh_row["vector"])),
                        "seventh_dominant_labels": list(seventh_row["labels"]),
                        "prescribed_root_type": seventh_row["prescribed_root_type"],
                        "prescribed_support_count": seventh_row["prescribed_support_count"],
                        "auxiliary_basis_in_ambient": rows(auxiliary_basis),
                        "complement_basis_in_ambient": rows(complement_basis),
                        "source_basis_rows_in_complement_basis": rows(frame_change),
                    }
                    if match is None:
                        match = {
                            "frame": frame,
                            "source": source_payload(frame, minimized),
                            "provenance": [],
                        }
                        ns_classes.append(match)
                    ns_classes_by_reduced_digest[reduced_digest] = match
                    match["provenance"].append(provenance)
            anchor_rows.append(
                {
                    "niemeier": ambient_label,
                    "anchor_index": int(anchor["anchor_index"]),
                    "sixth_candidates_total": total_sixths,
                    "sixth_candidates_visited": len(sixths),
                    **{key: int(value) for key, value in sorted(anchor_accounting.items())},
                }
            )
        ns_classes.sort(
            key=lambda row: (
                row["source"]["mw_rank_for_rho_19"],
                row["source"]["root_type"],
                row["source"]["gram_sha256"],
            )
        )
        target_rows = sorted(
            [
                {
                    "frame_id": row["frame_id"],
                    "root_type": row["root_type"],
                    "root_rank": int(row["root_rank"]),
                    "mw_rank_for_rho_19": int(row["mw_rank_for_rho_19"]),
                    "gram_sha256": row["gram_sha256"],
                }
                for row in targets
            ],
            key=lambda row: row["frame_id"],
        )
        for index, row in enumerate(ns_classes, 1):
            source = row["source"]
            candidate_height = minimize_child_frame(row["frame"])["mw_height"]
            matches_h3_height_control = bool(
                int(ns["determinant"]) == 948
                and source["root_type"] == "E7+E8"
                and candidate_height is not None
                and rational_form_isometric(candidate_height, h3_height)
            )
            source_record = {
                "source_id": f"{ns['ns_id']}-S{index:03d}",
                "ns_id": ns["ns_id"],
                "determinant": int(ns["determinant"]),
                "source": source,
                "same_ns_high_rank_targets": target_rows,
                "embedding_count_in_enumerated_cover": len(row["provenance"]),
                "niemeier_provenance": row["provenance"],
                "matches_h3_e7e8_mw2_height_control": matches_h3_height_control,
            }
            source_classes.append(source_record)
        accounting["ns_classes_searched"] += 1
        accounting["source_reduced_gram_classes"] += len(ns_classes)
        accounting["mw_at_most_two_sources"] += len(ns_classes)
        accounting["h3_positive_control_hits"] += sum(
            row["source"]["root_type"] == "E7+E8"
            and row["source"]["mw_height_gram"] is not None
            and rational_form_isometric(
                matrix(QQ, row["source"]["mw_height_gram"]), h3_height
            )
            for row in ns_classes
        )
        per_ns.append(
            {
                "ns_id": ns["ns_id"],
                "determinant": int(ns["determinant"]),
                "target_count": len(target_rows),
                "target_frame_ids": [row["frame_id"] for row in target_rows],
                "standardized_auxiliary_gram": rows(standardized["gram"]),
                "standardizing_basis_rows": rows(standardized["basis"]),
                "sixth_norm": standardized["sixth_norm"],
                "seventh_norm": standardized["seventh_norm"],
                "source_reduced_gram_classes": len(ns_classes),
                "anchor_accounting": anchor_rows,
            }
        )
        print(
            "FOUNDRYPRESCRIBED|"
            f"ns={ns['ns_id']}|targets={len(target_rows)}|"
            f"sources={len(ns_classes)}|status=PASS_NS",
            flush=True,
        )

    source_classes.sort(key=lambda row: (row["ns_id"], row["source_id"]))
    if arguments.require_hit and not source_classes:
        raise RuntimeError("prescribed-root enumeration produced no source")
    if arguments.require_h3_control and not any(
        row["matches_h3_e7e8_mw2_height_control"] for row in source_classes
    ):
        raise RuntimeError("H3 E7+E8/MW2 positive control was not recovered")
    completeness = (
        "DETERMINISTIC_PREFIX_SEARCH_EXACT_RETAINED_SOURCES"
        if truncated
        else "EXACT_WITHIN_DECLARED_NS_AMBIENT_NORM_AND_ROOT_FILTERS"
    )
    result = {
        "schema": "elkies-k3.lattice-foundry-prescribed-root-sources.v1",
        "status": f"PASS_{completeness}",
        "objective": (
            "Directly enumerate MW0--2 Niemeier source complements for exact "
            "foundry NS classes already carrying catalogued MW15--17 targets."
        ),
        "proof_boundary": {
            "proved": (
                "Every emitted row is a primitive embedding of the displayed exact "
                "rank-seven auxiliary in a hash-pinned Niemeier lattice. Its rank-17 "
                "orthogonal complement is saturated, has the displayed exact ADE/MW "
                "data, and is in the same positive-definite genus as every attached "
                "same-NS high-rank target. The H3 regression compares the exact "
                "E7+E8 root type and integral-isometry class of the binary MW height "
                "form; it is not promoted to a full-frame isometry certificate."
            ),
            "search_completeness": (
                "Enumeration is complete only for the selected NS classes, D5 anchor "
                "orbits, Niemeier ambients, sixth norm bound, prescribed source-root "
                "and support filters, fixed-space dimension bound, and any declared "
                "candidate prefixes."
            ),
            "not_proved": (
                "Distinct deterministic reduced Grams are not claimed to be distinct "
                "integral-isometry or J2 classes. A same-auxiliary source/target pair "
                "is not a marked elliptic-neighbour "
                "corridor. No rational source marking, Galois descent, nef U embedding, "
                "equation, or specialization-rank behavior is inferred. A miss outside "
                "the declared finite window is not a non-existence theorem."
            ),
        },
        "search": {
            "ns_ids": [row[0]["ns_id"] for row in selected_ns],
            "ambient_labels": sorted({row["niemeier"] for row in anchors}),
            "target_mw_range": [arguments.target_mw_min, arguments.target_mw_max],
            "source_root_rank_range": [
                arguments.source_root_rank_min,
                arguments.source_root_rank_max,
            ],
            "source_mw_range_for_rho_19": [
                17 - arguments.source_root_rank_max,
                17 - arguments.source_root_rank_min,
            ],
            "source_support_range": [
                arguments.source_support_min,
                arguments.source_support_max,
            ],
            "all_a_only": arguments.all_a_only,
            "sixth_norm_max": arguments.sixth_norm_max,
            "fixed_dimension_max": arguments.fixed_dimension_max,
            "ns_limit": arguments.ns_limit,
            "anchor_limit": arguments.anchor_limit,
            "sixth_candidate_limit": arguments.sixth_candidate_limit,
            "seventh_label_limit": arguments.seventh_label_limit,
            "truncated_by_prefix_limit": truncated,
        },
        "accounting": {
            **{key: int(value) for key, value in sorted(accounting.items())},
            "selected_ns_classes": len(selected_ns),
            "selected_d5_anchors": len(anchors),
            "skipped_ns_classes": len(skipped_ns),
        },
        "skipped_ns_classes": skipped_ns,
        "ns_accounting": per_ns,
        "sources": source_classes,
        "inputs": {
            relative(path): digest(path)
            for path in (database_path, catalog_path, anchors_path, H3_FRAME)
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/enumerate_lattice_foundry_prescribed_root_sources.sage"
        ),
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--database", type=Path, default=DATABASE)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--anchors", type=Path, default=ANCHORS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ns-id", action="append", default=[])
    parser.add_argument("--all-ns", action="store_true")
    parser.add_argument("--ambient-label", action="append", default=[])
    parser.add_argument("--all-ambients", action="store_true")
    parser.add_argument("--target-mw-min", type=int, default=15)
    parser.add_argument("--target-mw-max", type=int, default=17)
    parser.add_argument("--source-root-rank-min", type=int, default=15)
    parser.add_argument("--source-root-rank-max", type=int, default=17)
    parser.add_argument("--source-support-min", type=int, default=1)
    parser.add_argument("--source-support-max", type=int, default=17)
    parser.add_argument("--all-a-only", action="store_true")
    parser.add_argument("--sixth-norm-max", type=int, default=24)
    parser.add_argument("--fixed-dimension-max", type=int, default=6)
    parser.add_argument("--ns-limit", type=int, default=0)
    parser.add_argument("--anchor-limit", type=int, default=0)
    parser.add_argument("--sixth-candidate-limit", type=int, default=0)
    parser.add_argument("--seventh-label-limit", type=int, default=0)
    parser.add_argument("--require-hit", action="store_true")
    parser.add_argument("--require-h3-control", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.all_ns and arguments.ns_id:
        parser.error("use either --all-ns or --ns-id, not both")
    if arguments.all_ambients and arguments.ambient_label:
        parser.error("use either --all-ambients or --ambient-label, not both")
    if not (0 <= arguments.target_mw_min <= arguments.target_mw_max <= 17):
        parser.error("invalid target MW range")
    if not (
        0
        <= arguments.source_root_rank_min
        <= arguments.source_root_rank_max
        <= 17
    ):
        parser.error("invalid source root-rank range")
    if arguments.sixth_norm_max <= 0 or arguments.sixth_norm_max % 2:
        parser.error("--sixth-norm-max must be a positive even integer")
    result = build(arguments)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.is_file() or output.read_text() != serialized:
            raise SystemExit("prescribed-root source ledger is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(serialized)
    print(
        "FOUNDRYPRESCRIBED|"
        f"ns={result['accounting']['ns_classes_searched']}|"
        f"sources={result['accounting']['source_reduced_gram_classes']}|"
        f"h3={result['accounting']['h3_positive_control_hits']}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
