#!/usr/bin/env sage
"""Complete rootless-J2 frame census for the determinant-1092 X1092 lattice.

Every embedding of the certified rank-seven Nishiyama auxiliary sends its D5
root sublattice to one of the already complete D5 anchor orbits.  The sixth
and seventh auxiliary vectors are enumerated in exact residual Weyl chambers;
strictly positive final Dynkin labels are equivalent to a rootless complement.
"""

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from sage.all import QQ, ZZ, ceil, diagonal_matrix, floor, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
ANCHORS = ROOT / "artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json"
AUXILIARY = ROOT / "artifacts/generated-results/elliptic-curves/det1092_nishiyama_auxiliary_v1.json"
PARENT = ROOT / "artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/det1092_rootless_j2_niemeier_census_v1.json"

D5 = matrix(ZZ, [
    [2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, -1, -1],
    [0, 0, -1, 2, 0], [0, 0, -1, 0, 2],
])
SIXTH_PAIRINGS = vector(ZZ, [1, 1, 0, 0, -1])
SIXTH_NORM = ZZ(20)
SEVENTH_PAIRINGS = vector(ZZ, [-1, 0, 1, 1, 0, -2])
SEVENTH_NORM = ZZ(22)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def gram_hash(value):
    return hashlib.sha256(("\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n").encode()).hexdigest()


def signed_roots(gram):
    data = pari(gram).qfminim(2)
    reps = matrix(ZZ, data[2].sage()).transpose()
    result = [vector(ZZ, row) for row in reps.rows()]
    return result + [-root for root in result]


def root_components(gram, roots):
    unseen, result = set(range(len(roots))), []
    while unseen:
        component, frontier = {min(unseen)}, [min(unseen)]
        unseen.difference_update(component)
        while frontier:
            current = frontier.pop()
            neighbours = {index for index in unseen if roots[current] * gram * roots[index] != 0}
            component.update(neighbours)
            unseen.difference_update(neighbours)
            frontier.extend(sorted(neighbours))
        result.append([roots[index] for index in sorted(component)])
    return result


def simple_roots(gram, roots):
    if not roots:
        return matrix(ZZ, 0, gram.nrows())
    for trial in range(1, 100):
        chamber = vector(ZZ, [(index + 1) ** 2 + trial * (index + 1) + trial**2 for index in range(gram.nrows())])
        values = [root * gram * chamber for root in roots]
        if all(values):
            break
    else:
        raise ArithmeticError("could not choose a regular root chamber")
    positive = [root for root, value in zip(roots, values) if value > 0]
    positive_set = {tuple(map(int, root)) for root in positive}
    result = matrix(ZZ, [root for root in positive if not any(tuple(map(int, root - other)) in positive_set for other in positive)])
    assert result.rank() == result.nrows()
    return result


def labels_up_to(cartan_inverse, bound, strict):
    rank, current, result = cartan_inverse.nrows(), [], []
    lower = 1 if strict else 0

    def visit(index, norm):
        if index == rank:
            result.append((tuple(current), norm))
            return
        cross = sum(current[previous] * cartan_inverse[previous, index] for previous in range(index))
        coefficient = lower
        while True:
            candidate = norm + 2 * coefficient * cross + coefficient**2 * cartan_inverse[index, index]
            if candidate > bound:
                break
            current.append(coefficient)
            visit(index + 1, candidate)
            current.pop()
            coefficient += 1

    visit(0, QQ(0))
    return result


def ldl(value):
    size, lower, diagonal = value.nrows(), matrix(QQ, value.nrows()), []
    for row in range(size):
        lower[row, row] = 1
        diagonal.append(value[row, row] - sum(lower[row, index] ** 2 * diagonal[index] for index in range(row)))
        assert diagonal[-1] > 0
        for later in range(row + 1, size):
            lower[later, row] = (value[later, row] - sum(lower[later, index] * diagonal[index] * lower[row, index] for index in range(row))) / diagonal[index]
    assert lower * diagonal_matrix(diagonal) * lower.transpose() == value
    return lower, diagonal


def ceil_sqrt(value):
    numerator, denominator = int(value.numerator()), int(value.denominator())
    answer = math.isqrt(numerator // denominator)
    while answer * answer * denominator < numerator:
        answer += 1
    return answer


def ellipsoid_shell(quadratic, centre, norm):
    dimension = quadratic.nrows()
    if not dimension:
        return [tuple()] if norm == 0 else []
    lower, diagonal = ldl(quadratic)
    coordinates, shifted, result = [ZZ(0)] * dimension, [QQ(0)] * dimension, []

    def visit(index, remaining):
        if index < 0:
            if remaining == 0:
                result.append(tuple(map(int, coordinates)))
            return
        tail = sum(lower[later, index] * shifted[later] for later in range(index + 1, dimension))
        local_centre = centre[index] - tail
        radius = ceil_sqrt(remaining / diagonal[index])
        for coordinate in range(int(floor(local_centre)) - radius - 1, int(ceil(local_centre)) + radius + 2):
            local = QQ(coordinate) - local_centre
            term = diagonal[index] * local**2
            if term <= remaining:
                coordinates[index], shifted[index] = coordinate, QQ(coordinate) - centre[index]
                visit(index - 1, remaining - term)

    visit(dimension - 1, norm)
    for candidate in result:
        difference = vector(QQ, candidate) - centre
        assert difference * quadratic * difference == norm
    return result


def primitive(value):
    diagonal = value.smith_form()[0]
    return all(abs(diagonal[index, index]) == 1 for index in range(value.nrows()))


def chamber_vectors(gram, base, pairings, norm, strict):
    """All primitive next vectors in one residual Weyl chamber, exactly."""
    roots = [root for root in signed_roots(gram) if root * gram * base.transpose() == 0]
    components = root_components(gram, roots)
    simple = [simple_roots(gram, component) for component in components]
    root_basis = matrix(ZZ, [row for component in simple for row in component.rows()]) if simple else matrix(ZZ, 0, gram.nrows())
    constraints = (gram * base.transpose()).augment(gram * root_basis.transpose())
    extras = []
    for coordinate in range(gram.nrows()):
        if constraints.ncols() == gram.nrows():
            break
        unit = matrix(ZZ, gram.nrows(), 1, [int(index == coordinate) for index in range(gram.nrows())])
        if constraints.augment(unit).rank() > constraints.rank():
            constraints = constraints.augment(unit)
            extras.append(coordinate)
    assert constraints.nrows() == constraints.ncols() == gram.nrows()
    inverse = constraints.inverse()
    coordinate_norm = inverse * gram * inverse.transpose()
    prefix_size = base.nrows() + root_basis.nrows()
    fixed_dimension = len(extras)
    fixed_quadratic = coordinate_norm[prefix_size:, prefix_size:]
    fixed_inverse = fixed_quadratic.inverse() if fixed_dimension else matrix(QQ, 0, 0)
    component_labels = [labels_up_to((component * gram * component.transpose()).inverse(), QQ(norm), strict) for component in simple]
    candidates, accounting = [], Counter()
    for choice in itertools.product(*component_labels) if component_labels else [tuple()]:
        labels = vector(ZZ, [entry for component, unused in choice for entry in component])
        prefix = vector(QQ, list(pairings) + list(labels))
        prefix_norm = prefix * coordinate_norm[:prefix_size, :prefix_size] * prefix
        if fixed_dimension:
            linear = prefix * coordinate_norm[:prefix_size, prefix_size:]
            centre = -linear * fixed_inverse
            minimum = prefix_norm - linear * fixed_inverse * linear
        else:
            centre, minimum = vector(QQ, []), prefix_norm
        if minimum > norm:
            continue
        accounting["label_choices"] += 1
        for fixed in ellipsoid_shell(fixed_quadratic, centre, QQ(norm) - minimum):
            coordinates = vector(QQ, list(prefix) + list(fixed))
            candidate = coordinates * inverse
            if not all(value.denominator() == 1 for value in candidate):
                accounting["nonintegral"] += 1
                continue
            candidate = vector(ZZ, candidate)
            assert candidate * gram * candidate == norm
            assert candidate * gram * base.transpose() == pairings
            enlarged = base.stack(matrix(ZZ, [candidate]))
            if not primitive(enlarged):
                accounting["nonprimitive"] += 1
                continue
            candidates.append({"vector": candidate, "labels": list(map(int, labels)), "fixed_dimension": fixed_dimension})
    return candidates, {"residual_root_rank": root_basis.nrows(), "fixed_dimension": fixed_dimension, **{key: int(value) for key, value in accounting.items()}}


def build():
    auxiliary = json.loads(AUXILIARY.read_text())
    parent = json.loads(PARENT.read_text())
    catalog, anchors = json.loads(CATALOG.read_text()), json.loads(ANCHORS.read_text())
    assert auxiliary["status"] == "PASS_EXACT_AUXILIARY_AND_UNIMODULAR_GLUE"
    auxiliary_gram = matrix(ZZ, auxiliary["auxiliary"]["gram"])
    assert auxiliary_gram[:5, :5] == D5
    assert vector(ZZ, auxiliary_gram.row(5)[:5]) == SIXTH_PAIRINGS
    assert vector(ZZ, auxiliary_gram.row(6)[:6]) == SEVENTH_PAIRINGS
    frame = matrix(QQ, parent["generic_height_gram"])
    assert frame.det() == 1092 and int(pari(frame).qfminim(2)[0]) == 0
    grams = {row["label"]: matrix(ZZ, row["gram"]) for row in catalog["rooted_niemeier_lattices"]}
    rootless, anchor_accounting = [], []
    for anchor in anchors["anchors"]:
        gram, d5 = grams[anchor["niemeier"]], matrix(ZZ, anchor["D5_basis_in_ambient"])
        assert d5 * gram * d5.transpose() == D5
        sixths, sixth_data = chamber_vectors(gram, d5, SIXTH_PAIRINGS, SIXTH_NORM, False)
        seventh_count = 0
        for sixth in sixths:
            first_six = d5.stack(matrix(ZZ, [sixth["vector"]]))
            sevenths, seventh_data = chamber_vectors(gram, first_six, SEVENTH_PAIRINGS, SEVENTH_NORM, True)
            for seventh in sevenths:
                seventh_count += 1
                basis = first_six.stack(matrix(ZZ, [seventh["vector"]]))
                assert basis * gram * basis.transpose() == matrix(ZZ, auxiliary["auxiliary"]["gram"])
                complement = (basis * gram).right_kernel_matrix()
                complement_gram = complement * gram * complement.transpose()
                assert complement_gram.nrows() == 17 and complement_gram.det() == 1092
                assert int(pari(complement_gram).qfminim(2)[0]) == 0
                assert pari(frame).qfisom(pari(complement_gram)) != 0
                rootless.append({"anchor": f'{anchor["niemeier"]}:{anchor["anchor_index"]}', "niemeier": anchor["niemeier"], "basis": basis, "complement": complement_gram, "sixth": sixth, "seventh": seventh, "seventh_data": seventh_data})
        anchor_accounting.append({"anchor": f'{anchor["niemeier"]}:{anchor["anchor_index"]}', "niemeier": anchor["niemeier"], "sixth_candidates": len(sixths), "seventh_rootless_embeddings": seventh_count, "sixth_data": sixth_data})
        print(f'DET1092J2|anchor={anchor_accounting[-1]["anchor"]}|sixth={len(sixths)}|rootless={seventh_count}', flush=True)
    classes = []
    for embedding in rootless:
        hit = next((item for item in classes if pari(item["gram"]).qfisom(pari(embedding["complement"])) != 0), None)
        if hit is None:
            hit = {"gram": embedding["complement"], "representative": embedding, "count": 0, "anchors": Counter()}
            classes.append(hit)
        hit["count"] += 1
        hit["anchors"][embedding["anchor"]] += 1
    class_rows = []
    for index, item in enumerate(classes, 1):
        gram, representative = item["gram"], item["representative"]
        class_rows.append({"class_index": index, "gram": rows(gram), "gram_sha256": gram_hash(gram), "determinant": int(gram.det()), "minimum": 4, "norm4_unoriented_pairs": int(pari(gram).qfminim(4)[0]) // 2, "automorphism_group_order": int(pari(gram).qfauto()[0]), "matches_recovered_curve302_frame": bool(pari(frame).qfisom(pari(gram)) != 0), "embedding_count_in_cover": item["count"], "anchor_counts": dict(sorted(item["anchors"].items())), "representative_embedding": {"anchor": representative["anchor"], "niemeier": representative["niemeier"], "auxiliary_basis_in_ambient": rows(representative["basis"]), "sixth_labels": representative["sixth"]["labels"], "seventh_labels": representative["seventh"]["labels"]}})
    assert sum(row["matches_recovered_curve302_frame"] for row in class_rows) == 1
    return {"schema": "elkies-k3.det1092-rootless-j2-niemeier-census.v1", "status": "PASS_COMPLETE_ROOTLESS_J2_CLASSIFICATION", "inputs": {str(path.relative_to(ROOT)): digest(path) for path in (CATALOG, ANCHORS, AUXILIARY, PARENT)}, "classification_scope": {"proved": "Every embedding of the certified D5-containing auxiliary is covered by a complete D5 anchor orbit, a residual-Weyl dominant sixth vector, and a strictly dominant seventh vector. Exact integral ellipsoid enumeration and primitive checks retain precisely rootless rank17 complements, which are deduplicated by integral isometry.", "boundary": "This is an O(NS)/J2 frame classification. It is not a J1 surface-automorphism classification, a rational marked-U realization, or an equation catalogue."}, "accounting": {"rooted_niemeier_classes": 23, "D5_anchor_orbits": len(anchors["anchors"]), "primitive_rootless_embeddings_in_cover": len(rootless), "rootless_complement_isometry_classes": len(class_rows), "anchors": anchor_accounting}, "rootless_classes": class_rows}


def canonical(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    payload = canonical(build())
    if arguments.build:
        assert not arguments.output.exists()
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(payload)
        print(f"WROTE {arguments.output.relative_to(ROOT)}")
    assert arguments.output.read_bytes() == payload
    print("PASS complete determinant1092 rootless J2 Niemeier classification")
