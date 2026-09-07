#!/usr/bin/env sage
"""Classify the two CM extensions of the H2 Kumar frame.

At a CM specialization a primitive negative vector v in the generic
transcendental lattice becomes algebraic.  With the elliptic U fixed, the
positive-definite fibration frame is the primitive closure of

    K_H2 + <abs(v^2)>.

The closure indices are the divisibilities of v: 316 for Delta=-3 and 79
for Delta=-24.  This script constructs the corresponding discriminant-form
glues and reports the enhanced root systems.
"""

from pathlib import Path
from itertools import combinations
import argparse

from sage.all import CartanMatrix, Matrix, QQ, QuadraticForm, ZZ, block_diagonal_matrix, gcd, matrix, pari, prod, vector, xgcd


BASE = Path(__file__).resolve().parents[1]

parser = argparse.ArgumentParser()
parser.add_argument(
    "--search-max-q",
    type=int,
    default=60,
    help="search Weyl-orbit representatives using the level-79 section up to this q",
)
parser.add_argument(
    "--orbit-counts-only",
    action="store_true",
    help="enumerate exact Weyl-orbit counts without constructing their neighbor frames",
)
parser.add_argument(
    "--print-q80-frame",
    action="store_true",
    help="print the selected q=80 E6+D5+A3 generic and CM frame Grams",
)
parser.add_argument(
    "--print-cm43-frame",
    action="store_true",
    help="print the discriminant-43 E7+E8/MW3 closure frame Gram",
)
parser.add_argument(
    "--print-marked-cm43-frame",
    action="store_true",
    help="print the Q79-marked glue-211 CM43 frame and q60 fiber",
)
parser.add_argument(
    "--print-q60-in-q80",
    action="store_true",
    help="print every q=60 presentation fiber in the selected q=80 NS basis",
)
parser.add_argument(
    "--print-markings-in-q80",
    action="store_true",
    help="print the old height-4 and level-79 frame directions in q80 NS coordinates",
)
args = parser.parse_args()


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


def matrix_rows(value):
    return tuple(tuple(row) for row in value.rows())


def ade_minimal_correction(invariant, dual_coordinates, component_gram):
    """Identify the minimal ADE discriminant-class norm exactly."""
    rank, root_count, determinant = invariant
    dual_coordinates = vector(QQ, dual_coordinates)
    order = ZZ(1)
    for value in dual_coordinates:
        order = order.lcm(ZZ(value.denominator()))
    raw_norm = dual_coordinates * component_gram * dual_coordinates

    def mod_two(value):
        value = QQ(value)
        return value - 2 * (value / 2).floor()

    candidates = []
    if root_count == rank * (rank + 1) and determinant == rank + 1:
        # A_rank; k and rank+1-k have the same correction.
        modulus = rank + 1
        for k in range(modulus):
            class_order = ZZ(1) if k == 0 else ZZ(modulus // gcd(k, modulus))
            correction = QQ(k * (modulus - k)) / modulus
            candidates.append((class_order, correction))
    elif root_count == 2 * rank * (rank - 1) and determinant == 4:
        # D_rank: zero, vector, and the two spinor classes.
        candidates.extend(((ZZ(1), QQ(0)), (ZZ(2), QQ(1))))
        spin_order = ZZ(4 if rank % 2 else 2)
        candidates.extend(
            ((spin_order, QQ(rank) / 4), (spin_order, QQ(rank) / 4))
        )
    elif invariant == (6, 72, 3):
        candidates.extend(((ZZ(1), QQ(0)), (ZZ(3), QQ(4) / 3)))
    elif invariant == (7, 126, 2):
        candidates.extend(((ZZ(1), QQ(0)), (ZZ(2), QQ(3) / 2)))
    elif invariant == (8, 240, 1):
        candidates.append((ZZ(1), QQ(0)))
    else:
        raise RuntimeError(f"unrecognized ADE component {invariant}")

    matches = {
        correction
        for class_order, correction in candidates
        if class_order == order and mod_two(correction) == mod_two(raw_norm)
    }
    assert len(matches) == 1, (invariant, order, raw_norm, matches)
    return matches.pop()


def section_data_for_lift(gram, section_height, lift, components=None):
    """Return exact component corrections and P.O for one section lift."""
    components = root_component_data(gram) if components is None else components
    corrections = []
    for invariant, component_basis in components:
        component_gram = component_basis * gram * component_basis.transpose()
        pairing = vector(QQ, lift) * gram * component_basis.transpose()
        dual_coordinates = pairing * component_gram.inverse()
        correction = ade_minimal_correction(
            invariant, dual_coordinates, component_gram
        )
        corrections.append((invariant, correction))
    total_correction = sum(value for _, value in corrections)
    pole = (section_height + total_correction - 4) / 2
    assert pole in ZZ and pole >= 0
    return tuple(corrections), ZZ(pole)


def section_pole_data(gram, height, lifts):
    """Return exact component corrections and P.O for MW basis lifts."""
    components = root_component_data(gram)
    rows = []
    for index, lift in enumerate(lifts.rows()):
        rows.append(
            section_data_for_lift(
                gram, height[index, index], lift, components=components
            )
        )
    return tuple(rows)


def optimal_section_pole_basis(gram, height, lifts):
    """Certify the least possible maximum P.O for an MW basis (rank <= 3)."""
    rank = height.nrows()
    if rank == 0 or rank > 3:
        return None
    components = root_component_data(gram)

    def pole(row):
        row = vector(ZZ, row)
        section_height = vector(QQ, row) * height * vector(QQ, row)
        lift = row * lifts
        return section_data_for_lift(
            gram, section_height, lift, components=components
        )[1]

    coordinate_profile = tuple(pole(row) for row in Matrix.identity(ZZ, rank).rows())
    upper = max(coordinate_profile)
    scale = ZZ(1)
    for value in height.list():
        scale = scale.lcm(ZZ(QQ(value).denominator()))
    scaled = (scale * height).change_ring(ZZ)
    if any(value % 2 for value in scaled.diagonal()):
        scale *= 2
        scaled *= 2
    form = QuadraticForm(ZZ, scaled)
    for bound in range(upper + 1):
        q_bound = (scale * (4 + 2 * bound)) // 2 + 1
        shells = form.short_vector_list_up_to_length(
            q_bound, up_to_sign_flag=True
        )
        short = []
        for shell in shells:
            for row in shell:
                row = vector(ZZ, row)
                if row and 0 <= pole(row) <= bound:
                    short.append(row)
        for rows in combinations(short, rank):
            basis = matrix(ZZ, [list(row) for row in rows])
            if abs(basis.det()) == 1:
                return (
                    bound,
                    tuple(pole(row) for row in basis.rows()),
                    matrix_rows(basis),
                )
    raise RuntimeError("coordinate basis failed to furnish a pole bound")


def fractional_root_class(gram, lift, component_basis):
    """Return the discriminant class of a frame lift on one ADE factor."""
    component_gram = component_basis * gram * component_basis.transpose()
    pairing = vector(QQ, lift) * gram * component_basis.transpose()
    coordinates = pairing * component_gram.inverse()
    return tuple(value - value.floor() for value in coordinates)


def class_multiple(multiplier, point):
    return tuple(
        multiplier * value - (multiplier * value).floor() for value in point
    )


def class_order(point, modulus):
    zero = tuple(QQ(0) for _ in point)
    return next(
        order
        for order in range(1, modulus + 1)
        if class_multiple(order, point) == zero
    )


def ade_component_name(invariant):
    rank, root_count, determinant = invariant
    if root_count == rank * (rank + 1) and determinant == rank + 1:
        return f"A{rank}"
    if root_count == 2 * rank * (rank - 1) and determinant == 4:
        return f"D{rank}"
    if invariant == (6, 72, 3):
        return "E6"
    if invariant == (7, 126, 2):
        return "E7"
    if invariant == (8, 240, 1):
        return "E8"
    raise RuntimeError(f"unrecognized ADE component {invariant}")


def cyclic_component_labels(classes, modulus):
    """Label a cyclic discriminant group, up to inversion of its generator."""
    zero = tuple(QQ(0) for _ in classes[0])
    generators = [
        point for point in classes if class_order(point, modulus) == modulus
    ]
    if generators:
        generator = generators[0]
        return tuple(
            next(
                multiplier
                for multiplier in range(modulus)
                if class_multiple(multiplier, generator) == point
            )
            for point in classes
        )
    # A set of section classes need not generate the whole component group.
    # For Z/4 the only remaining nonzero possibility is the unique order-two
    # class, whose label is independent of a choice of generator.
    assert modulus == 4
    return tuple(0 if point == zero else 2 for point in classes)


def ade_pair_correction(name, left, right):
    """Canonical local height pairing for cyclic ADE component labels."""
    if left == 0 or right == 0:
        return QQ(0)
    if name.startswith("A"):
        rank = ZZ(name[1:])
        order = rank + 1
        return QQ(min(left, right) * (order - max(left, right))) / order
    if name == "D5":
        table = (
            (QQ(5) / 4, QQ(1) / 2, QQ(3) / 4),
            (QQ(1) / 2, QQ(1), QQ(1) / 2),
            (QQ(3) / 4, QQ(1) / 2, QQ(5) / 4),
        )
        return table[left - 1][right - 1]
    if name == "D4":
        # The three nonzero classes of D4 form (Z/2)^2.  Each has norm one;
        # distinct nonzero triality classes pair to 1/2.
        return QQ(1) if left == right else QQ(1) / 2
    if name == "E6":
        return QQ(4) / 3 if left == right else QQ(2) / 3
    if name == "E7":
        assert left == right == 1
        return QQ(3) / 2
    if name == "E8":
        raise AssertionError("E8 has no nonzero component class")
    raise RuntimeError(f"missing local correction table for {name}")


def exact_section_profiles(gram, height, lifts, basis):
    """Recover component labels and all Shioda intersections for an MW basis."""
    basis = matrix(ZZ, basis)
    transformed_height = basis * height * basis.transpose()
    transformed_lifts = basis * lifts
    components = root_component_data(gram)
    names = tuple(ade_component_name(invariant) for invariant, _ in components)
    labels_by_component = []
    for (invariant, component_basis), name in zip(components, names):
        modulus = invariant[2]
        if modulus == 1:
            labels_by_component.append((0,) * transformed_lifts.nrows())
            continue
        classes = tuple(
            fractional_root_class(gram, lift, component_basis)
            for lift in transformed_lifts.rows()
        )
        if name == "D4":
            zero = tuple(QQ(0) for _ in classes[0])
            nonzero_classes = []
            labels = []
            for point in classes:
                if point == zero:
                    labels.append(0)
                    continue
                if point not in nonzero_classes:
                    nonzero_classes.append(point)
                labels.append(nonzero_classes.index(point) + 1)
            assert len(nonzero_classes) <= 3
            labels_by_component.append(tuple(labels))
        else:
            labels_by_component.append(cyclic_component_labels(classes, modulus))
    profiles = tuple(
        tuple(labels_by_component[column][row] for column in range(len(names)))
        for row in range(transformed_lifts.nrows())
    )

    def local(left, right):
        return sum(
            ade_pair_correction(name, left[index], right[index])
            for index, name in enumerate(names)
        )

    poles = tuple(
        (transformed_height[row, row] + local(profiles[row], profiles[row]) - 4)
        / 2
        for row in range(transformed_height.nrows())
    )
    assert all(value in ZZ and value >= 0 for value in poles)
    pair_intersections = []
    for row in range(transformed_height.nrows()):
        for column in range(row + 1, transformed_height.nrows()):
            intersection = (
                2
                + poles[row]
                + poles[column]
                - local(profiles[row], profiles[column])
                - transformed_height[row, column]
            )
            assert intersection in ZZ and intersection >= 0
            pair_intersections.append((row + 1, column + 1, ZZ(intersection)))
    return (
        names,
        transformed_height,
        profiles,
        tuple(ZZ(value) for value in poles),
        tuple(pair_intersections),
    )


def discriminant_generator(gram):
    smith, left, _ = gram.smith_form()
    diagonal = list(smith.diagonal())
    assert diagonal[:-1] == [1] * (gram.nrows() - 1)
    order = abs(ZZ(diagonal[-1]))
    last = vector(ZZ, [0] * (gram.nrows() - 1) + [1])
    integral_covector = left.inverse() * last
    dual_vector = gram.inverse() * integral_covector
    assert order * dual_vector in ZZ**gram.nrows()
    return order, dual_vector


def primitive_closure(gram, norm, disc_coefficient, rank1_coefficient):
    order, generator = discriminant_generator(gram)
    ambient = block_diagonal_matrix(gram, Matrix(ZZ, [[norm]]))
    glue = vector(
        QQ,
        list(disc_coefficient * generator)
        + [QQ(rank1_coefficient) / norm],
    )
    assert (glue * ambient * glue) % 2 == 0

    denominator = ZZ(glue.denominator())
    integer_generators = denominator * Matrix.identity(gram.nrows() + 1)
    integer_generators = integer_generators.stack(
        Matrix(ZZ, 1, gram.nrows() + 1, [ZZ(denominator * x) for x in glue])
    )
    basis = integer_generators.row_module().basis_matrix().change_ring(QQ) / denominator
    extension = basis * ambient * basis.transpose()
    assert all(x.denominator() == 1 for x in extension.list())
    extension = extension.change_ring(ZZ)
    index = ZZ(1 / abs(basis.det()))
    return extension, ZZ(index), order, glue, basis


def root_invariants(gram):
    form = QuadraticForm(ZZ, gram)
    short = form.short_vector_list_up_to_length(2, up_to_sign_flag=True)
    half_roots = short[1] if len(short) > 1 else []
    root_basis = Matrix(ZZ, [list(row) for row in half_roots]).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return (
        root_basis.rank(),
        2 * len(half_roots),
        abs(ZZ(root_gram.det())),
        root_basis,
    )


def root_component_data(gram):
    form = QuadraticForm(ZZ, gram)
    half_roots = [vector(ZZ, row) for row in form.short_vector_list_up_to_length(
        2, up_to_sign_flag=True
    )[1]]
    roots = half_roots + [-row for row in half_roots]
    parents = list(range(len(roots)))

    def find(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left, right):
        left = find(left)
        right = find(right)
        if left != right:
            parents[right] = left

    for i, root in enumerate(roots):
        for j in range(i):
            if root * gram * roots[j] != 0:
                union(i, j)

    groups = {}
    for index, root in enumerate(roots):
        groups.setdefault(find(index), []).append(root)

    components = []
    for group in groups.values():
        basis = matrix(ZZ, [list(root) for root in group]).row_module().basis_matrix()
        root_gram = basis * gram * basis.transpose()
        invariant = (basis.rank(), len(group), abs(ZZ(root_gram.det())))
        components.append((invariant, basis))
    return sorted(components, key=lambda item: item[0])


def root_components(gram):
    return [invariant for invariant, _ in root_component_data(gram)]


def root_torsion_order(root_basis):
    smith, _, _ = root_basis.smith_form()
    diagonal = [abs(ZZ(smith[i, i])) for i in range(root_basis.nrows())]
    return prod(value for value in diagonal if value)


def isotropic_mate(ns, fiber):
    pairing = list(ns * fiber)
    current = ZZ(0)
    coefficients = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(pairing):
        if value == 0:
            continue
        new_gcd, old_scale, new_scale = xgcd(current, ZZ(value))
        coefficients = [old_scale * x for x in coefficients]
        coefficients[index] += new_scale
        current = new_gcd
    assert abs(current) == 1
    if current == -1:
        coefficients = [-x for x in coefficients]
    mate = vector(ZZ, coefficients)
    assert fiber * ns * mate == 1
    square = ZZ(mate * ns * mate)
    assert square % 2 == 0
    mate -= (square // 2) * fiber
    assert mate * ns * mate == 0 and fiber * ns * mate == 1
    return mate


def build_neighbor_basis(frame, witness, a, b):
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    ns = block_diagonal_matrix(hyperbolic, -frame)
    witness = vector(ZZ, witness)
    q = ZZ(witness * frame * witness) // 2
    assert witness * frame * witness == 2 * q
    assert ZZ(a) * ZZ(b) == q
    fiber = vector(ZZ, [ZZ(a), ZZ(b)] + list(witness))
    assert fiber * ns * fiber == 0
    mate = isotropic_mate(ns, fiber)
    complement = matrix(ZZ, [list(fiber * ns), list(mate * ns)]).right_kernel_matrix()
    basis = matrix(ZZ, [list(fiber), list(mate)] + [list(row) for row in complement])
    assert abs(basis.det()) == 1
    child = -(complement * ns * complement.transpose())
    return ns, basis, child


def q60_neighbor_basis(frame, a=5, b=12):
    witness = vector(ZZ, [0, 0, -1, -1, -1, -1, -1,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
    assert witness * frame * witness == 120
    return build_neighbor_basis(frame, witness, a, b)


def mw_height_gram(gram, root_basis=None, return_lifts=False):
    """Recover the saturated Mordell--Weil height lattice from a frame."""
    if root_basis is None:
        root_basis = root_invariants(gram)[3]
    rank = gram.nrows()
    mw_rank = rank - root_basis.rank()
    if mw_rank == 0:
        empty = Matrix(QQ, 0, 0)
        if return_lifts:
            return empty, empty, Matrix(ZZ, 0, rank)
        return empty, empty

    orthogonal = (root_basis * gram).right_kernel_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    orthogonal_gram = orthogonal * gram * orthogonal.transpose()
    combined = root_basis.stack(orthogonal)
    index = abs(ZZ(combined.det()))
    inverse_combined = combined.inverse()

    representatives = {tuple([QQ(0)] * rank): vector(ZZ, [0] * rank)}
    queue = [vector(ZZ, [0] * rank)]

    def fractional_key(row):
        coordinates = vector(QQ, row) * inverse_combined
        return tuple(value - value.floor() for value in coordinates)

    representatives = {fractional_key(queue[0]): queue[0]}
    head = 0
    while head < len(queue) and len(representatives) < index:
        row = queue[head]
        head += 1
        for coordinate in range(rank):
            unit = vector(ZZ, [1 if i == coordinate else 0 for i in range(rank)])
            for sign in (1, -1):
                candidate = row + sign * unit
                key = fractional_key(candidate)
                if key not in representatives:
                    representatives[key] = candidate
                    queue.append(candidate)
                    if len(representatives) == index:
                        break
            if len(representatives) == index:
                break
    assert len(representatives) == index

    inverse_root_gram = root_gram.inverse()
    inverse_orthogonal_gram = orthogonal_gram.inverse()

    def project_to_mw(row):
        row = vector(QQ, row)
        root_coordinates = row * gram * root_basis.transpose() * inverse_root_gram
        projection = row - root_coordinates * root_basis
        return projection * gram * orthogonal.transpose() * inverse_orthogonal_gram

    representative_rows = list(representatives.values())
    projected = [project_to_mw(row) for row in representative_rows]
    generators = list(Matrix.identity(QQ, mw_rank).rows()) + projected
    denominator = ZZ(1)
    for row in generators:
        denominator = denominator.lcm(ZZ(vector(QQ, row).denominator()))
    integer_rows = Matrix(
        ZZ, [[ZZ(denominator * value) for value in row] for row in generators]
    )
    mw_basis = integer_rows.row_module().basis_matrix().change_ring(QQ) / denominator
    mw_lifts = []
    for basis_row in mw_basis.rows():
        for projected_row, representative_row in zip(
            projected, representative_rows
        ):
            difference = vector(QQ, basis_row) - vector(QQ, projected_row)
            if all(value.denominator() == 1 for value in difference):
                lift = vector(QQ, representative_row) + difference * orthogonal
                assert all(value.denominator() == 1 for value in lift)
                lift = vector(ZZ, lift)
                assert project_to_mw(lift) == vector(QQ, basis_row)
                mw_lifts.append(lift)
                break
        else:
            raise RuntimeError("failed to lift saturated MW basis")
    mw_lifts = matrix(ZZ, [list(row) for row in mw_lifts])
    height = mw_basis * orthogonal_gram * mw_basis.transpose()
    torsion = root_torsion_order(root_basis)
    assert height.det() == QQ(gram.det() * torsion**2) / QQ(root_gram.det())

    scale = ZZ(1)
    for value in height.list():
        scale = scale.lcm(ZZ(QQ(value).denominator()))
    integral_height = (scale * height).change_ring(ZZ)
    lll = Matrix(ZZ, integral_height.nrows(), integral_height.ncols(), pari(integral_height).qflllgram())
    candidates = []
    for transform in (lll, lll.transpose()):
        if abs(transform.det()) != 1:
            continue
        reduced = transform * height * transform.transpose()
        candidates.append(
            (max(abs(value) for value in reduced.list()), reduced, transform)
        )
    if candidates:
        _, reduced_height, reduced_transform = min(candidates, key=lambda item: item[0])
    else:
        reduced_height = height
        reduced_transform = Matrix.identity(ZZ, mw_rank)
    reduced_lifts = reduced_transform * mw_lifts
    if return_lifts:
        return height, reduced_height, reduced_lifts
    return height, reduced_height


def optimal_e6_pole_basis(height):
    """Find and certify the smallest possible maximum P.O for E6+E8/MW3."""
    assert height.nrows() == 3
    scaled = (3 * height).change_ring(ZZ)
    form = QuadraticForm(ZZ, scaled)

    def pole(row):
        row = vector(ZZ, row)
        section_height = vector(QQ, row) * height * vector(QQ, row)
        fractional = section_height - section_height.floor()
        if fractional == 0:
            local_correction = QQ(0)
        else:
            assert fractional == QQ(2) / 3
            local_correction = QQ(4) / 3
        value = (section_height - 4 + local_correction) / 2
        assert value in ZZ
        return ZZ(value)

    coordinate_profile = tuple(pole(row) for row in Matrix.identity(ZZ, 3).rows())
    upper = max(coordinate_profile)
    for bound in range(upper + 1):
        # q_form(row)=3*h(row)/2 and h<=4+2*bound.  The extra one makes
        # Sage's shell endpoint convention explicit.
        q_bound = ZZ(3 * (4 + 2 * bound) / 2) + 1
        shells = form.short_vector_list_up_to_length(
            q_bound, up_to_sign_flag=True
        )
        short = []
        for shell in shells:
            for row in shell:
                row = vector(ZZ, row)
                if row and 0 <= pole(row) <= bound:
                    short.append(row)
        for rows in combinations(short, 3):
            basis = matrix(ZZ, [list(row) for row in rows])
            if abs(basis.det()) == 1:
                profile = tuple(pole(row) for row in basis.rows())
                return bound, profile, basis
    raise RuntimeError("coordinate basis failed to furnish an E6 pole bound")


def dominant_weight_representatives(cartan, norm_bound, shift=None):
    """Enumerate dominant weights of bounded norm in one root-lattice coset."""
    inverse = cartan.inverse()
    rank = cartan.nrows()
    shift = vector(QQ, [0] * rank) if shift is None else vector(QQ, shift)
    representatives = []

    def recurse(prefix):
        padded = vector(ZZ, prefix + [0] * (rank - len(prefix)))
        partial_norm = padded * inverse * padded
        if partial_norm > norm_bound:
            return
        if len(prefix) == rank:
            root_coordinates = vector(QQ, prefix) * inverse - shift
            if all(value.denominator() == 1 for value in root_coordinates):
                representatives.append(
                    (vector(ZZ, prefix), partial_norm, vector(ZZ, root_coordinates))
                )
            return
        value = 0
        while True:
            trial = prefix + [value]
            padded_trial = vector(ZZ, trial + [0] * (rank - len(trial)))
            if padded_trial * inverse * padded_trial > norm_bound:
                break
            recurse(trial)
            value += 1

    recurse([])
    return representatives


def level79_weyl_orbits(max_q):
    """Exact frame-vector orbits with nonzero level-79 coordinate and q<237."""
    assert 60 <= max_q < 237
    e7 = CartanMatrix(["E", 7])
    e8 = CartanMatrix(["E", 8])
    minuscule = vector(QQ, [0, 0, 0, 0, 0, 0, 1]) * e7.inverse()
    height79 = QQ(237) / 2
    budget = 2 * ZZ(max_q) - height79
    e7_weights = dominant_weight_representatives(e7, budget, minuscule)
    e8_weights = dominant_weight_representatives(e8, budget)
    by_q = {}
    seen = set()
    for _, norm7, root7 in e7_weights:
        for _, norm8, root8 in e8_weights:
            root_norm = norm7 + norm8
            if root_norm > budget:
                continue
            section4 = ZZ(0)
            while root_norm + 4 * section4**2 <= budget:
                norm = height79 + root_norm + 4 * section4**2
                assert norm / 2 in ZZ
                q = ZZ(norm / 2)
                witness = vector(
                    ZZ, list(root7) + list(root8) + [section4, 1]
                )
                key = tuple(witness)
                if key not in seen:
                    seen.add(key)
                    by_q.setdefault(q, []).append(witness)
                section4 += 1
    return by_q


def proper_factor_pairs(q):
    return tuple(
        (ZZ(a), ZZ(q // a))
        for a in range(2, ZZ(q).isqrt() + 1)
        if q % a == 0 and q // a > 1
    )


def transport_neighbor_to_cm(
    generic_neighbor_basis, closure_basis, cm_frame, return_basis=False
):
    """Keep the q=60 fiber class after adjoining the CM algebraic class."""
    inverse_closure_basis = closure_basis.inverse()

    def lift(row):
        ambient_frame = vector(QQ, list(row[2:]) + [0])
        cm_coordinates = ambient_frame * inverse_closure_basis
        lifted = vector(QQ, list(row[:2]) + list(cm_coordinates))
        assert all(value.denominator() == 1 for value in lifted)
        return vector(ZZ, lifted)

    lifted_generic_basis = [lift(row) for row in generic_neighbor_basis.rows()]
    fiber = lifted_generic_basis[0]
    mate = lifted_generic_basis[1]
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    cm_ns = block_diagonal_matrix(hyperbolic, -cm_frame)
    assert fiber * cm_ns * fiber == 0
    assert mate * cm_ns * mate == 0
    assert fiber * cm_ns * mate == 1
    complement = matrix(
        ZZ, [list(fiber * cm_ns), list(mate * cm_ns)]
    ).right_kernel_matrix()
    enhanced = -(complement * cm_ns * complement.transpose())
    assert enhanced.is_positive_definite()
    cm_neighbor_basis = matrix(
        ZZ, [list(fiber), list(mate)] + [list(row) for row in complement]
    )
    assert abs(cm_neighbor_basis.det()) == 1
    inverse_cm_neighbor_basis = cm_neighbor_basis.inverse()
    generic_child_embedding_rows = []
    for row in lifted_generic_basis[2:]:
        coordinates = row * inverse_cm_neighbor_basis
        assert coordinates[0] == coordinates[1] == 0
        generic_child_embedding_rows.append(list(coordinates[2:]))
    generic_child_embedding = matrix(ZZ, generic_child_embedding_rows)
    if return_basis:
        return enhanced, generic_child_embedding, cm_neighbor_basis
    return enhanced, generic_child_embedding


frame = load_gram(BASE / "data/fibrations/kumar_e7e8_mw2_frame_2.txt")
assert frame.det() == 948
order, generator = discriminant_generator(frame)
assert order == 948
assert generator * frame * generator == QQ(245) / 948

# Delta=-3: glue (3*55 in Z/948, 1 in Z/316), of order 316.
# The other normalized isotropic choices are equivalent under sign and the
# H2 discriminant action 475.
cm3, index3, _, glue3, closure_basis3 = primitive_closure(frame, 316, 3 * 55, 1)
assert index3 == 316
assert cm3.det() == 3
roots3 = root_invariants(cm3)
assert roots3[:3] == (18, 486, 3)  # E8 + E8 + A2
print(
    f"KUMARCM|Delta=-3|index={index3}|det={cm3.det()}|"
    f"root_rank={roots3[0]}|roots={roots3[1]}|root_det={roots3[2]}",
    flush=True,
)

# Delta=-24: after normalizing the rank-one component to 2 mod 158, the
# isotropic ratio is +/-29 in Z/79; the two signs give isometric closures.
cm24, index24, _, glue24, closure_basis24 = primitive_closure(frame, 158, 12 * 29, 2)
assert index24 == 79
assert cm24.det() == 24
roots24 = root_invariants(cm24)
assert roots24[:3] == (16, 480, 1)  # E8 + E8
mw24_basis = (roots24[3] * cm24).right_kernel_matrix()
mw24 = mw24_basis * cm24 * mw24_basis.transpose()
assert mw24.det() == 24
print(
    f"KUMARCM|Delta=-24|index={index24}|det={cm24.det()}|"
    f"root_rank={roots24[0]}|roots={roots24[1]}|root_det={roots24[2]}|"
    f"mw_gram={mw24}",
    flush=True,
)

# Delta=-43: the Gross vector -11*i-j-5*k transports to the primitive
# generic-transcendental vector (169,167,-128), of norm -40764 and
# divisibility 948.  Hence the algebraic rank-one summand has norm 40764
# and the primitive closure has index 948.  The normalized isotropic glue
# below is one of the eight Atkin--Lehner-related solutions of
#
#     245*k^2 + 43 == 0 (mod 2*948).
cm43, index43, _, glue43, closure_basis43 = primitive_closure(
    frame, 40764, 53, 43
)
assert index43 == 948
assert cm43.det() == 43
roots43 = root_invariants(cm43)
components43 = root_components(cm43)
height43, reduced_height43, reduced_lifts43 = mw_height_gram(
    cm43, roots43[3], return_lifts=True
)
cm43_profile_data = exact_section_profiles(
    cm43,
    reduced_height43,
    reduced_lifts43,
    Matrix.identity(ZZ, 3),
)
assert cm43_profile_data[3] == (0, 0, 0)
print(
    f"KUMARCM|Delta=-43|index={index43}|det={cm43.det()}|"
    f"root_rank={roots43[0]}|roots={roots43[1]}|root_det={roots43[2]}|"
    f"components={components43}|MW={18-roots43[0]}|"
    f"mw_reduced={reduced_height43}",
    flush=True,
)

# The same norm-120 witness has five proper factor presentations a*b=60.
# They define different U embeddings even though every generic child has the
# same E8+E6 root data.  Compare them before privileging an equation chart.
factor_presentations = ((2, 30), (3, 20), (4, 15), (5, 12), (6, 10))
presentation_data = {}
for factor_a, factor_b in factor_presentations:
    generic_ns_i, neighbor_basis_i, child_i = q60_neighbor_basis(
        frame, factor_a, factor_b
    )
    roots_i = root_invariants(child_i)
    components_i = root_components(child_i)
    torsion_i = root_torsion_order(roots_i[3])
    _, reduced_height_i, reduced_lifts_i = mw_height_gram(
        child_i, roots_i[3], return_lifts=True
    )
    pole_bound_i, pole_profile_i, pole_basis_i = optimal_e6_pole_basis(
        reduced_height_i
    )
    assert child_i.det() == 948
    assert components_i == [(6, 72, 3), (8, 240, 1)]
    assert torsion_i == 1

    cm3_i, embedding3_i = transport_neighbor_to_cm(
        neighbor_basis_i, closure_basis3, cm3
    )
    cm24_i, embedding24_i = transport_neighbor_to_cm(
        neighbor_basis_i, closure_basis24, cm24
    )
    cm43_i, embedding43_i = transport_neighbor_to_cm(
        neighbor_basis_i, closure_basis43, cm43
    )
    assert child_i == embedding3_i * cm3_i * embedding3_i.transpose()
    assert child_i == embedding24_i * cm24_i * embedding24_i.transpose()
    assert child_i == embedding43_i * cm43_i * embedding43_i.transpose()

    roots3_i = root_invariants(cm3_i)
    roots24_i = root_invariants(cm24_i)
    roots43_i = root_invariants(cm43_i)
    components3_i = root_components(cm3_i)
    components24_i = root_components(cm24_i)
    components43_i = root_components(cm43_i)
    _, reduced_height3_i = mw_height_gram(cm3_i, roots3_i[3])
    _, reduced_height24_i = mw_height_gram(cm24_i, roots24_i[3])
    _, height43_i, lifts43_i = mw_height_gram(
        cm43_i, roots43_i[3], return_lifts=True
    )
    optimal43_i = optimal_section_pole_basis(cm43_i, height43_i, lifts43_i)
    generic_optimal_lifts_i = matrix(ZZ, pole_basis_i) * reduced_lifts_i
    cm43_root_gram_inverse_i = (
        roots43_i[3] * cm43_i * roots43_i[3].transpose()
    ).inverse()
    specialization43_rows_i = []
    for generic_lift in generic_optimal_lifts_i.rows():
        cm_lift = generic_lift * embedding43_i
        root_coordinates = (
            cm_lift * cm43_i * roots43_i[3].transpose()
            * cm43_root_gram_inverse_i
        )
        cm_projection = cm_lift - root_coordinates * roots43_i[3]
        pairings = cm_projection * cm43_i * lifts43_i.transpose()
        coordinates = pairings * height43_i.inverse()
        assert all(value in ZZ for value in coordinates)
        specialization43_rows_i.append(list(map(ZZ, coordinates)))
    specialization43_i = matrix(ZZ, specialization43_rows_i)
    generic_at_cm43_i = exact_section_profiles(
        cm43_i, height43_i, lifts43_i, specialization43_i
    )

    presentation_data[(factor_a, factor_b)] = (
        generic_ns_i,
        neighbor_basis_i,
        child_i,
        cm3_i,
        embedding3_i,
        cm24_i,
        embedding24_i,
        cm43_i,
        embedding43_i,
    )
    print(
        f"KUMARCMQ60PRESENTATION|ab={factor_a},{factor_b}|stage=generic|"
        f"components={components_i}|MW={child_i.nrows()-roots_i[0]}|"
        f"torsion={torsion_i}|mw_reduced={reduced_height_i}|"
        f"optimal_max_PO={pole_bound_i}|pole_profile={pole_profile_i}|"
        f"pole_basis={pole_basis_i}",
        flush=True,
    )
    print(
        f"KUMARCMQ60PRESENTATION|ab={factor_a},{factor_b}|Delta=-3|"
        f"components={components3_i}|MW={18-roots3_i[0]}|"
        f"root_jump={roots3_i[0]-roots_i[0]}|mw_reduced={reduced_height3_i}",
        flush=True,
    )
    print(
        f"KUMARCMQ60PRESENTATION|ab={factor_a},{factor_b}|Delta=-24|"
        f"components={components24_i}|MW={18-roots24_i[0]}|"
        f"root_jump={roots24_i[0]-roots_i[0]}|mw_reduced={reduced_height24_i}",
        flush=True,
    )
    print(
        f"KUMARCMQ60PRESENTATION|ab={factor_a},{factor_b}|Delta=-43|"
        f"components={components43_i}|MW={18-roots43_i[0]}|"
        f"root_jump={roots43_i[0]-roots_i[0]}|mw_reduced={height43_i}|"
        f"optimal_poles={optimal43_i}|generic_at_cm43={generic_at_cm43_i}",
        flush=True,
    )

# Retain the detailed standard-frame certificate for the originally pinned
# (5,12) presentation.
(
    generic_ns,
    neighbor_basis,
    child,
    child_cm3,
    child_embedding3,
    child_cm24,
    child_embedding24,
    child_cm43,
    child_embedding43,
) = presentation_data[(5, 12)]
pinned_child = load_gram(BASE / "data/fibrations/kumar_q60_e8_e6_mw3_frame.txt")
assert child == pinned_child
assert child_cm3.det() == 3
assert child == child_embedding3 * child_cm3 * child_embedding3.transpose()
child_roots3 = root_invariants(child_cm3)
child_components3 = root_components(child_cm3)
assert child_components3 == [(2, 6, 3), (8, 240, 1), (8, 240, 1)]
print(
    f"KUMARCMQ60|Delta=-3|det={child_cm3.det()}|"
    f"root_rank={child_roots3[0]}|roots={child_roots3[1]}|"
    f"root_det={child_roots3[2]}|components={child_components3}|"
    f"MW={18-child_roots3[0]}",
    flush=True,
)

assert child_cm24.det() == 24
assert child == child_embedding24 * child_cm24 * child_embedding24.transpose()
child_roots24 = root_invariants(child_cm24)
child_components24 = root_components(child_cm24)
assert child_components24 == [(2, 6, 3), (7, 126, 2), (8, 240, 1)]
child_torsion24 = root_torsion_order(child_roots24[3])
assert child_torsion24 == 1
child_mw_height24 = QQ(child_cm24.det() * child_torsion24**2) / child_roots24[2]
assert child_mw_height24 == 4
standard_child24 = block_diagonal_matrix(
    CartanMatrix(["E", 8]),
    CartanMatrix(["E", 7]),
    CartanMatrix(["A", 2]),
    Matrix(ZZ, [[4]]),
)
assert standard_child24.det() == 24
standard_isometry24 = QuadraticForm(ZZ, child_cm24).is_globally_equivalent_to(
    QuadraticForm(ZZ, standard_child24), return_matrix=True
)
assert standard_isometry24 is not False
assert (
    standard_isometry24.transpose() * child_cm24 * standard_isometry24
    == standard_child24
)

# Express the original Kumar fiber in the new q=60 U plus the concrete
# E8+E7+A2+<4> frame at the discriminant-24 endpoint.
old_fiber = vector(ZZ, [1] + [0] * 18)
old_fiber_in_generic_neighbor = old_fiber * neighbor_basis.inverse()
old_fiber_in_cm_child = vector(
    ZZ,
    list(old_fiber_in_generic_neighbor[:2])
    + list(old_fiber_in_generic_neighbor[2:] * child_embedding24),
)
standard_frame_coordinates = (
    old_fiber_in_cm_child[2:] * standard_isometry24.inverse().transpose()
)
old_fiber_in_standard_child = vector(
    ZZ,
    list(old_fiber_in_cm_child[:2]) + list(standard_frame_coordinates),
)
standard_cm_ns = block_diagonal_matrix(
    matrix(ZZ, [[0, 1], [1, 0]]), -standard_child24
)
assert old_fiber_in_standard_child * standard_cm_ns * old_fiber_in_standard_child == 0
print(
    f"KUMARCMQ60|Delta=-24|det={child_cm24.det()}|"
    f"root_rank={child_roots24[0]}|roots={child_roots24[1]}|"
    f"root_det={child_roots24[2]}|components={child_components24}|"
    f"MW={18-child_roots24[0]}|mw_height={child_mw_height24}|"
    f"torsion={child_torsion24}|standard_frame_isometric=1",
    flush=True,
)
print(
    "KUMARCMQ60|Delta=-24|inverse_kumar_fiber_standard_coordinates={}".format(
        tuple(old_fiber_in_standard_child)
    ),
    flush=True,
)

# At the discriminant-43 point the selected q=60 child returns to an
# E7+E8/MW3 frame in the same integral isometry class as the Kumar closure.
# This supplies a precise lattice gate for testing whether the geometric
# neighbor becomes a genuine surface automorphism at CM43.
assert child_cm43.det() == 43
assert child == child_embedding43 * child_cm43 * child_embedding43.transpose()
(
    child_cm43_gate,
    child_embedding43_gate,
    child_cm43_neighbor_basis,
) = transport_neighbor_to_cm(
    neighbor_basis, closure_basis43, cm43, return_basis=True
)
assert child_cm43_gate == child_cm43
assert child_embedding43_gate == child_embedding43
child_roots43 = root_invariants(child_cm43)
assert root_components(child_cm43) == [(7, 126, 2), (8, 240, 1)]
cm43_return_isometry = QuadraticForm(
    ZZ, child_cm43
).is_globally_equivalent_to(
    QuadraticForm(ZZ, cm43), return_matrix=True
)
assert cm43_return_isometry is not False
assert abs(cm43_return_isometry.det()) == 1
assert (
    cm43_return_isometry.transpose()
    * child_cm43
    * cm43_return_isometry
    == cm43
)
print(
    "KUMARCMQ60|Delta=-43|ab=5,12|components=E7+E8|MW=3|"
    "returns_to_Kumar_frame_isometry=1|isometry_det={}|"
    "q60_fiber_in_Kumar_NS={}".format(
        cm43_return_isometry.det(),
        tuple(child_cm43_neighbor_basis[0]),
    ),
    flush=True,
)

# Express the q=60 fiber in a saturated explicit divisor basis
# [F,O,root-basis(15),P1,P2,P3].  The three sections are the reduced CM43 MW
# basis and all have P.O=0, so their NS representatives are (1,1,lift).
cm43_ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -cm43)
explicit_rows43 = [
    vector(ZZ, [1, 0] + [0]*18),
    vector(ZZ, [-1, 1] + [0]*18),
]
explicit_rows43.extend(
    vector(ZZ, [0, 0] + list(root)) for root in roots43[3].rows()
)
explicit_rows43.extend(
    vector(ZZ, [1, 1] + list(lift)) for lift in reduced_lifts43.rows()
)
explicit_basis43 = matrix(ZZ, [list(row) for row in explicit_rows43])
assert abs(explicit_basis43.det()) == 1
explicit_gram43 = explicit_basis43 * cm43_ns * explicit_basis43.transpose()
q60_fiber_explicit43 = (
    child_cm43_neighbor_basis[0] * explicit_basis43.inverse()
)
assert all(value in ZZ for value in q60_fiber_explicit43)
q60_fiber_explicit43 = vector(ZZ, q60_fiber_explicit43)
assert q60_fiber_explicit43 * explicit_gram43 * q60_fiber_explicit43 == 0
print(
    "KUMARCMQ60|Delta=-43|explicit_basis=F,O,roots15,P1,P2,P3|"
    "fiber={}|horizontal={}".format(
        tuple(q60_fiber_explicit43),
        tuple(q60_fiber_explicit43[-3:]),
    ),
    flush=True,
)

# The explicit CM-43 equation fixes more marking than the abstract closure:
# in its geometric basis the height-four and level-79 sections are
#
#     P3=(0,0,1),  Q79=(4,-5,1).
#
# There are eight normalized glue signs for the rank-one CM extension.  Audit
# them all and retain the one(s) whose reduced height basis carries the q=60
# horizontal projection to this exact Q79 vector.  This prevents an
# unmarked isometry of the discriminant-43 frame from being used as an
# equation-level identification.
marked_height43 = matrix(QQ, (
    (QQ(5)/2, -QQ(1)/2, -1),
    (-QQ(1)/2, QQ(5)/2, 0),
    (-1, 0, 4),
))
marked_height4 = vector(ZZ, (0, 0, 1))
marked_q79 = vector(ZZ, (4, -5, 1))
marked_q8_witness = vector(ZZ, (
    156, -78, 0, 0, -78, 0, -78, 0, 0,
    0, 0, 0, 0, 0, 0, -1, -155, -32,
))
glue43_candidates = tuple(
    k for k in range(1, 948)
    if (245*k*k + 43) % (2*948) == 0
)
assert len(glue43_candidates) == 8
marked_glue_hits = []
marked_glue_frames = {}
for glue_candidate in glue43_candidates:
    cm_candidate, candidate_index, _, _, candidate_closure_basis = (
        primitive_closure(frame, 40764, glue_candidate, 43)
    )
    assert candidate_index == 948 and cm_candidate.det() == 43
    candidate_roots = root_invariants(cm_candidate)
    assert root_components(cm_candidate) == [(7, 126, 2), (8, 240, 1)]
    _, candidate_height, candidate_lifts = mw_height_gram(
        cm_candidate, candidate_roots[3], return_lifts=True
    )
    _, _, candidate_neighbor_basis = transport_neighbor_to_cm(
        neighbor_basis,
        candidate_closure_basis,
        cm_candidate,
        return_basis=True,
    )
    candidate_ns = block_diagonal_matrix(
        matrix(ZZ, [[0, 1], [1, 0]]), -cm_candidate
    )
    candidate_explicit_rows = [
        vector(ZZ, [1, 0] + [0]*18),
        vector(ZZ, [-1, 1] + [0]*18),
    ]
    candidate_explicit_rows.extend(
        vector(ZZ, [0, 0] + list(root))
        for root in candidate_roots[3].rows()
    )
    candidate_explicit_rows.extend(
        vector(ZZ, [1, 1] + list(lift))
        for lift in candidate_lifts.rows()
    )
    candidate_explicit_basis = matrix(
        ZZ, [list(row) for row in candidate_explicit_rows]
    )
    assert abs(candidate_explicit_basis.det()) == 1
    candidate_fiber = vector(
        ZZ, candidate_neighbor_basis[0]*candidate_explicit_basis.inverse()
    )
    candidate_horizontal = vector(ZZ, candidate_fiber[-3:])
    assert candidate_horizontal*candidate_height*candidate_horizontal == QQ(237)/2
    candidate_root_gram_inverse = (
        candidate_roots[3]*cm_candidate*candidate_roots[3].transpose()
    ).inverse()

    def candidate_cm_mw_coordinates(cm_vector):
        cm_vector = vector(QQ, cm_vector)
        root_coordinates = (
            cm_vector*cm_candidate*candidate_roots[3].transpose()
            * candidate_root_gram_inverse
        )
        projection = cm_vector-root_coordinates*candidate_roots[3]
        pairings = projection*cm_candidate*candidate_lifts.transpose()
        coordinates = pairings*candidate_height.inverse()
        assert all(value in ZZ for value in coordinates)
        return vector(ZZ, coordinates)

    def candidate_mw_coordinates(generic_vector):
        ambient_vector = vector(QQ, list(generic_vector)+[0])
        cm_vector = ambient_vector*candidate_closure_basis.inverse()
        assert all(value in ZZ for value in cm_vector)
        return candidate_cm_mw_coordinates(vector(ZZ, cm_vector))

    def candidate_frame_vector(generic_vector):
        ambient_vector = vector(QQ, list(generic_vector)+[0])
        cm_vector = ambient_vector*candidate_closure_basis.inverse()
        assert all(value in ZZ for value in cm_vector)
        return vector(ZZ, cm_vector)

    generic_height4 = vector(ZZ, [0]*15+[1, 0])
    generic_q79 = vector(ZZ, [0]*16+[1])
    candidate_height4 = candidate_mw_coordinates(generic_height4)
    candidate_q79 = candidate_mw_coordinates(generic_q79)
    assert candidate_q79 == candidate_horizontal

    # PARI returns C with C^t*H_candidate*C=H_marked.  Thus row
    # coordinates map by (C^t)^(-1).  The target height lattice has only
    # global sign as an automorphism, so the ordered pair (P3,Q79) fixes the
    # remaining sign completely.
    marking_isometry = matrix(
        ZZ,
        pari(2*candidate_height).qfisom(pari(2*marked_height43)),
    )
    assert (
        marking_isometry.transpose()*(2*candidate_height)*marking_isometry
        == 2*marked_height43
    )
    mapped_height4 = vector(
        ZZ, candidate_height4*marking_isometry.inverse().transpose()
    )
    mapped_q79 = vector(
        ZZ, candidate_q79*marking_isometry.inverse().transpose()
    )
    if mapped_height4 == -marked_height4 and mapped_q79 == -marked_q79:
        mapped_height4 = -mapped_height4
        mapped_q79 = -mapped_q79
        marking_isometry = -marking_isometry
    mapped_q8_horizontal = None
    if glue_candidate == 211:
        candidate_q8_horizontal = candidate_cm_mw_coordinates(marked_q8_witness)
        mapped_q8_horizontal = vector(
            ZZ,
            candidate_q8_horizontal*marking_isometry.inverse().transpose(),
        )
        assert mapped_q8_horizontal*marked_height43*mapped_q8_horizontal <= 16
    marked = mapped_height4 == marked_height4 and mapped_q79 == marked_q79
    if marked:
        marked_glue_hits.append((glue_candidate, tuple(candidate_fiber)))
        marked_glue_frames[glue_candidate] = (
            cm_candidate,
            vector(ZZ, candidate_neighbor_basis[0]),
        )
    print(
        "KUMARCM43MARKING|glue={}|height={}|height4={}|horizontal={}"
        "|mapped_height4={}|mapped_Q79={}|marked_Q79={}".format(
            glue_candidate,
            candidate_height,
            tuple(candidate_height4),
            tuple(candidate_horizontal),
            tuple(mapped_height4),
            tuple(mapped_q79),
            int(marked),
        ),
        flush=True,
    )
    if mapped_q8_horizontal is not None:
        marked_height4_frame = candidate_frame_vector(generic_height4)
        marked_q79_frame = candidate_frame_vector(generic_q79)
        print(
            f"KUMARCM43MARKING|glue=211|q8_witness={tuple(marked_q8_witness)}"
            f"|q8_marked_horizontal={tuple(mapped_q8_horizontal)}"
            f"|q8_horizontal_height="
            f"{mapped_q8_horizontal*marked_height43*mapped_q8_horizontal}",
            flush=True,
        )
        print(
            f"KUMARCM43MARKING|glue=211"
            f"|height4_frame={tuple(marked_height4_frame)}"
            f"|q79_frame={tuple(marked_q79_frame)}",
            flush=True,
        )
assert len(marked_glue_hits) in (2, 4)
print(
    f"KUMARCM43MARKING|Q79_glues={tuple(hit[0] for hit in marked_glue_hits)}"
    f"|representative_fiber={marked_glue_hits[0][1]}|status=PASS",
    flush=True,
)

# The least-degenerate compact chart found in the exact bounded search is the
# q=80, (4,20) presentation below.  Generically it is E6+D5+A3/MW3.  At the
# discriminant-24 anchor it only acquires one A1 root, retains all three MW
# directions, and every specialized basis section has P.O=0.
q80_witness = vector(
    ZZ, (2, 3, 4, 6, 5, 4, 3, 6, 9, 12, 18, 15, 12, 8, 4, 2, 1)
)
_, q80_basis, q80_child = build_neighbor_basis(frame, q80_witness, 4, 20)
q80_pinned = load_gram(
    BASE / "data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt"
)
assert q80_child == q80_pinned

if args.print_markings_in_q80:
    for marking, old_frame_vector in (
        ("height4", vector(ZZ, [0]*15+[1, 0])),
        ("Q79", vector(ZZ, [0]*15+[0, 1])),
    ):
        old_ns_vector = vector(ZZ, [0, 0]+list(old_frame_vector))
        marking_in_q80 = vector(ZZ, old_ns_vector*q80_basis.inverse())
        print(
            f"KUMARCMQ80NSMARKING|marking={marking}|"
            f"coordinates={tuple(marking_in_q80)}",
            flush=True,
        )

if args.print_q60_in_q80:
    q80_ns = block_diagonal_matrix(
        matrix(ZZ, [[0, 1], [1, 0]]), -q80_child
    )
    for factor_presentation in factor_presentations:
        q60_basis_i = presentation_data[factor_presentation][1]
        q60_fiber_in_q80 = vector(
            ZZ, q60_basis_i.row(0) * q80_basis.inverse()
        )
        if q60_fiber_in_q80[0] < 0 and q60_fiber_in_q80[1] < 0:
            q60_fiber_in_q80 = -q60_fiber_in_q80
        assert q60_fiber_in_q80*q80_ns*q60_fiber_in_q80 == 0
        direct_q = q60_fiber_in_q80[0]*q60_fiber_in_q80[1]
        assert direct_q >= 0
        print(
            f"KUMARCMQ60INQ80|source_ab={factor_presentation[0]},"
            f"{factor_presentation[1]}|direct_ab="
            f"{q60_fiber_in_q80[0]},{q60_fiber_in_q80[1]}|q={direct_q}|"
            f"coordinates={tuple(q60_fiber_in_q80)}",
            flush=True,
        )
q80_roots = root_invariants(q80_child)
q80_components = root_components(q80_child)
assert q80_components == [(3, 12, 4), (5, 40, 4), (6, 72, 3)]
_, q80_height, q80_lifts = mw_height_gram(
    q80_child, q80_roots[3], return_lifts=True
)
q80_optimal = optimal_section_pole_basis(q80_child, q80_height, q80_lifts)
assert q80_optimal[:2] == (3, (0, 0, 3))
q80_profiles = exact_section_profiles(
    q80_child, q80_height, q80_lifts, q80_optimal[2]
)
assert q80_profiles[3] == q80_optimal[1]

q80_cm24, q80_embedding24 = transport_neighbor_to_cm(
    q80_basis, closure_basis24, cm24
)
assert q80_child == q80_embedding24 * q80_cm24 * q80_embedding24.transpose()
q80_cm24_roots = root_invariants(q80_cm24)
q80_cm24_components = root_components(q80_cm24)
assert q80_cm24_components == [
    (1, 2, 2), (3, 12, 4), (5, 40, 4), (6, 72, 3)
]
_, q80_cm24_height, q80_cm24_lifts = mw_height_gram(
    q80_cm24, q80_cm24_roots[3], return_lifts=True
)
q80_cm24_optimal = optimal_section_pole_basis(
    q80_cm24, q80_cm24_height, q80_cm24_lifts
)
assert q80_cm24_optimal[:2] == (0, (0, 0, 0))
q80_cm24_profiles = exact_section_profiles(
    q80_cm24, q80_cm24_height, q80_cm24_lifts, q80_cm24_optimal[2]
)
assert q80_cm24_profiles[3] == q80_cm24_optimal[1]

q80_cm43, q80_embedding43 = transport_neighbor_to_cm(
    q80_basis, closure_basis43, cm43
)
assert q80_child == q80_embedding43 * q80_cm43 * q80_embedding43.transpose()
q80_cm43_roots = root_invariants(q80_cm43)
q80_cm43_components = root_components(q80_cm43)
_, q80_cm43_height, q80_cm43_lifts = mw_height_gram(
    q80_cm43, q80_cm43_roots[3], return_lifts=True
)
q80_cm43_optimal = optimal_section_pole_basis(
    q80_cm43, q80_cm43_height, q80_cm43_lifts
)
q80_cm43_profiles = exact_section_profiles(
    q80_cm43, q80_cm43_height, q80_cm43_lifts, q80_cm43_optimal[2]
)
assert q80_cm43_profiles[3] == q80_cm43_optimal[1]

# Track the generic optimal basis into the CM24 frame.  Pairing with a
# saturated CM MW basis recovers its exact MW coordinates without making any
# choice of root representatives.  This decides which polynomial CM sections
# actually continue away from the CM divisor.
q80_generic_optimal_lifts = matrix(ZZ, q80_optimal[2]) * q80_lifts

# Track the two geometrically marked Kumar MW directions through the q=80
# neighbor into the optimal q80 MW basis.  The old NS basis is U plus the
# positive-definite frame; q80_basis expresses the new U and frame in it.
# Pairing in the new frame quotient removes all dependence on root lifts.
def old_kumar_mw_to_q80_optimal(old_frame_vector):
    old_ns_vector = vector(QQ, [0, 0]+list(old_frame_vector))
    new_ns_coordinates = old_ns_vector*q80_basis.inverse()
    new_frame_vector = vector(QQ, new_ns_coordinates[2:])
    q80_root_gram_inverse = (
        q80_roots[3]*q80_child*q80_roots[3].transpose()
    ).inverse()
    root_coordinates = (
        new_frame_vector*q80_child*q80_roots[3].transpose()
        * q80_root_gram_inverse
    )
    projection = new_frame_vector-root_coordinates*q80_roots[3]
    reduced_pairings = projection*q80_child*q80_lifts.transpose()
    reduced_coordinates = reduced_pairings*q80_height.inverse()
    optimal_coordinates = (
        reduced_coordinates*matrix(ZZ, q80_optimal[2]).inverse()
    )
    assert all(value in ZZ for value in optimal_coordinates)
    return vector(ZZ, optimal_coordinates)


q80_marked_height4 = old_kumar_mw_to_q80_optimal(
    vector(ZZ, [0]*15+[1, 0])
)
q80_marked_q79 = old_kumar_mw_to_q80_optimal(
    vector(ZZ, [0]*15+[0, 1])
)
assert q80_marked_height4*q80_profiles[1]*q80_marked_height4 == 4
# Section height is fibration-dependent.  The old Kumar Q79 class has height
# 237/2 in E7+E8, but height 120 after the q80 neighbor.
assert q80_marked_q79*q80_profiles[1]*q80_marked_q79 == 120
assert q80_marked_height4*q80_profiles[1]*q80_marked_q79 == 0
def q80_profile_of_combination(coordinates):
    labels = []
    for component_index, name in enumerate(q80_profiles[0]):
        if name.startswith("A"):
            modulus = ZZ(name[1:])+1
        elif name == "D5":
            modulus = 4
        elif name == "E6":
            modulus = 3
        else:
            raise AssertionError(f"unexpected q80 component {name}")
        labels.append(ZZ(sum(
            coordinates[index]*q80_profiles[2][index][component_index]
            for index in range(3)
        )) % modulus)
    return tuple(labels)


q80_marked_component_profiles = tuple(map(
    q80_profile_of_combination,
    (q80_marked_height4, q80_marked_q79),
))
q80_marked_local = lambda left, right: sum(
    ade_pair_correction(name, left[index], right[index])
    for index, name in enumerate(q80_profiles[0])
)
q80_marked_poles = tuple(
    ZZ((height+q80_marked_local(profile, profile)-4)/2)
    for height, profile in zip((QQ(4), QQ(120)), q80_marked_component_profiles)
)
q80_marked_intersection = ZZ(
    2+sum(q80_marked_poles)
    - q80_marked_local(*q80_marked_component_profiles)
)
print(
    f"KUMARCMQ80MARKING|basis=generic_optimal|"
    f"height4={tuple(q80_marked_height4)}|Q79={tuple(q80_marked_q79)}|"
    f"q80_heights=4,120|profiles={q80_marked_component_profiles}|"
    f"P.O={q80_marked_poles}|height4.Q79={q80_marked_intersection}|"
    "pair=0|status=PASS",
    flush=True,
)
q80_specialization_rows = []
q80_cm24_root_basis = q80_cm24_roots[3]
q80_cm24_root_gram_inverse = (
    q80_cm24_root_basis
    * q80_cm24
    * q80_cm24_root_basis.transpose()
).inverse()
for generic_lift in q80_generic_optimal_lifts.rows():
    cm_lift = generic_lift * q80_embedding24
    root_coordinates = (
        cm_lift
        * q80_cm24
        * q80_cm24_root_basis.transpose()
        * q80_cm24_root_gram_inverse
    )
    cm_projection = cm_lift - root_coordinates * q80_cm24_root_basis
    pairings = cm_projection * q80_cm24 * q80_cm24_lifts.transpose()
    coordinates = pairings * q80_cm24_height.inverse()
    assert all(value in ZZ for value in coordinates)
    q80_specialization_rows.append(list(map(ZZ, coordinates)))
q80_specialization = matrix(ZZ, q80_specialization_rows)
assert abs(q80_specialization.det()) == 1
q80_specialized_profiles = exact_section_profiles(
    q80_cm24,
    q80_cm24_height,
    q80_cm24_lifts,
    q80_specialization,
)

# Specialize the already-reduced first q4 class.  The unique root orthogonal
# to the transported generic A3+D5+E6 root lattice is the new CM24 A1.  Its
# effective sign is fixed by the transported generic P1, whose A1 label is
# one.  This distinguishes the live generic first-neighbor class from compact
# degree-two pencils that exist only after CM specialization.
q80_first_q4_reduced = vector(
    ZZ, (-6, -8, -11, -16, -13, -10, -13, -20, -26, -39, -8, -8, 20, 8, -2, -4, 2)
)
q80_first_q4_cm24 = q80_first_q4_reduced * q80_embedding24
q80_generic_roots_cm24 = q80_roots[3] * q80_embedding24
q80_cm24_half_roots = QuadraticForm(ZZ, q80_cm24).short_vector_list_up_to_length(
    2, up_to_sign_flag=True
)[1]
q80_cm24_all_roots = [vector(ZZ, row) for row in q80_cm24_half_roots]
q80_cm24_all_roots += [-row for row in q80_cm24_all_roots]
q80_extra_a1_candidates = [
    row for row in q80_cm24_all_roots
    if row * q80_cm24 * q80_generic_roots_cm24.transpose()
       == vector(ZZ, [0] * q80_generic_roots_cm24.nrows())
]
assert len(q80_extra_a1_candidates) == 2
q80_generic_p1_cm24 = q80_generic_optimal_lifts[0] * q80_embedding24
q80_effective_extra_a1 = next(
    row for row in q80_extra_a1_candidates
    if -q80_generic_p1_cm24 * q80_cm24 * row == 1
)
q80_first_q4_extra_pairing = -q80_first_q4_cm24 * q80_cm24 * q80_effective_extra_a1
print(
    f"KUMARCMQ80FIRSTQ4|cm24_extra_A1={tuple(q80_effective_extra_a1)}|"
    f"P1.A1=1|D.A1={q80_first_q4_extra_pairing}",
    flush=True,
)

q80_specialization43_rows = []
q80_cm43_root_basis = q80_cm43_roots[3]
q80_cm43_root_gram_inverse = (
    q80_cm43_root_basis
    * q80_cm43
    * q80_cm43_root_basis.transpose()
).inverse()
for generic_lift in q80_generic_optimal_lifts.rows():
    cm_lift = generic_lift * q80_embedding43
    root_coordinates = (
        cm_lift
        * q80_cm43
        * q80_cm43_root_basis.transpose()
        * q80_cm43_root_gram_inverse
    )
    cm_projection = cm_lift - root_coordinates * q80_cm43_root_basis
    pairings = cm_projection * q80_cm43 * q80_cm43_lifts.transpose()
    coordinates = pairings * q80_cm43_height.inverse()
    assert all(value in ZZ for value in coordinates)
    q80_specialization43_rows.append(list(map(ZZ, coordinates)))
q80_specialization43 = matrix(ZZ, q80_specialization43_rows)
q80_specialized43_profiles = exact_section_profiles(
    q80_cm43,
    q80_cm43_height,
    q80_cm43_lifts,
    q80_specialization43,
)
q80_generic_to_cm_polynomial = (
    q80_specialization * matrix(ZZ, q80_cm24_optimal[2]).inverse()
)
q80_cm_polynomial_to_generic = q80_generic_to_cm_polynomial.inverse()
assert all(value in ZZ for value in q80_cm_polynomial_to_generic.list())
q80_cm_basis_generic_profiles = exact_section_profiles(
    q80_child,
    q80_profiles[1],
    q80_generic_optimal_lifts,
    q80_cm_polynomial_to_generic,
)
print(
    "KUMARCMQ80|q=80|ab=4,20|generic_components={}|generic_MW={}|"
    "generic_optimal_poles={}|cm24_components={}|cm24_MW={}|"
    "cm24_optimal_poles={}|cm43_components={}|cm43_MW={}|"
    "cm43_optimal_poles={}|root_jump=1|all_MW_directions_survive=1".format(
        q80_components,
        matrix_rows(q80_height),
        q80_optimal,
        q80_cm24_components,
        matrix_rows(q80_cm24_height),
        q80_cm24_optimal,
        q80_cm43_components,
        matrix_rows(q80_cm43_height),
        q80_cm43_optimal,
    ),
    flush=True,
)
for stage, data in (
    ("generic", q80_profiles),
    ("cm24", q80_cm24_profiles),
    ("cm43", q80_cm43_profiles),
    ("generic_at_cm24", q80_specialized_profiles),
    ("generic_at_cm43", q80_specialized43_profiles),
    ("cm24_basis_generic", q80_cm_basis_generic_profiles),
):
    names, profile_height, profiles, poles, pair_intersections = data
    print(
        "KUMARCMQ80PROFILE|stage={}|components={}|height={}|profiles={}|"
        "P.O={}|pairs={}".format(
            stage,
            ",".join(names),
            matrix_rows(profile_height),
            ";".join(
                f"P{index + 1}:" + ",".join(map(str, profile))
                for index, profile in enumerate(profiles)
            ),
            ",".join(map(str, poles)),
            ";".join(
                f"P{left}.P{right}={intersection}"
                for left, right, intersection in pair_intersections
            ),
        ),
        flush=True,
    )
print(
    "KUMARCMQ80SPECIALIZATION|generic_optimal_to_cm_raw={}|"
    "generic_optimal_to_cm_polynomial={}|generic_optimal_to_cm43={}".format(
        matrix_rows(q80_specialization),
        matrix_rows(q80_generic_to_cm_polynomial),
        matrix_rows(q80_specialization43),
    ),
    flush=True,
)

# A q<60 neighbor cannot use the level-79 section: its height alone is
# 237/2.  For 60<=q<237 its coefficient is necessarily +/-1.  Weyl reduction
# in E7+E8 and the independent sign of the height-4 section therefore give a
# complete, very small orbit enumeration without traversing tens of thousands
# of root-lattice vectors.
orbit_representatives = level79_weyl_orbits(args.search_max_q)
assert len(orbit_representatives.get(60, ())) == 1
for orbit_q in sorted(orbit_representatives):
    print(
        f"KUMARCMORBIT|q={orbit_q}|weyl_orbits={len(orbit_representatives[orbit_q])}|"
        f"proper_factor_presentations={len(proper_factor_pairs(orbit_q))}",
        flush=True,
    )

if args.print_q80_frame:
    print("KUMARCMQ80FRAME|stage=generic", flush=True)
    print(q80_child, flush=True)
    print(
        f"KUMARCMQ80FRAME|generic_root_basis={matrix_rows(q80_roots[3])}",
        flush=True,
    )
    print(
        f"KUMARCMQ80FRAME|generic_reduced_MW_lifts={matrix_rows(q80_lifts)}|"
        f"optimal_change={matrix_rows(matrix(ZZ, q80_optimal[2]))}|"
        f"generic_optimal_MW_lifts={matrix_rows(q80_generic_optimal_lifts)}",
        flush=True,
    )
    print(
        f"KUMARCMQ80FRAME|generic_to_cm24_frame={matrix_rows(q80_embedding24)}",
        flush=True,
    )
    print("KUMARCMQ80FRAME|stage=cm24", flush=True)
    print(q80_cm24, flush=True)

if args.print_cm43_frame:
    print("KUMARCM43FRAME|stage=closure", flush=True)
    print(cm43, flush=True)
    print("KUMARCM43FRAME|stage=end", flush=True)

if args.print_marked_cm43_frame:
    marked_frame, marked_q60_fiber = marked_glue_frames[211]
    print("KUMARCM43MARKEDFRAME|glue=211|stage=closure", flush=True)
    print(marked_frame, flush=True)
    print(
        f"KUMARCM43MARKEDFRAME|q60_fiber={tuple(marked_q60_fiber)}",
        flush=True,
    )
    print("KUMARCM43MARKEDFRAME|stage=end", flush=True)

if not args.orbit_counts_only and args.search_max_q > 60:
    tested = 0
    stable = []
    jump_histogram = {}
    jump_one_types = {}
    best_jump = None
    for orbit_q in sorted(orbit_representatives):
        if orbit_q == 60:
            continue
        pairs = proper_factor_pairs(orbit_q)
        if not pairs:
            continue
        for orbit_index, witness in enumerate(orbit_representatives[orbit_q], 1):
            assert witness * frame * witness == 2 * orbit_q
            for factor_a, factor_b in pairs:
                _, basis_i, child_i = build_neighbor_basis(
                    frame, witness, factor_a, factor_b
                )
                generic_roots_i = root_invariants(child_i)
                cm24_i, _ = transport_neighbor_to_cm(
                    basis_i, closure_basis24, cm24
                )
                cm24_roots_i = root_invariants(cm24_i)
                root_jump = cm24_roots_i[0] - generic_roots_i[0]
                assert root_jump >= 0
                tested += 1
                jump_histogram[root_jump] = jump_histogram.get(root_jump, 0) + 1
                if best_jump is None or root_jump < best_jump:
                    best_jump = root_jump
                    print(
                        f"KUMARCMSEARCH|stage=new_best|q={orbit_q}|orbit={orbit_index}|"
                        f"ab={factor_a},{factor_b}|generic_root_rank={generic_roots_i[0]}|"
                        f"cm24_root_rank={cm24_roots_i[0]}|root_jump={root_jump}",
                        flush=True,
                    )
                if root_jump == 1:
                    type_key = tuple(generic_roots_i[:3] + cm24_roots_i[:3])
                    jump_one_types.setdefault(
                        type_key,
                        (
                            orbit_q,
                            orbit_index,
                            factor_a,
                            factor_b,
                            tuple(witness),
                            child_i,
                            cm24_i,
                            generic_roots_i[3],
                            cm24_roots_i[3],
                        ),
                    )
                if root_jump == 0:
                    generic_components_i = root_components(child_i)
                    cm24_components_i = root_components(cm24_i)
                    stable.append(
                        (
                            orbit_q,
                            orbit_index,
                            factor_a,
                            factor_b,
                            tuple(witness),
                            generic_components_i,
                            cm24_components_i,
                        )
                    )
                    print(
                        f"KUMARCMSEARCH|stage=stable_hit|q={orbit_q}|orbit={orbit_index}|"
                        f"ab={factor_a},{factor_b}|witness={tuple(witness)}|"
                        f"generic_components={generic_components_i}|"
                        f"cm24_components={cm24_components_i}",
                        flush=True,
                    )
    print(
        f"KUMARCMSEARCH|stage=summary|max_q={args.search_max_q}|tested={tested}|"
        f"jump_histogram={sorted(jump_histogram.items())}|stable_hits={len(stable)}|"
        f"best_jump={best_jump}",
        flush=True,
    )
    for type_index, (type_key, record) in enumerate(
        sorted(jump_one_types.items()), 1
    ):
        (
            orbit_q,
            orbit_index,
            factor_a,
            factor_b,
            witness,
            child_i,
            cm24_i,
            generic_root_basis_i,
            cm24_root_basis_i,
        ) = record
        generic_components_i = root_components(child_i)
        cm24_components_i = root_components(cm24_i)
        _, generic_height_i, generic_lifts_i = mw_height_gram(
            child_i, generic_root_basis_i, return_lifts=True
        )
        _, cm24_height_i, cm24_lifts_i = mw_height_gram(
            cm24_i, cm24_root_basis_i, return_lifts=True
        )
        generic_poles_i = section_pole_data(
            child_i, generic_height_i, generic_lifts_i
        )
        cm24_poles_i = section_pole_data(
            cm24_i, cm24_height_i, cm24_lifts_i
        )
        generic_optimal_poles_i = optimal_section_pole_basis(
            child_i, generic_height_i, generic_lifts_i
        )
        cm24_optimal_poles_i = optimal_section_pole_basis(
            cm24_i, cm24_height_i, cm24_lifts_i
        )
        print(
            f"KUMARCMSEARCH|stage=jump_one_type|type={type_index}|"
            f"root_invariants={type_key}|q={orbit_q}|orbit={orbit_index}|"
            f"ab={factor_a},{factor_b}|witness={witness}|"
            f"generic_components={generic_components_i}|"
            f"generic_mw={matrix_rows(generic_height_i)}|"
            f"generic_section_data={generic_poles_i}|"
            f"generic_optimal_poles={generic_optimal_poles_i}|"
            f"cm24_components={cm24_components_i}|"
            f"cm24_mw={matrix_rows(cm24_height_i)}|"
            f"cm24_section_data={cm24_poles_i}|"
            f"cm24_optimal_poles={cm24_optimal_poles_i}",
            flush=True,
        )

print("KUMARCM|status=PASS", flush=True)
