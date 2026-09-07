#!/usr/bin/env sage
"""Certify the two generic q=8 neighbors of the A13+A1 MW3 frame.

The script reduces each raw isotropic class against the old zero section and
all old-fiber components, then proves nefness. Section walls are reduced to a
finite rank-three MW closest-vector calculation; root primitivity and the
standard norm identity exclude a negative degree-two multisection. The child
root/MW-frame isometries are certified separately by analyze_mw3_branch.sage.
"""

from pathlib import Path

from sage.all import (
    QQ,
    ZZ,
    Matrix,
    QuadraticForm,
    block_diagonal_matrix,
    gcd,
    lcm,
    matrix,
    vector,
)
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
FRAME_PATH = ROOT / "elkies-k3/data/fibrations/mw3_a13_a1_frame.txt"
U = matrix(ZZ, [[0, 1], [1, 0]])
WITNESSES = (
    (0, -2, -3, -2, 2, 2, -4, 5, -4, 5, 2, 0, -2, 0, 0, 0, 0),
    (0, 0, -1, 3, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0),
)


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def qform(gram):
    coefficients = []
    for row in range(gram.nrows()):
        for column in range(row, gram.ncols()):
            coefficients.append(
                gram[row, row] // 2 if row == column else gram[row, column]
            )
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


def roots_of_norm_two(gram):
    half = qform(gram).short_vector_list_up_to_length(
        2, up_to_sign_flag=True
    )[1]
    half = [vector(ZZ, row) for row in half]
    return half + [-row for row in half]


def lex_positive(row):
    return next(value > 0 for value in row if value)


def deterministic_simple_roots(gram):
    roots = roots_of_norm_two(gram)
    positive = [row for row in roots if lex_positive(row)]
    positive_set = {tuple(row) for row in positive}
    simple = []
    for row in positive:
        if not any(
            tuple(row-left) in positive_set for left in positive if left != row
        ):
            simple.append(row)
    simple = matrix(ZZ, [list(row) for row in simple])
    assert simple.nrows() == simple.rank() == 14
    return simple, positive


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        todo = [min(unseen)]
        unseen.remove(todo[0])
        component = []
        while todo:
            index = todo.pop()
            component.append(index)
            adjacent = [other for other in unseen if cartan[index, other]]
            for other in adjacent:
                unseen.remove(other)
                todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(sorted(result, key=lambda component: (len(component), component)))


def highest_roots(gram, simple, positive):
    cartan = simple * gram * simple.transpose()
    inverse_simple = simple.pseudoinverse()
    result = []
    for component in connected_components(cartan):
        candidates = []
        for root in positive:
            coordinates = vector(QQ, root) * inverse_simple
            if not all(value in ZZ and value >= 0 for value in coordinates):
                continue
            support = tuple(index for index, value in enumerate(coordinates) if value)
            if support and all(index in component for index in support):
                candidates.append((sum(coordinates), vector(ZZ, root)))
        result.append(max(candidates, key=lambda item: item[0])[1])
    return tuple(result)


def intersection(left, right, ns):
    return ZZ(vector(ZZ, left) * ns * vector(ZZ, right))


def chamber_reduce(divisor, curves, ns):
    divisor = vector(ZZ, divisor)
    sequence = []
    while True:
        for name, curve in curves:
            pairing = intersection(divisor, curve, ns)
            if pairing < 0:
                divisor += pairing * curve
                sequence.append((name, pairing))
                assert divisor * ns * divisor == 0
                break
        else:
            return divisor, tuple(sequence)


def mod_one(value):
    value = QQ(value)
    return value - value.floor()


frame = load_matrix(FRAME_PATH)
ns = block_diagonal_matrix(U, -frame)
simple, positive = deterministic_simple_roots(frame)
cartan = simple * frame * simple.transpose()
components = connected_components(cartan)
assert sorted(map(len, components)) == [1, 13]

fiber = vector(ZZ, [1, 0] + [0] * 17)
zero = vector(ZZ, [-1, 1] + [0] * 17)
curves = [("O", zero)]
for index, root in enumerate(simple.rows(), start=1):
    curves.append((f"R{index}", vector(ZZ, [0, 0] + list(root))))
for index, root in enumerate(highest_roots(frame, simple, positive), start=1):
    theta_zero = fiber - vector(ZZ, [0, 0] + list(root))
    assert theta_zero * ns * theta_zero == -2
    curves.append((f"Theta0_{index}", theta_zero))

# Saturated rank-three MW quotient and integral representatives of all glue
# cosets. The A13+A1 root lattice is primitive, so the old MW torsion is one.
root_basis = matrix(ZZ, [list(root) for root in roots_of_norm_two(frame)]).row_module().basis_matrix()
root_gram = root_basis * frame * root_basis.transpose()
root_smith = root_basis.smith_form()[0]
assert tuple(abs(root_smith[index, index]) for index in range(14)) == (1,) * 14
essential_basis = (root_basis * frame).right_kernel_matrix()
essential_gram = essential_basis * frame * essential_basis.transpose()
combined_basis = matrix(ZZ, list(root_basis.rows()) + list(essential_basis.rows()))
combined_inverse = combined_basis.inverse()
glue_index = abs(combined_basis.det())


def glue_key(row):
    coordinates = vector(QQ, row) * combined_inverse
    return tuple(mod_one(value) for value in coordinates)


zero_frame = vector(ZZ, [0] * 17)
cosets = {glue_key(zero_frame): zero_frame}
queue = [zero_frame]
head = 0
while head < len(queue) and len(cosets) < glue_index:
    current = queue[head]
    head += 1
    for index in range(17):
        basis_vector = vector(ZZ, [1 if position == index else 0 for position in range(17)])
        for sign in (1, -1):
            candidate = current + sign*basis_vector
            key = glue_key(candidate)
            if key not in cosets:
                cosets[key] = candidate
                queue.append(candidate)
                if len(cosets) == glue_index:
                    break
        if len(cosets) == glue_index:
            break
assert len(cosets) == glue_index

root_gram_inverse = root_gram.inverse()
essential_gram_inverse = essential_gram.inverse()


def root_orthogonal_projection(row):
    row = vector(QQ, row)
    coefficients = row * frame * root_basis.transpose() * root_gram_inverse
    return row - coefficients * root_basis


def essential_coordinates(row):
    return vector(QQ, row) * frame * essential_basis.transpose() * essential_gram_inverse


projected_cosets = tuple(
    (representative, essential_coordinates(root_orthogonal_projection(representative)))
    for representative in cosets.values()
)
coordinate_rows = [
    vector(QQ, [1 if row == column else 0 for column in range(3)])
    for row in range(3)
] + [coordinates for _, coordinates in projected_cosets]
denominator = lcm(value.denominator() for row in coordinate_rows for value in row)
mw_basis_scaled = matrix(
    ZZ,
    [[ZZ(denominator*value) for value in row] for row in coordinate_rows],
).row_module().basis_matrix()
mw_basis = mw_basis_scaled.change_ring(QQ) / denominator
height = mw_basis * essential_gram * mw_basis.transpose()
assert height.det() == QQ(237) / 7
root_lattice = IntegralLattice(root_gram)


def fractional_key(coordinates):
    return tuple(mod_one(value) for value in coordinates)


coset_by_projection = {
    fractional_key(coordinates): (representative, coordinates)
    for representative, coordinates in projected_cosets
}


def shortest_section_lift(mw_coordinates):
    target_coordinates = vector(QQ, mw_coordinates) * mw_basis
    representative, representative_coordinates = coset_by_projection[
        fractional_key(target_coordinates)
    ]
    difference = target_coordinates - representative_coordinates
    assert all(value in ZZ for value in difference)
    raw_lift = vector(ZZ, representative) + vector(ZZ, difference) * essential_basis
    root_coordinates = (
        raw_lift * frame * root_basis.transpose() * root_gram_inverse
    )
    iterator = root_lattice.enumerate_close_vectors(-root_coordinates)
    shift = vector(ZZ, next(iterator))
    lift = raw_lift + shift * root_basis
    assert ZZ(lift * frame * lift) % 2 == 0
    return vector(ZZ, lift)


def nearby_mw_points(center_twice):
    # n=2*m-center_twice and dist(m,center/2)^2<2 iff n.H.n<8.
    scale = lcm([value.denominator() for value in height.list()] + [2])
    if scale % 2:
        scale *= 2
    integral_height = (scale*height).change_ring(ZZ)
    if any(value % 2 for value in integral_height.diagonal()):
        integral_height *= 2
        scale *= 2
    form = QuadraticForm(ZZ, integral_height)
    shells = form.short_vector_list_up_to_length(8*scale, up_to_sign_flag=True)
    vectors = {tuple([0, 0, 0])}
    for shell in shells[1:]:
        for row in shell:
            row = vector(ZZ, row)
            vectors.add(tuple(row))
            vectors.add(tuple(-row))
    result = []
    for row_tuple in vectors:
        row = vector(ZZ, row_tuple)
        if not all((row[index]+center_twice[index]) % 2 == 0 for index in range(3)):
            continue
        mw = vector(ZZ, [(row[index]+center_twice[index]) // 2 for index in range(3)])
        difference = vector(QQ, mw) - vector(QQ, center_twice)/2
        distance = difference * height * difference
        if distance < 2:
            result.append((tuple(mw), distance))
    return tuple(sorted(set(result)))


print(
    f"A13Q8|old_components={tuple(map(len, components))}|"
    f"root_primitive=1|MW_torsion=1|height_gram={height}|height_det={height.det()}",
    flush=True,
)

for witness_index, witness in enumerate(WITNESSES, start=1):
    raw = vector(ZZ, [2, 4] + list(witness))
    assert raw * ns * raw == 0 and gcd(list(raw)) == 1
    reduced, sequence = chamber_reduce(raw, curves, ns)
    assert intersection(reduced, fiber, ns) == 2
    old_pairings = tuple(
        (name, intersection(reduced, curve, ns)) for name, curve in curves
    )
    assert all(pairing >= 0 for _, pairing in old_pairings)

    frame_part = vector(ZZ, reduced[2:])
    mw_projection = root_orthogonal_projection(frame_part)
    mw_coordinates = essential_coordinates(mw_projection) * mw_basis.inverse()
    assert all(value in ZZ for value in mw_coordinates)
    mw_coordinates = vector(ZZ, mw_coordinates)
    root_part = vector(QQ, frame_part) - mw_projection
    mw_norm = mw_projection * frame * mw_projection
    root_norm = root_part * frame * root_part
    fiber_coefficient = reduced[0] + reduced[1]
    assert reduced[1] == 2
    assert mw_norm + root_norm == 4*fiber_coefficient-8

    nearby = nearby_mw_points(mw_coordinates)
    section_tests = []
    for section_mw, distance in nearby:
        if section_mw == (0, 0, 0):
            continue  # O was already included in the chamber reduction.
        lift = shortest_section_lift(section_mw)
        norm = ZZ(lift * frame * lift)
        section_height = vector(QQ, section_mw) * height * vector(QQ, section_mw)
        correction = QQ(norm) - section_height
        assert correction >= 0

        # Refined root Cauchy bound:
        # D.S >= distance-2+(sqrt(c)-sqrt(c_D)/2)^2.
        # Avoid radicals by checking L=c+c_D/4-(2-distance) >= sqrt(c*c_D).
        linear = correction + root_norm/4 - (2-distance)
        cauchy_certified = linear >= 0 and linear**2 >= correction*root_norm
        assert cauchy_certified
        pole = norm // 2 - 2
        section_tests.append(
            (section_mw, distance, section_height, correction, pole)
        )

    # The universal root Cauchy bound is D.S >= dist(MW(S),MW(D)/2)^2-2,
    # so `nearby` contains every section that could pair negatively. If a
    # negative degree-two curve existed, D-C would be vertical and C would
    # have the same MW projection as D. The identity
    # ||r-rho||^2=2(D.C+1) would force r=rho and a half-integral fiber
    # coefficient, impossible. Thus components+sections suffice for nefness.
    print(
        f"A13Q8|neighbor={witness_index}|raw={tuple(raw)}|"
        f"reflection_count={len(sequence)}|reduction={sequence}|reduced={tuple(reduced)}|"
        f"D.F=2|D.O={intersection(reduced, zero, ns)}|"
        f"fiber_coefficient={fiber_coefficient}|MW={tuple(mw_coordinates)}|"
        f"MW_norm={mw_norm}|root_norm={root_norm}",
        flush=True,
    )
    print(
        f"A13Q8|neighbor={witness_index}|old_curve_pairings={old_pairings}|"
        f"nearby_MW={nearby}|section_cauchy_tests={tuple(section_tests)}|"
        "all_sections_nonnegative=1|negative_bisection_impossible=1|"
        "primitive_isotropic=1|nef=1|status=PASS_NEF",
        flush=True,
    )

print("A13Q8|neighbors=2|status=PASS", flush=True)
