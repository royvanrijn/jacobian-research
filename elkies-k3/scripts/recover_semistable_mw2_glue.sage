from sage.all import *
from pathlib import Path


CASES = (
    {
        "name": "A6+A4+A3+A2",
        "frame": Path("elkies-k3/data/fibrations/mw2_a6_a4_a3_a2_frame.txt"),
        "ranks": (6, 4, 3, 2),
        "denominator": 105,
        "target_diagonal": (81, 326),
        "height_det": QQ(79) / 35,
        "optimized_change": ((1, 0), (0, 1)),
        "expected_optimized": (
            ((81, 39), (39, 326)),
            ((5, 4, 2, 0), (5, 1, 0, 2)),
            (0, 1), 1, (1, 1, 1),
        ),
    },
    {
        "name": "A5+A4+2A3",
        "frame": Path("elkies-k3/data/fibrations/mw2_a5_a4_a3a3_frame.txt"),
        "ranks": (5, 4, 3, 3),
        "denominator": 60,
        "target_diagonal": (67, 117),
        "height_det": QQ(79) / 40,
        "optimized_change": ((1, 0), (1, -1)),
        "expected_optimized": (
            ((67, 40), (40, 130)),
            ((2, 4, 0, 1), (5, 0, 2, 0)),
            (0, 0), 1, (0, 0, 1),
        ),
    },
    {
        "name": "A9+A5+A1",
        "frame": Path("elkies-k3/data/fibrations/mw2_a9_a5_a1_frame.txt"),
        "ranks": (9, 5, 1),
        "denominator": 30,
        "target_diagonal": (117, 248),
        "height_det": QQ(158) / 5,
        "optimized_change": ((1, 0), (0, 1)),
        "expected_optimized": (
            ((117, 24), (24, 248)),
            ((8, 0, 1), (1, 5, 0)),
            (1, 3), 5, (3, 4, 5),
        ),
    },
)


def qform_from_gram(gram):
    coefficients = []
    for i in range(gram.nrows()):
        for j in range(i, gram.ncols()):
            coefficients.append(gram[i, i] // 2 if i == j else gram[i, j])
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


def fractional_class(coordinates):
    return tuple(QQ(value) - floor(QQ(value)) for value in coordinates)


def add_classes(left, right):
    return fractional_class(tuple(a + b for a, b in zip(left, right)))


def class_multiple(multiplier, point):
    result = tuple(QQ(0) for _ in point)
    for _ in range(multiplier):
        result = add_classes(result, point)
    return result


def class_order(point, exponent):
    zero = tuple(QQ(0) for _ in point)
    return next(
        order for order in range(1, exponent + 1)
        if class_multiple(order, point) == zero
    )


def mod_two(value):
    value = QQ(value)
    return value - 2 * floor(value / 2)


def discriminant_generator(root_gram):
    """Return a quotient generator and its standard A_n component label."""
    modulus = abs(ZZ(root_gram.det()))
    inverse = root_gram.inverse()
    generators = [fractional_class(inverse.row(i)) for i in range(inverse.nrows())]
    zero = tuple(QQ(0) for _ in range(inverse.nrows()))
    classes = {zero}
    queue = [zero]
    head = 0
    while head < len(queue):
        point = queue[head]
        head += 1
        for generator in generators:
            candidate = add_classes(point, generator)
            if candidate not in classes:
                classes.add(candidate)
                queue.append(candidate)
    assert len(classes) == modulus
    full_order = sorted(
        point for point in classes if class_order(point, modulus) == modulus
    )
    assert full_order
    generator = full_order[0]

    discriminant_norm = mod_two(
        vector(QQ, generator) * root_gram * vector(QQ, generator)
    )
    possible_labels = sorted({
        min(label, modulus - label)
        for label in range(1, modulus)
        if gcd(label, modulus) == 1
        and mod_two(QQ(label * (modulus - label)) / modulus) == discriminant_norm
    })
    assert len(possible_labels) == 1
    return generator, possible_labels[0]


def component_label(point, generator, generator_label, modulus):
    multiplier = next(
        value for value in range(modulus)
        if class_multiple(value, generator) == point
    )
    return (multiplier * generator_label) % modulus


def a_correction(rank, left, right):
    if left == 0 or right == 0:
        return QQ(0)
    order = rank + 1
    return QQ(min(left, right) * (order - max(left, right))) / order


def local_correction(ranks, left, right):
    return sum(
        a_correction(rank, left[index], right[index])
        for index, rank in enumerate(ranks)
    )


def transform_profiles(ranks, profiles, transform):
    return tuple(
        tuple(
            sum(transform[row, column] * profiles[column][component]
                for column in range(2)) % (rank + 1)
            for component, rank in enumerate(ranks)
        )
        for row in range(2)
    )


def shioda_data(ranks, gram, profiles):
    zero_intersections = tuple(
        (gram[row, row] + local_correction(ranks, profiles[row], profiles[row]) - 4) / 2
        for row in range(2)
    )
    pair_intersection = (
        2 + sum(zero_intersections)
        - local_correction(ranks, profiles[0], profiles[1])
        - gram[0, 1]
    )
    assert all(value in ZZ and value >= 0 for value in zero_intersections)
    assert pair_intersection in ZZ and pair_intersection >= 0
    return tuple(ZZ(value) for value in zero_intersections), ZZ(pair_intersection)


def find_target_transforms(scaled, target):
    bound = max(target.diagonal())
    minimum = pari(scaled).qfminim(bound)
    representatives = list(matrix(ZZ, minimum[2]).columns())
    vectors = representatives + [-point for point in representatives]
    first = [point for point in vectors if point * scaled * point == target[0, 0]]
    second = [point for point in vectors if point * scaled * point == target[1, 1]]
    transforms = []
    for left in first:
        for right in second:
            transform = matrix(ZZ, [left, right])
            if abs(transform.det()) != 1:
                continue
            if transform * scaled * transform.transpose() == target:
                transforms.append(transform)
    unique = {tuple(transform.list()): transform for transform in transforms}
    return [unique[key] for key in sorted(unique)]


def recover_case(case):
    F = matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in case["frame"].read_text().splitlines()
        if line.strip()
    ])
    assert F.nrows() == 17 and F.det() == 948 and F.is_positive_definite()

    half_roots = [
        vector(ZZ, root)
        for root in qform_from_gram(F).short_vector_list_up_to_length(2, True)[1]
    ]
    roots = half_roots + [-root for root in half_roots]
    graph = Graph()
    graph.add_vertices(range(len(roots)))
    for i in range(len(roots)):
        for j in range(i):
            if roots[i] * F * roots[j] != 0:
                graph.add_edge(i, j)
    components = sorted(graph.connected_components(sort=False), key=len, reverse=True)
    component_bases = [
        matrix(ZZ, [roots[index] for index in component]).row_module().basis_matrix()
        for component in components
    ]
    component_ranks = tuple(basis.rank() for basis in component_bases)
    assert component_ranks == case["ranks"]
    component_grams = [basis * F * basis.transpose() for basis in component_bases]
    assert tuple(abs(gram.det()) for gram in component_grams) == tuple(
        rank + 1 for rank in component_ranks
    )

    R = block_matrix([[basis] for basis in component_bases], subdivide=False)
    GR = R * F * R.transpose()
    C = (R * F).right_kernel_matrix()
    GC = C * F * C.transpose()
    assert R.rank() == 15 and C.rank() == 2

    A = block_matrix([[R], [C]], subdivide=False)
    index = abs(A.det())
    A_inverse = A.inverse()

    def coset_key(point):
        return fractional_class(vector(QQ, point) * A_inverse)

    zero = vector(ZZ, [0] * 17)
    cosets = {coset_key(zero): zero}
    queue = [zero]
    head = 0
    while head < len(queue) and len(cosets) < index:
        point = queue[head]
        head += 1
        for coordinate in range(17):
            unit = vector(ZZ, [0] * 17)
            unit[coordinate] = 1
            for sign in (1, -1):
                candidate = point + sign * unit
                key = coset_key(candidate)
                if key not in cosets:
                    cosets[key] = candidate
                    queue.append(candidate)
                    if len(cosets) == index:
                        break
            if len(cosets) == index:
                break
    assert len(cosets) == index

    GR_inverse = GR.inverse()
    GC_inverse = GC.inverse()

    def project_mw(point):
        point = vector(QQ, point)
        return point - (point * F * R.transpose()) * GR_inverse * R

    def coordinates_in_C(point):
        return (vector(QQ, point) * F * C.transpose()) * GC_inverse

    projected_cosets = [
        (coordinates_in_C(project_mw(representative)), representative)
        for representative in cosets.values()
    ]
    generators = [vector(QQ, (1, 0)), vector(QQ, (0, 1))]
    generators += [coordinates for coordinates, _ in projected_cosets]
    denominator = lcm(
        QQ(value).denominator() for row in generators for value in row
    )
    MW_integer = matrix(ZZ, [
        [ZZ(denominator * value) for value in row]
        for row in generators
    ]).row_module().basis_matrix()
    MW_basis = MW_integer.change_ring(QQ) / denominator
    H = MW_basis * GC * MW_basis.transpose()
    assert H.det() == case["height_det"]

    component_bounds = []
    start = 0
    for basis in component_bases:
        component_bounds.append((start, start + basis.nrows()))
        start += basis.nrows()
    discriminant_generators = [discriminant_generator(gram) for gram in component_grams]

    def profiles_for_transform(transform):
        target_basis = transform * MW_basis
        target = target_basis * GC * target_basis.transpose()
        root_classes = []
        for target_vector in target_basis.rows():
            lift = None
            for projected, representative in projected_cosets:
                difference = target_vector - projected
                if all(QQ(value).denominator() == 1 for value in difference):
                    lift = vector(QQ, representative) + difference * C
                    break
            assert lift is not None
            assert all(QQ(value).denominator() == 1 for value in lift)
            lift = vector(ZZ, lift)
            assert project_mw(lift) == target_vector * C

            root_part = vector(QQ, lift) - target_vector * C
            root_coordinates = (root_part * F * R.transpose()) * GR_inverse
            assert root_part == root_coordinates * R
            root_classes.append([
                fractional_class(root_coordinates[left:right])
                for left, right in component_bounds
            ])

        profiles = []
        for row in range(2):
            profile = []
            for component, rank in enumerate(component_ranks):
                generator, generator_label = discriminant_generators[component]
                profile.append(component_label(
                    root_classes[row][component], generator,
                    generator_label, rank + 1,
                ))
            profiles.append(tuple(profile))

        # Give repeated A factors a deterministic order.
        repeated = {}
        for component, rank in enumerate(component_ranks):
            repeated.setdefault(rank, []).append(component)
        order = list(range(len(component_ranks)))
        for indices in repeated.values():
            if len(indices) > 1:
                sorted_indices = sorted(indices, key=lambda i: (profiles[0][i], profiles[1][i]))
                for destination, source in zip(indices, sorted_indices):
                    order[destination] = source
        profiles = [tuple(profile[index] for index in order) for profile in profiles]
        ordered_ranks = tuple(component_ranks[index] for index in order)

        zero_intersections, pair_intersection = shioda_data(
            ordered_ranks, target, profiles
        )
        return tuple(profiles), zero_intersections, pair_intersection

    denominator = case["denominator"]
    scaled = (denominator * H).change_ring(ZZ)
    first, second = case["target_diagonal"]
    orientations = []
    for cross_sign in (-1, 1):
        # The absolute reduced off-diagonal is determined by the determinant.
        cross_square = first * second - ZZ(case["height_det"] * denominator^2)
        cross = cross_sign * ZZ(sqrt(cross_square))
        assert cross * cross == cross_square
        target_scaled = matrix(ZZ, [[first, cross], [cross, second]])
        transforms = find_target_transforms(scaled, target_scaled)
        assert transforms
        profiles, zero_intersections, pair_intersection = profiles_for_transform(transforms[0])
        orientations.append((
            target_scaled, profiles, zero_intersections, pair_intersection,
            len(transforms), index,
        ))

    positive = next(item for item in orientations if item[0][0, 1] > 0)
    base_scaled, base_profiles = positive[0], positive[1]
    optimized_change = matrix(ZZ, case["optimized_change"])
    optimized_scaled = optimized_change * base_scaled * optimized_change.transpose()
    optimized_profiles = transform_profiles(
        component_ranks, base_profiles, optimized_change
    )
    optimized_gram = optimized_scaled / denominator
    optimized_O, optimized_pair = shioda_data(
        component_ranks, optimized_gram, optimized_profiles
    )

    # Exhaust all MW vectors whose pole count is at most that of the selected
    # basis.  Since h(P)<=4+2(P.O), qfminim at this bound contains every
    # competing vector.  Check all unimodular pairs and certify the optimal
    # (maximum pole, total poles, pair-intersection) triple.
    max_pole = max(optimized_O)
    norm_bound = denominator * (4 + 2 * max_pole)
    minimum = pari(optimized_scaled).qfminim(norm_bound)
    representatives = list(matrix(ZZ, minimum[2]).columns())
    vectors = representatives + [-point for point in representatives]
    section_vectors = []
    for coefficients in vectors:
        profile = tuple(
            sum(coefficients[row] * optimized_profiles[row][component]
                for row in range(2)) % (rank + 1)
            for component, rank in enumerate(component_ranks)
        )
        height = QQ(coefficients * optimized_scaled * coefficients) / denominator
        pole = (height + local_correction(component_ranks, profile, profile) - 4) / 2
        assert pole in ZZ and pole >= 0
        if pole <= max_pole:
            section_vectors.append((coefficients, profile, ZZ(pole)))

    basis_scores = []
    for left, left_profile, left_pole in section_vectors:
        for right, right_profile, right_pole in section_vectors:
            if abs(matrix(ZZ, [left, right]).det()) != 1:
                continue
            height_pair = QQ(left * optimized_scaled * right) / denominator
            section_pair = (
                2 + left_pole + right_pole
                - local_correction(component_ranks, left_profile, right_profile)
                - height_pair
            )
            assert section_pair in ZZ and section_pair >= 0
            basis_scores.append((
                max(left_pole, right_pole), left_pole + right_pole,
                ZZ(section_pair),
            ))
    assert basis_scores
    optimal_score = min(basis_scores)

    expected_scaled, expected_profiles, expected_O, expected_pair, expected_score = case["expected_optimized"]
    assert optimized_scaled == matrix(ZZ, expected_scaled)
    assert optimized_profiles == expected_profiles
    assert optimized_O == expected_O
    assert optimized_pair == expected_pair
    assert optimal_score == expected_score

    optimized = (
        optimized_scaled, optimized_profiles, optimized_O, optimized_pair,
        optimal_score, len(section_vectors), norm_bound,
    )
    return component_ranks, len(roots), abs(GR.det()), orientations, optimized


for case in CASES:
    ranks, root_count, root_det, orientations, optimized = recover_case(case)
    print(
        f"SEMI_MW2|case={case['name']}|roots={root_count}|rootdet={root_det}"
        f"|ranks={','.join(map(str, ranks))}",
        flush=True,
    )
    for target, profiles, zero_intersections, pair_intersection, transforms, index in orientations:
        print(
            f"SEMI_MW2|case={case['name']}|scaled_gram="
            f"{target[0,0]},{target[0,1]};{target[1,0]},{target[1,1]}"
            f"|profiles=P1:{','.join(map(str, profiles[0]))};P2:{','.join(map(str, profiles[1]))}"
            f"|O={zero_intersections[0]},{zero_intersections[1]}"
            f"|pair={pair_intersection}|transforms={transforms}|glue_index={index}",
            flush=True,
        )
    target, profiles, zero_intersections, pair_intersection, score, vector_count, norm_bound = optimized
    print(
        f"SEMI_MW2|case={case['name']}|optimized_scaled_gram="
        f"{target[0,0]},{target[0,1]};{target[1,0]},{target[1,1]}"
        f"|profiles=P1:{','.join(map(str, profiles[0]))};P2:{','.join(map(str, profiles[1]))}"
        f"|O={zero_intersections[0]},{zero_intersections[1]}|pair={pair_intersection}"
        f"|optimal_score={score}|vectors_checked={vector_count}|norm_bound={norm_bound}",
        flush=True,
    )
    print(f"SEMI_MW2|case={case['name']}|PASS", flush=True)
