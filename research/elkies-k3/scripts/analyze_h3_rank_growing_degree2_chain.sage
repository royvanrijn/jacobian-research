#!/usr/bin/env sage -python
"""Certify the degree-two H3 continuation from D12/MW5 to A3^3/MW8.

The root-adapted Weyl searches select the chain

    D12/MW5 --q6--> A11/MW6 --q8--> A5+A5/MW7
              --q4--> 3A3/MW8 --q4--> A3+2A2/MW10.

Every displayed isotropic class already lies in its source component chamber.
Full nefness is proved without a bounded section search.  Section
intersections are exact shifted closest-vector problems in the root lattice
and MW quotient; negative bisections are excluded by the universal
degree-two norm/parity identity.
"""

import itertools
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]

CASES = (
    {
        "name": "D12Q6",
        "frame": ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q24-degree2-frames/q24-o0085-r12-n264-d4-add367fba084.txt",
        "root_type": "D12",
        "root_rank": 12,
        "q": 6,
        "a": 3,
        "witness": (1, 4, -1, 0, 7, 4, 3, 7, 7, 7, 7, 7, -1, 0, -1, -1, 0),
        "component_pairings": (0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0),
        "affine_pairings": (0,),
        "section_distances": (3,) * 8,
        "child": "A11/MW6",
        "child_root_data": (11, 132, 12),
    },
    {
        "name": "A11Q8",
        "frame": ROOT / "artifacts/generated-results/elkies-k3-h3-d12-o85-q6-degree2-frames/q6-o0042-r11-n132-d12-e7e61e5dd4c2.txt",
        "root_type": "A11",
        "root_rank": 11,
        "q": 8,
        "a": 4,
        "witness": (0, -1, -5, 2, -2, -5, -4, -3, -3, -3, -3, -1, 0, 0, 1, 0, -1),
        "component_pairings": (0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0),
        "affine_pairings": (1,),
        "section_distances": (),
        "child": "A5+A5/MW7",
        "child_root_data": (10, 60, 36),
    },
    {
        "name": "A5A5Q4",
        "frame": ROOT / "artifacts/generated-results/elkies-k3-h3-a11-middle-q8-degree2-frames/q8-o0922-r10-n60-d36-c9cd5a498117.txt",
        "root_type": "A5+A5",
        "root_rank": 10,
        "q": 4,
        "a": 2,
        "witness": (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, -1, 0, 0, 0, 1),
        "component_pairings": (1, 1, 1, 0, 1, 0, 0, 0, 0, 0),
        "affine_pairings": (0, 0),
        "section_distances": (2, 2, 3, 3, 3, 3, 3, 3),
        "child": "3A3/MW8",
        "child_root_data": (9, 36, 64),
    },
    {
        "name": "3A3Q4",
        "frame": ROOT / "artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2-frames/q4-o0472-r9-n36-d64-4841c34fa442.txt",
        "root_type": "3A3",
        "root_rank": 9,
        "q": 4,
        "a": 2,
        "witness": (-3, 0, 2, 2, 3, 1, 2, 3, 1, -1, 0, 0, 1, 0, 1, 0, -1),
        "component_pairings": (1, 0, 1, 0, 0, 0, 1, 1, 0),
        "affine_pairings": (0, 1, 1),
        "section_distances": (2, 2) + (3,) * 16 + (4, 4),
        "child": "A3+2A2/MW10",
        "child_root_data": (7, 24, 36),
    },
)


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


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
            for other in tuple(unseen):
                if cartan[index, other]:
                    unseen.remove(other)
                    todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(result)


def highest_roots(cartan):
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-row for row in half)
    result = []
    for component in connected_components(cartan):
        candidates = [
            row for row in roots
            if all(value >= 0 for value in row)
            and all(index in component or row[index] == 0 for index in range(cartan.nrows()))
        ]
        result.append(max(candidates, key=lambda row: sum(row)))
    return tuple(result)


def bezout_vector_for_pairing(ns, fiber):
    current = ZZ(0)
    result = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fiber):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        result = [left * entry for entry in result]
        result[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        result = [-entry for entry in result]
    return vector(ZZ, result)


def neighbor_frame(ns, fiber):
    mate = bezout_vector_for_pairing(ns, fiber)
    mate -= (mate * ns * mate // 2) * fiber
    complement = matrix(
        ZZ, [list(fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    full_basis = matrix(
        ZZ, [list(fiber), list(mate)] + [list(row) for row in complement]
    )
    assert fiber * ns * mate == 1 and mate * ns * mate == 0
    assert abs(full_basis.det()) == 1
    return -(complement * ns * complement.transpose())


def standard_root_model(root_type):
    if root_type == "D12":
        cartan = CartanMatrix(["D", 12])
        ambient = matrix(ZZ, 12, 12)
        for index in range(11):
            ambient[index, index] = 1
            ambient[index, index + 1] = -1
        ambient[11, 10] = 1
        ambient[11, 11] = 1
        valid = lambda row: sum(row) % 2 == 0
    elif root_type == "A5+A5":
        a5 = CartanMatrix(["A", 5])
        cartan = block_diagonal_matrix(a5, a5)
        ambient = matrix(ZZ, 10, 12)
        for row_offset, column_offset in ((0, 0), (5, 6)):
            for index in range(5):
                ambient[row_offset + index, column_offset + index] = 1
                ambient[row_offset + index, column_offset + index + 1] = -1
        valid = lambda row: sum(row[:6]) == 0 and sum(row[6:]) == 0
    elif root_type == "3A3":
        a3 = CartanMatrix(["A", 3])
        cartan = block_diagonal_matrix(a3, a3, a3)
        ambient = matrix(ZZ, 9, 12)
        for row_offset, column_offset in ((0, 0), (3, 4), (6, 8)):
            for index in range(3):
                ambient[row_offset + index, column_offset + index] = 1
                ambient[row_offset + index, column_offset + index + 1] = -1
        valid = lambda row: all(sum(row[offset:offset + 4]) == 0 for offset in (0, 4, 8))
    else:
        raise ValueError(f"no closest-vector model requested for {root_type}")
    assert ambient * ambient.transpose() == cartan
    return cartan, ambient, valid


def closest_root_squared(root_type, root, center):
    standard, ambient_simple, valid = standard_root_model(root_type)
    isometry = matrix(ZZ, pari(standard).qfisom(pari(root)))
    assert abs(isometry.det()) == 1
    assert isometry.transpose() * root * isometry == standard
    target = ambient_simple.transpose() * isometry.inverse() * center
    choices = []
    for value in target:
        lower, upper = floor(value), ceil(value)
        choices.append((lower,) if lower == upper else (lower, upper))
    minimum = None
    for row in itertools.product(*choices):
        if not valid(row):
            continue
        difference = vector(QQ, row) - target
        value = difference * difference
        if minimum is None or value < minimum:
            minimum = value
    assert minimum is not None
    return minimum


def section_distance_profile(root_type, root, root_mw, height, witness, root_rank):
    z = vector(ZZ, witness[root_rank:])
    denominator = lcm(value.denominator() for value in height.list())
    scaled_height = (denominator * height).change_ring(ZZ)
    bound = 8 * denominator - 1
    short = pari(scaled_height).qfminim(bound)
    half = matrix(ZZ, short[2]).transpose().rows()
    candidates = [vector(ZZ, [0] * len(z))] + list(half) + [-row for row in half]
    candidates = [
        row for row in candidates
        if all((row[index] - z[index]) % 2 == 0 for index in range(len(z)))
        and QQ(row * scaled_height * row) / (4 * denominator) < 2
    ]
    if not candidates:
        return ()

    root_coordinate = vector(QQ, witness[:root_rank])
    distances = []
    for n in candidates:
        m = vector(QQ, (n + z) / 2)
        quotient_difference = m - vector(QQ, z) / 2
        center = (
            root_coordinate / 2
            - root.inverse() * root_mw * quotient_difference
        )
        root_distance = closest_root_squared(root_type, root, center)
        quotient_distance = QQ(n * scaled_height * n) / (4 * denominator)
        distances.append(root_distance + quotient_distance)
    return tuple(sorted(distances))


for case in CASES:
    frame = load_gram(case["frame"])
    assert frame.nrows() == 17 and frame.det() == 948
    root_rank = case["root_rank"]
    root = frame[:root_rank, :root_rank]
    root_mw = frame[:root_rank, root_rank:]
    height = (
        frame[root_rank:, root_rank:]
        - frame[root_rank:, :root_rank] * root.inverse() * root_mw
    )
    witness = vector(ZZ, case["witness"])
    assert witness * frame * witness == 4 * case["a"]

    ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -frame)
    old_fiber = vector(ZZ, [1, 0] + [0] * 17)
    old_zero = vector(ZZ, [-1, 1] + [0] * 17)
    divisor = vector(ZZ, [case["a"], 2] + list(witness))
    assert divisor * ns * divisor == 0
    assert divisor * ns * old_fiber == 2
    assert divisor * ns * old_zero == case["a"] - 2
    assert gcd(tuple(ns * divisor)) == 1

    effective_simple = tuple(
        vector(ZZ, [0, 0] + [
            -ZZ(index == node) for index in range(17)
        ])
        for node in range(root_rank)
    )
    component_pairings = tuple(
        divisor * ns * curve for curve in effective_simple
    )
    assert component_pairings == case["component_pairings"]
    affine_curves = tuple(
        old_fiber + vector(
            ZZ, [0, 0] + list(highest) + [0] * (17 - root_rank)
        )
        for highest in highest_roots(root)
    )
    affine_pairings = tuple(divisor * ns * curve for curve in affine_curves)
    assert affine_pairings == case["affine_pairings"]
    assert all(value >= 0 for value in component_pairings + affine_pairings)

    distances = section_distance_profile(
        case["root_type"], root, root_mw, height, witness, root_rank
    )
    assert distances == case["section_distances"]
    assert not distances or min(distances) >= 2

    # For C=[k,2,w], ||w-v||^2=2(D.C+1).  A negative intersection
    # forces w=v and D.C=-1, but then v^2=4k+2, contradicting parity.
    assert (witness * frame * witness - 2) % 4 != 0

    child = neighbor_frame(ns, divisor)
    minimum = pari(child).qfminim(2)
    child_roots = matrix(ZZ, minimum[2]).transpose()
    child_root_basis = child_roots.row_module().basis_matrix()
    child_root_gram = child_root_basis * child * child_root_basis.transpose()
    root_data = (
        child_root_basis.rank(),
        ZZ(minimum[0]),
        abs(child_root_gram.det()),
    )
    assert child.det() == 948 and root_data == case["child_root_data"]

    print(
        "H3CHAIN|step={}|q={}|ab={},2|old_degree=2|O={}|reflections=0|"
        "component_pairings={}|affine_pairings={}|section_profile={}|"
        "bisection_parity=1|nef=1|child={}|root_data={}|status=PASS".format(
            case["name"], case["q"], case["a"], case["a"] - 2,
            ",".join(map(str, component_pairings)),
            ",".join(map(str, affine_pairings)),
            "empty" if not distances else ",".join(map(str, distances)),
            case["child"], ",".join(map(str, root_data)),
        ),
        flush=True,
    )

print(
    "H3CHAIN|chain=D12/MW5-q6-A11/MW6-q8-A5+A5/MW7-q4-3A3/MW8-"
    "q4-A3+2A2/MW10|"
    "all_old_degree=2|status=PASS_H3_RANK_GROWING_DEGREE2_CHAIN",
    flush=True,
)
