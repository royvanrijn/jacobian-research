#!/usr/bin/env sage -python
"""Certify nefness of the degree-two H3 suffix from MW10 to MW17.

The root-adapted quotient searches select

    A3+2A2/MW10 --q4--> 5A1/MW12 --q4--> 4A1/MW13
      --q4--> 3A1/MW14 --q4--> 2A1/MW15 --q4--> A1/MW16
      --q6--> rootless/MW17.

Every displayed class has old-fiber degree two and already lies in the
component/zero chamber.  Full nefness is exact: sections are one closest-
vector problem in the full positive frame, while a negative bisection is
excluded by even-lattice parity.  No bounded section enumeration is used.
"""

import hashlib
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]

CASES = (
    {
        "name": "A3A2A2Q4",
        "frame": "artifacts/generated-results/elkies-k3-h3-a3x3-q4-degree2-frames/q4-o0323-r7-n24-d36-87b284dff2bc.txt",
        "frame_sha256": "826c8e31ecd04d37e50e3413e518c58660cca8bfb70c39c73d8330c4f8acbf08",
        "root_rank": 7,
        "source": "A3+2A2/MW10",
        "source_root_data": (7, 24, 36),
        "q": 4,
        "a": 2,
        "orbit": 207,
        "witness": (-8, 7, -5, 6, -1, 2, 2, -1, 1, 0, -1, 1, 1, -1, 0, 1, -1),
        "component_pairings": (0, 0, 0, 1, 1, 0, 1),
        "child": "5A1/MW12",
        "child_root_data": (5, 10, 32),
    },
    {
        "name": "5A1Q4",
        "frame": "artifacts/generated-results/elkies-k3-h3-mw10-a3a2a2-q4-degree2-frames/q4-o0207-r5-n10-d32-a462a553a1e9.txt",
        "frame_sha256": "3ab059a6d216f4f75942e954a9135af34be6308402a3b927b8b9718148ec13b5",
        "root_rank": 5,
        "source": "5A1/MW12",
        "source_root_data": (5, 10, 32),
        "q": 4,
        "a": 2,
        "orbit": 52,
        "witness": (-8, 11, 6, -4, -2, 0, -1, 0, 0, 0, 0, 1, -1, -1, 2, 0, -1),
        "component_pairings": (0, 0, 1, 1, 1),
        "child": "4A1/MW13",
        "child_root_data": (4, 8, 16),
    },
    {
        "name": "4A1Q4",
        "frame": "artifacts/generated-results/elkies-k3-h3-mw12-5a1-q4-degree2-first-hit-frames/q4-o0052-r4-n8-d16-066f47d7fff3.txt",
        "frame_sha256": "ab7cc81b2f55c477d7132f75132bd1916e6e32ebef21a7f32c789252e7f9ae99",
        "root_rank": 4,
        "source": "4A1/MW13",
        "source_root_data": (4, 8, 16),
        "q": 4,
        "a": 2,
        "orbit": 114,
        "witness": (-4, -1, 5, -1, -3, 0, 1, -2, 1, 0, 1, -1, 1, -1, 0, -1, -1),
        "component_pairings": (0, 1, 1, 0),
        "child": "3A1/MW14",
        "child_root_data": (3, 6, 8),
    },
    {
        "name": "3A1Q4",
        "frame": "artifacts/generated-results/elkies-k3-h3-mw13-4a1-q4-degree2-first-hit-frames/q4-o0114-r3-n6-d8-018b225c409b.txt",
        "frame_sha256": "bc03a4868c1c1f26c9283103076926e59046090270c071afef8b5f2901d2c57e",
        "root_rank": 3,
        "source": "3A1/MW14",
        "source_root_data": (3, 6, 8),
        "q": 4,
        "a": 2,
        "orbit": 498,
        "witness": (-11, 17, 0, -2, -1, 0, 0, 1, 1, 0, -1, 1, 1, 0, 1, -1, 1),
        "component_pairings": (1, 0, 1),
        "child": "2A1/MW15",
        "child_root_data": (2, 4, 4),
    },
    {
        "name": "2A1Q4",
        "frame": "artifacts/generated-results/elkies-k3-h3-mw14-3a1-q4-degree2-first-hit-frames/q4-o0498-r2-n4-d4-e86cf7c2d2f9.txt",
        "frame_sha256": "624dd2b2945288f29b8a548535bb47b269a6c9992894c01bddbdcd0170c87c8f",
        "root_rank": 2,
        "source": "2A1/MW15",
        "source_root_data": (2, 4, 4),
        "q": 4,
        "a": 2,
        "orbit": 981,
        "witness": (-5, -42, -1, 3, -1, -1, 2, 1, 2, 0, 1, 0, 1, 0, 1, 0, 1),
        "component_pairings": (1, 1),
        "child": "A1/MW16",
        "child_root_data": (1, 2, 2),
    },
    {
        "name": "A1Q6",
        "frame": "artifacts/generated-results/elkies-k3-h3-mw15-2a1-q4-degree2-first-hit-frames/q4-o0981-r1-n2-d2-4f02793cfc09.txt",
        "frame_sha256": "7391d3c8e8f1a3a7207724aca84d25d86f8cbce55144c82ac04a007bd5fc5bef",
        "root_rank": 1,
        "source": "A1/MW16",
        "source_root_data": (1, 2, 2),
        "q": 6,
        "a": 3,
        "orbit": 2247,
        "witness": (-2, -5, -1, 2, 3, -2, 1, 1, -1, 0, 2, 1, 0, -1, 1, -1, 3),
        "component_pairings": (1,),
        "child": "rootless/MW17",
        "child_root_data": (0, 0, 1),
    },
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    if cartan.nrows() == 0:
        return ()
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-row for row in half)
    result = []
    for component in connected_components(cartan):
        candidates = [
            row for row in roots
            if all(value >= 0 for value in row)
            and all(
                index in component or row[index] == 0
                for index in range(cartan.nrows())
            )
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


def root_data(frame):
    minimum = pari(frame).qfminim(2)
    root_count = ZZ(minimum[0])
    if root_count == 0:
        return (0, 0, 1)
    roots = matrix(ZZ, minimum[2]).transpose()
    basis = roots.row_module().basis_matrix()
    gram = basis * frame * basis.transpose()
    return (basis.rank(), root_count, abs(gram.det()))


for case in CASES:
    path = ROOT / case["frame"]
    assert digest(path) == case["frame_sha256"]
    frame = load_gram(path)
    assert frame.nrows() == 17 and frame.det() == 948
    assert all(frame[index, index] % 2 == 0 for index in range(17))
    assert frame.is_positive_definite()
    assert root_data(frame) == case["source_root_data"]

    root_rank = case["root_rank"]
    root = frame[:root_rank, :root_rank]
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
    assert all(value >= 0 for value in component_pairings + affine_pairings)
    assert divisor * ns * old_zero >= 0

    # A section is S_w=((w.M.w-2)/2,1,w).  Completing the square gives
    # D.S_w=(w-v/2).M.(w-v/2)-2.  The lattice CVP is exact.
    center = vector(QQ, witness) / 2
    closest = vector(ZZ, next(IntegralLattice(frame).enumerate_close_vectors(center)))
    closest_distance = (closest - center) * frame * (closest - center)
    minimum_section_pairing = closest_distance - 2
    assert closest_distance >= 2

    # For a bisection C=(k,2,w), one has
    # (w-v)^2/2=D.C+1.  A negative intersection forces w=v, but then
    # w^2=v^2 is 0 mod 4 instead of the required 2 mod 4.
    assert witness * frame * witness % 4 == 0

    child = neighbor_frame(ns, divisor)
    child_data = root_data(child)
    assert child.det() == 948 and child_data == case["child_root_data"]

    print(
        "H3MW17NEF|step={}|source={}|q={}|ab={},2|orbit={}|"
        "old_degree=2|O={}|reflections=0|component_pairings={}|"
        "affine_pairings={}|closest_section={}|closest_distance={}|"
        "minimum_section_pairing={}|bisection_parity=1|nef=1|child={}|"
        "root_data={}|status=PASS".format(
            case["name"], case["source"], case["q"], case["a"],
            case["orbit"], case["a"] - 2,
            ",".join(map(str, component_pairings)),
            ",".join(map(str, affine_pairings)),
            ",".join(map(str, closest)), closest_distance,
            minimum_section_pairing, case["child"],
            ",".join(map(str, child_data)),
        ),
        flush=True,
    )

print(
    "H3MW17NEF|chain=A3+2A2/MW10-q4-5A1/MW12-q4-4A1/MW13-"
    "q4-3A1/MW14-q4-2A1/MW15-q4-A1/MW16-q6-rootless/MW17|"
    "all_old_degree=2|all_nef=1|status=PASS_H3_MW10_TO_ROOTLESS_NEF",
    flush=True,
)
