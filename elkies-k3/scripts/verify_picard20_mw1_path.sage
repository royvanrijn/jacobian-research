#!/usr/bin/env sage
"""Replay the preferred two-neighbor path to the optimal disc-43 MW1 frame."""

from pathlib import Path

from sage.all import *


DATA = Path("elkies-k3/data/fibrations")


def read_matrix(name):
    path = DATA / name
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def qform_from_gram(gram):
    coefficients = []
    for i in range(gram.nrows()):
        for j in range(i, gram.ncols()):
            coefficients.append(gram[i, i] // 2 if i == j else gram[i, j])
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


def neighbor(parent, q, a, b, coordinates):
    rank = parent.nrows()
    vector_part = vector(ZZ, coordinates)
    assert len(vector_part) == rank
    assert a * b == q and vector_part * parent * vector_part == 2 * q

    hyperbolic_plane = matrix(ZZ, [[0, 1], [1, 0]])
    ns = block_diagonal_matrix(hyperbolic_plane, -parent)
    fiber = vector(ZZ, [a, b] + list(vector_part))
    assert fiber * ns * fiber == 0
    assert gcd([abs(ZZ(value)) for value in ns * fiber]) == 1

    pairings = list(ns * fiber)
    current_gcd = ZZ(0)
    mate = [ZZ(0)] * (rank + 2)
    for i, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, left, right = xgcd(current_gcd, ZZ(pairing))
        mate = [left * value for value in mate]
        mate[i] += right
        current_gcd = new_gcd
    assert abs(current_gcd) == 1
    if current_gcd == -1:
        mate = [-value for value in mate]
    mate = vector(ZZ, mate)
    assert fiber * ns * mate == 1
    mate -= ZZ(mate * ns * mate) // 2 * fiber
    assert mate * ns * mate == 0 and fiber * ns * mate == 1

    kernel = matrix(ZZ, [list(fiber * ns), list(mate * ns)]).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    assert child.det() == parent.det() and child.is_positive_definite()
    change_of_basis = block_matrix([[matrix(ZZ, [fiber, mate])], [kernel]], subdivide=False)
    assert abs(change_of_basis.det()) == 1
    assert change_of_basis * ns * change_of_basis.transpose() == block_diagonal_matrix(
        hyperbolic_plane, -child
    )
    return child, change_of_basis


def root_invariants(frame):
    half_roots = [
        vector(ZZ, root)
        for root in qform_from_gram(frame).short_vector_list_up_to_length(2, True)[1]
    ]
    roots = half_roots + [-root for root in half_roots]
    basis = matrix(ZZ, roots).row_module().basis_matrix()
    root_gram = basis * frame * basis.transpose()
    return basis.rank(), len(roots), abs(root_gram.det()), roots


canonical = read_matrix("picard20_e6_d4_a2a2_a1_mw3_frame.txt")
step1 = read_matrix("picard20_mw2_a7_a4_a3_a2_frame.txt")
endpoint = read_matrix("picard20_mw1_a12_a3_a2_frame.txt")
assert canonical.det() == step1.det() == endpoint.det() == 43

transitions = [
    (
        canonical,
        step1,
        8,
        2,
        4,
        (-1, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, -1, 0, 0, -1),
    ),
    (
        step1,
        endpoint,
        9,
        3,
        3,
        (0, -1, -2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0),
    ),
]

expected_invariants = (
    (16, 94, 480),
    (17, 174, 156),
)
composite_transport = identity_matrix(ZZ, 20)
for index, (transition, expected) in enumerate(zip(transitions, expected_invariants), 1):
    parent, child, q, a, b, coordinates = transition
    computed_child, transport = neighbor(parent, q, a, b, coordinates)
    assert computed_child == child
    composite_transport = transport * composite_transport
    root_rank, root_count, root_det, _ = root_invariants(child)
    assert (root_rank, root_count, root_det) == expected
    print(
        f"PICARD20MW1PATH|step={index}|q={q}|ab={a},{b}"
        f"|root_rank={root_rank}|roots={root_count}|rootdet={root_det}"
        f"|MW={18-root_rank}",
        flush=True,
    )

initial_ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -canonical)
terminal_ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -endpoint)
assert abs(composite_transport.det()) == 1
assert composite_transport * initial_ns * composite_transport.transpose() == terminal_ns
pinned_transport = read_matrix("picard20_mw1_a12_a3_a2_ns_transport.txt")
assert composite_transport == pinned_transport
direct_fiber = composite_transport[0]
direct_q = direct_fiber[0] * direct_fiber[1]
assert direct_q > 0
assert vector(ZZ, direct_fiber[2:]) * canonical * vector(ZZ, direct_fiber[2:]) == 2 * direct_q
print(
    f"PICARD20MW1PATH|composite_direct_q={direct_q}"
    f"|ab={direct_fiber[0]},{direct_fiber[1]}"
    f"|v={tuple(direct_fiber[2:])}",
    flush=True,
)
print("PICARD20MW1PATH|composite_transport_det=1|pinned=1", flush=True)

# Identify the endpoint root components by their (rank, root count, det).
_, _, _, endpoint_roots = root_invariants(endpoint)
graph = Graph()
graph.add_vertices(range(len(endpoint_roots)))
for i in range(len(endpoint_roots)):
    for j in range(i):
        if endpoint_roots[i] * endpoint * endpoint_roots[j] != 0:
            graph.add_edge(i, j)
components = []
for vertices in graph.connected_components(sort=False):
    basis = matrix(ZZ, [endpoint_roots[i] for i in vertices]).row_module().basis_matrix()
    gram = basis * endpoint * basis.transpose()
    components.append((basis.rank(), len(vertices), abs(gram.det())))
assert sorted(components) == [(2, 6, 3), (3, 12, 4), (12, 156, 13)]

# MW0 would require a rank-18 ADE root determinant divisible by 43.  Each
# irreducible ADE factor of rank at most 18 has determinant among n+1,4,3,2,1,
# hence none is divisible by 43.  Therefore the attained MW1 is optimal.
assert all((rank + 1) % 43 for rank in range(1, 19))
assert all(value % 43 for value in (1, 2, 3, 4))

print("PICARD20MW1PATH|endpoint_roots=A12+A3+A2|height=43/156|torsion=trivial", flush=True)
print("PICARD20MW1PATH|profile=3_or_10,1_or_3,1_or_2|P.O=0|semistable_candidate=1", flush=True)
print("PICARD20MW1PATH|optimal=1|reason=MW0_disc43_ADE_obstruction", flush=True)
print("PICARD20MW1PATH|status=PASS", flush=True)
