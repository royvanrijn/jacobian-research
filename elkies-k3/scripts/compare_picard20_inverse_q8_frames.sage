#!/usr/bin/env sage
"""Compare the transported first inverse pencil with the canonical disc-43 q=8 move.

The determinant-948 inverse class is first embedded in the saturated Picard-20
NS lattice obtained by adjoining the third section.  All comparisons are then
made in the pinned U + (-picard20 frame) coordinates.
"""

from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3" / "data" / "fibrations"
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def load_backward_classes(path):
    rows = [
        line.split("\t")
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    assert rows[0] == ["stage", "class", "coordinates"]
    return {
        (stage, class_name): vector(ZZ, coordinates.split(","))
        for stage, class_name, coordinates in rows[1:]
    }


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    if current == -1:
        coefficients = [-value for value in coefficients]
    return vector(ZZ, coefficients)


def split_isotropic(ns, fiber):
    """Return the negative-definite complement frame and its basis transport."""
    fiber = vector(ZZ, fiber)
    assert fiber * ns * fiber == 0
    assert gcd([abs(ZZ(value)) for value in ns * fiber]) == 1
    mate = bezout_vector(ns * fiber)
    assert fiber * ns * mate == 1
    mate_square = ZZ(mate * ns * mate)
    assert mate_square % 2 == 0
    mate -= (mate_square // 2) * fiber
    assert mate * ns * mate == 0 and fiber * ns * mate == 1
    kernel = matrix(ZZ, [list(fiber * ns), list(mate * ns)]).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    transport = matrix(ZZ, [list(fiber), list(mate)] + kernel.rows())
    assert abs(transport.det()) == 1
    assert transport * ns * transport.transpose() == block_diagonal_matrix(U, -child)
    return child, transport


def qform_from_gram(gram):
    coefficients = []
    for row in range(gram.nrows()):
        for column in range(row, gram.ncols()):
            coefficients.append(
                gram[row, row] // 2 if row == column else gram[row, column]
            )
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


def root_components(frame):
    half_roots = [
        vector(ZZ, root)
        for root in qform_from_gram(frame).short_vector_list_up_to_length(2, True)[1]
    ]
    roots = half_roots + [-root for root in half_roots]
    if not roots:
        return 0, 0, (), roots
    graph = Graph()
    graph.add_vertices(range(len(roots)))
    for left in range(len(roots)):
        for right in range(left):
            if roots[left] * frame * roots[right] != 0:
                graph.add_edge(left, right)
    components = []
    for vertices in graph.connected_components(sort=False):
        basis = matrix(ZZ, [roots[index] for index in vertices]).row_module().basis_matrix()
        gram = basis * frame * basis.transpose()
        components.append((basis.rank(), len(vertices), abs(gram.det())))
    components.sort()
    root_basis = matrix(ZZ, roots).row_module().basis_matrix()
    return root_basis.rank(), len(roots), tuple(components), roots


def ade_name(components):
    irreducible = {
        (1, 2, 2): "A1",
        (2, 6, 3): "A2",
        (3, 12, 4): "A3",
        (4, 20, 5): "A4",
        (4, 24, 4): "D4",
        (5, 30, 6): "A5",
        (7, 56, 8): "A7",
    }
    assert all(component in irreducible for component in components)
    names = [irreducible[component] for component in components]
    grouped = []
    for name in dict.fromkeys(names):
        multiplicity = names.count(name)
        grouped.append((str(multiplicity) if multiplicity > 1 else "") + name)
    return "+".join(grouped)


def reduce_against_explicit_curves(divisor, gram):
    reduced = vector(ZZ, list(divisor))
    reflections = []
    for _ in range(10000):
        pairings = reduced * gram
        negative = [index for index in range(1, 19) if pairings[index] < 0]
        if not negative:
            return reduced, tuple(reflections)
        index = min(negative, key=lambda candidate: pairings[candidate])
        multiplicity = -ZZ(pairings[index])
        reduced[index] -= multiplicity
        reflections.append((index, multiplicity))
    raise RuntimeError("explicit-curve reflection reduction did not terminate")


# Reconstruct the saturated Picard-20 lattice exactly as in
# verify_picard20_ns_extension.sage.
old_frame = load_matrix(DATA / "mw2_e6_d4_a2a2_a1_frame.txt")
old_basis = load_matrix(DATA / "mw2_e6_d4_a2a2_a1_explicit_basis.txt")
old_ns = block_diagonal_matrix(U, -old_frame)
old_gram = old_basis * old_ns * old_basis.transpose()
s_pairings = [1, 0] + [0] * 15 + [0, 2]
s_pairings[2 + 2] = 1
s_pairings[2 + 6] = 1
extended_gram = block_matrix([
    [old_gram, matrix(ZZ, 19, 1, s_pairings)],
    [matrix(ZZ, 1, 19, s_pairings), matrix(ZZ, [[-2]])],
])
assert extended_gram.det() == -43

split = identity_matrix(ZZ, 20)
split[1] = vector(ZZ, [1, 1] + [0] * 18)
split[17] = vector(ZZ, [-2, -1] + [0] * 15 + [1, 0, 0])
split[18] = vector(ZZ, [-3, -1] + [0] * 15 + [0, 1, 0])
split[19] = vector(ZZ, [-2, -1] + [0] * 15 + [0, 0, 1])
canonical_frame = load_matrix(DATA / "picard20_e6_d4_a2a2_a1_mw3_frame.txt")
canonical_ns = block_diagonal_matrix(U, -canonical_frame)
assert split * extended_gram * split.transpose() == canonical_ns
assert split.det() == 1

backward = load_backward_classes(
    DATA / "mw2_e6_d4_a2a2_a1_inverse_neighbor_classes.tsv"
)
raw_explicit = backward[("a5_d4_2a2_a1_mw3", "fiber")]
reduced_explicit, reflections = reduce_against_explicit_curves(raw_explicit, old_gram)
assert len(reflections) == 42


def into_canonical(explicit_vector):
    extended = vector(ZZ, list(explicit_vector) + [0])
    canonical = extended * split.inverse()
    assert canonical in ZZ ** 20
    assert canonical * canonical_ns * canonical == 0
    return vector(ZZ, canonical)


raw_fiber = into_canonical(raw_explicit)
reduced_fiber = into_canonical(reduced_explicit)

stage_names = (
    "a5_d4_2a2_a1_mw3",
    "q25_mw4",
    "q25_mw7",
    "rank17",
)
stage_results = {}
for stage_name in stage_names:
    stage_fiber = into_canonical(backward[(stage_name, "fiber")])
    stage_child, stage_transport = split_isotropic(canonical_ns, stage_fiber)
    stage_invariants = root_components(stage_child)[:3]
    stage_results[stage_name] = (
        stage_fiber,
        stage_child,
        stage_transport,
        stage_invariants,
    )

q8_witness = vector(ZZ, (
    -1, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, -1, 0, 0, -1,
))
q8_fiber = vector(ZZ, [2, 4] + list(q8_witness))
assert q8_fiber * canonical_ns * q8_fiber == 0

raw_child, raw_transport = split_isotropic(canonical_ns, raw_fiber)
reduced_child, reduced_transport = split_isotropic(canonical_ns, reduced_fiber)
q8_child, q8_transport = split_isotropic(canonical_ns, q8_fiber)
assert q8_child == load_matrix(DATA / "picard20_mw2_a7_a4_a3_a2_frame.txt")
assert raw_child == stage_results["a5_d4_2a2_a1_mw3"][1]

raw_invariants = root_components(raw_child)[:3]
reduced_invariants = root_components(reduced_child)[:3]
q8_invariants = root_components(q8_child)[:3]
assert raw_invariants == reduced_invariants
assert raw_invariants == (
    15,
    80,
    ((2, 6, 3), (4, 20, 5), (4, 24, 4), (5, 30, 6)),
)
assert q8_invariants[0] == 16
assert ade_name(raw_invariants[2]) == "A2+A4+D4+A5"
assert ade_name(q8_invariants[2]) == "A2+A3+A4+A7"
assert not pari(raw_child).qfisom(pari(q8_child))
assert pari(raw_child).qfisom(pari(reduced_child))

expected_stage_ade = {
    "a5_d4_2a2_a1_mw3": "A2+A4+D4+A5",
    "q25_mw4": "A1+A2+A3+A4+D4",
    "q25_mw7": "A1+A2+A3+A4+D4",
    "rank17": "4A1+3A2+A4",
}
assert {
    stage_name: ade_name(stage_results[stage_name][3][2])
    for stage_name in stage_names
} == expected_stage_ade

fiber_intersection = ZZ(raw_fiber * canonical_ns * q8_fiber)
reduced_intersection = ZZ(reduced_fiber * canonical_ns * q8_fiber)
q8_in_raw_child = vector(ZZ, q8_fiber * raw_transport.inverse())
raw_in_q8_child = vector(ZZ, raw_fiber * q8_transport.inverse())
assert q8_in_raw_child * block_diagonal_matrix(U, -raw_child) * q8_in_raw_child == 0
assert raw_in_q8_child * block_diagonal_matrix(U, -q8_child) * raw_in_q8_child == 0
assert q8_in_raw_child[1] == fiber_intersection
assert raw_in_q8_child[1] == fiber_intersection

print(
    f"PICARD20INVERSEQ8|raw_fiber={tuple(raw_fiber)}"
    f"|canonical_ab={raw_fiber[0]},{raw_fiber[1]}"
    f"|canonical_q={raw_fiber[0] * raw_fiber[1]}",
    flush=True,
)
print(
    f"PICARD20INVERSEQ8|reduced_fiber={tuple(reduced_fiber)}"
    f"|reflections={len(reflections)}"
    f"|canonical_ab={reduced_fiber[0]},{reduced_fiber[1]}"
    f"|canonical_q={reduced_fiber[0] * reduced_fiber[1]}",
    flush=True,
)
print(
    f"PICARD20INVERSEQ8|transported_child_root_rank={raw_invariants[0]}"
    f"|roots={raw_invariants[1]}|components={raw_invariants[2]}"
    f"|geometric_MW={18 - raw_invariants[0]}",
    flush=True,
)
print(
    f"PICARD20INVERSEQ8|q8_child_root_rank={q8_invariants[0]}"
    f"|roots={q8_invariants[1]}|components={q8_invariants[2]}"
    f"|geometric_MW={18 - q8_invariants[0]}",
    flush=True,
)
for stage_name in stage_names:
    stage_fiber, _, _, stage_invariants = stage_results[stage_name]
    print(
        f"PICARD20INVERSEPATH|stage={stage_name}"
        f"|canonical_ab={stage_fiber[0]},{stage_fiber[1]}"
        f"|canonical_q={stage_fiber[0] * stage_fiber[1]}"
        f"|root_rank={stage_invariants[0]}|roots={stage_invariants[1]}"
        f"|ADE={ade_name(stage_invariants[2])}"
        f"|geometric_MW={18 - stage_invariants[0]}",
        flush=True,
    )
print(
    f"PICARD20INVERSEQ8|fiber_intersection_raw={fiber_intersection}"
    f"|fiber_intersection_reduced={reduced_intersection}",
    flush=True,
)
print(
    f"PICARD20INVERSEQ8|q8_in_transported_child={tuple(q8_in_raw_child)}"
    f"|neighbor_q={q8_in_raw_child[0] * q8_in_raw_child[1]}",
    flush=True,
)
print(
    f"PICARD20INVERSEQ8|transported_in_q8_child={tuple(raw_in_q8_child)}"
    f"|neighbor_q={raw_in_q8_child[0] * raw_in_q8_child[1]}",
    flush=True,
)
print(
    "PICARD20INVERSEPATH|terminal_rootless=0|terminal_ADE=4A1+3A2+A4"
    "|terminal_geometric_MW=4",
    flush=True,
)
print("PICARD20INVERSEQ8|same_frame=0|raw_reduced_isometric=1|status=PASS", flush=True)
