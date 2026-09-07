#!/usr/bin/env sage
"""Reduce the first q=4 class in the q80 fibration chamber.

This is an equation-planning certificate, not yet a Riemann--Roch execution.
It chooses a deterministic simple-root chamber for the pinned
E6+D5+A3 frame, reflects the first q80-to-rootless fiber against the zero
section and all old-fiber components, and proves nonnegativity against every
section by reducing the only remaining test to the rank-three MW shell of
height less than two. Primitivity of the old root lattice then excludes the
only other possible fixed curve, a degree-two multisection, and proves the
reduced class nef.
"""

import argparse
from pathlib import Path

from sage.all import (
    QQ,
    ZZ,
    Matrix,
    QuadraticForm,
    block_diagonal_matrix,
    matrix,
    vector,
)
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
FRAME_PATH = (
    ROOT / "elkies-k3/data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt"
)
PATH_TSV = ROOT / "elkies-k3/data/fibrations/kumar_q80_to_rootless_path.tsv"
U = matrix(ZZ, [[0, 1], [1, 0]])

parser = argparse.ArgumentParser()
parser.add_argument("--cvp-candidates", type=int, default=128)
args = parser.parse_args()


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def roots_of_norm_two(gram):
    half = QuadraticForm(ZZ, gram).short_vector_list_up_to_length(
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
        decomposable = any(
            tuple(row - left) in positive_set
            for left in positive
            if left != row
        )
        if not decomposable:
            simple.append(row)
    simple = matrix(ZZ, [list(row) for row in simple])
    assert simple.nrows() == simple.rank() == 14
    cartan = simple * gram * simple.transpose()
    assert set(cartan.diagonal()) == {2}
    assert all(
        cartan[i, j] in (0, -1)
        for i in range(cartan.nrows())
        for j in range(cartan.ncols())
        if i != j
    )
    return simple, positive


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    answer = []
    while unseen:
        todo = [min(unseen)]
        component = []
        unseen.remove(todo[0])
        while todo:
            index = todo.pop()
            component.append(index)
            adjacent = [
                other
                for other in list(unseen)
                if cartan[index, other] != 0
            ]
            for other in adjacent:
                unseen.remove(other)
                todo.append(other)
        answer.append(tuple(sorted(component)))
    return tuple(sorted(answer, key=lambda component: (len(component), component)))


def highest_roots(gram, simple, positive):
    cartan = simple * gram * simple.transpose()
    components = connected_components(cartan)
    inverse_simple = simple.pseudoinverse()
    result = []
    for component in components:
        candidates = []
        for root in positive:
            coordinates = vector(QQ, root) * inverse_simple
            if not all(value in ZZ and value >= 0 for value in coordinates):
                continue
            support = tuple(index for index, value in enumerate(coordinates) if value)
            if support and all(index in component for index in support):
                candidates.append((sum(coordinates), vector(ZZ, root), coordinates))
        height, root, coordinates = max(candidates, key=lambda item: item[0])
        result.append((component, root, vector(ZZ, coordinates), ZZ(height)))
    return result


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


frame = load_matrix(FRAME_PATH)
ns = block_diagonal_matrix(U, -frame)
simple, positive = deterministic_simple_roots(frame)
cartan = simple * frame * simple.transpose()
components = connected_components(cartan)
assert sorted(map(len, components)) == [3, 5, 6]

zero = vector(ZZ, [-1, 1] + [0] * 17)
fiber = vector(ZZ, [1, 0] + [0] * 17)
curves = [("O", zero)]
for index, root in enumerate(simple.rows()):
    curves.append((f"R{index + 1}", vector(ZZ, [0, 0] + list(root))))
for component_index, (_, root, coordinates, _) in enumerate(
    highest_roots(frame, simple, positive), start=1
):
    theta0 = fiber - vector(ZZ, [0, 0] + list(root))
    curves.append((f"Theta0_{component_index}", theta0))
    assert theta0 * ns * theta0 == -2

first_row = PATH_TSV.read_text().splitlines()[1].split("\t")
assert first_row[:4] == ["1", "4", "2", "2"]
witness = list(map(ZZ, first_row[4].split(",")))
raw = vector(ZZ, [2, 2] + witness)
assert raw * ns * raw == 0

reduced, root_sequence = chamber_reduce(raw, curves, ns)
assert intersection(reduced, fiber, ns) == 2
reduced_frame = vector(ZZ, reduced[2:])
root_coordinates = vector(QQ, reduced_frame) * simple.pseudoinverse()
assert all(value in ZZ for value in root_coordinates)
root_coordinates = vector(ZZ, root_coordinates)

print(
    "Q80FIRSTCHAMBER|components={}|simple_cartans={}".format(
        tuple(map(len, components)),
        tuple(
            tuple(tuple(cartan[i, j] for j in component) for i in component)
            for component in components
        ),
    ),
    flush=True,
)
print(
    "Q80FIRSTCHAMBER|affine_fiber_multiplicities={}".format(
        tuple(
            (
                tuple(index + 1 for index in component),
                tuple(coordinates[index] for index in component),
            )
            for component, _, coordinates, _ in highest_roots(
                frame, simple, positive
            )
        )
    ),
    flush=True,
)
print(f"Q80FIRSTCHAMBER|raw={tuple(raw)}|D.F={intersection(raw, fiber, ns)}|D.O={intersection(raw, zero, ns)}", flush=True)
print(f"Q80FIRSTCHAMBER|fiber_reduction={root_sequence}", flush=True)
print(f"Q80FIRSTCHAMBER|after_fibers={tuple(reduced)}", flush=True)
print(
    f"Q80FIRSTCHAMBER|trivial_expression=4F+2O+root({tuple(root_coordinates)})",
    flush=True,
)
print(
    "Q80FIRSTCHAMBER|fiber_pairings={}".format(
        tuple((name, intersection(reduced, curve, ns)) for name, curve in curves)
    ),
    flush=True,
)

# Work in the exact saturated generic MW basis used by the q80 profile
# certificate.  These integral representatives are deliberately not reduced
# against the root lattice; an exact closest-vector calculation supplies a
# short representative of each tested MW coset without enumerating the full
# 17-dimensional frame shell.
root_basis = matrix(
    ZZ, [list(row) for row in roots_of_norm_two(frame)]
).row_module().basis_matrix()
assert root_basis.nrows() == 14
root_smith = root_basis.smith_form()[0]
assert tuple(abs(root_smith[i, i]) for i in range(14)) == (1,) * 14
root_gram = root_basis * frame * root_basis.transpose()
optimal_lifts = matrix(
    ZZ,
    [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [-2, -3, -4, -6, -5, -4, 0, 0, 0, 0, 0, 0, 5, 14, -2, 1, 2],
        [6, 9, 12, 18, 15, 12, 0, 0, 0, 0, 0, 0, -15, -40, 7, -3, -6],
    ],
)
root_projection = (
    optimal_lifts
    * frame
    * root_basis.transpose()
    * root_gram.inverse()
)
mw_projection = optimal_lifts - root_projection * root_basis
height = mw_projection * frame * mw_projection.transpose()
assert height == matrix(
    QQ,
    [
        [QQ(2) / 3, 0, QQ(3) / 4],
        [0, 4, 2],
        [QQ(3) / 4, 2, QQ(37) / 4],
    ],
)
root_lattice = IntegralLattice(root_gram)


def short_coset_representative(mw_coordinates, expected_norm):
    raw_lift = vector(ZZ, mw_coordinates) * optimal_lifts
    root_coordinates = (
        raw_lift * frame * root_basis.transpose() * root_gram.inverse()
    )
    iterator = root_lattice.enumerate_close_vectors(-root_coordinates)
    for _ in range(args.cvp_candidates):
        shift = vector(ZZ, next(iterator))
        lift = raw_lift + shift * root_basis
        norm = ZZ(lift * frame * lift)
        if norm == expected_norm:
            return lift
    raise RuntimeError(
        f"failed to locate norm-{expected_norm} representative of MW coset "
        f"{tuple(mw_coordinates)} in {args.cvp_candidates} CVP candidates"
    )


# First certify the three published optimal representatives have norms
# 4,4,10, equivalently P.O=(0,0,3).
for index, expected_norm in enumerate((4, 4, 10)):
    coordinates = vector(ZZ, [1 if i == index else 0 for i in range(3)])
    lift = short_coset_representative(coordinates, expected_norm)
    assert lift * frame * lift == expected_norm

# Since the reduced class has zero MW projection and frame norm eight, a
# section with MW height h and root correction c has
#
#   D.S = h+c-<D_root,S_root> >= h+c-sqrt(8c) >= h-2.
#
# It therefore suffices to enumerate h<2 in the rank-three height lattice.
# With Q(x)=12*h(x), the strict bound is Q<24.
assert reduced_frame * frame * reduced_frame == 8
reduced_mw_coordinates = (
    reduced_frame * frame * mw_projection.transpose() * height.inverse()
)
assert reduced_mw_coordinates == vector(QQ, [0, 0, 0])
scaled_height = (24 * height).change_ring(ZZ)
mw_form = QuadraticForm(ZZ, scaled_height)
shells = mw_form.short_vector_list_up_to_length(24, up_to_sign_flag=True)
low_height = []
for shell in shells[1:]:
    for row in shell:
        row = vector(ZZ, row)
        low_height.extend((row, -row))
low_height = sorted(set(map(tuple, low_height)))
assert low_height == [(-1, 0, 0), (1, 0, 0)]


def local_self_correction(coordinates):
    left, _, right = map(ZZ, coordinates)
    a3 = (left + 3 * right) % 4
    d5 = left % 4
    e6 = left % 3
    correction_a3 = QQ(a3 * (4 - a3)) / 4
    correction_d5 = (
        QQ(0) if d5 == 0 else QQ(1) if d5 == 2 else QQ(5) / 4
    )
    correction_e6 = QQ(0) if e6 == 0 else QQ(4) / 3
    return correction_a3 + correction_d5 + correction_e6


sections = []
for coordinates_tuple in low_height:
    coordinates = vector(ZZ, coordinates_tuple)
    section_height = vector(QQ, coordinates) * height * vector(QQ, coordinates)
    correction = local_self_correction(coordinates)
    norm = section_height + correction
    assert norm in ZZ and norm >= 4 and norm % 2 == 0
    lift = short_coset_representative(coordinates, ZZ(norm))
    pole = ZZ(norm) // 2 - 2
    section = vector(ZZ, [pole + 1, 1] + list(lift))
    assert section * ns * section == -2
    pairing = intersection(reduced, section, ns)
    assert pairing >= 0
    sections.append((pairing, pole, tuple(coordinates), tuple(lift), section))
sections.sort(key=lambda item: (item[0], item[1], item[2], item[3]))

print(
    f"Q80FIRSTCHAMBER|section_height_lt_2={tuple(low_height)}|cvp_candidates={args.cvp_candidates}|tested={len(sections)}",
    flush=True,
)
for index, (pairing, pole, coordinates, lift, _) in enumerate(sections):
    print(
        f"Q80FIRSTCHAMBER|low_section={index}|pairing={pairing}|P.O={pole}|MW={coordinates}|lift={lift}",
        flush=True,
    )

print(
    "Q80FIRSTCHAMBER|all_sections_nonnegative=1|"
    "root_primitive=1|MW_torsion=1",
    flush=True,
)

# Complete the nefness argument.  If an irreducible (-2)-curve C had D.C<0,
# then C would be fixed in the effective class D.  Hence 1 <= C.F <= D.F=2.
# Degree one is a section, already handled above.  At degree two, D-C is
# vertical, so C has zero MW projection.  Root primitivity therefore gives
#
#   C = kF + 2O + r,  r in R=A3+D5+E6.
#
# From C^2=-2 and D=4F+2O+rho with rho^2=8,
#
#   r^2=4k-6,  D.C=2k-<rho,r>,
#   ||r-rho||^2 = 2(D.C+1).
#
# If D.C<0, integrality makes D.C<=-1.  Positive definiteness forces equality,
# r=rho, and D.C=-1; then 8=4k-6 gives the impossible k=7/2.
assert intersection(reduced, fiber, ns) == 2
assert reduced_frame * frame * reduced_frame == 8
assert QQ(8 + 6) / 4 == QQ(7) / 2
print(
    "Q80FIRSTCHAMBER|negative_bisection_impossible=1|"
    "all_effective_minus2_nonnegative=1|nef=1|status=PASS_NEF_CHAMBER",
    flush=True,
)
