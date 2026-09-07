#!/usr/bin/env sage
"""Reduce the second q=4 fiber in the explicit first-child lattice chamber.

This continues the equation planning after
``derive_q80_first_q4_pencil.sage``.  It reconstructs the first D9+A4/MW4
frame by exact neighbor transport, loads the second pinned q=4 witness, and
reflects it against the zero section and every old-fiber component.  The
result is a chamber-reduced divisor and an exact decomposition into root and
MW projections.  A saturated rank-four MW calculation checks the only four
section classes of height below two; a norm argument excludes a negative
bisection.  Thus the script proves full nefness.  It also reduces the next
pinned q=12 witness in the second-child D7+D5 frame, certifying that the next
equation is a degree-three divisor with nonzero MW projection.
"""

import csv
from pathlib import Path

from sage.all import (
    identity_matrix,
    lcm,
    QQ,
    ZZ,
    QuadraticForm,
    block_diagonal_matrix,
    gcd,
    matrix,
    pari,
    vector,
    xgcd,
)
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3/data/fibrations"
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)]*len(pairings)
    for index, pairing in enumerate(pairings):
        if not pairing:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left*value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    return vector(ZZ, coefficients if current == 1 else [-x for x in coefficients])


def neighbor(parent, qnorm, a, b, coordinates):
    ns = block_diagonal_matrix(U, -parent)
    fiber = vector(ZZ, [a, b]+list(coordinates))
    assert a*b == qnorm and coordinates*parent*coordinates == 2*qnorm
    assert fiber*ns*fiber == 0
    assert gcd([abs(ZZ(value)) for value in ns*fiber]) == 1
    mate = bezout_vector(list(ns*fiber))
    mate -= ZZ(mate*ns*mate)//2*fiber
    complement = matrix(
        ZZ, [list(fiber*ns), list(mate*ns)]
    ).right_kernel_matrix()
    child = -(complement*ns*complement.transpose())
    transport = matrix(ZZ, [list(fiber), list(mate)]+complement.rows())
    assert abs(transport.det()) == 1
    return child, transport


def roots_of_norm_two(gram):
    half = QuadraticForm(ZZ, gram).short_vector_list_up_to_length(
        2, up_to_sign_flag=True
    )[1]
    half = [vector(ZZ, row) for row in half]
    return half+[-row for row in half]


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
    assert simple.nrows() == simple.rank()
    return simple, positive


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    answer = []
    while unseen:
        todo = [min(unseen)]
        unseen.remove(todo[0])
        component = []
        while todo:
            index = todo.pop()
            component.append(index)
            for other in list(unseen):
                if cartan[index, other]:
                    unseen.remove(other)
                    todo.append(other)
        answer.append(tuple(sorted(component)))
    return tuple(sorted(answer, key=lambda component: (len(component), component)))


def highest_roots(gram, simple, positive):
    components = connected_components(simple*gram*simple.transpose())
    inverse_simple = simple.pseudoinverse()
    result = []
    for component in components:
        candidates = []
        for root in positive:
            coordinates = vector(QQ, root)*inverse_simple
            if not all(value in ZZ and value >= 0 for value in coordinates):
                continue
            support = tuple(index for index, value in enumerate(coordinates) if value)
            if support and all(index in component for index in support):
                candidates.append((sum(coordinates), root, coordinates))
        _, root, coordinates = max(candidates, key=lambda item: item[0])
        result.append((component, vector(ZZ, root), vector(ZZ, coordinates)))
    return result


def intersection(left, right, ns):
    return ZZ(vector(ZZ, left)*ns*vector(ZZ, right))


def chamber_reduce(divisor, curves, ns):
    divisor = vector(ZZ, divisor)
    sequence = []
    while True:
        for name, curve in curves:
            pairing = intersection(divisor, curve, ns)
            if pairing < 0:
                divisor += pairing*curve
                sequence.append((name, pairing))
                assert divisor*ns*divisor == 0
                break
        else:
            return divisor, tuple(sequence)


with (DATA / "kumar_q80_to_rootless_path.tsv").open() as handle:
    steps = list(csv.DictReader(handle, delimiter="\t"))
start = load_matrix(DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt")
first = steps[0]
first_frame, first_transport = neighbor(
    start,
    ZZ(first["q"]),
    ZZ(first["a"]),
    ZZ(first["b"]),
    vector(ZZ, map(ZZ, first["v"].split(","))),
)
assert first_frame.det() == 948

ns = block_diagonal_matrix(U, -first_frame)
simple, positive = deterministic_simple_roots(first_frame)
assert simple.nrows() == 13
cartan = simple*first_frame*simple.transpose()
components = connected_components(cartan)
assert sorted(map(len, components)) == [4, 9]

fiber = vector(ZZ, [1, 0]+[0]*17)
zero = vector(ZZ, [-1, 1]+[0]*17)
curves = [("O", zero)]
for index, root in enumerate(simple.rows(), 1):
    curves.append((f"R{index}", vector(ZZ, [0, 0]+list(root))))
affine_data = highest_roots(first_frame, simple, positive)
for component_index, (_, root, _) in enumerate(affine_data, 1):
    curves.append(
        (f"Theta0_{component_index}", fiber-vector(ZZ, [0, 0]+list(root)))
    )

second = steps[1]
raw = vector(
    ZZ,
    [ZZ(second["a"]), ZZ(second["b"])]
    + list(map(ZZ, second["v"].split(","))),
)
assert raw*ns*raw == 0
reduced, sequence = chamber_reduce(raw, curves, ns)
assert intersection(reduced, fiber, ns) == 2

frame_part = vector(QQ, reduced[2:])
root_gram = simple*first_frame*simple.transpose()
root_coordinates = frame_part*first_frame*simple.transpose()*root_gram.inverse()
root_projection = root_coordinates*simple
mw_projection = frame_part-root_projection
mw_norm = mw_projection*first_frame*mw_projection

# The projected image of the full integral frame is the saturated MW height
# lattice.  Compute it directly rather than relying on a version-dependent
# LLL basis or enumerating the full rank-17 shell.
projection = (
    identity_matrix(QQ, 17)
    - first_frame*simple.transpose()*root_gram.inverse()*simple
)
projection_denominator = lcm(value.denominator() for value in projection.list())
scaled_projection = (projection_denominator*projection).change_ring(ZZ)
mw_projected_integer = scaled_projection.row_module().basis_matrix()
assert mw_projected_integer.nrows() == 4
mw_basis = mw_projected_integer/projection_denominator
mw_height = mw_basis*first_frame*mw_basis.transpose()
assert mw_height.det() == QQ(237)/5
scaled_height = (20*mw_height).change_ring(ZZ)
lll_columns = scaled_height.LLL_gram()
mw_transform = lll_columns.transpose()
assert abs(mw_transform.det()) == 1
mw_basis = mw_transform*mw_basis
mw_height = mw_basis*first_frame*mw_basis.transpose()
assert 20*mw_height == matrix(
    ZZ,
    [
        (19, 5, 0, -7),
        (5, 35, 0, -5),
        (0, 0, 80, 40),
        (-7, -5, 40, 171),
    ],
)
mw_form = QuadraticForm(ZZ, (40*mw_height).change_ring(ZZ))
mw_shells = mw_form.short_vector_list_up_to_length(40, up_to_sign_flag=True)
low_mw = []
for shell in mw_shells[1:]:
    for row in shell:
        row = vector(ZZ, row)
        low_mw.extend((row, -row))
low_mw = tuple(sorted(set(map(tuple, low_mw))))


def integral_preimage(projected):
    """Return x in ZZ^17 with x*projection=projected."""
    target = vector(QQ, projected)*projection_denominator
    assert all(value in ZZ for value in target)
    target = vector(ZZ, target)
    linear_map = scaled_projection.transpose()
    diagonal, left, right = linear_map.smith_form()
    assert left*linear_map*right == diagonal
    transformed = left*target
    smith_coordinates = vector(ZZ, [0]*17)
    for index in range(17):
        value = diagonal[index, index]
        if value:
            assert transformed[index] % value == 0
            smith_coordinates[index] = transformed[index]//value
        else:
            assert transformed[index] == 0
    preimage = right*smith_coordinates
    assert preimage*scaled_projection == target
    assert vector(QQ, preimage)*projection == vector(QQ, projected)
    return vector(ZZ, preimage)


root_smith = simple.smith_form()[0]
assert tuple(abs(root_smith[index, index]) for index in range(13)) == (1,)*13
root_lattice = IntegralLattice(root_gram)


def shortest_lift(projected):
    raw_lift = integral_preimage(projected)
    root_coordinates = (
        vector(QQ, raw_lift)*first_frame*simple.transpose()*root_gram.inverse()
    )
    shift = vector(ZZ, next(root_lattice.enumerate_close_vectors(-root_coordinates)))
    lift = raw_lift+shift*simple
    assert vector(QQ, lift)*projection == vector(QQ, projected)
    return vector(ZZ, lift)


low_sections = []
for coordinates_tuple in low_mw:
    coordinates = vector(ZZ, coordinates_tuple)
    projected = coordinates*mw_basis
    height_value = projected*first_frame*projected
    assert height_value < 2
    lift = shortest_lift(projected)
    norm = ZZ(lift*first_frame*lift)
    assert norm >= 4 and norm % 2 == 0
    pole = norm//2-2
    section = vector(ZZ, [pole+1, 1]+list(lift))
    assert section*ns*section == -2
    pairing = intersection(reduced, section, ns)
    low_sections.append(
        (tuple(coordinates), height_value, norm-height_value, norm, pole, pairing, tuple(lift))
    )
assert all(item[5] >= 0 for item in low_sections)
assert vector(ZZ, reduced[2:])*first_frame*vector(ZZ, reduced[2:]) == 8

print(
    f"Q80SECONDCHAMBER|components={tuple(map(len, components))}|"
    f"cartans={tuple(tuple(tuple(cartan[i,j] for j in component) for i in component) for component in components)}|"
    f"affine_multiplicities={tuple(tuple(coordinates[index] for index in component) for component, _, coordinates in affine_data)}",
    flush=True,
)
print(f"Q80SECONDCHAMBER|MW_height_lt_2={low_mw}", flush=True)
for item in low_sections:
    coordinates, height_value, correction, norm, pole, pairing, lift = item
    print(
        f"Q80SECONDCHAMBER|low_section_MW={coordinates}|height={height_value}|"
        f"root_correction={correction}|frame_norm={norm}|P.O={pole}|"
        f"D.S={pairing}|lift={lift}",
        flush=True,
    )
print(
    "Q80SECONDCHAMBER|all_sections_nonnegative=1|"
    "root_primitive=1|MW_torsion=1",
    flush=True,
)
print(
    f"Q80SECONDCHAMBER|raw={tuple(raw)}|D.F={intersection(raw, fiber, ns)}|"
    f"D.O={intersection(raw, zero, ns)}",
    flush=True,
)
print(f"Q80SECONDCHAMBER|reduction={sequence}", flush=True)
print(
    f"Q80SECONDCHAMBER|after_fibers={tuple(reduced)}|"
    f"D.O={intersection(reduced, zero, ns)}|"
    f"pairings={tuple((name, intersection(reduced, curve, ns)) for name, curve in curves)}",
    flush=True,
)
print(
    f"Q80SECONDCHAMBER|root_coordinates={tuple(root_coordinates)}|"
    f"mw_projection={tuple(mw_projection)}|mw_norm={mw_norm}",
    flush=True,
)
print(
    f"Q80SECONDCHAMBER|MW_basis={tuple(tuple(row) for row in mw_basis.rows())}|"
    f"height={tuple(tuple(row) for row in mw_height.rows())}|"
    f"height_det={mw_height.det()}",
    flush=True,
)
# If a negative irreducible (-2)-curve C had degree two over the old base,
# D-C would be vertical, so C would have zero MW projection.  Root
# primitivity gives C=kF+2O+r.  With D=4F+2O+rho and rho^2=8,
#
#   ||r-rho||^2 = 2(D.C+1).
#
# Thus D.C<0 forces r=rho and D.C=-1, after which C^2=-2 gives k=7/2,
# impossible.  Degree-one curves are the sections checked above.
assert QQ(8+6)/4 == QQ(7)/2 and QQ(7)/2 not in ZZ
print(
    "Q80SECONDCHAMBER|negative_bisection_impossible=1|"
    "all_effective_minus2_nonnegative=1|nef=1|status=PASS_NEF_CHAMBER",
    flush=True,
)


# The next pinned step is qualitatively different.  Reconstruct the second
# child, reduce the q=12 witness in its D7+D5 chamber, and record its MW part.
second_frame, second_transport = neighbor(
    first_frame,
    ZZ(second["q"]),
    ZZ(second["a"]),
    ZZ(second["b"]),
    vector(ZZ, map(ZZ, second["v"].split(","))),
)
assert second_frame.det() == 948
third_ns = block_diagonal_matrix(U, -second_frame)
third_simple, third_positive = deterministic_simple_roots(second_frame)
assert third_simple.nrows() == 12
third_cartan = third_simple*second_frame*third_simple.transpose()
third_components = connected_components(third_cartan)
assert sorted(map(len, third_components)) == [5, 7]
third_fiber = vector(ZZ, [1, 0]+[0]*17)
third_zero = vector(ZZ, [-1, 1]+[0]*17)
third_curves = [("O", third_zero)]
for index, root in enumerate(third_simple.rows(), 1):
    third_curves.append((f"R{index}", vector(ZZ, [0, 0]+list(root))))
third_affine_data = highest_roots(second_frame, third_simple, third_positive)
for component_index, (_, root, _) in enumerate(third_affine_data, 1):
    third_curves.append(
        (f"Theta0_{component_index}", third_fiber-vector(ZZ, [0, 0]+list(root)))
    )

third = steps[2]
third_raw = vector(
    ZZ,
    [ZZ(third["a"]), ZZ(third["b"])]
    + list(map(ZZ, third["v"].split(","))),
)
assert tuple(third_raw) == (
    3, 4, -2, 1, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
)
third_reduced, third_sequence = chamber_reduce(third_raw, third_curves, third_ns)
assert tuple(third_reduced) == (
    3, 3, -18, -20, 18, 0, -24, -22, 9, 6, 5, 42, -71, -10, 6, 16,
    19, -8, -8,
)
assert intersection(third_reduced, third_fiber, third_ns) == 3
assert intersection(third_reduced, third_zero, third_ns) == 0
third_frame_part = vector(QQ, third_reduced[2:])
third_root_gram = third_simple*second_frame*third_simple.transpose()
third_root_coordinates = (
    third_frame_part*second_frame*third_simple.transpose()*third_root_gram.inverse()
)
assert tuple(third_root_coordinates) == (
    -6, -3, -3, -3, -2, -5, -1, -1, -2, -1, -2, -4,
)
third_mw_projection = third_frame_part-third_root_coordinates*third_simple
third_mw_norm = third_mw_projection*second_frame*third_mw_projection
assert third_mw_projection != 0 and third_mw_norm == 8
third_root_lattice = IntegralLattice(third_root_gram)
third_root_shift = vector(
    ZZ,
    next(third_root_lattice.enumerate_close_vectors(-third_root_coordinates)),
)
third_short_lift = vector(ZZ, third_frame_part)+third_root_shift*third_simple
assert vector(QQ, third_short_lift)-third_mw_projection in third_simple.row_module(QQ)
third_short_norm = ZZ(third_short_lift*second_frame*third_short_lift)
assert third_short_norm >= third_mw_norm and third_short_norm % 2 == 0
third_section_pole = third_short_norm//2-2
third_section = vector(
    ZZ, [third_section_pole+1, 1]+list(third_short_lift)
)
assert third_section*third_ns*third_section == -2
third_vertical_root = third_root_coordinates*third_simple
third_clean_decomposition = (
    third_section
    + 2*third_zero
    + 2*third_fiber
    + vector(ZZ, [0, 0]+list(third_vertical_root))
)
assert third_clean_decomposition == third_reduced
third_divisor_section_pairing = intersection(
    third_reduced, third_section, third_ns
)
third_section_component_pairings = tuple(
    (name, intersection(third_section, curve, third_ns))
    for name, curve in third_curves[1:]
)
assert all(pairing >= 0 for _, pairing in third_section_component_pairings)
print(
    f"Q80THIRDCHAMBER|raw={tuple(third_raw)}|D.F=4|D.O=-1|"
    f"reduction={third_sequence}",
    flush=True,
)
print(
    f"Q80THIRDCHAMBER|after_fibers={tuple(third_reduced)}|D.F=3|D.O=0|"
    f"root_coordinates={tuple(third_root_coordinates)}|"
    f"mw_projection={tuple(third_mw_projection)}|mw_norm={third_mw_norm}|"
    f"short_lift={tuple(third_short_lift)}|short_norm={third_short_norm}|"
    f"section_P.O={third_section_pole}|D.S={third_divisor_section_pairing}|"
    "decomposition=S+2O+2F+root_correction|"
    "equation_gate=MW_marked_trisection|status=PASS_CHAMBER",
    flush=True,
)
print(
    f"Q80THIRDCHAMBER|S_component_pairings={third_section_component_pairings}",
    flush=True,
)
third_vertical_components = tuple(
    (component, tuple(-third_root_coordinates[index] for index in component))
    for component in third_components
)
third_affine_multiplicities = tuple(
    (component, tuple(coordinates[index] for index in component))
    for component, _, coordinates in third_affine_data
)
print(
    f"Q80THIRDCHAMBER|vertical_components={third_vertical_components}|"
    f"affine_multiplicities={third_affine_multiplicities}",
    flush=True,
)


# Pin the best q=8 presentation found by the bounded q=4,...,12 comparison.
# It is a genuine alternate neighbor, but its root rank is thirteen and its
# MW rank only four, so it does not improve the rank-growing q=12 continuation.
alternate_q8 = vector(
    ZZ, (0, 0, -1, 0, 0, 0, -2, -1, -2, 0, -2, 0, 0, -1, 0, -1, 0)
)
alternate_raw = vector(ZZ, [2, 4]+list(alternate_q8))
alternate_reduced, alternate_sequence = chamber_reduce(
    alternate_raw, third_curves, third_ns
)
assert intersection(alternate_reduced, third_fiber, third_ns) == 2
assert intersection(alternate_reduced, third_zero, third_ns) == 0
assert all(
    intersection(alternate_reduced, curve, third_ns) >= 0
    for _, curve in third_curves
)
alternate_component_pairings = tuple(
    (name, intersection(alternate_reduced, curve, third_ns))
    for name, curve in third_curves[1:]
)
alternate_frame_part = vector(QQ, alternate_reduced[2:])
alternate_source_root_coordinates = (
    alternate_frame_part*second_frame*third_simple.transpose()
    * third_root_gram.inverse()
)
alternate_source_mw_projection = (
    alternate_frame_part-alternate_source_root_coordinates*third_simple
)
alternate_source_mw_norm = (
    alternate_source_mw_projection*second_frame*alternate_source_mw_projection
)

# A negative section is automatically very short.  For a section with frame
# lift l, D.C=||l||^2-<d,l>; Cauchy--Schwarz and ||d||^2=8 show that D.C<0
# forces ||l||^2<8.  Thus norms four and six are exhaustive.  This qfminim
# call is deliberately capped at six; unlike the obsolete norm-12 full-frame
# shell it is small on this D7+D5 frame.
assert vector(ZZ, alternate_reduced[2:])*second_frame*vector(
    ZZ, alternate_reduced[2:]
) == 8
alternate_short = pari(second_frame).qfminim(6)
alternate_short_vectors = matrix(ZZ, alternate_short[2]).transpose()
alternate_section_pairings = []
for lift in alternate_short_vectors.rows():
    norm = ZZ(lift*second_frame*lift)
    if norm not in (4, 6):
        continue
    for signed_lift in (vector(ZZ, lift), -vector(ZZ, lift)):
        pole = norm//2-2
        section = vector(ZZ, [pole+1, 1]+list(signed_lift))
        assert section*third_ns*section == -2
        alternate_section_pairings.append(
            intersection(alternate_reduced, section, third_ns)
        )
assert alternate_section_pairings
assert min(alternate_section_pairings) >= 0
# If a negative irreducible bisection C=(k,2,l) existed, the identity
# ||l-d||^2=2(D.C+1) would force l=d and D.C=-1.  Then C^2=-2 and d^2=8
# would give 8=4k+2, i.e. k=3/2, impossible in the integral NS basis.
assert QQ(8-2)/4 == QQ(3)/2 and QQ(3)/2 not in ZZ
alternate_frame, alternate_transport = neighbor(
    second_frame, 8, 2, 4, alternate_q8
)
alternate_roots = roots_of_norm_two(alternate_frame)
alternate_root_module = matrix(ZZ, [list(root) for root in alternate_roots]).row_module()
alternate_root_basis = alternate_root_module.basis_matrix()
alternate_root_rank = alternate_root_basis.nrows()
alternate_root_det = abs(
    (alternate_root_basis*alternate_frame*alternate_root_basis.transpose()).det()
)
assert (
    alternate_root_rank,
    len(alternate_roots),
    alternate_root_det,
    17-alternate_root_rank,
) == (13, 128, 24, 4)

# Recover a deterministic saturated height Gram for comparison with the
# source D7+D5/MW5 frame.  The projected image of the full integral frame is
# already the saturated Mordell--Weil lattice because the root lattice is
# primitive (and hence the fibration has no torsion).
alternate_root_gram = (
    alternate_root_basis*alternate_frame*alternate_root_basis.transpose()
)
alternate_root_smith = alternate_root_basis.smith_form()[0]
assert tuple(
    abs(alternate_root_smith[index, index])
    for index in range(alternate_root_basis.nrows())
) == (1,)*13
alternate_projection = (
    identity_matrix(QQ, 17)
    - alternate_frame*alternate_root_basis.transpose()
    * alternate_root_gram.inverse()*alternate_root_basis
)
alternate_denominator = lcm(
    value.denominator() for value in alternate_projection.list()
)
alternate_scaled_projection = (
    alternate_denominator*alternate_projection
).change_ring(ZZ)
alternate_mw_integer = (
    alternate_scaled_projection.row_module().basis_matrix()
)
assert alternate_mw_integer.nrows() == 4
alternate_mw_basis = alternate_mw_integer/alternate_denominator
alternate_mw_height = (
    alternate_mw_basis*alternate_frame*alternate_mw_basis.transpose()
)
assert alternate_mw_height.det() == QQ(79)/2
alternate_height_scale = lcm(
    value.denominator() for value in alternate_mw_height.list()
)
alternate_lll_columns = (
    alternate_height_scale*alternate_mw_height
).change_ring(ZZ).LLL_gram()
alternate_mw_change = alternate_lll_columns.transpose()
assert abs(alternate_mw_change.det()) == 1
alternate_mw_height = (
    alternate_mw_change*alternate_mw_height*alternate_mw_change.transpose()
)
print(
    f"Q80SECONDCHILDQ8|v={tuple(alternate_q8)}|ab=2,4|"
    "root_rank=13|roots=128|rootdet=24|ADE=E6+A7|MW=4|"
    f"reduction={alternate_sequence}|nef_against_O_and_components=1|"
    f"reduced={tuple(alternate_reduced)}|D.F=2|D.O=0|"
    f"zero={tuple(third_zero)}|component_pairings={alternate_component_pairings}|"
    f"source_MW_projection={tuple(alternate_source_mw_projection)}|"
    f"source_MW_norm={alternate_source_mw_norm}|"
    f"short_sections={len(alternate_section_pairings)}|"
    f"minimum_section_pairing={min(alternate_section_pairings)}|"
    "negative_bisection_impossible=1|nef=1|"
    f"root_primitive=1|torsion=1|height={tuple(tuple(row) for row in alternate_mw_height.rows())}|"
    "source=D7+D5/MW5|child=E6+A7/MW4|"
    "rho=19|trivial_rank=15|shioda_tate_MW=4|"
    "same_frame=0|status=PASS_ALTERNATE_HIGH_ROOT",
    flush=True,
)
