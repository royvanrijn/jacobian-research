#!/usr/bin/env sage
"""Transport the marked CM-43 q=8, q=60, and D9+E7 fibers to geometry.

The output basis is

    [F,O,E7 simple components(7),E8 simple components(8),P1,P2,P3].

The section intersections and component labels come from the explicit
CM-43 Kumar equation.  Splitting off U=<F,O+F> gives a positive frame
integrally isometric to the pinned glue-211 frame.  Matching the explicit
height-four section P3 fixes the remaining global sign and makes the three
fiber transports canonical up to Weyl relabeling of the displayed simple
components.
"""

from pathlib import Path

from sage.all import *


BASE = Path(__file__).resolve().parents[1]


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


# Geometric divisor Gram.  P1 and P2 meet the minuscule E7 component (node 7
# in Sage's ordering), while P3 and all three sections at E8 meet the identity
# component.  Exact section intersections are P1.P2=1, P1.P3=3, P2.P3=2.
geometric_ns = matrix(ZZ, 20, 20, 0)
geometric_ns[0, 1] = geometric_ns[1, 0] = 1
geometric_ns[1, 1] = -2
geometric_ns[2:9, 2:9] = -CartanMatrix(["E", 7])
geometric_ns[9:17, 9:17] = -CartanMatrix(["E", 8])
for section in (17, 18, 19):
    geometric_ns[section, section] = -2
    geometric_ns[0, section] = geometric_ns[section, 0] = 1
geometric_ns[17, 18] = geometric_ns[18, 17] = 1
geometric_ns[17, 19] = geometric_ns[19, 17] = 3
geometric_ns[18, 19] = geometric_ns[19, 18] = 2
for section in (17, 18):
    geometric_ns[8, section] = geometric_ns[section, 8] = 1
assert geometric_ns.det() == -43

# Split the old fibration's U and recover its positive frame.
old_fiber = vector(ZZ, [1]+[0]*19)
old_isotropic_mate = vector(ZZ, [1, 1]+[0]*18)  # O+F
orthogonal = matrix(
    ZZ,
    [
        list(old_fiber*geometric_ns),
        list(old_isotropic_mate*geometric_ns),
    ],
).right_kernel_matrix()
geometric_split = matrix(
    ZZ,
    [list(old_fiber), list(old_isotropic_mate)]
    + [list(row) for row in orthogonal.rows()],
)
assert geometric_split.det() == 1
split_gram = geometric_split*geometric_ns*geometric_split.transpose()
U = matrix(ZZ, ((0, 1), (1, 0)))
assert split_gram[:2, :2] == U
assert split_gram[:2, 2:].is_zero() and split_gram[2:, :2].is_zero()
geometric_frame = -split_gram[2:, 2:]
assert geometric_frame.det() == 43

marked_frame = load_gram(
    BASE / "data/fibrations/kumar_cm43_marked_e7e8_mw3_frame.txt"
)
frame_isometry = matrix(
    ZZ, pari(marked_frame).qfisom(pari(geometric_frame))
)
# PARI returns C with C^t*target*C=source.
assert frame_isometry.transpose()*geometric_frame*frame_isometry == marked_frame
assert abs(frame_isometry.det()) == 1

# Fix the sign by requiring the pinned height-four frame vector to be the
# geometric P3 divisor rather than its inverse section.
marked_p3_frame = vector(ZZ, [0]*15+[4, 0, -1])
geometric_p3 = vector(ZZ, [0]*19+[1])
geometric_p3_split = geometric_p3*geometric_split.inverse()
assert tuple(geometric_p3_split[:2]) == (1, 1)
if -marked_p3_frame*frame_isometry.transpose() == geometric_p3_split[2:]:
    frame_isometry = -frame_isometry
assert marked_p3_frame*frame_isometry.transpose() == geometric_p3_split[2:]


def to_geometric(marked_divisor):
    marked_divisor = vector(ZZ, marked_divisor)
    split_coordinates = vector(
        ZZ,
        list(marked_divisor[:2])
        + list(marked_divisor[2:]*frame_isometry.transpose()),
    )
    answer = vector(ZZ, split_coordinates*geometric_split)
    assert answer*geometric_ns*answer == 0
    return answer


old_fiber_class = vector(ZZ, [1]+[0]*19)
zero_section_class = vector(ZZ, [0, 1]+[0]*18)


q8_witness = vector(ZZ, (
    156, -78, 0, 0, -78, 0, -78, 0, 0,
    0, 0, 0, 0, 0, 0, -1, -155, -32,
))
q8_geometric = to_geometric(vector(ZZ, [2, 4]+list(q8_witness)))
assert q8_geometric == vector(ZZ, (
    8, 5,
    -1, -2, -2, -3, -2, -2, -1,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, -2, 0,
))

# The same fixed-component issue already occurs for q=8.  Its raw old-fiber
# degree is b=4, but subtracting 2O and the ensuing 14 E7 components leaves
# the nef a=2 pencil.
q8_nef = vector(ZZ, q8_geometric)
q8_reflection_sequence = []
q8_reduction_curves = [vector(ZZ, [0, 1]+[0]*18)]
q8_reduction_curves.extend(
    vector(ZZ, [0, 0]+[1 if column == row else 0 for column in range(18)])
    for row in range(7)
)
while True:
    pairings = [q8_nef*geometric_ns*curve for curve in q8_reduction_curves]
    negative = [index for index, pairing in enumerate(pairings) if pairing < 0]
    if not negative:
        break
    index = negative[0]
    pairing = ZZ(pairings[index])
    q8_reflection_sequence.append((index, pairing))
    q8_nef += pairing*q8_reduction_curves[index]
assert q8_reflection_sequence == [
    (0, -2),
    (5, -1), (4, -1), (3, -1), (1, -1), (7, -1), (6, -1),
    (5, -1), (4, -1), (2, -1), (3, -1), (4, -1), (5, -1),
    (6, -1), (7, -1),
]
assert q8_nef == vector(ZZ, (
    8, 3,
    -2, -3, -4, -6, -5, -4, -3,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, -2, 0,
))
q8_section = vector(ZZ, (
    10, 2,
    -2, -3, -4, -6, -5, -4, -3,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, -2, 0,
))
assert q8_section*geometric_ns*q8_section == -2
assert q8_nef == q8_section+zero_section_class-2*old_fiber_class
assert q8_nef*geometric_ns*old_fiber_class == 2

# At Picard rank 20 this is still not nef: the two CM-only sections below
# finish the chamber reduction and expose the old fiber itself.
inverse_p2_section = vector(ZZ, (
    4, 2,
    -2, -3, -4, -6, -5, -4, -3,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, -1, 0,
))
p1_minus_p2_section = vector(ZZ, (
    3, 1,
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, -1, 0,
))
assert inverse_p2_section*geometric_ns*inverse_p2_section == -2
assert p1_minus_p2_section*geometric_ns*p1_minus_p2_section == -2
assert q8_nef*geometric_ns*inverse_p2_section == -1
q8_after_inverse_p2 = q8_nef-inverse_p2_section
assert q8_after_inverse_p2 == old_fiber_class+p1_minus_p2_section
assert q8_after_inverse_p2*geometric_ns*p1_minus_p2_section == -1
q8_final = q8_after_inverse_p2-p1_minus_p2_section
assert q8_final == old_fiber_class
q8_initial_fixed = q8_geometric-q8_nef
assert q8_geometric == (
    old_fiber_class+q8_initial_fixed
    + inverse_p2_section+p1_minus_p2_section
)

q60_marked = vector(ZZ, (
    5, 12,
    0, 0, -1, -1, -1, -1, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0,
))
q60_geometric = to_geometric(q60_marked)
assert q60_geometric == vector(ZZ, (
    17, 12,
    0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    4, -5, 1,
))

# The transported vector is isotropic but not yet in the effective chamber:
# it pairs -7 with O and -1 with the minuscule E7 component.  Deterministic
# reflection/fixed-component reduction against [O,E7 simples] takes 28 steps.
displayed_curves = [vector(ZZ, [0, 1]+[0]*18)]
displayed_curves.extend(
    vector(ZZ, [0, 0]+[1 if column == row else 0 for column in range(18)])
    for row in range(7)
)
q60_nef = vector(ZZ, q60_geometric)
reflection_sequence = []
while True:
    pairings = [q60_nef*geometric_ns*curve for curve in displayed_curves]
    negative = [index for index, pairing in enumerate(pairings) if pairing < 0]
    if not negative:
        break
    index = negative[0]
    pairing = ZZ(pairings[index])
    reflection_sequence.append((index, pairing))
    q60_nef += pairing*displayed_curves[index]
assert len(reflection_sequence) == 28
assert reflection_sequence == [
    (0, -7), (7, -1), (6, -1), (5, -1), (4, -1), (2, -1),
    (3, -1), (1, -1), (4, -1), (3, -1), (5, -1), (4, -1),
    (2, -1), (6, -1), (5, -1), (4, -1), (3, -1), (1, -1),
    (7, -1), (6, -1), (5, -1), (4, -1), (2, -1), (3, -1),
    (4, -1), (5, -1), (6, -1), (7, -1),
]
assert q60_nef == vector(ZZ, (
    17, 5,
    -2, -3, -4, -6, -5, -4, -3,
    0, 0, 0, 0, 0, 0, 0, 0,
    4, -5, 1,
))
assert tuple(q60_nef*geometric_ns) == (
    5, 7, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
    4, 30, 17,
)

# The actual section Q79 has O-intersection 58 and the nonidentity E7 label.
# Its divisor class follows either from the Shioda map or from its exact
# intersections with the displayed basis.
q79_section = vector(ZZ, (
    60, 1,
    -2, -3, -4, -6, -5, -4, -3,
    0, 0, 0, 0, 0, 0, 0, 0,
    4, -5, 1,
))
assert q79_section*geometric_ns*q79_section == -2
assert q79_section*geometric_ns*zero_section_class == 58
assert q60_nef == q79_section+4*zero_section_class-43*old_fiber_class
assert q60_nef*geometric_ns*old_fiber_class == 5

# The remaining negative curves are CM-only sections.  They collapse the
# q=60 specialization all the way to the old Kumar fiber, explaining the
# E7+E8 return and why the complete CM RR space is larger than a pencil.
p3_minus_p2_section = vector(ZZ, (
    4, 1,
    -2, -3, -4, -6, -5, -4, -3,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, -1, 1,
))
assert p3_minus_p2_section*geometric_ns*p3_minus_p2_section == -2
assert q60_nef*geometric_ns*p1_minus_p2_section == -4
q60_after_p1_minus_p2 = q60_nef-4*p1_minus_p2_section
assert q60_after_p1_minus_p2 == old_fiber_class+p3_minus_p2_section
assert q60_after_p1_minus_p2*geometric_ns*p3_minus_p2_section == -1
q60_final = q60_after_p1_minus_p2-p3_minus_p2_section
assert q60_final == old_fiber_class
q60_initial_fixed = q60_geometric-q60_nef
assert q60_initial_fixed == vector(ZZ, (
    0, 7,
    2, 3, 4, 6, 5, 4, 3,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0,
))
assert q60_geometric == (
    old_fiber_class+q60_initial_fixed
    + 4*p1_minus_p2_section+p3_minus_p2_section
)

source_d9e7_witness = vector(ZZ, (
    0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 1, 1, 1, 2, 3, 4, 0, -1,
))
source_d9e7_geometric = to_geometric(
    vector(ZZ, [2, 4]+list(source_d9e7_witness))
)
assert source_d9e7_geometric == vector(ZZ, (
    4, 3,
    0, 0, 0, 0, 0, 0, 0,
    1, 2, 3, 4, 3, 5, 3, 1,
    0, 0, 1,
))

print(
    f"CM43DIVISOR|basis=F,O,E7x7,E8x8,P1,P2,P3"
    f"|q8_raw={tuple(q8_geometric)}|raw_old_degree=4",
    flush=True,
)
print(
    f"CM43DIVISOR|q8_nef={tuple(q8_nef)}|old_degree=2"
    f"|reflections={len(q8_reflection_sequence)}|identity=R+O-2F",
    flush=True,
)
print(
    "CM43DIVISOR|q8_CM_collapse=F+initial_fixed+(-P2)+(P1-P2)"
    "|full_chamber=old_F",
    flush=True,
)
print(
    f"CM43DIVISOR|q60_raw={tuple(q60_geometric)}"
    f"|raw_old_degree={q60_geometric*geometric_ns*old_fiber_class}",
    flush=True,
)
print(
    f"CM43DIVISOR|q60_nef={tuple(q60_nef)}|old_degree=5"
    f"|reflections={len(reflection_sequence)}"
    f"|identity=Q79+4O-43F",
    flush=True,
)
print(
    "CM43DIVISOR|q60_CM_collapse=F+initial_fixed+4(P1-P2)+(P3-P2)"
    "|full_chamber=old_F",
    flush=True,
)
print(f"CM43DIVISOR|source_D9E7={tuple(source_d9e7_geometric)}", flush=True)
print("CM43DIVISOR|status=PASS", flush=True)
