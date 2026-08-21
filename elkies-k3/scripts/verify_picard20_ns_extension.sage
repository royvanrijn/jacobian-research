#!/usr/bin/env sage
"""Verify the third section, its glue, and the discriminant-43 NS frame."""

from sage.all import *
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3" / "data" / "fibrations"


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


R = PolynomialRing(QQ, "t")
t = R.gen()
lam = QQ(9) / 25
mu = QQ(49) / 25
pole = QQ(16) / 25
q = t - pole

A = R([0, 0, -QQ(32447500)/583443, -QQ(906250)/194481,
       QQ(31250000)/194481, -QQ(19531250)/194481])
B = R([0, 0, 0, QQ(300827000000)/2315685267, QQ(340001171875)/1029193452,
       -QQ(498857421875)/257298363, QQ(29541015625)/10501974,
       -QQ(152587890625)/85766121, QQ(152587890625)/343064484])
P1X = R([1, -QQ(800)/1323, QQ(625)/147])
P1Y = R([1, -QQ(400)/441, -QQ(394375)/18522, QQ(484375)/9261,
         -QQ(390625)/18522])
N = R([0, QQ(77824)/33075, QQ(12400)/1323, -QQ(30500)/1323, QQ(5000)/441])
M = R([0, 0, QQ(4096)/343, QQ(281408)/9261, -QQ(1517000)/9261,
       QQ(1296875)/6174, -QQ(1015625)/9261, QQ(390625)/18522])

SX = -QQ(625)/441*t**2 + QQ(3550)/1323*t
SY = (QQ(390625)/18522*t**4 - QQ(359375)/9261*t**3
      + QQ(1875)/98*t**2)

assert P1Y**2 == P1X**3 + A*P1X + B
assert M**2 == N**3 + A*N*q**4 + B*q**6
assert SY**2 == SX**3 + A*SX + B

# The new polynomial section is disjoint from P1.  Its raw intersection with
# P2 has roots 0,24/25,32/25; the t=0 contact is at the unresolved I0* point
# and disappears because the sections choose distinct first exceptional
# centers.  Hence P1.S=0 and P2.S=2 on the resolved K3.
p1_gcd = gcd(P1X-SX, P1Y-SY)
p2_gcd = gcd(SX*q**2-N, SY*q**3-M).monic()
assert p1_gcd == 1
assert p2_gcd == t * (t-QQ(24)/25) * (t-QQ(32)/25)


def local_rational(numerator, power, precision=10):
    series_ring = PowerSeriesRing(QQ, "z", default_prec=precision)
    z = series_ring.gen()
    expansion = sum(series_ring(numerator[i])*z**i for i in range(numerator.degree()+1))
    expansion /= (z-pole)**power
    return R([expansion[i] for i in range(precision)])


def additive_path(X, Y):
    P = PolynomialRing(QQ, ("z", "xx", "yy"))
    z, xx, yy = P.gens()

    def embed(poly):
        return sum(P(poly[i])*z**i for i in range(poly.degree()+1))

    surface = yy**2 - xx**3 - embed(A)*xx - embed(B)
    section_x, section_y = R(X), R(Y)
    center_x = center_y = QQ(0)
    path = []
    for _ in range(8):
        section_x = R((section_x-center_x)//t)
        section_y = R((section_y-center_y)//t)
        surface = P(surface(z, center_x+z*xx, center_y+z*yy)//z**2)
        center_x, center_y = section_x(0), section_y(0)
        center = {z: 0, xx: center_x, yy: center_y}
        gradient = tuple(surface.derivative(variable).subs(center) for variable in (z, xx, yy))
        path.append((center_x, center_y, gradient))
        if any(gradient):
            return tuple(path)
    raise RuntimeError("I0* path did not reach a smooth chart")


p2_path = additive_path(local_rational(N, 2), local_rational(M, 3))
s_path = additive_path(SX, SY)
assert p2_path[0][:2] == (QQ(7600)/1323, 0)
assert s_path[0][:2] == (QQ(3550)/1323, 0)
assert p2_path[0][:2] != s_path[0][:2]

# The I3/I3/I2 nodes are avoided, so S meets the identity components there.
node_at_one = QQ(100) / 1323
for point, node in ((QQ(1), node_at_one), (lam, P1X(lam)), (mu, P1X(mu))):
    assert (SX(point), SY(point)) != (node, 0)

# At IV* infinity S and P2 have the same exceptional leading sign, opposite
# to P1.  Together with the distinct nonzero D4 class above, this gives
# S=(E6 class 2, D4 class d2, 0,0,0; S.O=0), up to D4 triality.
assert SX.degree() == 2 and SY.degree() == 4
assert M[7] == SY[4] == -P1Y[4]

# Shioda replay.  Profiles are
# P1=(1,0,0,1,1;0), P2=(2,d1,1,0,0;1), S=(2,d2,0,0,0;0).
local_11 = QQ(4)/3 + QQ(2)/3 + QQ(1)/2
local_22 = QQ(4)/3 + 1 + QQ(2)/3
local_ss = QQ(4)/3 + 1
local_12 = QQ(2)/3
local_1s = QQ(2)/3
local_2s = QQ(4)/3 + QQ(1)/2
H = matrix(QQ, [
    [4-local_11, 2+0+1-2-local_12, 2-0-local_1s],
    [2+0+1-2-local_12, 4+2-local_22, 2+1-2-local_2s],
    [2-0-local_1s, 2+1-2-local_2s, 4-local_ss],
])
expected_height = matrix(QQ, [
    [QQ(3)/2, QQ(1)/3, QQ(4)/3],
    [QQ(1)/3, 3, -QQ(5)/6],
    [QQ(4)/3, -QQ(5)/6, QQ(5)/3],
])
assert H == expected_height and H.det() == QQ(43)/216

# Extend the old explicit NS basis by S.  Choose D4[1] for d2; D4[4] is the
# triality-equivalent choice.  The row records intersections with
# F,O,E6,D4,A2,A2,A1,P1,P2.
U = matrix(ZZ, [[0, 1], [1, 0]])
old_frame = load_matrix(DATA / "mw2_e6_d4_a2a2_a1_frame.txt")
old_basis = load_matrix(DATA / "mw2_e6_d4_a2a2_a1_explicit_basis.txt")
old_ns = block_diagonal_matrix(U, -old_frame)
old_gram = old_basis * old_ns * old_basis.transpose()
s_pairings = [1, 0] + [0]*15 + [0, 2]
s_pairings[2+2] = 1       # E6[3], same E6 class as P2
s_pairings[2+6] = 1       # D4[1], a nonzero class distinct from P2's D4[3]
extended_gram = block_matrix([
    [old_gram, matrix(ZZ, 19, 1, s_pairings)],
    [matrix(ZZ, 1, 19, s_pairings), matrix(ZZ, [[-2]])],
])
assert extended_gram.det() == -43
eigenvalues = extended_gram.change_ring(RR).eigenvalues()
assert sum(value > 0 for value in eigenvalues) == 1
assert sum(value < 0 for value in eigenvalues) == 19
smith = tuple(abs(value) for value in extended_gram.smith_form()[0].diagonal())
assert smith == (1,)*19 + (43,)

# Split off U using G=F+O and orthogonalized section directions.
split = identity_matrix(ZZ, 20)
split[1] = vector(ZZ, [1, 1] + [0]*18)
split[17] = vector(ZZ, [-2, -1] + [0]*15 + [1, 0, 0])
split[18] = vector(ZZ, [-3, -1] + [0]*15 + [0, 1, 0])
split[19] = vector(ZZ, [-2, -1] + [0]*15 + [0, 0, 1])
split_gram = split * extended_gram * split.transpose()
assert split.det() == 1
assert split_gram[:2, :2] == U
assert split_gram[:2, 2:] == 0 and split_gram[2:, :2] == 0
new_frame = -split_gram[2:, 2:]
pinned_frame = load_matrix(DATA / "picard20_e6_d4_a2a2_a1_mw3_frame.txt")
assert new_frame == pinned_frame
assert new_frame.det() == 43 and new_frame.is_positive_definite()
root_data = pari(new_frame).qfminim(2)
assert ZZ(root_data[0]) == 110
root_basis = matrix(ZZ, root_data[2]).transpose().row_module().basis_matrix()
assert root_basis.rank() == 15

# Geometrize the first small neighbor in the optimal Picard-20 path.  In the
# split U+(-frame) basis its q=8 witness is the one replayed independently by
# verify_picard20_mw1_path.sage.  Translate it back to the explicit divisor
# basis F,O,E6(6),D4(4),A2(2),A2(2),A1,P1,P2,S.
q8_vector_part = vector(ZZ, (
    -1, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, -1, 0, 0, -1,
))
q8_split = vector(ZZ, [2, 4] + list(q8_vector_part))
q8_raw = q8_split * split
assert q8_raw == vector(ZZ, (
    8, 5, -1, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, -1, 0, 0, -1,
))
assert q8_raw * extended_gram * q8_raw == 0
assert gcd([abs(ZZ(value)) for value in extended_gram*q8_raw]) == 1

# The raw lattice witness is not in the explicit effective chamber.  Reflect
# it in the listed known effective (-2)-curves.  Indices are one-based in the
# explicit basis; index 1 is O, 2..16 are the displayed fiber components, and
# 17..19 are P1,P2,S.  The recorded second entry is the negative intersection
# immediately before that reflection, so this is a checkable Weyl certificate.
reflection_certificate = (
    (1, -2), (9, -2), (8, -3), (10, -2), (9, -3),
    (3, -1), (4, -1), (5, -1), (6, -2), (4, -1),
    (5, -1), (7, -2), (3, -1), (5, -2), (6, -2),
    (3, -1), (2, -1), (4, -1), (5, -1), (3, -1),
    (7, -1), (5, -1), (6, -1), (10, -1), (11, -1),
)
q8_nef = q8_raw
for one_based_index, expected_pairing in reflection_certificate:
    curve = vector(ZZ, [int(i == one_based_index) for i in range(20)])
    assert curve * extended_gram * curve == -2
    pairing = q8_nef * extended_gram * curve
    assert pairing == expected_pairing
    q8_nef += pairing * curve

assert q8_nef == vector(ZZ, (
    8, 3, -2, -4, -3, -6, -5, -3, -3, -5, -3, -3, 0, 0, 0, 0, -1, 0, 0, -1,
))
assert q8_nef * extended_gram * q8_nef == 0
q8_pairings = tuple(q8_nef * extended_gram)
assert q8_pairings == (2, 2, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 5, 3, 4)
assert all(value >= 0 for value in q8_pairings[1:])

# Restriction to the old generic elliptic fiber discards F and vertical roots:
# it is 3*O-S, of degree two.  Hence the q=8 neighbor is the line pencil
# through S, generated by x-SX and y-SY (up to vertical rescaling/translation).
assert q8_pairings[0] == 2

print("PICARD20NS|section=S|identity=PASS|P1.S=0|P2.S=2|S.O=0", flush=True)
print("PICARD20NS|profile=S:2,d2,0,0,0|D4_triality=up_to_d2_d3", flush=True)
print("PICARD20NS|height_gram=(1/6)*[9,2,8;2,18,-5;8,-5,10]|det=43/216", flush=True)
print("PICARD20NS|rank=20|signature=1,19|disc=43|disc_group=Z/43|saturated=1", flush=True)
print("PICARD20NS|frame_rank=18|frame_det=43|root_rank=15|roots=110|MW=3", flush=True)
print(
    f"PICARD20NS|q8_raw={tuple(q8_raw)}|q8_nef={tuple(q8_nef)}"
    f"|reflections={len(reflection_certificate)}|old_fiber_degree={q8_pairings[0]}"
    "|horizontal=3O-S|pencil=line_through_S",
    flush=True,
)
print("PICARD20NS|status=PASS", flush=True)
