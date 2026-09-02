#!/usr/bin/env sage-python
"""Certify the systematic E6 node--node linear-chord rank-four component.

status: ACTIVE_PROOF
claim: a one-dimensional QQ incidence component has two invariant and two
       independent anti-invariant directions; its saturated NS determinant is 78
inputs: none
outputs: elkies-k3-e6-rank4-linear-chord-incidence-v1.json

The base is a curve over QQ, not a bounded rational-height search.  This
checker does not assert that its genus-zero quotient has been parametrized by
QQ(k), nor does it construct a rootless MW17 neighbour.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import (
    CartanMatrix, Curve, EllipticCurve, GF, PolynomialRing, QQ, RR, ZZ,
    factor, gamma, matrix, pi, vector, zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results"
    / "elkies-k3-e6-rank4-linear-chord-incidence-v1.json"
)
FRAME_PATH = ROOT / "elkies-k3/data/lattice/e6_rank4_det78_frame.txt"


def integer_rows(value):
    return [[int(entry) for entry in row] for row in matrix(ZZ, value).rows()]


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


# ---------------------------------------------------------------------------
# The full split marked E6 rational-surface chart.
# ---------------------------------------------------------------------------
R = PolynomialRing(QQ, names=("v", "w", "M", "ell"))
v, w, M, ell = R.gens()
zv = v**2 + 2
zw = w**2 + 2
mv = v * (v**2 + 3)
mw = w * (w**2 + 3)
a = 2 * (v**2 + v*w + w**2 + 3) / (v+w)
c = -2 * (v**2*w**2 + 2*v**2 - v*w + 2*w**2 + 6) / (v+w)

U = PolynomialRing(R.fraction_field(), "u")
u = U.gen()
A = a*u - 3
B = u**2 + c*u - 2
if (u+mv)**2 != zv**3 + A*zv + B:
    raise ArithmeticError("first marked E6 section identity failed")
if (u+mw)**2 != zw**3 + A*zw + B:
    raise ArithmeticError("second marked E6 section identity failed")

reduced_discriminant = 4*A**3 + 27*B**2
if reduced_discriminant(u=0) != 0 or reduced_discriminant.derivative(u)(u=0) == 0:
    raise ArithmeticError("u=0 is not generically a simple I1 fibre")


# A line through -P_v with slope L=ell*u-v has residual discriminant
#
#   delta = L^4-6*zv*L^2-8*(u+mv)*L-3*zv^2-4*A.
#
# Requiring delta=u*(u-M)*(ell^2*u+q)^2, where
# q=-2*ell*v+M*ell^2/2 is forced by the cubic coefficient, leaves the two
# exact incidence equations below.  The same equations with v,w exchanged
# give the second anti-invariant direction.
node_A = (v+w)*ell**2*M*(ell*M-4*v)**2 - 32*(w**2+3)
node_B = ell*(3*ell**2*M**2 - 8*v*ell*M - 16*v**2 - 48) - 32
node_resultant = R(node_A.resultant(node_B, ell) / (32768*M**6))
if node_resultant == 0:
    raise ArithmeticError("node incidence resultant vanished")

swap = R.hom([w, v, M, ell], R)
swapped_resultant = swap(node_resultant)
difference = node_resultant - swapped_resultant
quotient = R(difference / (-9*(w-v)*(v+w)))
if swap(quotient) != quotient:
    raise ArithmeticError("unordered incidence quotient is not symmetric")


# Convert a symmetric polynomial in v,w to elementary coordinates
# S=v+w and P=v*w.  This keeps the exact component decomposition small.
SP = PolynomialRing(QQ, names=("S", "P", "M"))
S, P, MM = SP.gens()
power_sums = [SP(2), S]
for unused in range(2, 11):
    power_sums.append(S*power_sums[-1] - P*power_sums[-2])


def symmetric_to_elementary(poly):
    coefficients = {tuple(exponent): coefficient for exponent, coefficient in poly.dict().items()}
    answer = SP.zero()
    seen = set()
    for exponent, coefficient in coefficients.items():
        i, j, k, unused_ell = exponent
        if unused_ell:
            raise ArithmeticError("eliminated slope survived")
        if (i, j, k) in seen:
            continue
        if i == j:
            answer += coefficient * P**i * MM**k
            seen.add((i, j, k))
            continue
        partner = (j, i, k, 0)
        if coefficients.get(partner, 0) != coefficient:
            raise ArithmeticError("polynomial is not symmetric")
        answer += coefficient * P**min(i, j) * power_sums[abs(i-j)] * MM**k
        seen.add((i, j, k))
        seen.add((j, i, k))
    return answer


G = symmetric_to_elementary(quotient)
H = symmetric_to_elementary(node_resultant + swapped_resultant)
plane_resultant = G.resultant(H, P)
plane_factors = [entry for entry in plane_resultant.factor() if entry[0] not in (S, MM)]
if len(plane_factors) != 2 or any(exponent != 1 for unused, exponent in plane_factors):
    raise ArithmeticError("unexpected unordered incidence decomposition")

plane_factors = sorted((polynomial for polynomial, unused in plane_factors), key=lambda f: f.total_degree())
Plane = PolynomialRing(QQ, names=("S", "M"))
plane_S, plane_M = Plane.gens()
genus_zero_factor, genus_two_factor = (
    Plane(polynomial.subs({S: plane_S, P: 0, MM: plane_M}))
    for polynomial in plane_factors
)
genus_zero = Curve(genus_zero_factor)
genus_two = Curve(genus_two_factor)
if genus_zero.genus() != 0 or genus_two.genus() != 2:
    raise ArithmeticError("incidence component genera changed")

# The rational singular origin of the genus-zero quotient has a double
# tangent.  After S=M*(X-3/8), its next tangent cone is
# 16384*(256*X^2+9*M^2).  Thus this particular normalization branch splits
# over QQ(i), and no QQ(k)-parameterization is asserted by this checker.
BX = PolynomialRing(QQ, names=("X", "m"))
X, m = BX.gens()
strict = BX(genus_zero_factor.subs({plane_S: m*(X-QQ(3)/8), plane_M: m}) / m**2)
tangent_cone = sum(
    coefficient * X**exponent[0] * m**exponent[1]
    for exponent, coefficient in strict.dict().items()
    if sum(exponent) == 2
)
if tangent_cone != 16384*(256*X**2 + 9*m**2):
    raise ArithmeticError("genus-zero branch tangent changed")


# ---------------------------------------------------------------------------
# A complete good-reduction witness on the genus-zero component.
# ---------------------------------------------------------------------------
F = GF(11)
v0, w0, M0, ell1, ell2 = map(F, (2, 7, 4, 3, 4))
subs1 = {v: ZZ(2), w: ZZ(7), M: ZZ(4), ell: ZZ(3)}
subs2 = {v: ZZ(7), w: ZZ(2), M: ZZ(4), ell: ZZ(4)}
if F(node_A.subs(subs1)) or F(node_B.subs(subs1)):
    raise ArithmeticError("first modular incidence witness failed")
if F(node_A.subs(subs2)) or F(node_B.subs(subs2)):
    raise ArithmeticError("second modular incidence witness failed")
if F(genus_zero_factor(ZZ(9), ZZ(4))) != 0:
    raise ArithmeticError("modular witness is not on the genus-zero component")
if F(genus_two_factor(ZZ(9), ZZ(4))) == 0:
    raise ArithmeticError("modular witness lies on both components")

Rt = PolynomialRing(F, "t")
t = Rt.gen()
K = Rt.fraction_field()
base_H = 1-t**2
base_u = M0/base_H
sqrt_d = M0*t/base_H
a0 = F(2)*(v0**2+v0*w0+w0**2+F(3))/(v0+w0)
c0 = -F(2)*(v0**2*w0**2+F(2)*v0**2-v0*w0+F(2)*w0**2+F(6))/(v0+w0)
k3_A = base_H**3*(a0*M0-F(3)*base_H)
k3_B = base_H**4*(M0**2+c0*M0*base_H-F(2)*base_H**2)
curve = EllipticCurve(K, [0, 0, 0, K(k3_A), K(k3_B)])

invariant = []
anti_invariant = []
chord_sections = []
for parameter, slope in ((v0, ell1), (w0, ell2)):
    section_x = parameter**2 + F(2)
    section_m = parameter*(parameter**2+F(3))
    old_section = curve(
        K(section_x*base_H**2),
        K((M0+section_m*base_H)*base_H**2),
    )
    invariant.append(old_section)

    line = slope*base_u-parameter
    square_factor_constant = -F(2)*slope*parameter + M0*slope**2/F(2)
    square_factor = slope**2*base_u + square_factor_constant
    residual_x = (
        line**2-section_x+square_factor*sqrt_d
    )/F(2)
    residual_y = -(base_u+section_m)+line*(residual_x-section_x)
    chord = curve(K(residual_x*base_H**2), K(residual_y*base_H**3))
    anti = 2*chord-old_section
    chord_sections.append(chord)
    anti_invariant.append(anti)

    if anti[0](t=-t) != anti[0] or anti[1](t=-t) != -anti[1]:
        raise ArithmeticError("anti-invariant deck character failed")


def cleared_height(point):
    """Clear 2E6+A1 component groups and read the K3 canonical height."""
    cleared = 6*point
    x_coordinate = cleared[0]
    pole_degree = max(
        x_coordinate.denominator().degree(),
        x_coordinate.numerator().degree()-4,
    )
    if pole_degree < 0 or pole_degree % 2:
        raise ArithmeticError("invalid cleared pole degree")
    return QQ(4+pole_degree)/36, int(pole_degree)


h1, pole1 = cleared_height(anti_invariant[0])
h2, pole2 = cleared_height(anti_invariant[1])
hsum, pole_sum = cleared_height(anti_invariant[0]+anti_invariant[1])
hdiff, pole_diff = cleared_height(anti_invariant[0]-anti_invariant[1])
pairing = (hsum-h1-h2)/2
anti_gram = matrix(QQ, [[h1, pairing], [pairing, h2]])
if anti_gram != matrix(QQ, [[QQ(22)/3, QQ(4)/3], [QQ(4)/3, QQ(22)/3]]):
    raise ArithmeticError("anti-invariant height Gram changed")
if hdiff != 12 or anti_gram.det() != 52:
    raise ArithmeticError("anti-invariant independence audit failed")

invariant_gram = matrix(QQ, [[QQ(4)/3, -QQ(2)/3], [-QQ(2)/3, QQ(4)/3]])
pure_character_det = invariant_gram.det()*anti_gram.det()
if pure_character_det != QQ(208)/3:
    raise ArithmeticError("pure character determinant changed")

# The chord itself is R_i=(P_i+T_i)/2.  Passing from P,Q,T1,T2 to
# P,Q,R1,R2 therefore has index four.  In this saturated candidate basis the
# height determinant is 13/3 and the 2E6+A1 root determinant is 18.
saturated_gram = matrix(QQ, [
    [QQ(4)/3, -QQ(2)/3, QQ(2)/3, -QQ(1)/3],
    [-QQ(2)/3, QQ(4)/3, -QQ(1)/3, QQ(2)/3],
    [QQ(2)/3, -QQ(1)/3, QQ(13)/6, QQ(1)/6],
    [-QQ(1)/3, QQ(2)/3, QQ(1)/6, QQ(13)/6],
])
if saturated_gram.det() != QQ(13)/3:
    raise ArithmeticError("saturated MW determinant changed")
ns_determinant = ZZ(18*saturated_gram.det())
if ns_determinant != 78 or not ns_determinant.is_squarefree():
    raise ArithmeticError("NS determinant or squarefree saturation gate changed")

# Integral NS marking.  Up to swapping the two IV* fibres and applying the E6
# diagram involution, P and Q meet terminal component 1 in both IV* fibres.
# R1 and R2 meet terminal component 1 in the first IV* fibre, the identity
# component in the second, and the nonidentity I2 component.  The compact
# heights force these profiles: P,Q have correction 8/3, while R1,R2 have
# correction 4/3+1/2=11/6.  The smooth gcds at the good-reduction witness give
# P.R2=Q.R1=1 and all other distinct section intersections zero.
section_intersections = matrix(ZZ, 4, 4)
section_intersections[0, 3] = section_intersections[3, 0] = 1
section_intersections[1, 2] = section_intersections[2, 1] = 1

ns = zero_matrix(ZZ, 19)
ns[0, 0] = -2
ns[0, 1] = ns[1, 0] = 1
ns[2:8, 2:8] = -CartanMatrix(["E", 6])
ns[8:14, 8:14] = -CartanMatrix(["E", 6])
ns[14, 14] = -2
for section_index in range(15, 19):
    ns[section_index, section_index] = -2
    ns[1, section_index] = ns[section_index, 1] = 1
for section_index in (15, 16):
    for component_index in (2, 8):
        ns[section_index, component_index] = 1
        ns[component_index, section_index] = 1
for section_index in (17, 18):
    for component_index in (2, 14):
        ns[section_index, component_index] = 1
        ns[component_index, section_index] = 1
for left in range(4):
    for right in range(left):
        ns[15 + left, 15 + right] = section_intersections[left, right]
        ns[15 + right, 15 + left] = section_intersections[left, right]

if ns.det() != ns_determinant:
    raise ArithmeticError("integral NS determinant changed")
if ns.elementary_divisors() != [1] * 18 + [78]:
    raise ArithmeticError("integral NS Smith invariants changed")

old_fibre = vector(ZZ, [0, 1] + [0] * 17)
old_zero = vector(ZZ, [1] + [0] * 18)
constraints = matrix(ZZ, [ns * old_fibre, ns * (old_zero + old_fibre)])
frame_basis = constraints.right_kernel().basis_matrix()
frame_transport = matrix(ZZ, [
    old_fibre, old_zero + old_fibre, *frame_basis.rows()
])
if abs(frame_transport.det()) != 1:
    raise ArithmeticError("integral old-U split is not unimodular")
split_gram = frame_transport * ns * frame_transport.transpose()
positive_frame = -split_gram[2:, 2:]
if (
    split_gram[:2, :2] != matrix(ZZ, [[0, 1], [1, 0]])
    or any(split_gram[i, j] for i in range(2) for j in range(2, 19))
    or positive_frame.det() != 78
    or not positive_frame.is_positive_definite()
):
    raise ArithmeticError("integral U-plus-frame split changed")
if load_gram(FRAME_PATH) != positive_frame:
    raise ArithmeticError("stored determinant-78 positive frame is stale")

# At t=2 the good specialized elliptic curve has order 17.  Generic torsion
# injects into the component groups Z/3, Z/3 and Z/2, hence has order dividing
# six; specialization eliminates it.
t_control = F(2)
control = EllipticCurve(F, [0, 0, 0, k3_A(t_control), k3_B(t_control)])
if control.discriminant() == 0 or control.cardinality() != 17:
    raise ArithmeticError("torsion control fibre changed")

# Necessary determinant screen for a rootless rank-17 MW lattice.
required_hermite = RR(4)/RR(ns_determinant)**(RR(1)/17)
blichfeldt = RR(2/pi)*RR(gamma(2+RR(17)/2))**(RR(2)/17)
if not required_hermite < blichfeldt:
    raise ArithmeticError("rootless feasibility screen unexpectedly failed")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

payload = {
    "schema": "elkies-k3.e6-rank4-linear-chord-incidence.v1",
    "status": "PASS_EXACT_E6_RANK4_INCIDENCE_CURVE",
    "base_curve": {
        "unordered_plane_factor_genus_zero": str(genus_zero_factor),
        "unordered_plane_factor_genus_two": str(genus_two_factor),
        "geometric_genera": [0, 2],
        "genus_zero_origin_second_tangent_cone": str(tangent_cone),
        "rational_parameterization_status": "not_certified",
    },
    "surface": {
        "equation": "y^2=x^3+(a*u-3)*x+(u^2+c*u-2)",
        "a": str(a),
        "c": str(c),
        "marked_sections": ["(v^2+2,u+v*(v^2+3))", "(w^2+2,u+w*(w^2+3))"],
        "branch_squareclass": "u*(u-M)",
        "k3_fibres_generic": "2IV*+I2+6I1",
        "root_lattice": "2E6+A1",
    },
    "good_reduction_witness": {
        "prime": 11,
        "v_w_M_ell1_ell2": [2, 7, 4, 3, 4],
        "anti_height_gram": [[str(x) for x in row] for row in anti_gram.rows()],
        "anti_height_determinant": str(anti_gram.det()),
        "cleared_pole_degrees_T1_T2_sum_diff": [pole1, pole2, pole_sum, pole_diff],
        "control_t": 2,
        "control_curve_order": int(control.cardinality()),
    },
    "mordell_weil": {
        "rank_split": "2+2",
        "pure_character_height_determinant": str(pure_character_det),
        "saturation_relations": ["2*R1=P+T1", "2*R2=Q+T2"],
        "pure_to_saturated_index": 4,
        "saturated_basis": ["P", "Q", "R1", "R2"],
        "saturated_height_gram": [[str(x) for x in row] for row in saturated_gram.rows()],
        "saturated_height_determinant": str(saturated_gram.det()),
        "generic_torsion": "trivial",
    },
    "neron_severi": {
        "root_determinant": 18,
        "absolute_determinant": int(ns_determinant),
        "generic_picard_rank": 19,
        "basis": [
            "O", "F", *[f"E6a_{index}" for index in range(1, 7)],
            *[f"E6b_{index}" for index in range(1, 7)], "A1_1",
            "P", "Q", "R1", "R2",
        ],
        "integral_gram": integer_rows(ns),
        "smith_invariants": [1] * 18 + [78],
        "section_component_profiles": {
            "P": ["E6a_1", "E6b_1", "I2_identity"],
            "Q": ["E6a_1", "E6b_1", "I2_identity"],
            "R1": ["E6a_1", "E6b_identity", "A1_1"],
            "R2": ["E6a_1", "E6b_identity", "A1_1"],
        },
        "distinct_section_intersection_gram": integer_rows(section_intersections),
        "positive_frame": {
            "path": "elkies-k3/data/lattice/e6_rank4_det78_frame.txt",
            "gram": integer_rows(positive_frame),
            "determinant": int(positive_frame.det()),
        },
    },
    "rootless_mw17_screen": {
        "required_hermite_invariant": float(required_hermite),
        "blichfeldt_upper_bound": float(blichfeldt),
        "determinant_screen": "passes",
        "rootless_frame_constructed": False,
    },
    "proof_boundary": {
        "proved": (
            "A one-dimensional QQ incidence component, four independent arithmetic "
            "directions over its function field, the saturated determinant-78 NS, "
            "and passage of the necessary rootless rank-17 Hermite screen."
        ),
        "not_proved": (
            "A QQ(k)-parameterization or a rational point on the genus-zero quotient; "
            "a rootless MW17 fibration; or completeness beyond the linear-chord incidence."
        ),
    },
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output_path = arguments.output.resolve()
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "E6RANK4|base_genera=0,2|rank_split=2+2|rho=19|"
    "mw_det=13/3|NS_det=78|rootless_screen=PASS|status=PASS_EXACT",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
