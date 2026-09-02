#!/usr/bin/env sage-python
"""Certify the E6+II rational-surface rank-sum-three construction.

status: ACTIVE_PROOF
claim: one-parameter E6 quadratic base changes with exact rank split 2+1,
       generic Picard rank 19, and no rootless MW17 fibration in the same NS
inputs: none
outputs: elkies-k3-e6-ii-rank3-quadratic-base-change-v1.json

The rational elliptic surface has fibres IV*+II+2I1 and two marked sections.
The quadratic cover is ramified at the II fibre and one smooth fibre.  Its K3
pullback has roots 2E6+A2 and one additional anti-invariant section.  A
Blichfeldt--Hermite bound rules out every rootless rank-17 frame of determinant
24, so the requested same-NS rootless search has an exact negative answer.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import (
    CartanMatrix,
    EllipticCurve,
    Infinity,
    PolynomialRing,
    QQ,
    RealBallField,
    ZZ,
    ceil,
    gamma,
    matrix,
    vector,
    zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-e6-ii-rank3-quadratic-base-change-v1.json"
)


def compact_height(curve, point, multiplier, chi):
    """Recover a height after clearing every fibre component group."""

    coordinate = (multiplier * point)[0]
    numerator = coordinate.numerator()
    denominator = coordinate.denominator()
    if not denominator.is_square():
        raise ArithmeticError("x-coordinate denominator is not a square")
    finite_poles = denominator.degree() // 2
    infinity_poles = ceil((numerator.degree() - 2 * chi) / 2)
    intersection = max(finite_poles, infinity_poles)
    return {
        "multiplier": int(multiplier),
        "x_numerator_degree": int(numerator.degree()),
        "x_denominator_degree": int(denominator.degree()),
        "intersection_with_zero": int(intersection),
        "height": str(QQ(2 * chi + 2 * intersection) / multiplier**2),
    }


def height_gram(curve, points, multiplier, chi):
    diagonal = [
        QQ(compact_height(curve, point, multiplier, chi)["height"])
        for point in points
    ]
    gram = matrix(QQ, len(points))
    for index, value in enumerate(diagonal):
        gram[index, index] = value
    for left in range(len(points)):
        for right in range(left):
            height_sum = QQ(
                compact_height(
                    curve, points[left] + points[right], multiplier, chi
                )["height"]
            )
            pairing = (height_sum - diagonal[left] - diagonal[right]) / 2
            gram[left, right] = gram[right, left] = pairing
    return gram


# The one rational parameter is r.  Work first on the rational surface over
# QQ(r)(u).
R_RING = PolynomialRing(QQ, "r")
r = R_RING.gen()
K = R_RING.fraction_field()
U_RING = PolynomialRing(K, "u")
u = U_RING.gen()
FU = U_RING.fraction_field()

a = K(2 * (r**2 + r + 1) / (r + 1))
c = K(-2 * r**2 / (r + 1))
A_surface = a * u
B_surface = u * (u + c)
surface = EllipticCurve(FU, [A_surface, B_surface])
P = surface(K(1), u + 1)
Q = surface(K(r**2), u + r**3)

if P[1] ** 2 != P[0] ** 3 + A_surface * P[0] + B_surface:
    raise ArithmeticError("first invariant section identity failed")
if Q[1] ** 2 != Q[0] ** 3 + A_surface * Q[0] + B_surface:
    raise ArithmeticError("second invariant section identity failed")

surface_discriminant_reduced = 4 * A_surface**3 + 27 * B_surface**2
surface_residual = 4 * a**3 * u + 27 * (u + c) ** 2
if surface_discriminant_reduced != u**2 * surface_residual:
    raise ArithmeticError("IV*+II+2I1 discriminant factorization failed")

surface_residual_discriminant = K(surface_residual.discriminant())
expected_surface_residual_discriminant = K(
    1024
    * (r - 1) ** 2
    * (2 * r + 1) ** 2
    * (r + 2) ** 2
    * (r**2 + r + 1) ** 3
    / (4 * (r + 1) ** 6)
)
if surface_residual_discriminant != expected_surface_residual_discriminant:
    raise ArithmeticError("residual I1 discriminant changed")

surface_gram = height_gram(surface, (P, Q), 3, 1)
expected_surface_gram = matrix(QQ, ((QQ(2) / 3, -QQ(1) / 3),
                                    (-QQ(1) / 3, QQ(2) / 3)))
if surface_gram != expected_surface_gram:
    raise ArithmeticError("rational-surface height Gram changed")


# The twist d*y^2=f has the constant section (0,1).  The cover
# u=-c/(1-t^2) has d=(c*t/(1-t^2))^2.
d = u * (u + c)
if d != B_surface:
    raise ArithmeticError("twist section (0,1) identity failed")

T_RING = PolynomialRing(K, "t")
t = T_RING.gen()
FT = T_RING.fraction_field()
H = 1 - t**2
A_k3 = -a * c * H**3
B_k3 = c**2 * t**2 * H**4
k3 = EllipticCurve(FT, [A_k3, B_k3])
P_k3 = k3(H**2, H**2 * (H - c))
Q_k3 = k3(r**2 * H**2, H**2 * (r**3 * H - c))
S_k3 = k3(0, c * t * H**2)

for label, point in (("P", P_k3), ("Q", Q_k3), ("S", S_k3)):
    if point[1] ** 2 != point[0] ** 3 + A_k3 * point[0] + B_k3:
        raise ArithmeticError(f"K3 section {label} failed")
if P_k3[0](-t) != P_k3[0] or P_k3[1](-t) != P_k3[1]:
    raise ArithmeticError("P is not deck invariant")
if Q_k3[0](-t) != Q_k3[0] or Q_k3[1](-t) != Q_k3[1]:
    raise ArithmeticError("Q is not deck invariant")
if S_k3[0](-t) != S_k3[0] or S_k3[1](-t) != -S_k3[1]:
    raise ArithmeticError("S is not deck anti-invariant")

k3_gram = height_gram(k3, (P_k3, Q_k3, S_k3), 3, 2)
expected_k3_gram = matrix(
    QQ,
    (
        (QQ(4) / 3, -QQ(2) / 3, 0),
        (-QQ(2) / 3, QQ(4) / 3, 0),
        (0, 0, QQ(2) / 3),
    ),
)
if k3_gram != expected_k3_gram or k3_gram.det() != QQ(8) / 9:
    raise ArithmeticError("K3 height Gram changed")

k3_discriminant_reduced = 4 * A_k3**3 + 27 * B_k3**2
k3_residual = 27 * c * t**4 - 4 * a**3 * H
if k3_discriminant_reduced != c**3 * H**8 * k3_residual:
    raise ArithmeticError("K3 discriminant factorization failed")

configuration_invariant = K(-4 * a**3 / (27 * c))
expected_configuration_invariant = K(
    16 * (r**2 + r + 1) ** 3 / (27 * r**2 * (r + 1) ** 2)
)
if configuration_invariant != expected_configuration_invariant:
    raise ArithmeticError("marked configuration invariant changed")
if configuration_invariant.derivative() == 0:
    raise ArithmeticError("K3 family became isotrivial")


# Integral NS in the basis O,F,E6a[1..6],E6b[1..6],A2[1..2],P,Q,S.
# P and Q meet component 1 in both IV* fibres.  S meets component 1 in
# the first IV*, component 6 in the second, and component 1 in the IV fibre.
# The three section curves are pairwise disjoint.
ns = zero_matrix(ZZ, 19)
ns[0, 0] = -2
ns[0, 1] = ns[1, 0] = 1
ns[2:8, 2:8] = -CartanMatrix(["E", 6])
ns[8:14, 8:14] = -CartanMatrix(["E", 6])
ns[14:16, 14:16] = -CartanMatrix(["A", 2])
P_index, Q_index, S_index = 16, 17, 18
for section_index in (P_index, Q_index, S_index):
    ns[section_index, section_index] = -2
    ns[1, section_index] = ns[section_index, 1] = 1
for section_index in (P_index, Q_index):
    ns[section_index, 2] = ns[2, section_index] = 1
    ns[section_index, 8] = ns[8, section_index] = 1
for component_index in (2, 13, 14):
    ns[S_index, component_index] = ns[component_index, S_index] = 1

if ns.det() != 24:
    raise ArithmeticError("Neron--Severi determinant changed")
if ns.elementary_divisors()[-3:] != [2, 2, 6]:
    raise ArithmeticError("Neron--Severi Smith invariants changed")


# The only possible proper index in this determinant-24 lattice is two.
# Eliminate all seven projective mod-two section combinations on one good
# specialization.  Also eliminate generic 3-primary torsion: every generic
# torsion section injects into the component groups (exponent three), while
# the selected good fibre has torsion Z/2.
r_control = QQ(2)
t_control = QQ(2)
a_control = QQ(a(r_control))
c_control = QQ(c(r_control))
H_control = 1 - t_control**2
fibre = EllipticCurve(
    QQ,
    [
        -a_control * c_control * H_control**3,
        c_control**2 * t_control**2 * H_control**4,
    ],
)
control_points = (
    fibre(H_control**2, H_control**2 * (H_control - c_control)),
    fibre(
        r_control**2 * H_control**2,
        H_control**2 * (r_control**3 * H_control - c_control),
    ),
    fibre(0, c_control * t_control * H_control**2),
)
divisibility = {}
for mask in range(1, 8):
    point = fibre(0)
    labels = []
    for index, label in enumerate(("P", "Q", "S")):
        if (mask >> index) & 1:
            point += control_points[index]
            labels.append(label)
    key = "+".join(labels)
    divisible = bool(point.is_divisible_by(2))
    divisibility[key] = {
        "point": [str(point[0]), str(point[1])],
        "divisible_by_2": divisible,
    }
if any(row["divisible_by_2"] for row in divisibility.values()):
    raise ArithmeticError("possible index-two MW saturation survived")
if fibre.torsion_subgroup().order() != 2:
    raise ArithmeticError("control-fibre torsion changed")


# Same-NS rootless gate.  A rootless rank-17 frame would be even of minimum
# at least four and determinant 24.  Its Hermite invariant would exceed the
# Blichfeldt upper bound for gamma_17.
RBF = RealBallField(200)
n = ZZ(17)
blichfeldt = (RBF(2) / RBF.pi()) * gamma(RBF(2) + RBF(n) / 2) ** (RBF(2) / n)
required_hermite = RBF(4) / RBF(24) ** (RBF(1) / n)
determinant_lower_bound = (RBF(4) / blichfeldt) ** n
if not (blichfeldt < required_hermite and determinant_lower_bound > 24):
    raise ArithmeticError("Blichfeldt rootless obstruction did not separate")


payload = {
    "schema": "elkies-k3.e6-ii-rank3-quadratic-base-change.v1",
    "status": "PASS_EXACT_E6_II_RANK_SUM_3_RHO19_ROOTLESS_IMPOSSIBLE",
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_e6_ii_rank3_quadratic_base_change.sage"
    ),
    "parameter": {
        "name": "r",
        "base_field": "QQ(r)",
        "rational_open_exclusions": ["-2", "-1", "-1/2", "0", "1"],
    },
    "rational_elliptic_surface": {
        "equation": "y^2=x^3+a*u*x+u*(u+c)",
        "a": str(a),
        "c": str(c),
        "fibre_profile": "IV*+II+2I1",
        "root_lattice": "E6",
        "root_rank": 6,
        "sections": [
            {"label": "P", "x": "1", "y": "u+1"},
            {"label": "Q", "x": "r^2", "y": "u+r^3"},
        ],
        "height_gram": [[str(value) for value in row] for row in surface_gram.rows()],
        "mordell_weil_rank_exact": 2,
        "torsion": "trivial",
        "residual_I1_discriminant": str(surface_residual_discriminant),
    },
    "quadratic_twist": {
        "squareclass": "d(u)=u*(u+c)",
        "model": "d(u)*y^2=x^3+a*u*x+u*(u+c)",
        "section": {"x": "0", "y": "1"},
        "rank_exact": 1,
    },
    "quadratic_base_change": {
        "map": "u=-c/(1-t^2)",
        "square_identity": "u*(u+c)=(c*t/(1-t^2))^2",
        "branch_values": ["u=0 (II)", "u=-c (smooth)"],
    },
    "k3": {
        "equation": "Y^2=X^3-a*c*(1-t^2)^3*X+c^2*t^2*(1-t^2)^4",
        "fibre_profile": "2IV*+IV+4I1",
        "root_lattice": "2E6+A2",
        "root_rank": 14,
        "sections": [
            {"label": "P", "character": "invariant", "X": "H^2", "Y": "H^2*(H-c)"},
            {"label": "Q", "character": "invariant", "X": "r^2*H^2", "Y": "H^2*(r^3*H-c)"},
            {"label": "S", "character": "anti-invariant", "X": "0", "Y": "c*t*H^2"},
        ],
        "height_gram": [[str(value) for value in row] for row in k3_gram.rows()],
        "height_determinant": str(k3_gram.det()),
        "rank_sum": 3,
        "rank_decomposition": "2 invariant + 1 anti-invariant",
        "mordell_weil_rank_exact_generic": 3,
        "picard_rank_exact_generic": 19,
        "marked_configuration_invariant": str(configuration_invariant),
        "neron_severi_determinant_absolute": 24,
        "neron_severi_smith_nontrivial": [2, 2, 6],
        "neron_severi_gram": [[int(value) for value in row] for row in ns.rows()],
        "saturation_control": {
            "r": 2,
            "t": 2,
            "curve_a_invariants": [str(value) for value in fibre.a_invariants()],
            "curve_discriminant": str(fibre.discriminant()),
            "specialized_torsion_order": int(fibre.torsion_subgroup().order()),
            "mod_two_combinations": divisibility,
        },
    },
    "same_ns_rootless_mw17": {
        "exists": False,
        "proof": (
            "A rootless Jacobian fibration would have an even positive-definite "
            "rank-17 frame of determinant 24 and minimum at least 4.  Its Hermite "
            "invariant exceeds Blichfeldt's upper bound for gamma_17."
        ),
        "rank": 17,
        "determinant": 24,
        "minimum_required": 4,
        "blichfeldt_upper_gamma17": str(blichfeldt),
        "required_hermite_invariant": str(required_hermite),
        "forced_determinant_lower_bound": str(determinant_lower_bound),
    },
    "proof_boundary": {
        "proved": (
            "The section identities, discriminant and fibre profiles, generic height "
            "matrices, rank split 2+1, generic Picard rank 19, trivial torsion, "
            "MW/NS saturation, determinant 24, and same-NS rootless impossibility."
        ),
        "not_proved": (
            "No rank-sum-four family is asserted here; no claim is made about exact "
            "ranks at every rational specialization or about rootless fibrations on "
            "other Neron--Severi lattices."
        ),
    },
}


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
output_path = arguments.output.resolve()
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "E6IIR3|surface=IVstar+II+2I1_MW2|twist=MW1|"
    "k3=2IVstar+IV+4I1_MW3_rho19_disc24|rootlessMW17=IMPOSSIBLE|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
