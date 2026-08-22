#!/usr/bin/env sage -python
"""Audit the H92 q6-child q8 rational-point/MW marking bridge.

The q8 marking artifact constructs an exact rational section S on the globally
minimal E8+E6 elliptic K3 and assigns it MW coordinates (-2,-2,0) in the
certified height lattice

    [[8/3,1/3,-1],
     [1/3,8/3,1],
     [-1,1,46]].

Independently, the exact rational functions for S have x-denominator h^2 and
y-denominator h^3, where h is squarefree of degree 46 and coprime to the
Weierstrass discriminant.  The formal zero parameter z=-x/y has one zero at
each h-point, so S.O=46 on smooth fibres.

For an elliptic K3, Shioda's self-height formula is

    <S,S> = 4 + 2(S.O) - sum_v contr_v(S).

The only reducible fibres are E8 and E6.  Bounding each local correction by
the largest diagonal entry of the corresponding inverse Cartan matrix gives
an intentionally loose but sufficient upper bound 30+6=36.  Hence the actual
rational section has height at least 60, while (-2,-2,0) has lattice height
24.  The point-to-MW bridge is therefore inconsistent.

This script is a regression guard for that known contradiction.  It succeeds
only when the contradiction is reproduced; it does not repair the marking.
"""

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
MARKING = GEN / "elkies-k3-h92-q6-child-q8-marking.json"


def polynomial(ring, values):
    return ring([QQ(value) for value in values])


def monic_power_root(value, exponent):
    root = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        root *= factor.monic() ** (multiplicity // exponent)
    return root.monic()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--marking", type=Path, default=MARKING)
args = parser.parse_args()

child = json.loads(args.child.read_text())
marking = json.loads(args.marking.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"

ring = PolynomialRing(QQ, "T")
field = ring.fraction_field()
A = polynomial(ring, child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B = polynomial(ring, child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
Delta = polynomial(ring, child["minimal_short_weierstrass"]["Delta_coefficients_low_to_high"])
assert A.degree() <= 8 and B.degree() <= 12 and Delta.degree() <= 24
curve = EllipticCurve(field, [0, 0, 0, field(A), field(B)])

section = marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
sx = field(polynomial(ring, section["x_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, section["x_denominator_coefficients_low_to_high"])
)
sy = field(polynomial(ring, section["y_numerator_coefficients_low_to_high"])) / field(
    polynomial(ring, section["y_denominator_coefficients_low_to_high"])
)
S = curve(sx, sy)
assert S in curve

h_x = monic_power_root(ring(sx.denominator()), 2)
h_y = monic_power_root(ring(sy.denominator()), 3)
assert h_x == h_y
h = h_x
assert h.degree() == 46 and h.gcd(Delta) in QQ
z_over_h = (-sx / sy) / field(h)
assert ring(z_over_h.numerator()).gcd(h) in QQ
assert ring(z_over_h.denominator()).gcd(h) in QQ
smooth_O_intersections = int(h.degree())

height_gram = matrix(QQ, [
    [QQ(8) / 3, QQ(1) / 3, -1],
    [QQ(1) / 3, QQ(8) / 3, 1],
    [-1, 1, 46],
])
coords = vector(QQ, marking["selected_q8"]["relative_child_section_MW_coordinates"])
assert tuple(coords) == (-2, -2, 0)
claimed_height = QQ(coords * height_gram * coords)
assert claimed_height == 24

E8 = matrix(QQ, [
    [2,0,-1,0,0,0,0,0],
    [0,2,0,-1,0,0,0,0],
    [-1,0,2,-1,0,0,0,0],
    [0,-1,-1,2,-1,0,0,0],
    [0,0,0,-1,2,-1,0,0],
    [0,0,0,0,-1,2,-1,0],
    [0,0,0,0,0,-1,2,-1],
    [0,0,0,0,0,0,-1,2],
])
E6 = matrix(QQ, [
    [2,-1,0,0,0,0],
    [-1,2,-1,0,0,0],
    [0,-1,2,-1,0,-1],
    [0,0,-1,2,-1,0],
    [0,0,0,-1,2,0],
    [0,0,-1,0,0,2],
])
max_E8 = max(E8.inverse().diagonal())
max_E6 = max(E6.inverse().diagonal())
max_total_correction = max_E8 + max_E6
height_lower_bound = QQ(4 + 2 * smooth_O_intersections) - max_total_correction
required_correction = QQ(4 + 2 * smooth_O_intersections) - claimed_height
contradiction = claimed_height < height_lower_bound
assert contradiction

print(
    "Q8MARKHEIGHT|"
    f"minimal_A_degree={A.degree()}|minimal_B_degree={B.degree()}|"
    f"delta_degree={Delta.degree()}|smooth_O_intersections={smooth_O_intersections}|"
    f"claimed_coords={','.join(map(str, coords))}|claimed_height={claimed_height}",
    flush=True,
)
print(
    "Q8MARKHEIGHT_BOUND|"
    f"max_E8_correction={max_E8}|max_E6_correction={max_E6}|"
    f"max_total_correction={max_total_correction}|"
    f"required_correction_for_claim={required_correction}|"
    f"height_lower_bound={height_lower_bound}",
    flush=True,
)
print(
    "Q8MARKHEIGHT_RESULT|contradiction=1|claimed_height=24|"
    f"lower_bound={height_lower_bound}|gap={height_lower_bound-claimed_height}|"
    "status=PASS_CONFIRMED_MARKING_BRIDGE_CONTRADICTION",
    flush=True,
)
