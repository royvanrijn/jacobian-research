#!/usr/bin/env sage -python
"""Derive the exact generic-fibre ambient for the first H3 q=8 pencil.

The selected q=8 fibre class is recorded by the exhaustive lattice transport
``classify_h3_q6_child_q8_orbits.sage``.  Pulled back to the original H3
E7+E8/MW2 frame it differs from

    9 O + 9 (-P1)

only by a vertical divisor.  Thus its restriction to the old generic fibre is
the degree-18 divisor ``9(O)+9(-P1)``.  With

    m = (y-y(P1))/(x-x(P1)),

the elliptic function field is quadratic over ``QQ(t)(m)``.  This script
proves that the 18 functions

    1,m,...,m^9, x,x*m,...,x*m^7

form the complete generic-fibre Riemann--Roch space.  It intentionally does
not guess the vertical/resolution conditions: deriving those finite quotient
blocks is the remaining equation-level q=8 task.
"""

import argparse
import hashlib
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1 = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
Q8_ORBITS = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json"


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_frame(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
exec(compile(CORE.read_text(), str(CORE), "exec"))

q8 = json.loads(Q8_ORBITS.read_text())["q8"]
hits = q8["d13_mw4_hits"]
assert len(hits) == 2

# This is the representative independently used in
# analyze_h3_first_q6_chamber.sage: it has q6 MW projection (0,-2,0).
hit = next(hit for hit in hits if hit["mw_projection"] == [0, -2, 0])
raw_divisor = vector(ZZ, hit["fiber_source_h3_ns"])
assert len(raw_divisor) == 19

fiber = vector(ZZ, [1, 0] + [0] * 17)
zero = vector(ZZ, [-1, 1] + [0] * 17)
twice_minuscule = (2, 3, 4, 6, 5, 4, 3)
minus_p1 = vector(
    ZZ, [5, 1] + [-value for value in twice_minuscule] + [0] * 8 + [1, 0]
)

# The q=8 classifier works in a nef chamber of the q=6 child.  Its raw
# pullback has fixed old-fibre components on the original H3 equation.  Strip
# them in the old E7+E8 chamber before using the class as a source linear
# system.  The generic-fibre restriction is preserved, but the vertical
# divisor is not; this rules out extrapolating a q=6 local module by power.
source_ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -load_frame(FRAME))
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)
divisor = vector(ZZ, raw_divisor)
fixed_component_reflections = []
while True:
    negative = tuple(
        (node, ZZ(divisor * source_ns * curve))
        for node, curve in enumerate(simple)
        if divisor * source_ns * curve < 0
    )
    if not negative:
        break
    node, pairing = negative[0]
    fixed_component_reflections.append((node + 1, int(pairing)))
    divisor += pairing * simple[node]

assert len(fixed_component_reflections) == 122
assert tuple(divisor * source_ns * curve for curve in simple) == (
    0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0,
)
highest_e7 = (2, 2, 3, 4, 3, 2, 1)
highest_e8 = (2, 3, 4, 6, 5, 4, 3, 2)
affine_e7 = fiber - sum(
    (coefficient * simple[index] for index, coefficient in enumerate(highest_e7)),
    vector(ZZ, [0] * 19),
)
affine_e8 = fiber - sum(
    (coefficient * simple[7 + index] for index, coefficient in enumerate(highest_e8)),
    vector(ZZ, [0] * 19),
)
assert (divisor * source_ns * affine_e7, divisor * source_ns * affine_e8) == (11, 16)
assert divisor * source_ns * divisor == 0
assert divisor * source_ns * fiber == 18

# In the frame [U,E7,E8,P1,P2], the source-nef divisor differs from
# 9O+9(-P1) by a vertical class.  This is the exact target for a resolved
# E7/E8 quotient calculation.
vertical_difference = divisor - 9 * zero - 9 * minus_p1
assert vertical_difference[1] == 0
assert tuple(vertical_difference[-2:]) == (0, 0)
assert tuple(vertical_difference) == (
    -11, 0, 2, 3, 4, 6, 5, 5, 6,
    -4, -5, -7, -10, -8, -6, -4, -2, 0, 0,
)
# Retain the natural repeated-section representative exactly.  It is linearly
# equivalent on the generic fibre to a one-marked-point divisor, but that
# conversion would conceal the q8 chord-power basis and must not replace this
# literal divisor before the resolved vertical conditions are imposed.
generic_support = certify_generic_fibre_horizontal_support(
    source_ns,
    divisor,
    fiber,
    (("O", 9, zero), ("-P1", 9, minus_p1)),
    tuple(
        ("source_simple_{}".format(index + 1), int(vertical_difference[index + 2]), simple[index])
        for index in range(len(simple))
        if vertical_difference[index + 2]
    ),
    fiber_twist=int(vertical_difference[0]),
    expected_old_fiber_degree=18,
)
assert generic_support["reconstructed_divisor"] == tuple(divisor)
assert [(item["name"], item["multiplicity"]) for item in generic_support["horizontal_support"]] == [
    ("O", 9), ("-P1", 9),
]

section = json.loads(P1.read_text())
assert section["status"] == "PASS_EXACT_H92_P1"
assert section["exact_weierstrass_square"]

anchor = SourceFileLoader("h92_q8_ambient_anchor", str(ANCHOR)).load_module()
r, s = anchor.EXPECTED_H92
_, formulas = anchor.parse_h92(H92)
A1, A, B1, B, B2 = tuple(QQ(value(r, s)) for value in formulas)

base_ring = PolynomialRing(QQ, "u")
u0 = base_ring.gen()
base_field = base_ring.fraction_field()
u = base_field.gen()
T = 1 / u
old_a = A1 * T**3 + A * T**4
old_b = B1 * T**5 + B * T**6 + B2 * T**7

x_p = base_field(polynomial(base_ring, section["x_entrance_base"]["numerator_coefficients"]))
x_p /= base_field(polynomial(base_ring, section["x_entrance_base"]["denominator_coefficients"]))
y_p = base_field(polynomial(base_ring, section["y_entrance_base"]["numerator_coefficients"]))
y_p /= base_field(polynomial(base_ring, section["y_entrance_base"]["denominator_coefficients"]))
assert y_p**2 == x_p**3 + old_a * x_p + old_b
assert y_p

# Substitute y=y(P1)+m(x-x(P1)) into the Weierstrass equation.  The marked
# root x=x(P1) factors out; the remaining factor is monic quadratic in x.
chord_ring = PolynomialRing(base_field, "m")
m = chord_ring.gen()
relation_ring = PolynomialRing(chord_ring, "x")
x = relation_ring.gen()
relation = (
    (relation_ring(y_p) + relation_ring(m) * (x - relation_ring(x_p)))**2
    - x**3 - relation_ring(old_a) * x - relation_ring(old_b)
)
quadratic, remainder = relation.quo_rem(x - relation_ring(x_p))
assert not remainder
assert quadratic.degree(x) == 2
assert quadratic.leading_coefficient() == -1
quadratic = -quadratic
assert quadratic.leading_coefficient() == 1

# At O, m=y/x+regular has a simple pole and x has pole order two.  At -P1,
# x-x(P1) is a uniformizer because y(P1) is nonzero, so m has a simple pole
# while x is regular.  The following list is therefore contained in
# L(9O+9(-P1)); it has the required degree 18 and is independent because
# {1,x} is a QQ(t)(m)-basis by the monic quadratic above.
basis = [
    {
        "kind": entry["kind"],
        "m_power": entry["m_power"],
        "x_power": entry["x_power"],
        "pole_order_at_O": entry["pole_order_at_O"],
        "pole_order_at_minus_P1": entry["pole_order_at_marked_section"],
    }
    for entry in balanced_marked_chord_power_basis(9, "m", "x")
]
assert len(basis) == 18
assert max(row["pole_order_at_O"] for row in basis) == 9
assert max(row["pole_order_at_minus_P1"] for row in basis) == 9

payload = {
    "schema": "elkies-k3.h92-q8-generic-rr-ambient.v1",
    "status": "PASS_EXACT_Q8_GENERIC_RR_AMBIENT",
    "scope": (
        "Exact generic-fibre Riemann--Roch ambient only. The remaining "
        "vertical and resolved-chart conditions needed to cut this space to "
        "a q=8 pencil are not asserted."
    ),
    "source_q8_lattice_class": list(map(int, divisor)),
    "source_q8_lattice_selection": {
        "q6_child_mw_projection": hit["mw_projection"],
        "raw_source_class_in_q6_child_chamber": list(map(int, raw_divisor)),
        "old_fibre_degree_before_fixed_components": int(raw_divisor[1]),
        "source_fixed_component_reflections": [list(entry) for entry in fixed_component_reflections],
        "source_nef_old_fibre_degree": int(divisor[1]),
        "source_nef_simple_component_pairings": [
            int(divisor * source_ns * curve) for curve in simple
        ],
        "source_nef_affine_component_pairings": [
            int(divisor * source_ns * affine_e7),
            int(divisor * source_ns * affine_e8),
        ],
        "vertical_difference_D_minus_9O_minus_9minusP1": list(map(int, vertical_difference)),
    },
    "generic_fibre_support_certificate": generic_support,
    "generic_fibre_divisor": "9*O + 9*(-P1)",
    "generic_fibre_degree": 18,
    "chord": "m=(y-y(P1))/(x-x(P1))",
    "quadratic_relation": str(quadratic),
    "quadratic_extension": "QQ(t)(x,y)=QQ(t)(m)[x]/(quadratic_relation)",
    "basis": basis,
    "dimension": 18,
    "riemann_roch": "genus_one_degree_18_implies_h0=18",
    "next_required_step": (
        "Derive finite vertical and resolved E7/E8 quotient conditions for "
        "this fixed basis; certify that their exact common kernel has dimension 2."
    ),
    "inputs": {
        "compiler_core": {"path": str(CORE.relative_to(ROOT)), "sha256": digest(CORE)},
    },
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q8AMBIENT|degree=18|basis=18|quadratic_in_x=1|"
    "vertical_difference=PASS|status=PASS_EXACT_Q8_GENERIC_RR_AMBIENT",
    flush=True,
)
