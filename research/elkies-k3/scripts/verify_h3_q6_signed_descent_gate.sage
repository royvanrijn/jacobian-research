#!/usr/bin/env sage -python
"""Verify the lattice gate for descending the two signed H3 q=6 pencils.

The marked q=6 divisor in ``analyze_h3_first_q6_chamber.sage`` is

    D_minus = O + (-P1) - F.

This checker constructs the opposite section from the same integral frame,
rather than changing only the Mordell--Weil coordinate and accidentally
retaining the wrong E7 correction.  It then computes the intersection of the
two signed pencils and the size of their natural trace/norm linear systems.

The result is conditional at exactly one point: an equation-level comparison
has not yet proved that conjugation of the oriented H21 cover exchanges these
two section classes.  If it does, the natural descent is a degree-21 map to a
rational surface.  It is not a rational elliptic pencil.
"""

from sage.all import *

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
FRAME_SHA256 = "ba09ec834a7229e11e4ca687d187f663b6368c3e2fac9b5133bb1570e7031599"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-h3-q6-signed-descent-gate.json"
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


def coordinates(value):
    return [int(entry) for entry in value]


parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

assert digest(FRAME) == FRAME_SHA256
frame = load_gram(FRAME)
assert frame.nrows() == 17 and frame.det() == 948
NS = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -frame)

F = vector(ZZ, [1, 0] + [0] * 17)
O = vector(ZZ, [-1, 1] + [0] * 17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)

# In the integral frame, -P1 has MW coordinate +1 and the twice-minuscule
# E7 correction.  The inverse section has MW coordinate -1, but its integral
# E7 representative is zero.  Both meet the same nonidentity E7 component.
twice_minuscule = (2, 3, 4, 6, 5, 4, 3)
minus_P1 = vector(
    ZZ,
    [5, 1]
    + [-value for value in twice_minuscule]
    + [0] * 8
    + [1, 0],
)
plus_P1 = vector(ZZ, [5, 1] + [0] * 15 + [-1, 0])

expected_component_pairings = (0, 0, 0, 0, 0, 0, 1) + (0,) * 8
for section in (minus_P1, plus_P1):
    assert section * NS * section == -2
    assert section * NS * F == 1
    assert section * NS * O == 4
    assert tuple(section * NS * curve for curve in simple) == (
        expected_component_pairings
    )

# Certify that the opposite E7 correction is forced by these component
# pairings once the MW coordinate is -1.
root_gram = frame[:7, :7]
mw_cross = vector(QQ, list(frame.row(15)[:7]))
target = vector(QQ, [0, 0, 0, 0, 0, 0, -1])
opposite_root = (target + mw_cross) * root_gram.inverse()
assert opposite_root == 0

D_minus = O + minus_P1 - F
D_plus = O + plus_P1 - F
for divisor in (D_minus, D_plus):
    assert divisor * NS * divisor == 0
    assert divisor * NS * F == 2
    assert divisor * NS * O == 1
    assert all(divisor * NS * curve >= 0 for curve in simple)

signed_intersection = ZZ(D_minus * NS * D_plus)
H = D_minus + D_plus
assert signed_intersection == 21
assert H * NS * H == 42
assert H * NS * F == 4
assert gcd(tuple(H)) == 1

# If both signed divisors are nef, their sum is nef and big.  K3
# Riemann--Roch/Kawamata--Viehweg gives h0(H)=2+H^2/2=23.  The two pencils are
# independent because their intersection is positive.  Hence their product
# map to P1 x P1 is generically finite of degree 21, and the pullback of
# H0(O(1,1)) is a basepoint-free rank-four module.  The canonical unordered
# trace/norm coordinates form rank three, not rank two.  A member of |H| has
# arithmetic genus 1+H^2/2=22.
h0_H = ZZ(2 + (H * NS * H) / 2)
arithmetic_genus_H = ZZ(1 + (H * NS * H) / 2)
assert h0_H == 23 and arithmetic_genus_H == 22

payload = {
    "schema": "elkies-k3.h3-q6-signed-descent-gate.v1",
    "status": "PASS_CONDITIONAL_SIGNED_DESCENT_GATE",
    "input": {
        "frame": str(FRAME.relative_to(ROOT)),
        "sha256": FRAME_SHA256,
    },
    "classes": {
        "minus_P1": coordinates(minus_P1),
        "plus_P1": coordinates(plus_P1),
        "D_minus": coordinates(D_minus),
        "D_plus": coordinates(D_plus),
        "H": coordinates(H),
    },
    "intersections": {
        "D_minus_squared": int(D_minus * NS * D_minus),
        "D_plus_squared": int(D_plus * NS * D_plus),
        "D_minus_D_plus": int(signed_intersection),
        "H_squared": int(H * NS * H),
        "H_old_fiber_degree": int(H * NS * F),
    },
    "linear_system": {
        "H_primitive": True,
        "h0_H_if_nef": int(h0_H),
        "arithmetic_genus_H": int(arithmetic_genus_H),
        "tensor_product_rank": 4,
        "unordered_trace_norm_rank": 3,
        "product_map_degree": int(signed_intersection),
        "natural_descent_is_rank_two": False,
        "natural_descent_is_elliptic_pencil": False,
    },
    "proof_boundary": (
        "The lattice identities are unconditional.  Their descent "
        "interpretation assumes that oriented-cover conjugation exchanges "
        "the section classes P1 and -P1; that marking rule still requires "
        "signed equation-level RR generators.  No claim is made that a "
        "special genus-22 Jacobian cannot have an accidental elliptic factor."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "H3Q6DESCENT|Dminus.Dplus=21|H2=42|H.F=4|H_primitive=1|"
    "h0_H=23|genus_H=22|tensor_rank=4|trace_norm_rank=3|"
    "product_degree=21",
    flush=True,
)
print(
    "H3Q6DESCENT|galois_sign_rule=UNPROVED|"
    "status=PASS_CONDITIONAL_SIGNED_DESCENT_GATE",
    flush=True,
)
