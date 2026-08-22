#!/usr/bin/env sage -python
"""Compile the II* vertical ideal forced by the corrected child q=8 target.

The physical-root certificate puts the q=8 vertical divisor at the II* fibre
on the E8 cycle ``(4,5,7,10,8,6,4,2)``.  Its standard E8 chart order is
``(2,4,6,10,4,7,5,8)``, which is the valuation cycle of ``X`` in the
unit-normalized II* germ

    Y^2 = X^3 + u^4*a(u)*X + u^5*b(u),    a(0)b(0) != 0.

Consequently the *vertical* contribution is the complete ideal
``(u^2,X,Y)`` of colength two.  This is deliberately only an additive-fibre
ideal target: it does not yet determine the coefficient module of the
generic chord basis ``<1,m>`` or assemble a global q=8 pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-physical-root-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-iistar-vertical-ideal.json"

E8_CARTAN = matrix(ZZ, [
    [2, 0, -1, 0, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, 2],
])
# Physical E8 simple root i maps to this unit-normalized II* chart component.
# The order is fixed by the affine attachment: root eight maps to B1.
CHART_ORDER = ("B1", "B2", "B3", "B4", "N3", "N40", "N4B", "N4inf")
PHYSICAL_TO_CHART = (2, 7, 6, 4, 8, 3, 5, 1)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--child", type=Path, default=CHILD)
parser.add_argument("--target", type=Path, default=TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("child", "target", "output"):
    setattr(args, name, getattr(args, name).resolve())

child = json.loads(args.child.read_text())
target = json.loads(args.target.read_text())
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"

T_ring = PolynomialRing(QQ, "T")
T = T_ring.gen()
model = child["minimal_short_weierstrass"]
A = polynomial(T_ring, model["A_coefficients_low_to_high"])
B = polynomial(T_ring, model["B_coefficients_low_to_high"])
ii_star = next(fibre for fibre in child["finite_fibres"] if fibre["kodaira"] == "II*")
factor = T_ring(ii_star["factor"])
assert factor.degree() == 1 and tuple(ii_star["minimal_orders"]) == (4, 5, 10)
base_point = -factor[0] / factor[1]

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
A_local = u_ring(A(base_point + u))
B_local = u_ring(B(base_point + u))
a, remainder = A_local.quo_rem(u**4)
assert not remainder and a(0)
b, remainder = B_local.quo_rem(u**5)
assert not remainder and b(0)

e8 = target["selected_q8"]["E8"]
physical_cycle = vector(ZZ, e8["vertical_cycle"])
physical_degrees = vector(ZZ, e8["component_degrees"])
assert physical_cycle == vector(ZZ, (4, 5, 7, 10, 8, 6, 4, 2))
assert physical_degrees == vector(ZZ, (-1, 0, 0, 0, 0, 0, 0, 0))
assert matrix(ZZ, e8["cartan"]) == E8_CARTAN
assert physical_cycle * E8_CARTAN == -physical_degrees

permutation = matrix(ZZ, 8, 8, lambda physical, chart:
                     ZZ(PHYSICAL_TO_CHART[physical] - 1 == chart))
chart_cycle = permutation.transpose() * physical_cycle
chart_degrees = permutation.transpose() * physical_degrees
assert chart_cycle == vector(ZZ, (2, 4, 6, 10, 4, 7, 5, 8))
assert chart_degrees == vector(ZZ, (0, -1, 0, 0, 0, 0, 0, 0))

# Unit-normalized II* valuations on the eight components.  The first three
# coordinate divisors meet respectively B1, B2, and N4B; their valuation
# cycles are therefore the corresponding inverse-Cartan columns.
u_value = vector(ZZ, (2, 2, 4, 6, 3, 4, 3, 5))
x_value = vector(ZZ, (2, 4, 6, 10, 4, 7, 5, 8))
y_value = vector(ZZ, (3, 5, 9, 15, 6, 10, 8, 12))
chart_cartan = permutation.transpose() * E8_CARTAN * permutation
assert chart_cartan * u_value == vector(ZZ, (1, 0, 0, 0, 0, 0, 0, 0))
assert chart_cartan * x_value == vector(ZZ, (0, 1, 0, 0, 0, 0, 0, 0))
assert chart_cartan * y_value == vector(ZZ, (0, 0, 0, 0, 0, 0, 1, 0))
assert chart_cycle == x_value

# The valuation semigroup gives the complete ideal without needing a chosen
# global chord trivialization.  The two unit coefficients a(0), b(0) ensure
# the usual II* leading-term equalities in every exceptional valuation.
for index in range(8):
    orders = (
        2 * y_value[index], 3 * x_value[index],
        4 * u_value[index] + x_value[index], 5 * u_value[index],
    )
    assert orders.count(min(orders)) >= 2
minimal_generators = []
for y_exponent in range(2):
    for u_exponent in range(6):
        for x_exponent in range(6):
            valuation = u_exponent * u_value + x_exponent * x_value + y_exponent * y_value
            if min(valuation - chart_cycle) < 0:
                continue
            candidate = (u_exponent, x_exponent, y_exponent)
            if any(
                generator[0] <= u_exponent
                and generator[1] <= x_exponent
                and generator[2] <= y_exponent
                for generator in minimal_generators
            ):
                continue
            minimal_generators.append(candidate)
assert minimal_generators == [(0, 1, 0), (2, 0, 0), (0, 0, 1)]

ring = PolynomialRing(QQ, names=("u", "X", "Y"), order="degrevlex")
uu, X, Y = ring.gens()
relation = Y**2 - X**3 - ring(a(uu)) * uu**4 * X - ring(b(uu)) * uu**5
ideal = ring.ideal((relation, uu**2, X, Y))
assert ideal.vector_space_dimension() == 2
assert ZZ(chart_cycle * chart_cartan * chart_cycle) // 2 == 2

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-iistar-vertical-ideal.v1",
    "status": "PASS_EXACT_Q6_CHILD_Q8_IISTAR_VERTICAL_IDEAL",
    "inputs": {
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "physical_root_target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)},
    },
    "ii_star_local_model": {
        "base": "u=T-(-factor[0]/factor[1])",
        "equation": "Y^2=X^3+u^4*a(u)*X+u^5*b(u)",
        "a0_nonzero": True,
        "b0_nonzero": True,
    },
    "vertical_target": {
        "physical_E8_cycle": list(map(int, physical_cycle)),
        "physical_E8_component_degrees": list(map(int, physical_degrees)),
        "chart_component_order": list(CHART_ORDER),
        "physical_E8_i_to_chart_component": list(PHYSICAL_TO_CHART),
        "chart_cycle": list(map(int, chart_cycle)),
        "chart_component_degrees": list(map(int, chart_degrees)),
    },
    "complete_ideal": {
        "generators": ["Y", "X", "u^2"],
        "quotient_basis": ["1", "u"],
        "colength": 2,
        "cycle_colength": 2,
    },
    "boundary": (
        "This derives only the II* vertical complete ideal selected by the q8 "
        "divisor. A chord-basis trivialization, the IV* module, smooth collision "
        "assembly, a global q8 pencil, a rootless equation, bisection covers, "
        "extension collisions, and generic rank 18 or 19 remain unproved."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8II|ideal=(u2,X,Y)|colength=2|"
    "status=PASS_EXACT_Q6_CHILD_Q8_IISTAR_VERTICAL_IDEAL",
    flush=True,
)
