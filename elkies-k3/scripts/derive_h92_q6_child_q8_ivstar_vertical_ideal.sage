#!/usr/bin/env sage -python
"""Compile the two E6-orbit IV* vertical ideals for the corrected child q=8 target.

At the IV* fibre of the explicit q6 child, write the unit-normalized germ as

    Y^2 = X^3 + u^3*a(u)*X + u^4*b(u),

where ``a(0)b(0)`` is nonzero.  The fibre is split over QQ, so ``b(0)=c^2``.
The standard resolved E6 chart has two arms exchanged by its diagram
involution.  The physical lattice target has component degree
``(-1,-1,0,0,0,0)``; without a chart attachment it can land on either arm.
They give the conjugate complete ideals

    (Y-c*u^2, u*X, X^2, u^3),  and  (Y+c*u^2, u*X, X^2, u^3).

Both have colength four.  This is deliberately an *unoriented* IV* result:
it narrows the additive q8 condition to two explicit finite ideals, but does
not select one of them as the physical q8 module or construct a pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, matrix, prod, vector


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-physical-root-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-vertical-ideal.json"

PHYSICAL_CARTAN = matrix(ZZ, [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, -1],
    [0, 0, -1, 2, -1, 0],
    [0, 0, 0, -1, 2, 0],
    [0, 0, -1, 0, 0, 2],
])
# This is the ordinary-blow-up E6 order.  The first arm is 1--6, and the
# two length-two arms are 2--4--6 and 3--5--6.  Direct chart pullbacks give
# the coordinate valuation cycles below; their Cartan products record where
# the strict transforms of u, X, and Y attach.
CHART_CARTAN = matrix(ZZ, [
    [2, 0, 0, 0, 0, -1],
    [0, 2, 0, -1, 0, 0],
    [0, 0, 2, 0, -1, 0],
    [0, -1, 0, 2, 0, -1],
    [0, 0, -1, 0, 2, -1],
    [-1, 0, 0, -1, -1, 2],
])
CHART_ORDER = ("leaf", "plus_outer", "minus_outer", "plus_inner", "minus_inner", "central")
# The two graph isomorphisms from the pinned physical E6 simple roots to the
# chart order.  They differ by the E6 arm involution and are the only
# remaining ambiguity in the current equation-to-lattice comparison.  In the
# displayed ordinary u-charts, ``plus_outer`` is Y/u^2=c and ``minus_outer``
# is Y/u^2=-c.
PHYSICAL_TO_CHART_CHOICES = (
    (2, 4, 6, 5, 3, 1),
    (3, 5, 6, 4, 2, 1),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


def permutation(mapping):
    return matrix(ZZ, 6, 6, lambda physical, chart:
                  ZZ(mapping[physical] - 1 == chart))


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
iv_star = next(fibre for fibre in child["finite_fibres"] if fibre["kodaira"] == "IV*")
factor = T_ring(iv_star["factor"])
assert factor.degree() == 1 and tuple(iv_star["minimal_orders"]) == (3, 4, 8)
base_point = -factor[0] / factor[1]

u_ring = PolynomialRing(QQ, "u")
u = u_ring.gen()
A_local = u_ring(A(base_point + u))
B_local = u_ring(B(base_point + u))
a, remainder = A_local.quo_rem(u**3)
assert not remainder and a(0)
b, remainder = B_local.quo_rem(u**4)
assert not remainder and b(0)
c = QQ(b(0)).sqrt()
assert c and c**2 == b(0)

# Start the actual ordinary-blow-up resolution rather than reading the two
# arms from the Kodaira symbol.  After the u-chart blow-up and then the
# x-chart blow-up at its sole singular point, the tangent cone is exactly
# (Y-c*u)(Y+c*u).  These are the two arms whose later exceptional curves are
# exchanged by the E6 diagram involution.
local_ring = PolynomialRing(QQ, names=("u", "X", "Y"))
uu_local, X_local, Y_local = local_ring.gens()
a_local = local_ring(a(uu_local))
b_local = local_ring(b(uu_local))
f0 = Y_local**2 - X_local**3 - uu_local**3 * a_local * X_local - uu_local**4 * b_local


def strict(value, substitutions, exceptional, multiplicity):
    transformed = local_ring(value(*substitutions))
    quotient, remainder = transformed.quo_rem(exceptional**multiplicity)
    assert not remainder
    return quotient


f1 = strict(
    f0, (uu_local, uu_local * X_local, uu_local * Y_local), uu_local, 2
)
assert f1 == (
    Y_local**2 - uu_local * X_local**3
    - uu_local**2 * a_local * X_local - uu_local**2 * b_local
)
f2 = strict(
    f1, (X_local * uu_local, X_local, X_local * Y_local), X_local, 2
)
expected_f2 = (
    Y_local**2 - X_local**2 * uu_local
    - X_local * uu_local**2 * local_ring(a(X_local * uu_local))
    - uu_local**2 * local_ring(b(X_local * uu_local))
)
assert f2 == expected_f2
assert f2.degree() >= 2
tangent_cone = sum(
    coefficient * prod(generator**power for generator, power in zip(
        (uu_local, X_local, Y_local), exponent
    ))
    for exponent, coefficient in f2.dict().items()
    if sum(exponent) == 2
)
assert tangent_cone == Y_local**2 - c**2 * uu_local**2
assert tangent_cone == (Y_local - c * uu_local) * (Y_local + c * uu_local)

e6 = target["selected_q8"]["E6"]
physical_cycle = vector(ZZ, e6["vertical_cycle"])
physical_degrees = vector(ZZ, e6["component_degrees"])
assert matrix(ZZ, e6["cartan"]) == PHYSICAL_CARTAN
assert physical_cycle * PHYSICAL_CARTAN == -physical_degrees
nef_target = target.get("normalization", {}).get("representative") == "nef"
if not nef_target:
    assert physical_cycle == vector(ZZ, (3, 5, 6, 4, 2, 3))
    assert physical_degrees == vector(ZZ, (-1, -1, 0, 0, 0, 0))
else:
    assert physical_cycle == vector(ZZ, (2, 3, 4, 3, 2, 2))
    assert physical_degrees == vector(ZZ, (-1, 0, 0, 0, -1, 0))

# The actual ordinary E6 resolution has these coordinate valuations.  They
# can also be recognized intrinsically from the strict-transform attachments:
# C*u=e_leaf, C*X=e_left_outer+e_right_outer, C*Y=e_central.
u_value = vector(ZZ, (2, 1, 1, 2, 2, 3))
x_value = vector(ZZ, (2, 2, 2, 3, 3, 4))
y_value = vector(ZZ, (3, 2, 2, 4, 4, 6))
assert CHART_CARTAN * u_value == vector(ZZ, (1, 0, 0, 0, 0, 0))
assert CHART_CARTAN * x_value == vector(ZZ, (0, 1, 1, 0, 0, 0))
assert CHART_CARTAN * y_value == vector(ZZ, (0, 0, 0, 0, 0, 1))
for index in range(6):
    orders = (
        2 * y_value[index], 3 * x_value[index],
        3 * u_value[index] + x_value[index], 4 * u_value[index],
    )
    assert orders.count(min(orders)) >= 2

orientation_data = []
for sign, mapping in zip(("minus", "plus"), PHYSICAL_TO_CHART_CHOICES):
    change = permutation(mapping)
    chart_cycle = change.transpose() * physical_cycle
    chart_degrees = change.transpose() * physical_degrees
    assert change.transpose() * PHYSICAL_CARTAN * change == CHART_CARTAN
    assert chart_cycle * CHART_CARTAN == -chart_degrees
    assert ZZ(chart_cycle * CHART_CARTAN * chart_cycle) // 2 == (2 if nef_target else 4)

    if nef_target:
        # After the four q8 chamber reflections both arm choices give the
        # same cycle.  Its valuation semigroup is the colength-two ideal
        # (u^2,X,Y), so no residual E6-arm choice remains.
        assert chart_cycle == vector(ZZ, (2, 2, 2, 3, 3, 4))
        assert chart_degrees == vector(ZZ, (0, -1, -1, 0, 0, 0))
        generators = ["Y", "X", "u^2"]
        colength = 2
    else:
        # The first ordinary double-point chart separates the two tangent
        # arms.  Its branch lift is stable under higher u-terms because only
        # b(0)=c^2 enters the tangent cone.
        if sign == "minus":
            assert chart_cycle == vector(ZZ, (3, 3, 2, 5, 4, 6))
            branch_generator = "Y-c*u^2"
        else:
            assert chart_cycle == vector(ZZ, (3, 2, 3, 4, 5, 6))
            branch_generator = "Y+c*u^2"
        generators = [branch_generator, "u*X", "X^2", "u^3"]
        colength = 4
    orientation_data.append({
        "sign": sign,
        "physical_E6_i_to_chart_component": list(mapping),
        "chart_cycle": list(map(int, chart_cycle)),
        "chart_component_degrees": list(map(int, chart_degrees)),
        "generators": generators,
        "cycle_colength": colength,
    })

# Exact colength check in the actual child germ.  The sign does not affect
# it: modulo u^3, u*X, X^2 and Y-sign*c*u^2, the four classes are 1,u,u^2,X.
ring = PolynomialRing(QQ, names=("u", "X", "Y"), order="degrevlex")
uu, X, Y = ring.gens()
relation = Y**2 - X**3 - ring(a(uu)) * uu**3 * X - ring(b(uu)) * uu**4
if nef_target:
    ideal = ring.ideal((relation, uu**2, X, Y))
    assert ideal.vector_space_dimension() == 2
else:
    for sign in (-1, 1):
        ideal = ring.ideal((relation, Y + sign * c * uu**2, uu * X, X**2, uu**3))
        assert ideal.vector_space_dimension() == 4

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-ivstar-vertical-ideal.v1",
    "status": "PASS_EXACT_Q6_CHILD_Q8_IVSTAR_VERTICAL_IDEAL_PAIR",
    "inputs": {
        "child_jacobian": {"path": str(args.child.relative_to(ROOT)), "sha256": digest(args.child)},
        "physical_root_target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)},
    },
    "target_representative": "nef" if nef_target else "dominant-d13",
    "iv_star_local_model": {
        "base": "u=T-(-factor[0]/factor[1])",
        "equation": "Y^2=X^3+u^3*a(u)*X+u^4*b(u)",
        "a0_nonzero": True,
        "b0_square": True,
        "chosen_square_root_c": str(c),
        "first_two_ordinary_blowup_charts": {
            "first_u_chart": str(f1),
            "second_x_chart": str(f2),
            "second_x_chart_tangent_cone": "(Y-c*u)*(Y+c*u)",
        },
    },
    "chart_resolution": {
        "component_order": list(CHART_ORDER),
        "cartan": [[int(value) for value in row] for row in CHART_CARTAN.rows()],
        "coordinate_valuations": {
            "u": list(map(int, u_value)),
            "X": list(map(int, x_value)),
            "Y": list(map(int, y_value)),
        },
    },
    "orientation_candidates": orientation_data,
    "common_quotient_basis": ["1", "u"] if nef_target else ["1", "u", "u^2", "X"],
    "boundary": (
        "The physical E6 lattice roots have not yet been attached to one of the "
        "two resolved chart arms, so this exports an exact conjugate pair rather "
        "than selecting a q8 IV* coefficient module. It does not trivialize the "
        "generic chord, combine the II*, IV*, and smooth conditions, construct a "
        "q8 pencil or rootless equation, produce bisections or extension collisions, "
        "or prove generic rank 18 or 19."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8IV|representative={}|ideals=2|colength={}|orientation={}"
    "|status=PASS_EXACT_Q6_CHILD_Q8_IVSTAR_VERTICAL_IDEAL_PAIR".format(
        "nef" if nef_target else "dominant-d13", 2 if nef_target else 4,
        "none" if nef_target else "unresolved_E6_arm",
    ),
    flush=True,
)
